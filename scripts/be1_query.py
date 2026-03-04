#!/usr/bin/env python3
"""
BE1 Dataset Query Tool
Extract expression profiles from the BE1 single-cell RNA-seq benchmark dataset

Enhanced version with support for comma-separated lists of genes and cells,
proper argparse type validation, and advanced query options.
"""

import scanpy as sc
import pandas as pd
import numpy as np
import json
import os
import sys
import argparse
import requests
import gzip
import shutil
import subprocess
from pathlib import Path
from typing import Optional, Dict, List, Union, Tuple
import urllib.request
from urllib.error import URLError
import tempfile
import time
import hashlib
import logging

# ============================================================================
# Custom argparse type functions for validation
# ============================================================================

def positive_int(value):
    """Ensure value is a positive integer"""
    try:
        ivalue = int(value)
        if ivalue <= 0:
            raise argparse.ArgumentTypeError(f"{value} must be a positive integer")
        return ivalue
    except ValueError:
        raise argparse.ArgumentTypeError(f"{value} must be an integer")

def positive_float(value):
    """Ensure value is a positive float"""
    try:
        fvalue = float(value)
        if fvalue <= 0:
            raise argparse.ArgumentTypeError(f"{value} must be a positive number")
        return fvalue
    except ValueError:
        raise argparse.ArgumentTypeError(f"{value} must be a number")

def non_negative_float(value):
    """Ensure value is a non-negative float"""
    try:
        fvalue = float(value)
        if fvalue < 0:
            raise argparse.ArgumentTypeError(f"{value} must be non-negative")
        return fvalue
    except ValueError:
        raise argparse.ArgumentTypeError(f"{value} must be a number")

def existing_file(value):
    """Ensure file exists"""
    if not os.path.exists(value):
        raise argparse.ArgumentTypeError(f"File '{value}' does not exist")
    return value

def existing_directory(value):
    """Ensure directory exists (create if needed)"""
    path = Path(value)
    path.mkdir(parents=True, exist_ok=True)
    return str(path)

def cell_line_type(value):
    """Validate cell line names with fuzzy matching"""
    valid_cell_lines = ['PC9', 'A549', 'HCC78', 'DV90', 'PBMC', 'PBMCs', 
                        'HTB178', 'CRL5868', 'CCL-185-IG', 'NCI-H596', 
                        'NCI-H1395', 'HTB', 'ALL']
    
    # Allow 'all' as special value
    if value.upper() == 'ALL':
        return 'ALL'
    
    # Case-insensitive exact match
    for valid in valid_cell_lines:
        if value.upper() == valid.upper():
            return valid
    
    # Try partial match
    matches = [cl for cl in valid_cell_lines if cl.upper() in value.upper() or value.upper() in cl.upper()]
    if len(matches) == 1:
        print(f"⚠️  Assuming '{matches[0]}' for cell line '{value}'")
        return matches[0]
    elif len(matches) > 1:
        raise argparse.ArgumentTypeError(
            f"Ambiguous cell line '{value}'. Did you mean one of: {', '.join(matches)}"
        )
    
    raise argparse.ArgumentTypeError(
        f"Invalid cell line: '{value}'. Choose from: {', '.join(valid_cell_lines)}"
    )

def mutation_type(value):
    """Validate mutation type"""
    valid_mutations = [
        'EGFR Del19', 'KRAS G12S', 'BRAF G469A', 'MET Del14',
        'ERBB2 V842I', 'ROS1 fusion', 'EML4-ALK fusion', 'Healthy control'
    ]
    
    # Allow 'all' as special value
    if value.upper() == 'ALL':
        return 'ALL'
    
    # Case-insensitive partial matching
    value_lower = value.lower()
    for valid in valid_mutations:
        if value_lower in valid.lower() or valid.lower() in value_lower:
            return valid
    
    raise argparse.ArgumentTypeError(
        f"Unknown mutation: '{value}'. Valid options: {valid_mutations}"
    )

def file_size_limit(max_mb):
    """Factory function to create file size validators"""
    def validator(value):
        if not os.path.exists(value):
            return value  # Let other validators handle non-existence
        size_mb = os.path.getsize(value) / (1024 * 1024)
        if size_mb > max_mb:
            raise argparse.ArgumentTypeError(
                f"File '{value}' is {size_mb:.1f}MB, exceeds limit of {max_mb}MB"
            )
        return value
    return validator

def percentage(value):
    """Validate percentage value (0-100)"""
    try:
        fvalue = float(value)
        if fvalue < 0 or fvalue > 100:
            raise argparse.ArgumentTypeError(f"{value} must be between 0 and 100")
        return fvalue
    except ValueError:
        raise argparse.ArgumentTypeError(f"{value} must be a number")

def range_type(min_val, max_val):
    """Factory function to create range validators"""
    def range_checker(value):
        try:
            fvalue = float(value)
            if fvalue < min_val or fvalue > max_val:
                raise argparse.ArgumentTypeError(
                    f"Value must be between {min_val} and {max_val}"
                )
            return fvalue
        except ValueError:
            raise argparse.ArgumentTypeError(f"{value} must be a number")
    return range_checker

def comma_separated_list(value):
    """Parse comma-separated string into list"""
    if not value:
        return []
    return [item.strip() for item in value.split(',') if item.strip()]

def gene_symbols(value):
    """Validate and parse gene symbols"""
    genes = comma_separated_list(value)
    # Basic validation: ensure each gene contains only valid characters
    for gene in genes:
        if not gene.replace('-', '').replace('_', '').isalnum():
            print(f"⚠️  Warning: Gene '{gene}' contains unusual characters")
    return genes

def cell_barcodes(value):
    """Validate and parse cell barcodes"""
    return comma_separated_list(value)


class BE1QueryTool:
    """Main class for querying BE1 dataset"""
    
    def __init__(self, geo_id: str, metadata_file: str, data_dir: str = "./be1_data", 
                 verbose: int = 0, log_file: Optional[str] = None):
        """
        Initialize the query tool
        
        Parameters:
        -----------
        geo_id : str
            GEO accession ID (e.g., GSE243665)
        metadata_file : str
            Path to JSON metadata file
        data_dir : str
            Directory to store downloaded data
        verbose : int
            Verbosity level (0-3)
        log_file : str or None
            Path to log file
        """
        self.geo_id = geo_id
        self.metadata_file = metadata_file
        self.data_dir = Path(data_dir)
        self.data_dir.mkdir(parents=True, exist_ok=True)
        self.verbose = verbose
        
        # Setup logging
        self._setup_logging(log_file)
        
        # Load metadata
        self.metadata = self._load_metadata()
        
        # Cell line mapping with correct file prefixes
        self.cell_line_mapping = {
            'PC9': {
                'prefix': 'GSE243665_PC9', 
                'mutation': 'EGFR Del19',
                'cells': 4492
            },
            'A549': {
                'prefix': 'GSE243665_A549', 
                'mutation': 'KRAS G12S',
                'cells': 6898
            },
            'CCL-185-IG': {
                'prefix': 'GSE243665_CCL-185-IG', 
                'mutation': 'EML4-ALK fusion',
                'cells': 6354
            },
            'NCI-H1395': {
                'prefix': 'GSE243665_CRL5868', 
                'mutation': 'BRAF G469A',
                'cells': 2673
            },
            'CRL5868': {
                'prefix': 'GSE243665_CRL5868', 
                'mutation': 'BRAF G469A',
                'cells': 2673
            },
            'DV90': {
                'prefix': 'GSE243665_DV90', 
                'mutation': 'ERBB2 V842I',
                'cells': 2998
            },
            'HCC78': {
                'prefix': 'GSE243665_HCC78', 
                'mutation': 'ROS1 fusion',
                'cells': 2748
            },
            'NCI-H596': {
                'prefix': 'GSE243665_HTB178', 
                'mutation': 'MET Del14',
                'cells': 2965
            },
            'HTB178': {
                'prefix': 'GSE243665_HTB178', 
                'mutation': 'MET Del14',
                'cells': 2965
            },
            'HTB': {
                'prefix': 'GSE243665_HTB178', 
                'mutation': 'MET Del14',
                'cells': 2965
            },
            'PBMC': {
                'prefix': 'GSE243665_PBMCs', 
                'mutation': 'Healthy control',
                'cells': 500
            },
            'PBMCs': {
                'prefix': 'GSE243665_PBMCs', 
                'mutation': 'Healthy control',
                'cells': 500
            }
        }
        
        # File sizes from GEO (in MB) for verification
        self.expected_file_sizes = {
            'GSE243665_HTB178_matrix.mtx.gz': 23.7,
            'GSE243665_HTB178_barcodes.tsv.gz': 0.015,  # ~15KB
            'GSE243665_HTB178_features.tsv.gz': 0.327,  # ~327KB
            'GSE243665_A549_matrix.mtx.gz': 59.7,
            'GSE243665_PC9_matrix.mtx.gz': 45.0,
            'GSE243665_CCL-185-IG_matrix.mtx.gz': 54.0,
            'GSE243665_CRL5868_matrix.mtx.gz': 16.8,
            'GSE243665_DV90_matrix.mtx.gz': 24.7,
            'GSE243665_HCC78_matrix.mtx.gz': 27.1,
            'GSE243665_PBMCs_matrix.mtx.gz': 1.0
        }
        
        # Loaded data cache
        self.loaded_data = {}
    
    def _setup_logging(self, log_file: Optional[str] = None):
        """Setup logging configuration"""
        log_level = logging.WARNING
        if self.verbose == 1:
            log_level = logging.INFO
        elif self.verbose >= 2:
            log_level = logging.DEBUG
        
        logging.basicConfig(
            level=log_level,
            format='%(asctime)s - %(levelname)s - %(message)s',
            filename=log_file,
            filemode='a' if log_file else None
        )
        self.logger = logging.getLogger(__name__)
    
    def _log(self, message: str, level: str = 'info'):
        """Log message with appropriate level"""
        if level == 'debug' and self.verbose >= 2:
            print(f"🔍 {message}")
        elif level == 'info' and self.verbose >= 1:
            print(f"ℹ️  {message}")
        elif level == 'warning':
            print(f"⚠️  {message}")
        elif level == 'error':
            print(f"❌ {message}")
        
        # Also log to file
        getattr(self.logger, level)(message)
    
    def _load_metadata(self) -> dict:
        """Load and parse the JSON metadata file"""
        try:
            with open(self.metadata_file, 'r') as f:
                return json.load(f)
        except FileNotFoundError:
            self._log(f"Metadata file not found: {self.metadata_file}", 'warning')
            self._log("Using default metadata from GEO", 'info')
            return self._get_default_metadata()
        except json.JSONDecodeError:
            self._log(f"Invalid JSON in metadata file: {self.metadata_file}", 'error')
            sys.exit(1)
    
    def _get_default_metadata(self) -> dict:
        """Return default metadata from GEO"""
        return {
            "GSE243665": {
                "study": [{
                    "Study ID": "GSE243665",
                    "Title": "A single cell RNAseq benchmark experiment embedding 'controlled' cancer heterogeneity"
                }],
                "sample": [
                    {"Sample ID": "BE1_PC9", "Cell Line": "PC9", "Driver Mutation": "EGFR Del19"},
                    {"Sample ID": "BE1_A549", "Cell Line": "A549", "Driver Mutation": "KRAS G12S"},
                    {"Sample ID": "BE1_NCI-H596", "Cell Line": "NCI-H596", "Driver Mutation": "MET Del14"},
                    {"Sample ID": "BE1_NCI-H1395", "Cell Line": "NCI-H1395", "Driver Mutation": "BRAF G469A"},
                    {"Sample ID": "BE1_DV90", "Cell Line": "DV90", "Driver Mutation": "ERBB2 V842I"},
                    {"Sample ID": "BE1_HCC78", "Cell Line": "HCC78", "Driver Mutation": "ROS1 fusion"},
                    {"Sample ID": "BE1_CCL185-IG", "Cell Line": "CCL-185-IG", "Driver Mutation": "EML4-ALK fusion"},
                    {"Sample ID": "BE1_PBMC", "Cell Line": "PBMC", "Driver Mutation": "Healthy control"}
                ]
            }
        }
    
    def check_file_integrity(self, file_path: Path) -> tuple:
        """
        Check file integrity by size and gzip validity
        Returns: (is_valid, message)
        """
        if not file_path.exists():
            return False, "File does not exist"
        
        # Check file size
        size_mb = file_path.stat().st_size / (1024 * 1024)
        filename = file_path.name
        
        if filename in self.expected_file_sizes:
            expected = self.expected_file_sizes[filename]
            size_diff = abs(size_mb - expected) / expected
            if size_diff > 0.1:  # 10% tolerance
                return False, f"File size mismatch: got {size_mb:.1f}MB, expected {expected:.1f}MB"
        
        # Check gzip validity
        try:
            with gzip.open(file_path, 'rb') as f:
                # Try to read first 1KB
                f.read(1024)
            return True, f"Valid gzip file ({size_mb:.1f}MB)"
        except Exception as e:
            # Try with zcat -f
            try:
                result = subprocess.run(
                    ['zcat', '-f', str(file_path)],
                    capture_output=True,
                    timeout=10,
                    check=False
                )
                if result.returncode == 0 and len(result.stdout) > 0:
                    return True, f"Repairable with zcat -f ({size_mb:.1f}MB)"
                else:
                    return False, f"Corrupted: {e}"
            except:
                return False, f"Corrupted: {e}"
    
    def check_data_exists(self, cell_line: str) -> bool:
        """Check if data for a specific cell line exists locally"""
        if cell_line not in self.cell_line_mapping:
            return False
        
        prefix = self.cell_line_mapping[cell_line]['prefix']
        
        # Check for all three required files with proper naming
        required_files = [
            f"{prefix}_barcodes.tsv.gz",
            f"{prefix}_features.tsv.gz",
            f"{prefix}_matrix.mtx.gz"
        ]
        
        self._log(f"Checking files for {cell_line}:", 'info')
        all_valid = True
        for file in required_files:
            file_path = self.data_dir / file
            is_valid, message = self.check_file_integrity(file_path)
            
            if is_valid:
                self._log(f"  ✅ {file}: {message}", 'debug')
            else:
                self._log(f"  ❌ {file}: {message}", 'warning')
                all_valid = False
        
        return all_valid
    
    def aggressive_matrix_repair(self, file_path: Path) -> bool:
        """
        Aggressive repair for matrix files using multiple methods
        """
        self._log(f"Attempting aggressive repair for {file_path.name}...", 'info')
        
        backup_path = file_path.with_suffix('.gz.corrupted')
        shutil.copy2(file_path, backup_path)
        
        try:
            # Method 1: Try to fix with gzip -d and recompress
            temp_file = file_path.with_suffix('.tmp')
            
            # Try different decompression methods
            methods = [
                (['gzip', '-d', '-c', '-f', str(file_path)], "gzip -d -c -f"),
                (['zcat', '-f', str(file_path)], "zcat -f"),
                (['gunzip', '-c', '-f', str(file_path)], "gunzip -c -f")
            ]
            
            for cmd, method_name in methods:
                try:
                    result = subprocess.run(
                        cmd,
                        capture_output=True,
                        timeout=30,
                        check=False
                    )
                    
                    if result.returncode == 0 and len(result.stdout) > 0:
                        # Successfully decompressed
                        with open(temp_file, 'wb') as f_out:
                            f_out.write(result.stdout)
                        
                        # Now recompress with proper gzip
                        with gzip.open(file_path, 'wb', compresslevel=6) as f_out:
                            with open(temp_file, 'rb') as f_in:
                                shutil.copyfileobj(f_in, f_out)
                        
                        # Verify repair
                        try:
                            with gzip.open(file_path, 'rb') as f_test:
                                f_test.read(1024)
                            self._log(f"✅ Successfully repaired with {method_name}", 'info')
                            temp_file.unlink(missing_ok=True)
                            backup_path.unlink(missing_ok=True)
                            return True
                        except:
                            continue
                except:
                    continue
            
            # Method 2: Try to read with Python's gzip ignoring errors
            try:
                with open(temp_file, 'wb') as f_out:
                    with gzip.open(file_path, 'rb') as f_in:
                        try:
                            data = f_in.read()
                            f_out.write(data)
                        except:
                            # If error, try to recover partial data
                            pass
                
                # Check if we got any data
                if temp_file.exists() and temp_file.stat().st_size > 0:
                    with gzip.open(file_path, 'wb') as f_out:
                        with open(temp_file, 'rb') as f_in:
                            shutil.copyfileobj(f_in, f_out)
                    
                    # Verify
                    try:
                        with gzip.open(file_path, 'rb') as f_test:
                            f_test.read(1024)
                        self._log(f"✅ Successfully repaired with partial data recovery", 'info')
                        temp_file.unlink(missing_ok=True)
                        backup_path.unlink(missing_ok=True)
                        return True
                    except:
                        pass
            except:
                pass
            
            self._log(f"❌ All repair methods failed for {file_path.name}", 'error')
            # Restore from backup
            shutil.copy2(backup_path, file_path)
            return False
            
        except Exception as e:
            self._log(f"Error during repair: {e}", 'error')
            if backup_path.exists():
                shutil.copy2(backup_path, file_path)
            return False
        finally:
            temp_file.unlink(missing_ok=True)
    
    def verify_and_repair_files(self, cell_line: str) -> bool:
        """
        Verify all files for a cell line and repair if corrupted
        """
        if cell_line not in self.cell_line_mapping:
            return False
        
        prefix = self.cell_line_mapping[cell_line]['prefix']
        files = [
            self.data_dir / f"{prefix}_barcodes.tsv.gz",
            self.data_dir / f"{prefix}_features.tsv.gz",
            self.data_dir / f"{prefix}_matrix.mtx.gz"
        ]
        
        all_valid = True
        for file_path in files:
            if not file_path.exists():
                self._log(f"❌ Missing file: {file_path.name}", 'error')
                all_valid = False
                continue
            
            is_valid, message = self.check_file_integrity(file_path)
            
            if is_valid:
                self._log(f"✅ {file_path.name}: {message}", 'info')
            else:
                self._log(f"⚠️ {file_path.name}: {message}", 'warning')
                
                # Try different repair methods based on file type
                if 'matrix' in file_path.name:
                    if self.aggressive_matrix_repair(file_path):
                        # Verify repair worked
                        is_valid, message = self.check_file_integrity(file_path)
                        if is_valid:
                            self._log(f"✅ {file_path.name}: Repaired successfully", 'info')
                        else:
                            self._log(f"❌ {file_path.name}: Repair failed", 'error')
                            all_valid = False
                    else:
                        all_valid = False
                else:
                    # For non-matrix files, use standard repair
                    if self.repair_gzip_with_zcat(file_path):
                        is_valid, message = self.check_file_integrity(file_path)
                        if is_valid:
                            self._log(f"✅ {file_path.name}: Repaired successfully", 'info')
                        else:
                            self._log(f"❌ {file_path.name}: Repair failed", 'error')
                            all_valid = False
                    else:
                        all_valid = False
        
        return all_valid
    
    def repair_gzip_with_zcat(self, file_path: Path) -> bool:
        """
        Repair corrupted gzip file using zcat -f
        """
        self._log(f"Repairing {file_path.name} with zcat -f...", 'info')
        
        backup_path = file_path.with_suffix('.gz.backup')
        shutil.copy2(file_path, backup_path)
        
        try:
            temp_file = file_path.with_suffix('.tmp')
            
            result = subprocess.run(
                ['zcat', '-f', str(file_path)],
                capture_output=True,
                timeout=30,
                check=False
            )
            
            if result.returncode == 0 and len(result.stdout) > 0:
                with gzip.open(file_path, 'wb', compresslevel=6) as f_out:
                    f_out.write(result.stdout)
                
                # Verify
                try:
                    with gzip.open(file_path, 'rb') as f_test:
                        f_test.read(100)
                    self._log(f"✅ Successfully repaired {file_path.name}", 'info')
                    temp_file.unlink(missing_ok=True)
                    backup_path.unlink(missing_ok=True)
                    return True
                except:
                    pass
            
            self._log(f"❌ Failed to repair {file_path.name}", 'error')
            shutil.copy2(backup_path, file_path)
            return False
            
        except Exception as e:
            self._log(f"Error during repair: {e}", 'error')
            if backup_path.exists():
                shutil.copy2(backup_path, file_path)
            return False
        finally:
            temp_file.unlink(missing_ok=True)
    
    def force_decompress_files(self, cell_line: str, output_dir: Path) -> bool:
        """
        Force decompress all files for a cell line using multiple methods
        """
        if cell_line not in self.cell_line_mapping:
            return False
        
        prefix = self.cell_line_mapping[cell_line]['prefix']
        files = [
            ('barcodes.tsv.gz', 'barcodes.tsv'),
            ('features.tsv.gz', 'features.tsv'),
            ('matrix.mtx.gz', 'matrix.mtx')
        ]
        
        success = True
        for gz_name, out_name in files:
            gz_file = self.data_dir / f"{prefix}_{gz_name}"
            out_file = output_dir / out_name
            
            if not gz_file.exists():
                self._log(f"❌ Missing file: {gz_file.name}", 'error')
                success = False
                continue
            
            self._log(f"Decompressing {gz_file.name} -> {out_name}", 'debug')
            
            # Try multiple decompression methods
            decompressed = False
            
            # Method 1: zcat -f
            try:
                with open(out_file, 'w') as f_out:
                    result = subprocess.run(
                        ['zcat', '-f', str(gz_file)],
                        capture_output=True,
                        text=True,
                        timeout=30,
                        check=False
                    )
                    
                    if result.returncode == 0 and len(result.stdout) > 0:
                        f_out.write(result.stdout)
                        self._log(f"✅ Decompressed with zcat -f", 'debug')
                        decompressed = True
            except:
                pass
            
            # Method 2: Python gzip
            if not decompressed:
                try:
                    with gzip.open(gz_file, 'rt') as f_in:
                        with open(out_file, 'w') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                    self._log(f"✅ Decompressed with Python gzip", 'debug')
                    decompressed = True
                except Exception as e:
                    self._log(f"❌ Python gzip failed: {e}", 'debug')
            
            # Method 3: For matrix files, try binary mode
            if not decompressed and 'matrix' in gz_name:
                try:
                    result = subprocess.run(
                        ['zcat', '-f', str(gz_file)],
                        capture_output=True,
                        timeout=30,
                        check=False
                    )
                    if result.returncode == 0 and len(result.stdout) > 0:
                        with open(out_file, 'wb') as f_out:
                            f_out.write(result.stdout)
                        self._log(f"✅ Decompressed with zcat -f (binary)", 'debug')
                        decompressed = True
                except:
                    pass
            
            if not decompressed:
                self._log(f"❌ Failed to decompress {gz_file.name}", 'error')
                success = False
        
        return success
    
    def list_available_files(self):
        """List all downloaded files in the data directory with detailed status"""
        print(f"\n📁 Files in {self.data_dir}:")
        files = list(self.data_dir.glob("*.gz"))
        if files:
            for f in sorted(files):
                size_mb = f.stat().st_size / (1024 * 1024)
                is_valid, message = self.check_file_integrity(f)
                
                if is_valid:
                    status = "✅"
                else:
                    if "Repairable" in message:
                        status = "🔄"
                    else:
                        status = "❌"
                
                print(f"  {status} {f.name} ({size_mb:.1f} MB) - {message}")
        else:
            print("  No files found")
    
    def download_cell_line_data(self, cell_line: str, force: bool = False, 
                               timeout: int = 60) -> bool:
        """
        Download data for a specific cell line from GEO
        """
        if cell_line not in self.cell_line_mapping:
            self._log(f"Unknown cell line: {cell_line}", 'error')
            print(f"Available: {list(self.cell_line_mapping.keys())}")
            return False
        
        if not force and self.check_data_exists(cell_line):
            self._log(f"Data for {cell_line} already exists locally", 'info')
            return True
        
        prefix = self.cell_line_mapping[cell_line]['prefix']
        
        # Try multiple download sources
        download_sources = [
            {
                'name': 'GEO HTTPS',
                'base_url': "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE243nnn/GSE243665/suppl/",
                'files': [
                    f"{prefix}_barcodes.tsv.gz",
                    f"{prefix}_features.tsv.gz",
                    f"{prefix}_matrix.mtx.gz"
                ]
            },
            {
                'name': 'GEO FTP',
                'base_url': "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE243nnn/GSE243665/suppl/",
                'files': [
                    f"{prefix}_barcodes.tsv.gz",
                    f"{prefix}_features.tsv.gz",
                    f"{prefix}_matrix.mtx.gz"
                ]
            }
        ]
        
        print(f"\n📥 Downloading {cell_line} data...")
        success = True
        
        for source in download_sources:
            print(f"\nTrying {source['name']}...")
            source_success = True
            
            for file in source['files']:
                url = source['base_url'] + file
                output_file = self.data_dir / file
                
                if output_file.exists() and not force:
                    print(f"  ⏭️  {file} already exists, skipping")
                    continue
                
                print(f"  Downloading {file}...")
                try:
                    # Try with requests
                    response = requests.get(url.replace('ftp://', 'https://'), 
                                          stream=True, 
                                          timeout=timeout,
                                          headers={'User-Agent': 'Mozilla/5.0'})
                    response.raise_for_status()
                    
                    # Download with progress
                    total_size = int(response.headers.get('content-length', 0))
                    downloaded = 0
                    
                    with open(output_file, 'wb') as f:
                        for chunk in response.iter_content(chunk_size=8192):
                            if chunk:
                                f.write(chunk)
                                downloaded += len(chunk)
                                if total_size and downloaded % (1024*1024) < 8192:
                                    progress = (downloaded / total_size) * 100
                                    print(f"\r    Downloaded {downloaded/(1024*1024):.1f} MB / {total_size/(1024*1024):.1f} MB ({progress:.1f}%)", end='')
                    
                    print(f"\n    ✅ Downloaded ({self._get_file_size(output_file):.1f} MB)")
                    
                except Exception as e:
                    print(f"    ❌ Failed: {e}")
                    source_success = False
            
            if source_success:
                print(f"✅ Successfully downloaded from {source['name']}")
                success = True
                break
            else:
                print(f"⚠️  Failed to download from {source['name']}")
        
        if success:
            print(f"✅ Successfully downloaded {cell_line} data")
            print(f"\n🔧 Verifying downloaded files...")
            self.verify_and_repair_files(cell_line)
        else:
            print(f"❌ Failed to download {cell_line} data from all sources")
            print("\nPlease download manually from GEO:")
            print("  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243665")
            print(f"and place files in: {self.data_dir}")
        
        return success
    
    def _get_file_size(self, filepath: Path) -> float:
        """Get file size in MB"""
        return filepath.stat().st_size / (1024 * 1024)
    
    def load_cell_line(self, cell_line: str, force_reload: bool = False) -> Optional[sc.AnnData]:
        """
        Load data for a specific cell line
        """
        if cell_line not in self.cell_line_mapping:
            matches = [cl for cl in self.cell_line_mapping.keys() 
                      if cell_line.upper() in cl.upper()]
            if matches:
                cell_line = matches[0]
                self._log(f"Using cell line: {cell_line}", 'info')
            else:
                self._log(f"Unknown cell line: {cell_line}", 'error')
                print(f"Available: {list(self.cell_line_mapping.keys())}")
                return None
        
        if not force_reload and cell_line in self.loaded_data:
            self._log(f"Using cached data for {cell_line}", 'info')
            return self.loaded_data[cell_line]
        
        # Check if data exists and is valid
        if not self.check_data_exists(cell_line):
            self._log(f"Data for {cell_line} not found or corrupted.", 'warning')
            self.list_available_files()
            
            response = input(f"\nDownload/repair data for {cell_line}? (y/n): ")
            if response.lower() == 'y':
                if not self.download_cell_line_data(cell_line):
                    return None
            else:
                return None
        
        # Verify and repair files before loading
        self._log(f"Verifying files for {cell_line}...", 'info')
        if not self.verify_and_repair_files(cell_line):
            self._log(f"Files for {cell_line} are corrupted and could not be repaired", 'error')
            response = input("Try downloading again? (y/n): ")
            if response.lower() == 'y':
                if not self.download_cell_line_data(cell_line, force=True):
                    return None
            else:
                return None
        
        prefix = self.cell_line_mapping[cell_line]['prefix']
        mutation = self.cell_line_mapping[cell_line]['mutation']
        
        print(f"\n📊 Loading {cell_line}...")
        
        # Try multiple loading methods
        loading_methods = [
            self._load_with_scanpy,
            self._load_with_decompressed_files,
            self._load_with_manual_construction
        ]
        
        for method in loading_methods:
            adata = method(cell_line, prefix, mutation)
            if adata is not None:
                self.loaded_data[cell_line] = adata
                return adata
        
        self._log(f"All loading methods failed for {cell_line}", 'error')
        return None
    
    def _load_with_scanpy(self, cell_line: str, prefix: str, mutation: str) -> Optional[sc.AnnData]:
        """Load using scanpy's built-in reader"""
        try:
            print(f"  Method 1: Using scanpy read_10x_mtx...")
            adata = sc.read_10x_mtx(
                str(self.data_dir),
                gex_only=False,
                prefix=prefix + '_',
                var_names='gene_symbols'
            )
            
            adata.obs['cell_line'] = cell_line
            adata.obs['driver_mutation'] = mutation
            adata.obs['geo_accession'] = self.geo_id
            
            print(f"  ✅ Loaded {adata.n_obs} cells, {adata.n_vars} genes")
            
            expected = self.cell_line_mapping[cell_line]['cells']
            if abs(adata.n_obs - expected) / expected > 0.1:
                print(f"  ⚠️  Warning: Expected ~{expected} cells, found {adata.n_obs}")
            
            return adata
            
        except Exception as e:
            print(f"  ❌ Scanpy loading failed: {e}")
            return None
    
    def _load_with_decompressed_files(self, cell_line: str, prefix: str, mutation: str) -> Optional[sc.AnnData]:
        """Load by decompressing files first"""
        try:
            print(f"  Method 2: Loading from decompressed files...")
            
            with tempfile.TemporaryDirectory() as temp_dir:
                temp_path = Path(temp_dir)
                
                if not self.force_decompress_files(cell_line, temp_path):
                    print(f"  ❌ Failed to decompress files")
                    return None
                
                # Check if matrix file was decompressed
                matrix_file = temp_path / "matrix.mtx"
                if not matrix_file.exists():
                    print(f"  ❌ Matrix file not decompressed")
                    return None
                
                # Try to read with scanpy
                adata = sc.read_10x_mtx(
                    str(temp_path),
                    gex_only=False,
                    var_names='gene_symbols'
                )
                
                adata.obs['cell_line'] = cell_line
                adata.obs['driver_mutation'] = mutation
                adata.obs['geo_accession'] = self.geo_id
                
                print(f"  ✅ Loaded {adata.n_obs} cells from decompressed files")
                return adata
                
        except Exception as e:
            print(f"  ❌ Decompressed loading failed: {e}")
            return None
    
    def _load_with_manual_construction(self, cell_line: str, prefix: str, mutation: str) -> Optional[sc.AnnData]:
        """Manually construct AnnData from files"""
        try:
            print(f"  Method 3: Manual construction...")
            
            # Create temp dir for decompressed files
            with tempfile.TemporaryDirectory() as temp_dir:
                temp_path = Path(temp_dir)
                
                # Decompress files
                files_to_decompress = [
                    (f"{prefix}_barcodes.tsv.gz", "barcodes.tsv"),
                    (f"{prefix}_features.tsv.gz", "features.tsv"),
                    (f"{prefix}_matrix.mtx.gz", "matrix.mtx")
                ]
                
                for gz_name, out_name in files_to_decompress:
                    gz_file = self.data_dir / gz_name
                    out_file = temp_path / out_name
                    
                    # Try multiple methods to decompress matrix
                    if 'matrix' in gz_name:
                        # For matrix file, try binary decompression
                        try:
                            result = subprocess.run(
                                ['zcat', '-f', str(gz_file)],
                                capture_output=True,
                                timeout=60,
                                check=False
                            )
                            if result.returncode == 0 and len(result.stdout) > 0:
                                with open(out_file, 'wb') as f:
                                    f.write(result.stdout)
                                continue
                        except:
                            pass
                    
                    # For other files or if binary method failed
                    try:
                        with gzip.open(gz_file, 'rt') as f_in:
                            with open(out_file, 'w') as f_out:
                                shutil.copyfileobj(f_in, f_out)
                    except Exception as e:
                        print(f"  ⚠️  Failed to decompress {gz_name}: {e}")
                        return None
                
                # Now try to read with scanpy
                adata = sc.read_10x_mtx(
                    str(temp_path),
                    gex_only=False,
                    var_names='gene_symbols'
                )
                
                adata.obs['cell_line'] = cell_line
                adata.obs['driver_mutation'] = mutation
                adata.obs['geo_accession'] = self.geo_id
                
                print(f"  ✅ Loaded {adata.n_obs} cells via manual construction")
                return adata
                
        except Exception as e:
            print(f"  ❌ Manual construction failed: {e}")
            return None
    
    def parse_comma_list(self, input_str: Optional[str]) -> List[str]:
        """
        Parse a comma-separated string into a list
        Handles spaces and empty strings
        """
        if not input_str:
            return []
        
        # Split by comma and strip whitespace
        items = [item.strip() for item in input_str.split(',')]
        # Remove empty strings
        return [item for item in items if item]
    
    def query_genes(self, genes: List[str], cell_line: Optional[str] = None,
                   mutation: Optional[str] = None, min_umi: int = 1,
                   include_zeros: bool = True, normalize: str = 'log1p') -> Optional[pd.DataFrame]:
        """
        Query expression for multiple genes
        
        Parameters:
        -----------
        genes : List[str]
            List of gene symbols to query
        cell_line : str or None
            Filter by cell line
        mutation : str or None
            Filter by mutation type
        min_umi : int
            Minimum UMI count to consider gene expressed
        include_zeros : bool
            Include cells with zero expression
        normalize : str
            Normalization method ('none', 'log1p', 'cpm')
            
        Returns:
        --------
        DataFrame with expression values for all genes (wide format)
        """
        # Load appropriate data
        if cell_line and cell_line != 'ALL':
            adata = self.load_cell_line(cell_line)
            if adata is None:
                return None
        else:
            adata = self.load_all_cell_lines()
            if adata is None:
                return None
        
        # Filter by mutation if specified
        if mutation and mutation != 'ALL' and not cell_line:
            mutation_cells = adata.obs['driver_mutation'] == mutation
            if mutation_cells.sum() == 0:
                print(f"❌ No cells found with mutation: {mutation}")
                return None
            adata = adata[mutation_cells].copy()
            print(f"  Filtered to {adata.n_obs} cells with mutation {mutation}")
        
        # Find each gene
        found_genes = []
        missing_genes = []
        gene_indices = []
        
        for gene in genes:
            if gene in adata.var_names:
                found_genes.append(gene)
                gene_indices.append(list(adata.var_names).index(gene))
            else:
                # Try partial match
                matches = [g for g in adata.var_names if gene.upper() in g.upper()]
                if matches:
                    found_genes.append(matches[0])
                    gene_indices.append(list(adata.var_names).index(matches[0]))
                    print(f"ℹ️  Using {matches[0]} for gene '{gene}'")
                else:
                    missing_genes.append(gene)
        
        if missing_genes:
            print(f"⚠️  Genes not found: {', '.join(missing_genes)}")
        
        if not found_genes:
            print("❌ No valid genes found")
            return None
        
        # Extract expression for all genes
        expression_data = {}
        for gene, idx in zip(found_genes, gene_indices):
            expr = adata.X[:, idx]
            if hasattr(expr, 'toarray'):
                expr = expr.toarray().flatten()
            expression_data[gene] = expr
        
        # Create results DataFrame
        results = pd.DataFrame({
            'cell_barcode': adata.obs_names,
            'cell_line': adata.obs['cell_line'],
            'driver_mutation': adata.obs['driver_mutation']
        })
        
        # Add expression columns with normalization
        for gene in found_genes:
            expr = expression_data[gene]
            results[gene] = expr
            
            if normalize == 'log1p':
                results[f'{gene}_log1p'] = np.log1p(expr)
            elif normalize == 'cpm':
                # Counts per million (crude approximation)
                total_counts = adata.X.sum(axis=1)
                if hasattr(total_counts, 'toarray'):
                    total_counts = total_counts.toarray().flatten()
                cpm = (expr / total_counts) * 1e6
                cpm = np.nan_to_num(cpm, nan=0.0, posinf=0.0, neginf=0.0)
                results[f'{gene}_cpm'] = cpm
            
            results[f'{gene}_expressed'] = expr >= min_umi
        
        # Filter zeros if requested
        if not include_zeros:
            # Keep cells where at least one gene is expressed
            expressed_mask = results[[f'{gene}_expressed' for gene in found_genes]].any(axis=1)
            results = results[expressed_mask].copy()
        
        # Print summary
        print(f"\n{'='*60}")
        print(f"📈 Expression profiles for {len(found_genes)} genes")
        print(f"{'='*60}")
        print(f"Total cells: {len(results)}")
        
        for gene in found_genes:
            expressed = (results[gene] >= min_umi).sum()
            pct_expressed = (expressed / len(results)) * 100 if len(results) > 0 else 0
            mean_expr = results[gene].mean()
            print(f"\n{gene}:")
            print(f"  Cells expressing (≥{min_umi} UMIs): {expressed} ({pct_expressed:.1f}%)")
            print(f"  Mean expression: {mean_expr:.3f}")
            print(f"  Max expression: {results[gene].max():.3f}")
        
        if cell_line and cell_line != 'ALL':
            print(f"\nCell line: {cell_line}")
        if mutation and mutation != 'ALL':
            print(f"Mutation: {mutation}")
        
        return results
    
    def query_cells(self, cell_barcodes: List[str], cell_line: Optional[str] = None,
                   top_n: int = 20) -> Optional[Dict]:
        """
        Query expression for multiple cells
        
        Parameters:
        -----------
        cell_barcodes : List[str]
            List of cell barcodes to query (can be partial)
        cell_line : str or None
            Specific cell line to search in
        top_n : int
            Number of top genes to return per cell
            
        Returns:
        --------
        Dictionary with cell information for all matching cells
        """
        # Determine which cell lines to search
        if cell_line and cell_line != 'ALL':
            cell_lines_to_search = [cell_line]
        else:
            cell_lines_to_search = list(self.cell_line_mapping.keys())
        
        all_results = {}
        
        for cl in cell_lines_to_search:
            adata = self.load_cell_line(cl)
            if adata is None:
                continue
            
            for barcode_pattern in cell_barcodes:
                matching_cells = [cell for cell in adata.obs_names if barcode_pattern in cell]
                
                if matching_cells:
                    for cell in matching_cells:
                        cell_idx = list(adata.obs_names).index(cell)
                        
                        expression = adata.X[cell_idx, :]
                        if hasattr(expression, 'toarray'):
                            expression = expression.toarray().flatten()
                        
                        gene_expr = pd.DataFrame({
                            'gene': adata.var_names,
                            'expression': expression
                        }).sort_values('expression', ascending=False)
                        
                        top_genes = gene_expr[gene_expr['expression'] > 0].head(top_n)
                        
                        all_results[cell] = {
                            'cell_line': cl,
                            'mutation': adata.obs['driver_mutation'].iloc[cell_idx],
                            'total_umis': float(expression.sum()),
                            'genes_detected': int((expression > 0).sum()),
                            'top_genes': top_genes[['gene', 'expression']].to_dict('records'),
                            'search_pattern': barcode_pattern
                        }
        
        if not all_results:
            print(f"❌ No cells found matching any of the patterns: {', '.join(cell_barcodes)}")
            return None
        
        # Group results by search pattern
        print(f"\n{'='*60}")
        print(f"🔬 Found {len(all_results)} matching cells for {len(cell_barcodes)} patterns")
        print(f"{'='*60}")
        
        for pattern in cell_barcodes:
            pattern_cells = {k: v for k, v in all_results.items() if v['search_pattern'] == pattern}
            if pattern_cells:
                print(f"\nPattern '{pattern}': {len(pattern_cells)} cells")
                for cell, info in list(pattern_cells.items())[:3]:  # Show first 3
                    print(f"  {cell}: {info['cell_line']}, {info['total_umis']:.0f} UMIs")
                if len(pattern_cells) > 3:
                    print(f"  ... and {len(pattern_cells)-3} more")
        
        return all_results
    
    def query_genes_in_cells(self, genes: List[str], cell_barcodes: List[str],
                            cell_line: Optional[str] = None, min_umi: int = 1) -> Optional[pd.DataFrame]:
        """
        Get expression of multiple genes in multiple cells
        
        Parameters:
        -----------
        genes : List[str]
            List of gene symbols
        cell_barcodes : List[str]
            List of cell barcodes
        cell_line : str or None
            Specific cell line to search in
        min_umi : int
            Minimum UMI count to consider gene expressed
            
        Returns:
        --------
        DataFrame with expression values for each gene in each cell
        """
        # First find all matching cells
        cells_info = self.query_cells(cell_barcodes, cell_line)
        
        if not cells_info:
            return None
        
        # Load data for each unique cell line
        cell_line_data = {}
        for cell, info in cells_info.items():
            cl = info['cell_line']
            if cl not in cell_line_data:
                adata = self.load_cell_line(cl)
                if adata is not None:
                    cell_line_data[cl] = adata
        
        # Collect results
        results = []
        found_genes = {}
        
        for cell, info in cells_info.items():
            adata = cell_line_data.get(info['cell_line'])
            if adata is None:
                continue
            
            if cell not in adata.obs_names:
                continue
            
            cell_idx = list(adata.obs_names).index(cell)
            
            row = {
                'cell_barcode': cell,
                'cell_line': info['cell_line'],
                'mutation': info['mutation'],
                'total_umis': info['total_umis'],
                'genes_detected': info['genes_detected']
            }
            
            for gene in genes:
                # Find gene if not already mapped
                if gene not in found_genes:
                    if gene in adata.var_names:
                        found_genes[gene] = gene
                    else:
                        matches = [g for g in adata.var_names if gene.upper() in g.upper()]
                        if matches:
                            found_genes[gene] = matches[0]
                            if gene != matches[0]:
                                print(f"ℹ️  Using {matches[0]} for gene '{gene}'")
                        else:
                            found_genes[gene] = None
                
                target_gene = found_genes.get(gene)
                if target_gene and target_gene in adata.var_names:
                    gene_idx = list(adata.var_names).index(target_gene)
                    expr = adata.X[cell_idx, gene_idx]
                    if hasattr(expr, 'toarray'):
                        expr = expr.toarray().flatten()[0]
                    
                    row[f'{gene}_expression'] = float(expr)
                    row[f'{gene}_log1p'] = float(np.log1p(expr))
                    row[f'{gene}_present'] = expr >= min_umi
                else:
                    row[f'{gene}_expression'] = 0.0
                    row[f'{gene}_log1p'] = 0.0
                    row[f'{gene}_present'] = False
            
            results.append(row)
        
        if not results:
            print("❌ No results found")
            return None
        
        results_df = pd.DataFrame(results)
        
        # Print summary
        print(f"\n{'='*60}")
        print(f"🎯 Expression query: {len(genes)} genes in {len(results)} cells")
        print(f"{'='*60}")
        
        for gene in genes:
            present_count = results_df[f'{gene}_present'].sum()
            print(f"{gene}: present in {present_count}/{len(results_df)} cells (≥{min_umi} UMIs)")
        
        return results_df
    
    def load_all_cell_lines(self, cell_lines: Optional[List[str]] = None, 
                           batch_size: int = 5000) -> Optional[sc.AnnData]:
        """Load multiple cell lines and combine them"""
        if cell_lines is None:
            cell_lines = list(self.cell_line_mapping.keys())
        
        cell_lines = list(dict.fromkeys(cell_lines))
        
        print(f"\n📊 Loading {len(cell_lines)} cell lines...")
        adatas = []
        
        for cell_line in cell_lines:
            adata = self.load_cell_line(cell_line)
            if adata is not None:
                adatas.append(adata)
        
        if not adatas:
            print("❌ No data loaded")
            return None
        
        if len(adatas) == 1:
            return adatas[0]
        
        print(f"\n📊 Combining {len(adatas)} cell lines (batch size: {batch_size})...")
        combined = adatas[0].concatenate(
            adatas[1:],
            batch_key='batch',
            batch_categories=[adata.obs['cell_line'].iloc[0] for adata in adatas]
        )
        print(f"  ✅ Total cells: {combined.n_obs}")
        
        return combined


def main():
    parser = argparse.ArgumentParser(
        description='Query BE1 single-cell RNA-seq dataset',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic queries
  python be1_query.py --geo GSE243665 --metadata be1_metadata.json --gene EGFR,KRAS,BRAF
  python be1_query.py --geo GSE243665 --metadata be1_metadata.json --cell AAACCCAAGAGG-1
  
  # With validation and thresholds
  python be1_query.py --geo GSE243665 --metadata be1_metadata.json --gene EGFR --min-umi 2 --threshold-type absolute
  
  # With filtering
  python be1_query.py --geo GSE243665 --metadata be1_metadata.json --gene EGFR,KRAS --cell-line A549 --mutation "EGFR Del19"
  
  # Advanced options
  python be1_query.py --geo GSE243665 --metadata be1_metadata.json --gene EGFR --top-n 50 --verbose --format json --output results.json
  
  # Batch processing
  python be1_query.py --geo GSE243665 --metadata be1_metadata.json --gene EGFR,KRAS --output results.csv --batch-size 1000 --no-cache
  
  # Multiple genes in multiple cells
  python be1_query.py --geo GSE243665 --metadata be1_metadata.json --genes-in-cells EGFR,KRAS --cells AAACCCA,AAACCCG --format tsv
  
Utility commands:
  python be1_query.py --geo GSE243665 --metadata be1_metadata.json --list-files
  python be1_query.py --geo GSE243665 --metadata be1_metadata.json --check-data A549
  python be1_query.py --geo GSE243665 --metadata be1_metadata.json --repair-files HTB178
  python be1_query.py --geo GSE243665 --metadata be1_metadata.json --verify-all
  python be1_query.py --geo GSE243665 --metadata be1_metadata.json --force-download A549 --timeout 120
        """
    )
    
    # Required arguments with validation
    parser.add_argument('--geo', default='GSE243665', 
                       help='GEO accession ID (e.g., GSE243665)')
    
    parser.add_argument('--metadata', default='be1_metadata.json', 
                       type=existing_file,
                       help='JSON metadata file path')
    
    parser.add_argument('--data-dir', default='./be1_data',
                       type=existing_directory,
                       help='Directory for data storage (will be created if needed)')
    
    # Query options with comma-separated lists
    parser.add_argument('--gene', type=gene_symbols,
                       help='Gene symbol(s) to query (comma-separated, e.g., EGFR,KRAS,BRAF)')
    
    parser.add_argument('--cell', type=cell_barcodes,
                       help='Cell barcode(s) to query (comma-separated, can be partial)')
    
    parser.add_argument('--genes-in-cells', nargs=2, metavar=('GENES', 'CELLS'),
                       type=lambda x: gene_symbols(x) if x else x,
                       help='Get expression of GENES (comma-separated) in CELLS (comma-separated)')
    
    parser.add_argument('--cell-line', type=cell_line_type, default='ALL',
                       help='Filter by cell line (PC9, A549, etc., or ALL)')
    
    parser.add_argument('--mutation', type=mutation_type, default='ALL',
                       help='Filter by mutation type')
    
    # Numerical parameters with validation
    parser.add_argument('--min-umi', type=positive_int, default=1,
                       help='Minimum UMI count to consider gene expressed (default: 1)')
    
    parser.add_argument('--threshold', type=non_negative_float, default=0.0,
                       help='Expression threshold value (default: 0.0)')
    
    parser.add_argument('--threshold-type', 
                       choices=['absolute', 'percentile', 'relative', 'zscore'],
                       default='absolute',
                       help='Type of expression threshold (default: absolute)')
    
    parser.add_argument('--percentile', type=range_type(0, 100), default=90,
                       help='Percentile threshold for expression (0-100, default: 90)')
    
    parser.add_argument('--top-n', type=positive_int, default=20,
                       help='Number of top genes to show per cell (default: 20)')
    
    parser.add_argument('--batch-size', type=positive_int, default=5000,
                       help='Batch size for processing large datasets (default: 5000)')
    
    parser.add_argument('--timeout', type=positive_int, default=60,
                       help='Download timeout in seconds (default: 60)')
    
    # Advanced options
    parser.add_argument('--normalization', 
                       choices=['none', 'log1p', 'cpm'],
                       default='log1p',
                       help='Normalization method (default: log1p)')
    
    parser.add_argument('--format',
                       choices=['csv', 'tsv', 'json'],
                       default='csv',
                       help='Output format (default: csv)')
    
    parser.add_argument('--compression',
                       choices=['none', 'gzip', 'bz2'],
                       default='none',
                       help='Output compression (default: none)')
    
    # Flags
    parser.add_argument('--verbose', '-v', action='count', default=0,
                       help='Increase verbosity (-v, -vv, -vvv)')
    
    parser.add_argument('--include-zeros', action='store_true',
                       help='Include zero expression values in output')
    
    parser.add_argument('--force', action='store_true',
                       help='Force operations (download, repair, etc.)')
    
    parser.add_argument('--no-cache', action='store_true',
                       help='Disable data caching')
    
    parser.add_argument('--dry-run', action='store_true',
                       help='Show what would be done without actually doing it')
    
    # Utility options
    parser.add_argument('--list-files', action='store_true',
                       help='List downloaded files with status')
    
    parser.add_argument('--check-data', type=cell_line_type,
                       help='Check if data exists for a cell line')
    
    parser.add_argument('--repair-files', type=cell_line_type,
                       help='Attempt to repair corrupted files for a cell line')
    
    parser.add_argument('--force-download', type=cell_line_type,
                       help='Force re-download data for a cell line')
    
    parser.add_argument('--verify-all', action='store_true',
                       help='Verify all downloaded files')
    
    # Output options
    parser.add_argument('--output', '-o',
                       help='Save results to file')
    
    parser.add_argument('--log-file',
                       help='Log file path')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Validate that at least one query option is provided
    if not any([args.gene, args.cell, args.genes_in_cells, 
                args.list_files, args.check_data, args.repair_files, 
                args.force_download, args.verify_all]):
        parser.print_help()
        print("\n❌ Please specify a query type (--gene, --cell, --genes-in-cells, or a utility command)")
        sys.exit(1)
    
    # Initialize query tool with validation
    try:
        qt = BE1QueryTool(args.geo, args.metadata, args.data_dir, 
                         verbose=args.verbose, log_file=args.log_file)
        
        # Override cache setting if specified
        if args.no_cache:
            qt.loaded_data = {}
        
    except Exception as e:
        print(f"❌ Failed to initialize query tool: {e}")
        sys.exit(1)
    
    # Utility commands
    if args.list_files:
        qt.list_available_files()
        return
    
    if args.check_data:
        if qt.check_data_exists(args.check_data):
            print(f"✅ Data for {args.check_data} exists and is valid")
        else:
            print(f"❌ Data for {args.check_data} is missing or corrupted")
        return
    
    if args.repair_files:
        if args.dry_run:
            print(f"🔍 [DRY RUN] Would attempt to repair files for {args.repair_files}")
            return
        print(f"🔧 Attempting to repair files for {args.repair_files}...")
        if qt.verify_and_repair_files(args.repair_files):
            print(f"✅ Successfully repaired files for {args.repair_files}")
        else:
            print(f"❌ Could not repair all files for {args.repair_files}")
        return
    
    if args.force_download:
        if args.dry_run:
            print(f"🔍 [DRY RUN] Would force download data for {args.force_download}")
            return
        print(f"📥 Force downloading data for {args.force_download}...")
        if qt.download_cell_line_data(args.force_download, force=True, timeout=args.timeout):
            print(f"✅ Successfully downloaded {args.force_download}")
        else:
            print(f"❌ Failed to download {args.force_download}")
        return
    
    if args.verify_all:
        print("🔍 Verifying all cell lines...")
        all_valid = True
        for cell_line in qt.cell_line_mapping.keys():
            if qt.check_data_exists(cell_line):
                print(f"  ✅ {cell_line}")
            else:
                print(f"  ❌ {cell_line}")
                all_valid = False
        print(f"\nOverall: {'✅ All valid' if all_valid else '❌ Some files need repair'}")
        return
    
    # Execute queries based on arguments
    results = None
    start_time = time.time()
    
    if args.genes_in_cells:
        genes_str, cells_str = args.genes_in_cells
        # Parse the comma-separated strings
        genes = qt.parse_comma_list(genes_str)
        cells = qt.parse_comma_list(cells_str)
        
        if not genes or not cells:
            print("❌ Please specify both genes and cells")
            return
        
        if args.verbose:
            print(f"📊 Querying {len(genes)} genes in {len(cells)} cells")
        
        results = qt.query_genes_in_cells(genes, cells, args.cell_line, args.min_umi)
    
    elif args.gene:
        genes = args.gene  # Already parsed by gene_symbols
        
        if args.verbose:
            print(f"📊 Querying {len(genes)} genes with min_umi={args.min_umi}")
        
        results = qt.query_genes(genes, args.cell_line, args.mutation, 
                                 args.min_umi, args.include_zeros, args.normalization)
    
    elif args.cell:
        cells = args.cell  # Already parsed by cell_barcodes
        
        if args.verbose:
            print(f"📊 Querying {len(cells)} cells")
        
        if len(cells) == 1:
            # Single cell query
            cell_info = qt.query_cell(cells[0], args.cell_line, args.top_n)
            if cell_info:
                # Convert to DataFrame for consistent output
                rows = []
                for cell, info in cell_info.items():
                    row = {
                        'cell_barcode': cell,
                        'cell_line': info['cell_line'],
                        'mutation': info['mutation'],
                        'total_umis': info['total_umis'],
                        'genes_detected': info['genes_detected']
                    }
                    for i, gene_info in enumerate(info['top_genes'][:args.top_n]):
                        row[f'top_gene_{i+1}'] = gene_info['gene']
                        row[f'top_gene_{i+1}_expr'] = gene_info['expression']
                    rows.append(row)
                results = pd.DataFrame(rows)
        else:
            results = qt.query_cells(cells, args.cell_line, args.top_n)
    
    # Handle results
    if results is not None:
        # Apply threshold if specified
        if args.threshold > 0 and isinstance(results, pd.DataFrame):
            if args.threshold_type == 'absolute':
                # Apply to gene expression columns
                for col in results.columns:
                    if col.endswith('_expression') or col in (args.gene or []):
                        threshold_col = f"{col}_thresholded"
                        results[threshold_col] = results[col] >= args.threshold
        
        # Save results
        if args.output:
            # Handle different formats
            compression = None if args.compression == 'none' else args.compression
            
            if args.format == 'csv':
                results.to_csv(args.output, index=False, compression=compression)
            elif args.format == 'tsv':
                results.to_csv(args.output, sep='\t', index=False, compression=compression)
            elif args.format == 'json':
                results.to_json(args.output, orient='records', indent=2)
            
            elapsed = time.time() - start_time
            print(f"\n💾 Results saved to: {args.output} ({elapsed:.1f} seconds)")
        
        # Print summary if verbose
        if args.verbose:
            print(f"\n📊 Query completed in {time.time() - start_time:.1f} seconds")
            if isinstance(results, pd.DataFrame):
                print(f"   Rows: {len(results)}, Columns: {len(results.columns)}")
    
    else:
        print("❌ No results returned")
        sys.exit(1)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\n⚠️  Interrupted by user")
        sys.exit(130)
    except Exception as e:
        print(f"\n❌ Unexpected error: {e}")
        if 'args' in locals() and hasattr(args, 'verbose') and args.verbose >= 2:
            import traceback
            traceback.print_exc()
        sys.exit(1)
