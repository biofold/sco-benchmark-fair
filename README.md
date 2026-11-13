# **ELIXIR - Single Cell Omics Benchmark - FAIR Metadata** 

**[→ View the FAIR Data Portal](https://biofold.github.io/sco-benchmark-fair)**

**Dataset:** [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24196083.svg)](https://doi.org/10.5281/zenodo.24196083) &nbsp; **Project report:** [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24150233.svg)](https://doi.org/10.5281/zenodo.24150233)


## **ELIXIR Single Cell Omics Benchmark Dataset 1 (GSE243665)** 

### **Executive Summary**

This document summarizes the FAIRification activities performed on the BE1 dataset, a single-cell RNA sequencing benchmark experiment designed to provide controlled cancer heterogeneity using lung cancer cell lines with seven distinct driver mutations. The dataset, deposited under GEO accession GSE243665, has been systematically transformed to comply with the Single Cell Schema (version 0.4) and enhanced with comprehensive metadata following FAIR (Findable, Accessible, Interoperable, Reusable) principles.

---

## **1\. Dataset Overview**

### **1.1 Original Dataset Description**

Title: A single cell RNAseq benchmark experiment embedding "controlled" cancer heterogeneity

Purpose: Provide a controlled heterogeneous environment for developing and validating computational methodologies for analysing cancer heterogeneity using scRNA-seq data.

Key Characteristics:

* Seven lung cancer cell lines with distinct driver mutations:  
  * PC9 (EGFR Del19, activating mutation)  
  * A549 (KRAS p.G12S, growth and proliferation)  
  * NCI-H596/HTB178 (MET Del14, enhanced apoptosis protection)  
  * NCI-H1395/CRL5868 (BRAF p.G469A, gain of function)  
  * DV90 (ERBB2 p.V842I, increased kinase activity)  
  * HCC78 (SLC34A2-ROS1 Fusion)  
  * CCL-185-IG (EML4-ALK Fusion in A549 isogenic background)  
* PBMC population (CMO 308) included in the GEO series; it is not described in the article  
* Total cells: 29,128 (seven cell lines) \+ 478 (PBMC) = 29,606  
* Platform: 10X Genomics (NovaSeq X Plus 10B flow-cell according to the article; the GEO platform record GPL24676 is Illumina NovaSeq 6000\)

### **1.2 Original Data Locations**

| Repository | Accession/DOI | Content |
| :---- | :---- | :---- |
| GEO | GSE243665 | Count matrices, configuration files |
| SRA | SRP462078 / PRJNA1019356 | Raw FASTQ files |
| figshare | 10.6084/m9.figshare.23939481.v1 | Count matrices (additional format) |
| figshare | 10.6084/m9.figshare.23284748.v1 | Supporting data (CCLE, driver gene targets) |
| figshare | 10.6084/m9.figshare.24744996.v1 | Supporting code and documentation |

---

## **2\. FAIRification Methodology**

### **2.1 Schema Selection**

The dataset was FAIRified using the Single Cell Schema version 0.4 (Darwin Core extension for single-cell RNA-seq), which provides:

* 15 standardized metadata sheets covering all aspects of single-cell experiments  
* Controlled vocabularies for technology, library preparation, and sequencing parameters  
* Cross-referencing capabilities between experimental components  
* Compatibility with existing bioinformatics standards

### **2.2 Information Sources Integrated**

Multiple sources were harmonized to create a comprehensive metadata record:

1. GEO Repository (GSE243665) \- Accession numbers, sample information, file listings  
2. Published Paper (Scientific Data, 2024\) \- Experimental design, validation results  
3. SRA Database \- Raw data access, BioProject linkage  
4. figshare Records \- Supplementary data, code, and documentation  
5. 10X Genomics Documentation \- Technical parameters, kit versions

### **2.3 FAIRification Process**

The following steps were performed:

Findable:

* Assigned persistent identifiers (DOI, GEO accession, SRA accession, BioProject ID)  
* Created comprehensive metadata with searchable keywords  
* Documented all data locations with direct URLs  
* Included bibliographic citation with PMID

Accessible:

* Provided direct FTP links to all supplementary files  
* Documented embargo status (now public as of Dec 2023\)  
* Included authentication/authorization information where applicable  
* Listed file sizes and formats for each data file

Interoperable:

* Used standardized schema with controlled vocabularies  
* Included all parameters for reproducing the analysis (Cell Ranger commands, configuration files)  
* Provided cross-references between related data objects  
* Used community-standard file formats (MTX, HDF5, FASTQ, CSV)

Reusable:

* Comprehensive experimental protocol documentation  
* Quality control metrics for each sample  
* Validation results against CCLE data  
* Usage notes including R Shiny App for generating controlled mixtures  
* Licensing: metadata package CC BY-NC 4.0; source article and figshare records CC BY 4.0

---

## **3\. Schema Implementation Details**

### **3.1 Completed Metadata Sheets**

The 15 sheets of the Single Cell Schema were populated:

| Sheet Name | Key Information Added |
| :---- | :---- |
| study | Title, description, citation, created date, technology, licence |
| person | Corresponding author of the data submission, with affiliation, funder and grant |
| sample | 8 populations (7 cell lines \+ PBMC) with taxon IDs, PMID references, source information |
| dissociation | Protocol ID, dissociation method, viability criteria, kit information |
| cell\_suspension | Per-population cell counts, suspension details |
| lib\_prep | 10X Genomics v3.1 kit details, amplification method, library construction parameters |
| sequencing | Platform, instrument model, read layout, barcode structure, Q30 metrics |
| analysis\_derived\_data | Cell type inference, QC filters, derived attributes |
| raw\_data\_processing | Cell Ranger v7.1.0 parameters, reference genome (GRCh38-2020-A) |
| downstream\_processing | Validation clustering parameters (Louvain, resolution 0.1) |
| data\_availability\_checklist | All data locations with DOIs and accessions |
| file | FASTQ files (SRA accession SRP462078); the 25 GEO files are listed in section 4.3 |
| expression\_data\_process\_setting | Complete processing parameters including intron mode |
| expression\_data\_file | Matrix file details (figshare) |
| geo\_samples | GSM accessions of the 8 GEO samples (4 gene-expression and 4 cell-multiplexing libraries) |

### **3.2 Controlled Vocabularies Used**

* Library strategy: RNA-Seq  
* Library source: TRANSCRIPTOMIC  
* Library selection: cDNA  
* Platform: ILLUMINA  
* Instrument model: Illumina NovaSeq 6000 (GEO platform record; the article reports a NovaSeq X Plus 10B flow-cell)  
* Suspension type: Cell  
* Spike-in: No  
* Primeness: 3'

### **3.3 Cross-References Established**

The schema maintains referential integrity through:

* Study ID (GSE243665) links all sheets  
* Sample IDs map to Cell Suspension IDs  
* Dissociation Protocol ID referenced by all suspensions  
* Library Preparation ID links to Sequencing ID  
* File IDs connected to Library Prep and Sequencing

---

## **4\. Enhanced Metadata Elements**

### **4.1 Technical Parameters**

| Parameter | Value |
| :---- | :---- |
| Cell Ranger version | 7.1.0 |
| Reference genome | GRCh38-2020-A |
| Gene annotation (QC plots) | Homo\_sapiens.GRCh38.99.gtf (rCASC) |
| Sequencing instrument | NovaSeq X Plus 10B flow-cell |
| Total reads | 2.46 billion |
| Q30 score | ≥71.26% |
| UMI barcode | 12 bp (Read 1, offset 16; 10x Genomics 3' v3.1 specification) |
| Cell barcode | 16 bp (Read 1, offset 0\) |
| cDNA read | 150 bp (Read 2\) |

### **4.2 Per-Population Statistics**

Cell Ranger statistics of the article (Table 1 of Arigoni et al., 2024); PBMC are only available in GEO.

| Cell Line | Cells | Median reads per cell | Median genes per cell | Total genes detected | Median UMI counts per cell |
| :---- | :---- | :---- | :---- | :---- | :---- |
| A549 | 6,898 | 21,728 | 2,315 | 27,584 | 5,961 |
| CCL-185-IG | 6,354 | 21,493 | 2,294 | 26,855 | 5,935 |
| NCI-H1395 (CRL5868) | 2,673 | 17,240 | 1,607 | 25,705 | 4,532 |
| DV90 | 2,998 | 21,241 | 2,160 | 25,786 | 5,628 |
| HCC78 | 2,748 | 28,423 | 2,665 | 25,425 | 7,808 |
| NCI-H596 (HTB178) | 2,965 | 19,652 | 2,113 | 25,281 | 5,300 |
| PC9 | 4,492 | 27,026 | 2,746 | 26,881 | 7,437 |
| PBMC | 478 | not in the article | | | |

Total: 29,128 cells in the seven cell lines; 29,606 including PBMC.

### **4.3 File Inventory**

GEO Supplementary Files (24 files):

* 8 samples × 3 files (barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz)  
* 1 configuration file (multi\_gex.csv.gz)

Total size: 267 MB compressed

figshare Resources:

* Count matrices (10X Genomics format)  
* CCLE supporting data with driver gene relationships  
* Code for extracting CCLE information

---

## **5\. Validation and Quality Control**

### **5.1 Technical Validation (from publication)**

1. BE1-500 dataset: 500 cells randomly selected from each cell line, clustered using Louvain modularity (resolution 0.1) → 6 clusters, each predominantly from specific cell line  
2. CCLE comparison: Pseudo-bulks compared with bulk transcriptomes from Cancer Cell Lines Encyclopedia → hierarchical clustering confirmed alignment between single-cell pseudo-bulks and respective cell line transcriptomes  
3. External dataset integration: Comparison with Tian et al. (2019), Aissa et al. (2021), and Clark et al. (2023) datasets confirmed expected biological relationships  
4. Quality filtering: rCASC mitoRiboUmi plot identified low-quality cells (NCI-H1395 reduced from 2673 to 1939 after QC)

### **5.2 Cell Ranger Quality Metrics**

No alerts generated for any cell line during processing. Complete metrics available in web\_summary.html files for each sample.

---

## **6\. FAIR Assessment Results**

### **6.1 Findability (self-assessment)**

* F1: Dataset assigned globally unique persistent identifiers  
  * GEO: GSE243665  
  * DOI: 10.1038/s41597-024-03002-y  
  * figshare DOIs: 10.6084/m9.figshare.23939481.v1, 10.6084/m9.figshare.23284748.v1  
  * SRA: SRP462078  
  * BioProject: PRJNA1019356  
* F2: Rich metadata describes context, purpose, and experimental design  
* F3: Metadata includes clear linkage to all data identifiers  
* F4: Dataset deposited in GEO, SRA and figshare; the article is indexed in PubMed

### **6.2 Accessibility (self-assessment)**

* A1: Standard protocol (HTTP, FTP) for data retrieval  
  * FTP links to all supplementary files  
  * SRA Run Selector for raw data  
  * figshare download options  
* A1.1: Open access (metadata package CC BY-NC 4.0; source article and figshare records CC BY 4.0)  
* A1.2: Authentication not required for public data  
* A2: Metadata remains accessible even if data becomes unavailable

### **6.3 Interoperability (self-assessment)**

* I1: Standard formal language (JSON-LD compatible structure)  
* I2: Community standards (Single Cell Schema v0.4, Darwin Core)  
* I3: Qualified references to other data (PMID links, cell line identifiers)

### **6.4 Reusability (self-assessment)**

* R1: Rich provenance documentation  
  * Complete experimental protocols  
  * Cell Ranger processing commands  
  * Quality control thresholds  
* R1.1: Clear usage licence (metadata package CC BY-NC 4.0; source article and figshare records CC BY 4.0)  
* R1.2: Detailed provenance including culture conditions, dissociation methods  
* R1.3: Community standards compliance ensures domain-relevant reuse

### **6.5 Machine-readable FAIR report**

The automated FAIR assessment of this repository is available as a machine-readable JSON file, [`fair-report/fair-report.json`](fair-report/fair-report.json), validated by the JSON Schema [`fair-report/fair-report.schema.json`](fair-report/fair-report.schema.json). It records the assessed repository, commit, release and DOI, the total score and the four FAIR principle scores with the criteria that are not met, the suggested improvements, and the limitations of the checklist. The file is produced by `scripts/test_fair_gitrepo.py` followed by `scripts/export_fair_report_json.py`; it can be regenerated at any release (see the usage in the header of the exporter script). Current values: total 81.25/100 (Findable 85, Accessible 80, Interoperable 75, Reusable 85); the file records the assessed commit and release.

### **6.6 Automation (GitHub Actions)**

Two workflows in `.github/workflows/` automate the checks:

* **Validate metadata** (`validate.yml`) runs on every push or pull request that changes the metadata. It validates `CITATION.cff` (with `cffconvert`), `.zenodo.json`, `fair-report/fair-report.json` against its JSON Schema and the Single Cell Schema record (`scrna-schema/`).
* **Regenerate FAIR report** (`regenerate-fair-report.yml`) is started on demand from the *Actions* tab (*Run workflow*, entering the release tag; the tag must match the version in `CITATION.cff`). It re-assesses the repository against the comparison repositories listed in `fair-report/comparison_repos.txt`, regenerates `fair-report.json`, the FAIR dashboards and the score numbers shown in the README and pages, validates the result and opens a pull request. Nothing is committed to the default branch: after review, merge the pull request and publish a new tagged release, which Zenodo archives automatically.

The benchmark data are a fixed community benchmark and are not refreshed on a schedule; they change only through deliberate curation followed by a new release.

---

## **7\. Usage Notes and Derived Resources**

### **7.1 R Shiny Application**

A dedicated R Shiny App ([http://aisc.hpc4ai.unito.it:3838/](http://aisc.hpc4ai.unito.it:3838/)) enables:

* Generation of sparse matrices by blending the seven cell lines at various ratios  
* Output in 10X Genomics format with cell barcodes containing cell line names  
* Creation of user-defined heterogeneity levels for benchmarking

### **7.2 Docker Container**

Cell Ranger v7.1.0 available as Docker container:

`text`

`docker.io/repbioinfo/cellranger.2023.7.1.0`

### **7.3 Configuration File**

The `multi_gex.csv` file (available in GEO) contains complete sample multiplexing configuration.

---

## **8\. Limitations and Considerations**

1. NCI-H1395 quality: Higher proportion of stressed cells (after QC reduces from 2673 to 1939 cells)  
2. Sequencing instrument: the article reports a NovaSeq X Plus 10B flow-cell, while the GEO platform record is Illumina NovaSeq 6000  
3. PBMC population: present in GEO (CMO 308, 478 cells) but not described in the article  
4. Per-population Cell Ranger statistics: use Table 1 of the article (section 4.2)

---

## **9\. Conclusions**

The BE1 dataset has been successfully FAIRified according to ELIXIR Single Cell Omics Community standards. The comprehensive metadata now enables:

1. Reproducibility: All experimental and computational parameters documented  
2. Reusability: Clear licensing and usage guidelines  
3. Interoperability: Standard schema enables integration with other datasets  
4. Findability: Multiple persistent identifiers and searchable metadata

The FAIRified dataset serves as an exemplar for controlled heterogeneity benchmark experiments and provides a valuable resource for developing and validating computational methods in cancer single-cell genomics.

---

## **10\. References**

1. Arigoni M, Ratto ML, Riccardo F, Balmas E et al. A single cell RNAseq benchmark experiment embedding "controlled" cancer heterogeneity. Sci Data 2024 Feb 2;11(1):159. doi:10.1038/s41597-024-03002-y. PMID: 38307867  
2. Single Cell Schema v0.4. [https://singlecellschemas.org](https://singlecellschemas.org/)  
3. Wilkinson, M., Dumontier, M., Aalbersberg, I. et al. The FAIR Guiding Principles for scientific data management and stewardship. Sci Data 3, 160018 (2016)  
4. 10X Genomics. Chromium Next GEM Single Cell 3' Kit v3.1 User Guide (CG00390 Rev C)

---

Prepared by: Emidio Capriotti, Pier Luigi Martelli, Castrense Savojardo (University of Bologna, ELIXIR Italy node)  
Date: 12 November 2026  
Dataset version: 1.0 (public since Dec 2023\)  
Schema version: dwc\_sc\_rnaseq v0.4

---

## **11\. How to cite**

Please cite both the original publication and this package.

1. Arigoni M, Ratto ML, Riccardo F, et al. A single cell RNAseq benchmark experiment embedding "controlled" cancer heterogeneity. Sci Data 2024;11:159. doi:[10.1038/s41597-024-03002-y](https://doi.org/10.1038/s41597-024-03002-y)
2. Capriotti E, Martelli PL, Savojardo C. *ELIXIR - Single Cell Omics Benchmark - FAIR Metadata*. Version 1.2.0. Zenodo. doi:[10.5281/zenodo.24196083](https://doi.org/10.5281/zenodo.24196083) (concept DOI, all versions).
3. Project report: Capriotti E, Martelli PL, Savojardo C. *FAIRification of the ELIXIR Single Cell Omics Benchmark*. Zenodo. doi:[10.5281/zenodo.24150233](https://doi.org/10.5281/zenodo.24150233) (concept DOI, all versions).

*This project is part of the ELIXIR Single Cell Omics Community and Tools Platform efforts to promote FAIR data practices.*
