# Gemella morbillorum and colorectal cancer

This repository contains the computational analyses associated with our study investigating the potential role of *Gemella morbillorum* in colorectal cancer progression.

## Overview

The repository integrates three computational analyses:

1. Public fecal metagenomic analysis and machine-learning-based feature selection
2. RNA-seq analysis of HCT116 colorectal cancer cells following *Gemella morbillorum* infection
3. TCGA COAD/READ survival analysis of candidate host genes identified from the RNA-seq analysis

## Study workflow

Public metagenomic analysis  
→ Identification of colorectal cancer progression-associated bacterial species  
→ Selection of *Gemella morbillorum*  
→ Infection of HCT116 colorectal cancer cells  
→ RNA-seq and pathway analysis  
→ Identification of candidate migration-associated host genes  
→ TCGA survival analysis

## Repository structure

### `01_public_metagenome/`

Public fecal metagenomic analysis and machine-learning-based identification of bacterial species associated with colorectal cancer progression.

Main analyses include:

- Random Forest
- Extra Trees
- Mann–Whitney U test
- *Gemella morbillorum* detection-rate analysis
- *Gemella morbillorum* relative-abundance analysis across CRC stages

See `01_public_metagenome/README.md` for details.

### `02_hct116_rnaseq/`

RNA-seq analysis of HCT116 colorectal cancer cells infected with *Gemella morbillorum*.

Main analyses include:

- fastp quality control
- Salmon transcript quantification
- tximport gene-level summarization
- DESeq2 differential expression analysis
- Hallmark GSEA
- GO over-representation analysis
- KEGG over-representation analysis
- Visualization of differentially expressed genes and enriched pathways

See `02_hct116_rnaseq/README.md` for details.

### `03_tcga_survival/`

TCGA COAD/READ survival analysis of candidate host genes identified from the HCT116 RNA-seq analysis.

Candidate genes analyzed include:

- SERPINE1
- PLAUR
- SDC4
- CCN1

Main analyses include:

- Kaplan–Meier overall survival analysis
- Univariate Cox proportional hazards regression
- Multivariable Cox regression
- Proportional hazards assumption testing

See `03_tcga_survival/README.md` for details.

## Data availability

Large raw datasets are not redistributed through this repository.

Public metagenomic datasets should be obtained from their original repositories or publications as described in `01_public_metagenome/README.md`.

Raw RNA-seq data generated in this study are deposited in the DNA Data Bank of Japan (DDBJ). Accession information will be added upon public release.

TCGA clinical and gene-expression data should be obtained from the corresponding public data repositories as described in `03_tcga_survival/README.md`.

## Reproducibility

Each analysis directory contains the scripts and instructions required to reproduce the corresponding computational analysis.

Software environments and version information for the RNA-seq analysis are provided in:

`02_hct116_rnaseq/environment/`

## Citation

If you use this repository, please cite the associated publication.

Citation information will be added after publication.

## License

License information for this repository is provided in the repository license file.
