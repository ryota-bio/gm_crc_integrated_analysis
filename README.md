# Gemella morbillorum and colorectal cancer

This repository contains the computational analyses associated with our study investigating the potential role of *Gemella morbillorum* in colorectal cancer progression.

## Associated manuscript

This repository contains the analysis code associated with the manuscript:

**Gemella morbillorum promotes migration and induces migration- and cell adhesion-associated transcriptional responses in HCT 116 cells**

Authors: Ryota Mori, Toshifumi Hara, Ryo Kutsuna, Junko Tomida, Tomoharu Takeuchi, and Yoshiaki Kawamura

Citation information will be updated upon publication.

## Overview

This repository integrates three computational analyses performed in the study:

1. Public fecal metagenomic analysis and machine-learning-based feature selection
2. RNA-seq analysis of HCT116 colorectal cancer cells following *Gemella morbillorum* infection
3. TCGA COAD/READ survival analysis of candidate host genes identified from the RNA-seq analysis

The analyses are organized into separate subdirectories to facilitate reproducibility.

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
- Identification of bacterial species commonly ranked among the top features
- *Gemella morbillorum* detection-rate analysis
- *Gemella morbillorum* relative-abundance analysis across colorectal cancer stages

Main contents:

- `scripts/` — analysis scripts
- `results/` — tabulated analysis results
- `figures/` — figures generated from the analyses
- `requirements.txt` — Python package requirements
- `run_all.sh` — script for running the analysis workflow

See `01_public_metagenome/README.md` for detailed instructions and data-source information.

---

### `02_hct116_rnaseq/`

RNA-seq analysis of HCT116 colorectal cancer cells infected with *Gemella morbillorum*.

The RNA-seq experiment consisted of three experimental groups:

- Uninfected control (Ctr)
- Cells harvested immediately after 4 h of *G. morbillorum* infection (Inf4h)
- Cells cultured for an additional 24 h after the 4 h infection period following bacterial removal (Inf24h)

Three biological replicates were analyzed for each group.

Main analyses include:

- fastp quality control
- Salmon transcript quantification
- tximport gene-level summarization
- DESeq2 differential expression analysis
- Volcano plot generation
- Hallmark gene set enrichment analysis (GSEA)
- GSEA enrichment plots
- Leading-edge gene heatmaps
- Gene Ontology over-representation analysis
- KEGG pathway over-representation analysis
- Visualization of enriched pathways

Main contents:

- `scripts/` — RNA-seq analysis scripts
- `metadata/` — sample metadata
- `environment/` — software environment and version information

See `02_hct116_rnaseq/README.md` for detailed execution instructions.

---

### `03_tcga_survival/`

TCGA COAD/READ survival analysis of candidate host genes identified from the HCT116 RNA-seq analysis.

Candidate genes analyzed include:

- `SERPINE1`
- `PLAUR`
- `SDC4`
- `CCN1`

These genes were selected based on the RNA-seq and functional enrichment analyses as candidate genes associated with cell migration and adhesion following *G. morbillorum* infection.

Main analyses include:

- Kaplan–Meier overall survival analysis
- Univariate Cox proportional hazards regression
- Multivariable Cox regression
- Proportional hazards assumption testing

Main contents:

- `scripts/` — survival analysis scripts
- `results/` — statistical results
- `figures/` — Kaplan–Meier plots and Cox regression forest plots

See `03_tcga_survival/README.md` for detailed instructions.

## Data availability

Large raw datasets are not redistributed through this repository.

### Public metagenomic data

The public fecal metagenomic datasets used in this study were obtained from previously published studies. These datasets should be downloaded from their original repositories or publications as described in `01_public_metagenome/README.md`.

### RNA-seq data

Raw RNA-seq data generated in this study have been deposited in the DNA Data Bank of Japan (DDBJ) under the following accession numbers:

**DRR1069701–DRR1069709**

Gene-level RNA-seq data have been deposited in the Genomic Expression Archive (GEA) under accession number:

**E-GEAD-1294**

### TCGA data

TCGA COAD/READ gene-expression and clinical data should be obtained from the corresponding public data repositories as described in `03_tcga_survival/README.md`.

## Reproducibility

Each analysis directory contains the scripts and instructions required to reproduce the corresponding computational analysis.

The RNA-seq software environment and version information are provided in:

`02_hct116_rnaseq/environment/`

The public metagenomic analysis includes a Python requirements file:

`01_public_metagenome/requirements.txt`

Because several analyses rely on publicly available datasets that are not redistributed in this repository, users should obtain the required input datasets from the original repositories described in the corresponding README files.

## Code availability

All custom scripts used for the public metagenomic and machine-learning analyses, HCT116 RNA-seq analysis, and TCGA COAD/READ survival analysis are provided in this repository.

Repository URL:

https://github.com/ryota-bio/gm_crc_integrated_analysis

A versioned release associated with the manuscript will be created upon finalization of the analysis code.

## Citation

If you use this repository, please cite the associated publication.

Citation information will be added upon publication.

## License

This repository is distributed under the MIT License.

See `LICENSE` for details.
