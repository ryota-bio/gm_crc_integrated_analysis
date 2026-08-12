# TCGA-COAD/READ candidate gene survival analysis

This repository contains scripts and summary outputs for overall survival analysis of selected candidate genes in TCGA-COAD/READ colorectal cancer samples.

## Candidate genes

- SERPINE1
- PLAUR
- CCN1
- SDC4

## Analysis overview

For each gene, TCGA-COAD/READ patients were divided into high- and low-expression groups using the median expression level.

Overall survival was analyzed using Kaplan-Meier curves, log-rank tests, Cox proportional hazards models, stage-adjusted Cox models, early/advanced stage-adjusted Cox models, early/advanced stage plus sex-adjusted Cox models, and cox.zph proportional hazards assumption tests.

## Input data

Input data are not redistributed in this repository.

Required input files:

- data/expression/tcga_coad_read_vst_protein_coding.rds
- data/clinical/cbioportal_coadread_tcga_pub_clinical_data.tsv

## Run

Rscript scripts/01_candidate_gene_survival_analysis_tcga_coad_read.R

## Outputs

Survival result tables are saved in results/survival/.

Figures are saved in figures/.

## Summary

High expression of SERPINE1 and CCN1 showed nominal associations with poorer overall survival in TCGA-COAD/READ in univariate analyses. PLAUR and SDC4 were not significantly associated with overall survival.

These analyses are exploratory association analyses and do not demonstrate causality.
