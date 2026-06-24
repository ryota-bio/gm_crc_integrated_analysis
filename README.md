# RNA-seq analysis of HCT116 cells infected with *Gemella morbillorum*

This repository contains the final scripts used to reproduce the RNA-seq analyses reported in the manuscript.

## Experimental design

- Cell line: HCT116
- Conditions:
  - Control
  - *Gemella morbillorum* infection for 4 h
  - *Gemella morbillorum* infection for 24 h
- Replicates: n = 3 per group
- Sequencing: paired-end RNA-seq
- Reference: GENCODE release 45

## Workflow

1. Quality control and trimming using fastp
2. Salmon index construction using GENCODE v45 transcript sequences
3. Transcript quantification using Salmon
4. Gene-level count and TPM matrix generation using tximport
5. Differential expression analysis using DESeq2
6. Hallmark gene set enrichment analysis using clusterProfiler and msigdbr
7. Figure generation in R

## Repository structure

- scripts/: analysis scripts
- metadata/: sample metadata
- environment/: conda environment and software versions
- README.md: repository description

## Input files

Raw FASTQ files are not included in this repository.

Place paired-end FASTQ files in:

    data/raw_fastq/

Reference files are not included in this repository.

Place the following GENCODE release 45 files in:

    reference/
    gencode.v45.transcripts.fa.gz
    gencode.v45.annotation.gtf.gz

## Running the analysis

Create and activate the conda environment:

    conda env create -f environment/conda_environment.yml
    conda activate gm_rnaseq

Run the pipeline:

    bash scripts/01_fastp_qc.sh
    bash scripts/02_salmon_index.sh
    bash scripts/03_salmon_quant.sh
    Rscript scripts/04_tximport_make_count_tpm.R
    Rscript scripts/05_deseq2_deg_analysis.R
    Rscript scripts/06_clusterprofiler_gsea_analysis.R
    Rscript scripts/07_make_figures.R

## Output files

The main processed expression matrices are generated in:

    results/gene_matrices/Gm_HCT116_salmon_gene_counts.csv
    results/gene_matrices/Gm_HCT116_salmon_gene_tpm.csv

The count matrix represents Salmon/tximport-derived gene-level estimated counts generated with countsFromAbundance = "lengthScaledTPM".

DESeq2 results are generated in:

    results/deseq2_deg/Gm_HCT116_DESeq2_Inf4h_vs_Ctr.csv
    results/deseq2_deg/Gm_HCT116_DESeq2_Inf24h_vs_Ctr.csv

Hallmark GSEA results are generated in:

    results/gsea/Gm_HCT116_GSEA_HALLMARK_Inf4h_vs_Ctr.csv
    results/gsea/Gm_HCT116_GSEA_HALLMARK_Inf24h_vs_Ctr.csv

## Software versions

Key software versions are listed in:

    environment/software_versions.txt

## Data availability

Raw sequencing data and processed gene expression matrices will be deposited in GEO/DDBJ under accession numbers to be assigned.

## Code availability

The scripts used for RNA-seq quality control, transcript quantification, gene-level matrix generation, differential expression analysis, gene set enrichment analysis, and figure generation are provided in this repository.
