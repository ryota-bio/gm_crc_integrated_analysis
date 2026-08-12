# RNA-seq analysis of HCT116 cells infected with *Gemella morbillorum*

This repository contains scripts used for RNA-seq analysis of HCT116 colorectal cancer cells infected with *Gemella morbillorum*.

The analysis includes read quality control, transcript quantification, gene-level count and TPM matrix generation, differential expression analysis, Hallmark GSEA, GO Biological Process over-representation analysis, KEGG pathway over-representation analysis, and figure generation.

## Experimental design

HCT116 cells were analyzed under the following conditions:

- Control
- *Gemella morbillorum* infection for 4 h
- *Gemella morbillorum* infection for 24 h

Each condition included three biological replicates.

RNA-seq libraries were generated as paired-end reads. Transcript quantification was performed using Salmon with a GENCODE release 45 transcriptome reference.

## Workflow

1. Quality control and read trimming using fastp
2. Salmon index construction
3. Transcript quantification using Salmon
4. Gene-level count and TPM matrix generation using tximport
5. Differential expression analysis using DESeq2
6. Hallmark GSEA using clusterProfiler
7. PCA figure generation
8. Hallmark GSEA NES bar plot generation
9. GO Biological Process ORA dot plot generation
10. KEGG pathway ORA dot plot generation
11. Hallmark GSEA enrichment plot generation
12. Hallmark GSEA leading edge heatmap generation

## Repository structure

```text
gm_hct116_rnaseq/
├── README.md
├── LICENSE
├── environment/
│   ├── conda_environment.yml
│   └── software_versions.txt
├── metadata/
│   └── sample_metadata.csv
└── scripts/
    ├── 01_fastp_qc.sh
    ├── 02_salmon_index.sh
    ├── 03_salmon_quant.sh
    ├── 04_tximport_make_count_tpm.R
    ├── 05_deseq2_deg_analysis.R
    ├── 08_clusterprofiler_gsea_analysis.R
    ├── 06_make_basic_figures.R
    ├── 09_make_gsea_nes_barplot.R
    ├── 14_make_go_ora_dotplots.R
    └── 15_make_kegg_ora_dotplots.R
```

## Input files

Raw FASTQ files are not included in this repository.

The analysis assumes the following local directory structure:

```text
data/raw_fastq/
reference/
results/
```

The sample information is provided in:

```text
metadata/sample_metadata.csv
```

## Running the analysis

Run the scripts from the project root directory.

```bash
bash scripts/01_fastp_qc.sh
bash scripts/02_salmon_index.sh
bash scripts/03_salmon_quant.sh

Rscript scripts/04_tximport_make_count_tpm.R
Rscript scripts/05_deseq2_deg_analysis.R
Rscript scripts/06_make_basic_figures.R
Rscript scripts/07_make_volcano_plots.R

Rscript scripts/08_clusterprofiler_gsea_analysis.R
Rscript scripts/09_make_gsea_nes_barplot.R
Rscript scripts/10_make_gsea_enrichment_plots.R
Rscript scripts/11_make_gsea_leading_edge_heatmap.R

Rscript scripts/12_run_go_ora_analysis.R
Rscript scripts/13_run_kegg_ora_analysis.R
Rscript scripts/14_make_go_ora_dotplots.R
Rscript scripts/15_make_kegg_ora_dotplots.R
```

## Output files

Gene-level expression matrices:

- `results/gene_matrices/Gm_HCT116_salmon_gene_counts.csv`
- `results/gene_matrices/Gm_HCT116_salmon_gene_tpm.csv`

Differential expression results:

- `results/deseq2_deg/Gm_HCT116_DESeq2_Inf4h_vs_Ctr.csv`
- `results/deseq2_deg/Gm_HCT116_DESeq2_Inf24h_vs_Ctr.csv`

Hallmark GSEA results:

- `results/gsea/Gm_HCT116_GSEA_HALLMARK_Inf4h_vs_Ctr.csv`
- `results/gsea/Gm_HCT116_GSEA_HALLMARK_Inf24h_vs_Ctr.csv`

GO-BP and KEGG ORA results:

- `results/ora/*GO_BP*_ORA.csv`
- `results/ora/*KEGG*_ORA.csv`

Figures:

- `results/figures/Figure_PCA_RNAseq.pdf`
- `results/figures/Inf4h_Hallmark_GSEA_Top20_by_padj.pdf`
- `results/figures/Inf4h_Hallmark_GSEA_Top20_by_padj.png`
- `results/figures/Inf24h_Hallmark_GSEA_Top20_by_padj.pdf`
- `results/figures/Inf24h_Hallmark_GSEA_Top20_by_padj.png`
- `results/figures/volcano_plots/*.pdf/png`
- `results/figures/volcano_plots/*_highlighted_genes.csv`
- `results/figures/*GO_BP*_byCount_dotplot.pdf`
- `results/figures/*GO_BP*_byCount_dotplot.png`
- `results/figures/*GO_BP*_byCount_dotplot_data.csv`
- `results/figures/*KEGG*_byCount_dotplot.pdf`
- `results/figures/*KEGG*_byCount_dotplot.png`
- `results/figures/*KEGG*_byCount_dotplot_data.csv`
- `results/figures/gsea_enrichment_plots/*.pdf/png`
- `results/figures/gsea_leading_edge_heatmaps/*.pdf/png`
- `results/figures/gsea_leading_edge_heatmaps/*_genes.csv`

## Notes on ORA

GO-BP and KEGG over-representation analyses were performed using differentially expressed genes defined by:

```text
adjusted P value < 0.05
log2 fold change > 1 for upregulated genes
log2 fold change < -1 for downregulated genes
```

Dot plots show the top enriched terms or pathways ranked by Count. The x-axis represents GeneRatio, dot size represents gene count, and dot fill color represents adjusted P value.

## Notes on GSEA enrichment and leading edge heatmaps

Hallmark GSEA enrichment plots are generated for selected pathways to visualize the running enrichment score and leading-edge position.

Leading edge heatmaps show the top 30 leading-edge genes for selected Hallmark pathways. Control, 4 h infection, and 24 h infection samples are shown together in the same heatmap. Gene expression values are visualized as row-wise Z-scores using variance-stabilized expression values.

## Software versions

Software versions and package information are provided in:

```text
environment/software_versions.txt
```

The conda environment file is provided in:

```text
environment/conda_environment.yml
```

## Data availability

Raw and processed RNA-seq data will be deposited in a public nucleotide sequence database under an accession number to be provided upon acceptance.

## Code availability

All scripts used for RNA-seq processing, differential expression analysis, enrichment analysis, and figure generation are provided in this repository.

