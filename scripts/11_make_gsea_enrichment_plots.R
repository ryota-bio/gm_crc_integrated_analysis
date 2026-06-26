#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(clusterProfiler)
  library(msigdbr)
  library(enrichplot)
})

# ============================================================
# Hallmark GSEA enrichment plots
# HCT116 cells infected with Gemella morbillorum
#
# This script generates enrichment plots for selected Hallmark
# pathways using ranked gene lists from DESeq2 results.
#
# Input:
#   results/deseq2_deg/Gm_HCT116_DESeq2_Inf4h_vs_Ctr.csv
#   results/deseq2_deg/Gm_HCT116_DESeq2_Inf24h_vs_Ctr.csv
#
# Output:
#   results/figures/gsea_enrichment_plots/*.pdf
#   results/figures/gsea_enrichment_plots/*.png
# ============================================================

# -----------------------------
# Project paths
# -----------------------------

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])

if (length(script_path) == 0) {
  project_dir <- getwd()
} else {
  project_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
}

if (!dir.exists(file.path(project_dir, "results"))) {
  project_dir <- getwd()
}

deg_dir <- file.path(project_dir, "results", "deseq2_deg")
fig_dir <- file.path(project_dir, "results", "figures")
gsea_plot_dir <- file.path(fig_dir, "gsea_enrichment_plots")

dir.create(gsea_plot_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Parameters
# -----------------------------

min_gs_size <- 10
max_gs_size <- 500
pvalue_cutoff <- 1
eps_value <- 0

selected_pathways <- list(
  Inf4h_vs_Ctr = c(
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_HYPOXIA"
  ),
  Inf24h_vs_Ctr = c(
    "HALLMARK_HYPOXIA",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_GLYCOLYSIS",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
  )
)

deg_files <- list(
  Inf4h_vs_Ctr = file.path(deg_dir, "Gm_HCT116_DESeq2_Inf4h_vs_Ctr.csv"),
  Inf24h_vs_Ctr = file.path(deg_dir, "Gm_HCT116_DESeq2_Inf24h_vs_Ctr.csv")
)

for (f in deg_files) {
  if (!file.exists(f)) {
    stop("File not found: ", f)
  }
}

# -----------------------------
# Hallmark gene sets
# -----------------------------

get_hallmark_msigdb <- function() {
  hallmark_df <- tryCatch(
    {
      msigdbr(species = "Homo sapiens", collection = "H")
    },
    error = function(e) {
      msigdbr(species = "Homo sapiens", category = "H")
    }
  )

  hallmark_df %>%
    dplyr::select(gs_name, gene_symbol) %>%
    dplyr::distinct()
}

hallmark_t2g <- get_hallmark_msigdb()

# -----------------------------
# Helper functions
# -----------------------------

make_gene_rank <- function(deg_df) {
  required_cols <- c("gene_name", "log2FoldChange")
  missing_cols <- setdiff(required_cols, colnames(deg_df))

  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if ("stat" %in% colnames(deg_df)) {
    rank_col <- "stat"
  } else {
    rank_col <- "log2FoldChange"
  }

  rank_df <- deg_df %>%
    filter(
      !is.na(gene_name),
      !is.na(.data[[rank_col]])
    ) %>%
    group_by(gene_name) %>%
    slice_max(order_by = abs(.data[[rank_col]]), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(.data[[rank_col]]))

  gene_list <- rank_df[[rank_col]]
  names(gene_list) <- rank_df$gene_name
  gene_list <- sort(gene_list, decreasing = TRUE)

  gene_list
}

clean_pathway_name <- function(x) {
  x %>%
    str_replace("^HALLMARK_", "") %>%
    str_replace_all("_", " ") %>%
    str_to_title()
}

safe_filename <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9_]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace_all("_$", "")
}

run_gsea <- function(gene_list) {
  clusterProfiler::GSEA(
    geneList = gene_list,
    TERM2GENE = hallmark_t2g,
    minGSSize = min_gs_size,
    maxGSSize = max_gs_size,
    pvalueCutoff = pvalue_cutoff,
    pAdjustMethod = "BH",
    eps = eps_value,
    verbose = FALSE
  )
}

save_plot_pdf_png <- function(p, pdf_file, png_file, width = 7.2, height = 5.4) {
  ggsave(
    filename = pdf_file,
    plot = p,
    width = width,
    height = height
  )

  png(
    filename = png_file,
    width = width,
    height = height,
    units = "in",
    res = 300
  )
  print(p)
  dev.off()
}

plot_one_enrichment <- function(gsea_res, comparison_name, pathway_id) {
  gsea_df <- as.data.frame(gsea_res)

  if (!(pathway_id %in% gsea_df$ID)) {
    message("[INFO] Pathway not found in GSEA result: ", pathway_id)
    return(invisible(NULL))
  }

  pathway_label <- clean_pathway_name(pathway_id)

  pathway_info <- gsea_df %>%
    filter(ID == pathway_id) %>%
    slice(1)

  nes_value <- round(pathway_info$NES, 2)
  padj_value <- signif(pathway_info$p.adjust, 3)

  plot_title <- paste0(
    comparison_name,
    " - ",
    pathway_label,
    "\nNES = ",
    nes_value,
    ", adjusted P = ",
    padj_value
  )

  p <- enrichplot::gseaplot2(
    gsea_res,
    geneSetID = pathway_id,
    title = plot_title,
    base_size = 12,
    rel_heights = c(1.5, 0.4, 0.8)
  )

  output_prefix <- paste0(
    comparison_name,
    "_GSEA_enrichment_",
    safe_filename(pathway_id)
  )

  pdf_file <- file.path(gsea_plot_dir, paste0(output_prefix, ".pdf"))
  png_file <- file.path(gsea_plot_dir, paste0(output_prefix, ".png"))

  save_plot_pdf_png(
    p = p,
    pdf_file = pdf_file,
    png_file = png_file,
    width = 7.2,
    height = 5.4
  )

  message("[INFO] Saved: ", pdf_file)
  message("[INFO] Saved: ", png_file)
}

# -----------------------------
# Main analysis
# -----------------------------

for (comparison_name in names(deg_files)) {
  message("[INFO] Processing: ", comparison_name)

  deg_df <- readr::read_csv(deg_files[[comparison_name]], show_col_types = FALSE)

  gene_list <- make_gene_rank(deg_df)

  message("[INFO] Ranked genes: ", length(gene_list))

  gsea_res <- run_gsea(gene_list)

  gsea_csv <- file.path(
    gsea_plot_dir,
    paste0(comparison_name, "_Hallmark_GSEA_recomputed_for_enrichment_plots.csv")
  )

  readr::write_csv(as.data.frame(gsea_res), gsea_csv)
  message("[INFO] Saved: ", gsea_csv)

  for (pathway_id in selected_pathways[[comparison_name]]) {
    plot_one_enrichment(
      gsea_res = gsea_res,
      comparison_name = comparison_name,
      pathway_id = pathway_id
    )
  }
}

message("[INFO] Hallmark GSEA enrichment plot generation completed.")
