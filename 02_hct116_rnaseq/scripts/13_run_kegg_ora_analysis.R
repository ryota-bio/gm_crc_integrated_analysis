#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

# ============================================================
# KEGG ORA only
# HCT116 cells infected with Gemella morbillorum
#
# This script performs KEGG pathway over-representation analysis
# for upregulated and downregulated DEGs in:
#   - Inf4h_vs_Ctr
#   - Inf24h_vs_Ctr
#
# This script only writes ORA result CSV files.
# It does not create dot plots.
# ============================================================

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

message("[INFO] Project directory: ", project_dir)

deg_dir <- file.path(project_dir, "results", "deseq2_deg")
ora_dir <- file.path(project_dir, "results", "ora")

dir.create(ora_dir, recursive = TRUE, showWarnings = FALSE)

deg_files <- list(
  Inf4h_vs_Ctr = file.path(deg_dir, "Gm_HCT116_DESeq2_Inf4h_vs_Ctr.csv"),
  Inf24h_vs_Ctr = file.path(deg_dir, "Gm_HCT116_DESeq2_Inf24h_vs_Ctr.csv")
)

for (f in deg_files) {
  if (!file.exists(f)) {
    stop("File not found: ", f)
  }
}

get_symbol_column <- function(df) {
  if ("gene_name" %in% colnames(df)) {
    return("gene_name")
  }
  if ("symbol" %in% colnames(df)) {
    return("symbol")
  }
  if ("gene_symbol" %in% colnames(df)) {
    return("gene_symbol")
  }
  stop("No gene symbol column found. Expected gene_name, symbol, or gene_symbol.")
}

run_kegg_ora <- function(deg_file, contrast_name,
                         padj_cutoff = 0.05,
                         log2fc_cutoff = 1) {
  message("[INFO] Reading: ", deg_file)

  deg <- readr::read_csv(deg_file, show_col_types = FALSE)

  required_cols <- c("log2FoldChange", "padj")
  missing_cols <- setdiff(required_cols, colnames(deg))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ", deg_file, ": ",
         paste(missing_cols, collapse = ", "))
  }

  symbol_col <- get_symbol_column(deg)

  deg <- deg %>%
    filter(!is.na(.data[[symbol_col]]), .data[[symbol_col]] != "") %>%
    filter(!is.na(log2FoldChange), !is.na(padj))

  gene_map <- bitr(
    unique(deg[[symbol_col]]),
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  ) %>%
    distinct(SYMBOL, ENTREZID, .keep_all = TRUE)

  universe_entrez <- unique(gene_map$ENTREZID)

  up_symbols <- deg %>%
    filter(padj < padj_cutoff, log2FoldChange >= log2fc_cutoff) %>%
    pull(.data[[symbol_col]]) %>%
    unique()

  down_symbols <- deg %>%
    filter(padj < padj_cutoff, log2FoldChange <= -log2fc_cutoff) %>%
    pull(.data[[symbol_col]]) %>%
    unique()

  gene_sets <- list(
    up = up_symbols,
    down = down_symbols
  )

  for (direction in names(gene_sets)) {
    symbols <- gene_sets[[direction]]

    entrez <- gene_map %>%
      filter(SYMBOL %in% symbols) %>%
      pull(ENTREZID) %>%
      unique()

    out_file <- file.path(
      ora_dir,
      paste0(contrast_name, "_KEGG_", direction, "_ORA.csv")
    )

    if (length(entrez) == 0) {
      warning("[WARN] No ENTREZ genes for: ", contrast_name, " ", direction)
      readr::write_csv(tibble(), out_file)
      next
    }

    ekegg <- enrichKEGG(
      gene = entrez,
      universe = universe_entrez,
      organism = "hsa",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.20
    )

    ekegg <- setReadable(
      ekegg,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID"
    )

    ekegg_df <- as.data.frame(ekegg)

    readr::write_csv(ekegg_df, out_file)

    message(
      "[INFO] Saved: ", out_file,
      " (", nrow(ekegg_df), " terms)"
    )
  }
}

for (contrast_name in names(deg_files)) {
  run_kegg_ora(
    deg_file = deg_files[[contrast_name]],
    contrast_name = contrast_name
  )
}

message("[INFO] KEGG ORA analysis completed.")
