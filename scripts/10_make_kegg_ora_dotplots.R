#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

# ============================================================
# KEGG ORA dot plots
# HCT116 cells infected with Gemella morbillorum
#
# This script performs KEGG pathway ORA for
# upregulated and downregulated DEGs in:
#   - Inf4h vs Ctr
#   - Inf24h vs Ctr
#
# Input:
#   results/deseq2_deg/Gm_HCT116_DESeq2_Inf4h_vs_Ctr.csv
#   results/deseq2_deg/Gm_HCT116_DESeq2_Inf24h_vs_Ctr.csv
#
# Output:
#   results/ora/*KEGG*_ORA.csv
#   results/figures/*KEGG*_dotplot.pdf
#   results/figures/*KEGG*_dotplot.png
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
ora_dir <- file.path(project_dir, "results", "ora")
fig_dir <- file.path(project_dir, "results", "figures")

dir.create(ora_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Parameters
# -----------------------------

padj_cutoff <- 0.05
log2fc_up_cutoff <- 1
log2fc_down_cutoff <- -1
top_n <- 20

# Use the same x-axis range for all KEGG dot plots
x_axis_max <- 0.20

# -----------------------------
# Input DEG files
# -----------------------------

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
# Helper functions
# -----------------------------

map_symbols_to_entrez <- function(symbols) {
  symbols <- symbols[!is.na(symbols)]
  symbols <- unique(symbols)

  if (length(symbols) == 0) {
    return(character(0))
  }

  mapped <- suppressMessages(
    clusterProfiler::bitr(
      symbols,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
  )

  unique(mapped$ENTREZID)
}

parse_gene_ratio <- function(x) {
  numerator <- as.numeric(sub("/.*", "", x))
  denominator <- as.numeric(sub(".*/", "", x))
  numerator / denominator
}

make_kegg_dotplot <- function(
  kegg_df,
  output_prefix,
  plot_title,
  top_n = 20,
  x_axis_max = 0.20
) {
  if (is.null(kegg_df) || nrow(kegg_df) == 0) {
    message("[INFO] No KEGG pathways to plot: ", output_prefix)
    return(invisible(NULL))
  }

  required_cols <- c("ID", "Description", "GeneRatio", "Count", "p.adjust")
  missing_cols <- setdiff(required_cols, colnames(kegg_df))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns for plotting: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  plot_df <- kegg_df %>%
    filter(!is.na(p.adjust)) %>%
    arrange(p.adjust) %>%
    slice_head(n = top_n) %>%
    mutate(
      GeneRatioNumeric = parse_gene_ratio(GeneRatio),
      Description_wrapped = stringr::str_wrap(Description, width = 36),
      Description_wrapped = factor(
        Description_wrapped,
        levels = rev(Description_wrapped)
      )
    )

  if (nrow(plot_df) == 0) {
    message("[INFO] No KEGG pathways after filtering: ", output_prefix)
    return(invisible(NULL))
  }

  # Increase plot height according to the number of displayed pathways
  fig_height <- max(6.5, 0.36 * nrow(plot_df) + 1.8)

  p <- ggplot(
    plot_df,
    aes(
      x = GeneRatioNumeric,
      y = Description_wrapped,
      size = Count,
      fill = p.adjust
    )
  ) +
    geom_point(
      shape = 21,
      color = "black",
      stroke = 0.35,
      alpha = 0.9
    ) +
    scale_fill_gradient(
      low = "#F8766D",
      high = "#619CFF",
      name = "p.adjust"
    ) +
    scale_size_continuous(name = "Count") +
    scale_x_continuous(
      breaks = seq(0, x_axis_max, by = 0.05)
    ) +
    coord_cartesian(xlim = c(0, x_axis_max)) +
    labs(
      title = plot_title,
      x = "GeneRatio",
      y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(size = 8.5, lineheight = 0.85),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 11),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      plot.margin = margin(t = 10, r = 20, b = 10, l = 10)
    )

  pdf_file <- file.path(fig_dir, paste0(output_prefix, ".pdf"))
  png_file <- file.path(fig_dir, paste0(output_prefix, ".png"))
  plot_data_file <- file.path(fig_dir, paste0(output_prefix, "_data.csv"))

  ggsave(
    filename = pdf_file,
    plot = p,
    width = 8,
    height = fig_height
  )

  ggsave(
    filename = png_file,
    plot = p,
    width = 8,
    height = fig_height,
    dpi = 300
  )

  readr::write_csv(plot_df, plot_data_file)

  message("[INFO] Saved: ", pdf_file)
  message("[INFO] Saved: ", png_file)
  message("[INFO] Saved: ", plot_data_file)
}

run_kegg_ora_for_one_set <- function(
  entrez_genes,
  universe_entrez,
  comparison_name,
  direction_name
) {
  if (length(entrez_genes) < 5) {
    message(
      "[INFO] Skipped KEGG ORA because fewer than 5 mapped genes: ",
      comparison_name,
      " ",
      direction_name
    )
    return(invisible(NULL))
  }

  kegg_res <- enrichKEGG(
    gene = entrez_genes,
    universe = universe_entrez,
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20
  )

  kegg_df <- as.data.frame(kegg_res)

  # Convert ENTREZ IDs to gene symbols when possible
  if (nrow(kegg_df) > 0) {
    kegg_df <- tryCatch(
      {
        readable_kegg <- setReadable(
          kegg_res,
          OrgDb = org.Hs.eg.db,
          keyType = "ENTREZID"
        )
        as.data.frame(readable_kegg)
      },
      error = function(e) {
        message("[INFO] KEGG setReadable failed; using ENTREZ IDs.")
        kegg_df
      }
    )
  }

  kegg_csv <- file.path(
    ora_dir,
    paste0(comparison_name, "_KEGG_", direction_name, "_ORA.csv")
  )

  readr::write_csv(kegg_df, kegg_csv)
  message("[INFO] Saved: ", kegg_csv)

  make_kegg_dotplot(
    kegg_df = kegg_df,
    output_prefix = paste0(comparison_name, "_KEGG_", direction_name, "_byPadj_dotplot"),
    plot_title = paste0(comparison_name, " KEGG ORA (", direction_name, ", by p.adjust)"),
    top_n = top_n,
    x_axis_max = x_axis_max
  )
}

# -----------------------------
# Main analysis
# -----------------------------

for (comparison_name in names(deg_files)) {
  message("[INFO] Processing: ", comparison_name)

  deg_df <- readr::read_csv(deg_files[[comparison_name]], show_col_types = FALSE)

  required_cols <- c("gene_name", "log2FoldChange", "padj")
  missing_cols <- setdiff(required_cols, colnames(deg_df))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in DEG file: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  deg_df <- deg_df %>%
    filter(!is.na(gene_name), !is.na(log2FoldChange), !is.na(padj))

  universe_symbols <- deg_df %>%
    pull(gene_name) %>%
    unique()

  up_symbols <- deg_df %>%
    filter(padj < padj_cutoff, log2FoldChange > log2fc_up_cutoff) %>%
    pull(gene_name) %>%
    unique()

  down_symbols <- deg_df %>%
    filter(padj < padj_cutoff, log2FoldChange < log2fc_down_cutoff) %>%
    pull(gene_name) %>%
    unique()

  universe_entrez <- map_symbols_to_entrez(universe_symbols)
  up_entrez <- map_symbols_to_entrez(up_symbols)
  down_entrez <- map_symbols_to_entrez(down_symbols)

  message("[INFO] Universe mapped genes: ", length(universe_entrez))
  message("[INFO] Up genes mapped: ", length(up_entrez))
  message("[INFO] Down genes mapped: ", length(down_entrez))

  run_kegg_ora_for_one_set(
    entrez_genes = up_entrez,
    universe_entrez = universe_entrez,
    comparison_name = comparison_name,
    direction_name = "up"
  )

  run_kegg_ora_for_one_set(
    entrez_genes = down_entrez,
    universe_entrez = universe_entrez,
    comparison_name = comparison_name,
    direction_name = "down"
  )
}

message("[INFO] KEGG ORA dot plot generation completed.")
