#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
})

# ============================================================
# GSEA NES bar plots
# HCT116 cells infected with Gemella morbillorum
#
# This script creates separate Hallmark GSEA NES bar plots
# for Inf4h vs Ctr and Inf24h vs Ctr.
#
# Input:
#   results/gsea/Gm_HCT116_GSEA_HALLMARK_Inf4h_vs_Ctr.csv
#   results/gsea/Gm_HCT116_GSEA_HALLMARK_Inf24h_vs_Ctr.csv
#
# Output:
#   results/figures/Inf4h_Hallmark_GSEA_Top20_by_padj.pdf
#   results/figures/Inf4h_Hallmark_GSEA_Top20_by_padj.png
#   results/figures/Inf4h_Hallmark_GSEA_Top20_by_padj_data.csv
#
#   results/figures/Inf24h_Hallmark_GSEA_Top20_by_padj.pdf
#   results/figures/Inf24h_Hallmark_GSEA_Top20_by_padj.png
#   results/figures/Inf24h_Hallmark_GSEA_Top20_by_padj_data.csv
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

gsea_dir <- file.path(project_dir, "results", "gsea")
fig_dir  <- file.path(project_dir, "results", "figures")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Function to create one plot
# -----------------------------

make_gsea_barplot <- function(
  input_file,
  output_prefix,
  plot_title,
  top_n = 20
) {
  if (!file.exists(input_file)) {
    stop("File not found: ", input_file)
  }

  gsea_df <- readr::read_csv(input_file, show_col_types = FALSE)

  required_cols <- c("ID", "NES", "p.adjust")
  missing_cols <- setdiff(required_cols, colnames(gsea_df))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in ",
      input_file,
      ": ",
      paste(missing_cols, collapse = ", ")
    )
  }

  plot_df <- gsea_df %>%
    filter(!is.na(NES), !is.na(p.adjust)) %>%
    arrange(p.adjust) %>%
    slice_head(n = top_n) %>%
    mutate(
      pathway_label = ID,
      pathway_label = factor(
        pathway_label,
        levels = pathway_label[order(NES)]
      )
    )

  if (nrow(plot_df) == 0) {
    stop("No pathways available for plotting: ", input_file)
  }

  # Symmetric x-axis around zero
  max_abs_nes <- max(abs(plot_df$NES), na.rm = TRUE)
  x_lim <- c(-1.1 * max_abs_nes, 1.1 * max_abs_nes)

  p <- ggplot(plot_df, aes(x = pathway_label, y = NES, fill = p.adjust)) +
    geom_col(width = 0.8, color = "black", linewidth = 0.25) +
    coord_flip() +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    scale_y_continuous(limits = x_lim) +
    scale_fill_gradient(
      low = "red",
      high = "blue",
      name = "p.adjust"
    ) +
    labs(
      title = plot_title,
      x = NULL,
      y = "NES"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 12),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )

  pdf_file <- file.path(fig_dir, paste0(output_prefix, ".pdf"))
  png_file <- file.path(fig_dir, paste0(output_prefix, ".png"))
  csv_file <- file.path(fig_dir, paste0(output_prefix, "_data.csv"))

  ggsave(
    filename = pdf_file,
    plot = p,
    width = 8,
    height = 5.2
  )

  ggsave(
    filename = png_file,
    plot = p,
    width = 8,
    height = 5.2,
    dpi = 300
  )

  readr::write_csv(plot_df, csv_file)

  message("[INFO] Saved: ", pdf_file)
  message("[INFO] Saved: ", png_file)
  message("[INFO] Saved: ", csv_file)
}

# -----------------------------
# Input files
# -----------------------------

gsea_4h_file <- file.path(
  gsea_dir,
  "Gm_HCT116_GSEA_HALLMARK_Inf4h_vs_Ctr.csv"
)

gsea_24h_file <- file.path(
  gsea_dir,
  "Gm_HCT116_GSEA_HALLMARK_Inf24h_vs_Ctr.csv"
)

# -----------------------------
# Create plots
# -----------------------------

make_gsea_barplot(
  input_file = gsea_4h_file,
  output_prefix = "Inf4h_Hallmark_GSEA_Top20_by_padj",
  plot_title = "Inf4h_Hallmark_GSEA_Top20_by_padj",
  top_n = 20
)

make_gsea_barplot(
  input_file = gsea_24h_file,
  output_prefix = "Inf24h_Hallmark_GSEA_Top20_by_padj",
  plot_title = "Inf24h_Hallmark_GSEA_Top20_by_padj",
  top_n = 20
)

message("[INFO] GSEA NES bar plot generation completed.")
