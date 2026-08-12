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
# This script creates Hallmark GSEA NES bar plots for
# Inf4h_vs_Ctr and Inf24h_vs_Ctr using common NES axis limits
# and a common p.adjust color scale.
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

gsea_dir <- file.path(project_dir, "results", "gsea")
fig_dir  <- file.path(project_dir, "results", "figures")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

gsea_4h_file <- file.path(
  gsea_dir,
  "Gm_HCT116_GSEA_HALLMARK_Inf4h_vs_Ctr.csv"
)

gsea_24h_file <- file.path(
  gsea_dir,
  "Gm_HCT116_GSEA_HALLMARK_Inf24h_vs_Ctr.csv"
)

for (f in c(gsea_4h_file, gsea_24h_file)) {
  if (!file.exists(f)) {
    stop("File not found: ", f)
  }
}

clean_pathway_label <- function(x) {
  x %>%
    str_replace("^HALLMARK_", "") %>%
    str_replace_all("_", " ") %>%
    str_to_title()
}

prepare_plot_df <- function(input_file, top_n = 20) {
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

  gsea_df %>%
    filter(!is.na(NES), !is.na(p.adjust)) %>%
    arrange(p.adjust) %>%
    slice_head(n = top_n) %>%
    mutate(
      pathway_label = clean_pathway_label(ID),
      pathway_label = factor(
        pathway_label,
        levels = pathway_label[order(NES)]
      )
    )
}

plot_df_4h <- prepare_plot_df(gsea_4h_file, top_n = 20)
plot_df_24h <- prepare_plot_df(gsea_24h_file, top_n = 20)

all_plot_df <- bind_rows(plot_df_4h, plot_df_24h)

common_max_abs_nes <- max(abs(all_plot_df$NES), na.rm = TRUE)
common_nes_lim <- c(-1.1 * common_max_abs_nes, 1.1 * common_max_abs_nes)

common_padj_lim <- range(all_plot_df$p.adjust, na.rm = TRUE)

message(
  "[INFO] Common NES axis limits: ",
  paste(round(common_nes_lim, 3), collapse = " to ")
)

message(
  "[INFO] Common p.adjust color limits: ",
  paste(signif(common_padj_lim, 3), collapse = " to ")
)

make_gsea_barplot <- function(
  plot_df,
  output_prefix,
  plot_title,
  nes_lim,
  padj_lim
) {
  if (nrow(plot_df) == 0) {
    stop("No pathways available for plotting: ", output_prefix)
  }

  p <- ggplot(plot_df, aes(x = pathway_label, y = NES, fill = p.adjust)) +
    geom_col(width = 0.8, color = "black", linewidth = 0.25) +
    coord_flip() +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    scale_y_continuous(limits = nes_lim) +
    scale_fill_gradient(
      low = "red",
      high = "blue",
      limits = padj_lim,
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

make_gsea_barplot(
  plot_df = plot_df_4h,
  output_prefix = "Inf4h_Hallmark_GSEA_Top20_by_padj",
  plot_title = "Inf4h Hallmark GSEA Top20 by padj",
  nes_lim = common_nes_lim,
  padj_lim = common_padj_lim
)

make_gsea_barplot(
  plot_df = plot_df_24h,
  output_prefix = "Inf24h_Hallmark_GSEA_Top20_by_padj",
  plot_title = "Inf24h Hallmark GSEA Top20 by padj",
  nes_lim = common_nes_lim,
  padj_lim = common_padj_lim
)

message("[INFO] GSEA NES bar plot generation completed.")
