#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

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
out_dir <- file.path(project_dir, "results", "figures", "volcano_plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

deg_files <- list(
  Inf4h_vs_Ctr  = file.path(deg_dir, "Gm_HCT116_DESeq2_Inf4h_vs_Ctr.csv"),
  Inf24h_vs_Ctr = file.path(deg_dir, "Gm_HCT116_DESeq2_Inf24h_vs_Ctr.csv")
)

for (f in deg_files) {
  if (!file.exists(f)) {
    stop("File not found: ", f)
  }
}

padj_cutoff <- 0.05
log2fc_cutoff <- 1
y_axis_cap <- 300
highlight_symbols <- c("SERPINE1", "PLAUR", "CCN1", "CYR61", "SDC4")

prepare_volcano_df <- function(input_file) {
  df_raw <- read_csv(input_file, show_col_types = FALSE)

  required_cols <- c("gene_name", "log2FoldChange", "padj")
  missing_cols <- setdiff(required_cols, colnames(df_raw))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ", input_file, ": ",
         paste(missing_cols, collapse = ", "))
  }

  df <- df_raw %>%
    filter(!is.na(log2FoldChange), !is.na(padj)) %>%
    mutate(
      padj_plot = ifelse(padj <= 0, 1e-300, padj),
      neglog10_padj_raw = -log10(padj_plot),
      neglog10_padj = pmin(neglog10_padj_raw, y_axis_cap),
      regulation = case_when(
        padj < padj_cutoff & log2FoldChange >= log2fc_cutoff  ~ "Up",
        padj < padj_cutoff & log2FoldChange <= -log2fc_cutoff ~ "Down",
        TRUE ~ "NS"
      )
    )

  return(df)
}

# ---- 共通x軸範囲を先に計算 ----
all_df <- lapply(deg_files, prepare_volcano_df)
all_fc <- unlist(lapply(all_df, function(x) x$log2FoldChange))
x_limit <- ceiling(max(abs(all_fc), na.rm = TRUE))
x_limit <- x_limit + 0.5

message("[INFO] Common x-axis limit: ", -x_limit, " to ", x_limit)

make_volcano_plot <- function(df, comparison_name) {
  label_df <- df %>%
    filter(gene_name %in% highlight_symbols) %>%
    mutate(
      label = case_when(
        gene_name == "CYR61" ~ "CCN1",
        TRUE ~ gene_name
      )
    ) %>%
    distinct(label, .keep_all = TRUE)

  message("[INFO] ", comparison_name, ": genes plotted = ", nrow(df))
  message("[INFO] ", comparison_name, ": highlighted genes found = ",
          paste(label_df$label, collapse = ", "))

  p <- ggplot(df, aes(x = log2FoldChange, y = neglog10_padj)) +
    geom_point(aes(color = regulation), size = 1.2, alpha = 0.75) +
    scale_color_manual(
      values = c(
        "Down" = "blue",
        "NS"   = "grey75",
        "Up"   = "red"
      ),
      breaks = c("Down", "NS", "Up")
    ) +
    geom_vline(
      xintercept = c(-log2fc_cutoff, log2fc_cutoff),
      linetype = "dashed",
      color = "black",
      linewidth = 0.5
    ) +
    geom_hline(
      yintercept = -log10(padj_cutoff),
      linetype = "dashed",
      color = "black",
      linewidth = 0.5
    ) +
    # 強調表示：黒枠
    geom_point(
      data = label_df,
      aes(fill = regulation),
      shape = 21,
      color = "black",
      size = 3.8,
      stroke = 0.9,
      inherit.aes = TRUE
    ) +
    scale_fill_manual(
      values = c(
        "Down" = "blue",
        "NS"   = "grey75",
        "Up"   = "red"
      ),
      guide = "none"
    ) +
    geom_text_repel(
      data = label_df,
      aes(label = label),
      size = 5,
      box.padding = 0.35,
      point.padding = 0.25,
      segment.color = "black",
      segment.size = 0.3,
      max.overlaps = Inf,
      min.segment.length = 0
    ) +
    coord_cartesian(
      xlim = c(-x_limit, x_limit),
      ylim = c(0, y_axis_cap)
    ) +
    labs(
      title = paste0(comparison_name, " Volcano plot"),
      x = "log2 fold change",
      y = expression(-log[10]("adjusted P value")),
      color = NULL
    ) +
    theme_bw(base_size = 16) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 22),
      legend.position = "top",
      legend.text = element_text(size = 14),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 14),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank()
    )

  file_stub <- paste0("Volcano_", comparison_name)

  pdf_file <- file.path(out_dir, paste0(file_stub, ".pdf"))
  png_file <- file.path(out_dir, paste0(file_stub, ".png"))
  label_file <- file.path(out_dir, paste0(file_stub, "_highlighted_genes.csv"))

  ggsave(pdf_file, p, width = 9, height = 7)
  ggsave(png_file, p, width = 9, height = 7, dpi = 300)
  write_csv(label_df, label_file)

  message("[INFO] Saved: ", pdf_file)
  message("[INFO] Saved: ", png_file)
  message("[INFO] Saved: ", label_file)
}

make_volcano_plot(all_df[["Inf4h_vs_Ctr"]], "Inf4h_vs_Ctr")
make_volcano_plot(all_df[["Inf24h_vs_Ctr"]], "Inf24h_vs_Ctr")

message("[INFO] Volcano plot generation completed.")
