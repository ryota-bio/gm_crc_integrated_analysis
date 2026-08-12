#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
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

ora_dir <- file.path(project_dir, "results", "ora")
fig_dir <- file.path(project_dir, "results", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

parse_gene_ratio <- function(x) {
  x <- as.character(x)
  ifelse(
    grepl("/", x),
    as.numeric(sub("/.*", "", x)) / as.numeric(sub(".*/", "", x)),
    as.numeric(x)
  )
}

make_go_dotplot <- function(input_file,
                            output_prefix,
                            plot_title,
                            top_n = 15,
                            x_axis_max = 0.20) {
  if (!file.exists(input_file)) {
    warning("[WARN] File not found: ", input_file)
    return(NULL)
  }

  df <- readr::read_csv(input_file, show_col_types = FALSE)

  required_cols <- c("Description", "GeneRatio", "Count", "p.adjust")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    warning("[WARN] Missing required columns in ", input_file, ": ",
            paste(missing_cols, collapse = ", "))
    return(NULL)
  }

  plot_df <- df %>%
    filter(!is.na(GeneRatio), !is.na(Count), !is.na(p.adjust)) %>%
    mutate(
      GeneRatio_numeric = parse_gene_ratio(GeneRatio),
      Description = str_to_sentence(Description)
    ) %>%
    arrange(desc(Count), p.adjust) %>%
    slice_head(n = top_n)

  plot_df$Description <- factor(
    plot_df$Description,
    levels = plot_df$Description[order(plot_df$Count)]
  )

  p <- ggplot(
    plot_df,
    aes(
      x = GeneRatio_numeric,
      y = Description,
      size = Count,
      fill = p.adjust
    )
  ) +
    geom_point(
      shape = 21,
      color = "black",
      stroke = 0.5,
      alpha = 1
    ) +
    scale_fill_gradient(
      low = "red",
      high = "blue",
      name = "p.adjust",
      limits = c(0, 0.05),
      breaks = c(0.01, 0.02, 0.03, 0.04)
    ) +
    scale_size_continuous(
      name = "Count",
      breaks = c(20, 40, 60),
      range = c(2, 8)
    ) +
    scale_x_continuous(
      limits = c(0, x_axis_max),
      breaks = seq(0, x_axis_max, by = 0.05)
    ) +
    labs(
      title = plot_title,
      x = "GeneRatio",
      y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 11),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    ) +
    guides(
      fill = guide_colorbar(order = 1),
      size = guide_legend(
        order = 2,
        override.aes = list(
          shape = 21,
          fill = "white",
          color = "black",
          stroke = 0.5
        )
      )
    )

  pdf_file <- file.path(fig_dir, paste0(output_prefix, ".pdf"))
  png_file <- file.path(fig_dir, paste0(output_prefix, ".png"))
  csv_file <- file.path(fig_dir, paste0(output_prefix, "_data.csv"))

  ggsave(pdf_file, p, width = 7.2, height = 5.8)
  ggsave(png_file, p, width = 7.2, height = 5.8, dpi = 300)
  readr::write_csv(plot_df, csv_file)

  message("[INFO] Saved: ", pdf_file)
  message("[INFO] Saved: ", png_file)
  message("[INFO] Saved: ", csv_file)
}

make_go_dotplot(file.path(ora_dir, "Inf4h_vs_Ctr_GO_BP_up_ORA.csv"),
                "Inf4h_vs_Ctr_GO_BP_up_byCount_dotplot",
                "Inf4h vs Ctr GO BP up")

make_go_dotplot(file.path(ora_dir, "Inf4h_vs_Ctr_GO_BP_down_ORA.csv"),
                "Inf4h_vs_Ctr_GO_BP_down_byCount_dotplot",
                "Inf4h vs Ctr GO BP down")

make_go_dotplot(file.path(ora_dir, "Inf24h_vs_Ctr_GO_BP_up_ORA.csv"),
                "Inf24h_vs_Ctr_GO_BP_up_byCount_dotplot",
                "Inf24h vs Ctr GO BP up")

make_go_dotplot(file.path(ora_dir, "Inf24h_vs_Ctr_GO_BP_down_ORA.csv"),
                "Inf24h_vs_Ctr_GO_BP_down_byCount_dotplot",
                "Inf24h vs Ctr GO BP down")

message("[INFO] GO BP ORA dot plot generation completed.")
