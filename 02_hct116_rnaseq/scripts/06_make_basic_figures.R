suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(DESeq2)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
project_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

metadata_file <- file.path(project_dir, "metadata", "sample_metadata.csv")
vsd_file <- file.path(project_dir, "results", "deseq2_deg", "Gm_HCT116_vst_object.rds")
fig_dir <- file.path(project_dir, "results", "figures")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- read_csv(metadata_file, show_col_types = FALSE)
vsd <- readRDS(vsd_file)

pca_data <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p <- ggplot(pca_data, aes(PC1, PC2, color = group)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_classic(base_size = 14)

ggsave(
  filename = file.path(fig_dir, "Figure_PCA_RNAseq.pdf"),
  plot = p,
  width = 5,
  height = 4
)

message("[INFO] Figure generation completed.")
