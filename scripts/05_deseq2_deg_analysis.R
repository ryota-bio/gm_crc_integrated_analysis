suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(DESeq2)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
project_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

metadata_file <- file.path(project_dir, "metadata", "sample_metadata.csv")
counts_file <- file.path(project_dir, "results", "gene_matrices", "Gm_HCT116_salmon_gene_counts.csv")
out_dir <- file.path(project_dir, "results", "deseq2_deg")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- read_csv(metadata_file, show_col_types = FALSE)
counts_df <- read_csv(counts_file, show_col_types = FALSE)

gene_info <- counts_df |> select(gene_id, gene_name)

count_mat <- counts_df |>
  select(-gene_name) |>
  column_to_rownames("gene_id") |>
  as.matrix()

count_mat <- round(count_mat)

metadata <- metadata |>
  mutate(
    group = case_when(
      condition == "Control" ~ "Ctr",
      condition == "Gm_infected" & time_point == "4h" ~ "Inf4h",
      condition == "Gm_infected" & time_point == "24h" ~ "Inf24h",
      TRUE ~ NA_character_
    )
  )

metadata <- metadata |> filter(sample_id %in% colnames(count_mat))
metadata$group <- factor(metadata$group, levels = c("Ctr", "Inf4h", "Inf24h"))
rownames(metadata) <- metadata$sample_id

count_mat <- count_mat[, metadata$sample_id]

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = metadata,
  design = ~ group
)

keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

dds <- DESeq(dds)

write_deg <- function(contrast_name, res_obj) {
  res_df <- as.data.frame(res_obj) |>
    rownames_to_column("gene_id") |>
    left_join(gene_info, by = "gene_id") |>
    select(gene_id, gene_name, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) |>
    arrange(padj)

  out_file <- file.path(out_dir, paste0("Gm_HCT116_DESeq2_", contrast_name, ".csv"))
  write_csv(res_df, out_file)
  message("[INFO] Saved: ", out_file)
}

res_4h <- results(dds, contrast = c("group", "Inf4h", "Ctr"))
res_24h <- results(dds, contrast = c("group", "Inf24h", "Ctr"))

write_deg("Inf4h_vs_Ctr", res_4h)
write_deg("Inf24h_vs_Ctr", res_24h)

vsd <- vst(dds, blind = FALSE)
saveRDS(vsd, file.path(out_dir, "Gm_HCT116_vst_object.rds"))

message("[INFO] DESeq2 analysis completed.")
