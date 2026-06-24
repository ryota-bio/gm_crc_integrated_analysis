suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(tximport)
  library(rtracklayer)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
project_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

metadata_file <- file.path(project_dir, "metadata", "sample_metadata.csv")
salmon_dir <- file.path(project_dir, "results", "salmon_quant")
ref_dir <- file.path(project_dir, "reference")
out_dir <- file.path(project_dir, "results", "gene_matrices")
env_dir <- file.path(project_dir, "environment")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(env_dir, recursive = TRUE, showWarnings = FALSE)

message("[INFO] Project directory: ", project_dir)

metadata <- read_csv(metadata_file, show_col_types = FALSE)

required_cols <- c("sample_id", "condition", "time_point", "replicate", "fastq_R1", "fastq_R2")
missing_cols <- setdiff(required_cols, colnames(metadata))
if (length(missing_cols) > 0) {
  stop("Missing columns in sample_metadata.csv: ", paste(missing_cols, collapse = ", "))
}

files <- file.path(salmon_dir, metadata$sample_id, "quant.sf")
names(files) <- metadata$sample_id

if (!all(file.exists(files))) {
  missing_files <- files[!file.exists(files)]
  stop("Missing Salmon quant.sf files:\n", paste(missing_files, collapse = "\n"))
}

gtf_file <- file.path(ref_dir, "gencode.v45.annotation.gtf.gz")
if (!file.exists(gtf_file)) {
  stop("GTF file not found: ", gtf_file,
       "\nPlease place gencode.v45.annotation.gtf.gz in the reference/ directory.")
}

message("[INFO] Importing GTF: ", gtf_file)
gtf <- rtracklayer::import(gtf_file)

tx2gene <- as.data.frame(gtf) |>
  filter(type == "transcript") |>
  select(tx_id = transcript_id, gene_id = gene_id, gene_name = gene_name) |>
  distinct(tx_id, gene_id, gene_name)

message("[INFO] Running tximport")
txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene[, c("tx_id", "gene_id")],
  countsFromAbundance = "lengthScaledTPM",
  ignoreAfterBar = TRUE
)

counts <- as.data.frame(txi$counts) |>
  rownames_to_column("gene_id")

tpm <- as.data.frame(txi$abundance) |>
  rownames_to_column("gene_id")

gene_annot <- tx2gene |>
  select(gene_id, gene_name) |>
  distinct(gene_id, .keep_all = TRUE)

counts_out <- gene_annot |>
  right_join(counts, by = "gene_id") |>
  arrange(gene_id)

tpm_out <- gene_annot |>
  right_join(tpm, by = "gene_id") |>
  arrange(gene_id)

counts_file <- file.path(out_dir, "Gm_HCT116_salmon_gene_counts.csv")
tpm_file <- file.path(out_dir, "Gm_HCT116_salmon_gene_tpm.csv")

write_csv(counts_out, counts_file)
write_csv(tpm_out, tpm_file)

writeLines(capture.output(sessionInfo()), file.path(env_dir, "R_sessionInfo.txt"))

message("[INFO] Saved: ", counts_file)
message("[INFO] Saved: ", tpm_file)
message("[INFO] Done.")
