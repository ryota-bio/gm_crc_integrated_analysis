suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(clusterProfiler)
  library(msigdbr)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
project_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

deg_dir <- file.path(project_dir, "results", "deseq2_deg")
out_dir <- file.path(project_dir, "results", "gsea")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("[INFO] Project directory: ", project_dir)

# MSigDB Hallmark gene sets
hallmark_raw <- msigdbr(
  species = "Homo sapiens",
  collection = "H"
)

# Handle possible column-name changes across msigdbr versions
if ("gene_symbol" %in% colnames(hallmark_raw)) {
  gene_col <- "gene_symbol"
} else if ("human_gene_symbol" %in% colnames(hallmark_raw)) {
  gene_col <- "human_gene_symbol"
} else {
  stop("Could not find gene symbol column in msigdbr output.")
}

hallmark <- hallmark_raw |>
  dplyr::select(gs_name, gene_symbol = dplyr::all_of(gene_col)) |>
  dplyr::distinct(gs_name, gene_symbol)

run_gsea_hallmark <- function(deg_file, contrast_name) {
  message("[INFO] Running GSEA: ", contrast_name)

  if (!file.exists(deg_file)) {
    stop("DESeq2 result file not found: ", deg_file)
  }

  deg <- readr::read_csv(deg_file, show_col_types = FALSE) |>
    dplyr::filter(!is.na(gene_name), !is.na(stat)) |>
    dplyr::group_by(gene_name) |>
    dplyr::slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) |>
    dplyr::ungroup()

  gene_list <- deg$stat
  names(gene_list) <- deg$gene_name

  gene_list <- sort(gene_list, decreasing = TRUE)

  gsea_res <- clusterProfiler::GSEA(
    geneList = gene_list,
    TERM2GENE = hallmark,
    pvalueCutoff = 1,
    minGSSize = 10,
    maxGSSize = 500,
    eps = 0,
    verbose = FALSE
  )

  out_file <- file.path(
    out_dir,
    paste0("Gm_HCT116_GSEA_HALLMARK_", contrast_name, ".csv")
  )

  readr::write_csv(as.data.frame(gsea_res), out_file)
  message("[INFO] Saved: ", out_file)
}

run_gsea_hallmark(
  file.path(deg_dir, "Gm_HCT116_DESeq2_Inf4h_vs_Ctr.csv"),
  "Inf4h_vs_Ctr"
)

run_gsea_hallmark(
  file.path(deg_dir, "Gm_HCT116_DESeq2_Inf24h_vs_Ctr.csv"),
  "Inf24h_vs_Ctr"
)

message("[INFO] GSEA analysis completed.")
