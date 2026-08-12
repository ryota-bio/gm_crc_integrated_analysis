#!/usr/bin/env Rscript

suppressPackageStartupMessages({
library(dplyr)
library(readr)
library(stringr)
library(DESeq2)
library(pheatmap)
})

# ============================================================

# Combined Hallmark GSEA leading edge heatmaps with dendrogram

#

# - One heatmap per Hallmark pathway

# - Ctr, Inf4h, and Inf24h are shown together

# - Leading edge genes are taken as the union of Inf4h_vs_Ctr

# and Inf24h_vs_Ctr

# - Row dendrogram is shown

# - Column order is fixed as Ctr -> Inf4h -> Inf24h

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

message("[INFO] Project directory: ", project_dir)

deg_dir <- file.path(project_dir, "results", "deseq2_deg")
gsea_dir <- file.path(project_dir, "results", "gsea")
out_dir <- file.path(project_dir, "results", "figures", "gsea_leading_edge_heatmaps")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------

# Input files

# -----------------------------

deg_files <- list(
Inf4h_vs_Ctr = file.path(deg_dir, "Gm_HCT116_DESeq2_Inf4h_vs_Ctr.csv"),
Inf24h_vs_Ctr = file.path(deg_dir, "Gm_HCT116_DESeq2_Inf24h_vs_Ctr.csv")
)

gsea_files <- list(
Inf4h_vs_Ctr = file.path(gsea_dir, "Gm_HCT116_GSEA_HALLMARK_Inf4h_vs_Ctr.csv"),
Inf24h_vs_Ctr = file.path(gsea_dir, "Gm_HCT116_GSEA_HALLMARK_Inf24h_vs_Ctr.csv")
)

vst_file <- file.path(deg_dir, "Gm_HCT116_vst_object.rds")

for (f in c(deg_files, gsea_files, vst_file)) {
if (!file.exists(f)) {
stop("File not found: ", f)
}
}

# -----------------------------

# Pathways to plot

# -----------------------------

selected_pathways <- c(
"HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
"HALLMARK_HYPOXIA",
"HALLMARK_TNFA_SIGNALING_VIA_NFKB",
"HALLMARK_GLYCOLYSIS"
)

# -----------------------------

# Plot parameters

# -----------------------------

cluster_rows_flag <- TRUE
cluster_cols_flag <- FALSE

show_rownames_flag <- TRUE

fontsize_row_value <- 8
fontsize_col_value <- 10

treeheight_row_value <- 60
treeheight_col_value <- ifelse(cluster_cols_flag, 50, 0)

border_color_value <- "grey75"

sample_group_order <- c("Ctr", "Inf4h", "Inf24h")

max_genes_per_pathway <- 30

heatmap_colors <- colorRampPalette(
c("#4575B4", "#F7F7F7", "#D73027")
)(100)

annotation_colors <- list(
Group = c(
Ctr = "#2ECC71",
Inf4h = "#F1948A",
Inf24h = "#85C1E9"
)
)

# -----------------------------

# Helper functions

# -----------------------------

row_zscore <- function(mat) {
z <- t(scale(t(mat)))
z[is.na(z)] <- 0
z
}

sanitize_filename <- function(x) {
x %>%
str_replace("^HALLMARK_", "") %>%
str_to_lower() %>%
str_replace_all("[^a-z0-9]+", "_") %>%
str_replace_all("^_|_$", "")
}

clean_pathway_name <- function(x) {
x %>%
str_replace("^HALLMARK_", "") %>%
str_replace_all("_", " ") %>%
str_to_title()
}

get_core_genes <- function(gsea_df, pathway_id) {
pathway_row <- gsea_df %>%
filter(ID == pathway_id)

if (nrow(pathway_row) == 0) {
return(character(0))
}

if (!"core_enrichment" %in% colnames(pathway_row)) {
return(character(0))
}

core_str <- pathway_row$core_enrichment[1]

if (is.na(core_str) || core_str == "") {
return(character(0))
}

genes <- str_split(core_str, "/", simplify = FALSE)[[1]]
genes <- unique(genes)
genes <- genes[genes != ""]

genes
}

detect_gene_id_column <- function(deg_df) {
candidates <- c("gene_id", "Geneid", "gene", "X", "X1", "...1")

for (candidate in candidates) {
if (candidate %in% colnames(deg_df)) {
return(candidate)
}
}

stop("Could not find a gene ID column in DEG table.")
}

make_gene_mapping <- function(deg_df_all, vst_mat) {
gene_id_col <- detect_gene_id_column(deg_df_all)

mapping <- deg_df_all %>%
filter(
!is.na(gene_name),
gene_name != "",
!is.na(.data[[gene_id_col]])
) %>%
transmute(
gene_name = as.character(gene_name),
gene_id = as.character(.data[[gene_id_col]])
) %>%
distinct() %>%
mutate(
gene_id_nover = sub("\\..*$", "", gene_id)
)

vst_map <- tibble(
vst_id = rownames(vst_mat),
gene_id_nover = sub("\\..*$", "", rownames(vst_mat))
)

mapping %>%
left_join(vst_map, by = "gene_id_nover") %>%
filter(!is.na(vst_id)) %>%
distinct(gene_name, vst_id, .keep_all = TRUE)
}

make_heatmap <- function(mat_plot, annotation_col, pathway_id) {
pathway_label <- clean_pathway_name(pathway_id)
pathway_file <- sanitize_filename(pathway_id)

file_stub <- paste0(
"CombinedTop30_Hallmark_LeadingEdge_",
pathway_file
)

pdf_file <- file.path(out_dir, paste0(file_stub, ".pdf"))
png_file <- file.path(out_dir, paste0(file_stub, ".png"))
gene_file <- file.path(out_dir, paste0(file_stub, "_genes.csv"))

write_csv(
data.frame(gene_name = rownames(mat_plot)),
gene_file
)

plot_width <- max(7.5, ncol(mat_plot) * 0.55 + 2.5)
plot_height <- max(6.5, nrow(mat_plot) * 0.22 + 2.5)

pheatmap_args <- list(
mat = mat_plot,
color = heatmap_colors,
cluster_rows = cluster_rows_flag,
cluster_cols = cluster_cols_flag,
annotation_col = annotation_col,
annotation_colors = annotation_colors,
show_rownames = show_rownames_flag,
show_colnames = TRUE,
fontsize_row = fontsize_row_value,
fontsize_col = fontsize_col_value,
border_color = border_color_value,
angle_col = 45,
treeheight_row = treeheight_row_value,
treeheight_col = treeheight_col_value,
main = paste0(pathway_label, " top 30 leading edge genes")
)

pdf(pdf_file, width = plot_width, height = plot_height)
do.call(pheatmap, pheatmap_args)
dev.off()

png(png_file, width = plot_width, height = plot_height, units = "in", res = 300)
do.call(pheatmap, pheatmap_args)
dev.off()

message("[INFO] Saved: ", pdf_file)
message("[INFO] Saved: ", png_file)
message("[INFO] Saved: ", gene_file)
}

# -----------------------------

# Load data

# -----------------------------

deg_list <- lapply(deg_files, function(x) {
read_csv(x, show_col_types = FALSE)
})

gsea_list <- lapply(gsea_files, function(x) {
read_csv(x, show_col_types = FALSE)
})

vst_obj <- readRDS(vst_file)
vst_mat <- assay(vst_obj)
sample_info <- as.data.frame(colData(vst_obj))
sample_info$sample <- rownames(sample_info)

if (!"group" %in% colnames(sample_info)) {
stop("The VST object colData must contain a 'group' column.")
}

sample_info <- sample_info %>%
mutate(
group = factor(group, levels = sample_group_order)
) %>%
arrange(group, sample)

sample_order <- sample_info$sample
sample_order <- sample_order[sample_order %in% colnames(vst_mat)]

vst_mat <- vst_mat[, sample_order, drop = FALSE]
sample_info <- sample_info %>%
filter(sample %in% sample_order)

annotation_col <- data.frame(
Group = sample_info$group
)
rownames(annotation_col) <- sample_info$sample

message("[INFO] Sample groups found:")
print(table(sample_info$group))

message("[INFO] Samples used:")
print(sample_order)

deg_df_all <- bind_rows(deg_list)
gene_mapping <- make_gene_mapping(deg_df_all, vst_mat)

# Rank genes for main-figure heatmaps.
# Priority is based on Inf24h_vs_Ctr |stat| when available.
# If not available, the maximum absolute score across contrasts is used.
priority_df <- bind_rows(lapply(names(deg_list), function(comparison_name) {
  df <- deg_list[[comparison_name]]

  if ("stat" %in% colnames(df)) {
    rank_col <- "stat"
  } else {
    rank_col <- "log2FoldChange"
  }

  df %>%
    filter(
      !is.na(gene_name),
      gene_name != "",
      !is.na(.data[[rank_col]])
    ) %>%
    transmute(
      gene_name = as.character(gene_name),
      comparison = comparison_name,
      rank_score = as.numeric(.data[[rank_col]]),
      abs_rank_score = abs(rank_score)
    )
})) %>%
  group_by(gene_name) %>%
  summarise(
    inf24h_score = suppressWarnings(max(abs_rank_score[comparison == "Inf24h_vs_Ctr"], na.rm = TRUE)),
    max_score = max(abs_rank_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    inf24h_score = ifelse(is.infinite(inf24h_score), NA_real_, inf24h_score),
    priority_score = ifelse(is.na(inf24h_score), max_score, inf24h_score)
  )

# -----------------------------

# Main loop

# -----------------------------

for (pathway_id in selected_pathways) {
message("[INFO] Processing pathway: ", pathway_id)

leading_genes <- unique(unlist(lapply(gsea_list, function(gsea_df) {
get_core_genes(gsea_df, pathway_id)
})))

leading_genes <- leading_genes[leading_genes != ""]

if (length(leading_genes) == 0) {
message("[INFO] No leading edge genes found for: ", pathway_id)
next
}

pathway_mapping <- gene_mapping %>%
filter(gene_name %in% leading_genes) %>%
distinct(gene_name, vst_id)

if (nrow(pathway_mapping) == 0) {
message("[INFO] No VST-matched genes found for: ", pathway_id)
next
}

pathway_mapping <- pathway_mapping %>%
filter(vst_id %in% rownames(vst_mat))

if (nrow(pathway_mapping) == 0) {
message("[INFO] No VST rownames matched for: ", pathway_id)
next
}

pathway_mapping <- pathway_mapping %>%
slice_head(n = max_genes_per_pathway)

message("[INFO] Number of genes plotted for ", pathway_id, ": ", nrow(pathway_mapping))

mat_sub <- vst_mat[pathway_mapping$vst_id, , drop = FALSE]

rownames(mat_sub) <- pathway_mapping$gene_name[
match(rownames(mat_sub), pathway_mapping$vst_id)
]

rownames(mat_sub) <- make.unique(rownames(mat_sub))

mat_plot <- row_zscore(mat_sub)

keep_rows <- apply(mat_plot, 1, function(x) {
all(!is.na(x)) && any(abs(x) > 1e-12)
})

mat_plot <- mat_plot[keep_rows, , drop = FALSE]

if (nrow(mat_plot) < 2) {
message("[INFO] Too few genes to draw heatmap for: ", pathway_id)
next
}

make_heatmap(
mat_plot = mat_plot,
annotation_col = annotation_col,
pathway_id = pathway_id
)
}

message("[INFO] Combined Hallmark leading edge heatmap generation completed.")
