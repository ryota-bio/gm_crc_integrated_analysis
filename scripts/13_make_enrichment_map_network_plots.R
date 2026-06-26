#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(igraph)
})

# ============================================================
# Enrichment map and cnet-like network plots
#
# This script generates:
#   1. Enrichment map:
#      enriched terms/pathways are connected based on gene overlap
#
#   2. Cnet-like plot:
#      enriched terms/pathways are connected to contributing genes
#
# Input:
#   results/gsea/Gm_HCT116_GSEA_HALLMARK_*.csv
#   results/ora/*_ORA.csv
#
# Output:
#   results/figures/enrichment_networks/*.pdf
#   results/figures/enrichment_networks/*.png
#   results/figures/enrichment_networks/*_nodes.csv
#   results/figures/enrichment_networks/*_edges.csv
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

gsea_dir <- file.path(project_dir, "results", "gsea")
ora_dir <- file.path(project_dir, "results", "ora")
out_dir <- file.path(project_dir, "results", "figures", "enrichment_networks")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Parameters
# -----------------------------

top_n_terms_emap <- 15
top_n_terms_cnet <- 8
max_genes_cnet <- 50

padj_cutoff <- 0.05
min_jaccard <- 0.10

label_cex_term <- 0.75
label_cex_gene <- 0.55

# -----------------------------
# Helper functions
# -----------------------------

safe_filename <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace_all("^_|_$", "") %>%
    str_to_lower()
}

wrap_label <- function(x, width = 28) {
  vapply(
    x,
    function(s) paste(strwrap(s, width = width), collapse = "\n"),
    character(1)
  )
}

clean_term_label <- function(x, width = 28) {
  x <- x %>%
    str_replace("^HALLMARK_", "") %>%
    str_replace_all("_", " ") %>%
    str_to_title()

  x <- x %>%
    str_replace("Tnfa Signaling Via Nfkb", "Tnfa Signaling
Via Nfkb") %>%
    str_replace("Epithelial Mesenchymal Transition", "Epithelial Mesenchymal
Transition") %>%
    str_replace("Cholesterol Homeostasis", "Cholesterol
Homeostasis") %>%
    str_replace("Unfolded Protein Response", "Unfolded Protein
Response") %>%
    str_replace("Inflammatory Response", "Inflammatory
Response") %>%
    str_replace("Myc Targets V1", "Myc Targets
V1") %>%
    str_replace("Myc Targets V2", "Myc Targets
V2") %>%
    str_replace("E2f Targets", "E2f
Targets") %>%
    str_replace("P53 Pathway", "P53
Pathway") %>%
    str_replace("Glycolysis", "Glycolysis") %>%
    str_replace("Hypoxia", "Hypoxia")

  wrap_label(x, width = width)
}

parse_gene_list <- function(x) {
  if (is.na(x) || x == "") {
    return(character(0))
  }
  genes <- unlist(str_split(x, "/"))
  genes <- unique(genes)
  genes <- genes[genes != ""]
  genes
}

get_gene_column <- function(df) {
  if ("core_enrichment" %in% colnames(df)) {
    return("core_enrichment")
  }
  if ("geneID" %in% colnames(df)) {
    return("geneID")
  }
  stop("No gene column found. Expected core_enrichment or geneID.")
}

get_description_column <- function(df) {
  if ("Description" %in% colnames(df)) {
    return("Description")
  }
  if ("ID" %in% colnames(df)) {
    return("ID")
  }
  stop("No Description or ID column found.")
}

get_padj_column <- function(df) {
  if ("p.adjust" %in% colnames(df)) {
    return("p.adjust")
  }
  if ("padj" %in% colnames(df)) {
    return("padj")
  }
  stop("No adjusted P-value column found.")
}

prepare_enrichment_df <- function(file_path, source_type) {
  df <- readr::read_csv(file_path, show_col_types = FALSE)

  gene_col <- get_gene_column(df)
  desc_col <- get_description_column(df)
  padj_col <- get_padj_column(df)

  if (!("ID" %in% colnames(df))) {
    df$ID <- df[[desc_col]]
  }

  if (!("NES" %in% colnames(df))) {
    df$NES <- NA_real_
  }

  out <- df %>%
    mutate(
      source_type = source_type,
      term_id = as.character(ID),
      term_name = as.character(.data[[desc_col]]),
      p_adjust = as.numeric(.data[[padj_col]]),
      gene_string = as.character(.data[[gene_col]])
    ) %>%
    filter(!is.na(p_adjust)) %>%
    mutate(
      gene_list = lapply(gene_string, parse_gene_list),
      gene_count = lengths(gene_list)
    ) %>%
    filter(gene_count > 0) %>%
    arrange(p_adjust)

  sig <- out %>%
    filter(p_adjust < padj_cutoff)

  if (nrow(sig) >= 3) {
    out <- sig
  }

  out
}

make_term_edges <- function(term_df, min_jaccard = 0.20) {
  if (nrow(term_df) < 2) {
    return(tibble())
  }

  edge_list <- list()
  k <- 1

  for (i in seq_len(nrow(term_df) - 1)) {
    for (j in seq((i + 1), nrow(term_df))) {
      genes_i <- term_df$gene_list[[i]]
      genes_j <- term_df$gene_list[[j]]

      inter <- length(intersect(genes_i, genes_j))
      union_n <- length(union(genes_i, genes_j))

      if (union_n == 0) {
        next
      }

      jaccard <- inter / union_n

      if (jaccard >= min_jaccard) {
        edge_list[[k]] <- tibble(
          from = term_df$term_id[i],
          to = term_df$term_id[j],
          overlap = inter,
          jaccard = jaccard
        )
        k <- k + 1
      }
    }
  }

  if (length(edge_list) == 0) {
    return(tibble())
  }

  bind_rows(edge_list)
}

plot_enrichment_map <- function(term_df, file_stub) {
  term_df <- term_df %>%
    slice_head(n = top_n_terms_emap)

  edges <- make_term_edges(term_df, min_jaccard = min_jaccard)

  nodes <- term_df %>%
    transmute(
      name = term_id,
      label = clean_term_label(term_name, width = 28),
      term_name = term_name,
      p_adjust = p_adjust,
      gene_count = gene_count,
      neg_log10_padj = -log10(p_adjust)
    )

  nodes_file <- file.path(out_dir, paste0(file_stub, "_emap_nodes.csv"))
  edges_file <- file.path(out_dir, paste0(file_stub, "_emap_edges.csv"))

  write_csv(nodes, nodes_file)
  write_csv(edges, edges_file)

  if (nrow(edges) == 0) {
    message("[INFO] No edges found for enrichment map: ", file_stub)
    return(invisible(NULL))
  }

  g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)

  set.seed(123)
  layout <- layout_with_fr(g)
  layout <- layout * 1.35

  v_size <- 8 + 2.5 * sqrt(V(g)$gene_count)
  v_color_value <- V(g)$neg_log10_padj
  v_color_scaled <- (v_color_value - min(v_color_value)) /
    (max(v_color_value) - min(v_color_value) + 1e-9)
  v_colors <- colorRampPalette(c("#91BFDB", "#FFFFBF", "#FC8D59"))(100)[
    pmax(1, round(v_color_scaled * 99) + 1)
  ]

  e_width <- 1 + 8 * E(g)$jaccard

  pdf_file <- file.path(out_dir, paste0(file_stub, "_emap.pdf"))
  png_file <- file.path(out_dir, paste0(file_stub, "_emap.png"))

  pdf(pdf_file, width = 14, height = 10)
  plot(
    g,
    layout = layout,
    vertex.size = v_size,
    vertex.color = v_colors,
    vertex.frame.color = "grey30",
    vertex.label = V(g)$label,
    vertex.label.cex = label_cex_term,
    vertex.label.color = "black",
    vertex.label.dist = 1.1,
    vertex.label.degree = -pi / 2,
    edge.width = e_width,
    edge.color = adjustcolor("grey40", alpha.f = 0.55),
    main = paste0(file_stub, " enrichment map")
  )
  dev.off()

  png(png_file, width = 14, height = 10, units = "in", res = 300)
  plot(
    g,
    layout = layout,
    vertex.size = v_size,
    vertex.color = v_colors,
    vertex.frame.color = "grey30",
    vertex.label = V(g)$label,
    vertex.label.cex = label_cex_term,
    vertex.label.color = "black",
    vertex.label.dist = 1.1,
    vertex.label.degree = -pi / 2,
    edge.width = e_width,
    edge.color = adjustcolor("grey40", alpha.f = 0.55),
    main = paste0(file_stub, " enrichment map")
  )
  dev.off()

  message("[INFO] Saved: ", pdf_file)
  message("[INFO] Saved: ", png_file)
  message("[INFO] Saved: ", nodes_file)
  message("[INFO] Saved: ", edges_file)
}

plot_cnet <- function(term_df, file_stub) {
  term_df <- term_df %>%
    slice_head(n = top_n_terms_cnet)

  term_gene <- term_df %>%
    select(term_id, term_name, p_adjust, gene_list) %>%
    tidyr::unnest(gene_list) %>%
    rename(gene = gene_list)

  if (nrow(term_gene) == 0) {
    message("[INFO] No term-gene edges found for cnet: ", file_stub)
    return(invisible(NULL))
  }

  gene_rank <- term_gene %>%
    count(gene, name = "degree") %>%
    arrange(desc(degree), gene) %>%
    slice_head(n = max_genes_cnet)

  term_gene <- term_gene %>%
    filter(gene %in% gene_rank$gene)

  term_nodes <- term_df %>%
    filter(term_id %in% unique(term_gene$term_id)) %>%
    transmute(
      name = term_id,
      label = clean_term_label(term_name, width = 25),
      type = "term",
      p_adjust = p_adjust,
      degree = NA_integer_
    )

  gene_nodes <- gene_rank %>%
    filter(gene %in% unique(term_gene$gene)) %>%
    transmute(
      name = gene,
      label = gene,
      type = "gene",
      p_adjust = NA_real_,
      degree = degree
    )

  nodes <- bind_rows(term_nodes, gene_nodes)

  edges <- term_gene %>%
    transmute(
      from = term_id,
      to = gene
    ) %>%
    distinct()

  nodes_file <- file.path(out_dir, paste0(file_stub, "_cnet_nodes.csv"))
  edges_file <- file.path(out_dir, paste0(file_stub, "_cnet_edges.csv"))

  write_csv(nodes, nodes_file)
  write_csv(edges, edges_file)

  g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)

  set.seed(123)
  layout <- layout_with_fr(g)
  layout <- layout * 1.35

  v_type <- V(g)$type
  v_size <- ifelse(v_type == "term", 18, 7)
  v_color <- ifelse(v_type == "term", "#FC8D59", "#91BFDB")
  v_label_cex <- ifelse(v_type == "term", label_cex_term, label_cex_gene)

  # Default:
  # Put all labels at the center of each node.
  v_label_dist <- rep(0, vcount(g))
  v_label_degree <- rep(0, vcount(g))

  # Fine adjustment for Inf4h cnet plot:
  # Move only EPITHELIAL_MESENCHYMAL_TRANSITION and Myc Targets V1 labels downward.
  if (str_detect(file_stub, "inf4h_vs_ctr")) {
    special_down_nodes <- V(g)$name %in% c(
      "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
      "HALLMARK_MYC_TARGETS_V1"
    )

    v_label_dist[special_down_nodes] <- 0.65
    v_label_degree[special_down_nodes] <- pi / 2
  }

  pdf_file <- file.path(out_dir, paste0(file_stub, "_cnet.pdf"))
  png_file <- file.path(out_dir, paste0(file_stub, "_cnet.png"))

  pdf(pdf_file, width = 18, height = 13)
  plot(
    g,
    layout = layout,
    vertex.size = v_size,
    vertex.color = v_color,
    vertex.frame.color = "grey30",
    vertex.label = V(g)$label,
    vertex.label.cex = v_label_cex,
    vertex.label.color = "black",
    vertex.label.dist = v_label_dist,
    vertex.label.degree = v_label_degree,
    edge.width = 0.8,
    edge.color = adjustcolor("grey50", alpha.f = 0.45),
    main = paste0(file_stub, " cnet-like plot")
  )
  dev.off()

  png(png_file, width = 18, height = 13, units = "in", res = 300)
  plot(
    g,
    layout = layout,
    vertex.size = v_size,
    vertex.color = v_color,
    vertex.frame.color = "grey30",
    vertex.label = V(g)$label,
    vertex.label.cex = v_label_cex,
    vertex.label.color = "black",
    vertex.label.dist = v_label_dist,
    vertex.label.degree = v_label_degree,
    edge.width = 0.8,
    edge.color = adjustcolor("grey50", alpha.f = 0.45),
    main = paste0(file_stub, " cnet-like plot")
  )
  dev.off()

  message("[INFO] Saved: ", pdf_file)
  message("[INFO] Saved: ", png_file)
  message("[INFO] Saved: ", nodes_file)
  message("[INFO] Saved: ", edges_file)
}

# -----------------------------
# Input file list
# -----------------------------

input_files <- list()

gsea_files <- list.files(
  gsea_dir,
  pattern = "GSEA_HALLMARK_.*\\.csv$",
  full.names = TRUE
)

for (f in gsea_files) {
  nm <- basename(f) %>%
    str_replace("^Gm_HCT116_", "") %>%
    str_replace("\\.csv$", "")
  input_files[[nm]] <- list(path = f, type = "GSEA_HALLMARK")
}

ora_files <- character(0)
# ORA files are intentionally excluded. This script now processes GSEA Hallmark results only.

if (length(input_files) == 0) {
  stop("No GSEA or ORA result files found.")
}

# -----------------------------
# Main
# -----------------------------

for (nm in names(input_files)) {
  message("[INFO] Processing: ", nm)

  file_path <- input_files[[nm]]$path
  source_type <- input_files[[nm]]$type

  enrich_df <- tryCatch(
    {
      prepare_enrichment_df(file_path, source_type)
    },
    error = function(e) {
      message("[INFO] Skipped: ", nm, " ; reason: ", e$message)
      return(NULL)
    }
  )

  if (is.null(enrich_df)) {
    next
  }

  if (nrow(enrich_df) < 2) {
    message("[INFO] Too few enriched terms for network: ", nm)
    next
  }

  file_stub <- safe_filename(nm)

  plot_enrichment_map(enrich_df, file_stub)
  plot_cnet(enrich_df, file_stub)
}

message("[INFO] Enrichment map and network plot generation completed.")

