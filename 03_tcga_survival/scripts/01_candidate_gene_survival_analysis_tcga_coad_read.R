#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(survival)
})
if (!requireNamespace("survminer", quietly = TRUE)) {
  stop("Package 'survminer' is required for Kaplan-Meier plots with risk tables.", call. = FALSE)
}

resolve_project_dir <- function() {
  cwd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (basename(cwd) == "tcga_crc_candidate_gene_survival") return(cwd)
  if (basename(cwd) == "scripts" && basename(dirname(cwd)) == "tcga_crc_candidate_gene_survival") {
    return(normalizePath(dirname(cwd), winslash = "/", mustWork = TRUE))
  }
  stop(
    paste0(
      "Cannot identify repository root from working directory: ", cwd,
      ". Run this script from the tcga_crc_candidate_gene_survival repository root ",
      "or its scripts/ directory."
    ),
    call. = FALSE
  )
}

project_dir <- resolve_project_dir()
expression_candidates <- file.path(
  project_dir, "data", "expression",
  c("tcga_coad_read_vst_protein_coding.rds", "tcga_coad_read_vst_protein_coding.csv")
)
annotation_path <- file.path(
  project_dir, "data", "expression", "tcga_coad_read_gene_annotation_protein_coding.csv"
)
sample_map_path <- file.path(project_dir, "data", "expression", "sample_file_map.csv")
clinical_path <- file.path(
  project_dir, "data", "clinical", "cbioportal_coadread_tcga_pub_clinical_data.tsv"
)
survival_dir <- file.path(project_dir, "results", "survival")
figure_dir <- file.path(project_dir, "figures")
dir.create(survival_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

sample_out <- file.path(survival_dir, "candidate_genes_tcga_coad_read_survival_sample_table.csv")
summary_out <- file.path(survival_dir, "candidate_genes_tcga_coad_read_survival_summary.csv")
tests_out <- file.path(survival_dir, "candidate_genes_tcga_coad_read_survival_tests.csv")
multivariable_out <- file.path(survival_dir, "candidate_genes_tcga_coad_read_multivariable_cox_results.csv")
ph_out <- file.path(survival_dir, "candidate_genes_tcga_coad_read_cox_ph_assumption.csv")
log_out <- file.path(survival_dir, "candidate_genes_tcga_coad_read_survival_log.txt")
forest_univariate_pdf <- file.path(figure_dir, "FigureS_candidate_gene_OS_univariate_Cox_forestplot.pdf")
forest_univariate_png <- file.path(figure_dir, "FigureS_candidate_gene_OS_univariate_Cox_forestplot.png")
forest_stage_binary_pdf <- file.path(figure_dir, "FigureS_candidate_gene_OS_stage_binary_adjusted_Cox_forestplot.pdf")
forest_stage_binary_png <- file.path(figure_dir, "FigureS_candidate_gene_OS_stage_binary_adjusted_Cox_forestplot.png")

target_genes <- c("SERPINE1", "PLAUR", "SDC4", "CCN1")
log_lines <- character()
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", sprintf(...))
  message(msg)
  log_lines <<- c(log_lines, msg)
}
write_log <- function() writeLines(log_lines, log_out)
stop_with_log <- function(...) {
  msg <- sprintf(...)
  log_msg("ERROR: %s", msg)
  write_log()
  stop(msg, call. = FALSE)
}

normalize_patient_id <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x == "" | is.na(x)] <- NA_character_
  ifelse(is.na(x), NA_character_, substr(x, 1, 12))
}
normalize_colname <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[ .()\\-]+", "_", x)
  x <- gsub("[^a-z0-9_]+", "", x)
  x <- gsub("_+", "_", x)
  gsub("^_+|_+$", "", x)
}
find_col <- function(dt, candidates) {
  actual <- normalize_colname(names(dt))
  for (candidate in normalize_colname(candidates)) {
    hit <- which(actual == candidate)
    if (length(hit)) return(names(dt)[hit[1]])
  }
  NA_character_
}
numeric_clean <- function(x) {
  suppressWarnings(as.numeric(gsub("[^0-9eE+.-]", "", as.character(x))))
}
parse_os_event <- function(x) {
  s <- tolower(trimws(as.character(x)))
  out <- rep(NA_integer_, length(s))
  out[grepl("deceased|dead", s) | grepl("^1([:_ ].*)?$", s)] <- 1L
  out[grepl("living|alive", s) | grepl("^0([:_ ].*)?$", s)] <- 0L
  out
}

load_annotation_targets <- function() {
  if (!file.exists(annotation_path)) {
    log_msg("Optional gene annotation not found; using expression row names/gene_name: %s", annotation_path)
    return(data.table(
      gene = target_genes,
      annotation_symbol = NA_character_,
      ensembl_clean = NA_character_
    ))
  }
  ann <- fread(annotation_path, data.table = TRUE, showProgress = FALSE)
  symbol_col <- find_col(ann, c("gene_name", "gene_symbol", "symbol"))
  id_col <- find_col(ann, c("gene_id_clean", "gene_id", "ensembl_gene_id"))
  if (is.na(symbol_col) || is.na(id_col)) stop_with_log("Annotation lacks gene symbol or Ensembl ID columns.")
  ann[, symbol_upper := toupper(trimws(as.character(get(symbol_col))))]
  ann[, ensembl_clean := sub("\\..*$", "", as.character(get(id_col)))]
  rbindlist(lapply(target_genes, function(gene_name) {
    aliases <- if (gene_name == "CCN1") c("CCN1", "CYR61") else gene_name
    hit <- ann[symbol_upper %in% aliases]
    if (!nrow(hit)) {
      return(data.table(gene = gene_name, annotation_symbol = NA_character_, ensembl_clean = NA_character_))
    }
    hit <- hit[match(aliases, symbol_upper, nomatch = length(aliases) + 1)][1]
    data.table(
      gene = gene_name,
      annotation_symbol = hit$symbol_upper,
      ensembl_clean = hit$ensembl_clean
    )
  }))
}

expression_to_matrix <- function(path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    obj <- readRDS(path)
    if (is.matrix(obj)) {
      return(list(matrix = obj, gene_symbols = rownames(obj)))
    }
    if (inherits(obj, "data.frame")) {
      dt <- as.data.table(obj, keep.rownames = "row_name")
      id_col <- find_col(dt, c("gene_id", "gene_id_clean", "row_name"))
      symbol_col <- find_col(dt, c("gene_name", "gene_symbol", "symbol"))
      if (is.na(id_col)) id_col <- "row_name"
      ids <- as.character(dt[[id_col]])
      gene_symbols <- if (is.na(symbol_col)) ids else as.character(dt[[symbol_col]])
      value_cols <- setdiff(names(dt), c(id_col, symbol_col, "row_name"))
      mat <- as.matrix(dt[, ..value_cols])
      storage.mode(mat) <- "numeric"
      rownames(mat) <- ids
      return(list(matrix = mat, gene_symbols = gene_symbols))
    }
    stop_with_log("Unsupported expression RDS class: %s", paste(class(obj), collapse = ", "))
  }
  dt <- fread(path, data.table = TRUE, showProgress = FALSE)
  id_col <- find_col(dt, c("gene_id", "gene_id_clean"))
  if (is.na(id_col)) id_col <- names(dt)[1]
  symbol_col <- find_col(dt, c("gene_name", "gene_symbol", "symbol"))
  ids <- as.character(dt[[id_col]])
  gene_symbols <- if (is.na(symbol_col)) ids else as.character(dt[[symbol_col]])
  value_cols <- setdiff(names(dt), c(id_col, symbol_col))
  mat <- as.matrix(dt[, ..value_cols])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- ids
  list(matrix = mat, gene_symbols = gene_symbols)
}

extract_candidate_expression <- function(path, targets) {
  expression_data <- expression_to_matrix(path)
  mat <- expression_data$matrix
  gene_symbols <- toupper(trimws(as.character(expression_data$gene_symbols)))
  ids_clean <- sub("\\..*$", "", rownames(mat))
  rows <- list()
  found <- character()
  missing <- character()
  for (gene_name in target_genes) {
    target <- targets[gene == gene_name]
    aliases <- if (gene_name == "CCN1") c("CCN1", "CYR61") else gene_name
    hit <- which(
      gene_symbols %in% aliases |
        toupper(rownames(mat)) %in% aliases |
        (!is.na(target$ensembl_clean[1]) & ids_clean == target$ensembl_clean[1])
    )
    if (!length(hit)) {
      missing <- c(missing, gene_name)
      next
    }
    found <- c(found, gene_name)
    values <- colMeans(mat[hit, , drop = FALSE], na.rm = TRUE)
    rows[[gene_name]] <- data.table(
      sample_id = colnames(mat),
      gene = gene_name,
      expression = as.numeric(values)
    )
  }
  if (!length(rows)) stop_with_log("None of the candidate genes were found in expression data.")
  out <- rbindlist(rows)
  out <- out[is.finite(expression)]
  out[, `:=`(
    patient_id = normalize_patient_id(sample_id),
    sample_type_code = substr(sample_id, 14, 15)
  )]
  attr(out, "genes_found") <- found
  attr(out, "genes_missing") <- missing
  out
}

filter_primary_first_sample <- function(expr) {
  expr <- expr[sample_type_code == "01"]
  if (file.exists(sample_map_path)) {
    map <- fread(sample_map_path, data.table = TRUE, showProgress = FALSE)
    map_sample <- find_col(map, c("sample_id", "tcga_barcode"))
    map_type <- find_col(map, "sample_type_code")
    if (!is.na(map_sample) && !is.na(map_type)) {
      primary_ids <- unique(as.character(map[
        suppressWarnings(as.integer(get(map_type))) == 1L,
        get(map_sample)
      ]))
      expr <- expr[sample_id %in% primary_ids]
    }
  }
  setorder(expr, gene, patient_id, sample_id)
  duplicate_n <- expr[, .N, by = .(gene, patient_id)][N > 1, .N]
  log_msg("Gene-patient combinations with multiple Primary Tumor samples: %d", duplicate_n)
  log_msg("Duplicate expression handling: first sample_id retained per gene and patient.")
  expr[, .SD[1], by = .(gene, patient_id)]
}

load_clinical <- function(path) {
  if (!file.exists(path)) stop_with_log("Clinical file not found: %s", path)
  dt <- fread(
    path, sep = "\t", data.table = TRUE, showProgress = FALSE,
    fill = TRUE, check.names = FALSE, encoding = "UTF-8"
  )
  log_msg("Clinical columns: %s", paste(names(dt), collapse = ", "))
  patient_col <- find_col(dt, c(
    "Patient ID", "patient_id", "PATIENT_ID", "patient", "case_id", "participant_id", "submitter_id"
  ))
  sample_col <- find_col(dt, c(
    "Sample ID", "sample_id", "SAMPLE_ID", "sample", "Tumor Sample Barcode", "tumor_sample_barcode"
  ))
  if (is.na(patient_col) && is.na(sample_col)) {
    stop_with_log("No patient/sample ID column. Available columns: %s", paste(names(dt), collapse = ", "))
  }
  dt[, patient_id := if (!is.na(patient_col)) {
    normalize_patient_id(get(patient_col))
  } else {
    normalize_patient_id(get(sample_col))
  }]
  time_col <- find_col(dt, c(
    "Overall Survival (Months)", "overall_survival_months", "OS_MONTHS", "os_months", "os_time"
  ))
  status_col <- find_col(dt, c(
    "Overall Survival Status", "overall_survival_status", "OS_STATUS", "os_status", "vital_status"
  ))
  if (is.na(time_col) || is.na(status_col)) {
    stop_with_log(
      "Required OS columns missing. time=%s; status=%s; available=%s",
      time_col, status_col, paste(names(dt), collapse = ", ")
    )
  }
  stage_col <- find_col(dt, c("Tumor Stage 2009", "stage", "stage_group", "ajcc_pathologic_stage"))
  sex_col <- find_col(dt, c("Sex", "sex", "gender"))
  log_msg("Patient ID column used: %s", ifelse(is.na(patient_col), "none", patient_col))
  log_msg("Sample ID column used: %s", ifelse(is.na(sample_col), "none", sample_col))
  log_msg("OS time/status columns used: %s / %s", time_col, status_col)
  log_msg("Stage/sex columns used: %s / %s", stage_col, sex_col)
  dt[, `:=`(
    OS_time_months = numeric_clean(get(time_col)),
    OS_event = parse_os_event(get(status_col)),
    stage_raw = if (is.na(stage_col)) NA_character_ else as.character(get(stage_col)),
    sex = if (is.na(sex_col)) NA_character_ else as.character(get(sex_col)),
    input_order = .I
  )]
  dt <- dt[!is.na(patient_id)]
  dup <- dt[, .N, by = patient_id][N > 1]
  log_msg("Duplicate clinical patients: %d; first input row retained.", nrow(dup))
  setorder(dt, patient_id, input_order)
  dt[, .SD[1], by = patient_id][, .(
    patient_id, OS_time_months, OS_event, stage_raw, sex
  )]
}

add_stage_groups <- function(dt) {
  x <- copy(dt)
  s <- toupper(trimws(as.character(x$stage_raw)))
  s[is.na(s)] <- ""
  x[, stage_group := fifelse(
    grepl("STAGE[ _-]*IV|(^|[^I])IV([ABC]|$)", s), "Stage IV",
    fifelse(
      grepl("STAGE[ _-]*III|(^|[^I])III([ABC]|$)", s), "Stage III",
      fifelse(
        grepl("STAGE[ _-]*II|(^|[^I])II([ABC]|$)", s), "Stage II",
        fifelse(grepl("STAGE[ _-]*I|(^|[^I])I([ABC]|$)", s), "Stage I", "Unknown")
      )
    )
  )]
  x[, stage_group := factor(stage_group, levels = c(
    "Stage I", "Stage II", "Stage III", "Stage IV", "Unknown"
  ))]
  x[, stage_binary := fifelse(
    as.character(stage_group) %in% c("Stage I", "Stage II"), "Early tumor",
    fifelse(
      as.character(stage_group) %in% c("Stage III", "Stage IV"), "Advanced tumor", NA_character_
    )
  )]
  x[, stage_binary := factor(stage_binary, levels = c("Early tumor", "Advanced tumor"))]
  x
}

median_survival <- function(x) {
  if (!nrow(x)) return(NA_real_)
  unname(summary(survfit(Surv(OS_time_months, OS_event) ~ 1, data = x))$table["median"])
}

fit_one_cox <- function(gene_name, model_name, formula, model_dt, min_n = 10L, min_events = 5L) {
  n_model <- nrow(model_dt)
  events <- sum(model_dt$OS_event)
  skip <- NULL
  if (n_model < min_n) skip <- sprintf("Skipped: n=%d is below %d.", n_model, min_n)
  if (is.null(skip) && events < min_events) skip <- sprintf("Skipped: events=%d is below %d.", events, min_events)
  if (is.null(skip) && uniqueN(model_dt$expression_group_median) < 2) {
    skip <- "Skipped: fewer than two expression groups."
  }
  if (!is.null(skip)) {
    return(list(
      result = data.table(
        gene = gene_name, model = model_name, variable = "gene_group_median",
        level_or_comparison = paste0(gene_name, "-high vs ", gene_name, "-low"),
        n = n_model, events = events, HR = NA_real_, lower_95CI = NA_real_,
        upper_95CI = NA_real_, p_value = NA_real_,
        fdr_within_model_for_gene_term = NA_real_,
        note = paste(skip, "Association analysis, not causation.")
      ),
      ph = data.table(
        gene = gene_name, model = model_name, variable = "GLOBAL",
        chisq = NA_real_, p_value = NA_real_, note = paste("cox.zph not run.", skip)
      ),
      fit = NULL, reason = skip
    ))
  }
  error_text <- NULL
  fit <- tryCatch(
    coxph(formula, data = model_dt, x = TRUE, model = TRUE),
    error = function(e) {
      error_text <<- conditionMessage(e)
      NULL
    }
  )
  if (is.null(fit)) {
    return(fit_one_cox(
      gene_name, model_name, formula, model_dt[0],
      min_n = 1L, min_events = 1L
    ))
  }
  sm <- summary(fit)
  coef_dt <- as.data.table(sm$coefficients, keep.rownames = "term")
  ci_dt <- as.data.table(sm$conf.int, keep.rownames = "term")
  result <- merge(
    coef_dt[, .(term, HR = `exp(coef)`, p_value = `Pr(>|z|)`)],
    ci_dt[, .(term, lower_95CI = `lower .95`, upper_95CI = `upper .95`)],
    by = "term", all = TRUE, sort = FALSE
  )
  result[, variable := fifelse(
    grepl("^expression_group_median", term), "gene_group_median",
    fifelse(grepl("^stage_group", term), "stage_group",
      fifelse(grepl("^stage_binary", term), "stage_binary",
        fifelse(grepl("^sex", term), "sex", term)
      )
    )
  )]
  result[, level_or_comparison := fifelse(
    variable == "gene_group_median",
    paste0(gene_name, "-high vs ", gene_name, "-low"),
    fifelse(
      variable == "stage_binary", "Advanced tumor vs Early tumor",
      fifelse(variable == "stage_group", paste0(sub("^stage_group", "", term), " vs Stage I"),
        fifelse(variable == "sex", paste0(sub("^sex", "", term), " vs reference sex"), term)
      )
    )
  )]
  result[, `:=`(
    gene = gene_name, model = model_name, n = n_model, events = events,
    fdr_within_model_for_gene_term = NA_real_,
    note = "Cox proportional hazards model; association analysis, not causation."
  )]
  result <- result[, .(
    gene, model, variable, level_or_comparison, n, events, HR,
    lower_95CI, upper_95CI, p_value, fdr_within_model_for_gene_term, note
  )]
  z <- tryCatch(cox.zph(fit), error = function(e) NULL)
  ph <- if (is.null(z)) {
    data.table(
      gene = gene_name, model = model_name, variable = "GLOBAL",
      chisq = NA_real_, p_value = NA_real_, note = "cox.zph failed."
    )
  } else {
    zdt <- as.data.table(z$table, keep.rownames = "variable")
    setnames(zdt, "p", "p_value")
    zdt[, .(
      gene = gene_name, model = model_name, variable, chisq, p_value,
      note = "cox.zph test; p < 0.05 suggests possible non-proportional hazards."
    )]
  }
  list(result = result, ph = ph, fit = fit, reason = error_text)
}

save_km <- function(gene_name, fit, dt, logrank_p, hr, lower, upper) {
  pdf_path <- file.path(
    figure_dir, sprintf("FigureS_%s_OS_KM_TCGA_COAD_READ_median_split.pdf", gene_name)
  )
  png_path <- file.path(
    figure_dir, sprintf("FigureS_%s_OS_KM_TCGA_COAD_READ_median_split.png", gene_name)
  )
  annotation <- sprintf(
    "Log-rank p=%s; HR=%.2f (95%% CI %.2f–%.2f)",
    format.pval(logrank_p, digits = 3, eps = 0.001), hr, lower, upper
  )
  g <- survminer::ggsurvplot(
    fit, data = dt, risk.table = TRUE, risk.table.height = 0.25,
    pval = annotation, conf.int = FALSE, palette = c("#4E79A7", "#D55E5E"),
    xlab = "Overall survival (months)", ylab = "Overall survival probability",
    title = sprintf("%s expression and overall survival in TCGA-COAD/READ", gene_name),
    legend.title = sprintf("%s expression", gene_name),
    legend.labs = c(paste0(gene_name, "-low"), paste0(gene_name, "-high")),
    ggtheme = theme_minimal(base_size = 11),
    tables.theme = survminer::theme_cleantable()
  )
  g$plot <- g$plot + labs(subtitle = "Median split")
  pdf(pdf_path, width = 7.5, height = 7.2)
  print(g)
  dev.off()
  png(png_path, width = 2250, height = 2160, res = 300)
  print(g)
  dev.off()
}

plot_forest <- function(plot_dt, pdf_path, png_path, title, fdr_col = "fdr") {
  x <- copy(plot_dt[is.finite(HR) & is.finite(lower_95CI) & is.finite(upper_95CI)])
  if (!nrow(x)) {
    log_msg("Forest plot skipped; no finite estimates: %s", title)
    return(invisible(NULL))
  }
  x[, gene := factor(gene, levels = rev(target_genes))]
  x[, label := sprintf("FDR=%s", format.pval(get(fdr_col), digits = 2, eps = 0.001))]
  p <- ggplot(x, aes(x = HR, y = gene)) +
    geom_vline(xintercept = 1, linetype = 2, color = "grey45") +
    geom_errorbarh(aes(xmin = lower_95CI, xmax = upper_95CI), height = 0.18, linewidth = 0.5) +
    geom_point(size = 2.5, color = "#B54A4A") +
    geom_text(aes(label = label), hjust = -0.08, size = 3.1) +
    scale_x_log10() +
    labs(x = "Hazard ratio (log scale)", y = NULL, title = title) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.y = element_blank(), plot.background = element_rect(fill = "white", color = NA))
  ggsave(pdf_path, p, width = 7.2, height = 4.8, device = "pdf")
  ggsave(png_path, p, width = 7.2, height = 4.8, dpi = 300)
}

log_msg("Project directory: %s", project_dir)
expression_path <- expression_candidates[file.exists(expression_candidates)][1]
if (is.na(expression_path)) {
  stop_with_log(
    "No VST expression input found. Required (either file): %s",
    paste(expression_candidates, collapse = " OR ")
  )
}
if (!file.exists(clinical_path)) {
  stop_with_log("Clinical input not found. Required file: %s", clinical_path)
}
log_msg("Expression file used: %s", expression_path)
log_msg("Clinical file used: %s", clinical_path)

targets <- load_annotation_targets()
expr_all <- extract_candidate_expression(expression_path, targets)
genes_found <- attr(expr_all, "genes_found")
genes_missing <- attr(expr_all, "genes_missing")
log_msg("Number of expression samples: %d", uniqueN(expr_all$sample_id))
log_msg("Genes analyzed: %s", paste(target_genes, collapse = ", "))
log_msg("Genes found: %s", paste(genes_found, collapse = ", "))
log_msg("Genes not found: %s", ifelse(length(genes_missing), paste(genes_missing, collapse = ", "), "none"))
expr <- filter_primary_first_sample(expr_all)
clinical <- add_stage_groups(load_clinical(clinical_path))
log_msg("Number of clinical samples: %d", nrow(clinical))

matched <- merge(expr[, .(patient_id, sample_id, gene, expression)], clinical, by = "patient_id")
matched <- matched[
  is.finite(expression) & is.finite(OS_time_months) &
    OS_time_months >= 0 & OS_event %in% c(0L, 1L)
]
log_msg("Matched patient number: %d", uniqueN(matched$patient_id))
if (!nrow(matched)) stop_with_log("No expression-clinical matched patients with valid OS data.")

sample_rows <- list()
summary_rows <- list()
test_rows <- list()
multi_objects <- list()
for (gene_name in genes_found) {
  x <- copy(matched[gene == gene_name])
  cutoff <- median(x$expression)
  x[, expression_group_median := factor(
    fifelse(expression >= cutoff, paste0(gene_name, "-high"), paste0(gene_name, "-low")),
    levels = c(paste0(gene_name, "-low"), paste0(gene_name, "-high"))
  )]
  sample_rows[[gene_name]] <- x[, .(
    patient_id, sample_id, gene, expression,
    expression_group_median = as.character(expression_group_median),
    OS_time_months, OS_event, stage_raw,
    stage_group = as.character(stage_group),
    stage_binary = as.character(stage_binary), sex
  )]
  for (group_name in levels(x$expression_group_median)) {
    gx <- x[expression_group_median == group_name]
    summary_rows[[paste(gene_name, group_name)]] <- data.table(
      gene = gene_name, cutoff_method = "median split", group = group_name,
      n = nrow(gx), events = sum(gx$OS_event),
      median_expression = median(gx$expression),
      median_survival_months = median_survival(gx)
    )
  }
  km_fit <- survfit(Surv(OS_time_months, OS_event) ~ expression_group_median, data = x)
  lr <- survdiff(Surv(OS_time_months, OS_event) ~ expression_group_median, data = x)
  lr_p <- pchisq(lr$chisq, df = length(lr$n) - 1, lower.tail = FALSE)
  uni <- coxph(Surv(OS_time_months, OS_event) ~ expression_group_median, data = x)
  us <- summary(uni)
  hr <- unname(us$coefficients[1, "exp(coef)"])
  lo <- unname(us$conf.int[1, "lower .95"])
  hi <- unname(us$conf.int[1, "upper .95"])
  up <- unname(us$coefficients[1, "Pr(>|z|)"])
  high <- x[expression_group_median == paste0(gene_name, "-high")]
  low <- x[expression_group_median == paste0(gene_name, "-low")]
  test_rows[[gene_name]] <- data.table(
    gene = gene_name, endpoint = "OS", cutoff_method = "median split",
    n_total = nrow(x), n_high = nrow(high), n_low = nrow(low),
    events_high = sum(high$OS_event), events_low = sum(low$OS_event),
    median_cutoff = cutoff, logrank_p = lr_p, logrank_fdr = NA_real_,
    univariate_HR = hr, univariate_lower_95CI = lo, univariate_upper_95CI = hi,
    univariate_p = up, univariate_fdr = NA_real_,
    note = "Univariate association analysis, not causation."
  )
  save_km(gene_name, km_fit, x, lr_p, hr, lo, hi)

  stage_x <- x[as.character(stage_group) != "Unknown"]
  stage_x[, stage_group := droplevels(stage_group)]
  binary_x <- x[!is.na(stage_binary)]
  binary_x[, stage_binary := droplevels(stage_binary)]
  sex_x <- binary_x[!is.na(sex) & trimws(sex) != ""]
  sex_x[, sex := droplevels(factor(trimws(sex)))]
  models <- list(
    fit_one_cox(
      gene_name, "model_stage",
      Surv(OS_time_months, OS_event) ~ expression_group_median + stage_group,
      stage_x
    ),
    fit_one_cox(
      gene_name, "model_stage_binary",
      Surv(OS_time_months, OS_event) ~ expression_group_median + stage_binary,
      binary_x
    )
  )
  if (uniqueN(sex_x$sex) >= 2) {
    models[[3]] <- fit_one_cox(
      gene_name, "model_stage_binary_sex",
      Surv(OS_time_months, OS_event) ~ expression_group_median + stage_binary + sex,
      sex_x
    )
  } else {
    models[[3]] <- fit_one_cox(
      gene_name, "model_stage_binary_sex",
      Surv(OS_time_months, OS_event) ~ expression_group_median + stage_binary,
      sex_x[0], min_n = 1L, min_events = 1L
    )
  }
  multi_objects[[gene_name]] <- models
  log_msg(
    "%s: expression min/median/mean/max=%s/%s/%s/%s; cutoff=%s; n high/low=%d/%d; events high/low=%d/%d",
    gene_name, signif(min(x$expression), 5), signif(median(x$expression), 5),
    signif(mean(x$expression), 5), signif(max(x$expression), 5), signif(cutoff, 6),
    nrow(high), nrow(low), sum(high$OS_event), sum(low$OS_event)
  )
  log_msg(
    "%s: log-rank p=%s; univariate HR=%s (95%% CI %s-%s), p=%s",
    gene_name, signif(lr_p, 5), signif(hr, 5), signif(lo, 5), signif(hi, 5), signif(up, 5)
  )
}

sample_dt <- rbindlist(sample_rows, use.names = TRUE, fill = TRUE)
summary_dt <- rbindlist(summary_rows, use.names = TRUE, fill = TRUE)
tests_dt <- rbindlist(test_rows, use.names = TRUE, fill = TRUE)
tests_dt[, logrank_fdr := p.adjust(logrank_p, method = "BH")]
tests_dt[, univariate_fdr := p.adjust(univariate_p, method = "BH")]
all_models <- unlist(multi_objects, recursive = FALSE)
multi_dt <- rbindlist(lapply(all_models, `[[`, "result"), use.names = TRUE, fill = TRUE)
ph_dt <- rbindlist(lapply(all_models, `[[`, "ph"), use.names = TRUE, fill = TRUE)
multi_dt[
  variable == "gene_group_median",
  fdr_within_model_for_gene_term := p.adjust(p_value, method = "BH"),
  by = model
]

fwrite(sample_dt, sample_out)
fwrite(summary_dt, summary_out)
fwrite(tests_dt, tests_out)
fwrite(multi_dt, multivariable_out)
fwrite(ph_dt, ph_out)

plot_forest(
  tests_dt[, .(
    gene, HR = univariate_HR, lower_95CI = univariate_lower_95CI,
    upper_95CI = univariate_upper_95CI, fdr = univariate_fdr
  )],
  forest_univariate_pdf, forest_univariate_png,
  "Candidate gene expression and overall survival: univariate Cox"
)
stage_binary_gene <- multi_dt[
  model == "model_stage_binary" & variable == "gene_group_median",
  .(
    gene, HR, lower_95CI, upper_95CI,
    fdr = fdr_within_model_for_gene_term
  )
]
plot_forest(
  stage_binary_gene,
  forest_stage_binary_pdf, forest_stage_binary_png,
  "Candidate gene expression and overall survival: stage-binary adjusted Cox"
)

for (gene_name in genes_found) {
  gene_multi <- multi_dt[gene == gene_name & variable == "gene_group_median"]
  for (i in seq_len(nrow(gene_multi))) {
    row <- gene_multi[i]
    log_msg(
      "%s %s: adjusted HR=%s (95%% CI %s-%s), p=%s, FDR=%s | %s",
      gene_name, row$model, signif(row$HR, 5), signif(row$lower_95CI, 5),
      signif(row$upper_95CI, 5), signif(row$p_value, 5),
      signif(row$fdr_within_model_for_gene_term, 5), row$note
    )
  }
  gene_ph <- ph_dt[gene == gene_name]
  log_msg(
    "%s PH assumption p-values: %s",
    gene_name,
    paste(sprintf("%s/%s=%s", gene_ph$model, gene_ph$variable, signif(gene_ph$p_value, 5)), collapse = "; ")
  )
}
log_msg("This is an association analysis, not causation.")
log_msg("Saved: %s", sample_out)
log_msg("Saved: %s", summary_out)
log_msg("Saved: %s", tests_out)
log_msg("Saved: %s", multivariable_out)
log_msg("Saved: %s", ph_out)
log_msg("Saved forest plots: %s; %s", forest_univariate_pdf, forest_stage_binary_pdf)
write_log()
