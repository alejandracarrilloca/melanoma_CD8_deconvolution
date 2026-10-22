#!/usr/bin/env Rscript

# ==============================================================================
# Step 27: Biological adjudication of the SPN-associated tumor-CD8 state
# ==============================================================================
#
# Questions:
#   1. Which genes and programs distinguish the Step 25 dysfunctional and
#      memory/transitional base states when patients are the replicates?
#   2. Is the dysfunctional state exhaustion-like, or better explained by
#      activation, cytotoxicity, proliferation, interferon, stress, or lineage
#      contamination?
#   3. After accounting for base state, sequencing depth, and overlays, which
#      programs remain associated with SPN detection within patients?
#
# State labels and eligibility were frozen before this analysis. SPN is not
# used to select cells, construct states, or define any gene program.
# ==============================================================================

options(stringsAsFactors = FALSE)

required_packages <- c(
  "data.table", "Matrix", "SingleCellExperiment", "SummarizedExperiment",
  "edgeR", "ggplot2"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ------------------------------------------------------------------------------
# Paths and constants
# ------------------------------------------------------------------------------

project_root <- normalizePath(getwd(), mustWork = TRUE)
analysis_dir <- file.path(
  project_root, "data", "analysis_ready", "scRNA", "GSE123139"
)
results_dir <- file.path(project_root, "results", "27")
figures_dir <- file.path(results_dir, "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

input_sce_file <- file.path(
  analysis_dir, "26_GSE123139_tumor_cd8_spn_annotated_sce.rds"
)
output_sce_file <- file.path(
  analysis_dir, "27_GSE123139_tumor_cd8_biological_programs_sce.rds"
)

random_seed <- 123139L
bootstrap_iterations <- 5000L
minimum_cells_per_patient_state <- 10L
minimum_cells_per_patient_model <- 100L
minimum_spn_group_cells <- 10L
use_robust_edger <- requireNamespace("statmod", quietly = TRUE)
base_states <- c("CD8_memory_transitional", "CD8_dysfunctional")

state_labels <- c(
  CD8_memory_transitional = "Memory/transitional",
  CD8_dysfunctional = "Dysfunctional"
)
state_colors <- c(
  CD8_memory_transitional = "#3B82B8",
  CD8_dysfunctional = "#7651B9"
)

theme_project <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.20)),
      plot.subtitle = element_text(color = "#5A5A5A"),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )
}

save_png <- function(filename, plot, width = 10, height = 8) {
  ggsave(
    file.path(figures_dir, filename), plot = plot, width = width,
    height = height, units = "in", dpi = 180, bg = "white"
  )
}

as_flag <- function(x) !is.na(x) & as.logical(x)

safe_wilcox_p <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3L || all(x == 0)) return(NA_real_)
  tryCatch(
    stats::wilcox.test(x, mu = 0, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
}

bootstrap_median_ci <- function(x, seed_offset = 0L) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(lower = NA_real_, upper = NA_real_))
  if (length(x) == 1L) return(c(lower = x, upper = x))
  set.seed(random_seed + seed_offset)
  boot <- replicate(
    bootstrap_iterations,
    stats::median(sample(x, length(x), replace = TRUE))
  )
  setNames(
    as.numeric(stats::quantile(boot, c(0.025, 0.975), names = FALSE)),
    c("lower", "upper")
  )
}

summarize_differences <- function(x, label, seed_offset = 0L) {
  x <- x[is.finite(x)]
  ci <- bootstrap_median_ci(x, seed_offset)
  data.table(
    comparison = label,
    n_patients = as.integer(length(x)),
    median_difference = if (length(x)) as.numeric(stats::median(x)) else NA_real_,
    bootstrap_ci_lower = unname(ci[[1L]]),
    bootstrap_ci_upper = unname(ci[[2L]]),
    percent_positive = if (length(x)) 100 * mean(x > 0) else NA_real_,
    wilcoxon_signed_rank_p = safe_wilcox_p(x)
  )
}

# Programs are deliberately specified before inspecting Step 27 results.
programs <- list(
  Exhaustion_core = c(
    "PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "TOX", "ENTPD1",
    "CXCL13", "LAYN", "BATF"
  ),
  Terminal_exhaustion = c(
    "HAVCR2", "ENTPD1", "CXCL13", "LAYN", "CD101", "TOX", "PRDM1",
    "EOMES", "TIGIT"
  ),
  Progenitor_dysfunction = c(
    "TCF7", "SLAMF6", "LEF1", "IL7R", "CCR7", "SELL", "MAL", "LTB"
  ),
  TCR_activation = c(
    "CD69", "TNFRSF9", "ICOS", "IL2RA", "FOS", "JUN", "JUNB",
    "NR4A1", "NR4A2", "EGR1"
  ),
  Cytotoxic_effector = c(
    "NKG7", "CCL5", "PRF1", "GZMB", "GZMH", "GNLY", "CTSW", "FGFBP2"
  ),
  Memory_naive = c(
    "TCF7", "LEF1", "CCR7", "IL7R", "SELL", "MAL", "LTB", "NOSIP"
  ),
  Tissue_resident = c(
    "CD69", "ITGAE", "CXCR6", "ZNF683", "RUNX3", "RGS1", "FABP5"
  ),
  Interferon_response = c(
    "ISG15", "IFI6", "IFIT1", "IFIT2", "IFIT3", "MX1", "ISG20",
    "STAT1", "IRF7", "IFI44L"
  ),
  Proliferation = c(
    "MKI67", "TOP2A", "STMN1", "TUBA1B", "TYMS", "UBE2C", "HMGB2",
    "CENPF", "NUSAP1"
  ),
  Stress = c(
    "FOS", "JUN", "JUNB", "DUSP1", "HSPA1A", "HSPA1B", "DNAJB1",
    "ATF3", "IER2"
  ),
  NK_lineage = c(
    "KLRD1", "NCR1", "FCER1G", "TYROBP", "TRDC", "TRGC1", "TRGC2",
    "KLRF1"
  ),
  CD4_Treg = c(
    "FOXP3", "IL2RA", "IKZF2", "CCR8", "TNFRSF18", "CTLA4"
  )
)

marker_categories <- list(
  Exhaustion = c(
    "PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "TOX", "ENTPD1",
    "CXCL13", "LAYN", "CD101", "BATF", "PRDM1"
  ),
  Progenitor_memory = c(
    "TCF7", "SLAMF6", "LEF1", "IL7R", "CCR7", "SELL", "MAL"
  ),
  Activation = c("CD69", "TNFRSF9", "ICOS", "IL2RA", "NR4A1", "EGR1"),
  Cytotoxic = c("NKG7", "CCL5", "PRF1", "GZMB", "GZMH", "GNLY", "CTSW"),
  Residency = c("ITGAE", "CXCR6", "ZNF683", "RGS1"),
  Lineage_review = c("KLRD1", "NCR1", "FCER1G", "TYROBP", "FOXP3")
)

# ------------------------------------------------------------------------------
# Load and validate Step 26
# ------------------------------------------------------------------------------

if (!file.exists(input_sce_file)) stop("Missing Step 26 SCE: ", input_sce_file)
message("Loading completed Step 26 object")
sce <- readRDS(input_sce_file)

required_metadata <- c(
  "subject_id", "step25_retained_cd8", "step25_base_state",
  "step25_signature_training_eligible", "step25_interferon_overlay",
  "step25_cytotoxic_overlay", "step25_proliferation_overlay",
  "endogenous_n_counts", "endogenous_n_genes", "step26_spn_detected",
  "step26_spn_logcounts"
)
missing_metadata <- setdiff(
  required_metadata, colnames(SummarizedExperiment::colData(sce))
)
if (length(missing_metadata)) {
  stop("Step 26 metadata is missing: ", paste(missing_metadata, collapse = ", "))
}
if (!all(c("counts", "logcounts") %in%
         SummarizedExperiment::assayNames(sce))) {
  stop("Step 26 must contain counts and logcounts assays.")
}
if (anyDuplicated(colnames(sce))) stop("Duplicated cell IDs occur.")

counts <- SummarizedExperiment::assay(sce, "counts")
logcounts <- SummarizedExperiment::assay(sce, "logcounts")
counts <- tryCatch(methods::as(counts, "dgCMatrix"), error = function(e) counts)

cell <- as.data.table(SummarizedExperiment::colData(sce))
cell[, cell_id := colnames(sce)]
cell[, subject_id := as.character(subject_id)]
cell[, step25_base_state := as.character(step25_base_state)]
cell[, log1p_endogenous_umis := log1p(as.numeric(endogenous_n_counts))]
cell[, any_overlay :=
  as_flag(step25_interferon_overlay) |
  as_flag(step25_cytotoxic_overlay) |
  as_flag(step25_proliferation_overlay)]

strict_candidates <- c(
  "strict_direct_cd8", "step21_strict_direct_cd8", "direct_cd8_evidence"
)
strict_column <- strict_candidates[strict_candidates %in% names(cell)][1L]
if (is.na(strict_column)) {
  warning("No strict direct-CD8 column was found; strict sensitivity omitted.")
  cell[, step27_strict_cd8 := FALSE]
} else {
  cell[, step27_strict_cd8 := as_flag(get(strict_column))]
}

retained_index <- which(
  as_flag(cell$step25_retained_cd8) &
    cell$step25_base_state %in% base_states
)
if (!length(retained_index)) stop("No retained Step 25 CD8 cells were found.")

training_index <- retained_index[
  as_flag(cell$step25_signature_training_eligible[retained_index])
]
if (!length(training_index)) stop("No signature-training-eligible cells remain.")

# ------------------------------------------------------------------------------
# Patient-paired pseudobulk differential expression
# ------------------------------------------------------------------------------

message("Running patient-paired pseudobulk differential expression")

run_paired_pseudobulk <- function(cell_indices, analysis_name) {
  meta <- cell[cell_indices, .(
    cell_index = cell_indices,
    subject_id,
    state = step25_base_state
  )]
  group_counts <- meta[, .(n_cells = .N), by = .(subject_id, state)]
  eligible_subjects <- group_counts[
    n_cells >= minimum_cells_per_patient_state,
    .(n_states = uniqueN(state)), by = subject_id
  ][n_states == 2L, subject_id]
  meta <- meta[subject_id %in% eligible_subjects]
  if (uniqueN(meta$subject_id) < 3L) {
    return(list(de = data.table(), qc = data.table(
      analysis = analysis_name, n_cells = nrow(meta),
      n_patients = uniqueN(meta$subject_id), n_genes_tested = 0L
    )))
  }

  group_key <- paste(meta$subject_id, meta$state, sep = "||")
  group_levels <- unique(group_key)
  group_factor <- factor(group_key, levels = group_levels)
  aggregation <- Matrix::sparse.model.matrix(~ 0 + group_factor)
  colnames(aggregation) <- group_levels
  pb_counts <- counts[, meta$cell_index, drop = FALSE] %*% aggregation

  first_row <- match(group_levels, group_key)
  pb_meta <- meta[first_row, .(subject_id, state)]
  pb_meta[, state := factor(state, levels = base_states)]
  pb_meta[, subject_id := factor(subject_id)]
  rownames(pb_meta) <- group_levels
  colnames(pb_counts) <- group_levels

  design <- stats::model.matrix(~ subject_id + state, data = pb_meta)
  y <- edgeR::DGEList(counts = pb_counts)
  keep <- edgeR::filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y)
  y <- edgeR::estimateDisp(y, design, robust = use_robust_edger)
  fit <- edgeR::glmQLFit(y, design, robust = use_robust_edger)
  state_coef <- grep("^state", colnames(design))
  if (length(state_coef) != 1L) stop("Could not identify the state coefficient.")
  test <- edgeR::glmQLFTest(fit, coef = state_coef)
  result <- as.data.table(edgeR::topTags(test, n = Inf, sort.by = "none")$table,
                          keep.rownames = "gene")
  setnames(result, c("logFC", "logCPM", "F", "PValue", "FDR"),
           c("log2_fold_change", "mean_log_cpm", "quasi_likelihood_F",
             "p_value", "fdr"))
  result[, `:=`(
    analysis = analysis_name,
    direction = fifelse(
      log2_fold_change > 0, "Higher_in_dysfunctional",
      fifelse(log2_fold_change < 0, "Higher_in_memory_transitional", "No_change")
    )
  )]
  setcolorder(result, c(
    "analysis", "gene", "log2_fold_change", "mean_log_cpm",
    "quasi_likelihood_F", "p_value", "fdr", "direction"
  ))
  list(
    de = result,
    qc = data.table(
      analysis = analysis_name,
      n_cells = as.integer(nrow(meta)),
      n_patients = as.integer(uniqueN(meta$subject_id)),
      n_pseudobulk_samples = as.integer(ncol(pb_counts)),
      n_genes_tested = as.integer(nrow(result)),
      n_fdr_005 = as.integer(sum(result$fdr < 0.05, na.rm = TRUE))
    )
  )
}

primary_pb <- run_paired_pseudobulk(training_index, "Primary_training_eligible")
no_overlay_index <- training_index[!cell$any_overlay[training_index]]
no_overlay_pb <- run_paired_pseudobulk(no_overlay_index, "No_overlay_sensitivity")
strict_index <- training_index[cell$step27_strict_cd8[training_index]]
strict_pb <- run_paired_pseudobulk(strict_index, "Strict_direct_CD8_sensitivity")

de_all <- rbindlist(list(
  primary_pb$de, no_overlay_pb$de, strict_pb$de
), fill = TRUE)
pseudobulk_qc <- rbindlist(list(
  primary_pb$qc, no_overlay_pb$qc, strict_pb$qc
), fill = TRUE)
if (!nrow(primary_pb$de)) stop("Primary pseudobulk analysis failed.")

de_concordance <- merge(
  primary_pb$de[, .(
    gene, primary_log2_fc = log2_fold_change, primary_fdr = fdr
  )],
  no_overlay_pb$de[, .(
    gene, no_overlay_log2_fc = log2_fold_change, no_overlay_fdr = fdr
  )], by = "gene", all = TRUE
)
de_concordance <- merge(
  de_concordance,
  strict_pb$de[, .(
    gene, strict_cd8_log2_fc = log2_fold_change, strict_cd8_fdr = fdr
  )], by = "gene", all = TRUE
)

# ------------------------------------------------------------------------------
# Independent program scores and patient-paired comparisons
# ------------------------------------------------------------------------------

message("Scoring independent biological programs")
program_presence <- rbindlist(lapply(names(programs), function(program_name) {
  requested <- unique(setdiff(programs[[program_name]], "SPN"))
  present <- intersect(requested, rownames(sce))
  data.table(
    program = program_name,
    n_genes_requested = as.integer(length(requested)),
    n_genes_present = as.integer(length(present)),
    genes_present = paste(present, collapse = ";")
  )
}))
if (any(program_presence$n_genes_present < 3L)) {
  stop("At least one program has fewer than three genes in the dataset.")
}

all_program_genes <- unique(unlist(lapply(programs, function(x) {
  intersect(setdiff(x, "SPN"), rownames(sce))
})))
program_expression <- as.matrix(
  logcounts[all_program_genes, retained_index, drop = FALSE]
)
gene_means <- rowMeans(program_expression)
gene_sds <- apply(program_expression, 1L, stats::sd)
gene_sds[!is.finite(gene_sds) | gene_sds == 0] <- 1
program_expression <- sweep(program_expression, 1L, gene_means, "-")
program_expression <- sweep(program_expression, 1L, gene_sds, "/")
program_expression[program_expression > 5] <- 5
program_expression[program_expression < -5] <- -5

for (program_name in names(programs)) {
  present <- intersect(setdiff(programs[[program_name]], "SPN"), all_program_genes)
  score <- colMeans(program_expression[present, , drop = FALSE])
  column_name <- paste0("step27_program_", program_name)
  sce[[column_name]] <- rep(NA_real_, ncol(sce))
  sce[[column_name]][retained_index] <- score
  cell[[column_name]] <- as.numeric(sce[[column_name]])
}

program_columns <- paste0("step27_program_", names(programs))
program_long <- melt(
  cell[retained_index],
  id.vars = c("subject_id", "step25_base_state", "cell_id"),
  measure.vars = program_columns,
  variable.name = "program", value.name = "program_score"
)
program_long[, program := sub("^step27_program_", "", program)]

patient_state_program <- program_long[, .(
  n_cells = as.integer(.N),
  mean_program_score = mean(program_score, na.rm = TRUE)
), by = .(subject_id, step25_base_state, program)][
  n_cells >= minimum_cells_per_patient_state
]

program_paired_tables <- list()
program_difference_rows <- list()
for (program_name in names(programs)) {
  temp <- patient_state_program[program == program_name]
  wide <- dcast(
    temp, subject_id ~ step25_base_state,
    value.var = "mean_program_score"
  )
  if (!all(base_states %in% names(wide))) next
  wide <- wide[complete.cases(wide[, ..base_states])]
  wide[, difference_dysfunctional_minus_memory :=
    CD8_dysfunctional - CD8_memory_transitional]
  wide[, program := program_name]
  program_paired_tables[[length(program_paired_tables) + 1L]] <- wide
  summary_row <- summarize_differences(
    wide$difference_dysfunctional_minus_memory,
    paste0("Dysfunctional minus memory/transitional: ", program_name),
    100L + length(program_difference_rows)
  )
  summary_row[, program := program_name]
  program_difference_rows[[length(program_difference_rows) + 1L]] <- summary_row
}
patient_program_differences <- rbindlist(program_paired_tables, fill = TRUE)
program_paired_summary <- rbindlist(program_difference_rows, fill = TRUE)
program_paired_summary[, fdr := p.adjust(wilcoxon_signed_rank_p, method = "BH")]

patient_balanced_program_state <- patient_state_program[, .(
  n_patient_groups = as.integer(.N),
  median_patient_mean_score = as.numeric(stats::median(mean_program_score)),
  q25_patient_mean_score = as.numeric(stats::quantile(mean_program_score, 0.25)),
  q75_patient_mean_score = as.numeric(stats::quantile(mean_program_score, 0.75))
), by = .(step25_base_state, program)]

# ------------------------------------------------------------------------------
# Within-patient SPN association with programs, adjusted for state and overlays
# ------------------------------------------------------------------------------

message("Testing which programs SPN tracks after state and depth adjustment")
model_dt <- cell[retained_index]
model_dt[, dysfunctional := as.integer(
  step25_base_state == "CD8_dysfunctional"
)]
model_dt[, spn_detected := as.integer(as_flag(step26_spn_detected))]
model_dt[, interferon_overlay := as.integer(as_flag(step25_interferon_overlay))]
model_dt[, cytotoxic_overlay := as.integer(as_flag(step25_cytotoxic_overlay))]
model_dt[, proliferation_overlay := as.integer(as_flag(step25_proliferation_overlay))]

fit_program_spn_model <- function(patient_table, program_column) {
  if (nrow(patient_table) < minimum_cells_per_patient_model ||
      sum(patient_table$spn_detected == 1L) < minimum_spn_group_cells ||
      sum(patient_table$spn_detected == 0L) < minimum_spn_group_cells ||
      uniqueN(patient_table$dysfunctional) < 2L) {
    return(data.table())
  }
  model_formula <- stats::as.formula(paste0(
    program_column,
    " ~ spn_detected + dysfunctional + log1p_endogenous_umis + ",
    "interferon_overlay + cytotoxic_overlay + proliferation_overlay"
  ))
  fit <- tryCatch(stats::lm(model_formula, data = patient_table),
                  error = function(e) NULL)
  if (is.null(fit) || !"spn_detected" %in% names(stats::coef(fit))) {
    return(data.table())
  }
  data.table(
    n_cells = as.integer(nrow(patient_table)),
    n_spn_detected = as.integer(sum(patient_table$spn_detected == 1L)),
    spn_adjusted_coefficient = unname(stats::coef(fit)[["spn_detected"]])
  )
}

spn_program_model_rows <- list()
for (program_name in names(programs)) {
  column_name <- paste0("step27_program_", program_name)
  model_result <- model_dt[, fit_program_spn_model(.SD, column_name),
                           by = subject_id]
  if (!nrow(model_result)) next
  model_result[, program := program_name]
  spn_program_model_rows[[length(spn_program_model_rows) + 1L]] <- model_result
}
spn_program_patient_models <- rbindlist(spn_program_model_rows, fill = TRUE)
spn_program_summary <- spn_program_patient_models[
  is.finite(spn_adjusted_coefficient), {
  ci <- bootstrap_median_ci(spn_adjusted_coefficient, 500L + .GRP)
  .(
    n_patients = as.integer(.N),
    median_spn_adjusted_coefficient = as.numeric(
      stats::median(spn_adjusted_coefficient)
    ),
    bootstrap_ci_lower = unname(ci[[1L]]),
    bootstrap_ci_upper = unname(ci[[2L]]),
    percent_positive = 100 * mean(spn_adjusted_coefficient > 0),
    wilcoxon_signed_rank_p = safe_wilcox_p(spn_adjusted_coefficient)
  )
}, by = program]
spn_program_summary[, fdr := p.adjust(wilcoxon_signed_rank_p, method = "BH")]

# ------------------------------------------------------------------------------
# Marker evidence and conservative state-label adjudication
# ------------------------------------------------------------------------------

message("Adjudicating the dysfunctional-state label")
marker_annotation <- rbindlist(lapply(names(marker_categories), function(category) {
  data.table(category = category, gene = marker_categories[[category]])
}))
marker_evidence <- merge(
  marker_annotation,
  primary_pb$de,
  by = "gene", all.x = TRUE
)
marker_evidence[, present_and_tested := !is.na(log2_fold_change)]
marker_evidence[, higher_in_dysfunctional_fdr_005 :=
  present_and_tested & log2_fold_change > 0 & fdr < 0.05]

primary_marker_candidates <- primary_pb$de[
  is.finite(fdr) & fdr < 0.05 & abs(log2_fold_change) >= 0.25
]
primary_marker_candidates[, technical_gene := grepl(
  "^(MT-|MTRNR|RPL|RPS|RP[0-9]+-|MALAT1$)", gene
)]
top_state_markers <- rbindlist(list(
  head(primary_marker_candidates[
    direction == "Higher_in_dysfunctional" & technical_gene == FALSE
  ][order(fdr, -abs(log2_fold_change))], 50L),
  head(primary_marker_candidates[
    direction == "Higher_in_memory_transitional" & technical_gene == FALSE
  ][order(fdr, -abs(log2_fold_change))], 50L)
), fill = TRUE)

program_lookup <- program_paired_summary[, .(
  program, program_median_difference = median_difference,
  program_fdr = fdr
)]
exhaustion_core <- program_lookup[program == "Exhaustion_core"]
terminal_exhaustion <- program_lookup[program == "Terminal_exhaustion"]
activation <- program_lookup[program == "TCR_activation"]
nk_lineage <- program_lookup[program == "NK_lineage"]
exhaustion_markers_up <- marker_evidence[
  category == "Exhaustion" & higher_in_dysfunctional_fdr_005 == TRUE,
  uniqueN(gene)
]

program_positive <- function(row) {
  nrow(row) == 1L && is.finite(row$program_fdr) &&
    row$program_fdr < 0.05 && row$program_median_difference > 0
}
core_pass <- program_positive(exhaustion_core)
terminal_pass <- program_positive(terminal_exhaustion)
activation_pass <- program_positive(activation)
nk_specific_markers_up <- marker_evidence[
  gene %in% c("KLRD1", "NCR1", "FCER1G", "TYROBP") &
    higher_in_dysfunctional_fdr_005 == TRUE,
  uniqueN(gene)
]
nk_warning <- program_positive(nk_lineage) && nk_specific_markers_up >= 2L

recommended_label <- if (core_pass && terminal_pass &&
    exhaustion_markers_up >= 3L && !nk_warning) {
  "Exhaustion-like dysfunctional CD8 state supported"
} else if ((core_pass || terminal_pass || exhaustion_markers_up >= 2L) &&
           !nk_warning) {
  "Dysfunctional CD8 state with partial exhaustion evidence"
} else {
  "Retain dysfunctional CD8 label; exhaustion not sufficiently supported"
}

state_label_adjudication <- data.table(
  criterion = c(
    "Exhaustion core program higher in dysfunctional",
    "Terminal exhaustion program higher in dysfunctional",
    "TCR activation program higher in dysfunctional",
    "Exhaustion markers higher at FDR < 0.05",
    "NK-lineage contamination warning",
    "Recommended label"
  ),
  value = c(
    as.character(core_pass), as.character(terminal_pass),
    as.character(activation_pass), as.character(exhaustion_markers_up),
    as.character(nk_warning), recommended_label
  )
)

# ------------------------------------------------------------------------------
# Outputs and figures
# ------------------------------------------------------------------------------

message("Writing Step 27 tables")
fwrite(pseudobulk_qc, file.path(results_dir, "27_pseudobulk_qc.csv"))
fwrite(de_all, file.path(results_dir, "27_paired_pseudobulk_de_all.csv"))
fwrite(primary_pb$de, file.path(results_dir, "27_primary_paired_pseudobulk_de.csv"))
fwrite(top_state_markers, file.path(results_dir, "27_top_state_markers.csv"))
fwrite(de_concordance, file.path(results_dir, "27_pseudobulk_de_concordance.csv"))
fwrite(program_presence, file.path(results_dir, "27_program_gene_coverage.csv"))
fwrite(patient_state_program, file.path(results_dir, "27_patient_state_program_scores.csv"))
fwrite(patient_program_differences, file.path(results_dir, "27_patient_program_differences.csv"))
fwrite(program_paired_summary, file.path(results_dir, "27_program_paired_summary.csv"))
fwrite(patient_balanced_program_state,
       file.path(results_dir, "27_patient_balanced_program_state.csv"))
fwrite(spn_program_patient_models,
       file.path(results_dir, "27_spn_program_patient_models.csv"))
fwrite(spn_program_summary, file.path(results_dir, "27_spn_program_summary.csv"))
fwrite(marker_evidence, file.path(results_dir, "27_key_marker_evidence.csv"))
fwrite(state_label_adjudication,
       file.path(results_dir, "27_state_label_adjudication.csv"))

dataset_qc <- data.table(
  metric = c(
    "input_cells", "retained_cd8_cells", "training_eligible_cells",
    "primary_pseudobulk_patients", "primary_genes_tested",
    "primary_genes_fdr_005", "programs_tested", "spn_programs_tested",
    "state_label", "state_accounting_validation", "signature_exported",
    "output_sce_size_mb"
  ),
  value = as.character(c(
    ncol(sce), length(retained_index), length(training_index),
    primary_pb$qc$n_patients, nrow(primary_pb$de),
    sum(primary_pb$de$fdr < 0.05, na.rm = TRUE), length(programs),
    nrow(spn_program_summary), recommended_label,
    as.integer(length(retained_index) == sum(as_flag(cell$step25_retained_cd8))),
    0L, NA_real_
  ))
)
fwrite(dataset_qc, file.path(results_dir, "27_dataset_qc.csv"))

message("Rendering Step 27 figures")
heatmap_dt <- copy(patient_balanced_program_state)
heatmap_dt[, state_label := factor(
  state_labels[step25_base_state], levels = unname(state_labels[base_states])
)]
heatmap_dt[, program := factor(program, levels = rev(names(programs)))]
p_heatmap <- ggplot(
  heatmap_dt,
  aes(state_label, program, fill = median_patient_mean_score)
) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#4C78A8", mid = "white", high = "#7651B9",
                       midpoint = 0) +
  labs(
    title = "Patient-balanced programs in the two tumor-CD8 base states",
    subtitle = "States were frozen before Step 27 and SPN is absent from every program",
    x = NULL, y = NULL, fill = "Median patient\nmean z-score"
  ) +
  theme_project()
save_png("27A_patient_balanced_program_heatmap.png", p_heatmap, 10, 8)

forest_dt <- copy(program_paired_summary)
forest_dt[, significant := is.finite(fdr) & fdr < 0.05]
forest_dt[, program := factor(program, levels = program[
  order(median_difference)
])]
p_forest <- ggplot(
  forest_dt,
  aes(median_difference, program, color = significant)
) +
  geom_vline(xintercept = 0, color = "#888888") +
  geom_errorbarh(aes(
    xmin = bootstrap_ci_lower, xmax = bootstrap_ci_upper
  ), height = 0.18) +
  geom_point(size = 3) +
  scale_color_manual(values = c(`TRUE` = "#7651B9", `FALSE` = "#999999")) +
  labs(
    title = "Patient-paired biological differences between CD8 states",
    subtitle = "Positive values indicate higher program activity in dysfunctional cells",
    x = "Median patient difference", y = NULL, color = "FDR < 0.05"
  ) +
  theme_project()
save_png("27B_patient_paired_program_differences.png", p_forest, 10, 8)

volcano_dt <- copy(primary_pb$de)
volcano_dt[, minus_log10_fdr := -log10(pmax(fdr, 1e-300))]
volcano_dt[, significant := fdr < 0.05 & abs(log2_fold_change) >= 0.25]
label_genes <- unique(unlist(marker_categories))
volcano_labels <- volcano_dt[gene %in% label_genes & significant == TRUE]
volcano_labels <- head(volcano_labels[order(fdr)], 20L)
p_volcano <- ggplot(
  volcano_dt,
  aes(log2_fold_change, minus_log10_fdr, color = significant)
) +
  geom_point(size = 0.8, alpha = 0.65) +
  geom_vline(xintercept = 0, color = "#888888") +
  scale_color_manual(values = c(`TRUE` = "#7651B9", `FALSE` = "#C8C8C8")) +
  geom_text(
    data = volcano_labels, aes(label = gene), size = 3,
    check_overlap = TRUE, vjust = -0.5, color = "black"
  ) +
  labs(
    title = "Patient-paired pseudobulk differential expression",
    subtitle = "Positive log2 fold changes indicate higher expression in dysfunctional cells",
    x = "Dysfunctional minus memory/transitional log2 fold change",
    y = expression(-log[10](FDR)), color = "FDR < 0.05\nand |log2FC| >= 0.25"
  ) +
  theme_project()
save_png("27C_primary_pseudobulk_volcano.png", p_volcano, 10, 8)

marker_plot <- marker_evidence[present_and_tested == TRUE]
marker_plot[, gene := factor(gene, levels = rev(unique(gene)))]
p_markers <- ggplot(
  marker_plot,
  aes(log2_fold_change, gene, color = -log10(pmax(fdr, 1e-300)))
) +
  geom_vline(xintercept = 0, color = "#888888") +
  geom_point(size = 3) +
  facet_wrap(~category, scales = "free_y", ncol = 2) +
  scale_color_viridis_c(option = "magma") +
  labs(
    title = "Independent marker evidence for the dysfunctional CD8 state",
    subtitle = "Positive values indicate higher expression in dysfunctional cells",
    x = "Paired-pseudobulk log2 fold change", y = NULL,
    color = expression(-log[10](FDR))
  ) +
  theme_project()
save_png("27D_key_marker_evidence.png", p_markers, 12, 10)

spn_forest <- copy(spn_program_summary)
spn_forest[, significant := is.finite(fdr) & fdr < 0.05]
spn_forest[, program := factor(program, levels = program[
  order(median_spn_adjusted_coefficient)
])]
p_spn <- ggplot(
  spn_forest,
  aes(median_spn_adjusted_coefficient, program, color = significant)
) +
  geom_vline(xintercept = 0, color = "#888888") +
  geom_errorbarh(aes(
    xmin = bootstrap_ci_lower, xmax = bootstrap_ci_upper
  ), height = 0.18) +
  geom_point(size = 3) +
  scale_color_manual(values = c(`TRUE` = "#D95F59", `FALSE` = "#999999")) +
  labs(
    title = "Programs associated with SPN within patients",
    subtitle = paste(
      "Adjusted for CD8 base state, UMI depth, interferon, cytotoxic,",
      "and proliferation overlays"
    ),
    x = "Median patient-specific SPN-detection coefficient",
    y = NULL, color = "FDR < 0.05"
  ) +
  theme_project()
save_png("27E_state_adjusted_spn_program_associations.png", p_spn, 11, 8)

exhaustion_patient <- patient_program_differences[
  program %in% c("Exhaustion_core", "Terminal_exhaustion")
]
if (nrow(exhaustion_patient)) {
  p_exhaustion <- ggplot(
    exhaustion_patient,
    aes(program, difference_dysfunctional_minus_memory, color = program)
  ) +
    geom_hline(yintercept = 0, color = "#888888") +
    geom_boxplot(outlier.shape = NA, width = 0.5, color = "#333333") +
    geom_jitter(width = 0.10, size = 2.3, alpha = 0.85) +
    guides(color = "none") +
    labs(
      title = "Exhaustion-program evidence across individual patients",
      subtitle = "Positive values indicate higher activity in dysfunctional cells",
      x = NULL, y = "Patient-paired program difference"
    ) +
    theme_project()
  save_png("27F_patient_exhaustion_program_differences.png", p_exhaustion, 9, 7)
}

message("Saving the Step 27 program-annotated object")
saveRDS(sce, output_sce_file, compress = FALSE)
dataset_qc[metric == "output_sce_size_mb", value := as.character(
  file.info(output_sce_file)$size / 1024 ^ 2
)]
fwrite(dataset_qc, file.path(results_dir, "27_dataset_qc.csv"))

writeLines(
  c(
    "Step 27 completed successfully.",
    paste("Retained CD8 cells:", length(retained_index)),
    paste("Training-eligible cells:", length(training_index)),
    paste("Primary paired-pseudobulk patients:", primary_pb$qc$n_patients),
    paste("Primary genes at FDR < 0.05:",
          sum(primary_pb$de$fdr < 0.05, na.rm = TRUE)),
    paste("Recommended state label:", recommended_label),
    "SPN was not used to define cells, states, overlays, or gene programs.",
    "No deconvolution signature was exported."
  ),
  file.path(results_dir, "27_completion.txt")
)

message("Step 27 complete: CD8 dysfunction and SPN-associated biology adjudicated")