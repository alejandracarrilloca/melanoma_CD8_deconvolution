#!/usr/bin/env Rscript

# ==============================================================================
# Step 24: Strict lineage and state-robustness gate for GSE123139 tumor CD8
# ==============================================================================
#
# This step deliberately attempts to falsify the Step 23 labels before they are
# used to build a bulk-deconvolution signature. It tests lineage contamination,
# patient-balanced program support, strict-direct-CD8 sensitivity, marker
# reproducibility, centroid collinearity, and cytotoxic/proliferating lineage
# evidence. It does not export a signature matrix.
# ==============================================================================

options(stringsAsFactors = FALSE)

required_packages <- c(
  "data.table", "Matrix", "SingleCellExperiment", "SummarizedExperiment",
  "ggplot2"
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
# Paths and review thresholds
# ------------------------------------------------------------------------------

project_root <- normalizePath(getwd(), mustWork = TRUE)
analysis_dir <- file.path(
  project_root, "data", "analysis_ready", "scRNA", "GSE123139"
)
results_dir <- file.path(project_root, "results", "24")
figures_dir <- file.path(results_dir, "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

input_sce_file <- file.path(
  analysis_dir, "23_GSE123139_tumor_cd8_final_states_sce.rds"
)
output_sce_file <- file.path(
  analysis_dir, "24_GSE123139_tumor_cd8_robustness_audited_sce.rds"
)
step23_transfer_file <- file.path(
  project_root, "results", "23", "23_cross_patient_transfer_summary.csv"
)

minimum_cells_per_patient_state <- 10L
minimum_cells_per_patient_state_marker <- 5L
minimum_cells_in_patient_remainder <- 10L
minimum_supported_patients <- 10L
minimum_cross_patient_recall <- 60
minimum_t_supported_percent <- 50
maximum_potential_non_t_percent <- 10
minimum_marker_delta_correlation <- 0.50
minimum_strict_robust_markers <- 5L
random_seed <- 123139L
maximum_plot_cells <- 30000L

state_levels <- c(
  "CD8_memory_transitional", "CD8_dysfunctional",
  "CD8_cytotoxic", "CD8_proliferating"
)
state_colors <- c(
  CD8_memory_transitional = "#3B82B8",
  CD8_dysfunctional = "#7651B9",
  CD8_cytotoxic = "#D95F59",
  CD8_proliferating = "#2F9E73"
)

lineage_programs <- list(
  T_identity = c("CD3D", "CD3E", "CD3G", "TRAC", "CD247", "LCK"),
  CD8_identity = c("CD8A", "CD8B"),
  CD4_Treg_specific = c("CD4", "FOXP3", "IL2RA", "CTLA4", "IKZF2"),
  NK_specific = c(
    "KLRD1", "FCGR3A", "TYROBP", "FCER1G", "NCR1", "NCR3"
  ),
  Gamma_delta_T = c("TRDC", "TRGC1", "TRGC2", "CD3D", "CD3E"),
  Myeloid = c(
    "LST1", "CTSS", "AIF1", "LILRB1", "CSF1R", "S100A8", "S100A9"
  ),
  B_cell = c("MS4A1", "CD79A", "CD79B", "CD37", "CD74", "CD22", "CD19"),
  Melanoma = c("MLANA", "PMEL", "SOX10", "MITF", "TYR", "DCT"),
  Cytotoxic = c("NKG7", "GNLY", "PRF1", "GZMB", "GZMH", "FGFBP2", "CTSW"),
  Dysfunction = c("PDCD1", "HAVCR2", "LAG3", "TIGIT", "TOX", "CXCL13", "ENTPD1"),
  Memory_transitional = c("TCF7", "IL7R", "CCR7", "LTB", "MAL", "LEF1", "SELL", "SLAMF6"),
  Proliferation = c("MKI67", "TOP2A", "STMN1", "TUBA1B", "HMGB2"),
  Interferon_response = c("ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "STAT1")
)

theme_project <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.22)),
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

safe_wilcox <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3L || all(x == 0)) return(1)
  tryCatch(
    stats::wilcox.test(x, mu = 0, exact = FALSE)$p.value,
    error = function(e) 1
  )
}

# ------------------------------------------------------------------------------
# Load and validate
# ------------------------------------------------------------------------------

if (!file.exists(input_sce_file)) stop("Missing Step 23 SCE: ", input_sce_file)
if (!file.exists(step23_transfer_file)) {
  stop("Missing Step 23 transfer summary: ", step23_transfer_file)
}
message("Loading completed Step 23 object")
sce <- readRDS(input_sce_file)

required_metadata <- c(
  "subject_id", "specimen_id", "amplification_batch", "step20_direct_T",
  "step20_direct_CD8", "step21_strict_direct_cd8", "step21_cd8_cluster",
  "step23_cd8_state", "step23_interferon_high",
  "step23_signature_training_eligible", "tcr_detected",
  "tcr_productive_primary", "tcr_clonotype_key", "endogenous_n_counts",
  "endogenous_n_genes"
)
missing_metadata <- setdiff(
  required_metadata, colnames(SummarizedExperiment::colData(sce))
)
if (length(missing_metadata)) {
  stop("Step 23 metadata is missing: ", paste(missing_metadata, collapse = ", "))
}
if (!all(c("counts", "logcounts") %in%
         SummarizedExperiment::assayNames(sce))) {
  stop("Step 23 must contain counts and logcounts assays.")
}
if (!"UMAP_step21" %in% SingleCellExperiment::reducedDimNames(sce)) {
  stop("Step 21 UMAP coordinates are missing from the Step 23 object.")
}

state <- as.character(sce$step23_cd8_state)
subject <- as.character(sce$subject_id)
if (!setequal(unique(state), state_levels)) stop("Unexpected Step 23 states.")

counts <- SummarizedExperiment::assay(sce, "counts")
logcounts <- SummarizedExperiment::assay(sce, "logcounts")
counts <- tryCatch(methods::as(counts, "dgCMatrix"), error = function(e) counts)
logcounts <- tryCatch(
  methods::as(logcounts, "dgCMatrix"),
  error = function(e) stop("Could not coerce logcounts: ", conditionMessage(e))
)

# ------------------------------------------------------------------------------
# Independent lineage programs and conservative direct-evidence flags
# ------------------------------------------------------------------------------

message("Scoring independent lineage and state programs")
score_program <- function(genes) {
  present <- intersect(genes, rownames(sce))
  if (!length(present)) return(rep(NA_real_, ncol(sce)))
  x <- logcounts[present, , drop = FALSE]
  gene_mean <- Matrix::rowMeans(x)
  gene_variance <- pmax(Matrix::rowMeans(x ^ 2) - gene_mean ^ 2, 0)
  gene_sd <- sqrt(gene_variance)
  keep <- is.finite(gene_sd) & gene_sd > 1e-8
  x <- x[keep, , drop = FALSE]
  gene_mean <- gene_mean[keep]
  gene_sd <- gene_sd[keep]
  if (!length(gene_sd)) return(rep(0, ncol(sce)))
  scaled_sum <- Matrix::colSums(
    Matrix::Diagonal(x = 1 / gene_sd) %*% x
  ) - sum(gene_mean / gene_sd)
  as.numeric(scaled_sum / length(gene_sd))
}

direct_detected <- function(genes, minimum_genes = 2L) {
  present <- intersect(genes, rownames(sce))
  if (!length(present)) return(rep(FALSE, ncol(sce)))
  as.numeric(Matrix::colSums(counts[present, , drop = FALSE] > 0)) >=
    minimum_genes
}

program_scores <- vapply(
  lineage_programs, score_program, numeric(ncol(sce))
)
colnames(program_scores) <- paste0("step24_score_", names(lineage_programs))
for (column_name in colnames(program_scores)) {
  sce[[column_name]] <- program_scores[, column_name]
}

t_supported <- as.logical(sce$step20_direct_T) |
  as.logical(sce$tcr_detected) | as.logical(sce$step21_strict_direct_cd8)
strict_direct_cd8 <- as.logical(sce$step21_strict_direct_cd8)
nk_direct <- direct_detected(lineage_programs$NK_specific, 2L)
myeloid_direct <- direct_detected(lineage_programs$Myeloid, 2L)
b_direct <- direct_detected(lineage_programs$B_cell, 2L)
melanoma_direct <- direct_detected(lineage_programs$Melanoma, 2L)
cd4_treg_direct <- direct_detected(lineage_programs$CD4_Treg_specific, 2L)
potential_non_t <- (nk_direct | myeloid_direct | b_direct | melanoma_direct) &
  !t_supported & !strict_direct_cd8
potential_cd4_treg <- cd4_treg_direct & !strict_direct_cd8

sce$step24_t_supported <- t_supported
sce$step24_nk_specific_direct <- nk_direct
sce$step24_myeloid_direct <- myeloid_direct
sce$step24_b_cell_direct <- b_direct
sce$step24_melanoma_direct <- melanoma_direct
sce$step24_cd4_treg_direct <- cd4_treg_direct
sce$step24_potential_non_t <- potential_non_t
sce$step24_potential_cd4_treg <- potential_cd4_treg

cell_metadata <- as.data.table(as.data.frame(
  SummarizedExperiment::colData(sce)
))
cell_metadata[, cell_id := colnames(sce)]
setcolorder(cell_metadata, c("cell_id", setdiff(names(cell_metadata), "cell_id")))
if (anyDuplicated(names(cell_metadata))) {
  stop("Duplicated columns occur in Step 24 cell metadata.")
}

evidence_by_state <- cell_metadata[, .(
  n_cells = .N,
  n_subjects = uniqueN(subject_id),
  n_subjects_with_at_least_10_cells = sum(table(subject_id) >= 10L),
  n_t_supported = sum(step24_t_supported),
  percent_t_supported = 100 * mean(step24_t_supported),
  n_strict_direct_cd8 = sum(step21_strict_direct_cd8),
  percent_strict_direct_cd8 = 100 * mean(step21_strict_direct_cd8),
  n_tcr = sum(tcr_detected),
  percent_tcr = 100 * mean(tcr_detected),
  n_nk_specific_direct = sum(step24_nk_specific_direct),
  percent_nk_specific_direct = 100 * mean(step24_nk_specific_direct),
  n_myeloid_direct = sum(step24_myeloid_direct),
  n_b_cell_direct = sum(step24_b_cell_direct),
  n_melanoma_direct = sum(step24_melanoma_direct),
  n_cd4_treg_direct = sum(step24_cd4_treg_direct),
  percent_cd4_treg_direct = 100 * mean(step24_cd4_treg_direct),
  n_potential_non_t = sum(step24_potential_non_t),
  percent_potential_non_t = 100 * mean(step24_potential_non_t),
  n_interferon_high = sum(step23_interferon_high),
  n_signature_training_eligible = sum(step23_signature_training_eligible)
), by = step23_cd8_state][order(match(step23_cd8_state, state_levels))]
fwrite(evidence_by_state, file.path(results_dir, "24_lineage_evidence_by_state.csv"))

# Patient-balanced summaries prevent large patients from defining lineage calls.
score_columns <- colnames(program_scores)
patient_state_scores <- cell_metadata[, c(
  "subject_id", "step23_cd8_state", score_columns
), with = FALSE][,
  c(list(n_cells = .N), lapply(.SD, mean, na.rm = TRUE)),
  by = .(subject_id, step23_cd8_state), .SDcols = score_columns
]
patient_state_scores_long <- melt(
  patient_state_scores,
  id.vars = c("subject_id", "step23_cd8_state", "n_cells"),
  variable.name = "program", value.name = "patient_mean_z_score"
)
patient_state_scores_long[, program := sub("^step24_score_", "", program)]
patient_state_scores_long[, sufficient_cells :=
  n_cells >= minimum_cells_per_patient_state]

patient_balanced_scores <- patient_state_scores_long[
  sufficient_cells == TRUE,
  .(
    n_supported_patients = uniqueN(subject_id),
    median_patient_mean_z_score = as.numeric(median(patient_mean_z_score)),
    first_quartile = as.numeric(quantile(patient_mean_z_score, 0.25)),
    third_quartile = as.numeric(quantile(patient_mean_z_score, 0.75)),
    percent_patients_positive = 100 * mean(patient_mean_z_score > 0)
  ), by = .(step23_cd8_state, program)
][order(match(step23_cd8_state, state_levels), program)]
fwrite(patient_state_scores_long,
       file.path(results_dir, "24_patient_state_program_scores.csv"))
fwrite(patient_balanced_scores,
       file.path(results_dir, "24_patient_balanced_program_summary.csv"))

# ------------------------------------------------------------------------------
# Strict-direct-CD8 sensitivity of state programs
# ------------------------------------------------------------------------------

message("Testing state programs in the strict direct-CD8 subset")
state_program_names <- c(
  "Cytotoxic", "Dysfunction", "Memory_transitional", "Proliferation",
  "T_identity", "CD8_identity", "NK_specific", "CD4_Treg_specific"
)
state_program_columns <- paste0("step24_score_", state_program_names)

summarize_subset_programs <- function(use_cells, subset_name) {
  dt <- cell_metadata[use_cells, c(
    "subject_id", "step23_cd8_state", state_program_columns
  ), with = FALSE]
  subject_means <- dt[,
    c(list(n_cells = .N), lapply(.SD, mean, na.rm = TRUE)),
    by = .(subject_id, step23_cd8_state), .SDcols = state_program_columns
  ]
  long <- melt(
    subject_means,
    id.vars = c("subject_id", "step23_cd8_state", "n_cells"),
    variable.name = "program", value.name = "patient_mean_z_score"
  )
  long[, program := sub("^step24_score_", "", program)]
  long <- long[n_cells >= minimum_cells_per_patient_state_marker]
  long[, .(
    n_supported_patients = uniqueN(subject_id),
    median_patient_mean_z_score = as.numeric(median(patient_mean_z_score)),
    percent_patients_positive = 100 * mean(patient_mean_z_score > 0)
  ), by = .(step23_cd8_state, program)][, subset := subset_name][]
}

program_sensitivity <- rbindlist(list(
  summarize_subset_programs(rep(TRUE, ncol(sce)), "All selected CD8"),
  summarize_subset_programs(strict_direct_cd8, "Strict direct-CD8")
), use.names = TRUE)
setcolorder(program_sensitivity, c(
  "subset", "step23_cd8_state", "program", "n_supported_patients",
  "median_patient_mean_z_score", "percent_patients_positive"
))
fwrite(program_sensitivity,
       file.path(results_dir, "24_strict_cd8_program_sensitivity.csv"))

program_sensitivity_wide <- dcast(
  program_sensitivity,
  step23_cd8_state + program ~ subset,
  value.var = "median_patient_mean_z_score"
)
setnames(
  program_sensitivity_wide,
  c("All selected CD8", "Strict direct-CD8"),
  c("all_cd8_median", "strict_cd8_median")
)
program_sensitivity_wide[, strict_minus_all :=
  strict_cd8_median - all_cd8_median]
fwrite(program_sensitivity_wide,
       file.path(results_dir, "24_strict_cd8_program_comparison.csv"))

# ------------------------------------------------------------------------------
# Patient-balanced final-state marker replication
# ------------------------------------------------------------------------------

message("Testing final-state markers in full and strict-CD8 training subsets")
row_data <- SummarizedExperiment::rowData(sce)
if (!all(c("step21_highly_variable", "step21_excluded_from_hvg") %in%
         colnames(row_data))) {
  stop("Step 21 gene annotations are missing from rowData.")
}
program_genes <- unique(unlist(lineage_programs, use.names = FALSE))
candidate_genes <- unique(c(
  rownames(sce)[as.logical(row_data$step21_highly_variable)],
  intersect(program_genes, rownames(sce))
))
candidate_genes <- candidate_genes[
  candidate_genes %in% rownames(sce) &
    !as.logical(row_data$step21_excluded_from_hvg[
      match(candidate_genes, rownames(sce))
    ])
]
candidate_logcounts <- logcounts[candidate_genes, , drop = FALSE]

marker_test <- function(use_cells, subset_name) {
  results <- vector("list", length(state_levels))
  names(results) <- state_levels
  for (state_name in state_levels) {
    eligible_subjects <- unique(subject[use_cells & state == state_name])
    eligible_subjects <- eligible_subjects[vapply(
      eligible_subjects,
      function(subject_name) {
        sum(use_cells & subject == subject_name & state == state_name) >=
          minimum_cells_per_patient_state_marker &&
          sum(use_cells & subject == subject_name & state != state_name) >=
          minimum_cells_in_patient_remainder
      }, logical(1L)
    )]
    if (length(eligible_subjects) < 3L) {
      warning("Fewer than three patients support ", subset_name, " markers for ",
              state_name)
      next
    }
    delta <- vapply(
      eligible_subjects,
      function(subject_name) {
        in_state <- use_cells & subject == subject_name & state == state_name
        remainder <- use_cells & subject == subject_name & state != state_name
        as.numeric(
          Matrix::rowMeans(candidate_logcounts[, in_state, drop = FALSE]) -
            Matrix::rowMeans(candidate_logcounts[, remainder, drop = FALSE])
        )
      }, numeric(length(candidate_genes))
    )
    if (is.null(dim(delta))) delta <- matrix(delta, ncol = 1L)
    p_value <- apply(delta, 1L, safe_wilcox)
    results[[state_name]] <- data.table(
      subset = subset_name,
      step23_cd8_state = state_name,
      gene = candidate_genes,
      n_paired_patients = ncol(delta),
      mean_within_patient_delta = rowMeans(delta),
      median_within_patient_delta = apply(delta, 1L, median),
      percent_patients_positive = 100 * rowMeans(delta > 0),
      p_value = p_value,
      adjusted_p_value = p.adjust(p_value, method = "BH")
    )
  }
  rbindlist(results, use.names = TRUE, fill = TRUE)
}

training_eligible <- as.logical(sce$step23_signature_training_eligible)
marker_full <- marker_test(training_eligible, "Full training set")
marker_strict <- marker_test(
  training_eligible & strict_direct_cd8, "Strict direct-CD8 training set"
)
marker_table <- rbindlist(list(marker_full, marker_strict), use.names = TRUE)
marker_table[, robust_marker :=
  median_within_patient_delta > 0 & percent_patients_positive >= 60 &
    adjusted_p_value <= 0.10]
fwrite(marker_table,
       file.path(results_dir, "24_patient_balanced_state_markers.tsv.gz"),
       sep = "\t")

marker_wide <- merge(
  marker_full[, .(
    step23_cd8_state, gene,
    full_delta = median_within_patient_delta,
    full_robust = median_within_patient_delta > 0 &
      percent_patients_positive >= 60 & adjusted_p_value <= 0.10
  )],
  marker_strict[, .(
    step23_cd8_state, gene,
    strict_delta = median_within_patient_delta,
    strict_robust = median_within_patient_delta > 0 &
      percent_patients_positive >= 60 & adjusted_p_value <= 0.10
  )],
  by = c("step23_cd8_state", "gene"), all = TRUE
)
marker_robustness <- marker_wide[, .(
  n_genes_compared = sum(is.finite(full_delta) & is.finite(strict_delta)),
  delta_spearman = suppressWarnings(cor(
    full_delta, strict_delta, method = "spearman", use = "complete.obs"
  )),
  n_full_robust_markers = sum(full_robust %in% TRUE),
  n_strict_robust_markers = sum(strict_robust %in% TRUE),
  n_robust_in_both = sum(full_robust %in% TRUE & strict_robust %in% TRUE)
), by = step23_cd8_state][order(match(step23_cd8_state, state_levels))]
fwrite(marker_robustness,
       file.path(results_dir, "24_marker_robustness_summary.csv"))

# ------------------------------------------------------------------------------
# State centroid collinearity and special-state lineage evidence
# ------------------------------------------------------------------------------

message("Calculating state-centroid correlations")
centroid_genes <- candidate_genes
centroid_matrix <- vapply(
  state_levels,
  function(state_name) {
    Matrix::rowMeans(candidate_logcounts[
      , training_eligible & state == state_name, drop = FALSE
    ])
  }, numeric(length(centroid_genes))
)
rownames(centroid_matrix) <- centroid_genes
colnames(centroid_matrix) <- state_levels
centroid_correlation <- stats::cor(centroid_matrix, method = "spearman")
# Build the long-form table explicitly.  as.data.table(as.table(...)) relies on
# dimname dispatch that differs across data.table versions and can pass integer
# dimension codes to setnames(), producing "x is type 'integer'".
centroid_correlation_table <- data.table(
  state_1 = rep(
    as.character(colnames(centroid_correlation)),
    each = nrow(centroid_correlation)
  ),
  state_2 = rep(
    as.character(rownames(centroid_correlation)),
    times = ncol(centroid_correlation)
  ),
  spearman_rho = as.numeric(centroid_correlation)
)
stopifnot(
  nrow(centroid_correlation_table) == length(state_levels)^2L,
  is.character(centroid_correlation_table$state_1),
  is.character(centroid_correlation_table$state_2)
)
fwrite(centroid_correlation_table,
       file.path(results_dir, "24_state_centroid_correlations.csv"))
message("State-centroid correlation table complete")

special_state_evidence <- cell_metadata[
  as.character(step23_cd8_state) %in% c("CD8_cytotoxic", "CD8_proliferating"),
  .(
    n_cells = .N,
    n_t_supported = sum(step24_t_supported),
    percent_t_supported = 100 * mean(step24_t_supported),
    n_strict_direct_cd8 = sum(step21_strict_direct_cd8),
    percent_strict_direct_cd8 = 100 * mean(step21_strict_direct_cd8),
    n_tcr = sum(tcr_detected),
    percent_tcr = 100 * mean(tcr_detected),
    n_nk_specific_direct = sum(step24_nk_specific_direct),
    percent_nk_specific_direct = 100 * mean(step24_nk_specific_direct),
    n_potential_non_t = sum(step24_potential_non_t),
    percent_potential_non_t = 100 * mean(step24_potential_non_t)
  ), by = .(step23_cd8_state, subject_id)
][order(step23_cd8_state, subject_id)]
fwrite(special_state_evidence,
       file.path(results_dir, "24_cytotoxic_proliferating_lineage_by_patient.csv"))

# Clonotype sharing supplies lineage support for cycling cells with suppressed
# lineage genes.
base_clone_state <- unique(cell_metadata[
  step23_cd8_state != "CD8_proliferating" & tcr_productive_primary == TRUE &
    !is.na(tcr_clonotype_key) & nzchar(tcr_clonotype_key),
  .(tcr_clonotype_key, base_state = as.character(step23_cd8_state))
])
proliferating_clones <- unique(cell_metadata[
  step23_cd8_state == "CD8_proliferating" & tcr_productive_primary == TRUE &
    !is.na(tcr_clonotype_key) & nzchar(tcr_clonotype_key),
  .(tcr_clonotype_key)
])
proliferating_clone_links <- merge(
  proliferating_clones, base_clone_state, by = "tcr_clonotype_key", all.x = TRUE
)
fwrite(proliferating_clone_links,
       file.path(results_dir, "24_proliferating_clonotype_base_state_links.csv"))

# ------------------------------------------------------------------------------
# Transparent screening table
# ------------------------------------------------------------------------------

transfer_summary <- fread(step23_transfer_file)
setnames(transfer_summary, "reference_state", "step23_cd8_state")
screening <- merge(
  evidence_by_state[, .(
    step23_cd8_state, n_cells, n_subjects_with_at_least_10_cells,
    percent_t_supported, percent_potential_non_t
  )],
  marker_robustness, by = "step23_cd8_state", all.x = TRUE
)
screening <- merge(
  screening,
  transfer_summary[, .(
    step23_cd8_state, cross_patient_recall_percent = accuracy_percent
  )],
  by = "step23_cd8_state", all.x = TRUE
)
screening[, support_pass :=
  n_subjects_with_at_least_10_cells >= minimum_supported_patients]
screening[, t_lineage_screen_pass :=
  percent_t_supported >= minimum_t_supported_percent]
screening[, non_t_screen_pass :=
  percent_potential_non_t <= maximum_potential_non_t_percent]
screening[, strict_marker_screen_pass :=
  is.finite(delta_spearman) & delta_spearman >= minimum_marker_delta_correlation &
    n_strict_robust_markers >= minimum_strict_robust_markers]
screening[, transfer_screen_pass :=
  is.na(cross_patient_recall_percent) |
    cross_patient_recall_percent >= minimum_cross_patient_recall]
screening[, all_available_screens_pass :=
  support_pass & t_lineage_screen_pass & non_t_screen_pass &
    strict_marker_screen_pass & transfer_screen_pass]
screening[, recommended_role := fifelse(
  step23_cd8_state == "CD8_cytotoxic" & !transfer_screen_pass,
  "Do not export standalone; test as cytotoxic overlay or merge",
  fifelse(
    all_available_screens_pass,
    "Eligible for cross-dataset validation",
    "Hold for review before signature export"
  )
)]
screening[, state_order__ := match(as.character(step23_cd8_state), state_levels)]
setorder(screening, state_order__)
screening[, state_order__ := NULL]
fwrite(screening, file.path(results_dir, "24_state_robustness_screening.csv"))

# ------------------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------------------

p_lineage <- ggplot(
  patient_balanced_scores,
  aes(x = program, y = step23_cd8_state, fill = median_patient_mean_z_score)
) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "#3B6FB6", mid = "white", high = "#7651B9", midpoint = 0,
    name = "Median patient\nmean z-score"
  ) +
  labs(
    title = "Patient-balanced lineage and state programs",
    subtitle = "Positive non-CD8 lineage programs trigger review rather than automatic relabeling",
    x = NULL, y = NULL
  ) + theme_project() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
save_png("24A_patient_balanced_lineage_state_heatmap.png", p_lineage, 13, 7)

evidence_plot <- melt(
  evidence_by_state,
  id.vars = "step23_cd8_state",
  measure.vars = c(
    "percent_t_supported", "percent_strict_direct_cd8", "percent_tcr",
    "percent_nk_specific_direct", "percent_cd4_treg_direct",
    "percent_potential_non_t"
  ),
  variable.name = "evidence", value.name = "percent_cells"
)
evidence_plot[, evidence := factor(evidence, levels = unique(evidence))]
p_evidence <- ggplot(
  evidence_plot,
  aes(x = step23_cd8_state, y = percent_cells, fill = step23_cd8_state)
) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~evidence, scales = "free_y") +
  scale_fill_manual(values = state_colors) +
  labs(
    title = "Direct lineage evidence and contamination screens",
    subtitle = "TCR absence is not treated as evidence against T-cell identity",
    x = NULL, y = "% of state cells"
  ) + theme_project() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_png("24B_direct_lineage_evidence_by_state.png", p_evidence, 14, 8)

p_sensitivity <- ggplot(
  program_sensitivity,
  aes(x = subset, y = median_patient_mean_z_score,
      color = step23_cd8_state, group = step23_cd8_state)
) +
  geom_hline(yintercept = 0, color = "#999999", linewidth = 0.3) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.1) +
  facet_wrap(~program, scales = "free_y") +
  scale_color_manual(values = state_colors) +
  labs(
    title = "State programs under strict direct-CD8 sensitivity analysis",
    subtitle = "Large shifts indicate dependence on the broader cluster-aware inclusion rule",
    x = NULL, y = "Median patient mean z-score", color = "CD8 state"
  ) + theme_project() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
save_png("24C_strict_cd8_program_sensitivity.png", p_sensitivity, 13, 9)

centroid_plot <- copy(centroid_correlation_table)
centroid_plot[, state_1 := factor(state_1, levels = state_levels)]
centroid_plot[, state_2 := factor(state_2, levels = rev(state_levels))]
p_centroid <- ggplot(
  centroid_plot, aes(x = state_1, y = state_2, fill = spearman_rho)
) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", spearman_rho)), size = 3.4) +
  scale_fill_gradient2(
    low = "#3B6FB6", mid = "white", high = "#7651B9", midpoint = 0,
    limits = c(-1, 1), name = "Spearman rho"
  ) +
  labs(
    title = "Correlation among candidate signature-state centroids",
    subtitle = "High off-diagonal correlations warn against separate bulk deconvolution states",
    x = NULL, y = NULL
  ) + theme_project() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_png("24D_state_centroid_correlations.png", p_centroid, 9, 8)

marker_plot <- marker_wide[
  is.finite(full_delta) & is.finite(strict_delta) &
    (full_robust %in% TRUE | strict_robust %in% TRUE)
]
p_marker <- ggplot(
  marker_plot,
  aes(x = full_delta, y = strict_delta, color = step23_cd8_state)
) +
  geom_hline(yintercept = 0, color = "#AAAAAA", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "#AAAAAA", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "#777777") +
  geom_point(size = 0.75, alpha = 0.65) +
  facet_wrap(~step23_cd8_state, scales = "free") +
  scale_color_manual(values = state_colors, guide = "none") +
  labs(
    title = "Patient-balanced marker replication in strict CD8 cells",
    subtitle = "Only genes robust in either the full or strict analysis are displayed",
    x = "Full training-set median delta",
    y = "Strict direct-CD8 median delta"
  ) + theme_project()
save_png("24E_strict_cd8_marker_replication.png", p_marker, 12, 9)

umap <- SingleCellExperiment::reducedDim(sce, "UMAP_step21")
set.seed(random_seed)
plot_index <- if (ncol(sce) > maximum_plot_cells) {
  sort(sample.int(ncol(sce), maximum_plot_cells))
} else seq_len(ncol(sce))
contamination_plot <- data.table(
  UMAP1 = umap[plot_index, 1L], UMAP2 = umap[plot_index, 2L],
  category = fifelse(
    potential_non_t[plot_index], "Potential non-T evidence",
    fifelse(t_supported[plot_index], "T-supported", "No direct T evidence")
  )
)
contamination_plot[, category := factor(
  category,
  levels = c("No direct T evidence", "T-supported", "Potential non-T evidence")
)]
setorder(contamination_plot, category)
p_contamination <- ggplot(
  contamination_plot, aes(UMAP1, UMAP2, color = category)
) +
  geom_point(size = 0.42, alpha = 0.72) +
  scale_color_manual(values = c(
    "No direct T evidence" = "#DDDDDD", "T-supported" = "#3B82B8",
    "Potential non-T evidence" = "#D95F59"
  )) +
  labs(
    title = "Independent T-lineage support across the selected compartment",
    subtitle = "Potential non-T calls require multiple lineage markers and absence of T/CD8 support",
    color = NULL
  ) + theme_project()
save_png("24F_umap_lineage_support_screen.png", p_contamination, 10, 8)

# ------------------------------------------------------------------------------
# Save and report
# ------------------------------------------------------------------------------

if (sum(evidence_by_state$n_cells) != ncol(sce)) {
  stop("Step 24 state accounting failed.")
}
if (anyNA(sce$step24_t_supported) || anyNA(sce$step24_potential_non_t)) {
  stop("Step 24 evidence flags contain missing values.")
}

sce$step24_robustness_screen_complete <- TRUE
saveRDS(sce, output_sce_file, compress = FALSE)
output_size_mb <- unname(file.info(output_sce_file)$size / 1024^2)

cytotoxic_recall <- transfer_summary[
  step23_cd8_state == "CD8_cytotoxic", accuracy_percent
]
dataset_qc <- data.table(
  metric = c(
    "cells", "subjects", "states_tested", "lineage_programs",
    "candidate_marker_genes", "strict_direct_cd8_cells",
    "t_supported_cells", "potential_non_t_cells",
    "cytotoxic_cross_patient_recall_percent",
    "states_passing_all_available_screens", "state_accounting_validation",
    "robustness_screen_complete", "signature_exported", "output_sce_size_mb"
  ),
  value = as.character(c(
    ncol(sce), uniqueN(subject), length(state_levels), length(lineage_programs),
    length(candidate_genes), sum(strict_direct_cd8), sum(t_supported),
    sum(potential_non_t), cytotoxic_recall,
    sum(screening$all_available_screens_pass),
    sum(evidence_by_state$n_cells) == ncol(sce), TRUE, FALSE, output_size_mb
  ))
)
fwrite(dataset_qc, file.path(results_dir, "24_dataset_qc.csv"))

figure_manifest <- data.table(
  figure = c(
    "24A_patient_balanced_lineage_state_heatmap.png",
    "24B_direct_lineage_evidence_by_state.png",
    "24C_strict_cd8_program_sensitivity.png",
    "24D_state_centroid_correlations.png",
    "24E_strict_cd8_marker_replication.png",
    "24F_umap_lineage_support_screen.png"
  ),
  purpose = c(
    "Patient-balanced lineage and state program review",
    "Direct T/CD8 evidence and conservative contamination flags",
    "Sensitivity of state programs to strict CD8 inclusion",
    "Bulk-signature collinearity warning",
    "Replication of patient-balanced markers in strict CD8 cells",
    "Spatial distribution of lineage-support flags in the CD8 UMAP"
  )
)
fwrite(figure_manifest, file.path(results_dir, "24_figure_manifest.csv"))
writeLines(capture.output(sessionInfo()), file.path(results_dir, "24_sessionInfo.txt"))

message("Step 24 strict CD8 robustness audit complete.")
message("No signature was exported.")
message("Output SCE: ", output_sce_file)
print(dataset_qc)
print(evidence_by_state)
print(marker_robustness)
print(screening)