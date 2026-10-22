#!/usr/bin/env Rscript

# ==============================================================================
# Step 32: Consolidate SPN-continuum and clonotype evidence in GSE123139
# ==============================================================================
#
# Purpose
#   1. Estimate patient-specific, depth-adjusted SPN associations across the
#      continuous Step 31 CD8 axis, without treating cells as replicates.
#   2. Determine whether the Q4-Q5 enrichment survives adjustment for base
#      state, activation/cytotoxicity, proliferation, and interferon response.
#   3. Test clone expansion versus continuum position at the clone level.
#   4. Test whether expanded clones occupy more fine states than expected from
#      clone size alone using patient-stratified state-label permutations.
#   5. Estimate clone-centered SPN-continuum associations using only variation
#      among cells of the same productive clonotype.
#   6. Produce a compact evidence ledger for downstream interpretation.
#
# Interpretation
#   All cross-patient summaries use patients as independent units. The
#   continuum is cross-sectional and does not establish temporal direction.
#   The permutation test preserves clone sizes and patient-specific state
#   frequencies. SPN was not used to define the continuum or fine states.
# ==============================================================================

options(stringsAsFactors = FALSE)
options(error = function() {
  cat("\n--- STEP 32 TRACEBACK ---\n", file = stderr())
  traceback(20L)
  quit(save = "no", status = 1L, runLast = FALSE)
})

required_packages <- c(
  "data.table", "SingleCellExperiment", "SummarizedExperiment", "ggplot2"
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

project_root <- normalizePath(getwd(), mustWork = TRUE)
analysis_dir <- file.path(
  project_root, "data", "analysis_ready", "scRNA", "GSE123139"
)
results_dir <- file.path(project_root, "results", "32")
figures_dir <- file.path(results_dir, "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

input_sce_file <- file.path(
  analysis_dir, "31_GSE123139_cd8_continuum_clonal_states_sce.rds"
)
output_sce_file <- file.path(
  analysis_dir, "32_GSE123139_spn_continuum_clone_consolidated_sce.rds"
)

random_seed <- 123139L
bootstrap_iterations <- 5000L
permutation_iterations <- 1000L
minimum_patient_cells <- 100L
minimum_tcr_patient_cells <- 50L
minimum_informative_clones <- 10L
minimum_clone_size <- 2L

base_states <- c("CD8_memory_transitional", "CD8_dysfunctional")
program_names <- c(
  "TCR_activation", "Cytotoxic_effector", "Proliferation",
  "Interferon_response", "Exhaustion_core", "Terminal_exhaustion"
)

as_flag <- function(x) !is.na(x) & as.logical(x)

zscore_safe <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

safe_wilcox_p <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3L || all(x == 0)) return(NA_real_)
  tryCatch(
    suppressWarnings(stats::wilcox.test(x, mu = 0, exact = FALSE)$p.value),
    error = function(e) NA_real_
  )
}

bootstrap_median_ci <- function(x, seed_offset = 0L) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(lower = NA_real_, upper = NA_real_))
  if (length(x) == 1L) return(c(lower = x, upper = x))
  set.seed(random_seed + as.integer(seed_offset))
  boot <- replicate(
    bootstrap_iterations,
    stats::median(sample(x, length(x), replace = TRUE))
  )
  setNames(
    as.numeric(stats::quantile(boot, c(0.025, 0.975), names = FALSE)),
    c("lower", "upper")
  )
}

summarize_effect <- function(x, label, scale, seed_offset = 0L) {
  x <- x[is.finite(x)]
  ci <- bootstrap_median_ci(x, seed_offset)
  data.table(
    comparison = label,
    effect_scale = scale,
    n_patients = as.integer(length(x)),
    median_effect = if (length(x)) as.numeric(stats::median(x)) else NA_real_,
    bootstrap_ci_lower = unname(ci["lower"]),
    bootstrap_ci_upper = unname(ci["upper"]),
    percent_positive = if (length(x)) 100 * mean(x > 0) else NA_real_,
    wilcoxon_signed_rank_p = safe_wilcox_p(x)
  )
}

safe_coef <- function(fit, term) {
  output <- tryCatch(stats::coef(fit)[term], error = function(e) NA_real_)
  if (!length(output) || !is.finite(output)) NA_real_ else as.numeric(output)
}

safe_spearman <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L || uniqueN(x[keep]) < 2L || uniqueN(y[keep]) < 2L) {
    return(NA_real_)
  }
  suppressWarnings(stats::cor(x[keep], y[keep], method = "spearman"))
}

classification <- function(effect, fdr, expected_direction = "positive") {
  if (!is.finite(effect) || !is.finite(fdr)) return("Not estimable")
  direction_ok <- if (expected_direction == "positive") effect > 0 else effect < 0
  if (direction_ok && fdr < 0.05) return("Supported")
  if (direction_ok && fdr < 0.10) return("Suggestive")
  "Not supported"
}

message("Loading completed Step 31 object")
if (!file.exists(input_sce_file)) stop("Missing Step 31 SCE: ", input_sce_file)
sce <- readRDS(input_sce_file)

program_columns <- paste0("step27_program_", program_names)
required_metadata <- c(
  "subject_id", "step25_retained_cd8", "step25_base_state",
  "endogenous_n_counts", "step26_spn_detected", "step26_spn_logcounts",
  "step31_continuum_raw", "step31_continuum_percentile",
  "step31_fine_state", "tcr_productive_primary", "tcr_clonotype_key",
  program_columns
)
missing_metadata <- setdiff(
  required_metadata, colnames(SummarizedExperiment::colData(sce))
)
if (length(missing_metadata)) {
  stop("Required metadata are missing: ", paste(missing_metadata, collapse = ", "))
}

cell <- as.data.table(SummarizedExperiment::colData(sce))
cell[, cell_id := colnames(sce)]
cell[, subject_id := as.character(subject_id)]
cell[, step25_base_state := as.character(step25_base_state)]
cell[, step31_fine_state := as.character(step31_fine_state)]
cell[, tcr_clonotype_key := as.character(tcr_clonotype_key)]
cell[, retained_cd8 := as_flag(step25_retained_cd8)]
cell[, spn_detected := as.integer(as_flag(step26_spn_detected))]
cell[, spn_logcounts := as.numeric(step26_spn_logcounts)]
cell[, log1p_endogenous_umis := log1p(as.numeric(endogenous_n_counts))]

analysis_dt <- cell[
  retained_cd8 == TRUE & step25_base_state %chin% base_states &
    is.finite(step31_continuum_percentile) &
    is.finite(log1p_endogenous_umis)
]
if (!nrow(analysis_dt)) stop("No retained Step 31 CD8 cells were found.")
analysis_dt[, dysfunctional_indicator := as.integer(
  step25_base_state == "CD8_dysfunctional"
)]

# Patient-center all continuous predictors. The continuum percentile is already
# defined within patient, but centering makes coefficients and prediction grids
# more stable and keeps biological adjustment within patient.
continuous_columns <- c(
  "step31_continuum_percentile", "log1p_endogenous_umis", program_columns
)
for (column in continuous_columns) {
  analysis_dt[, (paste0(column, "_z")) := zscore_safe(get(column)),
              by = subject_id]
}

# ------------------------------------------------------------------------------
# Patient-specific adjusted linear and smooth continuum models
# ------------------------------------------------------------------------------

message("Fitting patient-specific adjusted continuum models")
fit_patient_continuum <- function(dt) {
  if (nrow(dt) < minimum_patient_cells || uniqueN(dt$spn_detected) < 2L) {
    return(data.table())
  }

  minimal_formula <- spn_detected ~ step31_continuum_percentile_z +
    log1p_endogenous_umis_z + dysfunctional_indicator
  adjusted_formula <- spn_detected ~ step31_continuum_percentile_z +
    log1p_endogenous_umis_z + dysfunctional_indicator +
    step27_program_TCR_activation_z + step27_program_Cytotoxic_effector_z +
    step27_program_Proliferation_z + step27_program_Interferon_response_z
  expression_formula <- spn_logcounts ~ step31_continuum_percentile_z +
    log1p_endogenous_umis_z + dysfunctional_indicator +
    step27_program_TCR_activation_z + step27_program_Cytotoxic_effector_z +
    step27_program_Proliferation_z + step27_program_Interferon_response_z

  minimal_fit <- tryCatch(
    stats::glm(minimal_formula, data = dt, family = stats::binomial()),
    error = function(e) NULL
  )
  adjusted_fit <- tryCatch(
    stats::glm(adjusted_formula, data = dt, family = stats::binomial()),
    error = function(e) NULL
  )
  expression_fit <- tryCatch(stats::lm(expression_formula, data = dt),
                             error = function(e) NULL)

  # Natural spline is used only for standardized marginal contrasts. For each
  # patient, predict every observed cell at the 10th and 90th continuum
  # percentiles while leaving its other covariates unchanged.
  smooth_formula <- spn_detected ~
    splines::ns(step31_continuum_percentile, df = 3L) +
    log1p_endogenous_umis_z + dysfunctional_indicator +
    step27_program_TCR_activation_z + step27_program_Cytotoxic_effector_z +
    step27_program_Proliferation_z + step27_program_Interferon_response_z
  smooth_fit <- tryCatch(
    stats::glm(smooth_formula, data = dt, family = stats::binomial()),
    error = function(e) NULL
  )
  smooth_difference <- NA_real_
  if (!is.null(smooth_fit)) {
    low <- copy(dt)
    high <- copy(dt)
    low[, step31_continuum_percentile := 0.10]
    high[, step31_continuum_percentile := 0.90]
    smooth_difference <- tryCatch(
      100 * mean(
        stats::predict(smooth_fit, newdata = high, type = "response") -
          stats::predict(smooth_fit, newdata = low, type = "response")
      ), error = function(e) NA_real_
    )
  }

  data.table(
    subject_id = dt$subject_id[1L],
    n_cells = as.integer(nrow(dt)),
    minimal_detection_log_odds = if (is.null(minimal_fit)) NA_real_ else
      safe_coef(minimal_fit, "step31_continuum_percentile_z"),
    fully_adjusted_detection_log_odds = if (is.null(adjusted_fit)) NA_real_ else
      safe_coef(adjusted_fit, "step31_continuum_percentile_z"),
    fully_adjusted_expression_coefficient = if (is.null(expression_fit)) NA_real_ else
      safe_coef(expression_fit, "step31_continuum_percentile_z"),
    adjusted_smooth_q90_minus_q10_detection_pp = smooth_difference
  )
}

patient_continuum_models <- rbindlist(
  lapply(split(analysis_dt, analysis_dt$subject_id), fit_patient_continuum),
  fill = TRUE
)

continuum_model_summary <- rbindlist(list(
  summarize_effect(
    patient_continuum_models$minimal_detection_log_odds,
    "SPN detection across continuum; depth and base-state adjusted",
    "Log-odds per patient-SD continuum", 100L
  ),
  summarize_effect(
    patient_continuum_models$fully_adjusted_detection_log_odds,
    "SPN detection across continuum; biological-process adjusted",
    "Log-odds per patient-SD continuum", 101L
  ),
  summarize_effect(
    patient_continuum_models$fully_adjusted_expression_coefficient,
    "SPN expression across continuum; biological-process adjusted",
    "SPN logcounts per patient-SD continuum", 102L
  ),
  summarize_effect(
    patient_continuum_models$adjusted_smooth_q90_minus_q10_detection_pp,
    "Adjusted smooth continuum contrast: 90th minus 10th percentile",
    "Standardized SPN-detection percentage points", 103L
  )
), fill = TRUE)
continuum_model_summary[, fdr := p.adjust(wilcoxon_signed_rank_p, method = "BH")]

# ------------------------------------------------------------------------------
# Productive clonotypes: expansion, clone-centered effects, and size-aware null
# ------------------------------------------------------------------------------

message("Constructing productive-clonotype summaries")
tcr_dt <- analysis_dt[
  as_flag(tcr_productive_primary) & !is.na(tcr_clonotype_key) &
    nzchar(tcr_clonotype_key)
]
tcr_dt[, clone_id := paste(subject_id, tcr_clonotype_key, sep = "::")]

clone_summary <- tcr_dt[, .(
  clone_size = as.integer(.N),
  median_continuum = as.numeric(stats::median(step31_continuum_raw)),
  mean_continuum = mean(step31_continuum_raw),
  n_fine_states = as.integer(uniqueN(step31_fine_state)),
  spn_detected_percent = 100 * mean(spn_detected),
  mean_spn_logcounts = mean(spn_logcounts)
), by = .(subject_id, clone_id)]
clone_summary[, log2_clone_size := log2(clone_size)]
clone_summary[, expansion_class := fifelse(
  clone_size >= 10L, "Large (>=10 cells)",
  fifelse(clone_size >= 3L, "Expanded (3-9 cells)", "Small (1-2 cells)")
)]

patient_clone_expansion <- clone_summary[, .(
  n_clones = as.integer(.N),
  spearman_clone_size_continuum = safe_spearman(log2_clone_size, median_continuum),
  spearman_clone_size_spn_detection = safe_spearman(
    log2_clone_size, spn_detected_percent
  )
), by = subject_id][n_clones >= minimum_informative_clones]

clone_expansion_summary <- rbindlist(list(
  summarize_effect(
    patient_clone_expansion$spearman_clone_size_continuum,
    "Clone size versus median CD8 continuum position",
    "Within-patient Spearman rho", 200L
  ),
  summarize_effect(
    patient_clone_expansion$spearman_clone_size_spn_detection,
    "Clone size versus SPN detection",
    "Within-patient Spearman rho", 201L
  )
), fill = TRUE)
clone_expansion_summary[, fdr := p.adjust(wilcoxon_signed_rank_p, method = "BH")]

# Clone-centered analysis removes all between-clone differences. Coefficients
# therefore ask whether cells positioned further along the continuum than other
# cells of the same clone also have higher SPN.
eligible_clone_ids <- clone_summary[clone_size >= minimum_clone_size, clone_id]
clone_centered <- tcr_dt[clone_id %chin% eligible_clone_ids]
center_columns <- c(
  "spn_detected", "spn_logcounts", "step31_continuum_raw",
  "log1p_endogenous_umis", "step27_program_TCR_activation",
  "step27_program_Cytotoxic_effector", "step27_program_Proliferation"
)
for (column in center_columns) {
  clone_centered[, (paste0(column, "_within_clone")) :=
                   get(column) - mean(get(column), na.rm = TRUE),
                 by = clone_id]
}

fit_clone_centered <- function(dt) {
  informative <- dt[, .(
    outcome_variation = uniqueN(spn_detected) > 1L,
    continuum_variation = stats::sd(step31_continuum_raw, na.rm = TRUE) > 0
  ), by = clone_id]
  informative_ids <- informative[outcome_variation & continuum_variation, clone_id]
  use <- dt[clone_id %chin% informative_ids]
  if (nrow(use) < minimum_tcr_patient_cells ||
      uniqueN(informative_ids) < minimum_informative_clones) return(data.table())

  detection_fit <- tryCatch(stats::lm(
    spn_detected_within_clone ~ step31_continuum_raw_within_clone +
      log1p_endogenous_umis_within_clone +
      step27_program_TCR_activation_within_clone +
      step27_program_Cytotoxic_effector_within_clone +
      step27_program_Proliferation_within_clone,
    data = use
  ), error = function(e) NULL)
  expression_fit <- tryCatch(stats::lm(
    spn_logcounts_within_clone ~ step31_continuum_raw_within_clone +
      log1p_endogenous_umis_within_clone +
      step27_program_TCR_activation_within_clone +
      step27_program_Cytotoxic_effector_within_clone +
      step27_program_Proliferation_within_clone,
    data = use
  ), error = function(e) NULL)

  data.table(
    subject_id = dt$subject_id[1L],
    n_cells = as.integer(nrow(use)),
    n_informative_clones = as.integer(uniqueN(use$clone_id)),
    clone_centered_detection_coefficient = if (is.null(detection_fit)) NA_real_ else
      safe_coef(detection_fit, "step31_continuum_raw_within_clone"),
    clone_centered_expression_coefficient = if (is.null(expression_fit)) NA_real_ else
      safe_coef(expression_fit, "step31_continuum_raw_within_clone")
  )
}

patient_clone_centered <- rbindlist(
  lapply(split(clone_centered, clone_centered$subject_id), fit_clone_centered),
  fill = TRUE
)
clone_centered_summary <- rbindlist(list(
  summarize_effect(
    patient_clone_centered$clone_centered_detection_coefficient,
    "Within-clonotype continuum versus SPN detection",
    "Linear-probability coefficient", 300L
  ),
  summarize_effect(
    patient_clone_centered$clone_centered_expression_coefficient,
    "Within-clonotype continuum versus SPN expression",
    "SPN-logcount coefficient", 301L
  )
), fill = TRUE)
clone_centered_summary[, fdr := p.adjust(wilcoxon_signed_rank_p, method = "BH")]

# Patient-stratified permutation: fine-state labels are shuffled among cells
# while clone membership and clone sizes remain fixed. The expected positive
# size-breadth correlation caused by sampling opportunity is therefore retained.
message("Running clone-size-aware state-sharing permutations")
set.seed(random_seed + 400L)
permutation_patient_results <- list()
null_matrices <- list()
eligible_patients <- unique(tcr_dt$subject_id)

for (patient in eligible_patients) {
  current <- tcr_dt[subject_id == patient]
  sizes <- current[, .N, by = clone_id]
  if (nrow(sizes) < minimum_informative_clones ||
      uniqueN(current$step31_fine_state) < 2L) next

  observed_breadth <- current[, .(
    n_fine_states = as.integer(uniqueN(step31_fine_state))
  ), by = clone_id][sizes, on = "clone_id"]
  observed_rho <- safe_spearman(log2(observed_breadth$N),
                                observed_breadth$n_fine_states)

  clone_vector <- current$clone_id
  state_vector <- current$step31_fine_state
  null_rho <- rep(NA_real_, permutation_iterations)
  for (iteration in seq_len(permutation_iterations)) {
    permuted <- data.table(
      clone_id = clone_vector,
      state = sample(state_vector, length(state_vector), replace = FALSE)
    )[, .(n_fine_states = as.integer(uniqueN(state))), by = clone_id][
      sizes, on = "clone_id"
    ]
    null_rho[iteration] <- safe_spearman(log2(permuted$N),
                                         permuted$n_fine_states)
  }
  null_rho <- null_rho[is.finite(null_rho)]
  if (!is.finite(observed_rho) || !length(null_rho)) next
  excess <- observed_rho - mean(null_rho)
  permutation_p <- (1 + sum(null_rho >= observed_rho)) / (1 + length(null_rho))
  permutation_patient_results[[patient]] <- data.table(
    subject_id = patient,
    n_productive_cells = as.integer(nrow(current)),
    n_clones = as.integer(nrow(sizes)),
    observed_size_breadth_rho = observed_rho,
    mean_null_size_breadth_rho = mean(null_rho),
    excess_size_breadth_rho = excess,
    one_sided_permutation_p = permutation_p
  )
  null_matrices[[patient]] <- null_rho
}

patient_state_sharing_null <- rbindlist(permutation_patient_results, fill = TRUE)
state_sharing_summary <- summarize_effect(
  patient_state_sharing_null$excess_size_breadth_rho,
  "Clone size-state breadth association beyond sampling opportunity",
  "Observed minus permutation-null Spearman rho", 400L
)
state_sharing_summary[, fdr := wilcoxon_signed_rank_p]

# ------------------------------------------------------------------------------
# Evidence ledger, outputs, and figures
# ------------------------------------------------------------------------------

message("Building the Step 32 evidence ledger")
all_summaries <- rbindlist(list(
  continuum_model_summary[, analysis_family := "Continuum"],
  clone_expansion_summary[, analysis_family := "Clone expansion"],
  clone_centered_summary[, analysis_family := "Within clonotype"],
  state_sharing_summary[, analysis_family := "State-sharing null"]
), fill = TRUE)

ledger_map <- data.table(
  hypothesis = c(
    "SPN detection increases along the CD8 dysfunction continuum after biological adjustment",
    "SPN expression increases along the CD8 dysfunction continuum after biological adjustment",
    "The adjusted nonlinear continuum contrast is positive",
    "Expanded clonotypes lie further along the dysfunction continuum",
    "SPN tracks continuum position within the same productive clonotype",
    "Large clones span more fine states than expected from clone size alone"
  ),
  comparison = c(
    "SPN detection across continuum; biological-process adjusted",
    "SPN expression across continuum; biological-process adjusted",
    "Adjusted smooth continuum contrast: 90th minus 10th percentile",
    "Clone size versus median CD8 continuum position",
    "Within-clonotype continuum versus SPN detection",
    "Clone size-state breadth association beyond sampling opportunity"
  )
)
evidence_ledger <- all_summaries[ledger_map, on = "comparison"]
evidence_ledger[, evidence_class := mapply(
  classification, median_effect, fdr, MoreArgs = list(expected_direction = "positive")
)]
setcolorder(evidence_ledger, c(
  "hypothesis", "evidence_class", "analysis_family", "comparison",
  "effect_scale", "n_patients", "median_effect", "bootstrap_ci_lower",
  "bootstrap_ci_upper", "percent_positive", "wilcoxon_signed_rank_p", "fdr"
))

dataset_qc <- data.table(
  metric = c(
    "retained_cd8_cells", "subjects", "patients_with_continuum_models",
    "productive_tcr_cells", "productive_clonotypes",
    "patients_with_clone_centered_models",
    "patients_with_state_sharing_permutations", "permutations_per_patient",
    "spn_used_to_define_continuum_or_states", "signature_exported"
  ),
  value = c(
    nrow(analysis_dt), uniqueN(analysis_dt$subject_id),
    nrow(patient_continuum_models), nrow(tcr_dt), nrow(clone_summary),
    nrow(patient_clone_centered), nrow(patient_state_sharing_null),
    permutation_iterations, 0L, 0L
  )
)

fwrite(dataset_qc, file.path(results_dir, "32_dataset_qc.csv"))
fwrite(patient_continuum_models,
       file.path(results_dir, "32_patient_continuum_models.csv"))
fwrite(continuum_model_summary,
       file.path(results_dir, "32_continuum_model_summary.csv"))
fwrite(clone_summary, file.path(results_dir, "32_productive_clone_summary.csv"))
fwrite(patient_clone_expansion,
       file.path(results_dir, "32_patient_clone_expansion_associations.csv"))
fwrite(clone_expansion_summary,
       file.path(results_dir, "32_clone_expansion_summary.csv"))
fwrite(patient_clone_centered,
       file.path(results_dir, "32_patient_clone_centered_models.csv"))
fwrite(clone_centered_summary,
       file.path(results_dir, "32_clone_centered_summary.csv"))
fwrite(patient_state_sharing_null,
       file.path(results_dir, "32_clone_state_sharing_permutation_by_patient.csv"))
fwrite(state_sharing_summary,
       file.path(results_dir, "32_clone_state_sharing_permutation_summary.csv"))
fwrite(evidence_ledger, file.path(results_dir, "32_evidence_ledger.csv"))

theme_step32 <- theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(color = "grey35"),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

model_long <- melt(
  patient_continuum_models,
  id.vars = c("subject_id", "n_cells"),
  measure.vars = c(
    "minimal_detection_log_odds", "fully_adjusted_detection_log_odds"
  ),
  variable.name = "model", value.name = "coefficient"
)
model_long[, model := factor(model, levels = c(
  "minimal_detection_log_odds", "fully_adjusted_detection_log_odds"
), labels = c("Depth + base state", "Biological-process adjusted"))]
p_coefficients <- ggplot(model_long, aes(model, coefficient, color = model)) +
  geom_hline(yintercept = 0, color = "grey55") +
  geom_boxplot(outlier.shape = NA, width = 0.45) +
  geom_jitter(width = 0.08, height = 0, size = 2, alpha = 0.8) +
  scale_color_manual(values = c("#3B82B8", "#7651B8")) +
  labs(
    title = "Patient-specific SPN association with the CD8 continuum",
    subtitle = "Positive coefficients indicate greater SPN detection toward dysfunction",
    x = NULL, y = "Continuum coefficient (SPN-detection log-odds)"
  ) + theme_step32 + theme(legend.position = "none")
ggsave(file.path(figures_dir, "32A_adjusted_continuum_coefficients.png"),
       p_coefficients, width = 8.5, height = 6, dpi = 300)

p_smooth <- ggplot(
  patient_continuum_models,
  aes(y = reorder(subject_id, adjusted_smooth_q90_minus_q10_detection_pp),
      x = adjusted_smooth_q90_minus_q10_detection_pp)
) +
  geom_vline(xintercept = 0, color = "grey55") +
  geom_point(color = "#7651B8", size = 2.5) +
  labs(
    title = "Adjusted SPN enrichment from low to high continuum position",
    subtitle = "Patient-specific standardized predictions: 90th minus 10th percentile",
    x = "SPN-detection difference (percentage points)", y = "Patient"
  ) + theme_step32
ggsave(file.path(figures_dir, "32B_smooth_continuum_contrasts.png"),
       p_smooth, width = 8.5, height = 6.5, dpi = 300)

if (nrow(patient_state_sharing_null)) {
  p_null <- ggplot(
    patient_state_sharing_null,
    aes(mean_null_size_breadth_rho, observed_size_breadth_rho)
  ) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey45") +
    geom_point(aes(size = n_clones), color = "#7651B8", alpha = 0.8) +
    coord_equal() +
    labs(
      title = "Does clone-state breadth exceed the clone-size null?",
      subtitle = "Each point is one patient; the null preserves clone sizes and state frequencies",
      x = "Mean permutation-null size-breadth rho",
      y = "Observed size-breadth rho", size = "Clones"
    ) + theme_step32
  ggsave(file.path(figures_dir, "32C_clone_state_breadth_null.png"),
         p_null, width = 8, height = 7, dpi = 300)
}

ledger_plot <- copy(evidence_ledger)
ledger_plot[, hypothesis_short := factor(
  c(
    "Adjusted detection-continuum", "Adjusted expression-continuum",
    "Nonlinear Q90-Q10 contrast", "Clone expansion-continuum",
    "Within-clone SPN-continuum", "State breadth beyond size null"
  ), levels = rev(c(
    "Adjusted detection-continuum", "Adjusted expression-continuum",
    "Nonlinear Q90-Q10 contrast", "Clone expansion-continuum",
    "Within-clone SPN-continuum", "State breadth beyond size null"
  ))
)]
p_ledger <- ggplot(ledger_plot, aes(evidence_class, hypothesis_short,
                                    color = evidence_class)) +
  geom_point(size = 4) +
  scale_color_manual(values = c(
    "Supported" = "#2A9D6F", "Suggestive" = "#E9A23B",
    "Not supported" = "#9A9A9A", "Not estimable" = "#D1D1D1"
  )) +
  labs(
    title = "Step 32 evidence ledger",
    subtitle = "Patient-aware consolidation of SPN, continuum, and clonotype hypotheses",
    x = NULL, y = NULL, color = "Evidence"
  ) + theme_step32 + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(file.path(figures_dir, "32D_evidence_ledger.png"),
       p_ledger, width = 10, height = 6, dpi = 300)

S4Vectors::metadata(sce)$step32_evidence_ledger <- evidence_ledger
S4Vectors::metadata(sce)$step32_continuum_model_summary <-
  continuum_model_summary
S4Vectors::metadata(sce)$step32_clone_centered_summary <-
  clone_centered_summary
saveRDS(sce, output_sce_file, compress = FALSE)

message("Step 32 complete: SPN-continuum and clone evidence consolidated")
