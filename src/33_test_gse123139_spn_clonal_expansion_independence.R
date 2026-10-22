#!/usr/bin/env Rscript

# ==============================================================================
# Step 33: Test whether SPN independently tracks CD8 clonal expansion
# ==============================================================================
#
# Primary question
#   Is SPN detection associated with productive TCR clone size after accounting
#   for sequencing depth, continuum position, activation, cytotoxicity,
#   proliferation, and interferon response?
#
# Design
#   - Productive clonotypes are defined within patient.
#   - Primary models are clone-level quasibinomial models, so each clone is one
#     biological unit with SPN-positive and SPN-negative cell counts.
#   - Patients are the independent replication units; patient-specific clone-
#     size coefficients are summarized across patients.
#   - Models are repeated within each inherited CD8 base state.
#   - A one-cell-per-clone resampling sensitivity prevents large clones from
#     dominating simply because they contribute more cells.
#   - SPN is not used to define clone size, continuum position, or CD8 state.
#
# This remains an observational association. Productive clonotype sharing does
# not prove antigen specificity, temporal direction, or a causal role for CD43.
# ==============================================================================

options(stringsAsFactors = FALSE)
options(error = function() {
  cat("\n--- STEP 33 TRACEBACK ---\n", file = stderr())
  traceback(20L)
  quit(save = "no", status = 1L, runLast = FALSE)
})

required_packages <- c(
  "data.table", "SingleCellExperiment", "SummarizedExperiment", "S4Vectors",
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

project_root <- normalizePath(getwd(), mustWork = TRUE)
analysis_dir <- file.path(
  project_root, "data", "analysis_ready", "scRNA", "GSE123139"
)
results_dir <- file.path(project_root, "results", "33")
figures_dir <- file.path(results_dir, "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

input_sce_file <- file.path(
  analysis_dir, "32_GSE123139_spn_continuum_clone_consolidated_sce.rds"
)
output_sce_file <- file.path(
  analysis_dir, "33_GSE123139_spn_clonal_expansion_sce.rds"
)

random_seed <- 123139L
bootstrap_iterations <- 5000L
resampling_iterations <- 250L
minimum_clones_per_patient_model <- 25L
minimum_clones_per_state_model <- 20L
minimum_positive_and_negative_cells <- 5L
minimum_resampling_success_fraction <- 0.50

base_states <- c("CD8_memory_transitional", "CD8_dysfunctional")
base_state_labels <- c(
  CD8_memory_transitional = "Memory/transitional",
  CD8_dysfunctional = "Dysfunctional"
)
program_names <- c(
  "TCR_activation", "Cytotoxic_effector", "Proliferation",
  "Interferon_response"
)
program_columns <- paste0("step27_program_", program_names)

as_flag <- function(x) !is.na(x) & as.logical(x)

zscore_safe <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

safe_coef <- function(fit, term) {
  output <- tryCatch(stats::coef(fit)[term], error = function(e) NA_real_)
  if (!length(output) || !is.finite(output)) NA_real_ else as.numeric(output)
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

evidence_class <- function(effect, fdr) {
  if (!is.finite(effect) || !is.finite(fdr)) return("Not estimable")
  if (effect > 0 && fdr < 0.05) return("Supported")
  if (effect > 0 && fdr < 0.10) return("Suggestive")
  "Not supported"
}

message("Loading completed Step 32 object")
if (!file.exists(input_sce_file)) stop("Missing Step 32 SCE: ", input_sce_file)
sce <- readRDS(input_sce_file)

required_metadata <- c(
  "subject_id", "step25_retained_cd8", "step25_base_state",
  "endogenous_n_counts", "step26_spn_detected", "step26_spn_logcounts",
  "step31_continuum_raw", "step31_fine_state",
  "tcr_productive_primary", "tcr_clonotype_key", program_columns
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
cell[, productive_tcr := as_flag(tcr_productive_primary)]
cell[, spn_detected := as.integer(as_flag(step26_spn_detected))]
cell[, spn_logcounts := as.numeric(step26_spn_logcounts)]
cell[, log1p_endogenous_umis := log1p(as.numeric(endogenous_n_counts))]

tcr_dt <- cell[
  retained_cd8 == TRUE & step25_base_state %chin% base_states &
    productive_tcr == TRUE & !is.na(tcr_clonotype_key) &
    nzchar(tcr_clonotype_key) & is.finite(step31_continuum_raw) &
    is.finite(log1p_endogenous_umis)
]
if (!nrow(tcr_dt)) stop("No retained productive-TCR CD8 cells were found.")
tcr_dt[, clone_id := paste(subject_id, tcr_clonotype_key, sep = "::")]

# ------------------------------------------------------------------------------
# Clone-level summaries and descriptive category contrasts
# ------------------------------------------------------------------------------

message("Constructing clone-level SPN and biological summaries")
clone_summary <- tcr_dt[, .(
  clone_size = as.integer(.N),
  n_spn_detected = as.integer(sum(spn_detected)),
  n_spn_undetected = as.integer(sum(spn_detected == 0L)),
  spn_detected_percent = 100 * mean(spn_detected),
  mean_spn_logcounts = mean(spn_logcounts),
  mean_log1p_endogenous_umis = mean(log1p_endogenous_umis),
  mean_continuum = mean(step31_continuum_raw),
  mean_tcr_activation = mean(step27_program_TCR_activation),
  mean_cytotoxicity = mean(step27_program_Cytotoxic_effector),
  mean_proliferation = mean(step27_program_Proliferation),
  mean_interferon = mean(step27_program_Interferon_response),
  n_base_states = as.integer(uniqueN(step25_base_state)),
  n_fine_states = as.integer(uniqueN(step31_fine_state))
), by = .(subject_id, clone_id)]
clone_summary[, log2_clone_size := log2(clone_size)]
clone_summary[, expansion_class := factor(fifelse(
  clone_size >= 10L, "Large (>=10 cells)",
  fifelse(clone_size >= 3L, "Expanded (3-9 cells)", "Small (1-2 cells)")
), levels = c("Small (1-2 cells)", "Expanded (3-9 cells)",
             "Large (>=10 cells)"))]

patient_category_summary <- clone_summary[, .(
  n_clones = as.integer(.N),
  n_cells = as.integer(sum(clone_size)),
  clone_weighted_spn_detected_percent = 100 * mean(n_spn_detected / clone_size),
  cell_weighted_spn_detected_percent = 100 * sum(n_spn_detected) / sum(clone_size),
  median_clone_continuum = as.numeric(stats::median(mean_continuum)),
  median_clone_depth = as.numeric(stats::median(mean_log1p_endogenous_umis))
), by = .(subject_id, expansion_class)]

# ------------------------------------------------------------------------------
# Patient-specific clone-level quasibinomial models
# ------------------------------------------------------------------------------

message("Fitting adjusted clone-level quasibinomial models")
fit_clone_model <- function(dt, minimum_clones, stratum = "All productive CD8") {
  if (nrow(dt) < minimum_clones ||
      sum(dt$n_spn_detected) < minimum_positive_and_negative_cells ||
      sum(dt$n_spn_undetected) < minimum_positive_and_negative_cells ||
      uniqueN(dt$clone_size) < 2L) return(data.table())

  model_columns <- c(
    "log2_clone_size", "mean_log1p_endogenous_umis", "mean_continuum",
    "mean_tcr_activation", "mean_cytotoxicity", "mean_proliferation",
    "mean_interferon"
  )
  for (column in model_columns) {
    dt[, (paste0(column, "_z")) := zscore_safe(get(column))]
  }

  depth_fit <- tryCatch(stats::glm(
    cbind(n_spn_detected, n_spn_undetected) ~
      log2_clone_size_z + mean_log1p_endogenous_umis_z,
    family = stats::quasibinomial(), data = dt
  ), error = function(e) NULL)
  adjusted_fit <- tryCatch(stats::glm(
    cbind(n_spn_detected, n_spn_undetected) ~
      log2_clone_size_z + mean_log1p_endogenous_umis_z + mean_continuum_z +
      mean_tcr_activation_z + mean_cytotoxicity_z + mean_proliferation_z +
      mean_interferon_z,
    family = stats::quasibinomial(), data = dt
  ), error = function(e) NULL)

  data.table(
    subject_id = dt$subject_id[1L],
    stratum = stratum,
    n_clones = as.integer(nrow(dt)),
    n_cells = as.integer(sum(dt$clone_size)),
    depth_adjusted_clone_size_log_odds = if (is.null(depth_fit)) NA_real_ else
      safe_coef(depth_fit, "log2_clone_size_z"),
    biologically_adjusted_clone_size_log_odds = if (is.null(adjusted_fit)) NA_real_ else
      safe_coef(adjusted_fit, "log2_clone_size_z")
  )
}

patient_primary_models <- rbindlist(lapply(
  split(clone_summary, clone_summary$subject_id),
  fit_clone_model, minimum_clones = minimum_clones_per_patient_model
), fill = TRUE)

# State-specific clone summaries use total clone size as the expansion measure,
# while SPN counts and covariates are calculated only from cells in that state.
clone_sizes <- clone_summary[, .(clone_id, total_clone_size = clone_size)]
clone_state_summary <- tcr_dt[, .(
  state_cells = as.integer(.N),
  n_spn_detected = as.integer(sum(spn_detected)),
  n_spn_undetected = as.integer(sum(spn_detected == 0L)),
  mean_log1p_endogenous_umis = mean(log1p_endogenous_umis),
  mean_continuum = mean(step31_continuum_raw),
  mean_tcr_activation = mean(step27_program_TCR_activation),
  mean_cytotoxicity = mean(step27_program_Cytotoxic_effector),
  mean_proliferation = mean(step27_program_Proliferation),
  mean_interferon = mean(step27_program_Interferon_response)
), by = .(subject_id, clone_id, step25_base_state)][
  clone_sizes, on = "clone_id", nomatch = 0L
]
clone_state_summary[, clone_size := total_clone_size]
clone_state_summary[, log2_clone_size := log2(clone_size)]

patient_state_models <- rbindlist(lapply(base_states, function(state) {
  current <- clone_state_summary[step25_base_state == state]
  rbindlist(lapply(
    split(current, current$subject_id), fit_clone_model,
    minimum_clones = minimum_clones_per_state_model,
    stratum = unname(base_state_labels[state])
  ), fill = TRUE)
}), fill = TRUE)

all_patient_models <- rbindlist(list(
  patient_primary_models, patient_state_models
), fill = TRUE)

model_summary <- rbindlist(lapply(unique(all_patient_models$stratum), function(s) {
  current <- all_patient_models[stratum == s]
  rbindlist(list(
    summarize_effect(
      current$depth_adjusted_clone_size_log_odds,
      paste0(s, ": clone size versus SPN detection; depth adjusted"),
      "SPN-detection log-odds per SD log2 clone size",
      100L + match(s, unique(all_patient_models$stratum))
    ),
    summarize_effect(
      current$biologically_adjusted_clone_size_log_odds,
      paste0(s, ": clone size versus SPN detection; biological-process adjusted"),
      "SPN-detection log-odds per SD log2 clone size",
      200L + match(s, unique(all_patient_models$stratum))
    )
  ))[, stratum := s]
}), fill = TRUE)
model_summary[, fdr := p.adjust(wilcoxon_signed_rank_p, method = "BH")]

# ------------------------------------------------------------------------------
# One-cell-per-clone resampling sensitivity
# ------------------------------------------------------------------------------

message("Running one-cell-per-clone resampling sensitivity")
tcr_with_size <- tcr_dt[clone_summary[, .(
  clone_id, clone_size, log2_clone_size, expansion_class
)], on = "clone_id", nomatch = 0L]

fit_one_cell_iteration <- function(dt) {
  sampled <- dt[, .SD[sample.int(.N, 1L)], by = clone_id]
  if (nrow(sampled) < minimum_clones_per_patient_model ||
      uniqueN(sampled$spn_detected) < 2L || uniqueN(sampled$clone_size) < 2L) {
    return(c(
      coefficient = NA_real_, expanded_minus_small_pp = NA_real_,
      expanded_minus_singleton_pp = NA_real_
    ))
  }
  sampled[, clone_size_z := zscore_safe(log2_clone_size)]
  sampled[, depth_z := zscore_safe(log1p_endogenous_umis)]
  sampled[, continuum_z := zscore_safe(step31_continuum_raw)]
  sampled[, activation_z := zscore_safe(step27_program_TCR_activation)]
  sampled[, cytotoxicity_z := zscore_safe(step27_program_Cytotoxic_effector)]
  sampled[, proliferation_z := zscore_safe(step27_program_Proliferation)]
  sampled[, interferon_z := zscore_safe(step27_program_Interferon_response)]
  fit <- tryCatch(stats::glm(
    spn_detected ~ clone_size_z + depth_z + continuum_z + activation_z +
      cytotoxicity_z + proliferation_z + interferon_z,
    family = stats::binomial(), data = sampled
  ), error = function(e) NULL)
  coefficient <- if (is.null(fit)) NA_real_ else safe_coef(fit, "clone_size_z")
  expanded <- sampled[clone_size >= 3L, mean(spn_detected)]
  small <- sampled[clone_size <= 2L, mean(spn_detected)]
  singleton <- sampled[clone_size == 1L, mean(spn_detected)]
  difference <- if (is.finite(expanded) && is.finite(small)) {
    100 * (expanded - small)
  } else NA_real_
  singleton_difference <- if (is.finite(expanded) && is.finite(singleton)) {
    100 * (expanded - singleton)
  } else NA_real_
  c(
    coefficient = coefficient,
    expanded_minus_small_pp = difference,
    expanded_minus_singleton_pp = singleton_difference
  )
}

set.seed(random_seed + 300L)
resampling_results <- rbindlist(lapply(
  split(tcr_with_size, tcr_with_size$subject_id),
  function(current) {
    iterations <- replicate(
      resampling_iterations, fit_one_cell_iteration(current), simplify = "matrix"
    )
    successful <- is.finite(iterations["coefficient", ])
    data.table(
      subject_id = current$subject_id[1L],
      n_clones = as.integer(uniqueN(current$clone_id)),
      successful_iterations = as.integer(sum(successful)),
      success_fraction = mean(successful),
      median_adjusted_clone_size_log_odds = if (any(successful))
        as.numeric(stats::median(iterations["coefficient", successful])) else NA_real_,
      median_expanded_minus_small_detection_pp = if (any(successful))
        as.numeric(stats::median(
          iterations["expanded_minus_small_pp", successful], na.rm = TRUE
        )) else NA_real_,
      median_expanded_minus_singleton_detection_pp = if (any(successful))
        as.numeric(stats::median(
          iterations["expanded_minus_singleton_pp", successful], na.rm = TRUE
        )) else NA_real_
    )
  }
), fill = TRUE)
resampling_eligible <- resampling_results[
  success_fraction >= minimum_resampling_success_fraction
]
resampling_summary <- rbindlist(list(
  summarize_effect(
    resampling_eligible$median_adjusted_clone_size_log_odds,
    "One-cell-per-clone: adjusted clone size versus SPN detection",
    "Median resampled SPN-detection log-odds per SD log2 clone size", 400L
  ),
  summarize_effect(
    resampling_eligible$median_expanded_minus_small_detection_pp,
    "One-cell-per-clone: expanded minus small clone SPN detection",
    "Median resampled SPN-detection percentage points", 401L
  ),
  summarize_effect(
    resampling_eligible$median_expanded_minus_singleton_detection_pp,
    "One-cell-per-clone: expanded minus singleton SPN detection",
    "Median resampled SPN-detection percentage points", 402L
  )
), fill = TRUE)
resampling_summary[, fdr := p.adjust(wilcoxon_signed_rank_p, method = "BH")]

# ------------------------------------------------------------------------------
# Evidence ledger, output, and figures
# ------------------------------------------------------------------------------

message("Building Step 33 evidence ledger")
primary_adjusted_label <- paste0(
  "All productive CD8: clone size versus SPN detection; biological-process adjusted"
)
primary_depth_label <- paste0(
  "All productive CD8: clone size versus SPN detection; depth adjusted"
)
ledger_sources <- rbindlist(list(model_summary, resampling_summary), fill = TRUE)
ledger_map <- data.table(
  hypothesis = c(
    "SPN detection is associated with clone size after depth adjustment",
    "SPN detection is associated with clone size independently of state and biological programs",
    "The adjusted clone-size association persists when each clone contributes one cell",
    "Expanded clones have higher SPN detection when each clone contributes one cell",
    "Expanded clones have higher SPN detection than productive TCR singletons",
    "The independent clone-size association is present within memory/transitional cells",
    "The independent clone-size association is present within dysfunctional cells"
  ),
  comparison = c(
    primary_depth_label,
    primary_adjusted_label,
    "One-cell-per-clone: adjusted clone size versus SPN detection",
    "One-cell-per-clone: expanded minus small clone SPN detection",
    "One-cell-per-clone: expanded minus singleton SPN detection",
    "Memory/transitional: clone size versus SPN detection; biological-process adjusted",
    "Dysfunctional: clone size versus SPN detection; biological-process adjusted"
  )
)
evidence_ledger <- ledger_sources[ledger_map, on = "comparison"]
evidence_ledger[, evidence_class := mapply(
  evidence_class, median_effect, fdr
)]
setcolorder(evidence_ledger, c(
  "hypothesis", "evidence_class", "comparison", "stratum", "effect_scale",
  "n_patients", "median_effect", "bootstrap_ci_lower", "bootstrap_ci_upper",
  "percent_positive", "wilcoxon_signed_rank_p", "fdr"
))

dataset_qc <- data.table(
  metric = c(
    "productive_tcr_cells", "productive_clonotypes", "subjects",
    "small_clonotypes", "expanded_clonotypes", "large_clonotypes",
    "patients_with_primary_clone_models", "patients_with_state_models",
    "patients_with_resampling_sensitivity", "resampling_iterations_per_patient",
    "spn_used_to_define_clone_size_or_state", "signature_exported"
  ),
  value = c(
    nrow(tcr_dt), nrow(clone_summary), uniqueN(tcr_dt$subject_id),
    clone_summary[clone_size <= 2L, .N],
    clone_summary[clone_size >= 3L & clone_size < 10L, .N],
    clone_summary[clone_size >= 10L, .N],
    nrow(patient_primary_models), nrow(patient_state_models),
    nrow(resampling_eligible), resampling_iterations, 0L, 0L
  )
)

fwrite(dataset_qc, file.path(results_dir, "33_dataset_qc.csv"))
fwrite(clone_summary, file.path(results_dir, "33_productive_clone_summary.csv"))
fwrite(patient_category_summary,
       file.path(results_dir, "33_patient_clone_category_summary.csv"))
fwrite(all_patient_models,
       file.path(results_dir, "33_patient_adjusted_clone_models.csv"))
fwrite(model_summary, file.path(results_dir, "33_adjusted_clone_model_summary.csv"))
fwrite(resampling_results,
       file.path(results_dir, "33_one_cell_per_clone_resampling_by_patient.csv"))
fwrite(resampling_summary,
       file.path(results_dir, "33_one_cell_per_clone_resampling_summary.csv"))
fwrite(evidence_ledger, file.path(results_dir, "33_evidence_ledger.csv"))

theme_step33 <- theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(color = "grey35"),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

category_plot_dt <- patient_category_summary[
  !is.na(expansion_class) & n_clones >= 2L
]
p_categories <- ggplot(
  category_plot_dt,
  aes(expansion_class, clone_weighted_spn_detected_percent,
      group = subject_id)
) +
  geom_line(color = "grey75", alpha = 0.7) +
  geom_point(aes(color = expansion_class), size = 2.5) +
  scale_color_manual(values = c(
    "Small (1-2 cells)" = "#A7A7A7", "Expanded (3-9 cells)" = "#3B82B8",
    "Large (>=10 cells)" = "#7651B8"
  )) +
  labs(
    title = "Patient-paired SPN detection across productive clone sizes",
    subtitle = "Each clone contributes equally within each patient and category",
    x = NULL, y = "Clone-weighted SPN-detected cells (%)", color = "Clone size"
  ) + theme_step33 + theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(figures_dir, "33A_spn_by_clone_size_category.png"),
       p_categories, width = 9, height = 6.5, dpi = 300)

primary_plot <- all_patient_models[stratum == "All productive CD8"]
primary_long <- melt(
  primary_plot,
  id.vars = c("subject_id", "n_clones", "n_cells"),
  measure.vars = c(
    "depth_adjusted_clone_size_log_odds",
    "biologically_adjusted_clone_size_log_odds"
  ),
  variable.name = "model", value.name = "coefficient"
)
primary_long[, model := factor(model, levels = c(
  "depth_adjusted_clone_size_log_odds",
  "biologically_adjusted_clone_size_log_odds"
), labels = c("Depth adjusted", "Biological-process adjusted"))]
p_coefficients <- ggplot(primary_long, aes(model, coefficient, color = model)) +
  geom_hline(yintercept = 0, color = "grey55") +
  geom_boxplot(outlier.shape = NA, width = 0.45) +
  geom_jitter(width = 0.08, height = 0, size = 2, alpha = 0.8) +
  scale_color_manual(values = c("#3B82B8", "#7651B8")) +
  labs(
    title = "Does SPN independently track productive clone size?",
    subtitle = "Patient-specific clone-level quasibinomial coefficients",
    x = NULL, y = "Clone-size coefficient (SPN-detection log-odds)"
  ) + theme_step33 + theme(legend.position = "none")
ggsave(file.path(figures_dir, "33B_adjusted_clone_size_coefficients.png"),
       p_coefficients, width = 8.5, height = 6, dpi = 300)

state_plot <- all_patient_models[stratum != "All productive CD8"]
p_states <- ggplot(
  state_plot,
  aes(stratum, biologically_adjusted_clone_size_log_odds, color = stratum)
) +
  geom_hline(yintercept = 0, color = "grey55") +
  geom_boxplot(outlier.shape = NA, width = 0.45) +
  geom_jitter(width = 0.08, height = 0, size = 2, alpha = 0.8) +
  scale_color_manual(values = c(
    "Memory/transitional" = "#3B82B8", "Dysfunctional" = "#7651B8"
  )) +
  labs(
    title = "Adjusted SPN-clone-size association within CD8 base states",
    subtitle = "Positive values indicate higher SPN detection in larger clones",
    x = NULL, y = "Biological-process-adjusted clone-size coefficient"
  ) + theme_step33 + theme(legend.position = "none")
ggsave(file.path(figures_dir, "33C_state_stratified_clone_coefficients.png"),
       p_states, width = 8.5, height = 6, dpi = 300)

p_resampling <- ggplot(
  resampling_eligible,
  aes(y = reorder(subject_id, median_adjusted_clone_size_log_odds),
      x = median_adjusted_clone_size_log_odds)
) +
  geom_vline(xintercept = 0, color = "grey55") +
  geom_point(color = "#7651B8", size = 2.5) +
  labs(
    title = "One-cell-per-clone sensitivity",
    subtitle = paste0("Patient medians across ", resampling_iterations,
                      " resamples; every clone contributes one cell"),
    x = "Adjusted clone-size coefficient", y = "Patient"
  ) + theme_step33
ggsave(file.path(figures_dir, "33D_one_cell_per_clone_sensitivity.png"),
       p_resampling, width = 8.5, height = 6.5, dpi = 300)

S4Vectors::metadata(sce)$step33_adjusted_clone_model_summary <- model_summary
S4Vectors::metadata(sce)$step33_resampling_summary <- resampling_summary
S4Vectors::metadata(sce)$step33_evidence_ledger <- evidence_ledger
saveRDS(sce, output_sce_file, compress = FALSE)

message("Step 33 complete: independent SPN-clonal expansion association tested")
