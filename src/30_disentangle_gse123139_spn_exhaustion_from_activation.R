#!/usr/bin/env Rscript

# ==============================================================================
# Step 30: Does SPN track CD8 exhaustion independently of activation?
# ==============================================================================
#
# This step distinguishes three biological explanations for the Step 26 result:
#   1. SPN specifically tracks exhaustion/dysfunction;
#   2. SPN tracks activation/cytotoxicity and only indirectly appears related to
#      exhaustion;
#   3. SPN marks state composition but not variation within a state.
#
# Patients, not cells, are the inferential replicates. Models are fitted within
# each patient, and the resulting SPN coefficients are summarized across
# patients. SPN was not used to define the states or any program score.
#
# Primary outcomes:
#   - Exhaustion_core
#   - Terminal_exhaustion
#
# Model sequence:
#   basic: SPN + depth (+ base state in the all-CD8 analysis)
#   activation_cytotoxicity_adjusted: basic + activation + cytotoxicity
#   process_adjusted: activation/cytotoxicity model + tissue residency +
#                     interferon + proliferation + stress
#   differentiation_adjusted: process_adjusted + progenitor dysfunction
#
# Scopes:
#   - all retained CD8 cells
#   - dysfunctional cells only
#   - memory/transitional cells only
#
# Primary exposure is SPN detection, appropriate for sparse MARS-seq data.
# Continuous SPN logcounts are analyzed as a sensitivity.
# ==============================================================================

options(stringsAsFactors = FALSE)

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
results_dir <- file.path(project_root, "results", "30")
figures_dir <- file.path(results_dir, "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

input_sce_file <- file.path(
  analysis_dir, "28_GSE123139_tumor_cd8_within_state_spn_biology_sce.rds"
)
output_rds_file <- file.path(
  analysis_dir, "30_GSE123139_spn_exhaustion_conditional_models.rds"
)

set.seed(123139L)
bootstrap_iterations <- 5000L
minimum_all_cd8_cells <- 80L
minimum_within_state_cells <- 50L
minimum_spn_group_cells <- 10L
maximum_vif_for_sensitivity <- 10

base_states <- c("CD8_memory_transitional", "CD8_dysfunctional")
outcomes <- c("Exhaustion_core", "Terminal_exhaustion")
process_adjusters <- c(
  "TCR_activation", "Cytotoxic_effector", "Tissue_resident",
  "Interferon_response", "Proliferation", "Stress"
)

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
  set.seed(123139L + as.integer(seed_offset))
  boot <- replicate(
    bootstrap_iterations,
    stats::median(sample(x, length(x), replace = TRUE))
  )
  setNames(
    as.numeric(stats::quantile(boot, c(0.025, 0.975), names = FALSE)),
    c("lower", "upper")
  )
}

zscore <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

calculate_max_vif <- function(model_matrix) {
  x <- model_matrix[, colnames(model_matrix) != "(Intercept)", drop = FALSE]
  if (ncol(x) < 2L) return(1)
  vifs <- vapply(seq_len(ncol(x)), function(j) {
    response <- x[, j]
    predictors <- x[, -j, drop = FALSE]
    if (stats::sd(response) == 0) return(Inf)
    fit <- tryCatch(stats::lm.fit(cbind(1, predictors), response),
                    error = function(e) NULL)
    if (is.null(fit)) return(Inf)
    residual_ss <- sum(fit$residuals^2)
    total_ss <- sum((response - mean(response))^2)
    if (!is.finite(total_ss) || total_ss <= 0) return(Inf)
    r_squared <- 1 - residual_ss / total_ss
    if (!is.finite(r_squared) || r_squared >= 1) Inf else 1 / (1 - r_squared)
  }, numeric(1L))
  max(vifs, na.rm = TRUE)
}

message("Loading completed Step 28 object")
if (!file.exists(input_sce_file)) stop("Missing Step 28 SCE: ", input_sce_file)
sce <- readRDS(input_sce_file)

program_columns <- paste0(
  "step27_program_",
  unique(c(outcomes, process_adjusters, "Progenitor_dysfunction", "Memory_naive"))
)
required_metadata <- c(
  "subject_id", "step25_retained_cd8", "step25_base_state",
  "endogenous_n_counts", "endogenous_n_genes", "step26_spn_detected",
  "step26_spn_logcounts", program_columns
)
missing_metadata <- setdiff(
  required_metadata, colnames(SummarizedExperiment::colData(sce))
)
if (length(missing_metadata)) {
  stop("Step 28 metadata is missing: ", paste(missing_metadata, collapse = ", "))
}

cell <- as.data.table(SummarizedExperiment::colData(sce))
cell[, cell_id := colnames(sce)]
cell[, subject_id := as.character(subject_id)]
cell[, step25_base_state := as.character(step25_base_state)]
cell[, retained_cd8 := !is.na(step25_retained_cd8) &
       as.logical(step25_retained_cd8)]
cell[, spn_detected := as.integer(
  !is.na(step26_spn_detected) & as.logical(step26_spn_detected)
)]
cell[, spn_logcounts := as.numeric(step26_spn_logcounts)]
cell[, log1p_endogenous_umis := log1p(as.numeric(endogenous_n_counts))]
cell[, log1p_endogenous_genes := log1p(as.numeric(endogenous_n_genes))]
cell[, dysfunctional := as.integer(step25_base_state == "CD8_dysfunctional")]

analysis_dt <- cell[
  retained_cd8 == TRUE & step25_base_state %chin% base_states
]
if (!nrow(analysis_dt)) stop("No retained Step 25 CD8 cells were found.")

model_specifications <- c(
  "basic", "activation_cytotoxicity_adjusted", "process_adjusted",
  "differentiation_adjusted"
)

adjusters_for_specification <- function(specification) {
  switch(
    specification,
    basic = character(0L),
    activation_cytotoxicity_adjusted = c(
      "TCR_activation", "Cytotoxic_effector"
    ),
    process_adjusted = process_adjusters,
    differentiation_adjusted = c(
      process_adjusters, "Progenitor_dysfunction"
    ),
    stop("Unknown model specification: ", specification)
  )
}

scopes <- list(
  all_retained_cd8 = list(state = NA_character_, minimum = minimum_all_cd8_cells),
  dysfunctional_only = list(
    state = "CD8_dysfunctional", minimum = minimum_within_state_cells
  ),
  memory_transitional_only = list(
    state = "CD8_memory_transitional", minimum = minimum_within_state_cells
  )
)

prepare_patient_table <- function(dt, outcome, exposure, adjusters, include_state) {
  outcome_column <- paste0("step27_program_", outcome)
  out <- copy(dt)
  out[, outcome_z := zscore(get(outcome_column))]
  if (exposure == "spn_detected") {
    out[, spn_exposure := as.numeric(spn_detected)]
  } else {
  out[, spn_exposure := zscore(spn_logcounts)]
  }
  out[, depth_umis_z := zscore(log1p_endogenous_umis)]
  # Gene complexity and UMI depth are strongly correlated. Use the component
  # of detected-gene complexity not explained by UMI depth to avoid entering
  # two nearly redundant technical covariates into every patient model.
  complexity_fit <- tryCatch(
    stats::lm(log1p_endogenous_genes ~ log1p_endogenous_umis, data = out),
    error = function(e) NULL
  )
  if (is.null(complexity_fit)) {
    out[, depth_complexity_z := NA_real_]
  } else {
    out[, depth_complexity_z := zscore(stats::residuals(complexity_fit))]
  }
  adjusted_columns <- character()
  for (program in adjusters) {
    source <- paste0("step27_program_", program)
    destination <- paste0("adjuster_", program)
    out[, (destination) := zscore(get(source))]
    adjusted_columns <- c(adjusted_columns, destination)
  }
  variables <- c(
    "outcome_z", "spn_exposure", "depth_umis_z", "depth_complexity_z",
    adjusted_columns
  )
  if (include_state) variables <- c(variables, "dysfunctional")
  out[stats::complete.cases(out[, ..variables])]
}

fit_patient_model <- function(dt, outcome, exposure, specification,
                              adjusters, include_state, minimum_cells) {
  n_positive <- sum(dt$spn_detected == 1L)
  n_negative <- sum(dt$spn_detected == 0L)
  if (nrow(dt) < minimum_cells ||
      n_positive < minimum_spn_group_cells ||
      n_negative < minimum_spn_group_cells ||
      (include_state && uniqueN(dt$dysfunctional) < 2L)) {
    return(data.table())
  }

  model_dt <- prepare_patient_table(
    dt, outcome, exposure, adjusters, include_state
  )
  if (nrow(model_dt) < minimum_cells) return(data.table())

  terms <- c(
    "spn_exposure", "depth_umis_z", "depth_complexity_z"
  )
  if (length(adjusters)) terms <- c(terms, paste0("adjuster_", adjusters))
  if (include_state) terms <- c(terms, "dysfunctional")
  formula <- stats::as.formula(paste("outcome_z ~", paste(terms, collapse = " + ")))
  design <- tryCatch(stats::model.matrix(formula, model_dt),
                     error = function(e) NULL)
  if (is.null(design) || qr(design)$rank < ncol(design)) return(data.table())

  fit <- tryCatch(stats::lm(formula, data = model_dt), error = function(e) NULL)
  if (is.null(fit) || !"spn_exposure" %in% names(stats::coef(fit))) {
    return(data.table())
  }
  coefficient <- unname(stats::coef(fit)[["spn_exposure"]])
  if (!is.finite(coefficient)) return(data.table())

  data.table(
    n_cells = as.integer(nrow(model_dt)),
    n_spn_detected = as.integer(sum(model_dt$spn_detected == 1L)),
    n_spn_not_detected = as.integer(sum(model_dt$spn_detected == 0L)),
    spn_coefficient = coefficient,
    model_r_squared = summary(fit)$r.squared,
    max_vif = calculate_max_vif(design),
    low_collinearity = calculate_max_vif(design) <= maximum_vif_for_sensitivity
  )
}

message("Fitting patient-specific conditional exhaustion models")
patient_model_rows <- list()
counter <- 0L
for (scope_name in names(scopes)) {
  scope_info <- scopes[[scope_name]]
  scope_dt <- if (is.na(scope_info$state)) {
    analysis_dt
  } else {
    analysis_dt[step25_base_state == scope_info$state]
  }
  include_state <- is.na(scope_info$state)

  for (outcome in outcomes) {
    for (exposure in c("spn_detected", "spn_logcounts")) {
      for (specification in model_specifications) {
        adjusters <- adjusters_for_specification(specification)
        temp <- scope_dt[, fit_patient_model(
          .SD, outcome, exposure, specification, adjusters,
          include_state, scope_info$minimum
        ), by = subject_id]
        if (!nrow(temp)) next
        temp[, `:=`(
          scope = scope_name,
          outcome = outcome,
          exposure = exposure,
          specification = specification
        )]
        counter <- counter + 1L
        patient_model_rows[[counter]] <- temp
      }
    }
  }
}
patient_models <- rbindlist(patient_model_rows, fill = TRUE)
if (!nrow(patient_models)) stop("No patient model passed the eligibility rules.")

summarize_coefficients <- function(dt, seed_offset) {
  values <- dt$spn_coefficient[is.finite(dt$spn_coefficient)]
  ci <- bootstrap_median_ci(values, seed_offset)
  data.table(
    n_patients = as.integer(length(values)),
    median_spn_coefficient = if (length(values)) median(values) else NA_real_,
    bootstrap_ci_lower = unname(ci[["lower"]]),
    bootstrap_ci_upper = unname(ci[["upper"]]),
    percent_positive = if (length(values)) 100 * mean(values > 0) else NA_real_,
    wilcoxon_signed_rank_p = safe_wilcox_p(values),
    median_max_vif = if (nrow(dt)) median(dt$max_vif, na.rm = TRUE) else NA_real_,
    percent_low_collinearity = if (nrow(dt)) 100 * mean(dt$low_collinearity) else NA_real_
  )
}

model_summary <- patient_models[,
  summarize_coefficients(.SD, 1000L + .GRP),
  by = .(scope, outcome, exposure, specification)
]
model_summary[, fdr := p.adjust(wilcoxon_signed_rank_p, method = "BH"),
              by = .(exposure, specification)]
model_summary[, evidence := fifelse(
  fdr < 0.05 & bootstrap_ci_lower > 0, "Independent positive association",
  fifelse(
    fdr < 0.05 & bootstrap_ci_upper < 0, "Independent negative association",
    "No independent association"
  )
)]

# A stricter sensitivity excludes patient fits with VIF > 10.
low_vif_summary <- patient_models[low_collinearity == TRUE,
  summarize_coefficients(.SD, 3000L + .GRP),
  by = .(scope, outcome, exposure, specification)
]
low_vif_summary[, fdr := p.adjust(wilcoxon_signed_rank_p, method = "BH"),
                by = .(exposure, specification)]

# Quantify how much the SPN coefficient changes after biological adjustment in
# the same patient. A large shift toward zero supports activation confounding.
basic_effects <- patient_models[specification == "basic", .(
  subject_id, scope, outcome, exposure,
  basic_spn_coefficient = spn_coefficient
)]
adjusted_effects <- patient_models[specification != "basic", .(
  subject_id, scope, outcome, exposure, specification,
  adjusted_spn_coefficient = spn_coefficient
)]
attenuation <- merge(
  adjusted_effects, basic_effects,
  by = c("subject_id", "scope", "outcome", "exposure"), all = FALSE
)
attenuation[, coefficient_change :=
  adjusted_spn_coefficient - basic_spn_coefficient]
attenuation[, absolute_effect_retained_percent := fifelse(
  abs(basic_spn_coefficient) > 1e-8,
  100 * abs(adjusted_spn_coefficient) / abs(basic_spn_coefficient),
  NA_real_
)]

attenuation_summary <- attenuation[, .(
  n_patients = .N,
  median_basic_coefficient = median(basic_spn_coefficient),
  median_adjusted_coefficient = median(adjusted_spn_coefficient),
  median_coefficient_change = median(coefficient_change),
  median_absolute_effect_retained_percent = median(
    absolute_effect_retained_percent, na.rm = TRUE
  )
), by = .(scope, outcome, exposure, specification)]

# Correlations show which biological processes are hardest to separate.
correlation_columns <- c(
  paste0("step27_program_", c(outcomes, process_adjusters,
                              "Progenitor_dysfunction", "Memory_naive")),
  "step26_spn_logcounts"
)
program_correlation <- stats::cor(
  as.matrix(analysis_dt[, ..correlation_columns]),
  method = "spearman", use = "pairwise.complete.obs"
)
program_correlation_long <- as.data.table(as.table(program_correlation))
setnames(program_correlation_long, c("variable_1", "variable_2", "spearman_rho"))

dataset_qc <- data.table(
  metric = c(
    "retained_cd8_cells", "subjects", "memory_transitional_cells",
    "dysfunctional_cells", "spn_detected_cells", "spn_not_detected_cells",
    "outcomes_tested", "model_specifications", "analysis_scopes",
    "patient_models_fit", "patients_in_primary_all_cd8_full_model",
    "spn_used_to_define_states_or_programs"
  ),
  value = as.character(c(
    nrow(analysis_dt), uniqueN(analysis_dt$subject_id),
    analysis_dt[step25_base_state == "CD8_memory_transitional", .N],
    analysis_dt[step25_base_state == "CD8_dysfunctional", .N],
    analysis_dt[spn_detected == 1L, .N], analysis_dt[spn_detected == 0L, .N],
    length(outcomes), length(model_specifications), length(scopes),
    nrow(patient_models),
    patient_models[
      scope == "all_retained_cd8" & outcome == "Exhaustion_core" &
        exposure == "spn_detected" &
        specification == "process_adjusted", uniqueN(subject_id)
    ],
    0L
  ))
)

fwrite(dataset_qc, file.path(results_dir, "30_dataset_qc.csv"))
fwrite(patient_models, file.path(results_dir, "30_patient_conditional_models.csv"))
fwrite(model_summary, file.path(results_dir, "30_conditional_model_summary.csv"))
fwrite(low_vif_summary, file.path(results_dir, "30_low_vif_sensitivity.csv"))
fwrite(attenuation, file.path(results_dir, "30_patient_effect_attenuation.csv"))
fwrite(attenuation_summary, file.path(results_dir, "30_effect_attenuation_summary.csv"))
fwrite(program_correlation_long,
       file.path(results_dir, "30_program_correlation_matrix.csv"))

theme_project <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.15)),
      plot.subtitle = element_text(color = "#5A5A5A"),
      axis.title = element_text(face = "bold"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

primary_plot <- model_summary[
  exposure == "spn_detected" & specification %chin% c(
    "basic", "activation_cytotoxicity_adjusted", "process_adjusted"
  )
]
primary_plot[, scope_label := factor(
  scope,
  levels = c("all_retained_cd8", "dysfunctional_only",
             "memory_transitional_only"),
  labels = c("All retained CD8", "Dysfunctional only",
             "Memory/transitional only")
)]
primary_plot[, specification_label := factor(
  specification,
  levels = c(
    "basic", "activation_cytotoxicity_adjusted", "process_adjusted"
  ),
  labels = c(
    "Depth/state adjusted", "Activation + cytotoxicity adjusted",
    "All biological processes adjusted"
  )
)]

p_forest <- ggplot(primary_plot, aes(
  x = median_spn_coefficient, y = scope_label,
  xmin = bootstrap_ci_lower, xmax = bootstrap_ci_upper,
  color = specification_label
)) +
  geom_vline(xintercept = 0, color = "grey60") +
  geom_errorbarh(height = 0.16, position = position_dodge(width = 0.42)) +
  geom_point(size = 2.6, position = position_dodge(width = 0.42)) +
  facet_wrap(~ outcome, scales = "free_x") +
  scale_color_manual(values = c("#3B82B8", "#D95F59", "#7651B9")) +
  labs(
    title = "Does SPN track exhaustion beyond activation and cytotoxicity?",
    subtitle = "Patient-specific coefficients; points are cross-patient medians with bootstrap 95% CIs",
    x = "Adjusted SPN-detection coefficient (outcome SD units)",
    y = NULL,
    color = NULL
  ) +
  theme_project()
ggsave(
  file.path(figures_dir, "30A_conditional_spn_exhaustion_models.png"),
  p_forest, width = 12, height = 7, dpi = 180, bg = "white"
)

heatmap_dt <- copy(program_correlation_long)
heatmap_dt[, variable_1 := sub("^step27_program_", "", variable_1)]
heatmap_dt[, variable_2 := sub("^step27_program_", "", variable_2)]
p_correlation <- ggplot(
  heatmap_dt, aes(variable_1, variable_2, fill = spearman_rho)
) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(
    low = "#3B82B8", mid = "white", high = "#7651B9",
    midpoint = 0, limits = c(-1, 1)
  ) +
  coord_equal() +
  labs(
    title = "Correlation among SPN and frozen CD8 programs",
    subtitle = "High correlations motivate the joint conditional models",
    x = NULL, y = NULL, fill = "Spearman rho"
  ) +
  theme_project(11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )
ggsave(
  file.path(figures_dir, "30B_program_correlation_heatmap.png"),
  p_correlation, width = 11, height = 9, dpi = 180, bg = "white"
)

saveRDS(
  list(
    dataset_qc = dataset_qc,
    patient_models = patient_models,
    model_summary = model_summary,
    low_vif_summary = low_vif_summary,
    attenuation = attenuation,
    attenuation_summary = attenuation_summary,
    program_correlation = program_correlation,
    model_specifications = model_specifications,
    adjusters_for_specification = lapply(
      model_specifications, adjusters_for_specification
    )
  ),
  output_rds_file,
  compress = FALSE
)

message("Step 30 complete: SPN-exhaustion specificity tested conditionally")