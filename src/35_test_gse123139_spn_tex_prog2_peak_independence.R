#!/usr/bin/env Rscript

# ==============================================================================
# Step 35: Is the SPN peak in Beltra Tex_prog2 independent of depth/proliferation?
# ==============================================================================
#
# Step 34 projected four SPN-independent exhaustion states and found the highest
# raw SPN-detection percentage in the cycling Tex_prog2 state. Because cycling
# cells also contain more RNA, this script performs a focused patient-aware
# sensitivity analysis.
#
# For each adjacent transition and for Tex_prog2 versus all other Tex states:
#   1. fit within-patient SPN-detection models adjusted for RNA depth;
#   2. add the frozen Step 27 proliferation score;
#   3. repeat for SPN logcounts;
#   4. convert detection models to standardized probability differences using
#      g-computation, while preserving each patient's observed covariates.
#
# SPN is not used to define the states or the proliferation program.
# Patients, rather than cells, are the inferential units.
# ============================================================================== 

options(stringsAsFactors = FALSE)
options(error = function() {
  cat("\n--- STEP 35 TRACEBACK ---\n", file = stderr())
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
results_dir <- file.path(project_root, "results", "35")
figures_dir <- file.path(results_dir, "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

input_sce_file <- file.path(
  analysis_dir, "34_GSE123139_beltra_tex_states_spn_sce.rds"
)
if (!file.exists(input_sce_file)) stop("Step 34 object not found: ", input_sce_file)

random_seed <- 123139L
bootstrap_iterations <- 5000L
minimum_cells_per_state <- 10L
minimum_total_pair_cells <- 40L
minimum_outcome_class_cells <- 3L

tex_states <- c("Tex_prog1", "Tex_prog2", "Tex_int", "Tex_term")
comparisons <- list(
  Prog2_minus_prog1 = c("Tex_prog2", "Tex_prog1"),
  Int_minus_prog2 = c("Tex_int", "Tex_prog2"),
  Term_minus_int = c("Tex_term", "Tex_int"),
  Prog2_minus_other_Tex = c("Tex_prog2", "Other_Tex")
)

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

summarize_effect <- function(x, comparison, model, metric, effect_scale,
                             seed_offset = 0L) {
  x <- x[is.finite(x)]
  ci <- bootstrap_median_ci(x, seed_offset)
  data.table(
    comparison = comparison,
    model = model,
    metric = metric,
    effect_scale = effect_scale,
    n_patients = as.integer(length(x)),
    median_effect = if (length(x)) as.numeric(stats::median(x)) else NA_real_,
    bootstrap_ci_lower = unname(ci["lower"]),
    bootstrap_ci_upper = unname(ci["upper"]),
    percent_positive = if (length(x)) 100 * mean(x > 0) else NA_real_,
    wilcoxon_signed_rank_p = safe_wilcox_p(x)
  )
}

theme_project <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.18)),
      plot.subtitle = element_text(color = "#5A5A5A"),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )
}

message("Loading completed Step 34 object")
sce <- readRDS(input_sce_file)

required_metadata <- c(
  "subject_id", "step25_retained_cd8", "step34_consensus_state",
  "endogenous_n_counts", "step26_spn_detected", "step26_spn_logcounts",
  "step27_program_Proliferation"
)
missing_metadata <- setdiff(
  required_metadata, colnames(SummarizedExperiment::colData(sce))
)
if (length(missing_metadata)) {
  stop("Required metadata are missing: ", paste(missing_metadata, collapse = ", "))
}

cell <- as.data.table(SummarizedExperiment::colData(sce))
cell[, `:=`(
  subject_id = as.character(subject_id),
  retained_cd8 = !is.na(step25_retained_cd8) & as.logical(step25_retained_cd8),
  tex_state = as.character(step34_consensus_state),
  spn_detected = as.integer(
    !is.na(step26_spn_detected) & as.logical(step26_spn_detected)
  ),
  spn_logcounts = as.numeric(step26_spn_logcounts),
  log1p_endogenous_umis = log1p(as.numeric(endogenous_n_counts)),
  proliferation_score = as.numeric(step27_program_Proliferation)
)]

analysis_dt <- cell[
  retained_cd8 == TRUE & tex_state %chin% tex_states &
    is.finite(log1p_endogenous_umis) & is.finite(proliferation_score) &
    is.finite(spn_logcounts)
]
if (!nrow(analysis_dt)) stop("No consensus Tex cells were available.")

fit_pair_patient <- function(patient_dt, high_state, low_state, comparison_name) {
  if (low_state == "Other_Tex") {
    patient_dt <- patient_dt[tex_state %chin% tex_states]
    patient_dt[, high__ := as.integer(tex_state == high_state)]
    n_high <- sum(patient_dt$high__ == 1L)
    n_low <- sum(patient_dt$high__ == 0L)
  } else {
    patient_dt <- patient_dt[tex_state %chin% c(high_state, low_state)]
    patient_dt[, high__ := as.integer(tex_state == high_state)]
    n_high <- sum(patient_dt$tex_state == high_state)
    n_low <- sum(patient_dt$tex_state == low_state)
  }

  empty <- data.table(
    comparison = comparison_name,
    high_state = high_state,
    low_state = low_state,
    n_cells = as.integer(nrow(patient_dt)),
    n_high = as.integer(n_high), n_low = as.integer(n_low),
    raw_detection_difference_pp = NA_real_,
    depth_adjusted_detection_log_odds = NA_real_,
    depth_adjusted_detection_difference_pp = NA_real_,
    full_adjusted_detection_log_odds = NA_real_,
    full_adjusted_detection_difference_pp = NA_real_,
    depth_adjusted_expression_coefficient = NA_real_,
    full_adjusted_expression_coefficient = NA_real_
  )

  if (nrow(patient_dt) < minimum_total_pair_cells ||
      n_high < minimum_cells_per_state || n_low < minimum_cells_per_state ||
      sum(patient_dt$spn_detected == 1L) < minimum_outcome_class_cells ||
      sum(patient_dt$spn_detected == 0L) < minimum_outcome_class_cells) {
    return(empty)
  }

  patient_dt[, depth_z__ := zscore_safe(log1p_endogenous_umis)]
  patient_dt[, proliferation_z__ := zscore_safe(proliferation_score)]

  detection_depth <- tryCatch(
    suppressWarnings(glm(
      spn_detected ~ high__ + depth_z__, data = patient_dt,
      family = quasibinomial()
    )), error = function(e) NULL
  )
  detection_full <- tryCatch(
    suppressWarnings(glm(
      spn_detected ~ high__ + depth_z__ + proliferation_z__,
      data = patient_dt, family = quasibinomial()
    )), error = function(e) NULL
  )
  expression_depth <- tryCatch(
    lm(spn_logcounts ~ high__ + depth_z__, data = patient_dt),
    error = function(e) NULL
  )
  expression_full <- tryCatch(
    lm(
      spn_logcounts ~ high__ + depth_z__ + proliferation_z__,
      data = patient_dt
    ), error = function(e) NULL
  )

  standardized_difference <- function(fit) {
    if (is.null(fit)) return(NA_real_)
    new_high <- copy(patient_dt)
    new_low <- copy(patient_dt)
    new_high[, high__ := 1L]
    new_low[, high__ := 0L]
    prediction_high <- tryCatch(
      predict(fit, newdata = new_high, type = "response"),
      error = function(e) rep(NA_real_, nrow(new_high))
    )
    prediction_low <- tryCatch(
      predict(fit, newdata = new_low, type = "response"),
      error = function(e) rep(NA_real_, nrow(new_low))
    )
    if (!all(is.finite(c(prediction_high, prediction_low)))) return(NA_real_)
    100 * mean(prediction_high - prediction_low)
  }

  coefficient <- function(fit) {
    if (is.null(fit) || !"high__" %in% names(stats::coef(fit))) return(NA_real_)
    as.numeric(stats::coef(fit)["high__"])
  }

  data.table(
    comparison = comparison_name,
    high_state = high_state,
    low_state = low_state,
    n_cells = as.integer(nrow(patient_dt)),
    n_high = as.integer(n_high), n_low = as.integer(n_low),
    raw_detection_difference_pp = 100 * (
      mean(patient_dt[high__ == 1L]$spn_detected) -
        mean(patient_dt[high__ == 0L]$spn_detected)
    ),
    depth_adjusted_detection_log_odds = coefficient(detection_depth),
    depth_adjusted_detection_difference_pp = standardized_difference(
      detection_depth
    ),
    full_adjusted_detection_log_odds = coefficient(detection_full),
    full_adjusted_detection_difference_pp = standardized_difference(
      detection_full
    ),
    depth_adjusted_expression_coefficient = coefficient(expression_depth),
    full_adjusted_expression_coefficient = coefficient(expression_full)
  )
}

message("Fitting patient-specific adjacent-state and Prog2-peak models")
patient_rows <- list()
for (comparison_name in names(comparisons)) {
  pair <- comparisons[[comparison_name]]
  temp <- analysis_dt[, fit_pair_patient(
    copy(.SD), pair[1L], pair[2L], comparison_name
  ), by = subject_id]
  patient_rows[[comparison_name]] <- temp
}
patient_models <- rbindlist(patient_rows, use.names = TRUE, fill = TRUE)

summary_specs <- list(
  list("Raw", "SPN detection", "raw_detection_difference_pp",
       "SPN-detection percentage points"),
  list("Depth adjusted", "SPN detection",
       "depth_adjusted_detection_difference_pp",
       "Standardized SPN-detection percentage points"),
  list("Depth + proliferation adjusted", "SPN detection",
       "full_adjusted_detection_difference_pp",
       "Standardized SPN-detection percentage points"),
  list("Depth adjusted", "SPN detection log-odds",
       "depth_adjusted_detection_log_odds", "SPN-detection log-odds"),
  list("Depth + proliferation adjusted", "SPN detection log-odds",
       "full_adjusted_detection_log_odds", "SPN-detection log-odds"),
  list("Depth adjusted", "SPN expression",
       "depth_adjusted_expression_coefficient", "SPN-logcount coefficient"),
  list("Depth + proliferation adjusted", "SPN expression",
       "full_adjusted_expression_coefficient", "SPN-logcount coefficient")
)

summary_rows <- list()
seed_counter <- 0L
for (comparison_name in names(comparisons)) {
  comparison_dt <- patient_models[comparison == comparison_name]
  for (spec in summary_specs) {
    seed_counter <- seed_counter + 1L
    summary_rows[[length(summary_rows) + 1L]] <- summarize_effect(
      comparison_dt[[spec[[3L]]]], comparison_name, spec[[1L]], spec[[2L]],
      spec[[4L]], seed_counter
    )
  }
}
adjusted_summary <- rbindlist(summary_rows)
adjusted_summary[, fdr := p.adjust(wilcoxon_signed_rank_p, method = "BH"),
                 by = .(metric, model)]

# A compact result focused on the biological question motivating Step 35.
prog2_peak_summary <- adjusted_summary[
  comparison %chin% c("Prog2_minus_prog1", "Int_minus_prog2",
                        "Prog2_minus_other_Tex") &
    metric == "SPN detection" &
    model %chin% c("Raw", "Depth adjusted",
                     "Depth + proliferation adjusted")
]

get_evidence <- function(comparison_name) {
  row <- adjusted_summary[
    comparison == comparison_name & metric == "SPN detection" &
      model == "Depth + proliferation adjusted"
  ]
  if (!nrow(row) || !is.finite(row$fdr)) return("Insufficient evidence")
  if (row$fdr < 0.05 && row$median_effect > 0) return("Supported")
  if (row$fdr < 0.05 && row$median_effect < 0) return("Supported opposite")
  "Not supported"
}

evidence_ledger <- rbindlist(lapply(names(comparisons), function(comparison_name) {
  row <- adjusted_summary[
    comparison == comparison_name & metric == "SPN detection" &
      model == "Depth + proliferation adjusted"
  ]
  data.table(
    hypothesis = switch(
      comparison_name,
      Prog2_minus_prog1 = "SPN detection is higher in Tex_prog2 than Tex_prog1 independently of RNA depth and proliferation",
      Int_minus_prog2 = "SPN detection continues to increase from Tex_prog2 to Tex_int independently of RNA depth and proliferation",
      Term_minus_int = "SPN detection increases from Tex_int to Tex_term independently of RNA depth and proliferation",
      Prog2_minus_other_Tex = "The apparent SPN peak is specific to Tex_prog2 independently of RNA depth and proliferation"
    ),
    evidence_class = get_evidence(comparison_name),
    comparison = comparison_name,
    n_patients = if (nrow(row)) row$n_patients else 0L,
    median_adjusted_difference_pp = if (nrow(row)) row$median_effect else NA_real_,
    bootstrap_ci_lower = if (nrow(row)) row$bootstrap_ci_lower else NA_real_,
    bootstrap_ci_upper = if (nrow(row)) row$bootstrap_ci_upper else NA_real_,
    percent_positive = if (nrow(row)) row$percent_positive else NA_real_,
    fdr = if (nrow(row)) row$fdr else NA_real_
  )
}))

dataset_qc <- data.table(
  metric = c(
    "consensus_tex_cells", "subjects", "comparisons_tested",
    "patients_with_prog2_prog1_models", "patients_with_int_prog2_models",
    "patients_with_term_int_models", "patients_with_prog2_other_models",
    "spn_used_to_define_states_or_proliferation", "signature_exported"
  ),
  value = c(
    nrow(analysis_dt), uniqueN(analysis_dt$subject_id), length(comparisons),
    sum(is.finite(patient_models[
      comparison == "Prog2_minus_prog1"
    ]$full_adjusted_detection_difference_pp)),
    sum(is.finite(patient_models[
      comparison == "Int_minus_prog2"
    ]$full_adjusted_detection_difference_pp)),
    sum(is.finite(patient_models[
      comparison == "Term_minus_int"
    ]$full_adjusted_detection_difference_pp)),
    sum(is.finite(patient_models[
      comparison == "Prog2_minus_other_Tex"
    ]$full_adjusted_detection_difference_pp)),
    0L, 0L
  )
)

message("Writing Step 35 tables")
fwrite(dataset_qc, file.path(results_dir, "35_dataset_qc.csv"))
fwrite(patient_models,
       file.path(results_dir, "35_patient_pair_models.csv"))
fwrite(adjusted_summary,
       file.path(results_dir, "35_pairwise_adjusted_summary.csv"))
fwrite(prog2_peak_summary,
       file.path(results_dir, "35_prog2_peak_summary.csv"))
fwrite(evidence_ledger,
       file.path(results_dir, "35_evidence_ledger.csv"))

message("Rendering Step 35 figures")
plot_dt <- adjusted_summary[
  metric == "SPN detection" &
    model %chin% c("Raw", "Depth adjusted",
                     "Depth + proliferation adjusted")
]
plot_dt[, comparison := factor(
  comparison,
  levels = c("Prog2_minus_prog1", "Int_minus_prog2", "Term_minus_int",
             "Prog2_minus_other_Tex")
)]
plot_dt[, model := factor(
  model,
  levels = c("Raw", "Depth adjusted", "Depth + proliferation adjusted")
)]

p_effect <- ggplot(
  plot_dt,
  aes(median_effect, comparison, color = model)
) +
  geom_vline(xintercept = 0, color = "#888888") +
  geom_errorbarh(
    aes(xmin = bootstrap_ci_lower, xmax = bootstrap_ci_upper),
    position = position_dodge(width = 0.55), height = 0.18
  ) +
  geom_point(position = position_dodge(width = 0.55), size = 2.7) +
  scale_color_manual(values = c(
    "Raw" = "#777777", "Depth adjusted" = "#3686B9",
    "Depth + proliferation adjusted" = "#7A4EBF"
  )) +
  labs(
    title = "Is the apparent SPN peak in Tex_prog2 independent?",
    subtitle = "Patient-specific effects summarized by cross-patient medians and bootstrap 95% CIs",
    x = "SPN-detection difference (percentage points)", y = NULL,
    color = "Model"
  ) + theme_project()
ggsave(
  file.path(figures_dir, "35A_adjusted_spn_tex_state_effects.png"),
  p_effect, width = 12, height = 7.5, units = "in", dpi = 180, bg = "white"
)

patient_plot <- patient_models[
  comparison == "Prog2_minus_other_Tex" &
    is.finite(full_adjusted_detection_difference_pp)
]
p_patient <- ggplot(
  patient_plot,
  aes(reorder(subject_id, full_adjusted_detection_difference_pp),
      full_adjusted_detection_difference_pp)
) +
  geom_hline(yintercept = 0, color = "#888888") +
  geom_col(fill = "#7A4EBF", width = 0.72) +
  coord_flip() +
  labs(
    title = "Adjusted Tex_prog2 SPN peak by patient",
    subtitle = "Tex_prog2 minus other Tex states; RNA depth and proliferation adjusted",
    x = "Patient", y = "Standardized SPN-detection difference (percentage points)"
  ) + theme_project()
ggsave(
  file.path(figures_dir, "35B_patient_adjusted_prog2_peak.png"),
  p_patient, width = 9, height = 7.5, units = "in", dpi = 180, bg = "white"
)

message("Step 35 complete: Tex_prog2 SPN peak tested after depth/proliferation adjustment")
