#!/usr/bin/env Rscript

# ==============================================================================
# Step 31: Map the CD8 differentiation continuum and clonal state sharing
# ==============================================================================
#
# Purpose
#   1. Resolve retained tumor CD8 cells into biologically anchored fine states
#      without using SPN.
#   2. Place cells on a cross-sectional memory/progenitor-to-terminal-
#      dysfunction continuum built from the frozen Step 27 programs.
#   3. Determine where SPN enrichment appears along that continuum using
#      patients as replicates.
#   4. Use productive TCR clonotypes to ask whether the same clones occupy
#      multiple fine states and whether clone expansion is associated with
#      state position or SPN.
#
# Important interpretation
#   The continuum is a transcriptional ordering, not temporal pseudotime.
#   Shared clonotypes support biological connectivity between states, but do
#   not establish direction or prove that one sampled cell became another.
#   SPN is not used to define programs, the continuum, or fine states.
# ==============================================================================

options(stringsAsFactors = FALSE)
options(error = function() {
  cat("\n--- STEP 31 TRACEBACK ---\n", file = stderr())
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
results_dir <- file.path(project_root, "results", "31")
figures_dir <- file.path(results_dir, "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

input_sce_file <- file.path(
  analysis_dir, "28_GSE123139_tumor_cd8_within_state_spn_biology_sce.rds"
)
output_sce_file <- file.path(
  analysis_dir, "31_GSE123139_cd8_continuum_clonal_states_sce.rds"
)

set.seed(123139L)
bootstrap_iterations <- 5000L
minimum_patient_bin_cells <- 10L
minimum_patient_state_cells <- 10L
minimum_clone_cells_expanded <- 3L
minimum_clone_cells_large <- 10L
fine_state_margin <- 0.20

base_states <- c("CD8_memory_transitional", "CD8_dysfunctional")
fine_state_levels <- c(
  "Memory/progenitor", "Activated/cytotoxic",
  "Terminal/dysfunctional", "Mixed/intermediate"
)
continuum_bin_levels <- paste0("Q", 1:5)

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

summarize_difference <- function(x, label, seed_offset = 0L) {
  x <- x[is.finite(x)]
  ci <- bootstrap_median_ci(x, seed_offset)
  data.table(
    comparison = label,
    n_independent_units = as.integer(length(x)),
    median_difference = if (length(x)) as.numeric(stats::median(x)) else NA_real_,
    bootstrap_ci_lower = unname(ci["lower"]),
    bootstrap_ci_upper = unname(ci["upper"]),
    percent_positive = if (length(x)) 100 * mean(x > 0) else NA_real_,
    wilcoxon_signed_rank_p = safe_wilcox_p(x)
  )
}

zscore_safe <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

entropy_safe <- function(counts) {
  counts <- counts[counts > 0]
  if (!length(counts)) return(NA_real_)
  p <- counts / sum(counts)
  -sum(p * log(p))
}

message("Loading the Step 28 program-annotated CD8 object")
if (!file.exists(input_sce_file)) stop("Missing Step 28 SCE: ", input_sce_file)
sce <- readRDS(input_sce_file)

programs_required <- c(
  "Memory_naive", "Progenitor_dysfunction", "TCR_activation",
  "Cytotoxic_effector", "Exhaustion_core", "Terminal_exhaustion",
  "Tissue_resident", "Interferon_response", "Proliferation"
)
program_columns <- paste0("step27_program_", programs_required)
required_metadata <- c(
  "subject_id", "step25_retained_cd8", "step25_base_state",
  "endogenous_n_counts", "endogenous_n_genes", "step26_spn_detected",
  "step26_spn_logcounts", "tcr_productive_primary", "tcr_clonotype_key",
  program_columns
)
missing_metadata <- setdiff(
  required_metadata, colnames(SummarizedExperiment::colData(sce))
)
if (length(missing_metadata)) {
  stop("Required metadata are missing: ", paste(missing_metadata, collapse = ", "))
}
if (!"UMAP_step21" %in% SingleCellExperiment::reducedDimNames(sce)) {
  stop("The Step 21 CD8 UMAP is missing.")
}

cell <- as.data.table(SummarizedExperiment::colData(sce))
cell[, cell_id := colnames(sce)]
cell[, subject_id := as.character(subject_id)]
cell[, step25_base_state := as.character(step25_base_state)]
cell[, retained_cd8 := !is.na(step25_retained_cd8) & as.logical(step25_retained_cd8)]
cell[, spn_detected := as.integer(
  !is.na(step26_spn_detected) & as.logical(step26_spn_detected)
)]
cell[, spn_logcounts := as.numeric(step26_spn_logcounts)]

analysis_dt <- cell[
  retained_cd8 == TRUE & step25_base_state %chin% base_states
]
if (!nrow(analysis_dt)) stop("No retained Step 25 CD8 cells were found.")

# ------------------------------------------------------------------------------
# Build axes and a cross-sectional continuum without SPN
# ------------------------------------------------------------------------------

message("Constructing patient-centered biological axes and continuum")
analysis_dt[, axis_memory_progenitor := rowMeans(cbind(
  step27_program_Memory_naive, step27_program_Progenitor_dysfunction
), na.rm = TRUE)]
analysis_dt[, axis_activated_cytotoxic := rowMeans(cbind(
  step27_program_TCR_activation, step27_program_Cytotoxic_effector
), na.rm = TRUE)]
analysis_dt[, axis_terminal_dysfunctional := rowMeans(cbind(
  step27_program_Exhaustion_core, step27_program_Terminal_exhaustion,
  step27_program_Tissue_resident
), na.rm = TRUE)]

axis_columns <- c(
  "axis_memory_progenitor", "axis_activated_cytotoxic",
  "axis_terminal_dysfunctional"
)
for (column in axis_columns) {
  centered_name <- paste0(column, "_patient_z")
  analysis_dt[, (centered_name) := zscore_safe(get(column)), by = subject_id]
}

analysis_dt[, step31_continuum_raw := (
  axis_terminal_dysfunctional_patient_z +
    0.50 * axis_activated_cytotoxic_patient_z -
    axis_memory_progenitor_patient_z
) / 2.50]
analysis_dt[, step31_continuum_percentile := {
  ranks <- frank(step31_continuum_raw, ties.method = "average")
  if (.N == 1L) 0.5 else (ranks - 1) / (.N - 1)
}, by = subject_id]
analysis_dt[, step31_continuum_bin := factor(
  paste0("Q", pmin(5L, floor(step31_continuum_percentile * 5) + 1L)),
  levels = continuum_bin_levels
)]

dominance_matrix <- as.matrix(analysis_dt[, .(
  axis_memory_progenitor_patient_z,
  axis_activated_cytotoxic_patient_z,
  axis_terminal_dysfunctional_patient_z
)])
dominant_index <- max.col(dominance_matrix, ties.method = "first")
sorted_axes <- t(apply(dominance_matrix, 1L, sort, decreasing = TRUE))
dominance_margin <- sorted_axes[, 1L] - sorted_axes[, 2L]
dominant_labels <- c(
  "Memory/progenitor", "Activated/cytotoxic", "Terminal/dysfunctional"
)[dominant_index]
dominant_labels[!is.finite(dominance_margin) |
                  dominance_margin < fine_state_margin] <- "Mixed/intermediate"
analysis_dt[, step31_fine_state := factor(
  dominant_labels, levels = fine_state_levels
)]
analysis_dt[, step31_state_margin := dominance_margin]

# Check that the oriented continuum behaves as intended across the inherited
# two-state solution. This is a validation, not a circular state definition.
patient_base_continuum <- analysis_dt[, .(
  n_cells = as.integer(.N),
  median_continuum = as.numeric(stats::median(step31_continuum_raw))
), by = .(subject_id, step25_base_state)][n_cells >= minimum_patient_state_cells]
base_continuum_wide <- dcast(
  patient_base_continuum, subject_id ~ step25_base_state,
  value.var = "median_continuum"
)
base_continuum_wide <- base_continuum_wide[
  is.finite(CD8_memory_transitional) & is.finite(CD8_dysfunctional)
]
base_continuum_wide[, difference :=
  CD8_dysfunctional - CD8_memory_transitional]
continuum_orientation_summary <- summarize_difference(
  base_continuum_wide$difference,
  "Dysfunctional minus memory/transitional continuum score", 100L
)

# ------------------------------------------------------------------------------
# Patient-balanced fine-state and continuum-bin summaries
# ------------------------------------------------------------------------------

message("Testing where SPN enrichment emerges along the continuum")
patient_bin <- analysis_dt[, .(
  n_cells = as.integer(.N),
  spn_detected_percent = 100 * mean(spn_detected),
  mean_spn_logcounts = mean(spn_logcounts),
  median_continuum = as.numeric(stats::median(step31_continuum_raw))
), by = .(subject_id, step31_continuum_bin)][
  n_cells >= minimum_patient_bin_cells
]

patient_bin_summary <- patient_bin[, .(
  n_patients = as.integer(.N),
  median_patient_spn_detected_percent = as.numeric(stats::median(spn_detected_percent)),
  q25_patient_spn_detected_percent = as.numeric(stats::quantile(spn_detected_percent, 0.25)),
  q75_patient_spn_detected_percent = as.numeric(stats::quantile(spn_detected_percent, 0.75)),
  median_patient_mean_spn_logcounts = as.numeric(stats::median(mean_spn_logcounts)),
  q25_patient_mean_spn_logcounts = as.numeric(stats::quantile(mean_spn_logcounts, 0.25)),
  q75_patient_mean_spn_logcounts = as.numeric(stats::quantile(mean_spn_logcounts, 0.75))
), by = step31_continuum_bin]

bin_detection_wide <- dcast(
  patient_bin, subject_id ~ step31_continuum_bin,
  value.var = "spn_detected_percent"
)
bin_expression_wide <- dcast(
  patient_bin, subject_id ~ step31_continuum_bin,
  value.var = "mean_spn_logcounts"
)
q5_q1_summary <- rbindlist(list(
  summarize_difference(
    bin_detection_wide[is.finite(Q1) & is.finite(Q5), Q5 - Q1],
    "Continuum Q5 minus Q1: SPN detection percentage points", 200L
  ),
  summarize_difference(
    bin_expression_wide[is.finite(Q1) & is.finite(Q5), Q5 - Q1],
    "Continuum Q5 minus Q1: mean SPN logcounts", 201L
  )
), fill = TRUE)
q5_q1_summary[, fdr := p.adjust(wilcoxon_signed_rank_p, method = "BH")]

patient_fine_state <- analysis_dt[, .(
  n_cells = as.integer(.N),
  spn_detected_percent = 100 * mean(spn_detected),
  mean_spn_logcounts = mean(spn_logcounts),
  median_continuum = as.numeric(stats::median(step31_continuum_raw))
), by = .(subject_id, step31_fine_state)][
  n_cells >= minimum_patient_state_cells
]
fine_state_summary <- patient_fine_state[, .(
  n_patients = as.integer(.N),
  median_patient_cells = as.numeric(stats::median(n_cells)),
  median_patient_spn_detected_percent = as.numeric(stats::median(spn_detected_percent)),
  q25_patient_spn_detected_percent = as.numeric(stats::quantile(spn_detected_percent, 0.25)),
  q75_patient_spn_detected_percent = as.numeric(stats::quantile(spn_detected_percent, 0.75)),
  median_patient_mean_spn_logcounts = as.numeric(stats::median(mean_spn_logcounts)),
  median_patient_continuum = as.numeric(stats::median(median_continuum))
), by = step31_fine_state]

fine_detection_wide <- dcast(
  patient_fine_state, subject_id ~ step31_fine_state,
  value.var = "spn_detected_percent"
)
fine_expression_wide <- dcast(
  patient_fine_state, subject_id ~ step31_fine_state,
  value.var = "mean_spn_logcounts"
)
terminal_memory_summary <- rbindlist(list(
  summarize_difference(
    fine_detection_wide[
      is.finite(`Memory/progenitor`) & is.finite(`Terminal/dysfunctional`),
      `Terminal/dysfunctional` - `Memory/progenitor`
    ],
    "Terminal/dysfunctional minus memory/progenitor: SPN detection percentage points",
    300L
  ),
  summarize_difference(
    fine_expression_wide[
      is.finite(`Memory/progenitor`) & is.finite(`Terminal/dysfunctional`),
      `Terminal/dysfunctional` - `Memory/progenitor`
    ],
    "Terminal/dysfunctional minus memory/progenitor: mean SPN logcounts",
    301L
  )
), fill = TRUE)
terminal_memory_summary[, fdr := p.adjust(wilcoxon_signed_rank_p, method = "BH")]

# ------------------------------------------------------------------------------
# Productive TCR clonotypes: expansion, state breadth, and matched contrasts
# ------------------------------------------------------------------------------

message("Quantifying clonal expansion and fine-state sharing")
# Normalize imported TCR metadata before character-only operations. Depending
# on how the original MARS-seq annotation was serialized, clonotype keys can
# arrive as integer, factor, or character vectors.
analysis_dt[, tcr_clonotype_key := as.character(tcr_clonotype_key)]
analysis_dt[, tcr_productive_flag := fifelse(
  is.na(tcr_productive_primary), FALSE,
  as.logical(tcr_productive_primary)
)]
clone_dt <- analysis_dt[
  tcr_productive_flag == TRUE &
    !is.na(tcr_clonotype_key) & nzchar(tcr_clonotype_key)
]
clone_dt[, clone_id := paste(subject_id, tcr_clonotype_key, sep = "::")]
# `%chin%` is intentionally character-only; factors are integer vectors
# internally and therefore must be normalized before clone-state matching.
clone_dt[, step31_fine_state := as.character(step31_fine_state)]
message("  Productive TCR cells selected: ", nrow(clone_dt))

clone_total <- clone_dt[, .(
  clone_size = as.integer(.N),
  mean_continuum = mean(step31_continuum_raw),
  median_continuum = as.numeric(stats::median(step31_continuum_raw)),
  spn_detected_percent = 100 * mean(spn_detected),
  mean_spn_logcounts = mean(spn_logcounts),
  n_fine_states = as.integer(uniqueN(step31_fine_state)),
  fine_state_entropy = entropy_safe(table(step31_fine_state)),
  memory_progenitor_cells = as.integer(sum(step31_fine_state == "Memory/progenitor")),
  activated_cytotoxic_cells = as.integer(sum(step31_fine_state == "Activated/cytotoxic")),
  terminal_dysfunctional_cells = as.integer(sum(step31_fine_state == "Terminal/dysfunctional")),
  mixed_intermediate_cells = as.integer(sum(step31_fine_state == "Mixed/intermediate"))
), by = .(clone_id, subject_id)]
message("  Clone-level summaries calculated: ", nrow(clone_total))
clone_total[, expansion_class := fifelse(
  clone_size >= minimum_clone_cells_large, "Large (>=10 cells)",
  fifelse(clone_size >= minimum_clone_cells_expanded,
          "Expanded (3-9 cells)", "Small (1-2 cells)")
)]

clone_state <- clone_dt[, .(
  n_cells = as.integer(.N),
  spn_detected_percent = 100 * mean(spn_detected),
  mean_spn_logcounts = mean(spn_logcounts),
  mean_continuum = mean(step31_continuum_raw)
), by = .(clone_id, subject_id, step31_fine_state)]
message("  Clone-by-state summaries calculated: ", nrow(clone_state))

clone_sharing_summary <- clone_total[, .(
  productive_clonotypes = as.integer(.N),
  productive_cells = as.integer(sum(clone_size)),
  median_clone_size = as.numeric(stats::median(clone_size)),
  clones_spanning_multiple_fine_states = as.integer(sum(n_fine_states >= 2L)),
  percent_clones_spanning_multiple_fine_states = 100 * mean(n_fine_states >= 2L),
  clones_spanning_memory_and_terminal = as.integer(sum(
    memory_progenitor_cells >= 1L & terminal_dysfunctional_cells >= 1L
  )),
  clones_spanning_memory_and_terminal_min_two_each = as.integer(sum(
    memory_progenitor_cells >= 2L & terminal_dysfunctional_cells >= 2L
  ))
)]
message("  Clone state-sharing summary calculated")

build_matched_clone_pairs <- function(minimum_per_state) {
  eligible <- clone_total[
    memory_progenitor_cells >= minimum_per_state &
      terminal_dysfunctional_cells >= minimum_per_state,
    .(clone_id, subject_id, clone_size)
  ]
  if (!nrow(eligible)) return(data.table())
  selected <- clone_state[
    step31_fine_state %chin% c("Memory/progenitor", "Terminal/dysfunctional")
  ][eligible, on = .(clone_id, subject_id), nomatch = 0L]

  # Aggregate directly instead of using dcast(). Fine-state labels contain
  # slashes, and some data.table/R combinations attempt to treat their factor
  # codes as character column names during reshaping.
  output <- selected[, .(
    memory_detection = spn_detected_percent[
      step31_fine_state == "Memory/progenitor"
    ][1L],
    terminal_detection = spn_detected_percent[
      step31_fine_state == "Terminal/dysfunctional"
    ][1L],
    memory_expression = mean_spn_logcounts[
      step31_fine_state == "Memory/progenitor"
    ][1L],
    terminal_expression = mean_spn_logcounts[
      step31_fine_state == "Terminal/dysfunctional"
    ][1L]
  ), by = .(clone_id, subject_id, clone_size)]
  output[, spn_detection_difference_terminal_minus_memory :=
    terminal_detection - memory_detection]
  output[, spn_expression_difference_terminal_minus_memory :=
    terminal_expression - memory_expression]
  output[, minimum_cells_per_state := as.integer(minimum_per_state)]
  output
}

clone_pairs <- rbindlist(list(
  build_matched_clone_pairs(1L), build_matched_clone_pairs(2L)
), fill = TRUE)
message("  Matched memory-terminal clone pairs calculated: ", nrow(clone_pairs))

clone_pair_summary <- if (nrow(clone_pairs)) {
  rbindlist(lapply(sort(unique(clone_pairs$minimum_cells_per_state)), function(k) {
    current <- clone_pairs[minimum_cells_per_state == k]
    rbindlist(list(
      summarize_difference(
        current$spn_detection_difference_terminal_minus_memory,
        paste0("Terminal minus memory within productive clones; >=", k,
               " cell(s) per state: SPN detection percentage points"),
        400L + k
      ),
      summarize_difference(
        current$spn_expression_difference_terminal_minus_memory,
        paste0("Terminal minus memory within productive clones; >=", k,
               " cell(s) per state: mean SPN logcounts"),
        500L + k
      )
    ))[, minimum_cells_per_state := as.integer(k)]
  }), fill = TRUE)
} else data.table()
if (nrow(clone_pair_summary)) {
  clone_pair_summary[, fdr := p.adjust(wilcoxon_signed_rank_p, method = "BH")]
}

# Compare expanded-clone cells with cells in small productive clones within
# each patient. This is descriptive of clonal expansion, not antigen specificity.
clone_dt <- clone_dt[clone_total[, .(clone_id, expansion_class)], on = "clone_id"]
message("  Expansion classes joined to productive TCR cells")
clone_dt[, expanded_clone := expansion_class != "Small (1-2 cells)"]
patient_expansion <- clone_dt[, .(
  n_cells = as.integer(.N),
  median_continuum = as.numeric(stats::median(step31_continuum_raw)),
  spn_detected_percent = 100 * mean(spn_detected),
  mean_spn_logcounts = mean(spn_logcounts)
), by = .(subject_id, expanded_clone)][n_cells >= minimum_patient_state_cells]

exp_continuum <- dcast(
  patient_expansion, subject_id ~ expanded_clone, value.var = "median_continuum"
)
exp_detection <- dcast(
  patient_expansion, subject_id ~ expanded_clone, value.var = "spn_detected_percent"
)
exp_expression <- dcast(
  patient_expansion, subject_id ~ expanded_clone, value.var = "mean_spn_logcounts"
)
expansion_summary <- rbindlist(list(
  summarize_difference(
    exp_continuum[is.finite(`TRUE`) & is.finite(`FALSE`), `TRUE` - `FALSE`],
    "Expanded minus small productive clones: median continuum", 600L
  ),
  summarize_difference(
    exp_detection[is.finite(`TRUE`) & is.finite(`FALSE`), `TRUE` - `FALSE`],
    "Expanded minus small productive clones: SPN detection percentage points", 601L
  ),
  summarize_difference(
    exp_expression[is.finite(`TRUE`) & is.finite(`FALSE`), `TRUE` - `FALSE`],
    "Expanded minus small productive clones: mean SPN logcounts", 602L
  )
), fill = TRUE)
expansion_summary[, fdr := p.adjust(wilcoxon_signed_rank_p, method = "BH")]

# ------------------------------------------------------------------------------
# Attach annotations, validate, and write outputs
# ------------------------------------------------------------------------------

message("Writing Step 31 tables and annotated object")
match_index <- match(colnames(sce), analysis_dt$cell_id)
SummarizedExperiment::colData(sce)$step31_continuum_raw <-
  analysis_dt$step31_continuum_raw[match_index]
SummarizedExperiment::colData(sce)$step31_continuum_percentile <-
  analysis_dt$step31_continuum_percentile[match_index]
SummarizedExperiment::colData(sce)$step31_continuum_bin <-
  as.character(analysis_dt$step31_continuum_bin[match_index])
SummarizedExperiment::colData(sce)$step31_fine_state <-
  as.character(analysis_dt$step31_fine_state[match_index])
SummarizedExperiment::colData(sce)$step31_state_margin <-
  analysis_dt$step31_state_margin[match_index]

state_accounting_validation <- as.integer(
  sum(table(analysis_dt$step31_fine_state)) == nrow(analysis_dt)
)
continuum_complete_validation <- as.integer(
  all(is.finite(analysis_dt$step31_continuum_raw))
)
spn_independent_definition_validation <- 1L

dataset_qc <- data.table(
  metric = c(
    "retained_cd8_cells", "subjects", "fine_states", "continuum_bins",
    "productive_tcr_cells", "productive_clonotypes", "expanded_clonotypes",
    "large_clonotypes", "clones_spanning_multiple_fine_states",
    "clones_spanning_memory_and_terminal",
    "clones_spanning_memory_and_terminal_min_two_each",
    "state_accounting_validation", "continuum_complete_validation",
    "spn_used_to_define_continuum_or_fine_states", "signature_exported"
  ),
  value = c(
    nrow(analysis_dt), uniqueN(analysis_dt$subject_id),
    uniqueN(analysis_dt$step31_fine_state),
    uniqueN(analysis_dt$step31_continuum_bin), nrow(clone_dt),
    nrow(clone_total),
    clone_total[clone_size >= minimum_clone_cells_expanded, .N],
    clone_total[clone_size >= minimum_clone_cells_large, .N],
    clone_total[n_fine_states >= 2L, .N],
    clone_total[memory_progenitor_cells >= 1L &
                  terminal_dysfunctional_cells >= 1L, .N],
    clone_total[memory_progenitor_cells >= 2L &
                  terminal_dysfunctional_cells >= 2L, .N],
    state_accounting_validation, continuum_complete_validation,
    1L - spn_independent_definition_validation, 0L
  )
)

fwrite(dataset_qc, file.path(results_dir, "31_dataset_qc.csv"))
fwrite(continuum_orientation_summary,
       file.path(results_dir, "31_continuum_orientation_validation.csv"))
fwrite(patient_bin, file.path(results_dir, "31_patient_continuum_bin_spn.csv"))
fwrite(patient_bin_summary,
       file.path(results_dir, "31_continuum_bin_spn_summary.csv"))
fwrite(q5_q1_summary,
       file.path(results_dir, "31_continuum_extreme_contrast_summary.csv"))
fwrite(patient_fine_state,
       file.path(results_dir, "31_patient_fine_state_summary.csv"))
fwrite(fine_state_summary,
       file.path(results_dir, "31_fine_state_spn_summary.csv"))
fwrite(terminal_memory_summary,
       file.path(results_dir, "31_terminal_memory_spn_contrast.csv"))
fwrite(clone_total, file.path(results_dir, "31_productive_clonotype_summary.csv"))
fwrite(clone_state, file.path(results_dir, "31_productive_clone_fine_states.csv"))
fwrite(clone_sharing_summary,
       file.path(results_dir, "31_clone_state_sharing_summary.csv"))
fwrite(clone_pairs, file.path(results_dir, "31_clone_matched_state_differences.csv"))
fwrite(clone_pair_summary,
       file.path(results_dir, "31_clone_matched_state_summary.csv"))
fwrite(patient_expansion,
       file.path(results_dir, "31_patient_clone_expansion_summary.csv"))
fwrite(expansion_summary,
       file.path(results_dir, "31_clone_expansion_association_summary.csv"))

# ------------------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------------------

theme_step31 <- theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(color = "grey35"),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

umap <- as.data.table(SingleCellExperiment::reducedDim(sce, "UMAP_step21"))
setnames(umap, c("UMAP1", "UMAP2"))
umap[, cell_id := colnames(sce)]
plot_dt <- umap[analysis_dt[, .(
  cell_id, step31_fine_state, step31_continuum_percentile, spn_detected
)], on = "cell_id", nomatch = 0L]

p_states <- ggplot(plot_dt, aes(UMAP1, UMAP2, color = step31_fine_state)) +
  geom_point(size = 0.45, alpha = 0.65) +
  scale_color_manual(values = c(
    "Memory/progenitor" = "#3B82B8", "Activated/cytotoxic" = "#E36A5C",
    "Terminal/dysfunctional" = "#7651B8", "Mixed/intermediate" = "#A7A7A7"
  )) +
  labs(
    title = "Program-defined fine states across retained tumor CD8 cells",
    subtitle = "States and continuum were defined without SPN",
    color = "Fine state"
  ) + theme_step31
ggsave(file.path(figures_dir, "31A_fine_states_umap.png"), p_states,
       width = 11, height = 7, dpi = 300)

p_continuum <- ggplot(plot_dt, aes(UMAP1, UMAP2, color = step31_continuum_percentile)) +
  geom_point(size = 0.45, alpha = 0.7) +
  scale_color_gradientn(
    colors = c("#2C7BB6", "#FFFFBF", "#D7191C"), limits = c(0, 1)
  ) +
  labs(
    title = "Cross-sectional CD8 transcriptional continuum",
    subtitle = "Memory/progenitor to activated and terminally dysfunctional",
    color = "Within-patient\ncontinuum percentile"
  ) + theme_step31
ggsave(file.path(figures_dir, "31B_continuum_umap.png"), p_continuum,
       width = 11, height = 7, dpi = 300)

p_spn_trend <- ggplot(
  patient_bin_summary,
  aes(step31_continuum_bin, median_patient_spn_detected_percent, group = 1)
) +
  geom_line(color = "#7651B8", linewidth = 0.9) +
  geom_point(color = "#7651B8", size = 3) +
  geom_errorbar(aes(
    ymin = q25_patient_spn_detected_percent,
    ymax = q75_patient_spn_detected_percent
  ), width = 0.12, color = "#7651B8") +
  labs(
    title = "Patient-balanced SPN detection along the CD8 continuum",
    subtitle = "Points are patient medians; bars are interquartile ranges",
    x = "Within-patient continuum quintile",
    y = "SPN-detected cells (%)"
  ) + theme_step31
ggsave(file.path(figures_dir, "31C_spn_along_continuum.png"), p_spn_trend,
       width = 9, height = 6, dpi = 300)

clone_breadth <- clone_total[, .N, by = .(expansion_class, n_fine_states)]
p_clone_breadth <- ggplot(
  clone_breadth,
  aes(factor(n_fine_states), N, fill = expansion_class)
) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c(
    "Small (1-2 cells)" = "#BDBDBD", "Expanded (3-9 cells)" = "#3B82B8",
    "Large (>=10 cells)" = "#7651B8"
  )) +
  labs(
    title = "Expanded productive clonotypes span more CD8 fine states",
    subtitle = "State sharing is evidence of clonal state heterogeneity, not directionality",
    x = "Number of fine states occupied", y = "Productive clonotypes",
    fill = "Clone size"
  ) + theme_step31
ggsave(file.path(figures_dir, "31D_clone_state_breadth.png"), p_clone_breadth,
       width = 9, height = 6, dpi = 300)

if (nrow(clone_pairs[minimum_cells_per_state == 1L])) {
  p_clone_spn <- ggplot(
    clone_pairs[minimum_cells_per_state == 1L],
    aes(spn_detection_difference_terminal_minus_memory)
  ) +
    geom_histogram(binwidth = 10, boundary = 0, fill = "#7651B8", color = "white") +
    geom_vline(xintercept = 0, color = "grey35") +
    labs(
      title = "SPN differences within shared productive clonotypes",
      subtitle = "Positive values indicate greater detection in terminal/dysfunctional cells",
      x = "Terminal minus memory/progenitor SPN detection (percentage points)",
      y = "Productive clonotypes"
    ) + theme_step31
  ggsave(file.path(figures_dir, "31E_clone_matched_spn_difference.png"),
         p_clone_spn, width = 9, height = 6, dpi = 300)
}

saveRDS(sce, output_sce_file, compress = FALSE)

message("Step 31 complete: CD8 continuum, fine states, and clonal sharing quantified")
