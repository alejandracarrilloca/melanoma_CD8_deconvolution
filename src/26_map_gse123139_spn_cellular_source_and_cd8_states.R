#!/usr/bin/env Rscript

# ==============================================================================
# Step 26: SPN/CD43 cellular source and tumor-CD8 state associations
# ==============================================================================
#
# This is the first single-cell-first biological analysis after the Step 25
# lineage gate. It answers two questions without using SPN to define any cell
# population or CD8 state:
#
#   1. Which broad cellular compartments in GSE123139 express SPN?
#   2. Within conservatively retained tumor CD8 cells, is SPN associated with
#      memory/transitional versus dysfunctional identity or with interferon,
#      cytotoxic, and proliferation overlays?
#
# Inference is patient-balanced. Cell-level models are fit separately within
# each patient and summarized across patients; cells are never treated as
# independent biological replicates. Productive TCR clonotypes spanning both
# base states provide a separate matched sensitivity analysis.
#
# This step does not build or export a deconvolution signature and does not use
# the label "exhausted" as a synonym for the Step 25 dysfunctional state.
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
# Paths and analysis constants
# ------------------------------------------------------------------------------

project_root <- normalizePath(getwd(), mustWork = TRUE)
analysis_dir <- file.path(
  project_root, "data", "analysis_ready", "scRNA", "GSE123139"
)
results_dir <- file.path(project_root, "results", "26")
figures_dir <- file.path(results_dir, "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

full_sce_file <- file.path(
  analysis_dir, "19_GSE123139_broad_clusters_sce.rds"
)
cd8_sce_file <- file.path(
  analysis_dir, "25_GSE123139_tumor_cd8_lineage_adjudicated_sce.rds"
)
output_cd8_sce_file <- file.path(
  analysis_dir, "26_GSE123139_tumor_cd8_spn_annotated_sce.rds"
)

random_seed <- 123139L
bootstrap_iterations <- 5000L
minimum_cells_per_patient_group <- 10L
minimum_overlay_positive_cells <- 5L
minimum_clone_cells_per_state_primary <- 1L
minimum_clone_cells_per_state_sensitivity <- 2L

base_state_levels <- c("CD8_memory_transitional", "CD8_dysfunctional")
base_state_labels <- c(
  CD8_memory_transitional = "Memory/transitional",
  CD8_dysfunctional = "Dysfunctional"
)
base_state_colors <- c(
  CD8_memory_transitional = "#3B82B8",
  CD8_dysfunctional = "#7651B9"
)

# These cluster-to-compartment mappings are inherited from the Step 19 program
# audit and the Step 20 T/NK inclusion decision. They are intentionally called
# broad expression compartments, not definitive fine cell types.
broad_compartment_map <- c(
  C01 = "T_or_NK", C02 = "T_or_NK", C03 = "Myeloid",
  C04 = "B_cell", C05 = "T_or_NK", C06 = "T_or_NK",
  C07 = "Myeloid", C08 = "Melanoma", C09 = "T_or_NK",
  C10 = "Other_low_information", C11 = "T_or_NK",
  C12 = "Other_low_information", C13 = "Other_low_information",
  C14 = "T_or_NK"
)
compartment_colors <- c(
  T_or_NK = "#3B82B8", B_cell = "#B5A642", Myeloid = "#59A14F",
  Melanoma = "#7651B9", Other_low_information = "#AFAFAF"
)

melanoma_anchor_genes <- c("MLANA", "PMEL", "SOX10", "MITF", "TYR", "DCT")

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

safe_wilcox_p <- function(difference) {
  difference <- difference[is.finite(difference)]
  if (length(difference) < 3L || all(difference == 0)) return(NA_real_)
  tryCatch(
    stats::wilcox.test(difference, mu = 0, exact = FALSE)$p.value,
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
    stats::median(sample(x, length(x), replace = TRUE), na.rm = TRUE)
  )
  ci <- stats::quantile(
    boot, c(0.025, 0.975), names = FALSE, na.rm = TRUE
  )
  setNames(as.numeric(ci), c("lower", "upper"))
}

summarize_differences <- function(x, metric, comparison, seed_offset = 0L) {
  x <- x[is.finite(x)]
  ci <- bootstrap_median_ci(x, seed_offset)
  data.table(
    comparison = comparison,
    metric = metric,
    n_independent_units = as.integer(length(x)),
    median_difference = if (length(x)) stats::median(x) else NA_real_,
    bootstrap_ci_lower = unname(ci[[1L]]),
    bootstrap_ci_upper = unname(ci[[2L]]),
    percent_positive = if (length(x)) 100 * mean(x > 0) else NA_real_,
    wilcoxon_signed_rank_p = safe_wilcox_p(x)
  )
}

# ------------------------------------------------------------------------------
# Load and validate the full Step 19 atlas
# ------------------------------------------------------------------------------

if (!file.exists(full_sce_file)) stop("Missing Step 19 SCE: ", full_sce_file)
if (!file.exists(cd8_sce_file)) stop("Missing Step 25 SCE: ", cd8_sce_file)

message("Loading the full Step 19 atlas")
full_sce <- readRDS(full_sce_file)
required_full_metadata <- c(
  "source_group", "subject_id", "specimen_id", "amplification_batch",
  "broad_cluster", "tcr_detected", "direct_t_cell_evidence",
  "direct_cd8_evidence"
)
missing_full_metadata <- setdiff(
  required_full_metadata,
  colnames(SummarizedExperiment::colData(full_sce))
)
if (length(missing_full_metadata)) {
  stop("Step 19 metadata is missing: ", paste(missing_full_metadata, collapse = ", "))
}
if (!all(c("counts", "logcounts") %in%
         SummarizedExperiment::assayNames(full_sce))) {
  stop("Step 19 must contain counts and logcounts assays.")
}
if (!"UMAP" %in% SingleCellExperiment::reducedDimNames(full_sce)) {
  stop("Step 19 UMAP coordinates are missing.")
}
if (!all(c("SPN", "PTPRC") %in% rownames(full_sce))) {
  stop("SPN or PTPRC is absent from the Step 19 gene set.")
}
if (anyDuplicated(colnames(full_sce))) stop("Duplicated Step 19 cell IDs occur.")

full_counts <- SummarizedExperiment::assay(full_sce, "counts")
full_logcounts <- SummarizedExperiment::assay(full_sce, "logcounts")
full_counts <- tryCatch(
  methods::as(full_counts, "dgCMatrix"), error = function(e) full_counts
)
full_logcounts <- tryCatch(
  methods::as(full_logcounts, "dgCMatrix"), error = function(e) full_logcounts
)

full_spn_count <- as.numeric(full_counts["SPN", ])
full_spn_log <- as.numeric(full_logcounts["SPN", ])
full_ptprc_count <- as.numeric(full_counts["PTPRC", ])
full_ptprc_log <- as.numeric(full_logcounts["PTPRC", ])
full_cluster <- as.character(full_sce$broad_cluster)
if (!all(full_cluster %in% names(broad_compartment_map))) {
  stop("Unexpected Step 19 broad clusters occur.")
}
full_compartment <- unname(broad_compartment_map[full_cluster])

present_melanoma_genes <- intersect(melanoma_anchor_genes, rownames(full_sce))
if (length(present_melanoma_genes) < 3L) {
  stop("Too few melanoma anchor genes are present.")
}
melanoma_detected_count <- as.integer(Matrix::colSums(
  full_counts[present_melanoma_genes, , drop = FALSE] > 0
))

full_sce$step26_broad_compartment <- factor(
  full_compartment, levels = names(compartment_colors)
)
full_sce$step26_spn_count <- full_spn_count
full_sce$step26_spn_logcounts <- full_spn_log
full_sce$step26_spn_detected <- full_spn_count > 0
full_sce$step26_ptprc_count <- full_ptprc_count
full_sce$step26_ptprc_logcounts <- full_ptprc_log
full_sce$step26_ptprc_detected <- full_ptprc_count > 0
full_sce$step26_melanoma_anchor_direct <- melanoma_detected_count >= 2L

full_dt <- as.data.table(SummarizedExperiment::colData(full_sce))
full_dt[, cell_id := colnames(full_sce)]
full_dt[, source_group := as.character(source_group)]
full_dt[, subject_id := as.character(subject_id)]
full_dt[, broad_cluster := as.character(broad_cluster)]
full_dt[, step26_broad_compartment := as.character(step26_broad_compartment)]

# ------------------------------------------------------------------------------
# Full-atlas SPN summaries and malignant-compartment audit
# ------------------------------------------------------------------------------

message("Quantifying SPN across broad cellular compartments")
spn_by_cluster <- full_dt[, .(
  n_cells = as.integer(.N),
  n_subjects = as.integer(uniqueN(subject_id)),
  n_spn_detected = as.integer(sum(step26_spn_detected)),
  percent_spn_detected = 100 * mean(step26_spn_detected),
  mean_spn_logcounts = mean(step26_spn_logcounts),
  median_spn_logcounts = as.numeric(stats::median(step26_spn_logcounts)),
  mean_spn_count = mean(step26_spn_count),
  n_ptprc_detected = as.integer(sum(step26_ptprc_detected)),
  percent_ptprc_detected = 100 * mean(step26_ptprc_detected)
), by = .(broad_cluster, step26_broad_compartment, source_group)]

spn_by_compartment_source <- full_dt[, .(
  n_cells = as.integer(.N),
  n_subjects = as.integer(uniqueN(subject_id)),
  n_spn_detected = as.integer(sum(step26_spn_detected)),
  percent_spn_detected = 100 * mean(step26_spn_detected),
  mean_spn_logcounts = mean(step26_spn_logcounts),
  median_spn_logcounts = as.numeric(stats::median(step26_spn_logcounts)),
  mean_spn_count = mean(step26_spn_count),
  percent_tcr = 100 * mean(as_flag(tcr_detected)),
  percent_direct_t = 100 * mean(as_flag(direct_t_cell_evidence)),
  percent_direct_cd8 = 100 * mean(as_flag(direct_cd8_evidence))
), by = .(step26_broad_compartment, source_group)]

patient_compartment_spn <- full_dt[, .(
  n_cells = as.integer(.N),
  spn_detected_percent = 100 * mean(step26_spn_detected),
  mean_spn_logcounts = mean(step26_spn_logcounts),
  median_spn_logcounts = as.numeric(stats::median(step26_spn_logcounts)),
  median_endogenous_umis = if ("endogenous_n_counts" %in% names(full_dt)) {
    as.numeric(stats::median(endogenous_n_counts, na.rm = TRUE))
  } else {
    NA_real_
  }
), by = .(subject_id, source_group, step26_broad_compartment)][
  n_cells >= minimum_cells_per_patient_group
]

patient_balanced_compartment_spn <- patient_compartment_spn[, .(
  n_patient_groups = as.integer(.N),
  median_patient_spn_detected_percent = as.numeric(
    stats::median(spn_detected_percent)
  ),
  q25_patient_spn_detected_percent = stats::quantile(spn_detected_percent, 0.25),
  q75_patient_spn_detected_percent = stats::quantile(spn_detected_percent, 0.75),
  median_patient_mean_spn_logcounts = as.numeric(
    stats::median(mean_spn_logcounts)
  ),
  q25_patient_mean_spn_logcounts = stats::quantile(mean_spn_logcounts, 0.25),
  q75_patient_mean_spn_logcounts = stats::quantile(mean_spn_logcounts, 0.75)
), by = .(source_group, step26_broad_compartment)]

melanoma_audit <- full_dt[step26_broad_compartment == "Melanoma", .(
  n_cells = as.integer(.N),
  n_melanoma_anchor_direct = as.integer(sum(step26_melanoma_anchor_direct)),
  n_spn_detected = as.integer(sum(step26_spn_detected)),
  percent_spn_detected = 100 * mean(step26_spn_detected),
  n_ptprc_detected = as.integer(sum(step26_ptprc_detected)),
  percent_ptprc_detected = 100 * mean(step26_ptprc_detected),
  n_spn_positive_ptprc_negative = as.integer(sum(
    step26_spn_detected & !step26_ptprc_detected
  )),
  percent_spn_positive_ptprc_negative = 100 * mean(
    step26_spn_detected & !step26_ptprc_detected
  ),
  n_anchor_supported_spn_positive_ptprc_negative = as.integer(sum(
    step26_melanoma_anchor_direct & step26_spn_detected &
      !step26_ptprc_detected
  )),
  mean_spn_logcounts = mean(step26_spn_logcounts)
), by = .(source_group, subject_id)]

melanoma_audit_summary <- full_dt[
  step26_broad_compartment == "Melanoma",
  .(
    n_cells = as.integer(.N),
    n_subjects = as.integer(uniqueN(subject_id)),
    n_spn_detected = as.integer(sum(step26_spn_detected)),
    percent_spn_detected = 100 * mean(step26_spn_detected),
    n_spn_positive_ptprc_negative = as.integer(sum(
      step26_spn_detected & !step26_ptprc_detected
    )),
    n_anchor_supported_spn_positive_ptprc_negative = as.integer(sum(
      step26_melanoma_anchor_direct & step26_spn_detected &
        !step26_ptprc_detected
    )),
    interpretive_warning = paste(
      "PTPRC non-detection does not prove malignant-cell origin because",
      "single-cell dropout and ambient RNA remain possible."
    )
  ),
  by = source_group
]

# ------------------------------------------------------------------------------
# Load Step 25 and add SPN measurements without redefining the CD8 states
# ------------------------------------------------------------------------------

message("Loading the Step 25 lineage-adjudicated tumor-CD8 object")
cd8_sce <- readRDS(cd8_sce_file)
required_cd8_metadata <- c(
  "subject_id", "specimen_id", "step25_retained_cd8", "step25_base_state",
  "step25_lineage_decision", "step25_signature_training_eligible",
  "step25_interferon_overlay", "step25_cytotoxic_overlay",
  "step25_proliferation_overlay", "tcr_productive_primary",
  "tcr_clonotype_key", "endogenous_n_counts", "endogenous_n_genes"
)
missing_cd8_metadata <- setdiff(
  required_cd8_metadata,
  colnames(SummarizedExperiment::colData(cd8_sce))
)
if (length(missing_cd8_metadata)) {
  stop("Step 25 metadata is missing: ", paste(missing_cd8_metadata, collapse = ", "))
}
if (!all(c("counts", "logcounts") %in%
         SummarizedExperiment::assayNames(cd8_sce))) {
  stop("Step 25 must contain counts and logcounts assays.")
}
if (!"UMAP_step21" %in% SingleCellExperiment::reducedDimNames(cd8_sce)) {
  stop("Step 21 UMAP is missing from Step 25.")
}
if (!"SPN" %in% rownames(cd8_sce)) stop("SPN is absent from Step 25.")

cd8_counts <- SummarizedExperiment::assay(cd8_sce, "counts")
cd8_logcounts <- SummarizedExperiment::assay(cd8_sce, "logcounts")
cd8_spn_count <- as.numeric(cd8_counts["SPN", ])
cd8_spn_log <- as.numeric(cd8_logcounts["SPN", ])
cd8_sce$step26_spn_count <- cd8_spn_count
cd8_sce$step26_spn_logcounts <- cd8_spn_log
cd8_sce$step26_spn_detected <- cd8_spn_count > 0

cd8_dt <- as.data.table(SummarizedExperiment::colData(cd8_sce))
cd8_dt[, cell_id := colnames(cd8_sce)]
cd8_dt[, subject_id := as.character(subject_id)]
cd8_dt[, step25_base_state := as.character(step25_base_state)]
cd8_dt[, log1p_endogenous_umis := log1p(as.numeric(endogenous_n_counts))]
retained_dt <- cd8_dt[step25_retained_cd8 == TRUE]
if (!setequal(unique(retained_dt$step25_base_state), base_state_levels)) {
  stop("Unexpected Step 25 base states occur among retained cells.")
}

# ------------------------------------------------------------------------------
# Patient-balanced SPN comparison between the two CD8 base states
# ------------------------------------------------------------------------------

message("Testing SPN across CD8 base states with patients as replicates")
cd8_spn_by_base_state <- retained_dt[, .(
  n_cells = as.integer(.N),
  n_subjects = as.integer(uniqueN(subject_id)),
  n_spn_detected = as.integer(sum(step26_spn_detected)),
  percent_spn_detected = 100 * mean(step26_spn_detected),
  mean_spn_logcounts = mean(step26_spn_logcounts),
  median_spn_logcounts = as.numeric(stats::median(step26_spn_logcounts)),
  mean_spn_count = mean(step26_spn_count),
  median_endogenous_umis = as.numeric(stats::median(endogenous_n_counts)),
  percent_productive_tcr = 100 * mean(as_flag(tcr_productive_primary))
), by = step25_base_state]

cd8_patient_state <- retained_dt[, .(
  n_cells = as.integer(.N),
  spn_detected_percent = 100 * mean(step26_spn_detected),
  mean_spn_logcounts = mean(step26_spn_logcounts),
  median_spn_logcounts = as.numeric(stats::median(step26_spn_logcounts)),
  mean_spn_count = mean(step26_spn_count),
  median_endogenous_umis = as.numeric(stats::median(endogenous_n_counts)),
  median_endogenous_genes = as.numeric(stats::median(endogenous_n_genes))
), by = .(subject_id, step25_base_state)][
  n_cells >= minimum_cells_per_patient_group
]

state_detection_wide <- dcast(
  cd8_patient_state,
  subject_id ~ step25_base_state,
  value.var = "spn_detected_percent"
)
state_expression_wide <- dcast(
  cd8_patient_state,
  subject_id ~ step25_base_state,
  value.var = "mean_spn_logcounts"
)
state_umi_wide <- dcast(
  cd8_patient_state,
  subject_id ~ step25_base_state,
  value.var = "median_endogenous_umis"
)

paired_state_differences <- merge(
  state_detection_wide, state_expression_wide,
  by = "subject_id", suffixes = c("_detection", "_expression")
)
paired_state_differences <- merge(
  paired_state_differences, state_umi_wide,
  by = "subject_id", suffixes = c("", "_umi")
)
required_paired_columns <- c(
  "CD8_memory_transitional_detection", "CD8_dysfunctional_detection",
  "CD8_memory_transitional_expression", "CD8_dysfunctional_expression",
  "CD8_memory_transitional", "CD8_dysfunctional"
)
if (!all(required_paired_columns %in% names(paired_state_differences))) {
  stop("Could not align the patient-level CD8 state summaries.")
}
paired_state_differences[, detection_difference_dysfunctional_minus_memory :=
  CD8_dysfunctional_detection - CD8_memory_transitional_detection]
paired_state_differences[, expression_difference_dysfunctional_minus_memory :=
  CD8_dysfunctional_expression - CD8_memory_transitional_expression]
paired_state_differences[, umi_difference_dysfunctional_minus_memory :=
  CD8_dysfunctional - CD8_memory_transitional]

paired_state_complete <- paired_state_differences[
  is.finite(detection_difference_dysfunctional_minus_memory) &
    is.finite(expression_difference_dysfunctional_minus_memory)
]

state_paired_summary <- rbindlist(list(
  summarize_differences(
    paired_state_complete$detection_difference_dysfunctional_minus_memory,
    "SPN detection percentage points",
    "Dysfunctional minus memory/transitional", 100L
  ),
  summarize_differences(
    paired_state_complete$expression_difference_dysfunctional_minus_memory,
    "Mean SPN logcounts",
    "Dysfunctional minus memory/transitional", 200L
  ),
  summarize_differences(
    paired_state_complete$umi_difference_dysfunctional_minus_memory,
    "Median endogenous UMIs",
    "Dysfunctional minus memory/transitional", 300L
  )
))
# The UMI comparison is a depth diagnostic rather than an SPN hypothesis, so it
# is deliberately excluded from the SPN multiple-testing family.
state_paired_summary[, fdr := NA_real_]
state_paired_summary[grepl("^SPN", metric), fdr := p.adjust(
  wilcoxon_signed_rank_p, method = "BH"
)]

# Each patient contributes one depth-adjusted state coefficient. The across-
# patient summary is descriptive and patient-balanced; it is not a cell-level
# pseudo-replication test.
fit_patient_models <- function(patient_table) {
  patient_table <- copy(patient_table)
  patient_table[, dysfunctional := as.integer(
    step25_base_state == "CD8_dysfunctional"
  )]
  if (uniqueN(patient_table$step25_base_state) < 2L ||
      min(table(patient_table$step25_base_state)) <
        minimum_cells_per_patient_group) {
    return(data.table())
  }

  linear_fit <- tryCatch(
    stats::lm(
      step26_spn_logcounts ~ dysfunctional + log1p_endogenous_umis,
      data = patient_table
    ),
    error = function(e) NULL
  )
  logistic_fit <- NULL
  if (uniqueN(patient_table$step26_spn_detected) == 2L) {
    logistic_fit <- tryCatch(
      suppressWarnings(stats::glm(
        step26_spn_detected ~ dysfunctional + log1p_endogenous_umis,
        family = stats::binomial(), data = patient_table
      )),
      error = function(e) NULL
    )
  }

  linear_coefficient <- if (!is.null(linear_fit) &&
      "dysfunctional" %in% names(stats::coef(linear_fit))) {
    unname(stats::coef(linear_fit)[["dysfunctional"]])
  } else {
    NA_real_
  }
  logistic_coefficient <- if (!is.null(logistic_fit) &&
      "dysfunctional" %in% names(stats::coef(logistic_fit))) {
    unname(stats::coef(logistic_fit)[["dysfunctional"]])
  } else {
    NA_real_
  }
  data.table(
    n_cells = as.integer(nrow(patient_table)),
    n_memory = as.integer(sum(patient_table$dysfunctional == 0L)),
    n_dysfunctional = as.integer(sum(patient_table$dysfunctional == 1L)),
    linear_state_coefficient = linear_coefficient,
    logistic_state_log_odds = logistic_coefficient
  )
}

patient_model_results <- retained_dt[, fit_patient_models(.SD), by = subject_id]
patient_model_summary <- rbindlist(list(
  summarize_differences(
    patient_model_results$linear_state_coefficient,
    "Within-patient depth-adjusted SPN logcount coefficient",
    "Dysfunctional versus memory/transitional", 400L
  ),
  summarize_differences(
    patient_model_results$logistic_state_log_odds,
    "Within-patient depth-adjusted SPN detection log-odds",
    "Dysfunctional versus memory/transitional", 500L
  )
))
patient_model_summary[, fdr := p.adjust(wilcoxon_signed_rank_p, method = "BH")]

# ------------------------------------------------------------------------------
# SPN association with non-exclusive CD8 overlays
# ------------------------------------------------------------------------------

message("Testing SPN across CD8 overlays within each base state")
overlay_columns <- c(
  Interferon = "step25_interferon_overlay",
  Cytotoxic = "step25_cytotoxic_overlay",
  Proliferation = "step25_proliferation_overlay"
)

overlay_patient_tables <- list()
overlay_summaries <- list()
overlay_counter <- 0L
for (overlay_name in names(overlay_columns)) {
  overlay_column <- overlay_columns[[overlay_name]]
  for (base_state in base_state_levels) {
    overlay_counter <- overlay_counter + 1L
    temp <- retained_dt[step25_base_state == base_state]
    temp[, overlay_present := as_flag(get(overlay_column))]
    patient_overlay <- temp[, .(
      n_cells = as.integer(.N),
      spn_detected_percent = 100 * mean(step26_spn_detected),
      mean_spn_logcounts = mean(step26_spn_logcounts),
      median_endogenous_umis = as.numeric(stats::median(endogenous_n_counts))
    ), by = .(subject_id, overlay_present)]

    patient_counts <- dcast(
      patient_overlay, subject_id ~ overlay_present,
      value.var = "n_cells", fill = 0L
    )
    if (!all(c("FALSE", "TRUE") %in% names(patient_counts))) next
    eligible_subjects <- patient_counts[
      get("FALSE") >= minimum_cells_per_patient_group &
        get("TRUE") >= minimum_overlay_positive_cells,
      subject_id
    ]
    patient_overlay <- patient_overlay[subject_id %in% eligible_subjects]
    if (!nrow(patient_overlay)) next
    patient_overlay[, `:=`(
      overlay = overlay_name,
      step25_base_state = base_state
    )]
    overlay_patient_tables[[length(overlay_patient_tables) + 1L]] <-
      patient_overlay

    detection_wide <- dcast(
      patient_overlay, subject_id ~ overlay_present,
      value.var = "spn_detected_percent"
    )
    expression_wide <- dcast(
      patient_overlay, subject_id ~ overlay_present,
      value.var = "mean_spn_logcounts"
    )
    if (!all(c("FALSE", "TRUE") %in% names(detection_wide)) ||
        !all(c("FALSE", "TRUE") %in% names(expression_wide))) next

    detection_difference <-
      detection_wide[["TRUE"]] - detection_wide[["FALSE"]]
    expression_difference <-
      expression_wide[["TRUE"]] - expression_wide[["FALSE"]]
    comparison_label <- paste(
      overlay_name, "positive minus negative within", base_state_labels[[base_state]]
    )
    summary_table <- rbindlist(list(
      summarize_differences(
        detection_difference, "SPN detection percentage points",
        comparison_label, 600L + overlay_counter
      ),
      summarize_differences(
        expression_difference, "Mean SPN logcounts",
        comparison_label, 700L + overlay_counter
      )
    ))
    summary_table[, `:=`(
      overlay = overlay_name,
      step25_base_state = base_state
    )]
    overlay_summaries[[length(overlay_summaries) + 1L]] <- summary_table
  }
}

overlay_patient_spn <- if (length(overlay_patient_tables)) {
  rbindlist(overlay_patient_tables, fill = TRUE)
} else {
  data.table()
}
overlay_paired_summary <- if (length(overlay_summaries)) {
  rbindlist(overlay_summaries, fill = TRUE)
} else {
  data.table()
}
if (nrow(overlay_paired_summary)) {
  overlay_paired_summary[, fdr := p.adjust(
    wilcoxon_signed_rank_p, method = "BH"
  )]
}

# ------------------------------------------------------------------------------
# Productive-clonotype matched state comparison
# ------------------------------------------------------------------------------

message("Testing SPN within productive clonotypes spanning both base states")
clone_dt <- retained_dt[
  as_flag(tcr_productive_primary) &
    !is.na(tcr_clonotype_key) & nzchar(as.character(tcr_clonotype_key))
]
clone_dt[, clone_id := paste(subject_id, tcr_clonotype_key, sep = "::")]
clone_state_summary <- clone_dt[, .(
  n_cells = as.integer(.N),
  spn_detected_percent = 100 * mean(step26_spn_detected),
  mean_spn_logcounts = mean(step26_spn_logcounts),
  median_endogenous_umis = as.numeric(stats::median(endogenous_n_counts))
), by = .(clone_id, subject_id, step25_base_state)]

build_clone_pairs <- function(minimum_per_state) {
  clone_counts <- dcast(
    clone_state_summary,
    clone_id + subject_id ~ step25_base_state,
    value.var = "n_cells", fill = 0L
  )
  if (!all(base_state_levels %in% names(clone_counts))) return(data.table())
  eligible <- clone_counts[
    CD8_memory_transitional >= minimum_per_state &
      CD8_dysfunctional >= minimum_per_state,
    .(clone_id, subject_id)
  ]
  if (!nrow(eligible)) return(data.table())

  eligible_summary <- merge(
    clone_state_summary, eligible,
    by = c("clone_id", "subject_id")
  )
  detection <- dcast(
    eligible_summary,
    clone_id + subject_id ~ step25_base_state,
    value.var = "spn_detected_percent"
  )
  expression <- dcast(
    eligible_summary,
    clone_id + subject_id ~ step25_base_state,
    value.var = "mean_spn_logcounts"
  )
  sizes <- dcast(
    eligible_summary,
    clone_id + subject_id ~ step25_base_state,
    value.var = "n_cells"
  )
  output <- merge(
    detection, expression, by = c("clone_id", "subject_id"),
    suffixes = c("_detection", "_expression")
  )
  output <- merge(
    output, sizes, by = c("clone_id", "subject_id"),
    suffixes = c("", "_n")
  )
  output[, detection_difference_dysfunctional_minus_memory :=
    CD8_dysfunctional_detection - CD8_memory_transitional_detection]
  output[, expression_difference_dysfunctional_minus_memory :=
    CD8_dysfunctional_expression - CD8_memory_transitional_expression]
  output[, minimum_cells_per_state := as.integer(minimum_per_state)]
  output
}

clone_pairs_primary <- build_clone_pairs(minimum_clone_cells_per_state_primary)
clone_pairs_sensitivity <- build_clone_pairs(
  minimum_clone_cells_per_state_sensitivity
)
clone_pair_results <- rbindlist(
  list(clone_pairs_primary, clone_pairs_sensitivity), fill = TRUE
)

clone_paired_summary <- if (nrow(clone_pair_results)) {
  rbindlist(lapply(
    sort(unique(clone_pair_results$minimum_cells_per_state)),
    function(minimum_required) {
      current <- clone_pair_results[
        minimum_cells_per_state == minimum_required
      ]
      rbindlist(list(
        summarize_differences(
          current$detection_difference_dysfunctional_minus_memory,
          "SPN detection percentage points",
          paste0(
            "Dysfunctional minus memory within productive clones; >=",
            minimum_required, " cell(s) per state"
          ), 800L + minimum_required
        ),
        summarize_differences(
          current$expression_difference_dysfunctional_minus_memory,
          "Mean SPN logcounts",
          paste0(
            "Dysfunctional minus memory within productive clones; >=",
            minimum_required, " cell(s) per state"
          ), 900L + minimum_required
        )
      ))[, minimum_cells_per_state := minimum_required]
    }
  ), fill = TRUE)
} else {
  data.table()
}
if (nrow(clone_paired_summary)) {
  clone_paired_summary[, fdr := p.adjust(
    wilcoxon_signed_rank_p, method = "BH"
  )]
}

# ------------------------------------------------------------------------------
# Output tables and validation
# ------------------------------------------------------------------------------

message("Writing Step 26 tables")
fwrite(spn_by_cluster, file.path(results_dir, "26_spn_by_broad_cluster.csv"))
fwrite(
  spn_by_compartment_source,
  file.path(results_dir, "26_spn_by_compartment_and_source.csv")
)
fwrite(
  patient_compartment_spn,
  file.path(results_dir, "26_patient_compartment_spn.csv")
)
fwrite(
  patient_balanced_compartment_spn,
  file.path(results_dir, "26_patient_balanced_compartment_spn.csv")
)
fwrite(melanoma_audit, file.path(results_dir, "26_melanoma_spn_ptprc_by_patient.csv"))
fwrite(
  melanoma_audit_summary,
  file.path(results_dir, "26_melanoma_spn_ptprc_summary.csv")
)
fwrite(cd8_spn_by_base_state, file.path(results_dir, "26_cd8_spn_by_base_state.csv"))
fwrite(cd8_patient_state, file.path(results_dir, "26_cd8_patient_state_spn.csv"))
fwrite(
  paired_state_complete,
  file.path(results_dir, "26_cd8_paired_patient_state_differences.csv")
)
fwrite(
  state_paired_summary,
  file.path(results_dir, "26_cd8_paired_state_summary.csv")
)
fwrite(
  patient_model_results,
  file.path(results_dir, "26_cd8_patient_depth_adjusted_models.csv")
)
fwrite(
  patient_model_summary,
  file.path(results_dir, "26_cd8_patient_depth_adjusted_summary.csv")
)
fwrite(
  overlay_patient_spn,
  file.path(results_dir, "26_cd8_overlay_patient_spn.csv")
)
fwrite(
  overlay_paired_summary,
  file.path(results_dir, "26_cd8_overlay_paired_summary.csv")
)
fwrite(
  clone_state_summary,
  file.path(results_dir, "26_cd8_productive_clone_state_spn.csv")
)
fwrite(
  clone_pair_results,
  file.path(results_dir, "26_cd8_clone_paired_differences.csv")
)
fwrite(
  clone_paired_summary,
  file.path(results_dir, "26_cd8_clone_paired_summary.csv")
)

full_accounting_validation <- sum(spn_by_cluster$n_cells) == ncol(full_sce)
cd8_accounting_validation <- sum(cd8_spn_by_base_state$n_cells) ==
  sum(as_flag(cd8_sce$step25_retained_cd8))
spn_measurement_alignment_validation <-
  identical(colnames(cd8_sce), cd8_dt$cell_id) &&
  identical(colnames(full_sce), full_dt$cell_id)
if (!full_accounting_validation) stop("Full-atlas accounting failed.")
if (!cd8_accounting_validation) stop("Retained-CD8 accounting failed.")
if (!spn_measurement_alignment_validation) stop("SPN cell alignment failed.")

dataset_qc <- data.table(
  metric = c(
    "full_atlas_cells", "full_atlas_subjects", "tumor_cells", "pbmc_cells",
    "full_atlas_spn_detected_cells", "full_atlas_spn_detected_percent",
    "melanoma_compartment_cells", "melanoma_spn_detected_cells",
    "melanoma_anchor_supported_spn_positive_ptprc_negative_cells",
    "retained_cd8_cells", "retained_cd8_spn_detected_cells",
    "retained_cd8_spn_detected_percent", "paired_state_patients",
    "depth_adjusted_model_patients", "productive_cross_state_clones",
    "productive_cross_state_clones_minimum_two_cells_per_state",
    "full_accounting_validation", "cd8_accounting_validation",
    "spn_measurement_alignment_validation", "signature_exported",
    "output_cd8_sce_size_mb"
  ),
  value = as.character(c(
    ncol(full_sce), uniqueN(full_dt$subject_id),
    sum(full_dt$source_group == "Tumor"), sum(full_dt$source_group == "PBMC"),
    sum(full_dt$step26_spn_detected), 100 * mean(full_dt$step26_spn_detected),
    sum(full_dt$step26_broad_compartment == "Melanoma"),
    sum(full_dt$step26_broad_compartment == "Melanoma" &
          full_dt$step26_spn_detected),
    sum(full_dt$step26_broad_compartment == "Melanoma" &
          full_dt$step26_melanoma_anchor_direct &
          full_dt$step26_spn_detected & !full_dt$step26_ptprc_detected),
    nrow(retained_dt), sum(retained_dt$step26_spn_detected),
    100 * mean(retained_dt$step26_spn_detected),
    nrow(paired_state_complete), nrow(patient_model_results),
    nrow(clone_pairs_primary), nrow(clone_pairs_sensitivity),
    as.integer(full_accounting_validation),
    as.integer(cd8_accounting_validation),
    as.integer(spn_measurement_alignment_validation), 0L, NA_real_
  ))
)
fwrite(dataset_qc, file.path(results_dir, "26_dataset_qc.csv"))

# ------------------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------------------

message("Rendering Step 26 figures")
full_umap <- as.data.table(SingleCellExperiment::reducedDim(full_sce, "UMAP"))
setnames(full_umap, c("UMAP1", "UMAP2"))
full_plot_dt <- cbind(
  full_umap,
  full_dt[, .(
    source_group, step26_broad_compartment, step26_spn_detected,
    step26_spn_logcounts
  )]
)

p_full_umap <- ggplot(full_plot_dt, aes(UMAP1, UMAP2)) +
  geom_point(color = "#E1E1E1", size = 0.16) +
  geom_point(
    data = full_plot_dt[step26_spn_detected == TRUE],
    aes(color = step26_spn_logcounts), size = 0.34, alpha = 0.85
  ) +
  scale_color_viridis_c(option = "magma", direction = -1) +
  labs(
    title = "SPN/CD43 expression across the full GSE123139 atlas",
    subtitle = "Colored cells have at least one SPN raw count",
    color = "SPN logcounts"
  ) +
  theme_project()
save_png("26A_full_atlas_umap_spn.png", p_full_umap, 11, 8)

compartment_plot_dt <- patient_balanced_compartment_spn[
  n_patient_groups >= 3L
]
p_compartment <- ggplot(
  compartment_plot_dt,
  aes(
    reorder(step26_broad_compartment, median_patient_spn_detected_percent),
    median_patient_spn_detected_percent,
    color = step26_broad_compartment
  )
) +
  geom_errorbar(aes(
    ymin = q25_patient_spn_detected_percent,
    ymax = q75_patient_spn_detected_percent
  ), width = 0.2) +
  geom_point(size = 3) +
  facet_wrap(~source_group, scales = "free_x") +
  scale_color_manual(values = compartment_colors, guide = "none") +
  labs(
    title = "Patient-balanced SPN detection by broad compartment",
    subtitle = "Points are medians; bars show the interquartile range across patient-compartment groups",
    x = NULL, y = "Median patient SPN-detected cells (%)"
  ) +
  theme_project() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_png("26B_patient_balanced_spn_by_compartment.png", p_compartment, 12, 7)

cd8_umap <- as.data.table(SingleCellExperiment::reducedDim(cd8_sce, "UMAP_step21"))
setnames(cd8_umap, c("UMAP1", "UMAP2"))
cd8_plot_dt <- cbind(
  cd8_umap,
  cd8_dt[, .(
    step25_retained_cd8, step25_base_state,
    step26_spn_detected, step26_spn_logcounts
  )]
)
cd8_retained_plot <- cd8_plot_dt[step25_retained_cd8 == TRUE]
p_cd8_umap <- ggplot(cd8_retained_plot, aes(UMAP1, UMAP2)) +
  geom_point(color = "#E1E1E1", size = 0.22) +
  geom_point(
    data = cd8_retained_plot[step26_spn_detected == TRUE],
    aes(color = step26_spn_logcounts), size = 0.48, alpha = 0.90
  ) +
  scale_color_viridis_c(option = "magma", direction = -1) +
  labs(
    title = "SPN/CD43 across the purified tumor-CD8 compartment",
    subtitle = "CD8 states were defined without SPN",
    color = "SPN logcounts"
  ) +
  theme_project()
save_png("26C_purified_cd8_umap_spn.png", p_cd8_umap, 11, 8)

paired_plot <- melt(
  cd8_patient_state,
  id.vars = c("subject_id", "step25_base_state", "n_cells"),
  measure.vars = c("spn_detected_percent", "mean_spn_logcounts"),
  variable.name = "metric", value.name = "value"
)
paired_plot[, state_label := factor(
  base_state_labels[step25_base_state],
  levels = unname(base_state_labels[base_state_levels])
)]
paired_plot[, metric := factor(
  metric,
  levels = c("spn_detected_percent", "mean_spn_logcounts"),
  labels = c("SPN detection (%)", "Mean SPN logcounts")
)]
p_paired <- ggplot(
  paired_plot,
  aes(state_label, value, group = subject_id)
) +
  geom_line(color = "#B6B6B6", alpha = 0.75) +
  geom_point(aes(color = step25_base_state), size = 2) +
  facet_wrap(~metric, scales = "free_y") +
  scale_color_manual(values = base_state_colors, guide = "none") +
  labs(
    title = "Patient-paired SPN across the two CD8 base states",
    subtitle = paste0(
      "Patient-state groups require at least ",
      minimum_cells_per_patient_group, " cells"
    ),
    x = NULL, y = NULL
  ) +
  theme_project()
save_png("26D_patient_paired_cd8_state_spn.png", p_paired, 11, 6)

model_plot <- melt(
  patient_model_results,
  id.vars = c("subject_id", "n_cells", "n_memory", "n_dysfunctional"),
  measure.vars = c("linear_state_coefficient", "logistic_state_log_odds"),
  variable.name = "model", value.name = "coefficient",
  na.rm = TRUE
)
model_plot[, model := factor(
  model,
  levels = c("linear_state_coefficient", "logistic_state_log_odds"),
  labels = c("SPN logcount coefficient", "SPN detection log-odds")
)]
p_models <- ggplot(model_plot, aes(model, coefficient)) +
  geom_hline(yintercept = 0, color = "#888888") +
  geom_boxplot(outlier.shape = NA, width = 0.5, fill = "white") +
  geom_jitter(width = 0.12, aes(color = model), size = 2, alpha = 0.85) +
  guides(color = "none") +
  labs(
    title = "Within-patient, depth-adjusted SPN association with dysfunction",
    subtitle = "Positive coefficients indicate higher SPN in dysfunctional cells",
    x = NULL, y = "Patient-specific dysfunctional-state coefficient"
  ) +
  theme_project()
save_png("26E_patient_depth_adjusted_state_coefficients.png", p_models, 9, 7)

if (nrow(overlay_paired_summary)) {
  overlay_plot_summary <- overlay_paired_summary[
    metric == "SPN detection percentage points"
  ]
  overlay_plot_summary[, base_state_label :=
    base_state_labels[step25_base_state]]
  p_overlay <- ggplot(
    overlay_plot_summary,
    aes(overlay, median_difference, color = base_state_label)
  ) +
    geom_hline(yintercept = 0, color = "#888888") +
    geom_errorbar(aes(
      ymin = bootstrap_ci_lower, ymax = bootstrap_ci_upper
    ), width = 0.15, position = position_dodge(width = 0.35)) +
    geom_point(size = 3, position = position_dodge(width = 0.35)) +
    scale_color_manual(
      values = setNames(base_state_colors, unname(base_state_labels))
    ) +
    labs(
      title = "Patient-balanced SPN association with CD8 overlays",
      subtitle = "Positive values indicate greater SPN detection in overlay-positive cells",
      x = NULL, y = "Median paired difference (percentage points)",
      color = "Base state"
    ) +
    theme_project()
  save_png("26F_spn_by_cd8_overlay.png", p_overlay, 10, 7)
}

if (nrow(clone_pairs_primary)) {
  clone_plot <- clone_pairs_primary[
    is.finite(expression_difference_dysfunctional_minus_memory)
  ]
  p_clone <- ggplot(
    clone_plot,
    aes(expression_difference_dysfunctional_minus_memory)
  ) +
    geom_vline(xintercept = 0, color = "#888888") +
    geom_histogram(bins = 30, fill = "#7651B9", color = "white") +
    labs(
      title = "SPN state differences within shared productive clonotypes",
      subtitle = "Positive values indicate higher SPN in dysfunctional cells of the same clone",
      x = "Dysfunctional minus memory/transitional mean SPN logcounts",
      y = "Productive clonotypes"
    ) +
    theme_project()
  save_png("26G_clone_matched_spn_state_difference.png", p_clone, 9, 6)
}

# ------------------------------------------------------------------------------
# Save annotated CD8 object and completion record
# ------------------------------------------------------------------------------

message("Saving the Step 26 SPN-annotated tumor-CD8 object")
saveRDS(cd8_sce, output_cd8_sce_file, compress = FALSE)
dataset_qc[metric == "output_cd8_sce_size_mb", value := as.character(
  file.info(output_cd8_sce_file)$size / 1024 ^ 2
)]
fwrite(dataset_qc, file.path(results_dir, "26_dataset_qc.csv"))

writeLines(
  c(
    "Step 26 completed successfully.",
    paste("Full atlas cells:", ncol(full_sce)),
    paste("SPN-detected atlas cells:", sum(full_dt$step26_spn_detected)),
    paste("Retained tumor-CD8 cells:", nrow(retained_dt)),
    paste("SPN-detected retained CD8 cells:", sum(retained_dt$step26_spn_detected)),
    paste("Paired patients for state comparison:", nrow(paired_state_complete)),
    paste("Depth-adjusted patient models:", nrow(patient_model_results)),
    paste("Productive clones spanning both states:", nrow(clone_pairs_primary)),
    "SPN was not used to define compartments, states, overlays, or eligibility.",
    "No deconvolution signature was exported.",
    paste(
      "Malignant SPN+/PTPRC- candidates remain provisional because",
      "dropout and ambient RNA can mimic this pattern."
    )
  ),
  file.path(results_dir, "26_completion.txt")
)

message("Step 26 complete: SPN source and patient-aware CD8 associations quantified")
