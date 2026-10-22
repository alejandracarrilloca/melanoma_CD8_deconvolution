#!/usr/bin/env Rscript

# ==============================================================================
# Step 28: Within-state biology associated with SPN/CD43 in tumor CD8 cells
# ==============================================================================
#
# SPN was not used to define the Step 25 base states or Step 27 programs.
# This step asks whether SPN identifies additional biology *within* each frozen
# CD8 base state. Patients, not cells, are the independent replicates.
#
# Primary analysis:
#   - paired patient pseudobulk SPN-detected versus SPN-undetected contrasts,
#     performed separately in each base state;
#   - library-size normalization plus patient blocking;
#   - adjustment for patient-group sequencing depth and overlay composition;
#   - SPN itself is excluded from discovery to avoid a tautological result.
#
# Sensitivities:
#   - cells without interferon/cytotoxic/proliferation overlays
#   - strict direct-CD8 cells;
#   - patient-wise program models adjusted for depth and overlays;
#   - productive-clonotype-controlled program models.
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

project_root <- normalizePath(getwd(), mustWork = TRUE)
analysis_dir <- file.path(
  project_root, "data", "analysis_ready", "scRNA", "GSE123139"
)
results_dir <- file.path(project_root, "results", "28")
figures_dir <- file.path(results_dir, "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

input_sce_file <- file.path(
  analysis_dir, "27_GSE123139_tumor_cd8_biological_programs_sce.rds"
)
output_sce_file <- file.path(
  analysis_dir, "28_GSE123139_tumor_cd8_within_state_spn_biology_sce.rds"
)

random_seed <- 123139L
bootstrap_iterations <- 5000L
minimum_cells_per_patient_group <- 10L
minimum_cells_per_patient_program_model <- 80L
minimum_informative_clones <- 3L
use_robust_edger <- requireNamespace("statmod", quietly = TRUE)

base_states <- c("CD8_memory_transitional", "CD8_dysfunctional")
state_labels <- c(
  CD8_memory_transitional = "Memory/transitional",
  CD8_dysfunctional = "Exhaustion-like dysfunctional"
)
state_colors <- c(
  CD8_memory_transitional = "#3B82B8",
  CD8_dysfunctional = "#7651B9"
)

theme_project <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.18)),
      plot.subtitle = element_text(color = "#5A5A5A"),
      axis.title = element_text(face = "bold"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
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
    stats::median(sample(x, length(x), replace = TRUE)
  ))
  setNames(
    as.numeric(stats::quantile(boot, c(0.025, 0.975), names = FALSE)),
    c("lower", "upper")
  )
}

message("Loading completed Step 27 object")
if (!file.exists(input_sce_file)) stop("Missing Step 27 SCE: ", input_sce_file)
sce <- readRDS(input_sce_file)

required_metadata <- c(
  "subject_id", "step25_retained_cd8", "step25_base_state",
  "step25_interferon_overlay", "step25_cytotoxic_overlay",
  "step25_proliferation_overlay", "endogenous_n_counts",
  "endogenous_n_genes", "step26_spn_detected", "step26_spn_logcounts",
  "tcr_productive_primary", "tcr_clonotype_key"
)
missing_metadata <- setdiff(
  required_metadata, colnames(SummarizedExperiment::colData(sce))
)
if (length(missing_metadata)) {
  stop("Step 27 metadata is missing: ", paste(missing_metadata, collapse = ", "))
}
if (!all(c("counts", "logcounts") %in%
         SummarizedExperiment::assayNames(sce))) {
  stop("Step 27 object must contain counts and logcounts assays.")
}
if (anyDuplicated(colnames(sce))) stop("Duplicated cell IDs occur.")

counts <- SummarizedExperiment::assay(sce, "counts")
logcounts <- SummarizedExperiment::assay(sce, "logcounts")
counts <- tryCatch(methods::as(counts, "dgCMatrix"), error = function(e) counts)

cell <- as.data.table(SummarizedExperiment::colData(sce))
cell[, cell_index := .I]
cell[, cell_id := colnames(sce)]
cell[, subject_id := as.character(subject_id)]
cell[, step25_base_state := as.character(step25_base_state)]
cell[, log1p_endogenous_umis := log1p(as.numeric(endogenous_n_counts))]
cell[, spn_detected := as.integer(as_flag(step26_spn_detected))]
cell[, spn_group := factor(
  fifelse(spn_detected == 1L, "SPN_detected", "SPN_not_detected"),
  levels = c("SPN_not_detected", "SPN_detected")
)]
cell[, interferon_overlay := as.integer(as_flag(step25_interferon_overlay))]
cell[, cytotoxic_overlay := as.integer(as_flag(step25_cytotoxic_overlay))]
cell[, proliferation_overlay := as.integer(as_flag(step25_proliferation_overlay))]
cell[, any_overlay := interferon_overlay == 1L |
  cytotoxic_overlay == 1L | proliferation_overlay == 1L]

strict_candidates <- c(
  "step27_strict_cd8", "strict_direct_cd8", "step21_strict_direct_cd8",
  "direct_cd8_evidence"
)
strict_column <- strict_candidates[strict_candidates %in% names(cell)][1L]
if (is.na(strict_column)) {
  warning("No strict direct-CD8 column found; strict sensitivity omitted.")
  cell[, strict_cd8 := FALSE]
} else {
  cell[, strict_cd8 := as_flag(get(strict_column))]
}

retained_index <- which(
  as_flag(cell$step25_retained_cd8) &
    cell$step25_base_state %in% base_states
)
if (!length(retained_index)) stop("No retained Step 25 CD8 cells were found.")

# ------------------------------------------------------------------------------
# State-specific paired pseudobulk discovery
# ------------------------------------------------------------------------------

message("Running within-state patient-paired SPN pseudobulk contrasts")

run_spn_pseudobulk <- function(cell_indices, state_name, analysis_name) {
  indices <- cell_indices[cell$step25_base_state[cell_indices] == state_name]
  meta <- cell[indices, .(
    cell_index, subject_id, spn_group, log1p_endogenous_umis,
    interferon_overlay, cytotoxic_overlay, proliferation_overlay
  )]

  group_counts <- meta[, .(n_cells = .N), by = .(subject_id, spn_group)]
  eligible_subjects <- group_counts[
    n_cells >= minimum_cells_per_patient_group,
    .(n_groups = uniqueN(spn_group)), by = subject_id
  ][n_groups == 2L, subject_id]
  meta <- meta[subject_id %in% eligible_subjects]

  empty_qc <- function() data.table(
    analysis = analysis_name, step25_base_state = state_name,
    n_cells = as.integer(nrow(meta)),
    n_patients = as.integer(uniqueN(meta$subject_id)),
    n_pseudobulk_samples = 0L, n_genes_tested = 0L, n_fdr_005 = 0L
  )
  if (uniqueN(meta$subject_id) < 3L) {
    return(list(de = data.table(), qc = empty_qc()))
  }

  meta[, group_key := paste(subject_id, spn_group, sep = "||")]
  group_levels <- unique(meta$group_key)
  group_factor <- factor(meta$group_key, levels = group_levels)
  aggregation <- Matrix::sparse.model.matrix(~ 0 + group_factor)
  colnames(aggregation) <- group_levels
  pb_counts <- counts[, meta$cell_index, drop = FALSE] %*% aggregation

  pb_meta <- meta[, .(
    subject_id = subject_id[1L],
    spn_group = spn_group[1L],
    mean_log1p_umis = mean(log1p_endogenous_umis),
    interferon_fraction = mean(interferon_overlay),
    cytotoxic_fraction = mean(cytotoxic_overlay),
    proliferation_fraction = mean(proliferation_overlay)
  ), by = group_key]
  pb_meta <- pb_meta[match(group_levels, pb_meta$group_key)]
  pb_meta[, subject_id := factor(subject_id)]
  pb_meta[, spn_group := factor(
    as.character(spn_group), levels = c("SPN_not_detected", "SPN_detected")
  )]
  rownames(pb_meta) <- group_levels
  colnames(pb_counts) <- group_levels

  candidate_covariates <- c(
    "mean_log1p_umis", "interferon_fraction", "cytotoxic_fraction",
    "proliferation_fraction"
  )
  variable_covariates <- candidate_covariates[vapply(
    pb_meta[, ..candidate_covariates],
    function(x) is.finite(stats::sd(x)) && stats::sd(x) > 1e-8,
    logical(1L)
  )]
  design_terms <- c("subject_id", variable_covariates, "spn_group")
  design <- stats::model.matrix(
    stats::reformulate(design_terms), data = pb_meta
  )
  if (qr(design)$rank < ncol(design)) {
    design <- stats::model.matrix(~ subject_id + spn_group, data = pb_meta)
  }

  y <- edgeR::DGEList(counts = pb_counts)
  keep <- edgeR::filterByExpr(y, design = design) & rownames(y) != "SPN"
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y)
  y <- edgeR::estimateDisp(y, design, robust = use_robust_edger)
  fit <- edgeR::glmQLFit(y, design, robust = use_robust_edger)
  spn_coef <- grep("^spn_groupSPN_detected$", colnames(design))
  if (length(spn_coef) != 1L) stop("Could not identify SPN-group coefficient.")
  test <- edgeR::glmQLFTest(fit, coef = spn_coef)
  result <- as.data.table(
    edgeR::topTags(test, n = Inf, sort.by = "none")$table,
    keep.rownames = "gene"
  )
  setnames(
    result, c("logFC", "logCPM", "F", "PValue", "FDR"),
    c("log2_fold_change", "mean_log_cpm", "quasi_likelihood_F",
      "p_value", "fdr")
  )
  result[, `:=`(
    analysis = analysis_name,
    step25_base_state = state_name,
    direction = fifelse(
      log2_fold_change > 0, "Higher_in_SPN_detected",
      fifelse(log2_fold_change < 0, "Higher_in_SPN_undetected", "No_change")
    )
  )]
  setcolorder(result, c(
    "analysis", "step25_base_state", "gene", "log2_fold_change",
    "mean_log_cpm", "quasi_likelihood_F", "p_value", "fdr", "direction"
  ))
  list(
    de = result,
    qc = data.table(
      analysis = analysis_name, step25_base_state = state_name,
      n_cells = as.integer(nrow(meta)),
      n_patients = as.integer(uniqueN(meta$subject_id)),
      n_pseudobulk_samples = as.integer(ncol(pb_counts)),
      n_genes_tested = as.integer(nrow(result)),
      n_fdr_005 = as.integer(sum(result$fdr < 0.05, na.rm = TRUE))
    )
  )
}

analysis_sets <- list(
  Primary_retained_CD8 = retained_index,
  No_overlay_sensitivity = retained_index[!cell$any_overlay[retained_index]],
  Strict_direct_CD8_sensitivity = retained_index[cell$strict_cd8[retained_index]]
)

pb_results <- list()
for (analysis_name in names(analysis_sets)) {
  for (state_name in base_states) {
    key <- paste(analysis_name, state_name, sep = "||")
    pb_results[[key]] <- run_spn_pseudobulk(
      analysis_sets[[analysis_name]], state_name, analysis_name
    )
  }
}
pseudobulk_qc <- rbindlist(lapply(pb_results, `[[`, "qc"), fill = TRUE)
de_all <- rbindlist(lapply(pb_results, `[[`, "de"), fill = TRUE)
primary_de <- de_all[analysis == "Primary_retained_CD8"]
if (!nrow(primary_de)) stop("Primary within-state pseudobulk analysis failed.")

# Robust genes require primary FDR < 0.05, same direction in both available
# sensitivities, and FDR < 0.10 in at least one sensitivity.
robust_rows <- list()
for (state_name in base_states) {
  primary <- de_all[
    analysis == "Primary_retained_CD8" & step25_base_state == state_name,
    .(gene, primary_log2_fc = log2_fold_change, primary_fdr = fdr)
  ]
  no_overlay <- de_all[
    analysis == "No_overlay_sensitivity" & step25_base_state == state_name,
    .(gene, no_overlay_log2_fc = log2_fold_change, no_overlay_fdr = fdr)
  ]
  strict <- de_all[
    analysis == "Strict_direct_CD8_sensitivity" &
      step25_base_state == state_name,
    .(gene, strict_log2_fc = log2_fold_change, strict_fdr = fdr)
  ]
  concordance <- Reduce(
    function(x, y) merge(x, y, by = "gene", all = TRUE),
    list(primary, no_overlay, strict)
  )
  concordance[, same_direction_no_overlay :=
    is.na(no_overlay_log2_fc) |
      sign(primary_log2_fc) == sign(no_overlay_log2_fc)]
  concordance[, same_direction_strict :=
    is.na(strict_log2_fc) | sign(primary_log2_fc) == sign(strict_log2_fc)]
  concordance[, sensitivity_support :=
    (!is.na(no_overlay_fdr) & no_overlay_fdr < 0.10) |
      (!is.na(strict_fdr) & strict_fdr < 0.10)]
  concordance[, robust_within_state_spn_gene :=
    primary_fdr < 0.05 & same_direction_no_overlay &
      same_direction_strict & sensitivity_support]
  concordance[, step25_base_state := state_name]
  robust_rows[[state_name]] <- concordance
}
de_robustness <- rbindlist(robust_rows, fill = TRUE)

robust_gene_summary <- de_robustness[, .(
  genes_primary_fdr_005 = as.integer(sum(primary_fdr < 0.05, na.rm = TRUE)),
  genes_robust = as.integer(sum(robust_within_state_spn_gene, na.rm = TRUE)),
  robust_higher_in_spn = as.integer(sum(
    robust_within_state_spn_gene & primary_log2_fc > 0, na.rm = TRUE
  )),
  robust_lower_in_spn = as.integer(sum(
    robust_within_state_spn_gene & primary_log2_fc < 0, na.rm = TRUE
  ))
), by = step25_base_state]

# ------------------------------------------------------------------------------
# Patient-wise program associations within each state
# ------------------------------------------------------------------------------

message("Testing depth- and overlay-adjusted SPN program associations")
program_columns <- grep(
  "^step27_program_", colnames(SummarizedExperiment::colData(sce)),
  value = TRUE
)
if (!length(program_columns)) stop("No Step 27 program scores were found.")

fit_patient_program <- function(dt, program_column) {
  if (nrow(dt) < minimum_cells_per_patient_program_model ||
      sum(dt$spn_detected == 1L) < minimum_cells_per_patient_group ||
      sum(dt$spn_detected == 0L) < minimum_cells_per_patient_group) {
    return(data.table())
  }
  formula <- stats::as.formula(paste0(
    program_column,
    " ~ spn_detected + log1p_endogenous_umis + interferon_overlay + ",
    "cytotoxic_overlay + proliferation_overlay"
  ))
  fit <- tryCatch(stats::lm(formula, data = dt), error = function(e) NULL)
  if (is.null(fit) || !"spn_detected" %in% names(stats::coef(fit))) {
    return(data.table())
  }
  data.table(
    n_cells = as.integer(nrow(dt)),
    n_spn_detected = as.integer(sum(dt$spn_detected == 1L)),
    adjusted_spn_coefficient = unname(stats::coef(fit)[["spn_detected"]])
  )
}

program_model_dt <- cell[retained_index]
program_patient_rows <- list()
for (state_name in base_states) {
  for (program_column in program_columns) {
    temp <- program_model_dt[step25_base_state == state_name,
      fit_patient_program(.SD, program_column), by = subject_id]
    if (!nrow(temp)) next
    temp[, `:=`(
      step25_base_state = state_name,
      program = sub("^step27_program_", "", program_column)
    )]
    program_patient_rows[[length(program_patient_rows) + 1L]] <- temp
  }
}
program_patient_models <- rbindlist(program_patient_rows, fill = TRUE)
program_summary <- program_patient_models[is.finite(adjusted_spn_coefficient), {
  ci <- bootstrap_median_ci(adjusted_spn_coefficient, 1000L + .GRP)
  .(
    n_patients = as.integer(.N),
    median_adjusted_spn_coefficient = median(adjusted_spn_coefficient),
    bootstrap_ci_lower = unname(ci[[1L]]),
    bootstrap_ci_upper = unname(ci[[2L]]),
    percent_positive = 100 * mean(adjusted_spn_coefficient > 0),
    wilcoxon_signed_rank_p = safe_wilcox_p(adjusted_spn_coefficient)
  )
}, by = .(step25_base_state, program)]
program_summary[, fdr_within_state := p.adjust(
  wilcoxon_signed_rank_p, method = "BH"
), by = step25_base_state]

# ------------------------------------------------------------------------------
# Productive-clonotype-controlled program sensitivity
# ------------------------------------------------------------------------------

message("Testing SPN program associations within productive clonotypes")
clone_dt <- cell[retained_index][
  as_flag(tcr_productive_primary) & !is.na(tcr_clonotype_key) &
    nzchar(as.character(tcr_clonotype_key))
]
clone_dt[, clone_id := paste(subject_id, tcr_clonotype_key, sep = "::")]

fit_clone_program <- function(dt, program_column) {
  informative <- dt[, .(
    n_positive = sum(spn_detected == 1L),
    n_negative = sum(spn_detected == 0L)
  ), by = clone_id][n_positive >= 1L & n_negative >= 1L, clone_id]
  dt <- dt[clone_id %in% informative]
  if (uniqueN(dt$clone_id) < minimum_informative_clones || nrow(dt) < 20L) {
    return(data.table())
  }
  formula <- stats::as.formula(paste0(
    program_column,
    " ~ spn_detected + log1p_endogenous_umis + interferon_overlay + ",
    "cytotoxic_overlay + proliferation_overlay + factor(clone_id)"
  ))
  fit <- tryCatch(stats::lm(formula, data = dt), error = function(e) NULL)
  if (is.null(fit) || !"spn_detected" %in% names(stats::coef(fit))) {
    return(data.table())
  }
  data.table(
    n_cells = as.integer(nrow(dt)),
    n_informative_clones = as.integer(uniqueN(dt$clone_id)),
    clone_adjusted_spn_coefficient = unname(stats::coef(fit)[["spn_detected"]])
  )
}

clone_program_rows <- list()
for (state_name in base_states) {
  for (program_column in program_columns) {
    temp <- clone_dt[step25_base_state == state_name,
      fit_clone_program(.SD, program_column), by = subject_id]
    if (!nrow(temp)) next
    temp[, `:=`(
      step25_base_state = state_name,
      program = sub("^step27_program_", "", program_column)
    )]
    clone_program_rows[[length(clone_program_rows) + 1L]] <- temp
  }
}
clone_program_patient_models <- rbindlist(clone_program_rows, fill = TRUE)
if (nrow(clone_program_patient_models)) {
  clone_program_summary <- clone_program_patient_models[
    is.finite(clone_adjusted_spn_coefficient), {
      ci <- bootstrap_median_ci(clone_adjusted_spn_coefficient, 2000L + .GRP)
      .(
        n_patients = as.integer(.N),
        median_clone_adjusted_spn_coefficient = median(
          clone_adjusted_spn_coefficient
        ),
        bootstrap_ci_lower = unname(ci[[1L]]),
        bootstrap_ci_upper = unname(ci[[2L]]),
        percent_positive = 100 * mean(clone_adjusted_spn_coefficient > 0),
        wilcoxon_signed_rank_p = safe_wilcox_p(
          clone_adjusted_spn_coefficient
        )
      )
    }, by = .(step25_base_state, program)
  ]
  clone_program_summary[, fdr_within_state := p.adjust(
    wilcoxon_signed_rank_p, method = "BH"
  ), by = step25_base_state]
} else {
  clone_program_summary <- data.table()
}

# ------------------------------------------------------------------------------
# Output tables and figures
# ------------------------------------------------------------------------------

message("Writing Step 28 tables")
fwrite(pseudobulk_qc, file.path(results_dir, "28_pseudobulk_qc.csv"))
fwrite(de_all, file.path(results_dir, "28_within_state_spn_pseudobulk_de.csv"))
fwrite(de_robustness, file.path(results_dir, "28_spn_gene_robustness.csv"))
fwrite(robust_gene_summary, file.path(results_dir, "28_robust_gene_summary.csv"))
fwrite(program_patient_models,
       file.path(results_dir, "28_spn_program_patient_models.csv"))
fwrite(program_summary,
       file.path(results_dir, "28_spn_program_summary.csv"))
fwrite(clone_program_patient_models,
       file.path(results_dir, "28_clone_program_patient_models.csv"))
fwrite(clone_program_summary,
       file.path(results_dir, "28_clone_program_summary.csv"))

primary_qc <- pseudobulk_qc[analysis == "Primary_retained_CD8"]
dataset_qc <- data.table(
  metric = c(
    "input_cells", "retained_cd8_cells", "subjects",
    "primary_memory_patients", "primary_dysfunctional_patients",
    "primary_memory_genes_fdr_005", "primary_dysfunctional_genes_fdr_005",
    "robust_memory_spn_genes", "robust_dysfunctional_spn_genes",
    "programs_tested", "clone_program_tests", "spn_excluded_from_discovery",
    "state_accounting_validation", "signature_exported"
  ),
  value = c(
    ncol(sce), length(retained_index), uniqueN(cell$subject_id[retained_index]),
    primary_qc[step25_base_state == "CD8_memory_transitional", n_patients],
    primary_qc[step25_base_state == "CD8_dysfunctional", n_patients],
    primary_qc[step25_base_state == "CD8_memory_transitional", n_fdr_005],
    primary_qc[step25_base_state == "CD8_dysfunctional", n_fdr_005],
    robust_gene_summary[
      step25_base_state == "CD8_memory_transitional", genes_robust
    ],
    robust_gene_summary[
      step25_base_state == "CD8_dysfunctional", genes_robust
    ],
    nrow(program_summary), nrow(clone_program_summary), 1L,
    as.integer(length(retained_index) == sum(
      cell$step25_base_state[retained_index] %in% base_states
    )),
    0L
  )
)
fwrite(dataset_qc, file.path(results_dir, "28_dataset_qc.csv"))

plot_de <- primary_de[
  is.finite(log2_fold_change) & is.finite(fdr)
]
plot_de[, minus_log10_fdr := -log10(pmax(fdr, 1e-300))]
plot_de[, significant := fdr < 0.05]
plot_de[, state_label := state_labels[step25_base_state]]
p_volcano <- ggplot(
  plot_de,
  aes(log2_fold_change, minus_log10_fdr, color = significant)
) +
  geom_point(alpha = 0.55, size = 1.0) +
  geom_vline(xintercept = 0, color = "grey55") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey55") +
  facet_wrap(~ state_label, scales = "free_y") +
  scale_color_manual(values = c(`FALSE` = "#BDBDBD", `TRUE` = "#7651B9")) +
  labs(
    title = "Within-state transcriptional associations with SPN/CD43",
    subtitle = "Positive effects indicate higher expression in SPN-detected cells; SPN itself is excluded",
    x = "SPN-detected versus undetected log2 fold-change",
    y = expression(-log[10](FDR)), color = "FDR < 0.05"
  ) + theme_project()
save_png("28A_within_state_spn_volcano.png", p_volcano, 13, 6)

program_plot <- copy(program_summary)
program_plot[, state_label := state_labels[step25_base_state]]
p_program <- ggplot(
  program_plot,
  aes(
    median_adjusted_spn_coefficient, reorder(program, median_adjusted_spn_coefficient),
    color = step25_base_state
  )
) +
  geom_vline(xintercept = 0, color = "grey55") +
  geom_errorbarh(
    aes(xmin = bootstrap_ci_lower, xmax = bootstrap_ci_upper), height = 0.15
  ) +
  geom_point(size = 2.6) +
  facet_wrap(~ state_label, scales = "free_y") +
  scale_color_manual(values = state_colors, guide = "none") +
  labs(
    title = "Programs associated with SPN within each CD8 state",
    subtitle = "Patient-wise coefficients adjust for sequencing depth and state overlays",
    x = "Median adjusted SPN coefficient", y = NULL
  ) + theme_project()
save_png("28B_within_state_spn_programs.png", p_program, 13, 8)

top_robust <- de_robustness[robust_within_state_spn_gene == TRUE]
top_robust <- top_robust[order(step25_base_state, primary_fdr, -abs(primary_log2_fc))]
top_robust <- top_robust[, head(.SD, 20L), by = step25_base_state]
if (nrow(top_robust)) {
  top_robust[, state_label := state_labels[step25_base_state]]
  top_robust[, gene_label := factor(
    gene, levels = rev(unique(gene[order(primary_log2_fc)]))
  )]
  p_genes <- ggplot(
    top_robust,
    aes(primary_log2_fc, gene_label, color = step25_base_state)
  ) +
    geom_vline(xintercept = 0, color = "grey55") +
    geom_point(size = 2.5) +
    facet_wrap(~ state_label, scales = "free_y") +
    scale_color_manual(values = state_colors, guide = "none") +
    labs(
      title = "Robust genes associated with SPN within CD8 states",
      subtitle = "Primary FDR < 0.05 with directionally concordant sensitivity support",
      x = "SPN-detected versus undetected log2 fold-change", y = NULL
    ) + theme_project()
  save_png("28C_robust_within_state_spn_genes.png", p_genes, 12, 9)
}

# Store compact cell-level annotations; do not promote any Step 28 result to a
# deconvolution signature without independent validation.
sce$step28_analysis_eligible <- rep(FALSE, ncol(sce))
sce$step28_analysis_eligible[retained_index] <- TRUE
sce$step28_spn_group <- rep(NA_character_, ncol(sce))
sce$step28_spn_group[retained_index] <- as.character(cell$spn_group[retained_index])

message("Saving Step 28 object")
saveRDS(sce, output_sce_file, compress = FALSE)
dataset_qc[metric == "signature_exported", value := 0L]
dataset_qc <- rbind(
  dataset_qc,
  data.table(
    metric = "output_sce_size_mb",
    value = as.numeric(file.info(output_sce_file)$size / 1024^2)
  )
)
fwrite(dataset_qc, file.path(results_dir, "28_dataset_qc.csv"))

message(
  "Step 28 complete: within-state SPN biology tested with patient, depth, ",
  "overlay, strict-CD8, and clonotype-aware controls"
)

