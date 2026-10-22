#!/usr/bin/env Rscript

# ==============================================================================
# Step 37: Validate the Step 36 reference with patient-stratified pseudo-mixtures
# ==============================================================================
# The fixed SPN-free reference is applied to synthetic mixtures assembled from
# individual patients. True proportions are known. Validation emphasizes total
# CD8 abundance and the dysfunctional fraction within CD8, which are the two
# quantities needed for the final TCGA-SKCM integration.
# ============================================================================== 

options(stringsAsFactors = FALSE)
options(error = function() {
  cat("\n--- STEP 37 TRACEBACK ---\n", file = stderr())
  traceback(20L)
  quit(save = "no", status = 1L, runLast = FALSE)
})

required_packages <- c(
  "data.table", "SingleCellExperiment", "SummarizedExperiment", "Matrix",
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
results_dir <- file.path(project_root, "results", "37")
figures_dir <- file.path(results_dir, "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

full_sce_file <- file.path(analysis_dir, "19_GSE123139_broad_clusters_sce.rds")
cd8_sce_file <- file.path(
  analysis_dir, "34_GSE123139_beltra_tex_states_spn_sce.rds"
)
signature_file <- file.path(
  project_root, "results", "36", "cibersortx",
  "36_cibersortx_signature_matrix.tsv"
)
for (path in c(full_sce_file, cd8_sce_file, signature_file)) {
  if (!file.exists(path)) stop("Required input not found: ", path)
}

random_seed <- 123139L
mixtures_per_patient <- 50L
cells_per_mixture <- 400L
minimum_cells_per_cd8_state <- 20L
minimum_reference_classes <- 3L

broad_compartment_map <- c(
  C01 = "T_or_NK", C02 = "T_or_NK", C03 = "Myeloid",
  C04 = "B_cell", C05 = "T_or_NK", C06 = "T_or_NK",
  C07 = "Myeloid", C08 = "Melanoma", C09 = "T_or_NK",
  C10 = "Other_low_information", C11 = "T_or_NK",
  C12 = "Other_low_information", C13 = "Other_low_information",
  C14 = "T_or_NK"
)
reference_levels <- c(
  "CD8_memory_transitional", "CD8_dysfunctional", "Other_T_NK",
  "B_cell", "Myeloid", "Melanoma"
)

safe_spearman <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 5L || length(unique(x[keep])) < 2L ||
      length(unique(y[keep])) < 2L) return(NA_real_)
  suppressWarnings(stats::cor(x[keep], y[keep], method = "spearman"))
}

fit_nonnegative <- function(reference, mixture) {
  objective <- function(beta) {
    residual <- as.numeric(mixture - reference %*% beta)
    sum(residual^2)
  }
  fit <- tryCatch(
    optim(
      rep(1 / ncol(reference), ncol(reference)), objective,
      method = "L-BFGS-B", lower = rep(0, ncol(reference)),
      control = list(maxit = 1000L)
    ), error = function(e) NULL
  )
  if (is.null(fit) || !is.finite(sum(fit$par)) || sum(fit$par) <= 0) {
    return(rep(NA_real_, ncol(reference)))
  }
  fit$par / sum(fit$par)
}

message("Loading atlas, CD8 labels, and Step 36 reference")
full_sce <- readRDS(full_sce_file)
cd8_sce <- readRDS(cd8_sce_file)
signature <- fread(signature_file)

if (!"GeneSymbol" %in% names(signature)) stop("Signature lacks GeneSymbol column.")
if (!all(reference_levels %in% names(signature))) {
  stop("Signature is missing expected reference classes.")
}
if ("SPN" %chin% signature$GeneSymbol) stop("SPN occurs in validation reference.")

genes <- signature$GeneSymbol
if (!all(genes %chin% rownames(full_sce))) stop("Signature genes missing from atlas.")
reference <- as.matrix(signature[, ..reference_levels])
storage.mode(reference) <- "double"

# Equalize gene influence while preserving between-class contrasts.
gene_scale <- apply(reference, 1L, max, na.rm = TRUE)
gene_scale[!is.finite(gene_scale) | gene_scale <= 0] <- 1
reference_scaled <- reference / gene_scale

full_dt <- as.data.table(SummarizedExperiment::colData(full_sce))
full_dt[, cell_id := colnames(full_sce)]
full_dt[, subject_id := as.character(subject_id)]
full_dt[, source_group := as.character(source_group)]
full_dt[, broad_cluster := as.character(broad_cluster)]
full_dt[, broad_compartment := unname(broad_compartment_map[broad_cluster])]
full_dt[, reference_class := fifelse(
  broad_compartment == "T_or_NK", "Other_T_NK", broad_compartment
)]

cd8_dt <- as.data.table(SummarizedExperiment::colData(cd8_sce))
cd8_dt[, cell_id := colnames(cd8_sce)]
cd8_dt[, base_state := as.character(step25_base_state)]
cd8_dt[, retained := !is.na(step25_retained_cd8) & as.logical(step25_retained_cd8)]
cd8_map <- cd8_dt[
  retained == TRUE & base_state %chin%
    c("CD8_memory_transitional", "CD8_dysfunctional"),
  .(cell_id, cd8_reference_class = base_state)
]
full_dt <- merge(full_dt, cd8_map, by = "cell_id", all.x = TRUE, sort = FALSE)
full_dt[!is.na(cd8_reference_class), reference_class := cd8_reference_class]
full_dt <- full_dt[
  source_group == "Tumor" & reference_class %chin% reference_levels
]
full_dt[, atlas_index := match(cell_id, colnames(full_sce))]
if (anyNA(full_dt$atlas_index)) stop("Atlas cell alignment failed.")

class_counts <- full_dt[, .N, by = .(subject_id, reference_class)]
eligible_patients <- class_counts[, .(
  n_classes = sum(N > 0),
  memory_cells = sum(N[reference_class == "CD8_memory_transitional"]),
  dysfunctional_cells = sum(N[reference_class == "CD8_dysfunctional"])
), by = subject_id][
  n_classes >= minimum_reference_classes &
    memory_cells >= minimum_cells_per_cd8_state &
    dysfunctional_cells >= minimum_cells_per_cd8_state
]
if (!nrow(eligible_patients)) stop("No patients qualify for pseudo-mixtures.")

logcounts <- SummarizedExperiment::assay(full_sce, "logcounts")
linear_expression <- pmax(
  2^as.matrix(logcounts[genes, full_dt$atlas_index, drop = FALSE]) - 1, 0
)

set.seed(random_seed)
mixture_rows <- list()
row_counter <- 0L
message("Generating patient-stratified pseudo-mixtures")
for (patient in eligible_patients$subject_id) {
  patient_rows <- which(full_dt$subject_id == patient)
  available <- reference_levels[
    reference_levels %chin% unique(full_dt$reference_class[patient_rows])
  ]
  class_cell_indices <- lapply(available, function(class_name) {
    patient_rows[full_dt$reference_class[patient_rows] == class_name]
  })
  names(class_cell_indices) <- available

  for (iteration in seq_len(mixtures_per_patient)) {
    raw_weights <- rgamma(length(available), shape = 1.5, rate = 1)
    weights <- raw_weights / sum(raw_weights)
    sampled_counts <- as.integer(rmultinom(1L, cells_per_mixture, weights))
    sampled_atlas_columns <- unlist(lapply(seq_along(available), function(j) {
      sample(class_cell_indices[[j]], sampled_counts[j], replace = TRUE)
    }), use.names = FALSE)
    mixture_expression <- rowMeans(
      linear_expression[, sampled_atlas_columns, drop = FALSE]
    )
    estimated <- fit_nonnegative(
      reference_scaled, mixture_expression / gene_scale
    )
    names(estimated) <- reference_levels
    truth <- setNames(rep(0, length(reference_levels)), reference_levels)
    truth[available] <- sampled_counts / sum(sampled_counts)

    row_counter <- row_counter + 1L
    result <- data.table(
      subject_id = patient, mixture_iteration = iteration,
      available_classes = length(available)
    )
    for (class_name in reference_levels) {
      result[, (paste0("true_", class_name)) := truth[class_name]]
      result[, (paste0("estimated_", class_name)) := estimated[class_name]]
    }
    true_total_cd8 <- truth["CD8_memory_transitional"] +
      truth["CD8_dysfunctional"]
    estimated_total_cd8 <- estimated["CD8_memory_transitional"] +
      estimated["CD8_dysfunctional"]
    result[, `:=`(
      true_total_cd8 = true_total_cd8,
      estimated_total_cd8 = estimated_total_cd8,
      true_dysfunctional_within_cd8 = if (true_total_cd8 > 0)
        truth["CD8_dysfunctional"] / true_total_cd8 else NA_real_,
      estimated_dysfunctional_within_cd8 = if (estimated_total_cd8 > 0)
        estimated["CD8_dysfunctional"] / estimated_total_cd8 else NA_real_
    )]
    mixture_rows[[row_counter]] <- result
  }
}
mixtures <- rbindlist(mixture_rows, fill = TRUE)

targets <- c(reference_levels, "total_cd8", "dysfunctional_within_cd8")
metric_rows <- rbindlist(lapply(targets, function(target) {
  true_column <- paste0("true_", target)
  estimated_column <- paste0("estimated_", target)
  mixtures[, .(
    spearman_rho = safe_spearman(get(true_column), get(estimated_column)),
    rmse = sqrt(mean((get(estimated_column) - get(true_column))^2,
                     na.rm = TRUE)),
    mae = mean(abs(get(estimated_column) - get(true_column)), na.rm = TRUE),
    median_bias = median(get(estimated_column) - get(true_column), na.rm = TRUE)
  ), by = subject_id][, target := target]
}), use.names = TRUE)

metric_summary <- metric_rows[, .(
  n_patients = sum(is.finite(spearman_rho)),
  median_patient_spearman_rho = median(spearman_rho, na.rm = TRUE),
  q25_patient_spearman_rho = quantile(spearman_rho, 0.25, na.rm = TRUE),
  q75_patient_spearman_rho = quantile(spearman_rho, 0.75, na.rm = TRUE),
  median_patient_rmse = median(rmse, na.rm = TRUE),
  median_patient_mae = median(mae, na.rm = TRUE),
  median_patient_bias = median(median_bias, na.rm = TRUE)
), by = target]

metric_summary[, validation_class := fifelse(
  target == "total_cd8" & median_patient_spearman_rho >= 0.70 &
    median_patient_mae <= 0.10,
  "Pass",
  fifelse(
    target == "dysfunctional_within_cd8" &
      median_patient_spearman_rho >= 0.60 & median_patient_mae <= 0.15,
    "Pass",
    fifelse(target %chin% c("total_cd8", "dysfunctional_within_cd8"),
            "Fail", "Descriptive")
  )
)]

dataset_qc <- data.table(
  metric = c(
    "signature_genes", "eligible_patients", "mixtures_per_patient",
    "total_pseudomixtures", "cells_per_mixture", "spn_in_reference",
    "total_cd8_validation_pass", "within_cd8_state_validation_pass"
  ),
  value = c(
    length(genes), nrow(eligible_patients), mixtures_per_patient,
    nrow(mixtures), cells_per_mixture, 0L,
    as.integer(metric_summary[target == "total_cd8"]$validation_class == "Pass"),
    as.integer(metric_summary[
      target == "dysfunctional_within_cd8"
    ]$validation_class == "Pass")
  )
)

message("Writing Step 37 validation outputs")
fwrite(dataset_qc, file.path(results_dir, "37_dataset_qc.csv"))
fwrite(eligible_patients, file.path(results_dir, "37_eligible_patients.csv"))
fwrite(mixtures, file.path(results_dir, "37_pseudomixture_predictions.tsv.gz"),
       sep = "\t")
fwrite(metric_rows, file.path(results_dir, "37_patient_validation_metrics.csv"))
fwrite(metric_summary, file.path(results_dir, "37_validation_summary.csv"))

plot_targets <- c("total_cd8", "dysfunctional_within_cd8")
plot_dt <- rbindlist(lapply(plot_targets, function(target) {
  mixtures[, .(
    subject_id,
    truth = get(paste0("true_", target)),
    estimate = get(paste0("estimated_", target)),
    target = target
  )]
}))
p <- ggplot(plot_dt, aes(truth, estimate)) +
  geom_abline(slope = 1, intercept = 0, color = "#777777") +
  geom_point(aes(color = subject_id), alpha = 0.35, size = 1.1) +
  facet_wrap(~target, scales = "free", ncol = 2) +
  coord_equal() +
  labs(
    title = "Patient-stratified pseudo-mixture validation",
    subtitle = "Known proportions compared with Step 36 reference estimates",
    x = "True fraction", y = "Estimated fraction", color = "Patient"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"), legend.position = "none",
    strip.background = element_blank(), strip.text = element_text(face = "bold")
  )
ggsave(
  file.path(figures_dir, "37A_pseudomixture_truth_vs_estimate.png"),
  p, width = 11, height = 6.5, units = "in", dpi = 180, bg = "white"
)

message("Step 37 complete: deconvolution reference validated on pseudo-mixtures")
