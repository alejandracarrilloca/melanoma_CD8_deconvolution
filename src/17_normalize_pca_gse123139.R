#!/usr/bin/env Rscript

# ==============================================================================
# Step 17: Normalize GSE123139 and assess uncorrected technical structure
# ==============================================================================
#
# This step preserves the Step 16 raw counts, estimates scran deconvolution
# size factors, adds log-normalized expression, selects 2,000 variable genes,
# and calculates an uncorrected PCA. No batch correction or cell-type labeling
# is performed here. The uncorrected representation is required to determine
# whether plate, amplification, FACS gate, source, or patient structure drives
# the leading components before choosing any correction strategy.
#
# Inputs:
#   data/analysis_ready/scRNA/GSE123139/
#     16_GSE123139_endogenous_qc_sce.rds
#     14_GSE123139_endogenous_cell_qc.tsv.gz
#
# Outputs:
#   data/analysis_ready/scRNA/GSE123139/
#     17_GSE123139_normalized_pca_sce.rds
#   results/17/
#     17_dataset_qc.csv
#     17_normalization_summary.csv
#     17_retention_by_specimen.csv
#     17_missing_specimens.csv
#     17_highly_variable_genes.csv
#     17_pca_variance.csv
#     17_pc_metadata_associations.csv
#     17_figure_manifest.csv
#     17_sessionInfo.txt
#     figures/*.png
#
# Run from the project root:
#   Rscript src/17_normalize_pca_gse123139.R
# ==============================================================================

options(stringsAsFactors = FALSE)

required_packages <- c(
  "data.table", "Matrix", "SingleCellExperiment", "SummarizedExperiment",
  "S4Vectors", "scuttle", "scran", "BiocSingular", "ggplot2", "scales"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Install them in the existing single-cell environment before rerunning."
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ------------------------------------------------------------------------------
# Paths and configuration
# ------------------------------------------------------------------------------

project_root <- normalizePath(getwd(), mustWork = TRUE)
analysis_dir <- file.path(
  project_root, "data", "analysis_ready", "scRNA", "GSE123139"
)
results_dir <- file.path(project_root, "results", "17")
figures_dir <- file.path(results_dir, "figures")

dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

input_sce_file <- file.path(
  analysis_dir, "16_GSE123139_endogenous_qc_sce.rds"
)
step14_qc_file <- file.path(
  analysis_dir, "14_GSE123139_endogenous_cell_qc.tsv.gz"
)
output_sce_file <- file.path(
  analysis_dir, "17_GSE123139_normalized_pca_sce.rds"
)

n_hvgs <- 2000L
n_pcs <- 50L
minimum_cells_for_hvg <- 50L
maximum_plot_cells <- 40000L

for (path in c(input_sce_file, step14_qc_file)) {
  if (!file.exists(path)) stop("Missing required input: ", path)
}

# ------------------------------------------------------------------------------
# Load and validate the final raw-count object
# ------------------------------------------------------------------------------

message("Loading Step 16 SingleCellExperiment")
sce <- readRDS(input_sce_file)
counts <- SummarizedExperiment::assay(sce, "counts")

if (nrow(sce) != 55765L || ncol(sce) != 52943L) {
  stop(
    "Unexpected Step 16 dimensions: ", nrow(sce), " genes x ",
    ncol(sce), " cells"
  )
}
if (anyDuplicated(rownames(sce)) || anyDuplicated(colnames(sce))) {
  stop("Duplicated gene or cell identifiers occur in the Step 16 object.")
}
if (any(grepl("^ERCC-", rownames(sce), ignore.case = TRUE))) {
  stop("ERCC rows unexpectedly occur in the Step 16 object.")
}

required_metadata <- c(
  "source_group", "subject_id", "specimen_id", "facs_gate",
  "amplification_batch", "plate_configuration", "endogenous_n_counts",
  "endogenous_n_genes"
)
missing_metadata <- setdiff(
  required_metadata,
  colnames(SummarizedExperiment::colData(sce))
)
if (length(missing_metadata) > 0L) {
  stop("Step 16 metadata is missing: ", paste(missing_metadata, collapse = ", "))
}

# ------------------------------------------------------------------------------
# Document specimen retention, including the specimen absent from Step 16
# ------------------------------------------------------------------------------

qc14 <- fread(step14_qc_file, check.names = FALSE)
qc14[, final_qc_pass_800 :=
  endogenous_n_counts >= 800L &
  endogenous_n_genes >= 200L &
  !is.na(endogenous_mitochondrial_percent) &
  endogenous_mitochondrial_percent <= 20 &
  !endogenous_high_library_outlier
]

retention_by_specimen <- qc14[, .(
  source_group = unique(source_group)[1L],
  subject_id = unique(subject_id)[1L],
  n_input_wells = .N,
  n_final_cells = sum(final_qc_pass_800),
  percent_final = 100 * mean(final_qc_pass_800)
), by = specimen_id][order(source_group, subject_id, specimen_id)]

missing_specimens <- retention_by_specimen[n_final_cells == 0L]

fwrite(
  retention_by_specimen,
  file.path(results_dir, "17_retention_by_specimen.csv")
)
fwrite(missing_specimens, file.path(results_dir, "17_missing_specimens.csv"))

# ------------------------------------------------------------------------------
# scran deconvolution normalization
# ------------------------------------------------------------------------------

message("Computing quick clusters for deconvolution normalization")
set.seed(123139)
normalization_clusters <- scran::quickCluster(sce, min.size = 100L)
sce$normalization_cluster <- as.character(normalization_clusters)

message("Estimating scran deconvolution size factors")
sce <- scran::computeSumFactors(sce, clusters = normalization_clusters)
size_factors <- as.numeric(
  SummarizedExperiment::colData(sce)$sizeFactor
)

if (any(!is.finite(size_factors)) || any(size_factors <= 0)) {
  stop("Invalid non-positive or non-finite scran size factors were produced.")
}

message("Calculating log-normalized expression")
sce <- scuttle::logNormCounts(sce)

normalization_summary <- as.data.table(as.data.frame(
  SummarizedExperiment::colData(sce)
))[, .(
  n_cells = .N,
  median_size_factor = median(sizeFactor),
  minimum_size_factor = min(sizeFactor),
  maximum_size_factor = max(sizeFactor),
  median_endogenous_umis = median(endogenous_n_counts),
  median_endogenous_genes = median(endogenous_n_genes)
), by = source_group][order(source_group)]

fwrite(
  normalization_summary,
  file.path(results_dir, "17_normalization_summary.csv")
)

# ------------------------------------------------------------------------------
# Model gene variance and select HVGs for an uncorrected PCA
# ------------------------------------------------------------------------------

message("Modeling gene-level variance")
variance_fit <- scran::modelGeneVar(sce)

detected_cells <- as.integer(Matrix::rowSums(counts > 0))
technical_gene <-
  grepl("^MT-", rownames(sce)) |
  grepl("^RP[SL][0-9]", rownames(sce)) |
  grepl("^ERCC-", rownames(sce), ignore.case = TRUE) |
  rownames(sce) %in% c("MALAT1", "HBA1", "HBA2", "HBB")

eligible_hvg <- detected_cells >= minimum_cells_for_hvg & !technical_gene

if (sum(eligible_hvg) < n_hvgs) {
  stop("Fewer than ", n_hvgs, " genes are eligible for HVG selection.")
}

hvg_genes <- scran::getTopHVGs(
  variance_fit[eligible_hvg, ],
  n = n_hvgs
)

if (length(hvg_genes) != n_hvgs || anyDuplicated(hvg_genes)) {
  stop("HVG selection did not return exactly ", n_hvgs, " unique genes.")
}

hvg_table <- data.table(
  gene_symbol = rownames(sce),
  n_cells_detected = detected_cells,
  mean = variance_fit$mean,
  total_variance = variance_fit$total,
  technical_variance = variance_fit$tech,
  biological_variance = variance_fit$bio,
  p_value = variance_fit$p.value,
  fdr = variance_fit$FDR,
  excluded_technical_gene = technical_gene,
  eligible_for_hvg = eligible_hvg,
  selected_hvg = rownames(sce) %chin% hvg_genes
)[order(-biological_variance)]

fwrite(
  hvg_table,
  file.path(results_dir, "17_highly_variable_genes.csv")
)

row_data <- SummarizedExperiment::rowData(sce)
row_data$n_cells_detected <- detected_cells
row_data$excluded_from_hvg_selection <- technical_gene
row_data$highly_variable <- rownames(sce) %chin% hvg_genes
row_data$biological_variance <- variance_fit$bio
SummarizedExperiment::rowData(sce) <- row_data

# ------------------------------------------------------------------------------
# Uncorrected PCA
# ------------------------------------------------------------------------------

message("Running uncorrected PCA on ", n_hvgs, " HVGs")
logcounts_hvg <- SummarizedExperiment::assay(sce, "logcounts")[
  hvg_genes, , drop = FALSE
]

message(
  "HVG logcounts class before PCA: ",
  paste(class(logcounts_hvg), collapse = ", ")
)

# Bioconductor assays can be returned as sparse or delayed matrix-like objects.
# Coerce explicitly to a compressed sparse matrix before transposing so that
# base::t.default() is never selected and the PCA input is not made dense.
logcounts_hvg_sparse <- tryCatch(
  methods::as(logcounts_hvg, "dgCMatrix"),
  error = function(e) {
    stop(
      "Could not coerce the HVG logcounts assay to dgCMatrix for PCA: ",
      conditionMessage(e)
    )
  }
)
pca_input <- Matrix::t(logcounts_hvg_sparse)

if (!inherits(pca_input, "Matrix") ||
    nrow(pca_input) != ncol(sce) ||
    ncol(pca_input) != n_hvgs) {
  stop("Unexpected sparse PCA input class or dimensions.")
}

set.seed(123139)
pca_result <- BiocSingular::runPCA(
  pca_input,
  rank = n_pcs,
  BSPARAM = BiocSingular::IrlbaParam()
)

rm(logcounts_hvg, logcounts_hvg_sparse, pca_input)
invisible(gc())

if (nrow(pca_result$x) != ncol(sce) || ncol(pca_result$x) != n_pcs) {
  stop("Unexpected PCA dimensions.")
}

rownames(pca_result$x) <- colnames(sce)
colnames(pca_result$x) <- paste0("PC", seq_len(n_pcs))
SingleCellExperiment::reducedDim(sce, "PCA") <- pca_result$x

total_hvg_variance <- sum(
  variance_fit$total[match(hvg_genes, rownames(variance_fit))],
  na.rm = TRUE
)
pca_variance <- data.table(
  pc = paste0("PC", seq_len(n_pcs)),
  pc_number = seq_len(n_pcs),
  singular_value = pca_result$d,
  variance = pca_result$d^2 / (ncol(sce) - 1),
  percent_hvg_variance = 100 *
    (pca_result$d^2 / (ncol(sce) - 1)) / total_hvg_variance
)
pca_variance[, cumulative_percent_hvg_variance :=
  cumsum(percent_hvg_variance)
]

fwrite(pca_variance, file.path(results_dir, "17_pca_variance.csv"))

# ------------------------------------------------------------------------------
# Descriptive PC associations with biological and technical metadata
# ------------------------------------------------------------------------------

metadata <- as.data.table(as.data.frame(
  SummarizedExperiment::colData(sce)
))
pca_scores <- as.data.table(pca_result$x)

categorical_variables <- c(
  "source_group", "subject_id", "specimen_id", "facs_gate",
  "amplification_batch", "plate_configuration", "normalization_cluster"
)
numeric_variables <- c(
  "endogenous_n_counts", "endogenous_n_genes", "sizeFactor"
)

categorical_eta_squared <- function(values, scores) {
  valid <- !is.na(values) & nzchar(as.character(values)) & is.finite(scores)
  values <- as.character(values[valid])
  scores <- scores[valid]
  if (length(scores) < 2L || length(unique(values)) < 2L) return(NA_real_)

  grand_mean <- mean(scores)
  total_ss <- sum((scores - grand_mean)^2)
  if (total_ss <= 0) return(NA_real_)

  group_summary <- data.table(group = values, score = scores)[, .(
    n = .N,
    group_mean = mean(score)
  ), by = group]
  between_ss <- group_summary[, sum(n * (group_mean - grand_mean)^2)]
  between_ss / total_ss
}

association_list <- list()
association_counter <- 0L

for (pc_index in seq_len(min(20L, n_pcs))) {
  pc_name <- paste0("PC", pc_index)
  scores <- pca_scores[[pc_name]]

  for (variable in categorical_variables) {
    association_counter <- association_counter + 1L
    association_list[[association_counter]] <- data.table(
      pc = pc_name,
      pc_number = pc_index,
      variable = variable,
      variable_type = "categorical_eta_squared",
      association_r2 = categorical_eta_squared(metadata[[variable]], scores),
      n_levels = uniqueN(metadata[[variable]], na.rm = TRUE)
    )
  }

  for (variable in numeric_variables) {
    values <- metadata[[variable]]
    valid <- is.finite(values) & is.finite(scores)
    correlation <- suppressWarnings(cor(
      values[valid], scores[valid], method = "spearman"
    ))

    association_counter <- association_counter + 1L
    association_list[[association_counter]] <- data.table(
      pc = pc_name,
      pc_number = pc_index,
      variable = variable,
      variable_type = "squared_spearman_correlation",
      association_r2 = correlation^2,
      n_levels = NA_integer_
    )
  }
}

pc_associations <- rbindlist(association_list)
fwrite(
  pc_associations,
  file.path(results_dir, "17_pc_metadata_associations.csv")
)

# ------------------------------------------------------------------------------
# Diagnostic figures
# ------------------------------------------------------------------------------

lavender <- "#B993F6"
lavender_dark <- "#7651B5"
gray_fill <- "#EEEEEE"
gray_outline <- "#949494"

theme_project <- theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(color = "#555555", size = 11),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

plot_metadata <- cbind(metadata, pca_scores[, .(PC1, PC2, PC3, PC4)])
set.seed(123139)
plot_metadata <- if (nrow(plot_metadata) > maximum_plot_cells) {
  plot_metadata[sample(.N, maximum_plot_cells)]
} else {
  copy(plot_metadata)
}

plot_17a <- ggplot(
  metadata,
  aes(x = source_group, y = sizeFactor, fill = source_group)
) +
  geom_violin(scale = "width", color = gray_outline, linewidth = 0.35) +
  geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white", linewidth = 0.4) +
  scale_fill_manual(values = c("Tumor" = lavender, "PBMC" = gray_fill)) +
  labs(
    title = "scran size-factor distributions",
    subtitle = "Deconvolution normalization was estimated across quick clusters",
    x = NULL,
    y = "Size factor",
    fill = "Sample source"
  ) +
  theme_project

plot_17b <- ggplot(
  plot_metadata,
  aes(x = PC1, y = PC2, fill = source_group)
) +
  geom_point(
    shape = 21, size = 1.3, stroke = 0.2,
    color = gray_outline, alpha = 0.55
  ) +
  scale_fill_manual(values = c("Tumor" = lavender, "PBMC" = gray_fill)) +
  labs(
    title = "Uncorrected PCA by sample source",
    subtitle = "No patient, plate, or batch correction has been applied",
    x = paste0("PC1 (", round(pca_variance$percent_hvg_variance[1L], 1), "%)"),
    y = paste0("PC2 (", round(pca_variance$percent_hvg_variance[2L], 1), "%)"),
    fill = "Sample source"
  ) +
  theme_project

plot_17c <- ggplot(
  plot_metadata,
  aes(x = PC1, y = PC2, fill = facs_gate)
) +
  geom_point(
    shape = 21, size = 1.3, stroke = 0.2,
    color = gray_outline, alpha = 0.55
  ) +
  labs(
    title = "Uncorrected PCA by FACS gate",
    subtitle = "FACS enrichment may reflect both technical design and real lineage composition",
    x = paste0("PC1 (", round(pca_variance$percent_hvg_variance[1L], 1), "%)"),
    y = paste0("PC2 (", round(pca_variance$percent_hvg_variance[2L], 1), "%)"),
    fill = "FACS gate"
  ) +
  theme_project

plot_17d <- ggplot(
  plot_metadata,
  aes(x = PC1, y = PC2, fill = plate_configuration)
) +
  geom_point(
    shape = 21, size = 1.3, stroke = 0.2,
    color = gray_outline, alpha = 0.55
  ) +
  labs(
    title = "Uncorrected PCA by plate configuration",
    subtitle = "ERCC rows were removed before normalization",
    x = paste0("PC1 (", round(pca_variance$percent_hvg_variance[1L], 1), "%)"),
    y = paste0("PC2 (", round(pca_variance$percent_hvg_variance[2L], 1), "%)"),
    fill = "Plate configuration"
  ) +
  theme_project

heatmap_data <- pc_associations[
  pc_number <= 10L & variable != "normalization_cluster"
]
heatmap_data[, variable := factor(
  variable,
  levels = rev(c(
    "source_group", "subject_id", "specimen_id", "facs_gate",
    "amplification_batch", "plate_configuration", "endogenous_n_counts",
    "endogenous_n_genes", "sizeFactor"
  ))
)]

plot_17e <- ggplot(
  heatmap_data,
  aes(x = factor(pc, levels = paste0("PC", 1:10)), y = variable, fill = association_r2)
) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient(low = gray_fill, high = lavender_dark, limits = c(0, 1)) +
  labs(
    title = "Descriptive associations between PCs and metadata",
    subtitle = "Categorical variables use eta-squared; numeric variables use squared Spearman correlation",
    x = NULL,
    y = NULL,
    fill = expression(R^2)
  ) +
  theme_project +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

figure_paths <- file.path(
  figures_dir,
  c(
    "17A_size_factor_distributions.png",
    "17B_uncorrected_pca_by_source.png",
    "17C_uncorrected_pca_by_facs_gate.png",
    "17D_uncorrected_pca_by_plate_configuration.png",
    "17E_pc_metadata_associations.png"
  )
)

ggsave(figure_paths[1L], plot_17a, width = 8, height = 6, dpi = 300, bg = "white")
ggsave(figure_paths[2L], plot_17b, width = 9, height = 7, dpi = 300, bg = "white")
ggsave(figure_paths[3L], plot_17c, width = 9, height = 7, dpi = 300, bg = "white")
ggsave(figure_paths[4L], plot_17d, width = 10, height = 7, dpi = 300, bg = "white")
ggsave(figure_paths[5L], plot_17e, width = 10, height = 6.5, dpi = 300, bg = "white")

fwrite(
  data.table(
    figure = paste0("17", LETTERS[1:5]),
    file = basename(figure_paths),
    purpose = c(
      "Inspect normalization size factors by source",
      "Inspect uncorrected PCA separation by source",
      "Inspect FACS-gate structure in uncorrected PCA",
      "Inspect residual plate-configuration structure after ERCC removal",
      "Compare biological and technical metadata associations across PCs"
    )
  ),
  file.path(results_dir, "17_figure_manifest.csv")
)

# ------------------------------------------------------------------------------
# Save the normalized object and final audit summary
# ------------------------------------------------------------------------------

sce_metadata <- S4Vectors::metadata(sce)
sce_metadata$step <- 17L
sce_metadata$normalization <- list(
  method = "scran deconvolution size factors followed by logNormCounts",
  normalization_clusters = uniqueN(normalization_clusters),
  logcounts_assay = TRUE
)
sce_metadata$feature_selection <- list(
  method = "scran modelGeneVar/getTopHVGs",
  n_hvgs = n_hvgs,
  minimum_cells_detected = minimum_cells_for_hvg,
  excluded_patterns = c("MT-", "RPS/RPL", "ERCC", "MALAT1", "HBA1/HBA2/HBB")
)
sce_metadata$dimensionality_reduction <- list(
  method = "BiocSingular IrlbaPCA",
  n_pcs = n_pcs,
  batch_corrected = FALSE
)
S4Vectors::metadata(sce) <- sce_metadata

saveRDS(sce, output_sce_file, compress = "gzip")

dataset_qc <- data.table(
  metric = c(
    "genes", "cells", "tumor_cells", "pbmc_cells", "subjects",
    "input_specimens", "retained_specimens", "missing_specimens",
    "normalization_clusters", "minimum_size_factor", "median_size_factor",
    "maximum_size_factor", "logcounts_assay_present", "selected_hvgs",
    "pca_components", "pca_batch_corrected", "sce_object_size_mb"
  ),
  value = as.character(c(
    nrow(sce), ncol(sce), sum(sce$source_group == "Tumor"),
    sum(sce$source_group == "PBMC"), uniqueN(sce$subject_id),
    nrow(retention_by_specimen), sum(retention_by_specimen$n_final_cells > 0L),
    nrow(missing_specimens), uniqueN(normalization_clusters),
    min(size_factors), median(size_factors), max(size_factors),
    "TRUE", length(hvg_genes), ncol(pca_result$x), "FALSE",
    as.numeric(object.size(sce)) / 1024^2
  ))
)

fwrite(dataset_qc, file.path(results_dir, "17_dataset_qc.csv"))
writeLines(capture.output(sessionInfo()), file.path(results_dir, "17_sessionInfo.txt"))

message("Step 17 complete.")
message("Normalized PCA object: ", output_sce_file)
message("Dimensions: ", nrow(sce), " genes x ", ncol(sce), " cells")
message("HVGs: ", length(hvg_genes), "; PCs: ", ncol(pca_result$x))
message("Missing specimen(s):")
print(missing_specimens)
message("Normalization summary:")
print(normalization_summary)
message("Dataset QC:")
print(dataset_qc)
