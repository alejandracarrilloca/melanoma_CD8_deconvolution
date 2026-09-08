#!/usr/bin/env Rscript

# Repair the Step 17 PCA variance table and the PCA axis labels without
# repeating normalization, variance modeling, or PCA. This is needed for
# BiocSingular versions whose runPCA() result does not contain a `d` element.

options(stringsAsFactors = FALSE)

required_packages <- c(
  "data.table", "ggplot2", "SingleCellExperiment", "SummarizedExperiment"
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

results_dir <- file.path("results", "17")
figures_dir <- file.path(results_dir, "figures")
sce_file <- file.path(
  "data", "analysis_ready", "scRNA", "GSE123139",
  "17_GSE123139_normalized_pca_sce.rds"
)
hvg_file <- file.path(results_dir, "17_highly_variable_genes.csv")

if (!file.exists(sce_file)) stop("Missing Step 17 SCE: ", sce_file)
if (!file.exists(hvg_file)) stop("Missing Step 17 HVG table: ", hvg_file)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading saved Step 17 object")
sce <- readRDS(sce_file)
pca_scores <- SingleCellExperiment::reducedDim(sce, "PCA")

if (is.null(pca_scores) || nrow(pca_scores) != ncol(sce) || ncol(pca_scores) < 2L) {
  stop("The saved Step 17 object does not contain a valid PCA matrix.")
}

hvg_table <- fread(hvg_file)
selected_hvgs <- hvg_table[selected_hvg == TRUE]
if (nrow(selected_hvgs) != 2000L) {
  stop("Expected 2,000 selected HVGs, found ", nrow(selected_hvgs), ".")
}

total_hvg_variance <- sum(selected_hvgs$total_variance, na.rm = TRUE)
if (!is.finite(total_hvg_variance) || total_hvg_variance <= 0) {
  stop("Total variance across selected HVGs is not positive and finite.")
}

component_variance <- apply(pca_scores, 2L, stats::var)
if (any(!is.finite(component_variance)) || any(component_variance <= 0)) {
  stop("PCA score variances are missing, non-finite, or non-positive.")
}

n_pcs <- ncol(pca_scores)
singular_values <- sqrt(component_variance * (nrow(pca_scores) - 1L))
pca_variance <- data.table(
  pc = colnames(pca_scores),
  pc_number = seq_len(n_pcs),
  singular_value = singular_values,
  variance = component_variance,
  percent_hvg_variance = 100 * component_variance / total_hvg_variance
)
pca_variance[, cumulative_percent_hvg_variance :=
  cumsum(percent_hvg_variance)
]

if (any(diff(pca_variance$variance) > .Machine$double.eps^0.5)) {
  stop("PCA component variances are not in non-increasing order.")
}

fwrite(pca_variance, file.path(results_dir, "17_pca_variance.csv"))

metadata <- as.data.table(as.data.frame(SummarizedExperiment::colData(sce)))
plot_metadata <- cbind(metadata, as.data.table(pca_scores[, 1:2, drop = FALSE]))
set.seed(123139)
if (nrow(plot_metadata) > 40000L) {
  plot_metadata <- plot_metadata[sample(.N, 40000L)]
}

lavender <- "#B993F6"
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

x_label <- paste0("PC1 (", round(pca_variance$percent_hvg_variance[1L], 1), "%)")
y_label <- paste0("PC2 (", round(pca_variance$percent_hvg_variance[2L], 1), "%)")

base_pca <- function(fill_variable, title, subtitle, fill_label) {
  ggplot(
    plot_metadata,
    aes(x = PC1, y = PC2, fill = .data[[fill_variable]])
  ) +
    geom_point(
      shape = 21, size = 1.3, stroke = 0.2,
      color = gray_outline, alpha = 0.55
    ) +
    labs(
      title = title, subtitle = subtitle,
      x = x_label, y = y_label, fill = fill_label
    ) +
    theme_project
}

plot_17b <- base_pca(
  "source_group", "Uncorrected PCA by sample source",
  "No patient, plate, or batch correction has been applied", "Sample source"
) + scale_fill_manual(values = c("Tumor" = lavender, "PBMC" = gray_fill))

plot_17c <- base_pca(
  "facs_gate", "Uncorrected PCA by FACS gate",
  "FACS enrichment may reflect both technical design and real lineage composition",
  "FACS gate"
)

plot_17d <- base_pca(
  "plate_configuration", "Uncorrected PCA by plate configuration",
  "ERCC rows were removed before normalization", "Plate configuration"
)

ggsave(
  file.path(figures_dir, "17B_uncorrected_pca_by_source.png"),
  plot_17b, width = 9, height = 7, dpi = 300, bg = "white"
)
ggsave(
  file.path(figures_dir, "17C_uncorrected_pca_by_facs_gate.png"),
  plot_17c, width = 9, height = 7, dpi = 300, bg = "white"
)
ggsave(
  file.path(figures_dir, "17D_uncorrected_pca_by_plate_configuration.png"),
  plot_17d, width = 10, height = 7, dpi = 300, bg = "white"
)

repair_qc <- data.table(
  metric = c(
    "pca_cells", "pca_components", "selected_hvgs",
    "total_hvg_variance", "pc1_percent_hvg_variance",
    "pc2_percent_hvg_variance", "cumulative_percent_first_50_pcs",
    "variance_table_rows", "pca_figures_rewritten"
  ),
  value = as.character(c(
    nrow(pca_scores), n_pcs, nrow(selected_hvgs), total_hvg_variance,
    pca_variance$percent_hvg_variance[1L],
    pca_variance$percent_hvg_variance[2L],
    pca_variance$cumulative_percent_hvg_variance[n_pcs],
    nrow(pca_variance), 3L
  ))
)
fwrite(repair_qc, file.path(results_dir, "17_pca_variance_repair_qc.csv"))

message("Step 17 PCA variance repair complete.")
print(repair_qc)
