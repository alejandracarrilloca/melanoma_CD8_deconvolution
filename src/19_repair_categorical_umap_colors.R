#!/usr/bin/env Rscript

# Repair Step 19 categorical UMAP figures without recomputing PCA, graph,
# Louvain clusters, UMAP, or lineage scores.

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

project_root <- normalizePath(getwd(), mustWork = TRUE)
input_sce_file <- file.path(
  project_root, "data", "analysis_ready", "scRNA", "GSE123139",
  "19_GSE123139_broad_clusters_sce.rds"
)
figures_dir <- file.path(project_root, "results", "19", "figures")
results_dir <- dirname(figures_dir)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_sce_file)) {
  stop(
    "Missing completed Step 19 object: ", input_sce_file,
    ". Run the full corrected Step 19 script first."
  )
}

message("Loading completed Step 19 object")
sce <- readRDS(input_sce_file)
if (!"UMAP" %in% SingleCellExperiment::reducedDimNames(sce)) {
  stop("The Step 19 object does not contain UMAP coordinates.")
}

required_metadata <- c(
  "source_group", "tcr_detected", "direct_t_cell_evidence",
  "direct_cd8_evidence"
)
missing_metadata <- setdiff(
  required_metadata,
  colnames(SummarizedExperiment::colData(sce))
)
if (length(missing_metadata)) {
  stop("Step 19 metadata is missing: ", paste(missing_metadata, collapse = ", "))
}

umap <- as.matrix(SingleCellExperiment::reducedDim(sce, "UMAP"))
if (nrow(umap) != ncol(sce) || ncol(umap) != 2L) {
  stop("Unexpected Step 19 UMAP dimensions.")
}

maximum_plot_cells <- 40000L
set.seed(123139L)
plot_index <- if (ncol(sce) > maximum_plot_cells) {
  sort(sample.int(ncol(sce), maximum_plot_cells))
} else {
  seq_len(ncol(sce))
}

metadata <- as.data.table(as.data.frame(SummarizedExperiment::colData(sce)))
plot_data <- metadata[plot_index, .(
  source_group = as.character(source_group),
  tcr_detected = as.logical(tcr_detected),
  direct_t_cell_evidence = as.logical(direct_t_cell_evidence),
  direct_cd8_evidence = as.logical(direct_cd8_evidence)
)]
plot_data[, `:=`(
  UMAP1 = umap[plot_index, 1L],
  UMAP2 = umap[plot_index, 2L],
  tcr_status = factor(
    ifelse(tcr_detected, "TCR mapped", "No mapped TCR"),
    levels = c("No mapped TCR", "TCR mapped")
  ),
  direct_evidence = factor(
    fifelse(
      direct_cd8_evidence,
      "Direct CD8 evidence",
      fifelse(
        direct_t_cell_evidence,
        "Direct T evidence",
        "No direct T evidence"
      )
    ),
    levels = c(
      "No direct T evidence", "Direct T evidence", "Direct CD8 evidence"
    )
  )
)]

theme_project <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.25)),
      plot.subtitle = element_text(color = "#5A5A5A"),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

save_png <- function(filename, plot) {
  ggsave(
    file.path(figures_dir, filename), plot = plot,
    width = 10, height = 8, units = "in", dpi = 180, bg = "white"
  )
}

p_source <- ggplot(plot_data, aes(UMAP1, UMAP2, color = source_group)) +
  geom_point(size = 0.30, alpha = 0.65) +
  scale_color_manual(values = c(PBMC = "#B8B8B8", Tumor = "#B48BE8")) +
  labs(
    title = "Sample source in broad-clustering UMAP space",
    subtitle = "PBMC and tumor are shown for interpretation, not used as graph inputs",
    color = "Sample source"
  ) +
  theme_project()
save_png("19B_umap_sample_source.png", p_source)

plot_tcr <- plot_data[order(tcr_detected)]
p_tcr <- ggplot(plot_tcr, aes(UMAP1, UMAP2, color = tcr_status)) +
  geom_point(size = 0.34, alpha = 0.65) +
  scale_color_manual(values = c(
    "No mapped TCR" = "#E4E4E4",
    "TCR mapped" = "#7651B9"
  )) +
  labs(
    title = "TCR-beta recovery in broad-clustering UMAP space",
    subtitle = "TCR is supporting evidence and is not required for T-cell assignment",
    color = NULL
  ) +
  theme_project()
save_png("19C_umap_tcr_recovery.png", p_tcr)

plot_evidence <- plot_data[order(as.integer(direct_evidence))]
p_evidence <- ggplot(
  plot_evidence,
  aes(UMAP1, UMAP2, color = direct_evidence)
) +
  geom_point(size = 0.34, alpha = 0.68) +
  scale_color_manual(values = c(
    "No direct T evidence" = "#E4E4E4",
    "Direct T evidence" = "#4C78A8",
    "Direct CD8 evidence" = "#7651B9"
  )) +
  labs(
    title = "Direct expression and TCR evidence for T-cell identity",
    subtitle = "Conservative evidence flags; cluster-aware CD8 assignment follows review",
    color = NULL
  ) +
  theme_project()
save_png("19F_umap_direct_t_cd8_evidence.png", p_evidence)

source_qc <- plot_data[, .(n_plotted = .N), by = .(category = source_group)]
source_qc[, figure := "19B_umap_sample_source.png"]
tcr_qc <- plot_data[, .(n_plotted = .N),
                    by = .(category = as.character(tcr_status))]
tcr_qc[, figure := "19C_umap_tcr_recovery.png"]
evidence_qc <- plot_data[, .(n_plotted = .N),
                         by = .(category = as.character(direct_evidence))]
evidence_qc[, figure := "19F_umap_direct_t_cd8_evidence.png"]
repair_qc <- rbindlist(
  list(source_qc, tcr_qc, evidence_qc),
  use.names = TRUE
)
setcolorder(repair_qc, c("figure", "category", "n_plotted"))
fwrite(repair_qc, file.path(results_dir, "19_categorical_figure_repair_qc.csv"))

message("Step 19 categorical figure repair complete.")
print(repair_qc)