#!/usr/bin/env Rscript

# ==============================================================================
# Step 19: Broad clustering and evidence-based T-lineage/CD8 identification
# ==============================================================================
#
# This step uses the normalized, PCA-reduced, TCR-augmented Step 18 object.
# It builds an SNN graph from the first 30 uncorrected PCs, derives broad
# Louvain clusters, calculates a reproducible UMAP, and scores lineage marker
# programs. TCR recovery and the published mc_grp labels are used only as
# supporting/validation evidence; neither is required for classification.
#
# The direct_t_cell_evidence and direct_cd8_evidence flags are intentionally
# conservative diagnostics, not final biological labels. Cluster-level marker
# evidence and representation across subjects must be reviewed before the CD8
# compartment is finalized in Step 20.
# ==============================================================================

options(stringsAsFactors = FALSE)

required_packages <- c(
  "data.table", "Matrix", "SingleCellExperiment", "SummarizedExperiment",
  "S4Vectors", "scran", "igraph", "BiocNeighbors", "uwot", "ggplot2",
  "scales"
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
# Paths and configuration
# ------------------------------------------------------------------------------

project_root <- normalizePath(getwd(), mustWork = TRUE)
analysis_dir <- file.path(
  project_root, "data", "analysis_ready", "scRNA", "GSE123139"
)
results_dir <- file.path(project_root, "results", "19")
figures_dir <- file.path(results_dir, "figures")
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

input_sce_file <- file.path(
  analysis_dir, "18_GSE123139_normalized_pca_tcr_sce.rds"
)
output_sce_file <- file.path(
  analysis_dir, "19_GSE123139_broad_clusters_sce.rds"
)

n_pcs_graph <- 30L
knn_k <- 30L
umap_neighbors <- 30L
umap_min_dist <- 0.30
maximum_plot_cells <- 40000L
random_seed <- 123139L

if (!file.exists(input_sce_file)) stop("Missing Step 18 SCE: ", input_sce_file)

# Broad programs. These are used as continuous evidence, not as labels.
marker_programs <- list(
  T_lineage = c("CD3D", "CD3E", "TRAC", "CD247", "LCK", "LAT", "TRBC1", "TRBC2"),
  CD8 = c("CD8A", "CD8B"),
  CD4_Treg = c("CD4", "IL7R", "CCR7", "LTB", "MAL", "FOXP3", "IL2RA", "CTLA4", "IKZF2"),
  NK = c("NKG7", "GNLY", "KLRD1", "FCGR3A", "TRDC"),
  B_cell = c("CD79A", "CD79B", "MS4A1", "CD19", "CD22", "CD74", "HLA-DRA"),
  Myeloid = c("LST1", "FCER1G", "TYROBP", "AIF1", "CTSS", "LILRB1", "CTSD"),
  Melanoma = c("PMEL", "MLANA", "MITF", "TYR", "DCT", "SOX10")
)

palette_main <- c(
  purple = "#7651B9", lavender = "#B48BE8", gray = "#B8B8B8",
  dark = "#222222", blue = "#4C78A8", green = "#59A14F",
  orange = "#F28E2B", red = "#E15759"
)

theme_project <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.25)),
      plot.subtitle = element_text(color = "#5A5A5A"),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )
}

save_png <- function(filename, plot, width, height) {
  ggsave(
    file.path(figures_dir, filename), plot = plot, width = width,
    height = height, units = "in", dpi = 180, bg = "white"
  )
}

# ------------------------------------------------------------------------------
# Load and validate Step 18
# ------------------------------------------------------------------------------

message("Loading Step 18 SingleCellExperiment")
sce <- readRDS(input_sce_file)

if (!"logcounts" %in% SummarizedExperiment::assayNames(sce)) {
  stop("Step 18 object does not contain a logcounts assay.")
}
if (!"PCA" %in% SingleCellExperiment::reducedDimNames(sce)) {
  stop("Step 18 object does not contain PCA coordinates.")
}
if (anyDuplicated(rownames(sce)) || anyDuplicated(colnames(sce))) {
  stop("Duplicated gene or cell identifiers occur in the Step 18 object.")
}

required_metadata <- c(
  "source_group", "subject_id", "specimen_id", "facs_gate",
  "amplification_batch", "tcr_detected", "tcr_productive_primary",
  "tcr_original_mc_grp"
)
missing_metadata <- setdiff(
  required_metadata,
  colnames(SummarizedExperiment::colData(sce))
)
if (length(missing_metadata)) {
  stop("Step 18 metadata is missing: ", paste(missing_metadata, collapse = ", "))
}

pca <- SingleCellExperiment::reducedDim(sce, "PCA")
if (!is.matrix(pca)) pca <- as.matrix(pca)
if (nrow(pca) != ncol(sce) || ncol(pca) < n_pcs_graph) {
  stop("PCA dimensions do not support the requested graph construction.")
}
if (any(!is.finite(pca[, seq_len(n_pcs_graph), drop = FALSE]))) {
  stop("Non-finite values occur in the PCA coordinates used for clustering.")
}
pca_used <- pca[, seq_len(n_pcs_graph), drop = FALSE]

# ------------------------------------------------------------------------------
# SNN graph, Louvain clustering, and UMAP
# ------------------------------------------------------------------------------

message("Building a ", knn_k, "-nearest-neighbor SNN graph from ",
        n_pcs_graph, " PCs")
set.seed(random_seed)
snn_graph <- scran::buildSNNGraph(
  pca_used,
  transposed = TRUE,
  k = knn_k,
  BNPARAM = BiocNeighbors::AnnoyParam()
)
if (igraph::vcount(snn_graph) != ncol(sce)) {
  stop("SNN graph vertex count does not equal the number of cells.")
}

message("Running Louvain broad clustering")
set.seed(random_seed)
louvain <- igraph::cluster_louvain(snn_graph, weights = igraph::E(snn_graph)$weight)
cluster_number <- as.integer(igraph::membership(louvain))

# Stable labels are assigned by decreasing cluster size rather than the
# arbitrary membership integer returned by igraph.
cluster_sizes_raw <- sort(table(cluster_number), decreasing = TRUE)
cluster_order <- as.integer(names(cluster_sizes_raw))
cluster_map <- setNames(sprintf("C%02d", seq_along(cluster_order)), cluster_order)
sce$broad_cluster <- unname(cluster_map[as.character(cluster_number)])
sce$broad_cluster <- factor(sce$broad_cluster, levels = unname(cluster_map))

message("Calculating deterministic UMAP from the same ", n_pcs_graph, " PCs")
set.seed(random_seed)
umap <- uwot::umap(
  pca_used,
  n_neighbors = umap_neighbors,
  min_dist = umap_min_dist,
  metric = "euclidean",
  n_components = 2L,
  nn_method = "annoy",
  n_threads = 1L,
  n_sgd_threads = 1L,
  verbose = TRUE,
  ret_model = FALSE,
  seed = random_seed
)
colnames(umap) <- c("UMAP1", "UMAP2")
rownames(umap) <- colnames(sce)
SingleCellExperiment::reducedDim(sce, "UMAP") <- umap

# ------------------------------------------------------------------------------
# Continuous lineage-program scores and direct marker evidence
# ------------------------------------------------------------------------------

all_markers <- unique(unlist(marker_programs, use.names = FALSE))
marker_inventory <- rbindlist(lapply(names(marker_programs), function(program) {
  data.table(
    program = program,
    gene = marker_programs[[program]],
    present = marker_programs[[program]] %chin% rownames(sce)
  )
}))
fwrite(marker_inventory, file.path(results_dir, "19_marker_inventory.csv"))

present_markers <- intersect(all_markers, rownames(sce))
if (!length(present_markers)) stop("None of the lineage markers is present.")

message("Scoring broad lineage programs")
marker_logcounts <- as.matrix(
  SummarizedExperiment::assay(sce, "logcounts")[present_markers, , drop = FALSE]
)
marker_means <- rowMeans(marker_logcounts)
marker_sds <- apply(marker_logcounts, 1L, stats::sd)
valid_scale <- is.finite(marker_sds) & marker_sds > 0
marker_z <- marker_logcounts
marker_z[valid_scale, ] <- sweep(
  sweep(marker_logcounts[valid_scale, , drop = FALSE], 1L,
        marker_means[valid_scale], "-"),
  1L, marker_sds[valid_scale], "/"
)
marker_z[!valid_scale, ] <- 0

program_score <- function(genes) {
  genes <- intersect(genes, rownames(marker_z))
  if (!length(genes)) return(rep(NA_real_, ncol(sce)))
  colMeans(marker_z[genes, , drop = FALSE])
}

score_matrix <- vapply(marker_programs, program_score, numeric(ncol(sce)))
if (!identical(dim(score_matrix), c(ncol(sce), length(marker_programs)))) {
  stop("Unexpected lineage program score dimensions.")
}
colnames(score_matrix) <- paste0("score_", names(marker_programs))
for (nm in colnames(score_matrix)) sce[[nm]] <- score_matrix[, nm]

counts <- SummarizedExperiment::assay(sce, "counts")
detected_count <- function(genes) {
  genes <- intersect(genes, rownames(sce))
  if (!length(genes)) return(integer(ncol(sce)))
  as.integer(Matrix::colSums(counts[genes, , drop = FALSE] > 0))
}

sce$n_detected_T_lineage_markers <- detected_count(marker_programs$T_lineage)
sce$n_detected_CD8_markers <- detected_count(marker_programs$CD8)

# These flags describe direct evidence only. They deliberately do not attempt
# to rescue dropout-negative members of a T-cell-enriched cluster.
non_t_max <- pmax(
  sce$score_B_cell, sce$score_Myeloid, sce$score_Melanoma,
  na.rm = TRUE
)
sce$direct_t_cell_evidence <-
  as.logical(sce$tcr_detected) |
  (sce$n_detected_T_lineage_markers >= 1L & sce$score_T_lineage > non_t_max)
sce$direct_cd8_evidence <-
  as.character(sce$source_group) == "Tumor" &
  sce$direct_t_cell_evidence &
  sce$n_detected_CD8_markers >= 1L

# ------------------------------------------------------------------------------
# Cluster-level audit tables
# ------------------------------------------------------------------------------

cell_metadata <- as.data.table(as.data.frame(
  SummarizedExperiment::colData(sce)
))

# Earlier objects may already contain a cell_id metadata column. Assigning
# row names through keep.rownames in that situation creates duplicate column
# names, which recent ggplot2 versions reject even when the duplicated column
# is not mapped to an aesthetic.
cell_metadata[, cell_id := colnames(sce)]
setcolorder(cell_metadata, c("cell_id", setdiff(names(cell_metadata), "cell_id")))
if (anyDuplicated(names(cell_metadata))) {
  stop(
    "Duplicated columns remain in the plotting metadata: ",
    paste(unique(names(cell_metadata)[duplicated(names(cell_metadata))]),
          collapse = ", ")
  )
}

dominant_fraction <- function(x) {
  x <- x[!is.na(x) & nzchar(as.character(x))]
  if (!length(x)) return(NA_real_)
  max(table(x)) / length(x)
}

cluster_summary <- cell_metadata[, .(
  n_cells = .N,
  percent_all_cells = 100 * .N / ncol(sce),
  n_tumor = sum(source_group == "Tumor"),
  percent_tumor = 100 * mean(source_group == "Tumor"),
  n_pbmc = sum(source_group == "PBMC"),
  n_subjects = uniqueN(subject_id),
  n_specimens = uniqueN(specimen_id),
  n_amplification_batches = uniqueN(amplification_batch),
  dominant_subject_fraction = dominant_fraction(subject_id),
  dominant_batch_fraction = dominant_fraction(amplification_batch),
  n_tcr = sum(tcr_detected),
  percent_tcr = 100 * mean(tcr_detected),
  n_direct_t = sum(direct_t_cell_evidence),
  percent_direct_t = 100 * mean(direct_t_cell_evidence),
  n_direct_cd8 = sum(direct_cd8_evidence),
  percent_direct_cd8 = 100 * mean(direct_cd8_evidence),
  median_endogenous_umis = as.numeric(median(endogenous_n_counts)),
  median_endogenous_genes = as.numeric(median(endogenous_n_genes))
), by = broad_cluster][order(broad_cluster)]

program_columns <- paste0("score_", names(marker_programs))
cluster_program_summary <- melt(
  cell_metadata[, c("broad_cluster", program_columns), with = FALSE][,
    lapply(.SD, mean, na.rm = TRUE),
    by = broad_cluster,
    .SDcols = program_columns
  ],
  id.vars = "broad_cluster",
  variable.name = "program",
  value.name = "mean_z_score"
)
cluster_program_summary[, program := sub("^score_", "", program)]
setorder(cluster_program_summary, broad_cluster, program)

cluster_source <- cell_metadata[, .N, by = .(broad_cluster, source_group)]
cluster_source[, percent_within_cluster := 100 * N / sum(N), by = broad_cluster]
setorder(cluster_source, broad_cluster, source_group)

cluster_published_state <- cell_metadata[
  tcr_detected == TRUE & !is.na(tcr_original_mc_grp) & nzchar(tcr_original_mc_grp),
  .N,
  by = .(broad_cluster, tcr_original_mc_grp)
]
cluster_published_state[, percent_within_cluster_tcr := 100 * N / sum(N),
                        by = broad_cluster]
setorder(cluster_published_state, broad_cluster, -N)

evidence_by_cluster <- cell_metadata[, .(
  n_cells = .N,
  n_direct_t = sum(direct_t_cell_evidence),
  percent_direct_t = 100 * mean(direct_t_cell_evidence),
  n_direct_cd8 = sum(direct_cd8_evidence),
  percent_direct_cd8 = 100 * mean(direct_cd8_evidence),
  n_tcr = sum(tcr_detected),
  percent_tcr = 100 * mean(tcr_detected)
), by = .(broad_cluster, source_group)][order(broad_cluster, source_group)]

fwrite(cluster_summary, file.path(results_dir, "19_cluster_summary.csv"))
fwrite(cluster_program_summary,
       file.path(results_dir, "19_cluster_program_scores.csv"))
fwrite(cluster_source, file.path(results_dir, "19_cluster_by_source.csv"))
fwrite(cluster_published_state,
       file.path(results_dir, "19_published_state_by_cluster.csv"))
fwrite(evidence_by_cluster,
       file.path(results_dir, "19_direct_lineage_evidence_by_cluster.csv"))

# ------------------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------------------

set.seed(random_seed)
plot_index <- if (ncol(sce) > maximum_plot_cells) {
  sort(sample.int(ncol(sce), maximum_plot_cells))
} else {
  seq_len(ncol(sce))
}

plot_data <- cell_metadata[plot_index, .(
  cell_id,
  broad_cluster,
  source_group,
  tcr_detected,
  direct_t_cell_evidence,
  direct_cd8_evidence
)]
plot_data[, `:=`(
  UMAP1 = umap[plot_index, 1L],
  UMAP2 = umap[plot_index, 2L],
  tcr_status = factor(
    ifelse(tcr_detected, "TCR mapped", "No mapped TCR"),
    levels = c("No mapped TCR", "TCR mapped")
  ),
  direct_evidence = factor(
    fifelse(direct_cd8_evidence, "Direct CD8 evidence",
      fifelse(direct_t_cell_evidence, "Direct T evidence", "No direct T evidence")
    ),
    levels = c("No direct T evidence", "Direct T evidence", "Direct CD8 evidence")
  )
)]

cluster_levels <- levels(sce$broad_cluster)
cluster_colors <- grDevices::hcl.colors(length(cluster_levels), "Dynamic")
names(cluster_colors) <- cluster_levels

p_cluster <- ggplot(plot_data, aes(UMAP1, UMAP2, color = broad_cluster)) +
  geom_point(size = 0.30, alpha = 0.65) +
  scale_color_manual(values = cluster_colors, drop = FALSE) +
  labs(
    title = "Broad SNN clusters in GSE123139",
    subtitle = paste0(n_pcs_graph, " uncorrected PCs; Louvain clustering; ",
                      scales::comma(nrow(plot_data)), " cells plotted at most"),
    color = "Broad cluster"
  ) +
  theme_project()
save_png("19A_umap_broad_clusters.png", p_cluster, 11, 8)

p_source <- ggplot(plot_data, aes(UMAP1, UMAP2, color = source_group)) +
  geom_point(size = 0.30, alpha = 0.65) +
  scale_color_manual(values = c(PBMC = palette_main["gray"],
                                Tumor = palette_main["lavender"])) +
  labs(
    title = "Sample source in broad-clustering UMAP space",
    subtitle = "PBMC and tumor are shown for interpretation, not used as graph inputs",
    color = "Sample source"
  ) +
  theme_project()
save_png("19B_umap_sample_source.png", p_source, 10, 8)

plot_tcr <- plot_data[order(tcr_detected)]
p_tcr <- ggplot(plot_tcr, aes(UMAP1, UMAP2, color = tcr_status)) +
  geom_point(size = 0.34, alpha = 0.65) +
  scale_color_manual(values = c(
    "No mapped TCR" = "#E4E4E4", "TCR mapped" = palette_main["purple"]
  )) +
  labs(
    title = "TCR-beta recovery in broad-clustering UMAP space",
    subtitle = "TCR is supporting evidence and is not required for T-cell assignment",
    color = NULL
  ) +
  theme_project()
save_png("19C_umap_tcr_recovery.png", p_tcr, 10, 8)

p_program <- ggplot(
  cluster_program_summary,
  aes(x = program, y = broad_cluster, fill = mean_z_score)
) +
  geom_tile(color = "white", linewidth = 0.35) +
  scale_fill_gradient2(
    low = "#3B6FB6", mid = "white", high = palette_main["purple"],
    midpoint = 0, name = "Mean gene-wise\nz-score"
  ) +
  labs(
    title = "Broad lineage programs by SNN cluster",
    subtitle = "Continuous marker evidence; values are not cell-type labels",
    x = NULL, y = "Broad cluster"
  ) +
  theme_project() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_png("19D_cluster_lineage_program_heatmap.png", p_program, 10, 8)

if (nrow(cluster_published_state)) {
  p_published <- ggplot(
    cluster_published_state,
    aes(x = tcr_original_mc_grp, y = broad_cluster,
        fill = percent_within_cluster_tcr)
  ) +
    geom_tile(color = "white", linewidth = 0.30) +
    scale_fill_gradient(low = "white", high = palette_main["purple"],
                        name = "% of mapped TCR\ncells in cluster") +
    labs(
      title = "Published T-cell states across new broad clusters",
      subtitle = "Published mc_grp labels are validation metadata only",
      x = "Published mc_grp", y = "Broad cluster"
    ) +
    theme_project() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  save_png("19E_published_states_by_cluster.png", p_published, 11, 8)
}

plot_evidence <- plot_data[order(as.integer(direct_evidence))]
p_evidence <- ggplot(plot_evidence, aes(UMAP1, UMAP2, color = direct_evidence)) +
  geom_point(size = 0.34, alpha = 0.68) +
  scale_color_manual(values = c(
    "No direct T evidence" = "#E4E4E4",
    "Direct T evidence" = palette_main["blue"],
    "Direct CD8 evidence" = palette_main["purple"]
  )) +
  labs(
    title = "Direct expression and TCR evidence for T-cell identity",
    subtitle = "Conservative evidence flags; cluster-aware CD8 assignment follows review",
    color = NULL
  ) +
  theme_project()
save_png("19F_umap_direct_t_cd8_evidence.png", p_evidence, 10, 8)

# ------------------------------------------------------------------------------
# Final validations and save
# ------------------------------------------------------------------------------

if (!identical(rownames(umap), colnames(sce))) {
  stop("UMAP row names do not align with SCE columns.")
}
if (anyNA(sce$broad_cluster)) stop("Missing broad-cluster assignments occur.")
if (sum(cluster_summary$n_cells) != ncol(sce)) {
  stop("Cluster summary does not account for every cell.")
}
if (!identical(
  SingleCellExperiment::reducedDim(sce, "PCA"), pca
)) {
  stop("Original Step 18 PCA coordinates changed unexpectedly.")
}

dataset_qc <- data.table(
  metric = c(
    "genes", "cells", "tumor_cells", "pbmc_cells", "subjects",
    "broad_clusters", "graph_pcs", "snn_k", "umap_neighbors",
    "umap_min_dist", "tcr_positive_cells", "direct_t_evidence_cells",
    "direct_cd8_evidence_cells", "pca_coordinates_preserved",
    "umap_coordinates_present", "cluster_accounting_validation"
  ),
  value = as.character(c(
    nrow(sce), ncol(sce), sum(sce$source_group == "Tumor"),
    sum(sce$source_group == "PBMC"), uniqueN(sce$subject_id),
    length(unique(sce$broad_cluster)), n_pcs_graph, knn_k, umap_neighbors,
    umap_min_dist, sum(sce$tcr_detected), sum(sce$direct_t_cell_evidence),
    sum(sce$direct_cd8_evidence), TRUE, TRUE, TRUE
  ))
)

fwrite(dataset_qc, file.path(results_dir, "19_dataset_qc.csv"))
saveRDS(sce, output_sce_file, compress = FALSE)

figure_manifest <- data.table(
  figure = c(
    "19A_umap_broad_clusters.png",
    "19B_umap_sample_source.png",
    "19C_umap_tcr_recovery.png",
    "19D_cluster_lineage_program_heatmap.png",
    "19F_umap_direct_t_cd8_evidence.png"
  ),
  purpose = c(
    "Broad Louvain clusters in UMAP space",
    "Tumor and PBMC distribution in UMAP space",
    "Mapped TCR-beta records in UMAP space",
    "Cluster-level broad lineage program scores",
    "Direct T-lineage and CD8 evidence in UMAP space"
  )
)
if (file.exists(file.path(figures_dir, "19E_published_states_by_cluster.png"))) {
  figure_manifest <- rbind(
    figure_manifest,
    data.table(
      figure = "19E_published_states_by_cluster.png",
      purpose = "Published mc_grp validation labels by cluster"
    )
  )
}
setorder(figure_manifest, figure)
fwrite(figure_manifest, file.path(results_dir, "19_figure_manifest.csv"))
writeLines(capture.output(sessionInfo()), file.path(results_dir, "19_sessionInfo.txt"))

message("Step 19 broad clustering complete.")
message("Broad clusters: ", length(unique(sce$broad_cluster)))
message("Direct T-cell evidence: ", sum(sce$direct_t_cell_evidence), " / ", ncol(sce))
message("Direct tumor CD8 evidence: ", sum(sce$direct_cd8_evidence), " / ", ncol(sce))
message("Output SCE: ", output_sce_file)
print(dataset_qc)
print(cluster_summary)
