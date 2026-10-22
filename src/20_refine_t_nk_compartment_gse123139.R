#!/usr/bin/env Rscript

# ==============================================================================
# Step 20: Refine the tumor T/NK compartment in GSE123139
# ==============================================================================
#
# Step 19 identified seven broad clusters containing essentially all credible#!/usr/bin/env Rscript

# ==============================================================================
# Step 20: Refine the tumor T/NK compartment in GSE123139
# ==============================================================================
#
# Step 19 identified seven broad clusters containing essentially all credible
# T-lineage, CD8, and TCR evidence: C01, C02, C05, C06, C09, C11, and C14.
# These clusters are intentionally carried together because several contain
# mixtures of CD8 T cells, CD4/Treg cells, and NK-like cytotoxic cells.
#
# This step keeps tumor cells from those clusters, re-estimates normalization,
# selects HVGs while excluding receptor/technical genes, and performs a fresh
# PCA/SNN/UMAP analysis. It produces refined lineage/state evidence and
# patient-representation diagnostics. It does not yet assign final CD8 states;
# the refined subclusters must first be reviewed for CD4/NK contamination and
# subject-specific structure.
# ==============================================================================

options(stringsAsFactors = FALSE)

required_packages <- c(
  "data.table", "Matrix", "SingleCellExperiment", "SummarizedExperiment",
  "S4Vectors", "scuttle", "scran", "BiocSingular", "BiocNeighbors",
  "igraph", "uwot", "ggplot2", "scales"
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
results_dir <- file.path(project_root, "results", "20")
figures_dir <- file.path(results_dir, "figures")
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

input_sce_file <- file.path(
  analysis_dir, "19_GSE123139_broad_clusters_sce.rds"
)
output_sce_file <- file.path(
  analysis_dir, "20_GSE123139_tumor_t_nk_refined_sce.rds"
)

included_broad_clusters <- c("C01", "C02", "C05", "C06", "C09", "C11", "C14")
n_hvgs <- 2000L
n_pcs <- 50L
n_pcs_graph <- 30L
knn_k <- 30L
minimum_cells_for_hvg <- 50L
maximum_plot_cells <- 30000L
random_seed <- 123139L

if (!file.exists(input_sce_file)) stop("Missing Step 19 SCE: ", input_sce_file)

marker_programs <- list(
  T_identity = c("CD3D", "CD3E", "CD3G", "TRAC", "CD247", "LCK", "LAT", "TRBC1", "TRBC2"),
  CD8_identity = c("CD8A", "CD8B"),
  CD4_memory = c("CD4", "IL7R", "CCR7", "LTB", "MAL", "LEF1", "SELL"),
  Treg = c("FOXP3", "IL2RA", "CTLA4", "TNFRSF18", "IKZF2"),
  NK = c("NKG7", "GNLY", "KLRD1", "FCGR3A", "TYROBP", "TRDC"),
  Cytotoxic = c("NKG7", "GNLY", "PRF1", "GZMB", "GZMH", "FGFBP2", "CTSW"),
  Dysfunction = c("PDCD1", "HAVCR2", "LAG3", "TIGIT", "TOX", "CXCL13", "ENTPD1"),
  Memory_transitional = c("TCF7", "IL7R", "CCR7", "LTB", "MAL", "LEF1", "SELL"),
  Proliferation = c("MKI67", "TOP2A", "STMN1", "TUBA1B", "HMGB2"),
  B_cell = c("CD79A", "CD79B", "MS4A1", "CD19", "CD22", "CD74"),
  Myeloid = c("LST1", "FCER1G", "TYROBP", "AIF1", "CTSS", "LILRB1")
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

save_png <- function(filename, plot, width = 10, height = 8) {
  ggsave(
    file.path(figures_dir, filename), plot = plot, width = width,
    height = height, units = "in", dpi = 180, bg = "white"
  )
}

# ------------------------------------------------------------------------------
# Load Step 19 and extract the inclusive tumor T/NK compartment
# ------------------------------------------------------------------------------

message("Loading completed Step 19 object")
sce19 <- readRDS(input_sce_file)
input_step19_cell_count <- ncol(sce19)

required_metadata <- c(
  "source_group", "subject_id", "specimen_id", "facs_gate",
  "amplification_batch", "broad_cluster", "tcr_detected",
  "tcr_productive_primary", "tcr_original_mc_grp", "direct_t_cell_evidence",
  "direct_cd8_evidence", "endogenous_n_counts", "endogenous_n_genes"
)
missing_metadata <- setdiff(
  required_metadata,
  colnames(SummarizedExperiment::colData(sce19))
)
if (length(missing_metadata)) {
  stop("Step 19 metadata is missing: ", paste(missing_metadata, collapse = ", "))
}
if (!"counts" %in% SummarizedExperiment::assayNames(sce19)) {
  stop("Step 19 object does not contain raw counts.")
}
if (!all(included_broad_clusters %in% unique(as.character(sce19$broad_cluster)))) {
  stop("One or more expected broad T/NK clusters are absent from Step 19.")
}

keep <- as.character(sce19$source_group) == "Tumor" &
  as.character(sce19$broad_cluster) %in% included_broad_clusters
sce <- sce19[, keep]
rm(sce19)
invisible(gc())

if (!ncol(sce)) stop("No cells remain in the focused tumor T/NK compartment.")
if (any(as.character(sce$source_group) != "Tumor")) {
  stop("Non-tumor cells unexpectedly remain in the Step 20 subset.")
}
message("Focused tumor T/NK input: ", ncol(sce), " cells")

input_by_broad_cluster <- as.data.table(as.data.frame(
  SummarizedExperiment::colData(sce)
))[, .(
  n_cells = .N,
  n_subjects = uniqueN(subject_id),
  n_specimens = uniqueN(specimen_id),
  n_tcr = sum(tcr_detected),
  percent_tcr = 100 * mean(tcr_detected),
  n_direct_cd8 = sum(direct_cd8_evidence),
  percent_direct_cd8 = 100 * mean(direct_cd8_evidence)
), by = broad_cluster][order(broad_cluster)]
fwrite(input_by_broad_cluster,
       file.path(results_dir, "20_input_by_step19_cluster.csv"))

# ------------------------------------------------------------------------------
# Re-normalize within the focused compartment
# ------------------------------------------------------------------------------

message("Computing normalization clusters within the focused compartment")
set.seed(random_seed)
normalization_clusters <- scran::quickCluster(sce, min.size = 100L)
sce$step20_normalization_cluster <- as.character(normalization_clusters)

message("Estimating focused scran size factors")
sce <- scran::computeSumFactors(sce, clusters = normalization_clusters)
size_factors <- as.numeric(SummarizedExperiment::colData(sce)$sizeFactor)
if (any(!is.finite(size_factors)) || any(size_factors <= 0)) {
  stop("Invalid focused-compartment size factors were produced.")
}

message("Calculating focused log-normalized expression")
sce <- scuttle::logNormCounts(sce)

normalization_summary <- data.table(
  n_cells = ncol(sce),
  normalization_clusters = uniqueN(sce$step20_normalization_cluster),
  minimum_size_factor = min(size_factors),
  median_size_factor = median(size_factors),
  maximum_size_factor = max(size_factors),
  median_endogenous_umis = median(as.numeric(sce$endogenous_n_counts)),
  median_endogenous_genes = median(as.numeric(sce$endogenous_n_genes))
)
fwrite(normalization_summary,
       file.path(results_dir, "20_normalization_summary.csv"))

# ------------------------------------------------------------------------------
# Focused HVG selection and PCA
# ------------------------------------------------------------------------------

message("Modeling gene-level variance in the focused compartment")
variance_fit <- scran::modelGeneVar(sce)
counts <- SummarizedExperiment::assay(sce, "counts")
detected_cells <- as.integer(Matrix::rowSums(counts > 0))

# Receptor-variable genes can cause expanded clonotypes to form clusters even
# when their transcriptional states are otherwise similar. They remain in the
# object for inspection but are excluded from the embedding features.
technical_or_receptor_gene <-
  grepl("^MT-", rownames(sce)) |
  grepl("^RP[SL][0-9]", rownames(sce)) |
  grepl("^ERCC-", rownames(sce), ignore.case = TRUE) |
  rownames(sce) %in% c("MALAT1", "HBA1", "HBA2", "HBB") |
  grepl("^TR[ABDG][VJCD]", rownames(sce)) |
  grepl("^IG[HKL][VJCD]", rownames(sce))

eligible_hvg <-
  detected_cells >= minimum_cells_for_hvg &
  !technical_or_receptor_gene &
  is.finite(variance_fit$bio)

hvg_candidates <- rownames(sce)[eligible_hvg]
if (length(hvg_candidates) < n_hvgs) {
  stop("Only ", length(hvg_candidates), " genes are eligible for focused HVG selection.")
}
hvg_order <- order(variance_fit$bio[eligible_hvg], decreasing = TRUE)
hvg_genes <- hvg_candidates[hvg_order[seq_len(n_hvgs)]]

hvg_table <- data.table(
  gene = rownames(sce),
  n_cells_detected = detected_cells,
  biological_variance = variance_fit$bio,
  excluded_technical_or_receptor = technical_or_receptor_gene,
  eligible_for_hvg = eligible_hvg,
  selected_hvg = rownames(sce) %chin% hvg_genes
)[order(-biological_variance)]
fwrite(hvg_table, file.path(results_dir, "20_highly_variable_genes.csv"))

row_data <- SummarizedExperiment::rowData(sce)
row_data$step20_n_cells_detected <- detected_cells
row_data$step20_excluded_from_hvg <- technical_or_receptor_gene
row_data$step20_highly_variable <- rownames(sce) %chin% hvg_genes
row_data$step20_biological_variance <- variance_fit$bio
SummarizedExperiment::rowData(sce) <- row_data

message("Running focused uncorrected PCA on ", n_hvgs, " HVGs")
logcounts_hvg <- SummarizedExperiment::assay(sce, "logcounts")[
  hvg_genes, , drop = FALSE
]
logcounts_hvg_sparse <- tryCatch(
  methods::as(logcounts_hvg, "dgCMatrix"),
  error = function(e) stop(
    "Could not coerce focused HVG logcounts to dgCMatrix: ",
    conditionMessage(e)
  )
)
pca_input <- Matrix::t(logcounts_hvg_sparse)
if (!inherits(pca_input, "Matrix") ||
    nrow(pca_input) != ncol(sce) || ncol(pca_input) != n_hvgs) {
  stop("Unexpected focused PCA input dimensions.")
}

set.seed(random_seed)
pca_result <- BiocSingular::runPCA(
  pca_input,
  rank = n_pcs,
  BSPARAM = BiocSingular::IrlbaParam()
)
rm(logcounts_hvg, logcounts_hvg_sparse, pca_input)
invisible(gc())

if (nrow(pca_result$x) != ncol(sce) || ncol(pca_result$x) != n_pcs) {
  stop("Unexpected focused PCA dimensions.")
}
rownames(pca_result$x) <- colnames(sce)
colnames(pca_result$x) <- paste0("PC", seq_len(n_pcs))
SingleCellExperiment::reducedDim(sce, "PCA_step20") <- pca_result$x

component_variance <- apply(pca_result$x, 2L, stats::var)
total_hvg_variance <- sum(
  variance_fit$total[match(hvg_genes, rownames(variance_fit))],
  na.rm = TRUE
)
pca_variance <- data.table(
  pc = paste0("PC", seq_len(n_pcs)),
  pc_number = seq_len(n_pcs),
  variance = component_variance,
  percent_hvg_variance = 100 * component_variance / total_hvg_variance,
  cumulative_percent_hvg_variance =
    100 * cumsum(component_variance) / total_hvg_variance
)
fwrite(pca_variance, file.path(results_dir, "20_pca_variance.csv"))

# ------------------------------------------------------------------------------
# Focused graph, subclusters, and UMAP
# ------------------------------------------------------------------------------

pca_used <- pca_result$x[, seq_len(n_pcs_graph), drop = FALSE]
message("Building focused SNN graph")
set.seed(random_seed)
snn_graph <- scran::buildSNNGraph(
  pca_used,
  transposed = TRUE,
  k = knn_k,
  BNPARAM = BiocNeighbors::AnnoyParam()
)
if (igraph::vcount(snn_graph) != ncol(sce)) {
  stop("Focused SNN graph vertex count does not match the cell count.")
}

message("Running focused Louvain clustering")
set.seed(random_seed)
louvain <- igraph::cluster_louvain(
  snn_graph,
  weights = igraph::E(snn_graph)$weight
)
membership_raw <- as.integer(igraph::membership(louvain))
sizes_raw <- sort(table(membership_raw), decreasing = TRUE)
membership_order <- as.integer(names(sizes_raw))
subcluster_map <- setNames(
  sprintf("TNC%02d", seq_along(membership_order)),
  membership_order
)
sce$step20_subcluster <- factor(
  unname(subcluster_map[as.character(membership_raw)]),
  levels = unname(subcluster_map)
)

message("Calculating focused deterministic UMAP")
set.seed(random_seed)
umap <- uwot::umap(
  pca_used,
  n_neighbors = knn_k,
  min_dist = 0.30,
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
SingleCellExperiment::reducedDim(sce, "UMAP_step20") <- umap

# ------------------------------------------------------------------------------
# Refined lineage and CD8-state program scores
# ------------------------------------------------------------------------------

marker_inventory <- rbindlist(lapply(names(marker_programs), function(program) {
  data.table(
    program = program,
    gene = marker_programs[[program]],
    present = marker_programs[[program]] %chin% rownames(sce)
  )
}))
fwrite(marker_inventory, file.path(results_dir, "20_marker_inventory.csv"))

present_markers <- intersect(
  unique(unlist(marker_programs, use.names = FALSE)),
  rownames(sce)
)
marker_logcounts <- as.matrix(
  SummarizedExperiment::assay(sce, "logcounts")[present_markers, , drop = FALSE]
)
marker_means <- rowMeans(marker_logcounts)
marker_sds <- apply(marker_logcounts, 1L, stats::sd)
valid_scale <- is.finite(marker_sds) & marker_sds > 0
marker_z <- matrix(
  0,
  nrow = nrow(marker_logcounts),
  ncol = ncol(marker_logcounts),
  dimnames = dimnames(marker_logcounts)
)
marker_z[valid_scale, ] <- sweep(
  sweep(marker_logcounts[valid_scale, , drop = FALSE], 1L,
        marker_means[valid_scale], "-"),
  1L, marker_sds[valid_scale], "/"
)

program_score <- function(genes) {
  genes <- intersect(genes, rownames(marker_z))
  if (!length(genes)) return(rep(NA_real_, ncol(sce)))
  colMeans(marker_z[genes, , drop = FALSE])
}
score_matrix <- vapply(marker_programs, program_score, numeric(ncol(sce)))
colnames(score_matrix) <- paste0("step20_score_", names(marker_programs))
for (nm in colnames(score_matrix)) sce[[nm]] <- score_matrix[, nm]

detected_count <- function(genes) {
  genes <- intersect(genes, rownames(sce))
  if (!length(genes)) return(integer(ncol(sce)))
  as.integer(Matrix::colSums(counts[genes, , drop = FALSE] > 0))
}
sce$step20_n_T_identity_markers <- detected_count(marker_programs$T_identity)
sce$step20_n_CD8_identity_markers <- detected_count(marker_programs$CD8_identity)
sce$step20_n_Treg_markers <- detected_count(marker_programs$Treg)

# Direct evidence remains conservative. Final assignment will additionally use
# the refined subcluster profiles and patient representation.
sce$step20_direct_T <-
  as.logical(sce$tcr_detected) |
  (sce$step20_n_T_identity_markers >= 1L &
     sce$step20_score_T_identity > sce$step20_score_B_cell &
     sce$step20_score_T_identity > sce$step20_score_Myeloid)
sce$step20_direct_CD8 <-
  sce$step20_direct_T & sce$step20_n_CD8_identity_markers >= 1L

# ------------------------------------------------------------------------------
# Refined subcluster audit outputs
# ------------------------------------------------------------------------------

cell_metadata <- as.data.table(as.data.frame(
  SummarizedExperiment::colData(sce)
))
cell_metadata[, cell_id := colnames(sce)]
setcolorder(cell_metadata, c("cell_id", setdiff(names(cell_metadata), "cell_id")))
if (anyDuplicated(names(cell_metadata))) {
  stop("Duplicated columns occur in the Step 20 cell metadata.")
}

dominant_fraction <- function(x) {
  x <- x[!is.na(x) & nzchar(as.character(x))]
  if (!length(x)) return(NA_real_)
  as.numeric(max(table(x)) / length(x))
}

subcluster_summary <- cell_metadata[, .(
  n_cells = .N,
  percent_focused_cells = 100 * .N / ncol(sce),
  n_subjects = uniqueN(subject_id),
  n_specimens = uniqueN(specimen_id),
  n_amplification_batches = uniqueN(amplification_batch),
  dominant_subject_fraction = dominant_fraction(subject_id),
  dominant_specimen_fraction = dominant_fraction(specimen_id),
  dominant_batch_fraction = dominant_fraction(amplification_batch),
  n_tcr = sum(tcr_detected),
  percent_tcr = 100 * mean(tcr_detected),
  n_direct_T = sum(step20_direct_T),
  percent_direct_T = 100 * mean(step20_direct_T),
  n_direct_CD8 = sum(step20_direct_CD8),
  percent_direct_CD8 = 100 * mean(step20_direct_CD8),
  median_endogenous_umis = as.numeric(median(endogenous_n_counts)),
  median_endogenous_genes = as.numeric(median(endogenous_n_genes))
), by = step20_subcluster][order(step20_subcluster)]

program_columns <- paste0("step20_score_", names(marker_programs))
subcluster_program_summary <- melt(
  cell_metadata[, c("step20_subcluster", program_columns), with = FALSE][,
                                                                         lapply(.SD, mean, na.rm = TRUE),
                                                                         by = step20_subcluster,
                                                                         .SDcols = program_columns
  ],
  id.vars = "step20_subcluster",
  variable.name = "program",
  value.name = "mean_z_score"
)
subcluster_program_summary[, program := sub("^step20_score_", "", program)]
setorder(subcluster_program_summary, step20_subcluster, program)

subcluster_origin <- cell_metadata[, .N,
                                   by = .(step20_subcluster, broad_cluster)
]
subcluster_origin[, percent_within_subcluster := 100 * N / sum(N),
                  by = step20_subcluster]
setorder(subcluster_origin, step20_subcluster, -N)

subcluster_published_state <- cell_metadata[
  tcr_detected == TRUE & !is.na(tcr_original_mc_grp) & nzchar(tcr_original_mc_grp),
  .N,
  by = .(step20_subcluster, tcr_original_mc_grp)
]
subcluster_published_state[, percent_within_subcluster_tcr := 100 * N / sum(N),
                           by = step20_subcluster]
setorder(subcluster_published_state, step20_subcluster, -N)

subcluster_subject <- cell_metadata[, .N,
                                    by = .(step20_subcluster, subject_id)
]
subcluster_subject[, percent_within_subcluster := 100 * N / sum(N),
                   by = step20_subcluster]
setorder(subcluster_subject, step20_subcluster, -N)

subcluster_facs_gate <- cell_metadata[, .N,
                                      by = .(step20_subcluster, facs_gate)
]
subcluster_facs_gate[, percent_within_subcluster := 100 * N / sum(N),
                     by = step20_subcluster]
setorder(subcluster_facs_gate, step20_subcluster, -N)

fwrite(subcluster_summary, file.path(results_dir, "20_subcluster_summary.csv"))
fwrite(subcluster_program_summary,
       file.path(results_dir, "20_subcluster_program_scores.csv"))
fwrite(subcluster_origin,
       file.path(results_dir, "20_subcluster_step19_origin.csv"))
fwrite(subcluster_published_state,
       file.path(results_dir, "20_published_state_by_subcluster.csv"))
fwrite(subcluster_subject,
       file.path(results_dir, "20_subject_by_subcluster.csv"))
fwrite(subcluster_facs_gate,
       file.path(results_dir, "20_facs_gate_by_subcluster.csv"))

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
  step20_subcluster,
  broad_cluster,
  subject_id,
  tcr_detected,
  step20_direct_T,
  step20_direct_CD8
)]
plot_data[, `:=`(
  UMAP1 = umap[plot_index, 1L],
  UMAP2 = umap[plot_index, 2L],
  evidence = factor(
    fifelse(step20_direct_CD8, "Direct CD8 evidence",
            fifelse(step20_direct_T, "Direct T evidence", "No direct T evidence")
    ),
    levels = c("No direct T evidence", "Direct T evidence", "Direct CD8 evidence")
  )
)]

subcluster_levels <- levels(sce$step20_subcluster)
subcluster_colors <- grDevices::hcl.colors(length(subcluster_levels), "Dynamic")
names(subcluster_colors) <- subcluster_levels

p_subcluster <- ggplot(
  plot_data,
  aes(UMAP1, UMAP2, color = step20_subcluster)
) +
  geom_point(size = 0.35, alpha = 0.70) +
  scale_color_manual(values = subcluster_colors, drop = FALSE) +
  labs(
    title = "Focused tumor T/NK subclusters",
    subtitle = "Fresh normalization, HVG selection, PCA, SNN graph, and UMAP",
    color = "Subcluster"
  ) +
  theme_project()
save_png("20A_focused_umap_subclusters.png", p_subcluster, 11, 8)

broad_levels <- sort(unique(as.character(plot_data$broad_cluster)))
broad_colors <- grDevices::hcl.colors(length(broad_levels), "Set 2")
names(broad_colors) <- broad_levels
p_origin <- ggplot(plot_data, aes(UMAP1, UMAP2, color = broad_cluster)) +
  geom_point(size = 0.35, alpha = 0.70) +
  scale_color_manual(values = broad_colors) +
  labs(
    title = "Step 19 origins of the focused tumor T/NK compartment",
    subtitle = "Broad clusters were used inclusively, not as final identities",
    color = "Step 19 cluster"
  ) +
  theme_project()
save_png("20B_focused_umap_step19_origin.png", p_origin, 11, 8)

plot_evidence <- plot_data[order(as.integer(evidence))]
p_evidence <- ggplot(plot_evidence, aes(UMAP1, UMAP2, color = evidence)) +
  geom_point(size = 0.36, alpha = 0.72) +
  scale_color_manual(values = c(
    "No direct T evidence" = "#E4E4E4",
    "Direct T evidence" = "#4C78A8",
    "Direct CD8 evidence" = "#7651B9"
  )) +
  labs(
    title = "Direct T-lineage and CD8 evidence after focused reclustering",
    subtitle = "Final CD8 inclusion will use subcluster profiles and patient support",
    color = NULL
  ) +
  theme_project()
save_png("20C_focused_umap_lineage_evidence.png", p_evidence, 10, 8)

p_program <- ggplot(
  subcluster_program_summary,
  aes(x = program, y = step20_subcluster, fill = mean_z_score)
) +
  geom_tile(color = "white", linewidth = 0.30) +
  scale_fill_gradient2(
    low = "#3B6FB6", mid = "white", high = "#7651B9",
    midpoint = 0, name = "Mean gene-wise\nz-score"
  ) +
  labs(
    title = "Refined lineage and CD8-state programs by subcluster",
    subtitle = "Programs guide compartment selection; state labels are not yet assigned",
    x = NULL, y = "Focused subcluster"
  ) +
  theme_project() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_png("20D_subcluster_program_heatmap.png", p_program, 12, 8)

if (nrow(subcluster_published_state)) {
  p_published <- ggplot(
    subcluster_published_state,
    aes(x = tcr_original_mc_grp, y = step20_subcluster,
        fill = percent_within_subcluster_tcr)
  ) +
    geom_tile(color = "white", linewidth = 0.30) +
    scale_fill_gradient(low = "white", high = "#7651B9",
                        name = "% of mapped TCR\ncells in subcluster") +
    labs(
      title = "Published T-cell states across focused subclusters",
      subtitle = "Published mc_grp remains validation metadata only",
      x = "Published mc_grp", y = "Focused subcluster"
    ) +
    theme_project() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  save_png("20E_published_states_by_subcluster.png", p_published, 12, 8)
}

p_subject <- ggplot(
  subcluster_summary,
  aes(x = step20_subcluster, y = 100 * dominant_subject_fraction)
) +
  geom_col(fill = "#B48BE8", color = "#7651B9", linewidth = 0.4) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "#555555") +
  labs(
    title = "Largest single-subject contribution to each subcluster",
    subtitle = "Subclusters above 50% require special caution before state assignment",
    x = "Focused subcluster", y = "Largest subject contribution (%)"
  ) +
  theme_project() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_png("20F_subcluster_subject_dominance.png", p_subject, 10, 7)

# ------------------------------------------------------------------------------
# Validate and save
# ------------------------------------------------------------------------------

if (anyNA(sce$step20_subcluster)) stop("Missing focused subcluster assignments.")
if (sum(subcluster_summary$n_cells) != ncol(sce)) {
  stop("Focused subcluster accounting failed.")
}
if (!identical(rownames(umap), colnames(sce))) {
  stop("Focused UMAP rows do not align with SCE columns.")
}

dataset_qc <- data.table(
  metric = c(
    "input_step19_cells", "focused_tumor_cells", "focused_subjects",
    "included_step19_clusters", "normalization_clusters", "selected_hvgs",
    "pca_components", "graph_pcs", "snn_k", "focused_subclusters",
    "tcr_positive_cells", "direct_T_evidence_cells",
    "direct_CD8_evidence_cells", "receptor_genes_excluded_from_hvg",
    "subcluster_accounting_validation", "output_sce_size_mb"
  ),
  value = as.character(c(
    input_step19_cell_count, ncol(sce), uniqueN(sce$subject_id),
    paste(included_broad_clusters, collapse = ";"),
    uniqueN(sce$step20_normalization_cluster), n_hvgs, n_pcs,
    n_pcs_graph, knn_k, uniqueN(sce$step20_subcluster),
    sum(sce$tcr_detected), sum(sce$step20_direct_T),
    sum(sce$step20_direct_CD8),
    sum(grepl("^TR[ABDG][VJCD]|^IG[HKL][VJCD]", rownames(sce))),
    TRUE, NA_real_
  ))
)

metadata_list <- S4Vectors::metadata(sce)
metadata_list$step20 <- list(
  included_step19_clusters = included_broad_clusters,
  normalization = "scran deconvolution within focused tumor T/NK compartment",
  hvg_method = "scran modelGeneVar; technical and receptor genes excluded",
  n_hvgs = n_hvgs,
  pca = "BiocSingular IrlbaPCA; uncorrected",
  n_pcs = n_pcs,
  graph_pcs = n_pcs_graph,
  snn_k = knn_k,
  final_cd8_labels_assigned = FALSE
)
S4Vectors::metadata(sce) <- metadata_list

saveRDS(sce, output_sce_file, compress = FALSE)
dataset_qc[metric == "output_sce_size_mb", value := as.character(
  file.info(output_sce_file)$size / 1024^2
)]
fwrite(dataset_qc, file.path(results_dir, "20_dataset_qc.csv"))

figure_manifest <- data.table(
  figure = c(
    "20A_focused_umap_subclusters.png",
    "20B_focused_umap_step19_origin.png",
    "20C_focused_umap_lineage_evidence.png",
    "20D_subcluster_program_heatmap.png",
    "20F_subcluster_subject_dominance.png"
  ),
  purpose = c(
    "Focused Louvain subclusters",
    "Step 19 broad-cluster origins",
    "Direct T-lineage and CD8 evidence",
    "Refined lineage and state-program scores",
    "Largest subject contribution by subcluster"
  )
)
if (file.exists(file.path(figures_dir, "20E_published_states_by_subcluster.png"))) {
  figure_manifest <- rbind(
    figure_manifest,
    data.table(
      figure = "20E_published_states_by_subcluster.png",
      purpose = "Published mc_grp validation labels"
    )
  )
}
setorder(figure_manifest, figure)
fwrite(figure_manifest, file.path(results_dir, "20_figure_manifest.csv"))
writeLines(capture.output(sessionInfo()), file.path(results_dir, "20_sessionInfo.txt"))

message("Step 20 focused tumor T/NK refinement complete.")
message("Focused cells: ", ncol(sce))
message("Focused subclusters: ", uniqueN(sce$step20_subcluster))
message("Output SCE: ", output_sce_file)
print(dataset_qc)
print(subcluster_summary)

# T-lineage, CD8, and TCR evidence: C01, C02, C05, C06, C09, C11, and C14.
# These clusters are intentionally carried together because several contain
# mixtures of CD8 T cells, CD4/Treg cells, and NK-like cytotoxic cells.
#
# This step keeps tumor cells from those clusters, re-estimates normalization,
# selects HVGs while excluding receptor/technical genes, and performs a fresh
# PCA/SNN/UMAP analysis. It produces refined lineage/state evidence and
# patient-representation diagnostics. It does not yet assign final CD8 states;
# the refined subclusters must first be reviewed for CD4/NK contamination and
# subject-specific structure.
# ==============================================================================

options(stringsAsFactors = FALSE)

required_packages <- c(
  "data.table", "Matrix", "SingleCellExperiment", "SummarizedExperiment",
  "S4Vectors", "scuttle", "scran", "BiocSingular", "BiocNeighbors",
  "igraph", "uwot", "ggplot2", "scales"
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
results_dir <- file.path(project_root, "results", "20")
figures_dir <- file.path(results_dir, "figures")
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

input_sce_file <- file.path(
  analysis_dir, "19_GSE123139_broad_clusters_sce.rds"
)
output_sce_file <- file.path(
  analysis_dir, "20_GSE123139_tumor_t_nk_refined_sce.rds"
)

included_broad_clusters <- c("C01", "C02", "C05", "C06", "C09", "C11", "C14")
n_hvgs <- 2000L
n_pcs <- 50L
n_pcs_graph <- 30L
knn_k <- 30L
minimum_cells_for_hvg <- 50L
maximum_plot_cells <- 30000L
random_seed <- 123139L

if (!file.exists(input_sce_file)) stop("Missing Step 19 SCE: ", input_sce_file)

marker_programs <- list(
  T_identity = c("CD3D", "CD3E", "CD3G", "TRAC", "CD247", "LCK", "LAT", "TRBC1", "TRBC2"),
  CD8_identity = c("CD8A", "CD8B"),
  CD4_memory = c("CD4", "IL7R", "CCR7", "LTB", "MAL", "LEF1", "SELL"),
  Treg = c("FOXP3", "IL2RA", "CTLA4", "TNFRSF18", "IKZF2"),
  NK = c("NKG7", "GNLY", "KLRD1", "FCGR3A", "TYROBP", "TRDC"),
  Cytotoxic = c("NKG7", "GNLY", "PRF1", "GZMB", "GZMH", "FGFBP2", "CTSW"),
  Dysfunction = c("PDCD1", "HAVCR2", "LAG3", "TIGIT", "TOX", "CXCL13", "ENTPD1"),
  Memory_transitional = c("TCF7", "IL7R", "CCR7", "LTB", "MAL", "LEF1", "SELL"),
  Proliferation = c("MKI67", "TOP2A", "STMN1", "TUBA1B", "HMGB2"),
  B_cell = c("CD79A", "CD79B", "MS4A1", "CD19", "CD22", "CD74"),
  Myeloid = c("LST1", "FCER1G", "TYROBP", "AIF1", "CTSS", "LILRB1")
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

save_png <- function(filename, plot, width = 10, height = 8) {
  ggsave(
    file.path(figures_dir, filename), plot = plot, width = width,
    height = height, units = "in", dpi = 180, bg = "white"
  )
}

# ------------------------------------------------------------------------------
# Load Step 19 and extract the inclusive tumor T/NK compartment
# ------------------------------------------------------------------------------

message("Loading completed Step 19 object")
sce19 <- readRDS(input_sce_file)
input_step19_cell_count <- ncol(sce19)

required_metadata <- c(
  "source_group", "subject_id", "specimen_id", "facs_gate",
  "amplification_batch", "broad_cluster", "tcr_detected",
  "tcr_productive_primary", "tcr_original_mc_grp", "direct_t_cell_evidence",
  "direct_cd8_evidence", "endogenous_n_counts", "endogenous_n_genes"
)
missing_metadata <- setdiff(
  required_metadata,
  colnames(SummarizedExperiment::colData(sce19))
)
if (length(missing_metadata)) {
  stop("Step 19 metadata is missing: ", paste(missing_metadata, collapse = ", "))
}
if (!"counts" %in% SummarizedExperiment::assayNames(sce19)) {
  stop("Step 19 object does not contain raw counts.")
}
if (!all(included_broad_clusters %in% unique(as.character(sce19$broad_cluster)))) {
  stop("One or more expected broad T/NK clusters are absent from Step 19.")
}

keep <- as.character(sce19$source_group) == "Tumor" &
  as.character(sce19$broad_cluster) %in% included_broad_clusters
sce <- sce19[, keep]
rm(sce19)
invisible(gc())

if (!ncol(sce)) stop("No cells remain in the focused tumor T/NK compartment.")
if (any(as.character(sce$source_group) != "Tumor")) {
  stop("Non-tumor cells unexpectedly remain in the Step 20 subset.")
}
message("Focused tumor T/NK input: ", ncol(sce), " cells")

input_by_broad_cluster <- as.data.table(as.data.frame(
  SummarizedExperiment::colData(sce)
))[, .(
  n_cells = .N,
  n_subjects = uniqueN(subject_id),
  n_specimens = uniqueN(specimen_id),
  n_tcr = sum(tcr_detected),
  percent_tcr = 100 * mean(tcr_detected),
  n_direct_cd8 = sum(direct_cd8_evidence),
  percent_direct_cd8 = 100 * mean(direct_cd8_evidence)
), by = broad_cluster][order(broad_cluster)]
fwrite(input_by_broad_cluster,
       file.path(results_dir, "20_input_by_step19_cluster.csv"))

# ------------------------------------------------------------------------------
# Re-normalize within the focused compartment
# ------------------------------------------------------------------------------

message("Computing normalization clusters within the focused compartment")
set.seed(random_seed)
normalization_clusters <- scran::quickCluster(sce, min.size = 100L)
sce$step20_normalization_cluster <- as.character(normalization_clusters)

message("Estimating focused scran size factors")
sce <- scran::computeSumFactors(sce, clusters = normalization_clusters)
size_factors <- as.numeric(SummarizedExperiment::colData(sce)$sizeFactor)
if (any(!is.finite(size_factors)) || any(size_factors <= 0)) {
  stop("Invalid focused-compartment size factors were produced.")
}

message("Calculating focused log-normalized expression")
sce <- scuttle::logNormCounts(sce)

normalization_summary <- data.table(
  n_cells = ncol(sce),
  normalization_clusters = uniqueN(sce$step20_normalization_cluster),
  minimum_size_factor = min(size_factors),
  median_size_factor = median(size_factors),
  maximum_size_factor = max(size_factors),
  median_endogenous_umis = median(as.numeric(sce$endogenous_n_counts)),
  median_endogenous_genes = median(as.numeric(sce$endogenous_n_genes))
)
fwrite(normalization_summary,
       file.path(results_dir, "20_normalization_summary.csv"))

# ------------------------------------------------------------------------------
# Focused HVG selection and PCA
# ------------------------------------------------------------------------------

message("Modeling gene-level variance in the focused compartment")
variance_fit <- scran::modelGeneVar(sce)
counts <- SummarizedExperiment::assay(sce, "counts")
detected_cells <- as.integer(Matrix::rowSums(counts > 0))

# Receptor-variable genes can cause expanded clonotypes to form clusters even
# when their transcriptional states are otherwise similar. They remain in the
# object for inspection but are excluded from the embedding features.
technical_or_receptor_gene <-
  grepl("^MT-", rownames(sce)) |
  grepl("^RP[SL][0-9]", rownames(sce)) |
  grepl("^ERCC-", rownames(sce), ignore.case = TRUE) |
  rownames(sce) %in% c("MALAT1", "HBA1", "HBA2", "HBB") |
  grepl("^TR[ABDG][VJCD]", rownames(sce)) |
  grepl("^IG[HKL][VJCD]", rownames(sce))

eligible_hvg <-
  detected_cells >= minimum_cells_for_hvg &
  !technical_or_receptor_gene &
  is.finite(variance_fit$bio)

hvg_candidates <- rownames(sce)[eligible_hvg]
if (length(hvg_candidates) < n_hvgs) {
  stop("Only ", length(hvg_candidates), " genes are eligible for focused HVG selection.")
}
hvg_order <- order(variance_fit$bio[eligible_hvg], decreasing = TRUE)
hvg_genes <- hvg_candidates[hvg_order[seq_len(n_hvgs)]]

hvg_table <- data.table(
  gene = rownames(sce),
  n_cells_detected = detected_cells,
  biological_variance = variance_fit$bio,
  excluded_technical_or_receptor = technical_or_receptor_gene,
  eligible_for_hvg = eligible_hvg,
  selected_hvg = rownames(sce) %chin% hvg_genes
)[order(-biological_variance)]
fwrite(hvg_table, file.path(results_dir, "20_highly_variable_genes.csv"))

row_data <- SummarizedExperiment::rowData(sce)
row_data$step20_n_cells_detected <- detected_cells
row_data$step20_excluded_from_hvg <- technical_or_receptor_gene
row_data$step20_highly_variable <- rownames(sce) %chin% hvg_genes
row_data$step20_biological_variance <- variance_fit$bio
SummarizedExperiment::rowData(sce) <- row_data

message("Running focused uncorrected PCA on ", n_hvgs, " HVGs")
logcounts_hvg <- SummarizedExperiment::assay(sce, "logcounts")[
  hvg_genes, , drop = FALSE
]
logcounts_hvg_sparse <- tryCatch(
  methods::as(logcounts_hvg, "dgCMatrix"),
  error = function(e) stop(
    "Could not coerce focused HVG logcounts to dgCMatrix: ",
    conditionMessage(e)
  )
)
pca_input <- Matrix::t(logcounts_hvg_sparse)
if (!inherits(pca_input, "Matrix") ||
    nrow(pca_input) != ncol(sce) || ncol(pca_input) != n_hvgs) {
  stop("Unexpected focused PCA input dimensions.")
}

set.seed(random_seed)
pca_result <- BiocSingular::runPCA(
  pca_input,
  rank = n_pcs,
  BSPARAM = BiocSingular::IrlbaParam()
)
rm(logcounts_hvg, logcounts_hvg_sparse, pca_input)
invisible(gc())

if (nrow(pca_result$x) != ncol(sce) || ncol(pca_result$x) != n_pcs) {
  stop("Unexpected focused PCA dimensions.")
}
rownames(pca_result$x) <- colnames(sce)
colnames(pca_result$x) <- paste0("PC", seq_len(n_pcs))
SingleCellExperiment::reducedDim(sce, "PCA_step20") <- pca_result$x

component_variance <- apply(pca_result$x, 2L, stats::var)
total_hvg_variance <- sum(
  variance_fit$total[match(hvg_genes, rownames(variance_fit))],
  na.rm = TRUE
)
pca_variance <- data.table(
  pc = paste0("PC", seq_len(n_pcs)),
  pc_number = seq_len(n_pcs),
  variance = component_variance,
  percent_hvg_variance = 100 * component_variance / total_hvg_variance,
  cumulative_percent_hvg_variance =
    100 * cumsum(component_variance) / total_hvg_variance
)
fwrite(pca_variance, file.path(results_dir, "20_pca_variance.csv"))

# ------------------------------------------------------------------------------
# Focused graph, subclusters, and UMAP
# ------------------------------------------------------------------------------

pca_used <- pca_result$x[, seq_len(n_pcs_graph), drop = FALSE]
message("Building focused SNN graph")
set.seed(random_seed)
snn_graph <- scran::buildSNNGraph(
  pca_used,
  transposed = TRUE,
  k = knn_k,
  BNPARAM = BiocNeighbors::AnnoyParam()
)
if (igraph::vcount(snn_graph) != ncol(sce)) {
  stop("Focused SNN graph vertex count does not match the cell count.")
}

message("Running focused Louvain clustering")
set.seed(random_seed)
louvain <- igraph::cluster_louvain(
  snn_graph,
  weights = igraph::E(snn_graph)$weight
)
membership_raw <- as.integer(igraph::membership(louvain))
sizes_raw <- sort(table(membership_raw), decreasing = TRUE)
membership_order <- as.integer(names(sizes_raw))
subcluster_map <- setNames(
  sprintf("TNC%02d", seq_along(membership_order)),
  membership_order
)
sce$step20_subcluster <- factor(
  unname(subcluster_map[as.character(membership_raw)]),
  levels = unname(subcluster_map)
)

message("Calculating focused deterministic UMAP")
set.seed(random_seed)
umap <- uwot::umap(
  pca_used,
  n_neighbors = knn_k,
  min_dist = 0.30,
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
SingleCellExperiment::reducedDim(sce, "UMAP_step20") <- umap

# ------------------------------------------------------------------------------
# Refined lineage and CD8-state program scores
# ------------------------------------------------------------------------------

marker_inventory <- rbindlist(lapply(names(marker_programs), function(program) {
  data.table(
    program = program,
    gene = marker_programs[[program]],
    present = marker_programs[[program]] %chin% rownames(sce)
  )
}))
fwrite(marker_inventory, file.path(results_dir, "20_marker_inventory.csv"))

present_markers <- intersect(
  unique(unlist(marker_programs, use.names = FALSE)),
  rownames(sce)
)
marker_logcounts <- as.matrix(
  SummarizedExperiment::assay(sce, "logcounts")[present_markers, , drop = FALSE]
)
marker_means <- rowMeans(marker_logcounts)
marker_sds <- apply(marker_logcounts, 1L, stats::sd)
valid_scale <- is.finite(marker_sds) & marker_sds > 0
marker_z <- matrix(
  0,
  nrow = nrow(marker_logcounts),
  ncol = ncol(marker_logcounts),
  dimnames = dimnames(marker_logcounts)
)
marker_z[valid_scale, ] <- sweep(
  sweep(marker_logcounts[valid_scale, , drop = FALSE], 1L,
        marker_means[valid_scale], "-"),
  1L, marker_sds[valid_scale], "/"
)

program_score <- function(genes) {
  genes <- intersect(genes, rownames(marker_z))
  if (!length(genes)) return(rep(NA_real_, ncol(sce)))
  colMeans(marker_z[genes, , drop = FALSE])
}
score_matrix <- vapply(marker_programs, program_score, numeric(ncol(sce)))
colnames(score_matrix) <- paste0("step20_score_", names(marker_programs))
for (nm in colnames(score_matrix)) sce[[nm]] <- score_matrix[, nm]

detected_count <- function(genes) {
  genes <- intersect(genes, rownames(sce))
  if (!length(genes)) return(integer(ncol(sce)))
  as.integer(Matrix::colSums(counts[genes, , drop = FALSE] > 0))
}
sce$step20_n_T_identity_markers <- detected_count(marker_programs$T_identity)
sce$step20_n_CD8_identity_markers <- detected_count(marker_programs$CD8_identity)
sce$step20_n_Treg_markers <- detected_count(marker_programs$Treg)

# Direct evidence remains conservative. Final assignment will additionally use
# the refined subcluster profiles and patient representation.
sce$step20_direct_T <-
  as.logical(sce$tcr_detected) |
  (sce$step20_n_T_identity_markers >= 1L &
     sce$step20_score_T_identity > sce$step20_score_B_cell &
     sce$step20_score_T_identity > sce$step20_score_Myeloid)
sce$step20_direct_CD8 <-
  sce$step20_direct_T & sce$step20_n_CD8_identity_markers >= 1L

# ------------------------------------------------------------------------------
# Refined subcluster audit outputs
# ------------------------------------------------------------------------------

cell_metadata <- as.data.table(as.data.frame(
  SummarizedExperiment::colData(sce)
))
cell_metadata[, cell_id := colnames(sce)]
setcolorder(cell_metadata, c("cell_id", setdiff(names(cell_metadata), "cell_id")))
if (anyDuplicated(names(cell_metadata))) {
  stop("Duplicated columns occur in the Step 20 cell metadata.")
}

dominant_fraction <- function(x) {
  x <- x[!is.na(x) & nzchar(as.character(x))]
  if (!length(x)) return(NA_real_)
  as.numeric(max(table(x)) / length(x))
}

subcluster_summary <- cell_metadata[, .(
  n_cells = .N,
  percent_focused_cells = 100 * .N / ncol(sce),
  n_subjects = uniqueN(subject_id),
  n_specimens = uniqueN(specimen_id),
  n_amplification_batches = uniqueN(amplification_batch),
  dominant_subject_fraction = dominant_fraction(subject_id),
  dominant_specimen_fraction = dominant_fraction(specimen_id),
  dominant_batch_fraction = dominant_fraction(amplification_batch),
  n_tcr = sum(tcr_detected),
  percent_tcr = 100 * mean(tcr_detected),
  n_direct_T = sum(step20_direct_T),
  percent_direct_T = 100 * mean(step20_direct_T),
  n_direct_CD8 = sum(step20_direct_CD8),
  percent_direct_CD8 = 100 * mean(step20_direct_CD8),
  median_endogenous_umis = as.numeric(median(endogenous_n_counts)),
  median_endogenous_genes = as.numeric(median(endogenous_n_genes))
), by = step20_subcluster][order(step20_subcluster)]

program_columns <- paste0("step20_score_", names(marker_programs))
subcluster_program_summary <- melt(
  cell_metadata[, c("step20_subcluster", program_columns), with = FALSE][,
    lapply(.SD, mean, na.rm = TRUE),
    by = step20_subcluster,
    .SDcols = program_columns
  ],
  id.vars = "step20_subcluster",
  variable.name = "program",
  value.name = "mean_z_score"
)
subcluster_program_summary[, program := sub("^step20_score_", "", program)]
setorder(subcluster_program_summary, step20_subcluster, program)

subcluster_origin <- cell_metadata[, .N,
  by = .(step20_subcluster, broad_cluster)
]
subcluster_origin[, percent_within_subcluster := 100 * N / sum(N),
                  by = step20_subcluster]
setorder(subcluster_origin, step20_subcluster, -N)

subcluster_published_state <- cell_metadata[
  tcr_detected == TRUE & !is.na(tcr_original_mc_grp) & nzchar(tcr_original_mc_grp),
  .N,
  by = .(step20_subcluster, tcr_original_mc_grp)
]
subcluster_published_state[, percent_within_subcluster_tcr := 100 * N / sum(N),
                           by = step20_subcluster]
setorder(subcluster_published_state, step20_subcluster, -N)

subcluster_subject <- cell_metadata[, .N,
  by = .(step20_subcluster, subject_id)
]
subcluster_subject[, percent_within_subcluster := 100 * N / sum(N),
                   by = step20_subcluster]
setorder(subcluster_subject, step20_subcluster, -N)

subcluster_facs_gate <- cell_metadata[, .N,
  by = .(step20_subcluster, facs_gate)
]
subcluster_facs_gate[, percent_within_subcluster := 100 * N / sum(N),
                     by = step20_subcluster]
setorder(subcluster_facs_gate, step20_subcluster, -N)

fwrite(subcluster_summary, file.path(results_dir, "20_subcluster_summary.csv"))
fwrite(subcluster_program_summary,
       file.path(results_dir, "20_subcluster_program_scores.csv"))
fwrite(subcluster_origin,
       file.path(results_dir, "20_subcluster_step19_origin.csv"))
fwrite(subcluster_published_state,
       file.path(results_dir, "20_published_state_by_subcluster.csv"))
fwrite(subcluster_subject,
       file.path(results_dir, "20_subject_by_subcluster.csv"))
fwrite(subcluster_facs_gate,
       file.path(results_dir, "20_facs_gate_by_subcluster.csv"))

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
  step20_subcluster,
  broad_cluster,
  subject_id,
  tcr_detected,
  step20_direct_T,
  step20_direct_CD8
)]
plot_data[, `:=`(
  UMAP1 = umap[plot_index, 1L],
  UMAP2 = umap[plot_index, 2L],
  evidence = factor(
    fifelse(step20_direct_CD8, "Direct CD8 evidence",
      fifelse(step20_direct_T, "Direct T evidence", "No direct T evidence")
    ),
    levels = c("No direct T evidence", "Direct T evidence", "Direct CD8 evidence")
  )
)]

subcluster_levels <- levels(sce$step20_subcluster)
subcluster_colors <- grDevices::hcl.colors(length(subcluster_levels), "Dynamic")
names(subcluster_colors) <- subcluster_levels

p_subcluster <- ggplot(
  plot_data,
  aes(UMAP1, UMAP2, color = step20_subcluster)
) +
  geom_point(size = 0.35, alpha = 0.70) +
  scale_color_manual(values = subcluster_colors, drop = FALSE) +
  labs(
    title = "Focused tumor T/NK subclusters",
    subtitle = "Fresh normalization, HVG selection, PCA, SNN graph, and UMAP",
    color = "Subcluster"
  ) +
  theme_project()
save_png("20A_focused_umap_subclusters.png", p_subcluster, 11, 8)

broad_levels <- sort(unique(as.character(plot_data$broad_cluster)))
broad_colors <- grDevices::hcl.colors(length(broad_levels), "Set 2")
names(broad_colors) <- broad_levels
p_origin <- ggplot(plot_data, aes(UMAP1, UMAP2, color = broad_cluster)) +
  geom_point(size = 0.35, alpha = 0.70) +
  scale_color_manual(values = broad_colors) +
  labs(
    title = "Step 19 origins of the focused tumor T/NK compartment",
    subtitle = "Broad clusters were used inclusively, not as final identities",
    color = "Step 19 cluster"
  ) +
  theme_project()
save_png("20B_focused_umap_step19_origin.png", p_origin, 11, 8)

plot_evidence <- plot_data[order(as.integer(evidence))]
p_evidence <- ggplot(plot_evidence, aes(UMAP1, UMAP2, color = evidence)) +
  geom_point(size = 0.36, alpha = 0.72) +
  scale_color_manual(values = c(
    "No direct T evidence" = "#E4E4E4",
    "Direct T evidence" = "#4C78A8",
    "Direct CD8 evidence" = "#7651B9"
  )) +
  labs(
    title = "Direct T-lineage and CD8 evidence after focused reclustering",
    subtitle = "Final CD8 inclusion will use subcluster profiles and patient support",
    color = NULL
  ) +
  theme_project()
save_png("20C_focused_umap_lineage_evidence.png", p_evidence, 10, 8)

p_program <- ggplot(
  subcluster_program_summary,
  aes(x = program, y = step20_subcluster, fill = mean_z_score)
) +
  geom_tile(color = "white", linewidth = 0.30) +
  scale_fill_gradient2(
    low = "#3B6FB6", mid = "white", high = "#7651B9",
    midpoint = 0, name = "Mean gene-wise\nz-score"
  ) +
  labs(
    title = "Refined lineage and CD8-state programs by subcluster",
    subtitle = "Programs guide compartment selection; state labels are not yet assigned",
    x = NULL, y = "Focused subcluster"
  ) +
  theme_project() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_png("20D_subcluster_program_heatmap.png", p_program, 12, 8)

if (nrow(subcluster_published_state)) {
  p_published <- ggplot(
    subcluster_published_state,
    aes(x = tcr_original_mc_grp, y = step20_subcluster,
        fill = percent_within_subcluster_tcr)
  ) +
    geom_tile(color = "white", linewidth = 0.30) +
    scale_fill_gradient(low = "white", high = "#7651B9",
                        name = "% of mapped TCR\ncells in subcluster") +
    labs(
      title = "Published T-cell states across focused subclusters",
      subtitle = "Published mc_grp remains validation metadata only",
      x = "Published mc_grp", y = "Focused subcluster"
    ) +
    theme_project() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  save_png("20E_published_states_by_subcluster.png", p_published, 12, 8)
}

p_subject <- ggplot(
  subcluster_summary,
  aes(x = step20_subcluster, y = 100 * dominant_subject_fraction)
) +
  geom_col(fill = "#B48BE8", color = "#7651B9", linewidth = 0.4) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "#555555") +
  labs(
    title = "Largest single-subject contribution to each subcluster",
    subtitle = "Subclusters above 50% require special caution before state assignment",
    x = "Focused subcluster", y = "Largest subject contribution (%)"
  ) +
  theme_project() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_png("20F_subcluster_subject_dominance.png", p_subject, 10, 7)

# ------------------------------------------------------------------------------
# Validate and save
# ------------------------------------------------------------------------------

if (anyNA(sce$step20_subcluster)) stop("Missing focused subcluster assignments.")
if (sum(subcluster_summary$n_cells) != ncol(sce)) {
  stop("Focused subcluster accounting failed.")
}
if (!identical(rownames(umap), colnames(sce))) {
  stop("Focused UMAP rows do not align with SCE columns.")
}

dataset_qc <- data.table(
  metric = c(
    "input_step19_cells", "focused_tumor_cells", "focused_subjects",
    "included_step19_clusters", "normalization_clusters", "selected_hvgs",
    "pca_components", "graph_pcs", "snn_k", "focused_subclusters",
    "tcr_positive_cells", "direct_T_evidence_cells",
    "direct_CD8_evidence_cells", "receptor_genes_excluded_from_hvg",
    "subcluster_accounting_validation", "output_sce_size_mb"
  ),
  value = as.character(c(
    input_step19_cell_count, ncol(sce), uniqueN(sce$subject_id),
    paste(included_broad_clusters, collapse = ";"),
    uniqueN(sce$step20_normalization_cluster), n_hvgs, n_pcs,
    n_pcs_graph, knn_k, uniqueN(sce$step20_subcluster),
    sum(sce$tcr_detected), sum(sce$step20_direct_T),
    sum(sce$step20_direct_CD8),
    sum(grepl("^TR[ABDG][VJCD]|^IG[HKL][VJCD]", rownames(sce))),
    TRUE, NA_real_
  ))
)

metadata_list <- S4Vectors::metadata(sce)
metadata_list$step20 <- list(
  included_step19_clusters = included_broad_clusters,
  normalization = "scran deconvolution within focused tumor T/NK compartment",
  hvg_method = "scran modelGeneVar; technical and receptor genes excluded",
  n_hvgs = n_hvgs,
  pca = "BiocSingular IrlbaPCA; uncorrected",
  n_pcs = n_pcs,
  graph_pcs = n_pcs_graph,
  snn_k = knn_k,
  final_cd8_labels_assigned = FALSE
)
S4Vectors::metadata(sce) <- metadata_list

saveRDS(sce, output_sce_file, compress = FALSE)
dataset_qc[metric == "output_sce_size_mb", value := as.character(
  file.info(output_sce_file)$size / 1024^2
)]
fwrite(dataset_qc, file.path(results_dir, "20_dataset_qc.csv"))

figure_manifest <- data.table(
  figure = c(
    "20A_focused_umap_subclusters.png",
    "20B_focused_umap_step19_origin.png",
    "20C_focused_umap_lineage_evidence.png",
    "20D_subcluster_program_heatmap.png",
    "20F_subcluster_subject_dominance.png"
  ),
  purpose = c(
    "Focused Louvain subclusters",
    "Step 19 broad-cluster origins",
    "Direct T-lineage and CD8 evidence",
    "Refined lineage and state-program scores",
    "Largest subject contribution by subcluster"
  )
)
if (file.exists(file.path(figures_dir, "20E_published_states_by_subcluster.png"))) {
  figure_manifest <- rbind(
    figure_manifest,
    data.table(
      figure = "20E_published_states_by_subcluster.png",
      purpose = "Published mc_grp validation labels"
    )
  )
}
setorder(figure_manifest, figure)
fwrite(figure_manifest, file.path(results_dir, "20_figure_manifest.csv"))
writeLines(capture.output(sessionInfo()), file.path(results_dir, "20_sessionInfo.txt"))

message("Step 20 focused tumor T/NK refinement complete.")
message("Focused cells: ", ncol(sce))
message("Focused subclusters: ", uniqueN(sce$step20_subcluster))
message("Output SCE: ", output_sce_file)
print(dataset_qc)
print(subcluster_summary)
