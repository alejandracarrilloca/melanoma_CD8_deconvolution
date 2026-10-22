#!/usr/bin/env Rscript

# ==============================================================================
# Step 21: Define and recluster the tumor CD8 compartment in GSE123139
# ==============================================================================
#
# Step 20 separated the inclusive tumor T/NK compartment into nine subclusters.
# This step applies the reviewed, cluster-aware CD8 inclusion rule and records a
# strict direct-marker-only sensitivity set. It then performs fresh normalization,
# HVG selection, PCA, SNN clustering, and UMAP within the primary CD8 compartment.
# Final CD8-state names are deliberately deferred until the new clusters have
# been reviewed for marker programs, patient support, and clonotype behavior.
#
# Primary inclusion rule:
#   * retain all cells in TNC01, TNC04, TNC06, TNC07, TNC08, and TNC09;
#   * rescue only direct-CD8 cells from mixed TNC02 and NK-dominant TNC05;
#   * exclude CD4/Treg-dominant TNC03.
#
# TNC08 and TNC09 are retained so that real cells are not discarded, but they
# are flagged as subject-dominated and cannot define a standalone state in the
# next step without support from other subjects and clusters.
# ============================================================================== 

options(stringsAsFactors = FALSE)

required_packages <- c(
  "data.table", "Matrix", "SingleCellExperiment", "SummarizedExperiment",
  "S4Vectors", "scuttle", "scran", "BiocSingular", "BiocNeighbors",
  "igraph", "uwot", "ggplot2"
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
results_dir <- file.path(project_root, "results", "21")
figures_dir <- file.path(results_dir, "figures")
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

input_sce_file <- file.path(
  analysis_dir, "20_GSE123139_tumor_t_nk_refined_sce.rds"
)
output_sce_file <- file.path(
  analysis_dir, "21_GSE123139_tumor_cd8_state_discovery_sce.rds"
)

whole_subcluster_inclusion <- c(
  "TNC01", "TNC04", "TNC06", "TNC07", "TNC08", "TNC09"
)
direct_cd8_rescue_subclusters <- c("TNC02", "TNC05")
excluded_subclusters <- "TNC03"
subject_dominated_step20_subclusters <- c("TNC08", "TNC09")

n_hvgs <- 2000L
n_pcs <- 50L
n_pcs_graph <- 30L
knn_k <- 30L
minimum_cells_for_hvg <- 30L
maximum_plot_cells <- 30000L
random_seed <- 123139L

if (!file.exists(input_sce_file)) stop("Missing Step 20 SCE: ", input_sce_file)

marker_programs <- list(
  T_identity = c("CD3D", "CD3E", "CD3G", "TRAC", "CD247", "LCK", "LAT", "TRBC1", "TRBC2"),
  CD8_identity = c("CD8A", "CD8B"),
  CD4_memory = c("CD4", "IL7R", "CCR7", "LTB", "MAL", "LEF1", "SELL"),
  Treg = c("FOXP3", "IL2RA", "CTLA4", "TNFRSF18", "IKZF2"),
  NK = c("NKG7", "GNLY", "KLRD1", "FCGR3A", "TYROBP", "TRDC"),
  Cytotoxic = c("NKG7", "GNLY", "PRF1", "GZMB", "GZMH", "FGFBP2", "CTSW"),
  Dysfunction = c("PDCD1", "HAVCR2", "LAG3", "TIGIT", "TOX", "CXCL13", "ENTPD1"),
  Memory_transitional = c("TCF7", "IL7R", "CCR7", "LTB", "MAL", "LEF1", "SELL", "SLAMF6"),
  Progenitor_dysfunction = c("TCF7", "SLAMF6", "CXCR5", "TOX", "TIGIT"),
  Proliferation = c("MKI67", "TOP2A", "STMN1", "TUBA1B", "HMGB2"),
  Interferon_response = c("ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "STAT1"),
  Stress = c("FOS", "JUN", "JUNB", "DUSP1", "HSPA1A", "HSPA1B")
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

dominant_fraction <- function(x) {
  x <- x[!is.na(x) & nzchar(as.character(x))]
  if (!length(x)) return(NA_real_)
  as.numeric(max(table(x)) / length(x))
}

# ------------------------------------------------------------------------------
# Load Step 20 and apply the reviewed CD8 inclusion rule
# ------------------------------------------------------------------------------

message("Loading completed Step 20 object")
sce20 <- readRDS(input_sce_file)

required_metadata <- c(
  "source_group", "subject_id", "specimen_id", "facs_gate",
  "amplification_batch", "broad_cluster", "step20_subcluster",
  "step20_direct_T", "step20_direct_CD8", "tcr_detected",
  "tcr_productive_primary", "tcr_original_mc_grp", "tcr_clonotype_key",
  "endogenous_n_counts", "endogenous_n_genes"
)
missing_metadata <- setdiff(
  required_metadata, colnames(SummarizedExperiment::colData(sce20))
)
if (length(missing_metadata)) {
  stop("Step 20 metadata is missing: ", paste(missing_metadata, collapse = ", "))
}
if (!"counts" %in% SummarizedExperiment::assayNames(sce20)) {
  stop("Step 20 object does not contain raw counts.")
}
if (any(as.character(sce20$source_group) != "Tumor")) {
  stop("The Step 20 input unexpectedly contains non-tumor cells.")
}
if (!"UMAP_step20" %in% SingleCellExperiment::reducedDimNames(sce20)) {
  stop("The Step 20 input is missing UMAP_step20 coordinates.")
}
umap20_all <- SingleCellExperiment::reducedDim(sce20, "UMAP_step20")
if (nrow(umap20_all) != ncol(sce20) || ncol(umap20_all) != 2L) {
  stop("Unexpected Step 20 UMAP dimensions.")
}

step20_cluster <- as.character(sce20$step20_subcluster)
strict_direct_cd8 <- as.logical(sce20$step20_direct_CD8)
whole_context <- step20_cluster %in% whole_subcluster_inclusion
rescued_direct_cd8 <- step20_cluster %in% direct_cd8_rescue_subclusters &
  strict_direct_cd8
primary_cd8 <- whole_context | rescued_direct_cd8

if (any(step20_cluster %in% excluded_subclusters & primary_cd8)) {
  stop("A cell from the excluded CD4/Treg subcluster entered the CD8 set.")
}

inclusion_basis_all <- rep("Excluded after Step 20 review", ncol(sce20))
inclusion_basis_all[whole_context] <- "CD8-context subcluster"
inclusion_basis_all[rescued_direct_cd8] <- "Direct-CD8 rescue from mixed subcluster"

selection_audit <- data.table(
  cell_id = colnames(sce20),
  step20_subcluster = step20_cluster,
  strict_direct_cd8 = strict_direct_cd8,
  primary_cd8 = primary_cd8,
  inclusion_basis = inclusion_basis_all,
  subject_dominated_step20_origin =
    step20_cluster %in% subject_dominated_step20_subclusters
)

selection_by_subcluster <- selection_audit[, .(
  n_step20_cells = .N,
  n_strict_direct_cd8 = sum(strict_direct_cd8),
  percent_strict_direct_cd8 = 100 * mean(strict_direct_cd8),
  n_primary_cd8 = sum(primary_cd8),
  percent_primary_cd8 = 100 * mean(primary_cd8),
  n_subject_dominated_origin = sum(subject_dominated_step20_origin)
), by = step20_subcluster][order(step20_subcluster)]

fwrite(selection_by_subcluster,
       file.path(results_dir, "21_cd8_inclusion_by_step20_subcluster.csv"))
fwrite(selection_audit,
       file.path(results_dir, "21_cd8_selection_audit.tsv.gz"), sep = "\t")

sce <- sce20[, primary_cd8]
sce$step21_inclusion_basis <- inclusion_basis_all[primary_cd8]
sce$step21_strict_direct_cd8 <- strict_direct_cd8[primary_cd8]
sce$step21_subject_dominated_step20_origin <-
  step20_cluster[primary_cd8] %in% subject_dominated_step20_subclusters
rm(sce20)
invisible(gc())

if (!ncol(sce)) stop("No cells satisfy the primary CD8 inclusion rule.")
message("Primary tumor CD8 compartment: ", ncol(sce), " cells")

inclusion_basis_summary <- as.data.table(as.data.frame(
  SummarizedExperiment::colData(sce)
))[, .(
  n_cells = .N,
  n_subjects = uniqueN(subject_id),
  n_specimens = uniqueN(specimen_id),
  n_tcr = sum(tcr_detected),
  n_strict_direct_cd8 = sum(step21_strict_direct_cd8)
), by = step21_inclusion_basis][order(step21_inclusion_basis)]
fwrite(inclusion_basis_summary,
       file.path(results_dir, "21_cd8_inclusion_basis_summary.csv"))

# ------------------------------------------------------------------------------
# Re-normalize within the selected tumor CD8 compartment
# ------------------------------------------------------------------------------

message("Computing normalization clusters within the selected CD8 compartment")
set.seed(random_seed)
normalization_clusters <- scran::quickCluster(sce, min.size = 100L)
sce$step21_normalization_cluster <- as.character(normalization_clusters)

message("Estimating CD8-compartment scran size factors")
sce <- scran::computeSumFactors(sce, clusters = normalization_clusters)
size_factors <- as.numeric(SummarizedExperiment::colData(sce)$sizeFactor)
if (any(!is.finite(size_factors)) || any(size_factors <= 0)) {
  stop("Invalid CD8-compartment size factors were produced.")
}

message("Calculating CD8-compartment log-normalized expression")
sce <- scuttle::logNormCounts(sce)

normalization_summary <- data.table(
  n_cells = ncol(sce),
  normalization_clusters = uniqueN(sce$step21_normalization_cluster),
  minimum_size_factor = min(size_factors),
  median_size_factor = median(size_factors),
  maximum_size_factor = max(size_factors),
  median_endogenous_umis = median(as.numeric(sce$endogenous_n_counts)),
  median_endogenous_genes = median(as.numeric(sce$endogenous_n_genes))
)
fwrite(normalization_summary,
       file.path(results_dir, "21_normalization_summary.csv"))

# ------------------------------------------------------------------------------
# CD8-compartment HVGs and PCA
# ------------------------------------------------------------------------------

message("Modeling gene-level variance in the selected CD8 compartment")
variance_fit <- scran::modelGeneVar(sce)
counts <- SummarizedExperiment::assay(sce, "counts")
detected_cells <- as.integer(Matrix::rowSums(counts > 0))

technical_or_receptor_gene <-
  grepl("^MT-", rownames(sce)) |
  grepl("^RP[SL][0-9]", rownames(sce)) |
  grepl("^ERCC-", rownames(sce), ignore.case = TRUE) |
  rownames(sce) %in% c("MALAT1", "HBA1", "HBA2", "HBB") |
  grepl("^TR[ABDG][VJCD]", rownames(sce)) |
  grepl("^IG[HKL][VJCD]", rownames(sce))

eligible_hvg <- detected_cells >= minimum_cells_for_hvg &
  !technical_or_receptor_gene & is.finite(variance_fit$bio)
hvg_candidates <- rownames(sce)[eligible_hvg]
if (length(hvg_candidates) < n_hvgs) {
  stop("Only ", length(hvg_candidates), " genes are eligible for CD8 HVGs.")
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
fwrite(hvg_table, file.path(results_dir, "21_highly_variable_genes.csv"))

row_data <- SummarizedExperiment::rowData(sce)
row_data$step21_n_cells_detected <- detected_cells
row_data$step21_excluded_from_hvg <- technical_or_receptor_gene
row_data$step21_highly_variable <- rownames(sce) %chin% hvg_genes
row_data$step21_biological_variance <- variance_fit$bio
SummarizedExperiment::rowData(sce) <- row_data

message("Running CD8-compartment uncorrected PCA on ", n_hvgs, " HVGs")
logcounts_hvg <- SummarizedExperiment::assay(sce, "logcounts")[
  hvg_genes, , drop = FALSE
]
logcounts_hvg_sparse <- tryCatch(
  methods::as(logcounts_hvg, "dgCMatrix"),
  error = function(e) stop("Could not coerce CD8 HVGs to dgCMatrix: ",
                           conditionMessage(e))
)
pca_input <- Matrix::t(logcounts_hvg_sparse)
if (!inherits(pca_input, "Matrix") || nrow(pca_input) != ncol(sce) ||
    ncol(pca_input) != n_hvgs) {
  stop("Unexpected CD8 PCA input dimensions.")
}

set.seed(random_seed)
pca_result <- BiocSingular::runPCA(
  pca_input, rank = n_pcs, BSPARAM = BiocSingular::IrlbaParam()
)
rm(logcounts_hvg, logcounts_hvg_sparse, pca_input)
invisible(gc())

if (nrow(pca_result$x) != ncol(sce) || ncol(pca_result$x) != n_pcs) {
  stop("Unexpected CD8 PCA dimensions.")
}
rownames(pca_result$x) <- colnames(sce)
colnames(pca_result$x) <- paste0("PC", seq_len(n_pcs))
SingleCellExperiment::reducedDim(sce, "PCA_step21") <- pca_result$x

component_variance <- apply(pca_result$x, 2L, stats::var)
total_hvg_variance <- sum(
  variance_fit$total[match(hvg_genes, rownames(variance_fit))], na.rm = TRUE
)
pca_variance <- data.table(
  pc = paste0("PC", seq_len(n_pcs)),
  pc_number = seq_len(n_pcs),
  variance = component_variance,
  percent_hvg_variance = 100 * component_variance / total_hvg_variance,
  cumulative_percent_hvg_variance =
    100 * cumsum(component_variance) / total_hvg_variance
)
fwrite(pca_variance, file.path(results_dir, "21_pca_variance.csv"))

# ------------------------------------------------------------------------------
# CD8-state discovery graph and UMAP
# ------------------------------------------------------------------------------

pca_used <- pca_result$x[, seq_len(n_pcs_graph), drop = FALSE]
message("Building CD8-state discovery SNN graph")
set.seed(random_seed)
snn_graph <- scran::buildSNNGraph(
  pca_used, transposed = TRUE, k = knn_k,
  BNPARAM = BiocNeighbors::AnnoyParam()
)
if (igraph::vcount(snn_graph) != ncol(sce)) {
  stop("CD8 SNN graph vertex count does not match the cell count.")
}

message("Running Louvain clustering for CD8-state discovery")
set.seed(random_seed)
louvain <- igraph::cluster_louvain(
  snn_graph, weights = igraph::E(snn_graph)$weight
)
membership_raw <- as.integer(igraph::membership(louvain))
sizes_raw <- sort(table(membership_raw), decreasing = TRUE)
membership_order <- as.integer(names(sizes_raw))
cluster_map <- setNames(
  sprintf("CD8C%02d", seq_along(membership_order)), membership_order
)
sce$step21_cd8_cluster <- factor(
  unname(cluster_map[as.character(membership_raw)]),
  levels = unname(cluster_map)
)

message("Calculating deterministic CD8-compartment UMAP")
set.seed(random_seed)
umap <- uwot::umap(
  pca_used, n_neighbors = knn_k, min_dist = 0.30,
  metric = "euclidean", n_components = 2L, nn_method = "annoy",
  n_threads = 1L, n_sgd_threads = 1L, verbose = TRUE,
  ret_model = FALSE, seed = random_seed
)
colnames(umap) <- c("UMAP1", "UMAP2")
rownames(umap) <- colnames(sce)
SingleCellExperiment::reducedDim(sce, "UMAP_step21") <- umap

# ------------------------------------------------------------------------------
# State, lineage, and contamination program scores
# ------------------------------------------------------------------------------

marker_inventory <- rbindlist(lapply(names(marker_programs), function(program) {
  data.table(
    program = program,
    gene = marker_programs[[program]],
    present = marker_programs[[program]] %chin% rownames(sce)
  )
}))
fwrite(marker_inventory, file.path(results_dir, "21_marker_inventory.csv"))

present_markers <- intersect(
  unique(unlist(marker_programs, use.names = FALSE)), rownames(sce)
)
marker_logcounts <- as.matrix(
  SummarizedExperiment::assay(sce, "logcounts")[
    present_markers, , drop = FALSE
  ]
)
marker_means <- rowMeans(marker_logcounts)
marker_sds <- apply(marker_logcounts, 1L, stats::sd)
valid_scale <- is.finite(marker_sds) & marker_sds > 0
marker_z <- matrix(
  0, nrow = nrow(marker_logcounts), ncol = ncol(marker_logcounts),
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
colnames(score_matrix) <- paste0("step21_score_", names(marker_programs))
for (nm in colnames(score_matrix)) sce[[nm]] <- score_matrix[, nm]

# ------------------------------------------------------------------------------
# Cluster, subject, state, and clonotype audit outputs
# ------------------------------------------------------------------------------

cell_metadata <- as.data.table(as.data.frame(
  SummarizedExperiment::colData(sce)
))
cell_metadata[, cell_id := colnames(sce)]
setcolorder(cell_metadata, c("cell_id", setdiff(names(cell_metadata), "cell_id")))
if (anyDuplicated(names(cell_metadata))) {
  stop("Duplicated columns occur in the Step 21 cell metadata.")
}

cluster_summary <- cell_metadata[, .(
  n_cells = .N,
  percent_cd8_cells = 100 * .N / ncol(sce),
  n_subjects = uniqueN(subject_id),
  n_specimens = uniqueN(specimen_id),
  n_amplification_batches = uniqueN(amplification_batch),
  dominant_subject_fraction = dominant_fraction(subject_id),
  dominant_specimen_fraction = dominant_fraction(specimen_id),
  dominant_batch_fraction = dominant_fraction(amplification_batch),
  n_tcr = sum(tcr_detected),
  percent_tcr = 100 * mean(tcr_detected),
  n_productive_tcr = sum(tcr_productive_primary),
  n_strict_direct_cd8 = sum(step21_strict_direct_cd8),
  percent_strict_direct_cd8 = 100 * mean(step21_strict_direct_cd8),
  n_subject_dominated_step20_origin =
    sum(step21_subject_dominated_step20_origin),
  percent_subject_dominated_step20_origin =
    100 * mean(step21_subject_dominated_step20_origin),
  median_endogenous_umis = as.numeric(median(endogenous_n_counts)),
  median_endogenous_genes = as.numeric(median(endogenous_n_genes))
), by = step21_cd8_cluster][order(step21_cd8_cluster)]

program_columns <- paste0("step21_score_", names(marker_programs))
cluster_program_summary <- melt(
  cell_metadata[, c("step21_cd8_cluster", program_columns), with = FALSE][,
    lapply(.SD, mean, na.rm = TRUE),
    by = step21_cd8_cluster, .SDcols = program_columns
  ],
  id.vars = "step21_cd8_cluster", variable.name = "program",
  value.name = "mean_z_score"
)
cluster_program_summary[, program := sub("^step21_score_", "", program)]
setorder(cluster_program_summary, step21_cd8_cluster, program)

origin_summary <- cell_metadata[, .N,
  by = .(step21_cd8_cluster, step20_subcluster, step21_inclusion_basis)
]
origin_summary[, percent_within_cd8_cluster := 100 * N / sum(N),
               by = step21_cd8_cluster]
setorder(origin_summary, step21_cd8_cluster, -N)

subject_summary <- cell_metadata[, .N,
  by = .(step21_cd8_cluster, subject_id)
]
subject_summary[, percent_within_cd8_cluster := 100 * N / sum(N),
                by = step21_cd8_cluster]
setorder(subject_summary, step21_cd8_cluster, -N)

published_summary <- cell_metadata[
  tcr_detected == TRUE & !is.na(tcr_original_mc_grp) &
    nzchar(tcr_original_mc_grp),
  .N, by = .(step21_cd8_cluster, tcr_original_mc_grp)
]
published_summary[, percent_within_cd8_cluster_tcr := 100 * N / sum(N),
                  by = step21_cd8_cluster]
setorder(published_summary, step21_cd8_cluster, -N)

clonotype_summary <- cell_metadata[
  tcr_productive_primary == TRUE & !is.na(tcr_clonotype_key) &
    nzchar(tcr_clonotype_key),
  .(
    n_cells = .N,
    n_cd8_clusters = uniqueN(step21_cd8_cluster),
    cd8_clusters = paste(sort(unique(as.character(step21_cd8_cluster))),
                         collapse = ";"),
    n_step20_subclusters = uniqueN(step20_subcluster),
    step20_subclusters = paste(sort(unique(as.character(step20_subcluster))),
                               collapse = ";")
  ), by = .(tcr_clonotype_key, subject_id)
][order(-n_cells, subject_id, tcr_clonotype_key)]

fwrite(cluster_summary, file.path(results_dir, "21_cd8_cluster_summary.csv"))
fwrite(cluster_program_summary,
       file.path(results_dir, "21_cd8_cluster_program_scores.csv"))
fwrite(origin_summary,
       file.path(results_dir, "21_step20_origin_by_cd8_cluster.csv"))
fwrite(subject_summary,
       file.path(results_dir, "21_subject_by_cd8_cluster.csv"))
fwrite(published_summary,
       file.path(results_dir, "21_published_state_by_cd8_cluster.csv"))
fwrite(clonotype_summary,
       file.path(results_dir, "21_clonotype_by_cd8_cluster.tsv.gz"), sep = "\t")

# ------------------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------------------

selection_plot <- selection_audit
selection_plot[, `:=`(
  UMAP1 = umap20_all[, 1L],
  UMAP2 = umap20_all[, 2L],
  selection = factor(
    fifelse(primary_cd8, inclusion_basis, "Excluded after Step 20 review"),
    levels = c(
      "Excluded after Step 20 review", "CD8-context subcluster",
      "Direct-CD8 rescue from mixed subcluster"
    )
  )
)]

set.seed(random_seed)
selection_plot_index <- if (nrow(selection_plot) > maximum_plot_cells) {
  sort(sample.int(nrow(selection_plot), maximum_plot_cells))
} else seq_len(nrow(selection_plot))
selection_plot <- selection_plot[selection_plot_index]
selection_plot <- selection_plot[order(as.integer(selection))]

p_selection <- ggplot(selection_plot, aes(UMAP1, UMAP2, color = selection)) +
  geom_point(size = 0.36, alpha = 0.72) +
  scale_color_manual(values = c(
    "Excluded after Step 20 review" = "#E4E4E4",
    "CD8-context subcluster" = "#7651B9",
    "Direct-CD8 rescue from mixed subcluster" = "#2F7FBF"
  )) +
  labs(
    title = "Reviewed tumor CD8 selection in Step 20 UMAP space",
    subtitle = "Purple: CD8-context subclusters; blue: direct-CD8 rescue from mixed subclusters",
    color = NULL
  ) + theme_project()
save_png("21A_reviewed_cd8_selection.png", p_selection, 11, 8)

set.seed(random_seed)
plot_index <- if (ncol(sce) > maximum_plot_cells) {
  sort(sample.int(ncol(sce), maximum_plot_cells))
} else seq_len(ncol(sce))
plot_data <- cell_metadata[plot_index, .(
  step21_cd8_cluster, step20_subcluster, step21_inclusion_basis,
  step21_subject_dominated_step20_origin
)]
plot_data[, `:=`(
  UMAP1 = umap[plot_index, 1L], UMAP2 = umap[plot_index, 2L]
)]

cluster_levels <- levels(sce$step21_cd8_cluster)
cluster_colors <- grDevices::hcl.colors(length(cluster_levels), "Dynamic")
names(cluster_colors) <- cluster_levels
p_cluster <- ggplot(plot_data, aes(UMAP1, UMAP2, color = step21_cd8_cluster)) +
  geom_point(size = 0.40, alpha = 0.72) +
  scale_color_manual(values = cluster_colors, drop = FALSE) +
  labs(
    title = "Tumor CD8 state-discovery clusters",
    subtitle = "Fresh CD8-compartment normalization, HVGs, PCA, SNN graph, and UMAP",
    color = "CD8 cluster"
  ) + theme_project()
save_png("21B_cd8_state_discovery_clusters.png", p_cluster, 11, 8)

p_program <- ggplot(
  cluster_program_summary,
  aes(x = program, y = step21_cd8_cluster, fill = mean_z_score)
) +
  geom_tile(color = "white", linewidth = 0.30) +
  scale_fill_gradient2(
    low = "#3B6FB6", mid = "white", high = "#7651B9",
    midpoint = 0, name = "Mean gene-wise\nz-score"
  ) +
  labs(
    title = "Lineage and CD8-state programs by discovery cluster",
    subtitle = "Final biological state names will be assigned only after patient-support review",
    x = NULL, y = "CD8 discovery cluster"
  ) + theme_project() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_png("21C_cd8_cluster_program_heatmap.png", p_program, 13, 8)

p_subject <- ggplot(
  cluster_summary,
  aes(x = step21_cd8_cluster, y = 100 * dominant_subject_fraction)
) +
  geom_col(fill = "#B48BE8", color = "#7651B9", linewidth = 0.4) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "#555555") +
  labs(
    title = "Largest single-subject contribution to each CD8 discovery cluster",
    subtitle = "Clusters above 50% cannot define a standalone signature state",
    x = "CD8 discovery cluster", y = "Largest subject contribution (%)"
  ) + theme_project() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_png("21D_cd8_cluster_subject_dominance.png", p_subject, 10, 7)

if (nrow(published_summary)) {
  p_published <- ggplot(
    published_summary,
    aes(x = tcr_original_mc_grp, y = step21_cd8_cluster,
        fill = percent_within_cd8_cluster_tcr)
  ) +
    geom_tile(color = "white", linewidth = 0.30) +
    scale_fill_gradient(low = "white", high = "#7651B9",
                        name = "% of mapped TCR\ncells in cluster") +
    labs(
      title = "Published T-cell states across CD8 discovery clusters",
      subtitle = "Published mc_grp is validation metadata, not a training label",
      x = "Published mc_grp", y = "CD8 discovery cluster"
    ) + theme_project() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  save_png("21E_published_states_by_cd8_cluster.png", p_published, 12, 8)
}

origin_plot <- plot_data
origin_plot[, origin_warning := factor(
  fifelse(
    step21_subject_dominated_step20_origin,
    "From subject-dominated TNC08/TNC09",
    "From patient-supported Step 20 origins"
  ),
  levels = c(
    "From patient-supported Step 20 origins",
    "From subject-dominated TNC08/TNC09"
  )
)]
origin_plot <- origin_plot[order(as.integer(origin_warning))]
p_origin_warning <- ggplot(
  origin_plot, aes(UMAP1, UMAP2, color = origin_warning)
) +
  geom_point(size = 0.38, alpha = 0.72) +
  scale_color_manual(values = c(
    "From patient-supported Step 20 origins" = "#DADADA",
    "From subject-dominated TNC08/TNC09" = "#7651B9"
  )) +
  labs(
    title = "Subject-dominated Step 20 origins in CD8 discovery space",
    subtitle = "These cells are retained, but cannot independently define a state",
    color = NULL
  ) + theme_project()
save_png("21F_subject_dominated_origins.png", p_origin_warning, 10, 8)

# ------------------------------------------------------------------------------
# Validate and save
# ------------------------------------------------------------------------------

if (anyNA(sce$step21_cd8_cluster)) stop("Missing CD8 cluster assignments.")
if (sum(cluster_summary$n_cells) != ncol(sce)) {
  stop("CD8 cluster accounting failed.")
}
if (!identical(rownames(umap), colnames(sce))) {
  stop("CD8 UMAP rows do not align with SCE columns.")
}
if (any(as.character(sce$step20_subcluster) %in% excluded_subclusters)) {
  stop("Excluded TNC03 cells remain in the final Step 21 object.")
}

dataset_qc <- data.table(
  metric = c(
    "step20_focused_cells", "strict_direct_cd8_cells",
    "primary_cd8_cells", "primary_cd8_subjects", "primary_cd8_specimens",
    "whole_context_cells", "rescued_direct_cd8_cells",
    "subject_dominated_origin_cells", "normalization_clusters",
    "selected_hvgs", "pca_components", "graph_pcs", "snn_k",
    "cd8_discovery_clusters", "tcr_positive_cells",
    "productive_tcr_cells", "unique_productive_clonotypes",
    "cluster_accounting_validation", "final_state_labels_assigned",
    "output_sce_size_mb"
  ),
  value = as.character(c(
    nrow(selection_audit), sum(selection_audit$strict_direct_cd8),
    ncol(sce), uniqueN(sce$subject_id), uniqueN(sce$specimen_id),
    sum(whole_context), sum(rescued_direct_cd8),
    sum(sce$step21_subject_dominated_step20_origin),
    uniqueN(sce$step21_normalization_cluster), n_hvgs, n_pcs,
    n_pcs_graph, knn_k, uniqueN(sce$step21_cd8_cluster),
    sum(sce$tcr_detected), sum(sce$tcr_productive_primary),
    uniqueN(sce$tcr_clonotype_key[
      sce$tcr_productive_primary & !is.na(sce$tcr_clonotype_key)
    ]),
    TRUE, FALSE, NA_real_
  ))
)

metadata_list <- S4Vectors::metadata(sce)
metadata_list$step21 <- list(
  primary_inclusion_rule = paste(
    "whole", paste(whole_subcluster_inclusion, collapse = ","),
    "; direct-CD8 rescue", paste(direct_cd8_rescue_subclusters, collapse = ","),
    "; exclude", paste(excluded_subclusters, collapse = ",")
  ),
  strict_sensitivity_definition = "step20_direct_CD8",
  subject_dominated_step20_origins = subject_dominated_step20_subclusters,
  normalization = "scran deconvolution within selected tumor CD8 compartment",
  hvg_method = "scran modelGeneVar; technical and receptor genes excluded",
  n_hvgs = n_hvgs,
  pca = "BiocSingular IrlbaPCA; uncorrected",
  n_pcs = n_pcs,
  graph_pcs = n_pcs_graph,
  snn_k = knn_k,
  final_cd8_state_labels_assigned = FALSE
)
S4Vectors::metadata(sce) <- metadata_list

saveRDS(sce, output_sce_file, compress = FALSE)
dataset_qc[metric == "output_sce_size_mb", value := as.character(
  file.info(output_sce_file)$size / 1024^2
)]
fwrite(dataset_qc, file.path(results_dir, "21_dataset_qc.csv"))

figure_manifest <- data.table(
  figure = c(
    "21A_reviewed_cd8_selection.png",
    "21B_cd8_state_discovery_clusters.png",
    "21C_cd8_cluster_program_heatmap.png",
    "21D_cd8_cluster_subject_dominance.png",
    "21F_subject_dominated_origins.png"
  ),
  purpose = c(
    "Reviewed primary CD8 selection in Step 20 space",
    "Fresh CD8-compartment discovery clusters",
    "Lineage and CD8-state program evidence",
    "Patient representation by CD8 discovery cluster",
    "Locations of retained subject-dominated Step 20 origins"
  )
)
if (file.exists(file.path(figures_dir,
                          "21E_published_states_by_cd8_cluster.png"))) {
  figure_manifest <- rbind(
    figure_manifest,
    data.table(
      figure = "21E_published_states_by_cd8_cluster.png",
      purpose = "Published mc_grp validation labels"
    )
  )
}
setorder(figure_manifest, figure)
fwrite(figure_manifest, file.path(results_dir, "21_figure_manifest.csv"))
writeLines(capture.output(sessionInfo()), file.path(results_dir, "21_sessionInfo.txt"))

message("Step 21 tumor CD8 selection and state-discovery clustering complete.")
message("Strict direct-CD8 sensitivity cells: ", sum(selection_audit$strict_direct_cd8))
message("Primary CD8 cells: ", ncol(sce))
message("CD8 discovery clusters: ", uniqueN(sce$step21_cd8_cluster))
message("Output SCE: ", output_sce_file)
print(dataset_qc)
print(selection_by_subcluster)
print(cluster_summary)
