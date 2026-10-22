#!/usr/bin/env Rscript

# ==============================================================================
# Step 22: Patient-balanced validation of GSE123139 tumor CD8 states
# ==============================================================================
#
# Step 21 produced seven CD8 state-discovery clusters. This step tests whether
# their expression programs recur across patients, identifies cluster markers
# using within-patient contrasts, and measures clonotype sharing between
# clusters. It does not yet turn subject-dominated interferon/transitional
# clusters into standalone signature states.
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
# Paths and configuration
# ------------------------------------------------------------------------------

project_root <- normalizePath(getwd(), mustWork = TRUE)
analysis_dir <- file.path(
  project_root, "data", "analysis_ready", "scRNA", "GSE123139"
)
results_dir <- file.path(project_root, "results", "22")
figures_dir <- file.path(results_dir, "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

input_sce_file <- file.path(
  analysis_dir, "21_GSE123139_tumor_cd8_state_discovery_sce.rds"
)

minimum_cells_per_subject_cluster <- 10L
minimum_cells_in_subject_remainder <- 10L
top_markers_per_cluster <- 8L
maximum_plot_cells <- 30000L
random_seed <- 123139L

if (!file.exists(input_sce_file)) stop("Missing Step 21 SCE: ", input_sce_file)

state_programs <- list(
  CD8_identity = c("CD8A", "CD8B"),
  Cytotoxic = c("NKG7", "GNLY", "PRF1", "GZMB", "GZMH", "FGFBP2", "CTSW"),
  Dysfunction = c("PDCD1", "HAVCR2", "LAG3", "TIGIT", "TOX", "CXCL13", "ENTPD1"),
  Memory_transitional = c("TCF7", "IL7R", "CCR7", "LTB", "MAL", "LEF1", "SELL", "SLAMF6"),
  Progenitor_dysfunction = c("TCF7", "SLAMF6", "CXCR5", "TOX", "TIGIT"),
  Proliferation = c("MKI67", "TOP2A", "STMN1", "TUBA1B", "HMGB2"),
  Interferon_response = c("ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "STAT1"),
  Stress = c("FOS", "JUN", "JUNB", "DUSP1", "HSPA1A", "HSPA1B"),
  NK = c("NKG7", "GNLY", "KLRD1", "FCGR3A", "TYROBP", "TRDC"),
  CD4_memory = c("CD4", "IL7R", "CCR7", "LTB", "MAL", "LEF1", "SELL"),
  Treg = c("FOXP3", "IL2RA", "CTLA4", "TNFRSF18", "IKZF2")
)

review_status <- data.table(
  step21_cd8_cluster = sprintf("CD8C%02d", 1:7),
  current_interpretation = c(
    "Unresolved low-RNA/intermediate",
    "Dysfunctional CD8 anchor",
    "Unresolved intermediate",
    "Proliferating CD8 anchor",
    "Interferon-response overlay; subject dominated",
    "Cytotoxic CD8 anchor",
    "Transitional-looking; subject dominated"
  ),
  may_define_standalone_state = c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE),
  reason = c(
    "No dominant state program in Step 21",
    "CD8 identity and dysfunction programs; 21 subjects",
    "No dominant state program in Step 21",
    "Strong proliferation program; 24 subjects",
    "Strong interferon program but 78.5% from one subject",
    "Strong cytotoxic program; 21 subjects",
    "67.6% from one subject and mostly subject-dominated Step 20 origins"
  )
)
fwrite(review_status, file.path(results_dir, "22_prevalidation_decisions.csv"))

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

safe_wilcox <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3L || all(x == 0)) return(1)
  tryCatch(
    stats::wilcox.test(x, mu = 0, exact = FALSE)$p.value,
    error = function(e) 1
  )
}

# ------------------------------------------------------------------------------
# Load and validate Step 21
# ------------------------------------------------------------------------------

message("Loading completed Step 21 object")
sce <- readRDS(input_sce_file)

required_metadata <- c(
  "subject_id", "specimen_id", "amplification_batch", "step20_subcluster",
  "step21_cd8_cluster", "step21_strict_direct_cd8",
  "step21_subject_dominated_step20_origin", "tcr_detected",
  "tcr_productive_primary", "tcr_clonotype_key", "tcr_original_mc_grp",
  "endogenous_n_counts", "endogenous_n_genes"
)
missing_metadata <- setdiff(
  required_metadata, colnames(SummarizedExperiment::colData(sce))
)
if (length(missing_metadata)) {
  stop("Step 21 metadata is missing: ", paste(missing_metadata, collapse = ", "))
}
if (!all(c("counts", "logcounts") %in%
         SummarizedExperiment::assayNames(sce))) {
  stop("Step 21 must contain counts and logcounts assays.")
}
if (!"UMAP_step21" %in% SingleCellExperiment::reducedDimNames(sce)) {
  stop("Step 21 UMAP coordinates are missing.")
}

cluster <- as.character(sce$step21_cd8_cluster)
subject <- as.character(sce$subject_id)
cluster_levels <- levels(sce$step21_cd8_cluster)
if (is.null(cluster_levels)) cluster_levels <- sort(unique(cluster))

cell_metadata <- as.data.table(as.data.frame(
  SummarizedExperiment::colData(sce)
))
cell_metadata[, cell_id := colnames(sce)]
setcolorder(cell_metadata, c("cell_id", setdiff(names(cell_metadata), "cell_id")))
if (anyDuplicated(names(cell_metadata))) {
  stop("Duplicated columns occur in Step 21 cell metadata.")
}

# ------------------------------------------------------------------------------
# Patient-balanced program summaries
# ------------------------------------------------------------------------------

message("Calculating patient-balanced state-program summaries")
program_columns <- paste0("step21_score_", names(state_programs))
missing_scores <- setdiff(
  program_columns, colnames(SummarizedExperiment::colData(sce))
)
if (length(missing_scores)) {
  stop("Step 21 is missing program scores: ", paste(missing_scores, collapse = ", "))
}

subject_cluster_programs <- cell_metadata[, c(
  "subject_id", "step21_cd8_cluster", program_columns
), with = FALSE][,
  c(list(n_cells = .N), lapply(.SD, mean, na.rm = TRUE)),
  by = .(subject_id, step21_cd8_cluster),
  .SDcols = program_columns
]

subject_cluster_programs_long <- melt(
  subject_cluster_programs,
  id.vars = c("subject_id", "step21_cd8_cluster", "n_cells"),
  variable.name = "program", value.name = "subject_mean_z_score"
)
subject_cluster_programs_long[, program := sub("^step21_score_", "", program)]
subject_cluster_programs_long[, sufficient_cells :=
  n_cells >= minimum_cells_per_subject_cluster]

patient_balanced_programs <- subject_cluster_programs_long[
  sufficient_cells == TRUE,
  .(
    n_supported_subjects = uniqueN(subject_id),
    total_cells_in_supported_subjects = sum(n_cells),
    median_subject_mean_z_score = median(subject_mean_z_score, na.rm = TRUE),
    mean_subject_mean_z_score = mean(subject_mean_z_score, na.rm = TRUE),
    first_quartile_subject_mean_z_score =
      as.numeric(quantile(subject_mean_z_score, 0.25, na.rm = TRUE)),
    third_quartile_subject_mean_z_score =
      as.numeric(quantile(subject_mean_z_score, 0.75, na.rm = TRUE)),
    n_subjects_positive = sum(subject_mean_z_score > 0, na.rm = TRUE),
    percent_subjects_positive = 100 * mean(subject_mean_z_score > 0, na.rm = TRUE)
  ), by = .(step21_cd8_cluster, program)
][order(step21_cd8_cluster, program)]

fwrite(subject_cluster_programs_long,
       file.path(results_dir, "22_subject_cluster_program_scores.csv"))
fwrite(patient_balanced_programs,
       file.path(results_dir, "22_patient_balanced_program_summary.csv"))

# ------------------------------------------------------------------------------
# Within-patient marker contrasts
# ------------------------------------------------------------------------------

message("Preparing candidate genes for within-patient marker contrasts")
row_data <- SummarizedExperiment::rowData(sce)
if (!"step21_highly_variable" %in% colnames(row_data)) {
  stop("Step 21 HVG annotations are missing from rowData.")
}

program_genes <- unique(unlist(state_programs, use.names = FALSE))
candidate_genes <- unique(c(
  rownames(sce)[as.logical(row_data$step21_highly_variable)],
  intersect(program_genes, rownames(sce))
))
candidate_genes <- candidate_genes[candidate_genes %in% rownames(sce)]
if (!length(candidate_genes)) stop("No genes available for marker validation.")

logcounts_candidates <- SummarizedExperiment::assay(
  sce, "logcounts"
)[candidate_genes, , drop = FALSE]
logcounts_candidates <- tryCatch(
  methods::as(logcounts_candidates, "dgCMatrix"),
  error = function(e) stop("Could not coerce candidate logcounts: ",
                           conditionMessage(e))
)
counts_candidates <- SummarizedExperiment::assay(
  sce, "counts"
)[candidate_genes, , drop = FALSE]

marker_results <- vector("list", length(cluster_levels))
names(marker_results) <- cluster_levels

for (cluster_name in cluster_levels) {
  message("Patient-balanced marker contrasts for ", cluster_name)
  eligible_subjects <- unique(subject[cluster == cluster_name])
  eligible_subjects <- eligible_subjects[vapply(
    eligible_subjects,
    function(subject_name) {
      in_cluster <- subject == subject_name & cluster == cluster_name
      in_remainder <- subject == subject_name & cluster != cluster_name
      sum(in_cluster) >= minimum_cells_per_subject_cluster &&
        sum(in_remainder) >= minimum_cells_in_subject_remainder
    }, logical(1L)
  )]

  if (length(eligible_subjects) < 3L) {
    stop("Fewer than three patients support within-patient contrasts for ",
         cluster_name)
  }

  delta_matrix <- vapply(
    eligible_subjects,
    function(subject_name) {
      in_cluster <- subject == subject_name & cluster == cluster_name
      in_remainder <- subject == subject_name & cluster != cluster_name
      as.numeric(Matrix::rowMeans(logcounts_candidates[, in_cluster, drop = FALSE]) -
                   Matrix::rowMeans(logcounts_candidates[, in_remainder, drop = FALSE]))
    },
    numeric(length(candidate_genes))
  )
  if (is.null(dim(delta_matrix))) {
    delta_matrix <- matrix(delta_matrix, ncol = 1L)
  }
  rownames(delta_matrix) <- candidate_genes
  colnames(delta_matrix) <- eligible_subjects

  in_cluster_all <- cluster == cluster_name
  in_remainder_all <- cluster != cluster_name
  p_values <- apply(delta_matrix, 1L, safe_wilcox)

  marker_results[[cluster_name]] <- data.table(
    step21_cd8_cluster = cluster_name,
    gene = candidate_genes,
    n_paired_subjects = ncol(delta_matrix),
    mean_within_subject_delta = rowMeans(delta_matrix, na.rm = TRUE),
    median_within_subject_delta = apply(delta_matrix, 1L, median, na.rm = TRUE),
    percent_subjects_positive = 100 * rowMeans(delta_matrix > 0, na.rm = TRUE),
    percent_subjects_negative = 100 * rowMeans(delta_matrix < 0, na.rm = TRUE),
    p_value = p_values,
    adjusted_p_value = p.adjust(p_values, method = "BH"),
    percent_cells_detected_in_cluster = 100 * as.numeric(
      Matrix::rowMeans(counts_candidates[, in_cluster_all, drop = FALSE] > 0)
    ),
    percent_cells_detected_in_remainder = 100 * as.numeric(
      Matrix::rowMeans(counts_candidates[, in_remainder_all, drop = FALSE] > 0)
    )
  )[order(-median_within_subject_delta, -percent_subjects_positive,
          adjusted_p_value)]
}

marker_table <- rbindlist(marker_results, use.names = TRUE)
marker_table[, marker_rank := frank(
  -median_within_subject_delta, ties.method = "first"
), by = step21_cd8_cluster]

robust_markers <- marker_table[
  median_within_subject_delta > 0 &
    percent_subjects_positive >= 60 & adjusted_p_value <= 0.10
]
robust_markers[, robust_marker_rank := frank(
  -median_within_subject_delta, ties.method = "first"
), by = step21_cd8_cluster]

top_marker_list <- rbindlist(lapply(cluster_levels, function(cluster_name) {
  robust <- robust_markers[
    step21_cd8_cluster == cluster_name
  ][order(robust_marker_rank)]
  if (nrow(robust) >= top_markers_per_cluster) {
    return(robust[seq_len(top_markers_per_cluster)])
  }
  fallback <- marker_table[
    step21_cd8_cluster == cluster_name & median_within_subject_delta > 0
  ][order(marker_rank)]
  fallback[seq_len(min(top_markers_per_cluster, nrow(fallback)))]
}), use.names = TRUE, fill = TRUE)

fwrite(marker_table,
       file.path(results_dir, "22_patient_balanced_gene_markers.tsv.gz"),
       sep = "\t")
fwrite(robust_markers,
       file.path(results_dir, "22_robust_patient_balanced_markers.csv"))
fwrite(top_marker_list,
       file.path(results_dir, "22_top_patient_balanced_markers.csv"))

# ------------------------------------------------------------------------------
# Clonotype sharing across discovery clusters
# ------------------------------------------------------------------------------

message("Calculating productive-clonotype sharing across CD8 clusters")
clone_membership <- unique(cell_metadata[
  tcr_productive_primary == TRUE & !is.na(tcr_clonotype_key) &
    nzchar(tcr_clonotype_key),
  .(tcr_clonotype_key, subject_id, step21_cd8_cluster)
])

clone_counts <- clone_membership[, .(
  n_productive_clonotypes = uniqueN(tcr_clonotype_key)
), by = step21_cd8_cluster]

cluster_pairs <- CJ(
  cluster_1 = cluster_levels,
  cluster_2 = cluster_levels,
  unique = TRUE
)
cluster_pairs <- cluster_pairs[cluster_1 <= cluster_2]
clonotype_sharing <- cluster_pairs[, {
  cluster_1_key <- as.character(cluster_1[1L])
  cluster_2_key <- as.character(cluster_2[1L])
  clones_1 <- clone_membership[
    step21_cd8_cluster == cluster_1_key, unique(tcr_clonotype_key)
  ]
  clones_2 <- clone_membership[
    step21_cd8_cluster == cluster_2_key, unique(tcr_clonotype_key)
  ]
  shared <- length(intersect(clones_1, clones_2))
  union_n <- length(union(clones_1, clones_2))
  .(
    n_clonotypes_cluster_1 = length(clones_1),
    n_clonotypes_cluster_2 = length(clones_2),
    n_shared_clonotypes = shared,
    clonotype_jaccard = if (union_n > 0) shared / union_n else NA_real_,
    percent_cluster_1_clonotypes_shared =
      if (length(clones_1)) 100 * shared / length(clones_1) else NA_real_,
    percent_cluster_2_clonotypes_shared =
      if (length(clones_2)) 100 * shared / length(clones_2) else NA_real_
  )
}, by = .(cluster_1, cluster_2)]

fwrite(clone_counts, file.path(results_dir, "22_clonotype_counts_by_cluster.csv"))
fwrite(clonotype_sharing,
       file.path(results_dir, "22_clonotype_sharing_between_clusters.csv"))

# ------------------------------------------------------------------------------
# Sensitivity and subject-support summaries
# ------------------------------------------------------------------------------

cluster_support <- cell_metadata[, .(
  n_cells = .N,
  n_subjects = uniqueN(subject_id),
  n_subjects_with_at_least_10_cells = uniqueN(
    subject_id[ave(seq_len(.N), subject_id, FUN = length) >= 10L]
  ),
  dominant_subject_fraction = max(table(subject_id)) / .N,
  n_strict_direct_cd8 = sum(step21_strict_direct_cd8),
  percent_strict_direct_cd8 = 100 * mean(step21_strict_direct_cd8),
  n_tcr = sum(tcr_detected),
  percent_tcr = 100 * mean(tcr_detected),
  n_subject_dominated_origin = sum(step21_subject_dominated_step20_origin),
  percent_subject_dominated_origin =
    100 * mean(step21_subject_dominated_step20_origin)
), by = step21_cd8_cluster][order(step21_cd8_cluster)]
fwrite(cluster_support, file.path(results_dir, "22_cluster_support_summary.csv"))

# ------------------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------------------

p_program <- ggplot(
  patient_balanced_programs,
  aes(x = program, y = step21_cd8_cluster,
      fill = median_subject_mean_z_score)
) +
  geom_tile(color = "white", linewidth = 0.30) +
  scale_fill_gradient2(
    low = "#3B6FB6", mid = "white", high = "#7651B9",
    midpoint = 0, name = "Median patient\nmean z-score"
  ) +
  labs(
    title = "Patient-balanced programs across CD8 discovery clusters",
    subtitle = paste0(
      "Each patient contributes equally; patient-cluster groups require at least ",
      minimum_cells_per_subject_cluster, " cells"
    ),
    x = NULL, y = "CD8 discovery cluster"
  ) + theme_project() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_png("22A_patient_balanced_program_heatmap.png", p_program, 13, 8)

top_genes <- unique(top_marker_list$gene)
top_marker_plot <- marker_table[gene %chin% top_genes]
gene_order <- unique(top_marker_list[
  order(step21_cd8_cluster, marker_rank), gene
])
top_marker_plot[, gene := factor(gene, levels = rev(gene_order))]
p_markers <- ggplot(
  top_marker_plot,
  aes(x = step21_cd8_cluster, y = gene, fill = median_within_subject_delta)
) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(
    low = "#3B6FB6", mid = "white", high = "#7651B9",
    midpoint = 0, name = "Median within-patient\nexpression difference"
  ) +
  labs(
    title = "Patient-balanced marker genes for CD8 discovery clusters",
    subtitle = "Markers are ranked using cluster-versus-remainder contrasts within the same patients",
    x = "CD8 discovery cluster", y = NULL
  ) + theme_project()
save_png("22B_patient_balanced_marker_heatmap.png", p_markers, 11, 12)

subject_program_plot <- subject_cluster_programs_long[
  sufficient_cells == TRUE & program %chin% c(
    "Cytotoxic", "Dysfunction", "Memory_transitional", "Proliferation",
    "Interferon_response"
  )
]
p_subject_program <- ggplot(
  subject_program_plot,
  aes(x = step21_cd8_cluster, y = subject_mean_z_score,
      color = step21_cd8_cluster)
) +
  geom_hline(yintercept = 0, color = "#777777", linewidth = 0.35) +
  geom_boxplot(outlier.shape = NA, alpha = 0.15, show.legend = FALSE) +
  geom_jitter(width = 0.15, height = 0, size = 1.1, alpha = 0.75,
              show.legend = FALSE) +
  facet_wrap(~program, scales = "free_y") +
  labs(
    title = "State-program evidence across individual patients",
    subtitle = "Each point is one patient with at least 10 cells in that discovery cluster",
    x = "CD8 discovery cluster", y = "Patient mean gene-wise z-score"
  ) + theme_project() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_png("22C_state_programs_across_patients.png", p_subject_program, 13, 9)

sharing_symmetric <- rbind(
  clonotype_sharing[, .(
    cluster_1, cluster_2, n_shared_clonotypes, clonotype_jaccard
  )],
  clonotype_sharing[cluster_1 != cluster_2, .(
    cluster_1 = cluster_2,
    cluster_2 = cluster_1,
    n_shared_clonotypes,
    clonotype_jaccard
  )]
)
p_clone <- ggplot(
  sharing_symmetric,
  aes(x = cluster_1, y = cluster_2, fill = clonotype_jaccard)
) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = n_shared_clonotypes), size = 3) +
  scale_fill_gradient(low = "white", high = "#7651B9",
                      name = "Clonotype\nJaccard index") +
  labs(
    title = "Productive TCR clonotype sharing across CD8 clusters",
    subtitle = "Numbers show shared clonotypes; diagonal values show within-cluster totals",
    x = NULL, y = NULL
  ) + theme_project() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_png("22D_clonotype_sharing_heatmap.png", p_clone, 9, 8)

umap <- SingleCellExperiment::reducedDim(sce, "UMAP_step21")
set.seed(random_seed)
plot_index <- if (ncol(sce) > maximum_plot_cells) {
  sort(sample.int(ncol(sce), maximum_plot_cells))
} else seq_len(ncol(sce))
overlay_plot <- data.table(
  UMAP1 = umap[plot_index, 1L],
  UMAP2 = umap[plot_index, 2L],
  cluster = cluster[plot_index],
  warning_origin = as.logical(
    sce$step21_subject_dominated_step20_origin[plot_index]
  )
)
overlay_plot[, category := factor(
  fifelse(
    cluster == "CD8C05", "Interferon cluster (CD8C05)",
    fifelse(cluster == "CD8C07", "Transitional-looking cluster (CD8C07)",
            "Patient-supported discovery clusters")
  ),
  levels = c(
    "Patient-supported discovery clusters",
    "Interferon cluster (CD8C05)",
    "Transitional-looking cluster (CD8C07)"
  )
)]
overlay_plot <- overlay_plot[order(as.integer(category))]
p_overlay <- ggplot(overlay_plot, aes(UMAP1, UMAP2, color = category)) +
  geom_point(size = 0.40, alpha = 0.72) +
  scale_color_manual(values = c(
    "Patient-supported discovery clusters" = "#DCDCDC",
    "Interferon cluster (CD8C05)" = "#3B6FB6",
    "Transitional-looking cluster (CD8C07)" = "#7651B9"
  )) +
  labs(
    title = "Subject-dominated overlays requiring state merging",
    subtitle = "CD8C05 and CD8C07 are tested as overlays, not standalone signature states",
    color = NULL
  ) + theme_project()
save_png("22E_subject_dominated_overlays.png", p_overlay, 10, 8)

# ------------------------------------------------------------------------------
# Validate and report
# ------------------------------------------------------------------------------

if (sum(cluster_support$n_cells) != ncol(sce)) {
  stop("Step 22 cluster accounting failed.")
}
if (!all(cluster_levels %in% patient_balanced_programs$step21_cd8_cluster)) {
  stop("One or more clusters lack patient-balanced program summaries.")
}
if (!all(cluster_levels %in% marker_table$step21_cd8_cluster)) {
  stop("One or more clusters lack within-patient marker results.")
}

dataset_qc <- data.table(
  metric = c(
    "cells", "subjects", "cd8_discovery_clusters", "candidate_marker_genes",
    "minimum_cells_per_subject_cluster", "patient_cluster_program_groups",
    "within_patient_marker_tests", "robust_markers_fdr10_support60",
    "productive_clonotypes", "cluster_pairs_tested",
    "standalone_state_anchors_before_validation",
    "subject_dominated_overlay_clusters", "cluster_accounting_validation",
    "final_state_labels_assigned"
  ),
  value = as.character(c(
    ncol(sce), uniqueN(subject), length(cluster_levels), length(candidate_genes),
    minimum_cells_per_subject_cluster,
    uniqueN(subject_cluster_programs_long[
      sufficient_cells == TRUE,
      paste(subject_id, step21_cd8_cluster, sep = "::")
    ]),
    nrow(marker_table), nrow(robust_markers),
    uniqueN(clone_membership$tcr_clonotype_key), nrow(clonotype_sharing),
    sum(review_status$may_define_standalone_state),
    sum(!review_status$may_define_standalone_state &
          grepl("subject dominated", review_status$current_interpretation)),
    TRUE, FALSE
  ))
)
fwrite(dataset_qc, file.path(results_dir, "22_dataset_qc.csv"))

figure_manifest <- data.table(
  figure = c(
    "22A_patient_balanced_program_heatmap.png",
    "22B_patient_balanced_marker_heatmap.png",
    "22C_state_programs_across_patients.png",
    "22D_clonotype_sharing_heatmap.png",
    "22E_subject_dominated_overlays.png"
  ),
  purpose = c(
    "Equal-patient-weight state program comparison",
    "Within-patient discovery-cluster markers",
    "Patient-to-patient recurrence of major state programs",
    "TCR clonotype sharing between discovery clusters",
    "Visualization of subject-dominated overlays"
  )
)
fwrite(figure_manifest, file.path(results_dir, "22_figure_manifest.csv"))
writeLines(capture.output(sessionInfo()), file.path(results_dir, "22_sessionInfo.txt"))

message("Step 22 patient-balanced CD8-state validation complete.")
message("Candidate genes tested: ", length(candidate_genes))
message("Robust patient-balanced markers: ", nrow(robust_markers))
message("Final state labels remain intentionally unassigned.")
print(dataset_qc)
print(patient_balanced_programs)
print(top_marker_list)
