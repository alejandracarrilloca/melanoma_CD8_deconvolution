#!/usr/bin/env Rscript

# ==============================================================================
# Step 23: Freeze final GSE123139 tumor-CD8 states
# ==============================================================================
#
# Step 22 established three reproducible state anchors (dysfunctional,
# proliferating, and cytotoxic), a broad memory/intermediate compartment, and
# two subject-dominated overlays. This step freezes four deconvolution-ready
# CD8 labels. CD8C05 is retained as an interferon overlay but receives a base
# state using productive-clonotype evidence first and cross-patient kNN second.
# CD8C07 is merged into memory/transitional and is never used as a standalone
# signature state.
# ==============================================================================

options(stringsAsFactors = FALSE)

required_packages <- c(
  "data.table", "Matrix", "SingleCellExperiment", "SummarizedExperiment",
  "BiocSingular", "BiocNeighbors", "ggplot2"
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
# Paths and fixed decisions
# ------------------------------------------------------------------------------

project_root <- normalizePath(getwd(), mustWork = TRUE)
analysis_dir <- file.path(
  project_root, "data", "analysis_ready", "scRNA", "GSE123139"
)
results_dir <- file.path(project_root, "results", "23")
figures_dir <- file.path(results_dir, "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

input_sce_file <- file.path(
  analysis_dir, "21_GSE123139_tumor_cd8_state_discovery_sce.rds"
)
output_sce_file <- file.path(
  analysis_dir, "23_GSE123139_tumor_cd8_final_states_sce.rds"
)

random_seed <- 123139L
n_transfer_hvgs <- 1500L
n_transfer_pcs <- 30L
k_vote <- 30L
k_candidates <- 500L
minimum_clone_anchor_cells <- 2L
minimum_clone_consensus <- 0.60
minimum_patient_cells_for_support <- 10L
maximum_plot_cells <- 30000L

final_state_levels <- c(
  "CD8_memory_transitional",
  "CD8_dysfunctional",
  "CD8_cytotoxic",
  "CD8_proliferating"
)

cluster_decisions <- data.table(
  step21_cd8_cluster = sprintf("CD8C%02d", 1:7),
  step23_action = c(
    "Merge with memory/transitional",
    "Dysfunctional anchor",
    "Merge with memory/transitional",
    "Proliferating state",
    "Interferon overlay; infer base state",
    "Cytotoxic anchor",
    "Merge with memory/transitional"
  ),
  fixed_state = c(
    "CD8_memory_transitional", "CD8_dysfunctional",
    "CD8_memory_transitional", "CD8_proliferating", NA,
    "CD8_cytotoxic", "CD8_memory_transitional"
  ),
  used_as_expression_anchor = c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE),
  rationale = c(
    "Modest recurrent memory signal; technical markers excluded downstream",
    "Recurrent CD8, dysfunction, and progenitor-dysfunction programs",
    "GZMK/GZMA intermediate compartment on the memory-effector continuum",
    "Strong recurrent proliferation program across 17 supported patients",
    "Strong interferon response but 78.5% of cells derive from one subject",
    "Recurrent cytotoxic program across all supported patients",
    "Only five well-supported patients; strongest clonotype sharing is with CD8C03"
  )
)
fwrite(cluster_decisions, file.path(results_dir, "23_cluster_state_decisions.csv"))

theme_project <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.22)),
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

state_colors <- c(
  CD8_memory_transitional = "#3B82B8",
  CD8_dysfunctional = "#7651B9",
  CD8_cytotoxic = "#D95F59",
  CD8_proliferating = "#2F9E73"
)

# ------------------------------------------------------------------------------
# Load and validate
# ------------------------------------------------------------------------------

if (!file.exists(input_sce_file)) stop("Missing Step 21 SCE: ", input_sce_file)
message("Loading completed Step 21 object")
sce <- readRDS(input_sce_file)

required_metadata <- c(
  "subject_id", "specimen_id", "amplification_batch", "step20_subcluster",
  "step21_cd8_cluster", "step21_strict_direct_cd8",
  "step21_subject_dominated_step20_origin", "tcr_detected",
  "tcr_productive_primary", "tcr_clonotype_key", "endogenous_n_counts",
  "endogenous_n_genes"
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
cell_id <- colnames(sce)
if (!setequal(unique(cluster), cluster_decisions$step21_cd8_cluster)) {
  stop("Unexpected Step 21 discovery-cluster labels.")
}

cell_metadata <- as.data.table(as.data.frame(
  SummarizedExperiment::colData(sce)
))
cell_metadata[, cell_id := cell_id]
setcolorder(cell_metadata, c("cell_id", setdiff(names(cell_metadata), "cell_id")))
if (anyDuplicated(names(cell_metadata)) || anyDuplicated(cell_metadata$cell_id)) {
  stop("Step 21 cell metadata contains duplicated names or cell IDs.")
}

# ------------------------------------------------------------------------------
# Build a state-focused transfer space that cannot be driven by overlays
# ------------------------------------------------------------------------------

message("Constructing an overlay-resistant expression space")
row_data <- SummarizedExperiment::rowData(sce)
required_row_data <- c(
  "step21_biological_variance", "step21_excluded_from_hvg"
)
missing_row_data <- setdiff(required_row_data, colnames(row_data))
if (length(missing_row_data)) {
  stop("Step 21 rowData is missing: ", paste(missing_row_data, collapse = ", "))
}

interferon_genes <- c(
  "ISG15", "IFIT1", "IFIT2", "IFIT3", "IFITM1", "IFI6", "IFI44L",
  "MX1", "MX2", "OAS1", "OAS2", "OAS3", "ISG20", "STAT1", "IRF7"
)
proliferation_genes <- c(
  "MKI67", "TOP2A", "STMN1", "TUBA1B", "TUBB", "HMGB2", "TYMS",
  "UBE2C", "CCNB1", "CCNB2", "CDK1", "PCNA", "MCM2", "MCM3",
  "MCM4", "MCM5", "MCM6", "MCM7"
)
stress_genes <- c(
  "FOS", "FOSB", "JUN", "JUNB", "JUND", "DUSP1", "HSPA1A",
  "HSPA1B", "DNAJB1", "IER2", "IER3"
)

gene_names <- rownames(sce)
extra_exclusion <- gene_names %chin% c(
  interferon_genes, proliferation_genes, stress_genes
) |
  grepl("^MT", gene_names) |
  grepl("^MTRNR", gene_names) |
  grepl("^RP[SL][0-9]", gene_names) |
  grepl("^HIST", gene_names) |
  gene_names %in% c("MALAT1", "GAPDH", "ACTB", "PTMA")

biological_variance <- as.numeric(row_data$step21_biological_variance)
eligible_transfer_gene <-
  !as.logical(row_data$step21_excluded_from_hvg) & !extra_exclusion &
  is.finite(biological_variance)
eligible_genes <- gene_names[eligible_transfer_gene]
if (length(eligible_genes) < n_transfer_hvgs) {
  stop("Only ", length(eligible_genes), " genes are eligible for transfer PCA.")
}
transfer_hvgs <- eligible_genes[
  order(biological_variance[eligible_transfer_gene], decreasing = TRUE)[
    seq_len(n_transfer_hvgs)
  ]
]

transfer_gene_table <- data.table(
  gene = gene_names,
  biological_variance = biological_variance,
  step21_excluded_from_hvg = as.logical(row_data$step21_excluded_from_hvg),
  excluded_overlay_or_technical = extra_exclusion,
  eligible_for_transfer = eligible_transfer_gene,
  selected_for_transfer = gene_names %chin% transfer_hvgs
)[order(-biological_variance)]
fwrite(transfer_gene_table, file.path(results_dir, "23_transfer_gene_selection.csv"))

transfer_cluster_levels <- c("CD8C01", "CD8C02", "CD8C03", "CD8C05", "CD8C06")
transfer_index <- which(cluster %in% transfer_cluster_levels)
anchor_index <- which(cluster %in% c("CD8C01", "CD8C02", "CD8C03", "CD8C06"))
overlay_index <- which(cluster == "CD8C05")
if (!length(anchor_index) || !length(overlay_index)) {
  stop("Anchor or interferon-overlay cells are missing.")
}

logcounts_transfer <- SummarizedExperiment::assay(
  sce, "logcounts"
)[transfer_hvgs, transfer_index, drop = FALSE]
logcounts_transfer <- tryCatch(
  methods::as(logcounts_transfer, "dgCMatrix"),
  error = function(e) stop("Could not coerce transfer expression: ",
                           conditionMessage(e))
)
pca_input <- Matrix::t(logcounts_transfer)
set.seed(random_seed)
pca_result <- BiocSingular::runPCA(
  pca_input, rank = n_transfer_pcs, BSPARAM = BiocSingular::IrlbaParam()
)
rm(logcounts_transfer, pca_input)
invisible(gc())

if (nrow(pca_result$x) != length(transfer_index) ||
    ncol(pca_result$x) != n_transfer_pcs) {
  stop("Unexpected Step 23 transfer-PCA dimensions.")
}
rownames(pca_result$x) <- cell_id[transfer_index]
colnames(pca_result$x) <- paste0("PC", seq_len(n_transfer_pcs))

transfer_position <- match(seq_len(ncol(sce)), transfer_index)
anchor_position <- transfer_position[anchor_index]
overlay_position <- transfer_position[overlay_index]
if (anyNA(anchor_position) || anyNA(overlay_position)) {
  stop("Could not align transfer-PCA positions to Step 21 cells.")
}

anchor_state <- fifelse(
  cluster[anchor_index] == "CD8C02", "CD8_dysfunctional",
  fifelse(cluster[anchor_index] == "CD8C06", "CD8_cytotoxic",
          "CD8_memory_transitional")
)
anchor_subject <- subject[anchor_index]
overlay_subject <- subject[overlay_index]
anchor_pca <- pca_result$x[anchor_position, , drop = FALSE]
overlay_pca <- pca_result$x[overlay_position, , drop = FALSE]

# Select the nearest k neighbors from other patients and use inverse-distance
# weighted voting. Cross-patient neighbors prevent a patient-specific effect
# from masquerading as a transferable state.
vote_neighbors <- function(index, distance, query_subject, reference_subject,
                           reference_label, k) {
  if (!identical(dim(index), dim(distance))) {
    stop("Neighbor index and distance dimensions differ.")
  }
  state_order <- c(
    "CD8_memory_transitional", "CD8_dysfunctional", "CD8_cytotoxic"
  )
  out <- vector("list", nrow(index))
  for (i in seq_len(nrow(index))) {
    candidate <- index[i, ]
    keep <- is.finite(candidate) & candidate > 0L &
      reference_subject[candidate] != query_subject[i]
    candidate <- candidate[keep]
    candidate_distance <- distance[i, ][keep]
    if (length(candidate) < k) {
      stop("Fewer than ", k, " cross-patient neighbors for query row ", i)
    }
    candidate <- candidate[seq_len(k)]
    candidate_distance <- candidate_distance[seq_len(k)]
    labels <- reference_label[candidate]
    weights <- 1 / pmax(candidate_distance, 1e-8)
    scores <- setNames(rep(0, length(state_order)), state_order)
    for (state_name in state_order) {
      scores[state_name] <- sum(weights[labels == state_name])
    }
    probabilities <- scores / sum(scores)
    winner <- state_order[which.max(probabilities[state_order])]
    out[[i]] <- data.table(
      knn_state = winner,
      knn_confidence = unname(probabilities[winner]),
      knn_memory_probability = probabilities["CD8_memory_transitional"],
      knn_dysfunction_probability = probabilities["CD8_dysfunctional"],
      knn_cytotoxic_probability = probabilities["CD8_cytotoxic"],
      kth_neighbor_distance = max(candidate_distance)
    )
  }
  rbindlist(out)
}

# ------------------------------------------------------------------------------
# Cross-patient validation of expression transfer
# ------------------------------------------------------------------------------

message("Validating anchor-state transfer across patients")
validation_candidates <- min(k_candidates + 1L, nrow(anchor_pca))
set.seed(random_seed)
anchor_knn <- BiocNeighbors::queryKNN(
  X = anchor_pca, query = anchor_pca, k = validation_candidates,
  BNPARAM = BiocNeighbors::AnnoyParam()
)
anchor_validation <- vote_neighbors(
  anchor_knn$index, anchor_knn$distance, anchor_subject, anchor_subject,
  anchor_state, k_vote
)
anchor_validation[, `:=`(
  cell_id = cell_id[anchor_index],
  subject_id = anchor_subject,
  discovery_cluster = cluster[anchor_index],
  reference_state = anchor_state,
  correct = knn_state == anchor_state
)]
setcolorder(anchor_validation, c(
  "cell_id", "subject_id", "discovery_cluster", "reference_state",
  "knn_state", "correct", setdiff(names(anchor_validation), c(
    "cell_id", "subject_id", "discovery_cluster", "reference_state",
    "knn_state", "correct"
  ))
))
fwrite(anchor_validation,
       file.path(results_dir, "23_cross_patient_transfer_validation.tsv.gz"),
       sep = "\t")

validation_summary <- anchor_validation[, .(
  n_cells = .N,
  n_correct = sum(correct),
  accuracy_percent = 100 * mean(correct),
  median_confidence = median(knn_confidence)
), by = reference_state][order(match(reference_state, final_state_levels))]
fwrite(validation_summary,
       file.path(results_dir, "23_cross_patient_transfer_summary.csv"))

validation_confusion <- anchor_validation[, .N,
  by = .(reference_state, predicted_state = knn_state)
][order(reference_state, predicted_state)]
fwrite(validation_confusion,
       file.path(results_dir, "23_cross_patient_transfer_confusion.csv"))

# ------------------------------------------------------------------------------
# Assign CD8C05 by clonotype first, then cross-patient expression neighbors
# ------------------------------------------------------------------------------

message("Inferring base states for the CD8C05 interferon overlay")
set.seed(random_seed)
overlay_knn_raw <- BiocNeighbors::queryKNN(
  X = anchor_pca, query = overlay_pca,
  k = min(k_candidates, nrow(anchor_pca)),
  BNPARAM = BiocNeighbors::AnnoyParam()
)
overlay_knn <- vote_neighbors(
  overlay_knn_raw$index, overlay_knn_raw$distance,
  overlay_subject, anchor_subject, anchor_state, k_vote
)
overlay_knn[, cell_id := cell_id[overlay_index]]

anchor_clone_cells <- data.table(
  clonotype_key = as.character(sce$tcr_clonotype_key[anchor_index]),
  productive = as.logical(sce$tcr_productive_primary[anchor_index]),
  state = anchor_state
)[productive == TRUE & !is.na(clonotype_key) & nzchar(clonotype_key)]

clone_state_counts <- anchor_clone_cells[, .N, by = .(clonotype_key, state)]
clone_totals <- clone_state_counts[, .(
  n_anchor_cells = sum(N),
  maximum_state_cells = max(N),
  n_states_observed = .N
), by = clonotype_key]
clone_state_counts <- merge(
  clone_state_counts, clone_totals, by = "clonotype_key", all.x = TRUE
)
setorder(clone_state_counts, clonotype_key, -N, state)
clone_consensus <- clone_state_counts[, .(
  clone_state = state[1L],
  clone_state_cells = N[1L],
  n_anchor_cells = n_anchor_cells[1L],
  n_states_observed = n_states_observed[1L],
  clone_consensus_fraction = N[1L] / n_anchor_cells[1L],
  unique_top_state = sum(N == max(N)) == 1L
), by = clonotype_key]
clone_consensus[, usable_clone_consensus :=
  n_anchor_cells >= minimum_clone_anchor_cells & unique_top_state &
    clone_consensus_fraction >= minimum_clone_consensus]
fwrite(clone_consensus,
       file.path(results_dir, "23_anchor_clonotype_state_consensus.csv"))

overlay_assignment <- data.table(
  cell_id = cell_id[overlay_index],
  subject_id = overlay_subject,
  clonotype_key = as.character(sce$tcr_clonotype_key[overlay_index]),
  productive_tcr = as.logical(sce$tcr_productive_primary[overlay_index])
)
overlay_assignment <- merge(
  overlay_assignment, overlay_knn, by = "cell_id", all.x = TRUE, sort = FALSE
)
overlay_assignment <- merge(
  overlay_assignment, clone_consensus, by = "clonotype_key",
  all.x = TRUE, sort = FALSE
)
overlay_assignment[, clone_supported :=
  productive_tcr == TRUE & !is.na(usable_clone_consensus) &
    usable_clone_consensus]
overlay_assignment[, inferred_base_state := fifelse(
  clone_supported, clone_state, knn_state
)]
overlay_assignment[, assignment_method := fifelse(
  clone_supported, "productive_clonotype_consensus", "cross_patient_knn"
)]
overlay_assignment[, assignment_confidence := fifelse(
  clone_supported, clone_consensus_fraction, knn_confidence
)]
overlay_assignment[, overlay_order__ := match(cell_id, cell_id[overlay_index])]
setorder(overlay_assignment, overlay_order__)
overlay_assignment[, overlay_order__ := NULL]
if (nrow(overlay_assignment) != length(overlay_index) ||
    anyDuplicated(overlay_assignment$cell_id) ||
    anyNA(overlay_assignment$inferred_base_state)) {
  stop("CD8C05 overlay assignment failed accounting or completeness checks.")
}
fwrite(overlay_assignment,
       file.path(results_dir, "23_interferon_overlay_assignments.tsv.gz"),
       sep = "\t")

overlay_assignment_summary <- overlay_assignment[, .(
  n_cells = .N,
  percent_overlay_cells = 100 * .N / nrow(overlay_assignment),
  median_assignment_confidence = median(assignment_confidence),
  n_subjects = uniqueN(subject_id)
), by = .(inferred_base_state, assignment_method)][
  order(match(inferred_base_state, final_state_levels), assignment_method)
]
fwrite(overlay_assignment_summary,
       file.path(results_dir, "23_interferon_overlay_assignment_summary.csv"))

# ------------------------------------------------------------------------------
# Freeze final cell labels
# ------------------------------------------------------------------------------

message("Freezing four final CD8 states")
final_state <- rep(NA_character_, ncol(sce))
assignment_method <- rep(NA_character_, ncol(sce))
assignment_confidence <- rep(1, ncol(sce))
interferon_high <- cluster == "CD8C05"

final_state[cluster %in% c("CD8C01", "CD8C03")] <-
  "CD8_memory_transitional"
assignment_method[cluster %in% c("CD8C01", "CD8C03")] <-
  "patient_supported_cluster_merge"

final_state[cluster == "CD8C02"] <- "CD8_dysfunctional"
assignment_method[cluster == "CD8C02"] <- "patient_supported_state_anchor"
final_state[cluster == "CD8C04"] <- "CD8_proliferating"
assignment_method[cluster == "CD8C04"] <- "patient_supported_state_anchor"
final_state[cluster == "CD8C06"] <- "CD8_cytotoxic"
assignment_method[cluster == "CD8C06"] <- "patient_supported_state_anchor"

final_state[cluster == "CD8C07"] <- "CD8_memory_transitional"
assignment_method[cluster == "CD8C07"] <- "subject_dominated_cluster_merge"

overlay_match <- match(cell_id[overlay_index], overlay_assignment$cell_id)
final_state[overlay_index] <- overlay_assignment$inferred_base_state[overlay_match]
assignment_method[overlay_index] <- overlay_assignment$assignment_method[overlay_match]
assignment_confidence[overlay_index] <-
  overlay_assignment$assignment_confidence[overlay_match]

if (anyNA(final_state) || any(!final_state %in% final_state_levels)) {
  stop("One or more cells lack a valid final CD8 state.")
}

sce$step23_cd8_state <- factor(final_state, levels = final_state_levels)
sce$step23_assignment_method <- assignment_method
sce$step23_assignment_confidence <- assignment_confidence
sce$step23_interferon_high <- interferon_high
sce$step23_state_anchor_cell <- cluster %in% c("CD8C01", "CD8C02", "CD8C03", "CD8C06")
sce$step23_signature_training_eligible <- !interferon_high & cluster != "CD8C07"

# ------------------------------------------------------------------------------
# Final-state audits
# ------------------------------------------------------------------------------

cell_metadata <- as.data.table(as.data.frame(
  SummarizedExperiment::colData(sce)
))
cell_metadata[, cell_id := colnames(sce)]
setcolorder(cell_metadata, c("cell_id", setdiff(names(cell_metadata), "cell_id")))

state_summary <- cell_metadata[, .(
  n_cells = .N,
  percent_cd8_cells = 100 * .N / ncol(sce),
  n_subjects = uniqueN(subject_id),
  n_subjects_with_at_least_10_cells = sum(table(subject_id) >=
                                           minimum_patient_cells_for_support),
  n_specimens = uniqueN(specimen_id),
  n_amplification_batches = uniqueN(amplification_batch),
  dominant_subject_fraction = max(table(subject_id)) / .N,
  n_tcr = sum(tcr_detected),
  percent_tcr = 100 * mean(tcr_detected),
  n_productive_tcr = sum(tcr_productive_primary),
  unique_productive_clonotypes = uniqueN(
    tcr_clonotype_key[tcr_productive_primary == TRUE &
                        !is.na(tcr_clonotype_key) & nzchar(tcr_clonotype_key)]
  ),
  n_interferon_high = sum(step23_interferon_high),
  percent_interferon_high = 100 * mean(step23_interferon_high),
  n_signature_training_eligible = sum(step23_signature_training_eligible),
  median_endogenous_umis = as.numeric(median(endogenous_n_counts)),
  median_endogenous_genes = as.numeric(median(endogenous_n_genes))
), by = step23_cd8_state][order(match(step23_cd8_state, final_state_levels))]
fwrite(state_summary, file.path(results_dir, "23_final_state_summary.csv"))

patient_state_counts <- cell_metadata[, .(n_cells = .N),
  by = .(subject_id, step23_cd8_state)
]
patient_totals <- patient_state_counts[, .(subject_cd8_cells = sum(n_cells)),
                                      by = subject_id]
patient_state_counts <- merge(
  patient_state_counts, patient_totals, by = "subject_id", all.x = TRUE
)
patient_state_counts[, percent_subject_cd8 := 100 * n_cells / subject_cd8_cells]
patient_state_counts[, state_order__ := match(
  as.character(step23_cd8_state), final_state_levels
)]
setorder(patient_state_counts, subject_id, state_order__)
patient_state_counts[, state_order__ := NULL]
fwrite(patient_state_counts,
       file.path(results_dir, "23_final_states_by_subject.csv"))

program_columns <- grep(
  "^step21_score_", colnames(SummarizedExperiment::colData(sce)), value = TRUE
)
if (!length(program_columns)) stop("No Step 21 program scores are available.")
state_program_summary <- cell_metadata[, lapply(.SD, mean, na.rm = TRUE),
  by = step23_cd8_state, .SDcols = program_columns
]
state_program_summary <- melt(
  state_program_summary, id.vars = "step23_cd8_state",
  variable.name = "program", value.name = "mean_z_score"
)
state_program_summary[, program := sub("^step21_score_", "", program)]
state_program_summary[, state_order__ := match(
  as.character(step23_cd8_state), final_state_levels
)]
setorder(state_program_summary, state_order__, program)
state_program_summary[, state_order__ := NULL]
fwrite(state_program_summary,
       file.path(results_dir, "23_final_state_program_scores.csv"))

final_clone_membership <- unique(cell_metadata[
  tcr_productive_primary == TRUE & !is.na(tcr_clonotype_key) &
    nzchar(tcr_clonotype_key),
  .(tcr_clonotype_key, step23_cd8_state)
])
final_clone_sharing <- CJ(
  state_1 = final_state_levels, state_2 = final_state_levels, unique = TRUE
)[state_1 <= state_2][, {
  s1 <- as.character(state_1[1L])
  s2 <- as.character(state_2[1L])
  c1 <- final_clone_membership[
    step23_cd8_state == s1, unique(tcr_clonotype_key)
  ]
  c2 <- final_clone_membership[
    step23_cd8_state == s2, unique(tcr_clonotype_key)
  ]
  shared <- length(intersect(c1, c2))
  union_n <- length(union(c1, c2))
  .(
    n_clonotypes_state_1 = length(c1),
    n_clonotypes_state_2 = length(c2),
    n_shared_clonotypes = shared,
    clonotype_jaccard = if (union_n) shared / union_n else NA_real_
  )
}, by = .(state_1, state_2)]
fwrite(final_clone_sharing,
       file.path(results_dir, "23_final_state_clonotype_sharing.csv"))

# ------------------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------------------

umap <- SingleCellExperiment::reducedDim(sce, "UMAP_step21")
set.seed(random_seed)
plot_index <- if (ncol(sce) > maximum_plot_cells) {
  sort(sample.int(ncol(sce), maximum_plot_cells))
} else seq_len(ncol(sce))
umap_plot <- data.table(
  UMAP1 = umap[plot_index, 1L],
  UMAP2 = umap[plot_index, 2L],
  final_state = factor(final_state[plot_index], levels = final_state_levels),
  interferon_high = interferon_high[plot_index]
)

p_state <- ggplot(umap_plot, aes(UMAP1, UMAP2, color = final_state)) +
  geom_point(size = 0.42, alpha = 0.72) +
  scale_color_manual(values = state_colors, drop = FALSE) +
  labs(
    title = "Final tumor-CD8 states in GSE123139",
    subtitle = "Four states frozen after patient-balanced and clonotype-aware validation",
    color = "Final CD8 state"
  ) + theme_project()
save_png("23A_umap_final_cd8_states.png", p_state, 11, 8)

p_patient <- ggplot(
  patient_state_counts,
  aes(x = subject_id, y = percent_subject_cd8, fill = step23_cd8_state)
) +
  geom_col(width = 0.85) +
  scale_fill_manual(values = state_colors, drop = FALSE) +
  labs(
    title = "Final CD8-state composition across patients",
    subtitle = "Percentages are calculated within each patient's retained tumor-CD8 compartment",
    x = "Patient", y = "% of patient CD8 cells", fill = "Final CD8 state"
  ) + theme_project() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_png("23B_final_state_composition_by_patient.png", p_patient, 13, 7)

p_program <- ggplot(
  state_program_summary,
  aes(x = program, y = step23_cd8_state, fill = mean_z_score)
) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "#3B6FB6", mid = "white", high = "#7651B9", midpoint = 0,
    name = "Mean gene-wise\nz-score"
  ) +
  labs(
    title = "Expression programs in the four final CD8 states",
    subtitle = "Interferon remains an overlay and is not a fifth state",
    x = NULL, y = NULL
  ) + theme_project() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_png("23C_final_state_program_heatmap.png", p_program, 12, 6)

overlay_umap <- data.table(
  UMAP1 = umap[overlay_index, 1L],
  UMAP2 = umap[overlay_index, 2L],
  inferred_state = factor(
    final_state[overlay_index], levels = final_state_levels[1:3]
  ),
  method = assignment_method[overlay_index],
  confidence = assignment_confidence[overlay_index]
)
p_overlay <- ggplot(
  overlay_umap, aes(UMAP1, UMAP2, color = inferred_state, shape = method)
) +
  geom_point(size = 1.15, alpha = 0.78) +
  scale_color_manual(values = state_colors, drop = FALSE) +
  labs(
    title = "Base-state assignments within the interferon overlay",
    subtitle = "Productive clonotype consensus is used before cross-patient expression neighbors",
    color = "Inferred base state", shape = "Assignment method"
  ) + theme_project()
save_png("23D_interferon_overlay_base_states.png", p_overlay, 10, 8)

p_validation <- ggplot(
  anchor_validation,
  aes(x = reference_state, fill = knn_state)
) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  scale_fill_manual(values = state_colors[final_state_levels[1:3]], drop = FALSE) +
  labs(
    title = "Cross-patient validation of base-state transfer",
    subtitle = paste0(k_vote, " nearest neighbors; all neighbors from other patients"),
    x = "Reference anchor state", y = "Predicted fraction",
    fill = "Predicted state"
  ) + theme_project() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
save_png("23E_cross_patient_transfer_validation.png", p_validation, 10, 7)

# ------------------------------------------------------------------------------
# Save and report
# ------------------------------------------------------------------------------

if (sum(state_summary$n_cells) != ncol(sce)) {
  stop("Final-state accounting does not match the SCE cell count.")
}
if (uniqueN(as.character(sce$step23_cd8_state)) != length(final_state_levels)) {
  stop("Not all four final states are represented.")
}
if (any(state_summary$n_subjects_with_at_least_10_cells < 10L)) {
  warning("At least one final state has fewer than 10 patients with >=10 cells.")
}
if (any(interferon_high & final_state == "CD8_proliferating")) {
  stop("The interferon overlay was unexpectedly assigned to proliferation.")
}

message("Saving final-state SingleCellExperiment")
saveRDS(sce, output_sce_file, compress = FALSE)
output_size_mb <- unname(file.info(output_sce_file)$size / 1024^2)

overall_validation_accuracy <- 100 * mean(anchor_validation$correct)
dataset_qc <- data.table(
  metric = c(
    "cells", "subjects", "final_cd8_states", "memory_transitional_cells",
    "dysfunctional_cells", "cytotoxic_cells", "proliferating_cells",
    "interferon_overlay_cells", "interferon_clone_assigned_cells",
    "interferon_knn_assigned_cells", "transfer_hvgs", "transfer_pcs",
    "cross_patient_knn_k", "cross_patient_anchor_validation_accuracy_percent",
    "all_cells_labeled", "final_state_accounting_validation",
    "output_sce_size_mb"
  ),
  value = as.character(c(
    ncol(sce), uniqueN(subject), length(final_state_levels),
    sum(final_state == "CD8_memory_transitional"),
    sum(final_state == "CD8_dysfunctional"),
    sum(final_state == "CD8_cytotoxic"),
    sum(final_state == "CD8_proliferating"), sum(interferon_high),
    sum(interferon_high & assignment_method == "productive_clonotype_consensus"),
    sum(interferon_high & assignment_method == "cross_patient_knn"),
    length(transfer_hvgs), n_transfer_pcs, k_vote,
    overall_validation_accuracy, !anyNA(final_state),
    sum(state_summary$n_cells) == ncol(sce), output_size_mb
  ))
)
fwrite(dataset_qc, file.path(results_dir, "23_dataset_qc.csv"))

figure_manifest <- data.table(
  figure = c(
    "23A_umap_final_cd8_states.png",
    "23B_final_state_composition_by_patient.png",
    "23C_final_state_program_heatmap.png",
    "23D_interferon_overlay_base_states.png",
    "23E_cross_patient_transfer_validation.png"
  ),
  purpose = c(
    "Final four-state annotation in the Step 21 embedding",
    "Patient support and within-patient state composition",
    "Major expression programs after state merging",
    "Audit of interferon-overlay base-state inference",
    "Cross-patient validation of the expression transfer rule"
  )
)
fwrite(figure_manifest, file.path(results_dir, "23_figure_manifest.csv"))
writeLines(capture.output(sessionInfo()), file.path(results_dir, "23_sessionInfo.txt"))

message("Step 23 final CD8-state assignment complete.")
message("Final states: ", paste(final_state_levels, collapse = ", "))
message("Cross-patient anchor validation accuracy: ",
        round(overall_validation_accuracy, 2), "%")
message("Output SCE: ", output_sce_file)
print(dataset_qc)
print(state_summary)
print(overlay_assignment_summary)