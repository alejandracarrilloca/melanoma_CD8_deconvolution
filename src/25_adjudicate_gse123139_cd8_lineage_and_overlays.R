#!/usr/bin/env Rscript

# ==============================================================================
# Step 25: Lineage adjudication and overlay-aware CD8 state consolidation
# ==============================================================================
#
# Step 24 showed that the selected tumor-CD8 compartment contains reproducible
# memory/transitional, dysfunctional, and proliferating biology, but that its
# broad inclusion rule also admits cells with non-T evidence. It also showed
# that the cytotoxic cluster is not a reproducible standalone base state.
#
# This step therefore:
#   1. replaces broad non-T flags with lineage-specific anchor evidence;
#   2. audits CD37/CD74 separately so they cannot define B-cell identity;
#   3. validates CD8-versus-alternative-lineage transfer across patients;
#   4. retains only high-confidence anchors or conservative cross-patient
#      rescues;
#   5. consolidates the retained compartment into memory/transitional and
#      dysfunctional base states; and
#   6. preserves cytotoxic, interferon, and proliferation as non-exclusive
#      overlays.
#
# No bulk-deconvolution signature is exported by this step.
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
# Paths and fixed analysis parameters
# ------------------------------------------------------------------------------

project_root <- normalizePath(getwd(), mustWork = TRUE)
analysis_dir <- file.path(
  project_root, "data", "analysis_ready", "scRNA", "GSE123139"
)
results_dir <- file.path(project_root, "results", "25")
figures_dir <- file.path(results_dir, "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

input_sce_file <- file.path(
  analysis_dir, "24_GSE123139_tumor_cd8_robustness_audited_sce.rds"
)
output_full_sce_file <- file.path(
  analysis_dir, "25_GSE123139_tumor_cd8_lineage_adjudicated_sce.rds"
)
output_purified_sce_file <- file.path(
  analysis_dir, "25_GSE123139_tumor_cd8_purified_sce.rds"
)

random_seed <- 123139L
lineage_knn_k <- 30L
state_knn_k <- 30L
state_hvg_count <- 1500L
minimum_anchor_cells <- 100L
minimum_anchor_subjects <- 5L
minimum_clone_cells <- 2L
minimum_clone_consensus <- 0.60
minimum_state_confidence <- 0.50

base_state_levels <- c("CD8_memory_transitional", "CD8_dysfunctional")
base_state_colors <- c(
  CD8_memory_transitional = "#3B82B8",
  CD8_dysfunctional = "#7651B9"
)
lineage_decision_colors <- c(
  high_confidence_cd8_anchor = "#2678B2",
  cross_patient_cd8_rescue = "#61A6D8",
  alternative_lineage_anchor = "#D95F59",
  conflicting_lineage_evidence = "#E69F00",
  unresolved_not_retained = "#C9C9C9"
)

# Genes used to define direct lineage evidence. CD37 and CD74 are intentionally
# excluded from the B-cell anchors and audited only as contextual genes.
lineage_genes <- list(
  T_core = c("CD3D", "CD3E", "CD3G", "TRAC", "CD247", "LCK", "LAT", "TRBC1", "TRBC2"),
  CD8 = c("CD8A", "CD8B"),
  B_specific = c("MS4A1", "CD79A", "CD79B", "CD19", "CD22"),
  B_context = c("CD37", "CD74", "HLA-DRA"),
  Myeloid_specific = c("LST1", "LILRB1", "CSF1R", "S100A8", "S100A9", "AIF1", "CTSS"),
  Melanoma_specific = c("MLANA", "PMEL", "SOX10", "MITF", "TYR", "DCT"),
  NK_specific = c("KLRD1", "FCGR3A", "FCER1G", "TYROBP", "NCR1", "NCR3"),
  CD4_Treg_specific = c("CD4", "FOXP3", "IL2RA", "CTLA4", "IKZF2")
)

overlay_genes <- list(
  Interferon = c("ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "STAT1", "IFI6", "IFI44L"),
  Cytotoxic = c("NKG7", "GNLY", "PRF1", "GZMB", "GZMH", "FGFBP2", "CTSW"),
  Proliferation = c("MKI67", "TOP2A", "STMN1", "TUBA1B", "HMGB2", "TYMS")
)

theme_project <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.20)),
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

as_flag <- function(x) !is.na(x) & as.logical(x)

# ------------------------------------------------------------------------------
# Load and validate the completed Step 24 object
# ------------------------------------------------------------------------------

if (!file.exists(input_sce_file)) stop("Missing Step 24 SCE: ", input_sce_file)
message("Loading completed Step 24 object")
sce <- readRDS(input_sce_file)

required_metadata <- c(
  "subject_id", "specimen_id", "amplification_batch", "step20_direct_T",
  "step20_direct_CD8", "step21_strict_direct_cd8", "step21_cd8_cluster",
  "step23_cd8_state", "step23_interferon_high",
  "step23_signature_training_eligible", "step24_t_supported",
  "step24_potential_non_t", "tcr_detected", "tcr_productive_primary",
  "tcr_clonotype_key", "endogenous_n_counts", "endogenous_n_genes"
)
missing_metadata <- setdiff(
  required_metadata, colnames(SummarizedExperiment::colData(sce))
)
if (length(missing_metadata)) {
  stop("Step 24 metadata is missing: ", paste(missing_metadata, collapse = ", "))
}
if (!all(c("counts", "logcounts") %in% SummarizedExperiment::assayNames(sce))) {
  stop("Step 24 must contain counts and logcounts assays.")
}
if (!all(c("PCA_step21", "UMAP_step21") %in%
         SingleCellExperiment::reducedDimNames(sce))) {
  stop("Step 21 PCA and UMAP coordinates must be preserved in Step 24.")
}

counts <- SummarizedExperiment::assay(sce, "counts")
logcounts <- SummarizedExperiment::assay(sce, "logcounts")
counts <- tryCatch(methods::as(counts, "dgCMatrix"), error = function(e) counts)
logcounts <- tryCatch(
  methods::as(logcounts, "dgCMatrix"),
  error = function(e) stop("Could not coerce logcounts: ", conditionMessage(e))
)

subject <- as.character(sce$subject_id)
step23_state <- as.character(sce$step23_cd8_state)
productive_tcr <- as_flag(sce$tcr_productive_primary)
old_strict_cd8 <- as_flag(sce$step21_strict_direct_cd8)

# ------------------------------------------------------------------------------
# Direct, lineage-specific evidence and the CD37/CD74 audit
# ------------------------------------------------------------------------------

message("Auditing exact genes behind Step 24 lineage flags")
detected_gene_count <- function(genes) {
  present <- intersect(genes, rownames(sce))
  if (!length(present)) return(rep.int(0L, ncol(sce)))
  as.integer(Matrix::colSums(counts[present, , drop = FALSE] > 0))
}

detected_counts <- lapply(lineage_genes, detected_gene_count)
t_core_direct <- detected_counts$T_core >= 2L
cd8_detected <- detected_counts$CD8 >= 1L
b_specific_direct <- detected_counts$B_specific >= 2L
b_context_any <- detected_counts$B_context >= 1L
myeloid_specific_direct <- detected_counts$Myeloid_specific >= 2L
melanoma_specific_direct <- detected_counts$Melanoma_specific >= 2L
nk_specific_direct <- detected_counts$NK_specific >= 2L
cd4_treg_specific_direct <- detected_counts$CD4_Treg_specific >= 2L

t_positive_evidence <- t_core_direct | productive_tcr | old_strict_cd8
specific_non_t_evidence <-
  b_specific_direct | myeloid_specific_direct | melanoma_specific_direct |
  nk_specific_direct
cd4_without_cd8 <- cd4_treg_specific_direct & !cd8_detected
lineage_conflict <-
  (specific_non_t_evidence | cd4_without_cd8) & t_positive_evidence

positive_anchor <-
  (old_strict_cd8 | ((productive_tcr | t_core_direct) & cd8_detected)) &
  !specific_non_t_evidence & !cd4_without_cd8
negative_anchor <-
  (specific_non_t_evidence | cd4_without_cd8) & !t_positive_evidence

if (sum(positive_anchor) < minimum_anchor_cells ||
    data.table::uniqueN(subject[positive_anchor]) < minimum_anchor_subjects) {
  stop("Too few cross-patient positive CD8 anchors.")
}
if (sum(negative_anchor) < minimum_anchor_cells ||
    data.table::uniqueN(subject[negative_anchor]) < minimum_anchor_subjects) {
  stop("Too few cross-patient alternative-lineage anchors.")
}

sce$step25_t_core_direct <- t_core_direct
sce$step25_cd8_detected <- cd8_detected
sce$step25_b_specific_direct <- b_specific_direct
sce$step25_b_context_any <- b_context_any
sce$step25_myeloid_specific_direct <- myeloid_specific_direct
sce$step25_melanoma_specific_direct <- melanoma_specific_direct
sce$step25_nk_specific_direct <- nk_specific_direct
sce$step25_cd4_treg_specific_direct <- cd4_treg_specific_direct
sce$step25_lineage_conflict <- lineage_conflict

anchor_class <- rep("ambiguous", ncol(sce))
anchor_class[positive_anchor] <- "positive_cd8_anchor"
anchor_class[negative_anchor] <- "alternative_lineage_anchor"
anchor_class[lineage_conflict] <- "conflicting_lineage_evidence"
sce$step25_lineage_anchor_class <- anchor_class

# ------------------------------------------------------------------------------
# Leave-one-patient-out lineage distances and conservative rescue threshold
# ------------------------------------------------------------------------------

message("Validating CD8 lineage transfer across patients")
pca <- SingleCellExperiment::reducedDim(sce, "PCA_step21")
if (ncol(pca) < 30L) stop("PCA_step21 contains fewer than 30 components.")
pca <- as.matrix(pca[, seq_len(30L), drop = FALSE])

cross_patient_mean_distance <- function(query_matrix, query_subject,
                                        reference_matrix, reference_subject,
                                        k, seed_offset = 0L) {
  output <- rep(NA_real_, nrow(query_matrix))
  subject_levels <- sort(unique(query_subject))
  for (i in seq_along(subject_levels)) {
    sid <- subject_levels[[i]]
    query_index <- which(query_subject == sid)
    reference_index <- which(reference_subject != sid)
    k_use <- min(as.integer(k), length(reference_index))
    if (k_use < 5L) stop("Insufficient cross-patient reference anchors for ", sid)
    set.seed(random_seed + seed_offset + i)
    neighbor_result <- BiocNeighbors::queryKNN(
      X = reference_matrix[reference_index, , drop = FALSE],
      query = query_matrix[query_index, , drop = FALSE],
      k = k_use,
      BNPARAM = BiocNeighbors::AnnoyParam()
    )
    output[query_index] <- rowMeans(neighbor_result$distance)
  }
  output
}

positive_distance <- cross_patient_mean_distance(
  pca, subject, pca[positive_anchor, , drop = FALSE],
  subject[positive_anchor], lineage_knn_k, 1000L
)
negative_distance <- cross_patient_mean_distance(
  pca, subject, pca[negative_anchor, , drop = FALSE],
  subject[negative_anchor], lineage_knn_k, 2000L
)
lineage_score <- log((negative_distance + 1e-8) / (positive_distance + 1e-8))

# At most five percent of the known alternative-lineage anchors should cross
# the rescue boundary. A score of zero is the equal-distance boundary.
lineage_rescue_threshold <- max(
  0, as.numeric(stats::quantile(
    lineage_score[negative_anchor], probs = 0.95, na.rm = TRUE, names = FALSE
  ))
)
cross_patient_rescue <-
  anchor_class == "ambiguous" & is.finite(lineage_score) &
  lineage_score >= lineage_rescue_threshold

lineage_decision <- rep("unresolved_not_retained", ncol(sce))
lineage_decision[positive_anchor] <- "high_confidence_cd8_anchor"
lineage_decision[cross_patient_rescue] <- "cross_patient_cd8_rescue"
lineage_decision[negative_anchor] <- "alternative_lineage_anchor"
lineage_decision[lineage_conflict] <- "conflicting_lineage_evidence"
retained_cd8 <- lineage_decision %in% c(
  "high_confidence_cd8_anchor", "cross_patient_cd8_rescue"
)

lineage_confidence <- abs(stats::plogis(lineage_score) - 0.5) * 2
sce$step25_lineage_positive_distance <- positive_distance
sce$step25_lineage_negative_distance <- negative_distance
sce$step25_lineage_score <- lineage_score
sce$step25_lineage_confidence <- lineage_confidence
sce$step25_lineage_decision <- lineage_decision
sce$step25_retained_cd8 <- retained_cd8

validation_index <- positive_anchor | negative_anchor
validation_reference <- ifelse(
  positive_anchor[validation_index], "CD8", "alternative_lineage"
)
validation_prediction <- ifelse(
  lineage_score[validation_index] >= lineage_rescue_threshold,
  "CD8", "alternative_lineage"
)
lineage_validation <- data.table(
  subject_id = subject[validation_index],
  reference_class = validation_reference,
  predicted_class = validation_prediction,
  lineage_score = lineage_score[validation_index]
)
lineage_validation_summary <- lineage_validation[, .(
  n_cells = as.integer(.N),
  n_correct = as.integer(sum(reference_class == predicted_class)),
  accuracy_percent = 100 * mean(reference_class == predicted_class),
  median_lineage_score = stats::median(lineage_score, na.rm = TRUE)
), by = reference_class]
lineage_sensitivity <- mean(validation_prediction[validation_reference == "CD8"] == "CD8")
lineage_specificity <- mean(
  validation_prediction[validation_reference == "alternative_lineage"] ==
    "alternative_lineage"
)
lineage_balanced_accuracy <- 100 * mean(c(lineage_sensitivity, lineage_specificity))

# ------------------------------------------------------------------------------
# Overlay-resistant expression space and two-state transfer
# ------------------------------------------------------------------------------

message("Consolidating retained cells into two base states")
retained_index <- which(retained_cd8)
if (length(retained_index) < 1000L) stop("Too few cells survived lineage adjudication.")

row_data_names <- colnames(SummarizedExperiment::rowData(sce))
if (all(c("step21_biological_variance", "step21_excluded_from_hvg") %in%
        row_data_names)) {
  biological_variance <- as.numeric(
    SummarizedExperiment::rowData(sce)$step21_biological_variance
  )
  eligible_gene <-
    is.finite(biological_variance) & biological_variance > 0 &
    !as_flag(SummarizedExperiment::rowData(sce)$step21_excluded_from_hvg)
} else {
  gene_mean <- Matrix::rowMeans(logcounts[, retained_index, drop = FALSE])
  biological_variance <- pmax(
    Matrix::rowMeans(logcounts[, retained_index, drop = FALSE] ^ 2) - gene_mean ^ 2,
    0
  )
  eligible_gene <- is.finite(biological_variance) & biological_variance > 0
}

excluded_state_space_genes <- unique(c(
  unlist(overlay_genes, use.names = FALSE),
  unlist(lineage_genes, use.names = FALSE),
  grep("^(MT-|RPL|RPS)", rownames(sce), value = TRUE),
  c("MALAT1", "NEAT1", "FOS", "JUN", "JUNB", "DUSP1", "HSPA1A", "HSPA1B")
))
eligible_gene <- eligible_gene & !(rownames(sce) %in% excluded_state_space_genes)
eligible_index <- which(eligible_gene)
if (length(eligible_index) < state_hvg_count) {
  stop("Too few eligible genes for overlay-resistant state transfer.")
}
hvg_order <- eligible_index[order(biological_variance[eligible_index], decreasing = TRUE)]
state_hvgs <- rownames(sce)[head(hvg_order, state_hvg_count)]

set.seed(random_seed + 3000L)
state_pca_result <- BiocSingular::runPCA(
  Matrix::t(logcounts[state_hvgs, retained_index, drop = FALSE]),
  rank = 30L, BSPARAM = BiocSingular::IrlbaParam()
)
state_pca <- state_pca_result$x
rownames(state_pca) <- colnames(sce)[retained_index]
colnames(state_pca) <- paste0("PC", seq_len(ncol(state_pca)))
state_pca_full <- matrix(
  NA_real_, nrow = ncol(sce), ncol = ncol(state_pca),
  dimnames = list(colnames(sce), colnames(state_pca))
)
state_pca_full[retained_index, ] <- state_pca
SingleCellExperiment::reducedDim(sce, "PCA_step25_overlay_resistant") <-
  state_pca_full

state_anchor <-
  retained_cd8 & positive_anchor &
  step23_state %in% base_state_levels &
  as_flag(sce$step23_signature_training_eligible) &
  !as_flag(sce$step23_interferon_high)
state_anchor_counts <- table(factor(
  step23_state[state_anchor], levels = base_state_levels
))
if (any(state_anchor_counts < minimum_anchor_cells)) {
  stop("A base state has too few high-confidence transfer anchors.")
}

cross_patient_state_vote <- function(query_matrix, query_subject,
                                     reference_matrix, reference_subject,
                                     reference_label, k) {
  predicted <- rep(NA_character_, nrow(query_matrix))
  confidence <- rep(NA_real_, nrow(query_matrix))
  subject_levels <- sort(unique(query_subject))
  for (i in seq_along(subject_levels)) {
    sid <- subject_levels[[i]]
    query_index <- which(query_subject == sid)
    reference_index <- which(reference_subject != sid)
    k_use <- min(as.integer(k), length(reference_index))
    if (k_use < 5L) stop("Insufficient cross-patient state anchors for ", sid)
    set.seed(random_seed + 4000L + i)
    neighbor_result <- BiocNeighbors::queryKNN(
      X = reference_matrix[reference_index, , drop = FALSE],
      query = query_matrix[query_index, , drop = FALSE],
      k = k_use,
      BNPARAM = BiocNeighbors::AnnoyParam()
    )
    neighbor_label <- matrix(
      reference_label[reference_index][neighbor_result$index],
      nrow = length(query_index), ncol = k_use
    )
    neighbor_weight <- 1 / (neighbor_result$distance + 1e-8)
    memory_vote <- rowSums(
      neighbor_weight * (neighbor_label == "CD8_memory_transitional")
    )
    dysfunctional_vote <- rowSums(
      neighbor_weight * (neighbor_label == "CD8_dysfunctional")
    )
    total_vote <- memory_vote + dysfunctional_vote
    predicted[query_index] <- ifelse(
      memory_vote >= dysfunctional_vote,
      "CD8_memory_transitional", "CD8_dysfunctional"
    )
    confidence[query_index] <- pmax(memory_vote, dysfunctional_vote) / total_vote
  }
  list(predicted = predicted, confidence = confidence)
}

retained_subject <- subject[retained_index]
anchor_position <- match(which(state_anchor), retained_index)
state_vote <- cross_patient_state_vote(
  state_pca, retained_subject,
  state_pca[anchor_position, , drop = FALSE],
  retained_subject[anchor_position], step23_state[state_anchor], state_knn_k
)

anchor_validation_position <- which(state_anchor[retained_index])
state_transfer_validation <- data.table(
  subject_id = retained_subject[anchor_validation_position],
  reference_state = step23_state[retained_index][anchor_validation_position],
  predicted_state = state_vote$predicted[anchor_validation_position],
  confidence = state_vote$confidence[anchor_validation_position]
)
state_transfer_summary <- state_transfer_validation[, .(
  n_cells = as.integer(.N),
  n_correct = as.integer(sum(reference_state == predicted_state)),
  accuracy_percent = 100 * mean(reference_state == predicted_state),
  median_confidence = stats::median(confidence, na.rm = TRUE)
), by = reference_state]

# Productive clonotypes spanning high-confidence base anchors provide the first
# assignment route for cells whose original Step 23 state is now an overlay.
anchor_clone_dt <- data.table(
  clonotype = as.character(sce$tcr_clonotype_key[state_anchor]),
  base_state = step23_state[state_anchor]
)[!is.na(clonotype) & nzchar(clonotype)]
clone_consensus <- anchor_clone_dt[, {
  state_counts <- table(base_state)
  top_state <- names(state_counts)[which.max(state_counts)]
  .(
    n_anchor_cells = as.integer(sum(state_counts)),
    consensus_state = as.character(top_state),
    consensus_fraction = as.numeric(max(state_counts) / sum(state_counts))
  )
}, by = clonotype]
clone_consensus <- clone_consensus[
  n_anchor_cells >= minimum_clone_cells &
    consensus_fraction >= minimum_clone_consensus
]

base_state <- rep(NA_character_, ncol(sce))
base_assignment_method <- rep(NA_character_, ncol(sce))
base_assignment_confidence <- rep(NA_real_, ncol(sce))

existing_base <- retained_cd8 & step23_state %in% base_state_levels
base_state[existing_base] <- step23_state[existing_base]
base_assignment_method[existing_base] <- "retained_step23_base_state"
base_assignment_confidence[existing_base] <- 1

overlay_target <- retained_cd8 & !existing_base
clone_match <- match(
  as.character(sce$tcr_clonotype_key), clone_consensus$clonotype
)
clone_assignable <- overlay_target & !is.na(clone_match)
base_state[clone_assignable] <-
  clone_consensus$consensus_state[clone_match[clone_assignable]]
base_assignment_method[clone_assignable] <- "productive_clonotype_consensus"
base_assignment_confidence[clone_assignable] <-
  clone_consensus$consensus_fraction[clone_match[clone_assignable]]

knn_target <- overlay_target & is.na(base_state)
knn_position <- match(which(knn_target), retained_index)
base_state[knn_target] <- state_vote$predicted[knn_position]
base_assignment_method[knn_target] <- "cross_patient_knn"
base_assignment_confidence[knn_target] <- state_vote$confidence[knn_position]

sce$step25_base_state <- base_state
sce$step25_base_assignment_method <- base_assignment_method
sce$step25_base_assignment_confidence <- base_assignment_confidence
sce$step25_interferon_overlay <- retained_cd8 & as_flag(sce$step23_interferon_high)
sce$step25_cytotoxic_overlay <- retained_cd8 & step23_state == "CD8_cytotoxic"
sce$step25_proliferation_overlay <- retained_cd8 &
  step23_state == "CD8_proliferating"
sce$step25_signature_training_eligible <-
  retained_cd8 & positive_anchor & step23_state %in% base_state_levels &
  !sce$step25_interferon_overlay & !sce$step25_cytotoxic_overlay &
  !sce$step25_proliferation_overlay &
  base_assignment_confidence >= minimum_state_confidence

# ------------------------------------------------------------------------------
# Tables and hard accounting checks
# ------------------------------------------------------------------------------

message("Writing lineage, transfer, and overlay QC")
cell_metadata <- as.data.table(SummarizedExperiment::colData(sce))
cell_metadata[, cell_id := colnames(sce)]

lineage_summary <- cell_metadata[, .(
  n_cells = as.integer(.N),
  percent_all_cells = 100 * .N / nrow(cell_metadata),
  n_subjects = as.integer(uniqueN(subject_id)),
  median_lineage_score = stats::median(step25_lineage_score, na.rm = TRUE),
  median_lineage_confidence = stats::median(step25_lineage_confidence, na.rm = TRUE)
), by = step25_lineage_decision]

decision_by_prior_state <- cell_metadata[, .(
  n_cells = as.integer(.N)
), by = .(step23_cd8_state, step25_lineage_decision)]
decision_by_prior_state[, percent_prior_state :=
  100 * n_cells / sum(n_cells), by = step23_cd8_state]

specific_evidence_by_state <- cell_metadata[, .(
  n_cells = as.integer(.N),
  percent_t_core = 100 * mean(step25_t_core_direct),
  percent_cd8_detected = 100 * mean(step25_cd8_detected),
  percent_productive_tcr = 100 * mean(as_flag(tcr_productive_primary)),
  percent_b_specific = 100 * mean(step25_b_specific_direct),
  percent_b_context_any = 100 * mean(step25_b_context_any),
  percent_myeloid_specific = 100 * mean(step25_myeloid_specific_direct),
  percent_melanoma_specific = 100 * mean(step25_melanoma_specific_direct),
  percent_nk_specific = 100 * mean(step25_nk_specific_direct),
  percent_cd4_treg_specific = 100 * mean(step25_cd4_treg_specific_direct),
  percent_lineage_conflict = 100 * mean(step25_lineage_conflict),
  percent_retained = 100 * mean(step25_retained_cd8)
), by = step23_cd8_state]

b_context_category <- ifelse(
  b_specific_direct & b_context_any, "B-specific plus context",
  ifelse(
    b_specific_direct, "B-specific without context",
    ifelse(b_context_any, "CD37/CD74/HLA-DRA context only", "Neither")
  )
)
b_context_audit <- data.table(
  step23_cd8_state = step23_state,
  category = b_context_category,
  t_positive_evidence = t_positive_evidence,
  retained_cd8 = retained_cd8
)[, .(
  n_cells = as.integer(.N),
  n_t_positive = as.integer(sum(t_positive_evidence)),
  percent_t_positive = 100 * mean(t_positive_evidence),
  n_retained = as.integer(sum(retained_cd8)),
  percent_retained = 100 * mean(retained_cd8)
), by = .(step23_cd8_state, category)]

final_base_state_summary <- cell_metadata[step25_retained_cd8 == TRUE, .(
  n_cells = as.integer(.N),
  percent_retained_cd8 = 100 * .N / sum(cell_metadata$step25_retained_cd8),
  n_subjects = as.integer(uniqueN(subject_id)),
  n_subjects_with_at_least_10_cells = as.integer(sum(table(subject_id) >= 10L)),
  n_tcr = as.integer(sum(as_flag(tcr_detected))),
  percent_tcr = 100 * mean(as_flag(tcr_detected)),
  n_strict_direct_cd8 = as.integer(sum(as_flag(step21_strict_direct_cd8))),
  percent_strict_direct_cd8 = 100 * mean(as_flag(step21_strict_direct_cd8)),
  n_signature_training_eligible = as.integer(sum(step25_signature_training_eligible)),
  median_endogenous_umis = stats::median(endogenous_n_counts, na.rm = TRUE),
  median_endogenous_genes = stats::median(endogenous_n_genes, na.rm = TRUE)
), by = step25_base_state]

overlay_summary <- rbindlist(list(
  cell_metadata[step25_retained_cd8 == TRUE, .(
    step25_base_state,
    overlay = "Interferon",
    present = step25_interferon_overlay
  )],
  cell_metadata[step25_retained_cd8 == TRUE, .(
    step25_base_state,
    overlay = "Cytotoxic",
    present = step25_cytotoxic_overlay
  )],
  cell_metadata[step25_retained_cd8 == TRUE, .(
    step25_base_state,
    overlay = "Proliferation",
    present = step25_proliferation_overlay
  )]
))[, .(
  n_cells = as.integer(sum(present)),
  percent_base_state = 100 * mean(present)
), by = .(step25_base_state, overlay)]

base_assignment_summary <- cell_metadata[step25_retained_cd8 == TRUE, .(
  n_cells = as.integer(.N),
  percent_retained_cd8 = 100 * .N / sum(cell_metadata$step25_retained_cd8),
  median_assignment_confidence = stats::median(
    step25_base_assignment_confidence, na.rm = TRUE
  ),
  n_subjects = as.integer(uniqueN(subject_id))
), by = .(step25_base_state, step25_base_assignment_method)]

lineage_accounting_validation <-
  sum(lineage_summary$n_cells) == ncol(sce) &&
  sum(retained_cd8) == sum(lineage_decision %in% c(
    "high_confidence_cd8_anchor", "cross_patient_cd8_rescue"
  ))
no_forbidden_cells_retained <- !any(retained_cd8 & (negative_anchor | lineage_conflict))
base_state_accounting_validation <-
  all(!is.na(base_state[retained_cd8])) &&
  all(is.na(base_state[!retained_cd8])) &&
  setequal(unique(base_state[retained_cd8]), base_state_levels)
if (!lineage_accounting_validation) stop("Lineage accounting validation failed.")
if (!no_forbidden_cells_retained) stop("A negative or conflicting anchor was retained.")
if (!base_state_accounting_validation) stop("Base-state accounting validation failed.")

dataset_qc <- data.table(
  metric = c(
    "input_cells", "subjects", "positive_cd8_anchor_cells",
    "alternative_lineage_anchor_cells", "lineage_conflict_cells",
    "cross_patient_rescued_cells", "retained_cd8_cells",
    "not_retained_cells", "lineage_rescue_threshold",
    "lineage_balanced_accuracy_percent", "base_states",
    "interferon_overlay_cells", "cytotoxic_overlay_cells",
    "proliferation_overlay_cells", "signature_training_eligible_cells",
    "lineage_accounting_validation", "no_forbidden_cells_retained",
    "base_state_accounting_validation", "signature_exported",
    "output_full_sce_size_mb", "output_purified_sce_size_mb"
  ),
  value = as.character(c(
    ncol(sce), uniqueN(subject), sum(positive_anchor), sum(negative_anchor),
    sum(lineage_conflict), sum(cross_patient_rescue), sum(retained_cd8),
    sum(!retained_cd8), lineage_rescue_threshold,
    lineage_balanced_accuracy, length(base_state_levels),
    sum(sce$step25_interferon_overlay), sum(sce$step25_cytotoxic_overlay),
    sum(sce$step25_proliferation_overlay),
    sum(sce$step25_signature_training_eligible),
    as.integer(lineage_accounting_validation),
    as.integer(no_forbidden_cells_retained),
    as.integer(base_state_accounting_validation), 0L, NA_real_, NA_real_
  ))
)

fwrite(dataset_qc, file.path(results_dir, "25_dataset_qc.csv"))
fwrite(lineage_summary, file.path(results_dir, "25_lineage_adjudication_summary.csv"))
fwrite(
  decision_by_prior_state,
  file.path(results_dir, "25_lineage_decision_by_step23_state.csv")
)
fwrite(
  specific_evidence_by_state,
  file.path(results_dir, "25_specific_lineage_evidence_by_step23_state.csv")
)
fwrite(b_context_audit, file.path(results_dir, "25_cd37_cd74_b_context_audit.csv"))
fwrite(
  lineage_validation,
  file.path(results_dir, "25_cross_patient_lineage_validation.csv")
)
fwrite(
  lineage_validation_summary,
  file.path(results_dir, "25_cross_patient_lineage_validation_summary.csv")
)
fwrite(
  state_transfer_validation,
  file.path(results_dir, "25_cross_patient_base_state_validation.csv")
)
fwrite(
  state_transfer_summary,
  file.path(results_dir, "25_cross_patient_base_state_validation_summary.csv")
)
fwrite(
  state_transfer_validation[, .N, by = .(reference_state, predicted_state)],
  file.path(results_dir, "25_cross_patient_base_state_confusion.csv")
)
fwrite(
  final_base_state_summary,
  file.path(results_dir, "25_final_base_state_summary.csv")
)
fwrite(base_assignment_summary, file.path(results_dir, "25_base_assignment_summary.csv"))
fwrite(overlay_summary, file.path(results_dir, "25_overlay_summary.csv"))
fwrite(
  data.table(gene = state_hvgs, rank = seq_along(state_hvgs)),
  file.path(results_dir, "25_overlay_resistant_hvgs.csv")
)

cell_export <- cell_metadata[, .(
  cell_id, subject_id, specimen_id, step21_cd8_cluster, step23_cd8_state,
  step25_lineage_anchor_class, step25_lineage_score,
  step25_lineage_confidence, step25_lineage_decision, step25_retained_cd8,
  step25_base_state, step25_base_assignment_method,
  step25_base_assignment_confidence, step25_interferon_overlay,
  step25_cytotoxic_overlay, step25_proliferation_overlay,
  step25_signature_training_eligible
)]
fwrite(
  cell_export,
  file.path(results_dir, "25_cell_lineage_and_state_assignments.tsv.gz"),
  sep = "\t"
)

# ------------------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------------------

message("Rendering Step 25 figures")
umap <- as.data.table(SingleCellExperiment::reducedDim(sce, "UMAP_step21"))
setnames(umap, c("UMAP1", "UMAP2"))
plot_dt <- cbind(
  umap,
  cell_metadata[, .(
    subject_id, step23_cd8_state, step25_lineage_anchor_class,
    step25_lineage_score, step25_lineage_decision, step25_retained_cd8,
    step25_base_state, step25_interferon_overlay,
    step25_cytotoxic_overlay, step25_proliferation_overlay
  )]
)

plot_dt[, step25_lineage_decision := factor(
  step25_lineage_decision,
  levels = names(lineage_decision_colors)
)]
p_lineage <- ggplot(plot_dt, aes(UMAP1, UMAP2, color = step25_lineage_decision)) +
  geom_point(size = 0.25, alpha = 0.72) +
  scale_color_manual(values = lineage_decision_colors, drop = FALSE) +
  labs(
    title = "Lineage adjudication of the Step 23 tumor-CD8 compartment",
    subtitle = "Only specific anchors or conservative cross-patient rescues are retained",
    color = "Step 25 decision"
  ) +
  theme_project()
save_png("25A_umap_lineage_adjudication.png", p_lineage, 11, 8)

evidence_long <- melt(
  specific_evidence_by_state,
  id.vars = c("step23_cd8_state", "n_cells"),
  measure.vars = grep("^percent_", names(specific_evidence_by_state), value = TRUE),
  variable.name = "evidence", value.name = "percent_cells"
)
p_evidence <- ggplot(
  evidence_long,
  aes(evidence, step23_cd8_state, fill = percent_cells)
) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "#F7F7FB", high = "#7651B9") +
  labs(
    title = "Specific lineage evidence by prior CD8 state",
    subtitle = "CD37/CD74 context is shown separately from B-specific anchors",
    x = NULL, y = "Step 23 state", fill = "% of cells"
  ) +
  theme_project() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_png("25B_specific_lineage_evidence.png", p_evidence, 12, 7)

score_plot_dt <- plot_dt[step25_lineage_anchor_class %in% c(
  "positive_cd8_anchor", "alternative_lineage_anchor", "ambiguous"
)]
p_score <- ggplot(
  score_plot_dt,
  aes(step25_lineage_anchor_class, step25_lineage_score,
      color = step25_lineage_anchor_class)
) +
  geom_hline(yintercept = lineage_rescue_threshold, linetype = 2) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.18, size = 0.18, alpha = 0.22) +
  guides(color = "none") +
  labs(
    title = "Cross-patient lineage scores",
    subtitle = "Dashed line is the conservative ambiguous-cell rescue threshold",
    x = NULL, y = "log(distance to alternative / distance to CD8)"
  ) +
  theme_project() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_png("25C_cross_patient_lineage_scores.png", p_score, 10, 7)

purified_plot_dt <- copy(plot_dt[step25_retained_cd8 == TRUE])
purified_plot_dt[, step25_base_state := factor(
  step25_base_state, levels = base_state_levels
)]
p_base <- ggplot(
  purified_plot_dt,
  aes(UMAP1, UMAP2, color = step25_base_state)
) +
  geom_point(size = 0.28, alpha = 0.75) +
  scale_color_manual(values = base_state_colors, drop = FALSE) +
  labs(
    title = "Purified tumor-CD8 compartment with two base states",
    subtitle = "Cytotoxic, interferon, and proliferation are retained as overlays",
    color = "Base state"
  ) +
  theme_project()
save_png("25D_umap_purified_base_states.png", p_base, 11, 8)

overlay_plot_dt <- rbindlist(list(
  purified_plot_dt[step25_interferon_overlay == TRUE, .(
    UMAP1, UMAP2, overlay = "Interferon"
  )],
  purified_plot_dt[step25_cytotoxic_overlay == TRUE, .(
    UMAP1, UMAP2, overlay = "Cytotoxic"
  )],
  purified_plot_dt[step25_proliferation_overlay == TRUE, .(
    UMAP1, UMAP2, overlay = "Proliferation"
  )]
))
p_overlay <- ggplot(purified_plot_dt, aes(UMAP1, UMAP2)) +
  geom_point(color = "#E2E2E2", size = 0.22) +
  geom_point(
    data = overlay_plot_dt,
    aes(color = overlay), size = 0.45, alpha = 0.85
  ) +
  facet_wrap(~overlay, nrow = 1) +
  scale_color_manual(values = c(
    Interferon = "#2C7FB8", Cytotoxic = "#D95F59", Proliferation = "#2F9E73"
  )) +
  guides(color = "none") +
  labs(
    title = "Non-exclusive programs over the two CD8 base states",
    subtitle = "Overlay labels do not create additional deconvolution states"
  ) +
  theme_project()
save_png("25E_umap_state_overlays.png", p_overlay, 14, 5.5)

patient_retention <- cell_metadata[, .(
  input_cells = as.integer(.N),
  retained_cells = as.integer(sum(step25_retained_cd8)),
  retained_percent = 100 * mean(step25_retained_cd8)
), by = subject_id]
p_patient <- ggplot(
  patient_retention,
  aes(reorder(as.character(subject_id), retained_percent), retained_percent)
) +
  geom_col(fill = "#3B82B8") +
  coord_flip() +
  labs(
    title = "CD8 retention after lineage adjudication by patient",
    x = "Patient", y = "% of selected cells retained"
  ) +
  theme_project()
save_png("25F_patient_retention.png", p_patient, 8, 9)
fwrite(patient_retention, file.path(results_dir, "25_patient_retention.csv"))

# ------------------------------------------------------------------------------
# Save full audit object and purified object, then finalize QC sizes
# ------------------------------------------------------------------------------

message("Saving Step 25 full and purified objects")
saveRDS(sce, output_full_sce_file, compress = FALSE)
purified_sce <- sce[, retained_cd8]
saveRDS(purified_sce, output_purified_sce_file, compress = FALSE)

dataset_qc[metric == "output_full_sce_size_mb", value := as.character(
  file.info(output_full_sce_file)$size / 1024 ^ 2
)]
dataset_qc[metric == "output_purified_sce_size_mb", value := as.character(
  file.info(output_purified_sce_file)$size / 1024 ^ 2
)]
fwrite(dataset_qc, file.path(results_dir, "25_dataset_qc.csv"))

writeLines(
  c(
    "Step 25 completed successfully.",
    paste("Input cells:", ncol(sce)),
    paste("Retained tumor-CD8 cells:", sum(retained_cd8)),
    paste("Positive CD8 anchors:", sum(positive_anchor)),
    paste("Cross-patient rescues:", sum(cross_patient_rescue)),
    paste("Alternative-lineage anchors:", sum(negative_anchor)),
    paste("Conflicting-lineage cells:", sum(lineage_conflict)),
    paste("Lineage balanced accuracy (%):", round(lineage_balanced_accuracy, 2)),
    paste("Memory/transitional base cells:", sum(base_state == "CD8_memory_transitional", na.rm = TRUE)),
    paste("Dysfunctional base cells:", sum(base_state == "CD8_dysfunctional", na.rm = TRUE)),
    paste("Interferon overlay cells:", sum(sce$step25_interferon_overlay)),
    paste("Cytotoxic overlay cells:", sum(sce$step25_cytotoxic_overlay)),
    paste("Proliferation overlay cells:", sum(sce$step25_proliferation_overlay)),
    "No signature matrix was exported."
  ),
  file.path(results_dir, "25_completion.txt")
)

message("Step 25 complete: lineage adjudicated, two base states retained, overlays preserved")
