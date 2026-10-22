#!/usr/bin/env Rscript

# ==============================================================================
# Step 36: Build an SPN-free melanoma CD8-state deconvolution reference
# ==============================================================================
# Produces a conservative reference separating memory/transitional and
# dysfunctional CD8 cells from other major melanoma compartments. Candidate
# genes must be significant in the Step 27 paired-patient pseudobulk analysis,
# show consistent direction across patients, and retain compartment specificity.
# SPN and technical/receptor genes are excluded.
# ============================================================================== 

options(stringsAsFactors = FALSE)
options(error = function() {
  cat("\n--- STEP 36 TRACEBACK ---\n", file = stderr())
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
results_dir <- file.path(project_root, "results", "36")
figures_dir <- file.path(results_dir, "figures")
export_dir <- file.path(results_dir, "cibersortx")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

full_sce_file <- file.path(
  analysis_dir, "19_GSE123139_broad_clusters_sce.rds"
)
cd8_sce_file <- file.path(
  analysis_dir, "34_GSE123139_beltra_tex_states_spn_sce.rds"
)
de_file <- file.path(results_dir, "..", "27",
                     "27_primary_paired_pseudobulk_de.csv")

for (path in c(full_sce_file, cd8_sce_file, de_file)) {
  if (!file.exists(path)) stop("Required input not found: ", path)
}

minimum_cells_per_patient_state <- 10L
minimum_direction_consistency <- 0.70
maximum_de_fdr <- 0.01
minimum_abs_log2_fc <- 0.50
maximum_markers_per_cd8_state <- 50L
minimum_markers_per_cd8_state <- 10L

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

# Frozen anchors ensure that the non-CD8 compartments are identifiable rather
# than being represented only by low expression of CD8-state markers.
compartment_anchor_genes <- list(
  Other_T_NK = c(
    "CD3D", "CD3E", "TRAC", "CD4", "IL7R", "KLRD1", "NCR1", "GNLY",
    "FCGR3A", "TYROBP"
  ),
  B_cell = c(
    "CD79A", "CD79B", "MS4A1", "CD19", "CD22", "CD74", "CD37",
    "HLA-DRA", "CD83"
  ),
  Myeloid = c(
    "LST1", "FCER1G", "TYROBP", "AIF1", "CTSS", "LILRB1", "S100A8",
    "S100A9", "CTSD", "LYZ"
  ),
  Melanoma = c("MLANA", "PMEL", "SOX10", "MITF", "TYR", "DCT")
)

technical_gene <- function(gene) {
  gene == "SPN" |
    grepl("^(MT-|MTRNR|MTATP|MTND|MRPL|MRPS)", gene) |
    grepl("^(RPL|RPS|RP[0-9]+-|SNOR|SNHG|MIR|hsa-mir)", gene) |
    grepl("^(TRAV|TRAJ|TRAC|TRBV|TRBJ|TRBC|IGHA|IGHG|IGHM|IGKC|IGLC)", gene)
}

message("Loading full atlas, finalized CD8 object, and paired pseudobulk DE")
full_sce <- readRDS(full_sce_file)
cd8_sce <- readRDS(cd8_sce_file)
de <- fread(de_file)

if (!all(c("counts", "logcounts") %in%
         SummarizedExperiment::assayNames(full_sce))) {
  stop("Step 19 atlas must contain counts and logcounts assays.")
}
required_full_metadata <- c("subject_id", "source_group", "broad_cluster")
missing_full <- setdiff(
  required_full_metadata, colnames(SummarizedExperiment::colData(full_sce))
)
if (length(missing_full)) stop("Missing Step 19 metadata: ", paste(missing_full, collapse = ", "))

required_cd8_metadata <- c(
  "subject_id", "step25_retained_cd8", "step25_base_state",
  "step25_signature_training_eligible"
)
missing_cd8 <- setdiff(
  required_cd8_metadata, colnames(SummarizedExperiment::colData(cd8_sce))
)
if (length(missing_cd8)) stop("Missing CD8 metadata: ", paste(missing_cd8, collapse = ", "))

required_de_columns <- c("gene", "log2_fold_change", "mean_log_cpm", "fdr")
if (!all(required_de_columns %in% names(de))) {
  stop("Step 27 DE table has unexpected columns.")
}
if (anyDuplicated(colnames(full_sce)) || anyDuplicated(colnames(cd8_sce))) {
  stop("Duplicated cell identifiers occur.")
}

common_genes <- intersect(rownames(full_sce), rownames(cd8_sce))
if (length(common_genes) < 1000L) stop("Too few common genes between objects.")

cd8_dt <- as.data.table(SummarizedExperiment::colData(cd8_sce))
cd8_dt[, cell_id := colnames(cd8_sce)]
cd8_dt[, subject_id := as.character(subject_id)]
cd8_dt[, base_state := as.character(step25_base_state)]
cd8_dt[, retained := !is.na(step25_retained_cd8) & as.logical(step25_retained_cd8)]
cd8_dt[, training_eligible := !is.na(step25_signature_training_eligible) &
          as.logical(step25_signature_training_eligible)]
training_index <- which(cd8_dt$retained & cd8_dt$training_eligible &
                          cd8_dt$base_state %chin%
                            c("CD8_memory_transitional", "CD8_dysfunctional"))
if (!length(training_index)) stop("No eligible CD8 signature-training cells.")

cd8_log <- SummarizedExperiment::assay(cd8_sce, "logcounts")
candidate_de <- de[
  gene %chin% common_genes & fdr <= maximum_de_fdr &
    abs(log2_fold_change) >= minimum_abs_log2_fc
]
candidate_de[, technical_gene := technical_gene(gene)]
candidate_de <- candidate_de[technical_gene == FALSE]
if (!nrow(candidate_de)) stop("No Step 27 genes passed the initial DE filters.")

# Cross-patient direction consistency among training-eligible CD8 cells.
candidate_genes <- unique(candidate_de$gene)
cd8_expr <- as.matrix(cd8_log[candidate_genes, training_index, drop = FALSE])
training_meta <- cd8_dt[training_index, .(subject_id, base_state)]

patient_state_means <- rbindlist(lapply(seq_along(candidate_genes), function(i) {
  data.table(
    subject_id = training_meta$subject_id,
    base_state = training_meta$base_state,
    gene = candidate_genes[i],
    expression = as.numeric(cd8_expr[i, ])
  )[, .(n_cells = .N, mean_expression = mean(expression)),
    by = .(gene, subject_id, base_state)][n_cells >= minimum_cells_per_patient_state]
}))

patient_wide <- dcast(
  patient_state_means,
  gene + subject_id ~ base_state,
  value.var = "mean_expression"
)
patient_wide <- patient_wide[
  is.finite(CD8_memory_transitional) & is.finite(CD8_dysfunctional)
]
patient_wide[, patient_difference :=
               CD8_dysfunctional - CD8_memory_transitional]
direction_summary <- patient_wide[, .(
  paired_patients = .N,
  percent_higher_in_dysfunctional = 100 * mean(patient_difference > 0),
  median_patient_difference = median(patient_difference)
), by = gene]

candidate_de <- merge(candidate_de, direction_summary, by = "gene", all.x = TRUE)
candidate_de[, expected_consistency := fifelse(
  log2_fold_change > 0,
  percent_higher_in_dysfunctional / 100,
  1 - percent_higher_in_dysfunctional / 100
)]

# Build atlas labels. Finalized CD8 labels replace the corresponding T/NK cells;
# all remaining tumor T/NK cells form a nuisance reference compartment.
full_dt <- as.data.table(SummarizedExperiment::colData(full_sce))
full_dt[, cell_id := colnames(full_sce)]
full_dt[, subject_id := as.character(subject_id)]
full_dt[, source_group := as.character(source_group)]
full_dt[, broad_cluster := as.character(broad_cluster)]
full_dt[, broad_compartment := unname(broad_compartment_map[broad_cluster])]
full_dt[, reference_class := fifelse(
  broad_compartment == "T_or_NK", "Other_T_NK", broad_compartment
)]

cd8_label_map <- cd8_dt[
  retained == TRUE & base_state %chin%
    c("CD8_memory_transitional", "CD8_dysfunctional"),
  .(cell_id, cd8_reference_class = base_state)
]
full_dt <- merge(full_dt, cd8_label_map, by = "cell_id", all.x = TRUE, sort = FALSE)
full_dt[!is.na(cd8_reference_class), reference_class := cd8_reference_class]

reference_cells <- full_dt[
  source_group == "Tumor" & reference_class %chin% reference_levels
]
full_index <- match(reference_cells$cell_id, colnames(full_sce))
if (anyNA(full_index)) stop("Failed to align reference cells to Step 19 atlas.")

# Compartment specificity is evaluated on patient-balanced mean log-expression.
full_log <- SummarizedExperiment::assay(full_sce, "logcounts")
specificity_expr <- as.matrix(
  full_log[candidate_genes, full_index, drop = FALSE]
)
specificity_rows <- rbindlist(lapply(seq_along(candidate_genes), function(i) {
  data.table(
    subject_id = reference_cells$subject_id,
    reference_class = reference_cells$reference_class,
    expression = as.numeric(specificity_expr[i, ])
  )[, .(mean_expression = mean(expression)),
    by = .(subject_id, reference_class)][, .(
      median_patient_expression = median(mean_expression)
    ), by = reference_class][, gene := candidate_genes[i]]
}))
specificity_wide <- dcast(
  specificity_rows, gene ~ reference_class,
  value.var = "median_patient_expression"
)
for (column in reference_levels) {
  if (!column %in% names(specificity_wide)) specificity_wide[, (column) := NA_real_]
}
non_cd8_columns <- c("Other_T_NK", "B_cell", "Myeloid", "Melanoma")
specificity_wide[, max_non_cd8_expression := do.call(
  pmax, c(.SD, list(na.rm = TRUE))
), .SDcols = non_cd8_columns]
specificity_wide[, dysfunctional_specificity_margin :=
                   CD8_dysfunctional - max_non_cd8_expression]
specificity_wide[, memory_specificity_margin :=
                   CD8_memory_transitional - max_non_cd8_expression]

candidate_de <- merge(candidate_de, specificity_wide, by = "gene", all.x = TRUE)
candidate_de[, target_state := fifelse(
  log2_fold_change > 0, "CD8_dysfunctional", "CD8_memory_transitional"
)]
candidate_de[, target_specificity_margin := fifelse(
  target_state == "CD8_dysfunctional", dysfunctional_specificity_margin,
  memory_specificity_margin
)]
candidate_de[, passes_direction :=
               expected_consistency >= minimum_direction_consistency]
candidate_de[, passes_specificity := target_specificity_margin > 0]
candidate_de[, eligible_marker := passes_direction & passes_specificity]

selected_markers <- candidate_de[eligible_marker == TRUE][
  order(target_state, fdr, -abs(log2_fold_change),
        -target_specificity_margin), head(.SD, maximum_markers_per_cd8_state),
  by = target_state
]
marker_counts <- selected_markers[, .N, by = target_state]
if (any(!c("CD8_memory_transitional", "CD8_dysfunctional") %chin%
        marker_counts$target_state) || any(marker_counts$N < minimum_markers_per_cd8_state)) {
  fwrite(candidate_de, file.path(results_dir, "36_candidate_marker_audit.csv"))
  stop("Too few robust markers for one or both CD8 states. Inspect candidate audit.")
}

present_anchor_genes <- intersect(
  unique(unlist(compartment_anchor_genes)), rownames(full_sce)
)
selected_genes <- unique(c(selected_markers$gene, present_anchor_genes))
gene_roles <- rbindlist(list(
  selected_markers[, .(
    gene, reference_role = target_state, source = "Paired-pseudobulk CD8 marker"
  )],
  rbindlist(lapply(names(compartment_anchor_genes), function(class_name) {
    data.table(
      gene = intersect(compartment_anchor_genes[[class_name]], rownames(full_sce)),
      reference_role = class_name,
      source = "Frozen compartment anchor"
    )
  }))
), use.names = TRUE)[!duplicated(gene)]

# Reference values are medians of patient-specific class means on the original
# linear normalized-expression scale. This prevents large patients dominating.
selected_log <- as.matrix(full_log[selected_genes, full_index, drop = FALSE])
selected_linear <- pmax(2^selected_log - 1, 0)
patient_class_reference <- rbindlist(lapply(seq_along(selected_genes), function(i) {
  data.table(
    gene = selected_genes[i],
    subject_id = reference_cells$subject_id,
    reference_class = reference_cells$reference_class,
    expression = as.numeric(selected_linear[i, ])
  )[, .(mean_expression = mean(expression)),
    by = .(gene, subject_id, reference_class)]
}))
reference_long <- patient_class_reference[, .(
  expression = median(mean_expression),
  contributing_patients = uniqueN(subject_id)
), by = .(gene, reference_class)]
reference_matrix <- dcast(
  reference_long, gene ~ reference_class, value.var = "expression", fill = 0
)
setcolorder(reference_matrix, c("gene", reference_levels))
setnames(reference_matrix, "gene", "GeneSymbol")

# A compact cell-level reference is also exported for CIBERSORTx signature
# creation. Every class is capped per patient to reduce imbalance.
set.seed(123139L)
reference_cells[, keep_for_export := FALSE]
reference_cells[, keep_for_export := seq_len(.N) %in%
                  sample.int(.N, min(.N, 100L)),
                by = .(subject_id, reference_class)]
export_cells <- reference_cells[keep_for_export == TRUE]
export_index <- match(export_cells$cell_id, colnames(full_sce))
export_matrix <- as.matrix(
  SummarizedExperiment::assay(full_sce, "counts")[
    selected_genes, export_index, drop = FALSE
  ]
)
export_table <- data.table(GeneSymbol = selected_genes)
export_table <- cbind(export_table, as.data.table(export_matrix))
setnames(export_table, c("GeneSymbol", export_cells$cell_id))
phenotype_table <- export_cells[, .(cell_id, reference_class, subject_id)]

message("Writing Step 36 QC and CIBERSORTx exports")
dataset_qc <- data.table(
  metric = c(
    "full_atlas_tumor_cells", "reference_cells", "reference_patients",
    "candidate_de_genes", "eligible_robust_markers",
    "memory_transitional_markers", "dysfunctional_markers",
    "compartment_anchor_genes", "total_reference_genes",
    "cell_level_export_cells", "spn_in_reference", "signature_exported"
  ),
  value = c(
    sum(full_dt$source_group == "Tumor"), nrow(reference_cells),
    uniqueN(reference_cells$subject_id), nrow(candidate_de),
    sum(candidate_de$eligible_marker),
    sum(selected_markers$target_state == "CD8_memory_transitional"),
    sum(selected_markers$target_state == "CD8_dysfunctional"),
    length(present_anchor_genes), length(selected_genes),
    nrow(export_cells), as.integer("SPN" %chin% selected_genes), 1L
  )
)
fwrite(dataset_qc, file.path(results_dir, "36_dataset_qc.csv"))
fwrite(candidate_de, file.path(results_dir, "36_candidate_marker_audit.csv"))
fwrite(selected_markers, file.path(results_dir, "36_selected_cd8_state_markers.csv"))
fwrite(gene_roles, file.path(results_dir, "36_reference_gene_roles.csv"))
fwrite(reference_long, file.path(results_dir, "36_reference_long.csv"))
fwrite(reference_matrix,
       file.path(export_dir, "36_cibersortx_signature_matrix.tsv"), sep = "\t")
fwrite(export_table,
       file.path(export_dir, "36_cibersortx_single_cell_counts.tsv.gz"), sep = "\t")
fwrite(phenotype_table,
       file.path(export_dir, "36_cibersortx_cell_phenotypes.tsv"), sep = "\t")

plot_dt <- selected_markers[, .(
  gene, target_state, log2_fold_change, expected_consistency,
  target_specificity_margin
)]
p <- ggplot(
  plot_dt,
  aes(log2_fold_change, reorder(gene, log2_fold_change), color = target_state)
) +
  geom_vline(xintercept = 0, color = "#888888") +
  geom_point(aes(size = expected_consistency), alpha = 0.85) +
  facet_wrap(~target_state, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c(
    CD8_memory_transitional = "#3B82B8", CD8_dysfunctional = "#7651B9"
  ), guide = "none") +
  labs(
    title = "Robust SPN-free markers for melanoma CD8-state deconvolution",
    subtitle = "Paired-patient DE, cross-patient direction, and compartment specificity required",
    x = "Dysfunctional minus memory/transitional log2 fold change",
    y = NULL, size = "Directional\nconsistency"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.background = element_blank(), strip.text = element_text(face = "bold")
  )
ggsave(
  file.path(figures_dir, "36A_selected_cd8_state_markers.png"),
  p, width = 12, height = 11, units = "in", dpi = 180, bg = "white"
)

if ("SPN" %chin% selected_genes) stop("SPN entered the Step 36 reference.")
message("Step 36 complete: SPN-free CD8-state deconvolution reference exported")
