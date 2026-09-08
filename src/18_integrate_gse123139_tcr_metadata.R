#!/usr/bin/env Rscript

# ==============================================================================
# Step 18: Integrate GSE123139 TCR-beta metadata with the normalized SCE
#
# The GEO TCR file has 47 tab-delimited header fields but 48 fields per data
# row. The unnamed first field duplicates Well_ID. This script assigns it the
# name record_id, validates record_id == Well_ID, and maps retained TCR records
# to the Step 17 SingleCellExperiment by Well_ID. Amp.Batch and Patient are
# independently cross-checked against the SCE metadata.
#
# Original mc_grp labels are retained for external comparison only. They are
# not used here to define new expression-based cell identities.
# ==============================================================================

options(stringsAsFactors = FALSE)

required_packages <- c(
  "data.table", "ggplot2", "SingleCellExperiment",
  "SummarizedExperiment", "S4Vectors"
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

input_sce_file <- file.path(
  "data", "analysis_ready", "scRNA", "GSE123139",
  "17_GSE123139_normalized_pca_sce.rds"
)
tcr_file <- file.path(
  "data", "scRNA", "GSE123139",
  "GSE123139_T_cells_tcrb_v2.txt.gz"
)
output_sce_file <- file.path(
  "data", "analysis_ready", "scRNA", "GSE123139",
  "18_GSE123139_normalized_pca_tcr_sce.rds"
)
results_dir <- file.path("results", "18")
figures_dir <- file.path(results_dir, "figures")
dir.create(dirname(output_sce_file), recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_sce_file)) stop("Missing Step 17 SCE: ", input_sce_file)
if (!file.exists(tcr_file)) stop("Missing TCR file: ", tcr_file)

# ------------------------------------------------------------------------------
# Read and validate the irregular TCR table
# ------------------------------------------------------------------------------

message("Reading TCR-beta metadata")
header_connection <- gzfile(tcr_file, open = "rt")
header_line <- readLines(header_connection, n = 1L, warn = FALSE)
close(header_connection)

header_fields <- strsplit(header_line, "\t", fixed = TRUE)[[1L]]
if (length(header_fields) != 47L) {
  stop("Expected 47 tab-delimited TCR header fields, found ", length(header_fields))
}

tcr <- fread(
  tcr_file,
  sep = "\t",
  header = FALSE,
  skip = 1L,
  quote = "",
  na.strings = c("NA", ""),
  fill = FALSE,
  check.names = FALSE,
  showProgress = FALSE
)

if (ncol(tcr) != 48L) {
  stop("Expected 48 tab-delimited fields per TCR data row, found ", ncol(tcr))
}
tcr_input_data_fields <- ncol(tcr)
setnames(tcr, c("record_id", header_fields))

required_tcr_columns <- c(
  "record_id", "Well_ID", "NKI.Plate.ID", "Patient", "clone_id",
  "reads_freq", "CDR3_first", "CDR3_translation_first", "Amp.Batch",
  "mc", "mc_grp", "seq_batch"
)
missing_tcr_columns <- setdiff(required_tcr_columns, names(tcr))
if (length(missing_tcr_columns)) {
  stop("Missing required TCR columns: ", paste(missing_tcr_columns, collapse = ", "))
}

tcr[, record_id := as.character(record_id)]
tcr[, Well_ID := as.character(Well_ID)]
tcr[, reported_patient_label := as.character(Patient)]
tcr[, derived_subject_id := sub(
  "^([pP][0-9]+).*",
  "\\1",
  reported_patient_label,
  perl = TRUE
)]
if (anyNA(tcr$record_id) || anyNA(tcr$Well_ID)) {
  stop("Missing record_id or Well_ID values occur in the TCR table.")
}
if (any(tcr$record_id != tcr$Well_ID)) {
  stop("The unnamed TCR record identifier does not always equal Well_ID.")
}
if (anyDuplicated(tcr$Well_ID)) {
  duplicate_ids <- unique(tcr$Well_ID[duplicated(tcr$Well_ID)])
  fwrite(
    tcr[Well_ID %chin% duplicate_ids],
    file.path(results_dir, "18_duplicated_tcr_well_ids.tsv.gz"),
    sep = "\t"
  )
  stop("Duplicated Well_ID values occur in the TCR table; see diagnostic output.")
}

# ------------------------------------------------------------------------------
# Map TCR records to the final normalized cells
# ------------------------------------------------------------------------------

message("Loading Step 17 SingleCellExperiment")
sce <- readRDS(input_sce_file)
if (anyDuplicated(colnames(sce))) stop("Duplicated cell IDs occur in the SCE.")

required_sce_metadata <- c(
  "source_group", "subject_id", "specimen_id", "facs_gate",
  "amplification_batch"
)
missing_sce_metadata <- setdiff(
  required_sce_metadata,
  colnames(SummarizedExperiment::colData(sce))
)
if (length(missing_sce_metadata)) {
  stop("Missing required SCE metadata: ", paste(missing_sce_metadata, collapse = ", "))
}

tcr[, retained_in_step17 := Well_ID %chin% colnames(sce)]
tcr_retained <- tcr[retained_in_step17 == TRUE]
if (!nrow(tcr_retained)) stop("No TCR Well_ID values map to the Step 17 SCE.")

sce_index <- match(tcr_retained$Well_ID, colnames(sce))
if (anyNA(sce_index)) stop("Internal TCR-to-SCE mapping failure.")

amp_mismatch <- as.character(tcr_retained$Amp.Batch) !=
  as.character(sce$amplification_batch[sce_index])
subject_mismatch <- as.character(tcr_retained$derived_subject_id) !=
  as.character(sce$subject_id[sce_index])
amp_mismatch[is.na(amp_mismatch)] <- TRUE
subject_mismatch[is.na(subject_mismatch)] <- TRUE

if (any(amp_mismatch) || any(subject_mismatch)) {
  mapping_mismatches <- tcr_retained[amp_mismatch | subject_mismatch]
  mapping_mismatches[, sce_amplification_batch :=
    as.character(sce$amplification_batch[sce_index[amp_mismatch | subject_mismatch]])
  ]
  mapping_mismatches[, sce_subject_id :=
    as.character(sce$subject_id[sce_index[amp_mismatch | subject_mismatch]])
  ]
  fwrite(
    mapping_mismatches,
    file.path(results_dir, "18_tcr_mapping_mismatches.tsv.gz"),
    sep = "\t"
  )
  stop("TCR mapping failed Amp.Batch or derived-subject cross-validation.")
}

# Replace any stale mismatch table from an earlier failed validation with an
# empty, schema-preserving table so downstream inspection cannot confuse it
# with the successful rerun.
successful_mapping_audit <- copy(tcr_retained[0L])
successful_mapping_audit[, sce_amplification_batch := character()]
successful_mapping_audit[, sce_subject_id := character()]
fwrite(
  successful_mapping_audit,
  file.path(results_dir, "18_tcr_mapping_mismatches.tsv.gz"),
  sep = "\t"
)

valid_text <- function(x) {
  x <- as.character(x)
  !is.na(x) & nzchar(x) & x != "NA"
}

tcr[, productive_primary :=
  valid_text(CDR3_translation_first) &
  !grepl("\\*", as.character(CDR3_translation_first))
]
tcr[, clonotype_key := NA_character_]
tcr[productive_primary == TRUE, clonotype_key := paste(
  derived_subject_id, CDR3_translation_first, sep = "::"
)]
tcr[productive_primary == FALSE & valid_text(CDR3_first), clonotype_key := paste(
  derived_subject_id, CDR3_first, sep = "::"
)]

tcr_retained <- tcr[retained_in_step17 == TRUE]
tcr_retained[, retained_clonotype_size := .N, by = clonotype_key]
tcr_retained[is.na(clonotype_key), retained_clonotype_size := NA_integer_]

mapping_index <- match(colnames(sce), tcr_retained$Well_ID)
sce$tcr_detected <- !is.na(mapping_index)
sce$tcr_productive_primary <- FALSE
sce$tcr_productive_primary[sce$tcr_detected] <-
  tcr_retained$productive_primary[mapping_index[sce$tcr_detected]]
sce$tcr_clone_id <- as.character(tcr_retained$clone_id[mapping_index])
sce$tcr_clonotype_key <- as.character(tcr_retained$clonotype_key[mapping_index])
sce$tcr_clone_size_retained <- as.integer(
  tcr_retained$retained_clonotype_size[mapping_index]
)
sce$tcr_reads_freq <- as.numeric(tcr_retained$reads_freq[mapping_index])
sce$tcr_reported_patient_label <- as.character(
  tcr_retained$reported_patient_label[mapping_index]
)
sce$tcr_derived_subject_id <- as.character(
  tcr_retained$derived_subject_id[mapping_index]
)
sce$tcr_cdr3_nt_primary <- as.character(tcr_retained$CDR3_first[mapping_index])
sce$tcr_cdr3_aa_primary <- as.character(
  tcr_retained$CDR3_translation_first[mapping_index]
)
sce$tcr_original_mc <- as.character(tcr_retained$mc[mapping_index])
sce$tcr_original_mc_grp <- as.character(tcr_retained$mc_grp[mapping_index])
sce$tcr_seq_batch <- as.character(tcr_retained$seq_batch[mapping_index])

if (sum(sce$tcr_detected) != nrow(tcr_retained)) {
  stop("The number of mapped SCE cells does not equal retained TCR records.")
}

# ------------------------------------------------------------------------------
# Audit outputs
# ------------------------------------------------------------------------------

cell_metadata <- as.data.table(as.data.frame(
  SummarizedExperiment::colData(sce)
), keep.rownames = "cell_id")

tcr_by_source <- cell_metadata[, .(
  n_cells = .N,
  n_tcr_detected = sum(tcr_detected),
  percent_tcr_detected = 100 * mean(tcr_detected),
  n_productive_primary = sum(tcr_productive_primary),
  percent_productive_primary = 100 * mean(tcr_productive_primary)
), by = source_group][order(source_group)]

tcr_by_facs_gate <- cell_metadata[, .(
  n_cells = .N,
  n_tcr_detected = sum(tcr_detected),
  percent_tcr_detected = 100 * mean(tcr_detected)
), by = .(source_group, facs_gate)][order(source_group, -percent_tcr_detected)]

original_state_summary <- cell_metadata[tcr_detected == TRUE, .(
  n_cells = .N,
  n_subjects = uniqueN(subject_id),
  n_specimens = uniqueN(specimen_id),
  n_amplification_batches = uniqueN(amplification_batch),
  n_productive_primary = sum(tcr_productive_primary),
  percent_productive_primary = 100 * mean(tcr_productive_primary)
), by = .(source_group, tcr_original_mc_grp)][order(
  source_group, -n_cells, tcr_original_mc_grp
)]

clone_summary <- tcr_retained[!is.na(clonotype_key), .(
  n_cells = .N,
  subject_id = first(as.character(derived_subject_id)),
  reported_patient_labels = paste(
    sort(unique(as.character(reported_patient_label))), collapse = ";"
  ),
  cdr3_aa_primary = first(as.character(CDR3_translation_first)),
  n_original_states = uniqueN(mc_grp),
  original_states = paste(sort(unique(mc_grp)), collapse = ";")
), by = clonotype_key][order(-n_cells, subject_id, clonotype_key)]

mapping_summary <- data.table(
  metric = c(
    "tcr_header_fields", "tcr_data_fields", "tcr_records",
    "unique_tcr_well_ids", "record_id_equals_well_id",
    "tcr_records_retained_in_step17", "tcr_records_excluded_by_qc",
    "step17_cells", "percent_step17_cells_with_tcr",
    "retained_primary_productive_tcr", "amp_batch_mismatches",
    "derived_subject_mismatches", "retained_unique_clonotypes"
  ),
  value = as.character(c(
    length(header_fields), tcr_input_data_fields, nrow(tcr), uniqueN(tcr$Well_ID),
    all(tcr$record_id == tcr$Well_ID), nrow(tcr_retained),
    sum(!tcr$retained_in_step17), ncol(sce), 100 * mean(sce$tcr_detected),
    sum(tcr_retained$productive_primary), sum(amp_mismatch),
    sum(subject_mismatch), uniqueN(tcr_retained$clonotype_key, na.rm = TRUE)
  ))
)

fwrite(mapping_summary, file.path(results_dir, "18_tcr_mapping_summary.csv"))
fwrite(tcr_by_source, file.path(results_dir, "18_tcr_by_source.csv"))
fwrite(tcr_by_facs_gate, file.path(results_dir, "18_tcr_by_facs_gate.csv"))
fwrite(
  original_state_summary,
  file.path(results_dir, "18_tcr_original_state_summary.csv")
)
fwrite(clone_summary, file.path(results_dir, "18_tcr_clonotype_summary.tsv.gz"), sep = "\t")
fwrite(
  tcr_retained,
  file.path(results_dir, "18_retained_tcr_records.tsv.gz"),
  sep = "\t"
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

plot_18a <- ggplot(
  tcr_by_source,
  aes(x = source_group, y = percent_tcr_detected, fill = source_group)
) +
  geom_col(color = gray_outline, linewidth = 0.4, width = 0.65) +
  geom_text(
    aes(label = paste0(n_tcr_detected, " / ", n_cells)),
    vjust = -0.4, size = 4
  ) +
  scale_fill_manual(values = c("Tumor" = lavender, "PBMC" = gray_fill)) +
  coord_cartesian(ylim = c(0, max(tcr_by_source$percent_tcr_detected) * 1.15)) +
  labs(
    title = "TCR-beta recovery among retained cells",
    subtitle = "TCR recovery supports T-cell identity but is not required for classification",
    x = NULL, y = "Cells with a TCR-beta record (%)", fill = "Sample source"
  ) +
  theme_project

pca_scores <- SingleCellExperiment::reducedDim(sce, "PCA")
if (is.null(pca_scores) || ncol(pca_scores) < 2L) {
  stop("Step 17 PCA coordinates are absent from the SCE.")
}
pca_plot <- data.table(
  cell_id = colnames(sce),
  PC1 = pca_scores[, 1L],
  PC2 = pca_scores[, 2L],
  tcr_detected = sce$tcr_detected,
  original_state = sce$tcr_original_mc_grp
)

set.seed(123139)
background <- pca_plot[tcr_detected == FALSE]
if (nrow(background) > 30000L) background <- background[sample(.N, 30000L)]
tcr_points <- pca_plot[tcr_detected == TRUE]

pca_variance_file <- file.path("results", "17", "17_pca_variance.csv")
x_label <- "PC1"
y_label <- "PC2"
if (file.exists(pca_variance_file)) {
  pca_variance <- fread(pca_variance_file)
  if (nrow(pca_variance) >= 2L) {
    x_label <- paste0("PC1 (", round(pca_variance$percent_hvg_variance[1L], 1), "%)")
    y_label <- paste0("PC2 (", round(pca_variance$percent_hvg_variance[2L], 1), "%)")
  }
}

plot_18b <- ggplot() +
  geom_point(
    data = background, aes(x = PC1, y = PC2),
    color = gray_fill, size = 0.8, alpha = 0.45
  ) +
  geom_point(
    data = tcr_points, aes(x = PC1, y = PC2),
    color = lavender_dark, size = 1.0, alpha = 0.65
  ) +
  labs(
    title = "TCR-beta recovery in uncorrected PCA space",
    subtitle = "Gray: no retained TCR record; purple: mapped TCR record",
    x = x_label, y = y_label
  ) +
  theme_project + theme(legend.position = "none")

plot_18c <- ggplot(
  tcr_points[!is.na(original_state)],
  aes(x = PC1, y = PC2, fill = original_state)
) +
  geom_point(
    shape = 21, color = gray_outline, stroke = 0.15,
    size = 1.2, alpha = 0.65
  ) +
  labs(
    title = "Published T-cell states in uncorrected PCA space",
    subtitle = "mc_grp labels are retained for comparison, not used as new classifications",
    x = x_label, y = y_label, fill = "Published mc_grp"
  ) +
  theme_project

figure_paths <- file.path(
  figures_dir,
  c(
    "18A_tcr_recovery_by_source.png",
    "18B_tcr_recovery_in_pca.png",
    "18C_published_tcr_states_in_pca.png"
  )
)
ggsave(figure_paths[1L], plot_18a, width = 8, height = 6, dpi = 300, bg = "white")
ggsave(figure_paths[2L], plot_18b, width = 9, height = 7, dpi = 300, bg = "white")
ggsave(figure_paths[3L], plot_18c, width = 10, height = 7, dpi = 300, bg = "white")

fwrite(
  data.table(
    figure = paste0("18", LETTERS[1:3]),
    file = basename(figure_paths),
    purpose = c(
      "Quantify retained TCR-beta recovery by sample source",
      "Locate TCR-positive cells in uncorrected PCA space",
      "Compare published T-cell state labels in uncorrected PCA space"
    )
  ),
  file.path(results_dir, "18_figure_manifest.csv")
)

# ------------------------------------------------------------------------------
# Save the augmented object and provenance
# ------------------------------------------------------------------------------

sce_metadata <- S4Vectors::metadata(sce)
sce_metadata$step <- 18L
sce_metadata$tcr_integration <- list(
  source_file = tcr_file,
  mapping_key = paste(
    "Well_ID to SCE column name, validated by Amp.Batch and a subject prefix",
    "derived from the TCR Patient label"
  ),
  unnamed_first_column = "record_id; required to equal Well_ID",
  patient_field_interpretation = paste(
    "The TCR Patient field can contain specimen-style labels such as p12-2;",
    "the ^p[0-9]+ prefix is compared with SCE subject_id"
  ),
  original_mc_grp_usage = "external comparison only",
  tcr_required_for_t_cell_classification = FALSE
)
S4Vectors::metadata(sce) <- sce_metadata

saveRDS(sce, output_sce_file, compress = "gzip")
writeLines(capture.output(sessionInfo()), file.path(results_dir, "18_sessionInfo.txt"))

dataset_qc <- data.table(
  metric = c(
    "genes", "cells", "tumor_cells", "pbmc_cells", "subjects",
    "tcr_records_total", "tcr_records_retained", "tcr_positive_cells",
    "primary_productive_tcr_cells", "unique_retained_clonotypes",
    "amp_batch_mapping_validation", "derived_subject_mapping_validation",
    "pca_coordinates_preserved", "output_sce_size_mb"
  ),
  value = as.character(c(
    nrow(sce), ncol(sce), sum(sce$source_group == "Tumor"),
    sum(sce$source_group == "PBMC"), uniqueN(sce$subject_id),
    nrow(tcr), nrow(tcr_retained), sum(sce$tcr_detected),
    sum(sce$tcr_productive_primary),
    uniqueN(sce$tcr_clonotype_key, na.rm = TRUE),
    ifelse(any(amp_mismatch), "failed", "passed"),
    ifelse(any(subject_mismatch), "failed", "passed"),
    identical(pca_scores, SingleCellExperiment::reducedDim(sce, "PCA")),
    as.numeric(object.size(sce)) / 1024^2
  ))
)
fwrite(dataset_qc, file.path(results_dir, "18_dataset_qc.csv"))

message("Step 18 complete.")
message("TCR records retained: ", nrow(tcr_retained), " / ", nrow(tcr))
message("TCR-positive fraction of Step 17 cells: ", round(100 * mean(sce$tcr_detected), 2), "%")
message("Augmented SCE: ", output_sce_file)
print(mapping_summary)
print(tcr_by_source)
