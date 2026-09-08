#!/usr/bin/env Rscript

# ==============================================================================
# Step 16: Build the final GSE123139 object after ERCC-aware QC
# ==============================================================================
#
# Step 14 established that 92 plate-specific features are ERCC spike-ins and
# that QC must be based on the 55,765 common endogenous features. This script
# creates the corrected SingleCellExperiment without overwriting the Step 13
# audit object.
#
# To avoid rebuilding all retained columns from raw files, the script:
#   1. Subsets the Step 13 matrix to cells that still pass endogenous QC.
#   2. Reads only final-QC cells that are absent from the Step 13 object.
#   3. Combines and reorders all cells to the Step 14 metadata.
#   4. Verifies every column sum and detected-gene count against Step 14.
#
# Inputs:
#   data/analysis_ready/scRNA/GSE123139/
#     13_GSE123139_filtered_sce.rds
#     13_GSE123139_common_gene_symbols.tsv.gz
#     14_GSE123139_endogenous_cell_qc.tsv.gz
#   data/scRNA/GSE123139/raw_counts/GSM*_*.txt.gz
#
# Outputs:
#   data/analysis_ready/scRNA/GSE123139/
#     16_GSE123139_endogenous_qc_sce.rds
#     16_GSE123139_endogenous_qc_cell_metadata.tsv.gz
#   results/16/
#     16_qc_transition_summary.csv
#     16_cell_retention_by_source.csv
#     16_added_cells.tsv.gz
#     16_removed_cells.tsv.gz
#     16_dataset_qc.csv
#     16_sessionInfo.txt
#
# Run from the project root:
#   Rscript src/16_build_gse123139_endogenous_qc_sce.R
# ==============================================================================

options(stringsAsFactors = FALSE)

required_packages <- c(
  "data.table", "Matrix", "SingleCellExperiment", "S4Vectors"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages(library(data.table))

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

project_root <- normalizePath(getwd(), mustWork = TRUE)
raw_dir <- file.path(project_root, "data", "scRNA", "GSE123139", "raw_counts")
analysis_dir <- file.path(
  project_root, "data", "analysis_ready", "scRNA", "GSE123139"
)
results_dir <- file.path(project_root, "results", "16")

dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

step13_sce_file <- file.path(
  analysis_dir, "13_GSE123139_filtered_sce.rds"
)
common_gene_file <- file.path(
  analysis_dir, "13_GSE123139_common_gene_symbols.tsv.gz"
)
corrected_qc_file <- file.path(
  analysis_dir, "14_GSE123139_endogenous_cell_qc.tsv.gz"
)

output_sce_file <- file.path(
  analysis_dir, "16_GSE123139_endogenous_qc_sce.rds"
)
output_metadata_file <- file.path(
  analysis_dir, "16_GSE123139_endogenous_qc_cell_metadata.tsv.gz"
)

minimum_endogenous_umi <- 800L
minimum_endogenous_genes <- 200L
maximum_endogenous_mitochondrial_percent <- 20

for (path in c(step13_sce_file, common_gene_file, corrected_qc_file)) {
  if (!file.exists(path)) stop("Missing required input: ", path)
}

# ------------------------------------------------------------------------------
# Load and validate the audited inputs
# ------------------------------------------------------------------------------

message("Loading the Step 13 SingleCellExperiment")
sce13 <- readRDS(step13_sce_file)
qc <- fread(corrected_qc_file, check.names = FALSE)
common_genes <- as.character(fread(common_gene_file)[[1L]])

required_qc_columns <- c(
  "cell_id", "count_file", "source_group", "subject_id", "specimen_id",
  "original_final_qc_pass",
  "endogenous_n_counts", "endogenous_n_genes",
  "endogenous_mitochondrial_percent", "endogenous_high_library_outlier"
)
missing_qc_columns <- setdiff(required_qc_columns, names(qc))

if (length(missing_qc_columns) > 0L) {
  stop(
    "Step 14 QC table is missing: ",
    paste(missing_qc_columns, collapse = ", ")
  )
}

if (anyDuplicated(qc$cell_id)) stop("Duplicated cell IDs occur in Step 14 QC.")
if (anyDuplicated(common_genes)) stop("Duplicated common genes are present.")
if (!identical(rownames(sce13), common_genes)) {
  stop("Step 13 matrix rows do not match the saved common-gene order.")
}
if (anyDuplicated(colnames(sce13))) {
  stop("Duplicated cell IDs occur in the Step 13 matrix.")
}

expected_step13_ids <- qc[original_final_qc_pass == TRUE, cell_id]
if (!setequal(colnames(sce13), expected_step13_ids)) {
  stop("Step 13 cells do not equal the original final-QC set in Step 14.")
}

qc[, final_qc_pass_800 :=
  endogenous_n_counts >= minimum_endogenous_umi &
  endogenous_n_genes >= minimum_endogenous_genes &
  !is.na(endogenous_mitochondrial_percent) &
  endogenous_mitochondrial_percent <=
    maximum_endogenous_mitochondrial_percent &
  !endogenous_high_library_outlier
]

final_metadata <- copy(qc[final_qc_pass_800 == TRUE])
final_ids <- final_metadata$cell_id
existing_ids <- final_ids[final_ids %chin% colnames(sce13)]
gained_ids <- final_ids[!final_ids %chin% colnames(sce13)]
removed_ids <- colnames(sce13)[!colnames(sce13) %chin% final_ids]

transition_summary <- data.table(
  transition = c(
    "Step 13 cells", "Retained from Step 13", "Removed after ERCC correction",
    "Added under final endogenous QC", "Final Step 16 cells"
  ),
  n_cells = c(
    ncol(sce13), length(existing_ids), length(removed_ids),
    length(gained_ids), length(final_ids)
  )
)

fwrite(
  transition_summary,
  file.path(results_dir, "16_qc_transition_summary.csv")
)
fwrite(
  qc[cell_id %chin% gained_ids],
  file.path(results_dir, "16_added_cells.tsv.gz"),
  sep = "\t"
)
fwrite(
  qc[cell_id %chin% removed_ids],
  file.path(results_dir, "16_removed_cells.tsv.gz"),
  sep = "\t"
)

message("Retaining ", length(existing_ids), " cells from Step 13")
message("Reading ", length(gained_ids), " recovered cells from raw plates")

# ------------------------------------------------------------------------------
# Retain existing columns, then release the full Step 13 object
# ------------------------------------------------------------------------------

counts_existing <- SummarizedExperiment::assay(sce13, "counts")[
  , existing_ids, drop = FALSE
]
row_data <- SummarizedExperiment::rowData(sce13)

rm(sce13)
invisible(gc())

# ------------------------------------------------------------------------------
# Read only cells absent from Step 13 but recovered by corrected QC
# ------------------------------------------------------------------------------

read_cell_header <- function(path) {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection))
  line <- readLines(connection, n = 1L, warn = FALSE)

  if (length(line) != 1L || !nzchar(line)) {
    stop("Could not read the cell header from: ", basename(path))
  }

  fields <- strsplit(line, "\t", fixed = TRUE)[[1L]]
  if (any(!nzchar(fields)) || anyDuplicated(fields)) {
    stop("Invalid cell IDs in: ", basename(path))
  }
  fields
}

gained_matrix_list <- list()
gained_counter <- 0L

if (length(gained_ids) > 0L) {
  gained_metadata <- final_metadata[cell_id %chin% gained_ids]
  gained_files <- unique(gained_metadata$count_file)

  for (i in seq_along(gained_files)) {
    file_name <- gained_files[[i]]
    path <- file.path(raw_dir, file_name)
    message(
      "Reading recovered cells from plate ", i, "/", length(gained_files),
      ": ", file_name
    )

    if (!file.exists(path)) stop("Missing raw count file: ", path)

    plate_metadata <- gained_metadata[count_file == file_name]
    cell_ids <- read_cell_header(path)
    retained_column_indices <- match(plate_metadata$cell_id, cell_ids)

    if (anyNA(retained_column_indices)) {
      stop("A recovered cell is absent from its raw plate: ", file_name)
    }

    plate_table <- fread(
      path,
      skip = 1L,
      header = FALSE,
      select = c(1L, retained_column_indices + 1L),
      col.names = c("gene_symbol", plate_metadata$cell_id),
      check.names = FALSE,
      showProgress = FALSE
    )

    plate_genes <- as.character(plate_table[["gene_symbol"]])
    gene_indices <- match(common_genes, plate_genes)

    if (anyNA(gene_indices)) {
      stop("A common gene is absent from: ", file_name)
    }

    dense_counts <- as.matrix(
      plate_table[gene_indices, -1L, with = FALSE]
    )
    storage.mode(dense_counts) <- "numeric"
    rownames(dense_counts) <- common_genes
    colnames(dense_counts) <- plate_metadata$cell_id

    gained_counter <- gained_counter + 1L
    gained_matrix_list[[gained_counter]] <- Matrix::Matrix(
      dense_counts, sparse = TRUE
    )

    rm(plate_table, dense_counts)
    invisible(gc())
  }

  counts_gained <- do.call(cbind, gained_matrix_list)
} else {
  counts_gained <- Matrix::Matrix(
    0,
    nrow = length(common_genes),
    ncol = 0L,
    sparse = TRUE,
    dimnames = list(common_genes, character())
  )
}

if (!setequal(colnames(counts_gained), gained_ids)) {
  stop("Recovered matrix columns do not equal the gained-cell set.")
}

# ------------------------------------------------------------------------------
# Combine, order, and validate against the endogenous QC metrics
# ------------------------------------------------------------------------------

counts_final <- cbind(counts_existing, counts_gained)
final_order <- match(final_ids, colnames(counts_final))

if (anyNA(final_order)) stop("A final cell is absent from the combined matrix.")
counts_final <- counts_final[, final_order, drop = FALSE]

if (!identical(colnames(counts_final), final_ids)) {
  stop("Final matrix columns and corrected metadata could not be aligned.")
}
if (anyDuplicated(colnames(counts_final))) {
  stop("Duplicated cell IDs occur in the final matrix.")
}
if (any(grepl("^ERCC-", rownames(counts_final), ignore.case = TRUE))) {
  stop("ERCC rows unexpectedly occur in the final matrix.")
}

matrix_n_counts <- as.numeric(Matrix::colSums(counts_final))
matrix_n_genes <- as.integer(Matrix::colSums(counts_final > 0))

if (!isTRUE(all.equal(
  matrix_n_counts,
  as.numeric(final_metadata$endogenous_n_counts),
  tolerance = 0,
  check.attributes = FALSE
))) {
  stop("Final matrix column sums do not equal Step 14 endogenous counts.")
}

if (!identical(matrix_n_genes, as.integer(final_metadata$endogenous_n_genes))) {
  stop("Final matrix detected-gene counts do not equal Step 14 QC.")
}

final_metadata[, final_object := TRUE]
final_metadata[, primary_tumor_analysis := source_group == "Tumor"]

col_data <- S4Vectors::DataFrame(
  as.data.frame(final_metadata),
  row.names = final_metadata$cell_id
)

sce16 <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = counts_final),
  rowData = row_data,
  colData = col_data,
  metadata = list(
    accession = "GSE123139",
    step = 16L,
    parent_object = basename(step13_sce_file),
    qc_definition = list(
      minimum_endogenous_umi = minimum_endogenous_umi,
      minimum_endogenous_detected_genes = minimum_endogenous_genes,
      maximum_endogenous_mitochondrial_percent =
        maximum_endogenous_mitochondrial_percent,
      excluded_plate_endogenous_high_library_outliers = TRUE,
      high_library_nmads = 4
    ),
    excluded_features = list(
      n_features = 92L,
      feature_type = "ERCC spike-ins",
      evidence = "Step 14 feature-name audit"
    ),
    assay_scale = "raw endogenous/common-gene UMI counts"
  )
)

saveRDS(sce16, output_sce_file, compress = "gzip")
fwrite(final_metadata, output_metadata_file, sep = "\t")

# ------------------------------------------------------------------------------
# Final summaries
# ------------------------------------------------------------------------------

retention_by_source <- qc[, .(
  n_input_wells = .N,
  n_final_cells = sum(final_qc_pass_800),
  percent_final = 100 * mean(final_qc_pass_800)
), by = source_group][order(source_group)]

dataset_qc <- data.table(
  metric = c(
    "genes", "final_cells", "final_tumor_cells", "final_pbmc_cells",
    "retained_from_step13", "removed_from_step13", "recovered_cells_added",
    "retained_subjects", "retained_specimens", "matrix_nonzero_entries",
    "ercc_rows_in_final_matrix", "column_sum_validation",
    "detected_gene_validation", "sce_object_size_mb"
  ),
  value = as.character(c(
    nrow(sce16), ncol(sce16), sum(sce16$source_group == "Tumor"),
    sum(sce16$source_group == "PBMC"), length(existing_ids),
    length(removed_ids), length(gained_ids), uniqueN(sce16$subject_id),
    uniqueN(sce16$specimen_id), Matrix::nnzero(counts_final),
    sum(grepl("^ERCC-", rownames(sce16), ignore.case = TRUE)),
    "passed", "passed", as.numeric(object.size(sce16)) / 1024^2
  ))
)

fwrite(
  retention_by_source,
  file.path(results_dir, "16_cell_retention_by_source.csv")
)
fwrite(dataset_qc, file.path(results_dir, "16_dataset_qc.csv"))
writeLines(capture.output(sessionInfo()), file.path(results_dir, "16_sessionInfo.txt"))

message("Step 16 complete.")
message("Final SingleCellExperiment: ", output_sce_file)
message("Dimensions: ", nrow(sce16), " genes x ", ncol(sce16), " cells")
message("QC transition summary:")
print(transition_summary)
message("Final retention by source:")
print(retention_by_source)
message("Dataset QC:")
print(dataset_qc)