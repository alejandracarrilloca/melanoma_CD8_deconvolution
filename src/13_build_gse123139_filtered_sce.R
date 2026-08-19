#!/usr/bin/env Rscript

# ==============================================================================
# Build the filtered and gene-harmonized GSE123139 SingleCellExperiment
# ==============================================================================
#
# This script finalizes cell-level QC after inspection in Step 12, harmonizes
# the two gene configurations found across plates by retaining their
# intersection, builds one sparse raw-count matrix, and saves a reusable
# SingleCellExperiment for downstream CD8 T-cell annotation.
#
# Final QC definition:
#   - at least 1,000 total UMIs
#   - at least 200 detected genes
#   - at most 20% mitochondrial UMIs
#   - not a plate-specific high-library outlier (> median + 4 MAD)
#
# Inputs:
#   data/scRNA/GSE123139/raw_counts/GSM*_*.txt.gz
#   data/analysis_ready/scRNA/GSE123139/12_GSE123139_cell_qc.tsv.gz
#
# Reusable outputs:
#   data/analysis_ready/scRNA/GSE123139/
#     13_GSE123139_filtered_sce.rds
#     13_GSE123139_filtered_cell_metadata.tsv.gz
#     13_GSE123139_common_gene_symbols.tsv.gz
#
# Run-specific outputs:
#   results/13/
#     13_final_qc_configuration.csv
#     13_gene_harmonization_qc.csv
#     13_plate_matrix_build_qc.csv
#     13_cell_retention_by_source.csv
#     13_cell_retention_by_subject.csv
#     13_dataset_qc.csv
#     13_sessionInfo.txt
#
# Run from the project root with:
#   Rscript src/13_build_gse123139_filtered_sce.R
# ==============================================================================

options(stringsAsFactors = FALSE)

required_packages <- c(
  "data.table",
  "Matrix",
  "SingleCellExperiment",
  "S4Vectors"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages(library(data.table))

# ------------------------------------------------------------------------------
# Paths and final QC thresholds
# ------------------------------------------------------------------------------

project_root <- getwd()

raw_dir <- file.path(
  project_root,
  "data", "scRNA", "GSE123139", "raw_counts"
)

analysis_dir <- file.path(
  project_root,
  "data", "analysis_ready", "scRNA", "GSE123139"
)

results_dir <- file.path(project_root, "results", "13")

dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

cell_qc_file <- file.path(
  analysis_dir,
  "12_GSE123139_cell_qc.tsv.gz"
)

sce_file <- file.path(
  analysis_dir,
  "13_GSE123139_filtered_sce.rds"
)

metadata_file <- file.path(
  analysis_dir,
  "13_GSE123139_filtered_cell_metadata.tsv.gz"
)

gene_file <- file.path(
  analysis_dir,
  "13_GSE123139_common_gene_symbols.tsv.gz"
)

minimum_umi <- 1000L
minimum_genes <- 200L
maximum_mitochondrial_percent <- 20
exclude_high_library_outliers <- TRUE

configuration <- data.table(
  parameter = c(
    "minimum_umi",
    "minimum_detected_genes",
    "maximum_mitochondrial_percent",
    "exclude_plate_high_library_outliers",
    "gene_harmonization",
    "saved_assay"
  ),
  value = c(
    as.character(minimum_umi),
    as.character(minimum_genes),
    as.character(maximum_mitochondrial_percent),
    as.character(exclude_high_library_outliers),
    "intersection_across_all_plates_in_first_plate_order",
    "raw_UMI_counts"
  )
)

fwrite(
  configuration,
  file.path(results_dir, "13_final_qc_configuration.csv")
)

# ------------------------------------------------------------------------------
# Validate inputs and define the retained cells
# ------------------------------------------------------------------------------

if (!file.exists(cell_qc_file)) {
  stop("Missing Step 12 cell-QC table: ", cell_qc_file)
}

count_files <- sort(list.files(
  raw_dir,
  pattern = "^GSM[0-9]+_.*\\.txt\\.gz$",
  full.names = TRUE
))

if (length(count_files) == 0L) {
  stop("No extracted GSE123139 count files found in: ", raw_dir)
}

cell_qc <- fread(cell_qc_file, check.names = FALSE)

required_qc_columns <- c(
  "cell_id",
  "count_file",
  "source_group",
  "subject_id",
  "specimen_id",
  "n_counts",
  "n_genes",
  "mitochondrial_percent",
  "high_library_outlier"
)

missing_qc_columns <- setdiff(required_qc_columns, names(cell_qc))

if (length(missing_qc_columns) > 0L) {
  stop(
    "Step 12 table is missing required column(s): ",
    paste(missing_qc_columns, collapse = ", ")
  )
}

if (anyDuplicated(cell_qc$cell_id)) {
  stop("Duplicated cell IDs are present in the Step 12 QC table.")
}

cell_qc[, final_qc_pass :=
  n_counts >= minimum_umi &
  n_genes >= minimum_genes &
  mitochondrial_percent <= maximum_mitochondrial_percent &
  !high_library_outlier
]

cell_qc[, primary_tumor_analysis :=
  final_qc_pass & source_group == "Tumor"
]

retained_metadata <- copy(cell_qc[final_qc_pass == TRUE])

if (nrow(retained_metadata) == 0L) {
  stop("No cells passed the final QC definition.")
}

retention_by_source <- cell_qc[, .(
  n_cell_columns = .N,
  n_final_qc_pass = sum(final_qc_pass),
  percent_final_qc_pass = 100 * mean(final_qc_pass)
), by = source_group][order(source_group)]

retention_by_subject <- cell_qc[, .(
  n_cell_columns = .N,
  n_final_qc_pass = sum(final_qc_pass),
  percent_final_qc_pass = 100 * mean(final_qc_pass),
  n_tumor_final_qc_pass = sum(primary_tumor_analysis),
  n_pbmc_final_qc_pass = sum(final_qc_pass & source_group == "PBMC")
), by = subject_id][order(subject_id)]

fwrite(
  retention_by_source,
  file.path(results_dir, "13_cell_retention_by_source.csv")
)

fwrite(
  retention_by_subject,
  file.path(results_dir, "13_cell_retention_by_subject.csv")
)

# ------------------------------------------------------------------------------
# Read the cell IDs stored in a plate-file first line
# ------------------------------------------------------------------------------

read_cell_header <- function(path) {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection))

  line <- readLines(connection, n = 1L, warn = FALSE)

  if (length(line) != 1L) {
    stop("Could not read the cell-ID header from: ", basename(path))
  }

  # These plate files have one more field in every data row than in the first
  # line: the first-line fields are all cell IDs, while the unlabelled first
  # data column contains gene symbols. Therefore, do not discard the first
  # header field (doing so shifts every selected cell column by one).
  fields <- strsplit(line, "\t", fixed = TRUE)[[1L]]

  if (any(!nzchar(fields)) || anyDuplicated(fields)) {
    stop("Invalid or duplicated cell IDs in the header of: ", basename(path))
  }

  fields
}

read_gene_column <- function(path) {
  as.character(fread(
    path,
    skip = 1L,
    header = FALSE,
    select = 1L,
    showProgress = FALSE
  )[[1L]])
}

# ------------------------------------------------------------------------------
# Determine a common, consistently ordered gene set across all plates
# ------------------------------------------------------------------------------

message("Determining the gene intersection across ", length(count_files), " plates")

reference_genes <- read_gene_column(count_files[[1L]])

if (anyDuplicated(reference_genes)) {
  stop("Duplicated gene symbols occur in the reference plate.")
}

common_genes <- reference_genes
union_genes <- reference_genes
gene_qc_list <- vector("list", length(count_files))

for (i in seq_along(count_files)) {
  path <- count_files[[i]]
  genes <- read_gene_column(path)

  if (anyDuplicated(genes)) {
    stop("Duplicated gene symbols occur in: ", basename(path))
  }

  identical_order <- identical(genes, reference_genes)
  identical_set <-
    length(genes) == length(reference_genes) &&
    all(reference_genes %chin% genes)

  gene_qc_list[[i]] <- data.table(
    count_file = basename(path),
    n_gene_rows = length(genes),
    identical_to_reference_order = identical_order,
    identical_to_reference_gene_set = identical_set
  )

  common_genes <- common_genes[common_genes %chin% genes]
  union_genes <- union(union_genes, genes)
}

gene_qc <- rbindlist(gene_qc_list)
gene_qc[, common_gene_count := length(common_genes)]
gene_qc[, union_gene_count := length(union_genes)]

if (length(common_genes) == 0L) {
  stop("The intersection of genes across plates is empty.")
}

fwrite(
  gene_qc,
  file.path(results_dir, "13_gene_harmonization_qc.csv")
)

fwrite(
  data.table(gene_symbol = common_genes),
  gene_file,
  sep = "\t"
)

rm(union_genes)
invisible(gc())

# ------------------------------------------------------------------------------
# Build one sparse filtered count matrix plate by plate
# ------------------------------------------------------------------------------

message(
  "Building sparse matrix for ",
  format(nrow(retained_metadata), big.mark = ","),
  " retained cells and ",
  format(length(common_genes), big.mark = ","),
  " common genes"
)

sparse_plate_matrices <- list()
plate_build_qc <- list()
matrix_counter <- 0L

for (i in seq_along(count_files)) {
  path <- count_files[[i]]
  file_name <- basename(path)
  message("Processing plate ", i, "/", length(count_files), ": ", file_name)

  cell_ids <- read_cell_header(path)
  plate_metadata <- retained_metadata[count_file == file_name]

  if (nrow(plate_metadata) == 0L) {
    plate_build_qc[[i]] <- data.table(
      count_file = file_name,
      n_cell_columns = length(cell_ids),
      n_retained_cells = 0L,
      n_common_genes = length(common_genes),
      harmonized_count_fraction_minimum = NA_real_,
      harmonized_count_fraction_median = NA_real_,
      harmonized_count_fraction_maximum = NA_real_
    )
    next
  }

  retained_column_indices <- match(plate_metadata$cell_id, cell_ids)

  if (anyNA(retained_column_indices)) {
    missing_cells <- plate_metadata$cell_id[is.na(retained_column_indices)]
    stop(
      "Retained cell(s) are missing from ", file_name, ": ",
      paste(head(missing_cells, 5L), collapse = ", ")
    )
  }

  selected_indices <- c(1L, retained_column_indices + 1L)
  selected_names <- c("gene_symbol", plate_metadata$cell_id)

  plate_table <- fread(
    path,
    skip = 1L,
    header = FALSE,
    select = selected_indices,
    col.names = selected_names,
    check.names = FALSE,
    showProgress = FALSE
  )

  plate_genes <- as.character(plate_table[["gene_symbol"]])
  gene_indices <- match(common_genes, plate_genes)

  if (anyNA(gene_indices)) {
    stop("A common gene is unexpectedly absent from: ", file_name)
  }

  dense_counts <- as.matrix(
    plate_table[gene_indices, -1L, with = FALSE]
  )
  storage.mode(dense_counts) <- "numeric"
  rownames(dense_counts) <- common_genes
  colnames(dense_counts) <- plate_metadata$cell_id

  sparse_counts <- Matrix::Matrix(dense_counts, sparse = TRUE)
  harmonized_counts <- as.numeric(Matrix::colSums(sparse_counts))
  original_counts <- plate_metadata$n_counts
  retained_fraction <- harmonized_counts / original_counts

  if (any(harmonized_counts > original_counts + 1e-8, na.rm = TRUE)) {
    stop("Harmonized counts exceed original QC counts in: ", file_name)
  }

  matrix_counter <- matrix_counter + 1L
  sparse_plate_matrices[[matrix_counter]] <- sparse_counts

  plate_build_qc[[i]] <- data.table(
    count_file = file_name,
    n_cell_columns = length(cell_ids),
    n_retained_cells = ncol(sparse_counts),
    n_common_genes = nrow(sparse_counts),
    harmonized_count_fraction_minimum = min(retained_fraction, na.rm = TRUE),
    harmonized_count_fraction_median = median(retained_fraction, na.rm = TRUE),
    harmonized_count_fraction_maximum = max(retained_fraction, na.rm = TRUE)
  )

  rm(
    plate_table,
    plate_genes,
    dense_counts,
    sparse_counts,
    harmonized_counts
  )
  invisible(gc())
}

plate_build_qc <- rbindlist(plate_build_qc, fill = TRUE)

fwrite(
  plate_build_qc,
  file.path(results_dir, "13_plate_matrix_build_qc.csv")
)

if (length(sparse_plate_matrices) == 0L) {
  stop("No plate matrices were constructed.")
}

counts <- do.call(cbind, sparse_plate_matrices)

if (anyDuplicated(colnames(counts))) {
  stop("Duplicated cell IDs occur in the combined count matrix.")
}

metadata_order <- match(colnames(counts), retained_metadata$cell_id)

if (anyNA(metadata_order)) {
  stop("Some matrix cells are absent from the retained metadata.")
}

retained_metadata <- retained_metadata[metadata_order]

if (!identical(colnames(counts), retained_metadata$cell_id)) {
  stop("Count-matrix columns and metadata rows could not be aligned.")
}

retained_metadata[, harmonized_n_counts :=
  as.numeric(Matrix::colSums(counts))
]

retained_metadata[, harmonized_count_fraction :=
  harmonized_n_counts / n_counts
]

# ------------------------------------------------------------------------------
# Construct and save the SingleCellExperiment
# ------------------------------------------------------------------------------

col_data <- S4Vectors::DataFrame(
  as.data.frame(retained_metadata),
  row.names = retained_metadata$cell_id
)

row_data <- S4Vectors::DataFrame(
  gene_symbol = common_genes,
  row.names = common_genes
)

sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = counts),
  rowData = row_data,
  colData = col_data,
  metadata = list(
    accession = "GSE123139",
    step = 13L,
    qc_definition = list(
      minimum_umi = minimum_umi,
      minimum_detected_genes = minimum_genes,
      maximum_mitochondrial_percent = maximum_mitochondrial_percent,
      excluded_plate_high_library_outliers =
        exclude_high_library_outliers
    ),
    gene_harmonization =
      "intersection across all plates, ordered as in the first plate",
    assay_scale = "raw UMI counts"
  )
)

saveRDS(sce, sce_file, compress = "gzip")
fwrite(retained_metadata, metadata_file, sep = "\t")

# ------------------------------------------------------------------------------
# Final validation summary
# ------------------------------------------------------------------------------

dataset_qc <- data.table(
  metric = c(
    "input_plates",
    "input_cell_columns",
    "retained_cells",
    "retained_tumor_cells",
    "retained_pbmc_cells",
    "retained_subjects",
    "retained_specimens",
    "reference_gene_count",
    "common_gene_count",
    "plates_with_reference_gene_order",
    "plates_with_alternative_gene_order_or_set",
    "matrix_nonzero_entries",
    "harmonized_count_fraction_minimum",
    "harmonized_count_fraction_median",
    "harmonized_count_fraction_maximum",
    "sce_object_size_mb"
  ),
  value = as.character(c(
    length(count_files),
    nrow(cell_qc),
    ncol(sce),
    sum(sce$source_group == "Tumor"),
    sum(sce$source_group == "PBMC"),
    uniqueN(sce$subject_id),
    uniqueN(sce$specimen_id),
    length(reference_genes),
    nrow(sce),
    sum(gene_qc$identical_to_reference_order),
    sum(!gene_qc$identical_to_reference_order),
    Matrix::nnzero(counts),
    min(sce$harmonized_count_fraction, na.rm = TRUE),
    median(sce$harmonized_count_fraction, na.rm = TRUE),
    max(sce$harmonized_count_fraction, na.rm = TRUE),
    as.numeric(object.size(sce)) / 1024^2
  ))
)

fwrite(dataset_qc, file.path(results_dir, "13_dataset_qc.csv"))

writeLines(
  capture.output(sessionInfo()),
  file.path(results_dir, "13_sessionInfo.txt")
)

message("Step 13 complete.")
message("SingleCellExperiment: ", sce_file)
message("Dimensions: ", nrow(sce), " genes x ", ncol(sce), " cells")
message("Retained cells by source:")
print(retention_by_source)
message("Dataset QC:")
print(dataset_qc)
