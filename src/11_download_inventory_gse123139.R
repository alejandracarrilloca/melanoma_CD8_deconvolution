#!/usr/bin/env Rscript

# ==============================================================================
# Step 11: Download and inventory GSE123139 melanoma single-cell RNA-seq data
# ==============================================================================
#
# Purpose:
#   Download the processed MARS-seq UMI-count files and GEO sample metadata,
#   verify the 204 plate-level files, and create a cell-level metadata table for
#   the downstream CD8 T-cell analysis.
#
# Inputs downloaded by this script:
#   data/scRNA/GSE123139/GSE123139_RAW.tar
#   data/scRNA/GSE123139/GSE123139_T_cells_tcrb_v2.txt.gz
#   GEO metadata for GSE123139 (through GEOquery)
#
# Outputs:
#   results/11/11_download_manifest.csv
#   results/11/11_geo_sample_metadata.csv
#   results/11/11_plate_file_inventory.csv
#   results/11/11_dataset_qc.csv
#   results/11/11_marker_availability.csv
#   results/11/11_sessionInfo.txt
#   data/analysis_ready/scRNA/GSE123139/
#     11_GSE123139_cell_metadata.csv.gz
#
# This step does not filter, normalize, cluster, or annotate CD8 states. Those
# operations should follow only after the metadata and count-file joins pass QC.
#
# Run from the project root:
#   Rscript src/11_download_inventory_gse123139.R
# ==============================================================================

options(stringsAsFactors = FALSE, timeout = max(3600, getOption("timeout")))

required_packages <- c("data.table", "GEOquery")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(GEOquery)
})

# ------------------------------------------------------------------------------
# Paths and URLs
# ------------------------------------------------------------------------------

project_root <- normalizePath(getwd(), mustWork = TRUE)

input_dir <- file.path(project_root, "data", "scRNA", "GSE123139")
raw_counts_dir <- file.path(input_dir, "raw_counts")
analysis_ready_dir <- file.path(
  project_root, "data", "analysis_ready", "scRNA", "GSE123139"
)
results_dir <- file.path(project_root, "results", "11")

dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(raw_counts_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(analysis_ready_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

raw_tar <- file.path(input_dir, "GSE123139_RAW.tar")
tcr_file <- file.path(input_dir, "GSE123139_T_cells_tcrb_v2.txt.gz")

raw_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/",
  "GSE123139/suppl/GSE123139_RAW.tar"
)

tcr_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/",
  "GSE123139/suppl/GSE123139_T_cells_tcrb_v2.txt.gz"
)

download_if_missing <- function(url, destination) {
  if (file.exists(destination) && file.info(destination)$size > 0) {
    message("Using existing file: ", destination)
    return("existing")
  }

  message("Downloading: ", basename(destination))
  temporary <- paste0(destination, ".partial")

  if (file.exists(temporary)) {
    unlink(temporary)
  }

  utils::download.file(
    url = url,
    destfile = temporary,
    mode = "wb",
    quiet = FALSE
  )

  if (!file.exists(temporary) || file.info(temporary)$size == 0) {
    stop("Download failed or produced an empty file: ", destination)
  }

  if (!file.rename(temporary, destination)) {
    stop("Could not move completed download to: ", destination)
  }

  "downloaded"
}

raw_download_status <- download_if_missing(raw_url, raw_tar)
tcr_download_status <- download_if_missing(tcr_url, tcr_file)

# ------------------------------------------------------------------------------
# Extract plate-level count files
# ------------------------------------------------------------------------------

existing_plate_files <- list.files(
  raw_counts_dir,
  pattern = "^GSM[0-9]+_.+\\.txt(\\.gz)?$",
  full.names = TRUE,
  recursive = TRUE
)

if (length(existing_plate_files) == 0L) {
  message("Extracting processed UMI-count files")
  utils::untar(raw_tar, exdir = raw_counts_dir)
}

plate_files <- sort(list.files(
  raw_counts_dir,
  pattern = "^GSM[0-9]+_.+\\.txt(\\.gz)?$",
  full.names = TRUE,
  recursive = TRUE
))

if (length(plate_files) == 0L) {
  stop("No plate-level count files were found after extraction.")
}

# ------------------------------------------------------------------------------
# Retrieve and standardize GEO sample metadata
# ------------------------------------------------------------------------------

geo_metadata_file <- file.path(
  results_dir,
  "11_geo_sample_metadata.csv"
)

collapse_meta <- function(meta, field) {
  values <- unlist(meta[names(meta) == field], use.names = FALSE)
  values <- trimws(as.character(values))
  values <- values[nzchar(values)]
  if (length(values) == 0L) NA_character_ else paste(unique(values), collapse = " | ")
}

all_characteristics <- function(meta) {
  characteristic_names <- grep(
    "^characteristics_ch1",
    names(meta),
    value = TRUE
  )
  values <- unlist(meta[characteristic_names], use.names = FALSE)
  trimws(as.character(values))
}

extract_characteristic <- function(meta, label) {
  values <- all_characteristics(meta)
  pattern <- paste0("^", label, "\\s*:\\s*")
  matched <- grep(pattern, values, ignore.case = TRUE, value = TRUE)

  if (length(matched) == 0L) {
    return(NA_character_)
  }

  result <- sub(pattern, "", matched, ignore.case = TRUE)
  paste(unique(trimws(result)), collapse = " | ")
}

if (file.exists(geo_metadata_file)) {
  message("Using existing GEO metadata: ", geo_metadata_file)
  geo_metadata <- fread(geo_metadata_file)
} else {
  message("Retrieving GEO metadata for GSE123139")
  gse <- GEOquery::getGEO(
    "GSE123139",
    GSEMatrix = FALSE,
    getGPL = FALSE
  )

  gsm_list <- GEOquery::GSMList(gse)

  geo_metadata <- rbindlist(lapply(
    names(gsm_list),
    function(current_gsm_accession) {
      meta <- GEOquery::Meta(gsm_list[[current_gsm_accession]])

      data.table(
        gsm_accession = current_gsm_accession,
        title = collapse_meta(meta, "title"),
        source_name = collapse_meta(meta, "source_name_ch1"),
        plate_id = extract_characteristic(meta, "plate id"),
        amplification_batch = extract_characteristic(
          meta,
          "amplification batch"
        ),
        sample_source = extract_characteristic(meta, "sample source"),
        patient_id = extract_characteristic(meta, "patient id"),
        facs_gate = extract_characteristic(meta, "facs gate"),
        supplementary_file = collapse_meta(meta, "supplementary_file_1")
      )
    }
  ), fill = TRUE)

  setorder(geo_metadata, gsm_accession)
  fwrite(geo_metadata, geo_metadata_file)
}

# ------------------------------------------------------------------------------
# Inventory plate files and build cell-level metadata
# ------------------------------------------------------------------------------

read_header <- function(path) {
  connection <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }

  on.exit(close(connection))
  header_line <- readLines(connection, n = 1L, warn = FALSE)

  if (length(header_line) != 1L || !nzchar(header_line)) {
    stop("Could not read a header from: ", path)
  }

  strsplit(header_line, "\t", fixed = TRUE)[[1L]]
}

plate_inventory_list <- vector("list", length(plate_files))
cell_metadata_list <- vector("list", length(plate_files))

for (index in seq_along(plate_files)) {
  path <- plate_files[index]
  file_name <- basename(path)
  gsm_accession <- sub("^(GSM[0-9]+)_.*$", "\\1", file_name)
  amplification_batch <- sub(
    "^GSM[0-9]+_(.+)\\.txt(\\.gz)?$",
    "\\1",
    file_name
  )

  header <- read_header(path)
  cell_ids <- header

  # The first header field is blank because row names contain gene symbols.
  if (length(cell_ids) > 0L && !nzchar(cell_ids[1L])) {
    cell_ids <- cell_ids[-1L]
  }

  gsm_key <- gsm_accession
  matched_metadata <- geo_metadata[gsm_accession == gsm_key]

  plate_inventory_list[[index]] <- data.table(
    file_name = file_name,
    gsm_accession = gsm_accession,
    amplification_batch_from_file = amplification_batch,
    file_size_bytes = file.info(path)$size,
    n_cell_columns = length(cell_ids),
    n_unique_cell_ids = uniqueN(cell_ids),
    n_duplicated_cell_ids_within_file = sum(duplicated(cell_ids)),
    first_cell_id = if (length(cell_ids) > 0L) cell_ids[1L] else NA_character_,
    last_cell_id = if (length(cell_ids) > 0L) tail(cell_ids, 1L) else NA_character_,
    matched_to_geo_metadata = nrow(matched_metadata) == 1L
  )

  if (nrow(matched_metadata) == 1L && length(cell_ids) > 0L) {
    cell_metadata_list[[index]] <- data.table(
      cell_id = cell_ids,
      gsm_accession = gsm_accession,
      amplification_batch = matched_metadata$amplification_batch,
      plate_id = matched_metadata$plate_id,
      patient_id = matched_metadata$patient_id,
      sample_source = matched_metadata$sample_source,
      facs_gate = matched_metadata$facs_gate,
      source_name = matched_metadata$source_name,
      count_file = file_name
    )
  } else {
    cell_metadata_list[[index]] <- data.table(
      cell_id = cell_ids,
      gsm_accession = gsm_accession,
      amplification_batch = amplification_batch,
      plate_id = NA_character_,
      patient_id = NA_character_,
      sample_source = NA_character_,
      facs_gate = NA_character_,
      source_name = NA_character_,
      count_file = file_name
    )
  }
}

plate_inventory <- rbindlist(plate_inventory_list, fill = TRUE)
cell_metadata <- rbindlist(cell_metadata_list, fill = TRUE)

cell_metadata[, duplicated_cell_id_across_files := duplicated(cell_id) |
  duplicated(cell_id, fromLast = TRUE)]

fwrite(
  plate_inventory,
  file.path(results_dir, "11_plate_file_inventory.csv")
)

fwrite(
  cell_metadata,
  file.path(
    analysis_ready_dir,
    "11_GSE123139_cell_metadata.csv.gz"
  )
)

# ------------------------------------------------------------------------------
# Confirm genes required for CD8-state annotation
# ------------------------------------------------------------------------------

representative_file <- plate_files[1L]
representative_genes <- fread(
  representative_file,
  select = 1L,
  header = TRUE,
  data.table = FALSE,
  check.names = FALSE
)[[1L]]

required_markers <- list(
  lineage = c("PTPRC", "CD3D", "CD3E", "TRAC", "CD8A", "CD8B"),
  cd43 = c("SPN"),
  dysfunctional = c("PDCD1", "HAVCR2", "LAG3", "TIGIT", "TOX", "CXCL13"),
  transitional_memory_like = c("TCF7", "IL7R", "CCR7", "LTB", "MAL"),
  cytotoxic_effector = c("NKG7", "GNLY", "PRF1", "GZMB", "GZMH", "FGFBP2"),
  proliferating = c("MKI67", "TOP2A", "STMN1", "TUBA1B")
)

marker_availability <- rbindlist(lapply(names(required_markers), function(group) {
  data.table(
    marker_group = group,
    gene = required_markers[[group]],
    present = required_markers[[group]] %in% representative_genes
  )
}))

fwrite(
  marker_availability,
  file.path(results_dir, "11_marker_availability.csv")
)

# ------------------------------------------------------------------------------
# QC summaries and provenance
# ------------------------------------------------------------------------------

dataset_qc <- data.table(
  metric = c(
    "geo_samples",
    "extracted_plate_files",
    "plate_files_matched_to_geo_metadata",
    "total_cell_columns",
    "unique_cell_ids",
    "duplicated_cell_ids_across_files",
    "unique_patients",
    "tumor_cell_columns",
    "pbmc_cell_columns",
    "required_markers_present",
    "required_markers_tested"
  ),
  value = as.character(c(
    nrow(geo_metadata),
    nrow(plate_inventory),
    sum(plate_inventory$matched_to_geo_metadata),
    nrow(cell_metadata),
    uniqueN(cell_metadata$cell_id),
    uniqueN(cell_metadata[duplicated_cell_id_across_files == TRUE, cell_id]),
    uniqueN(cell_metadata[!is.na(patient_id), patient_id]),
    nrow(cell_metadata[tolower(sample_source) == "tumor"]),
    nrow(cell_metadata[tolower(sample_source) == "pbmc"]),
    sum(marker_availability$present),
    nrow(marker_availability)
  ))
)

fwrite(
  dataset_qc,
  file.path(results_dir, "11_dataset_qc.csv")
)

download_manifest <- data.table(
  file = c(basename(raw_tar), basename(tcr_file)),
  source_url = c(raw_url, tcr_url),
  status = c(raw_download_status, tcr_download_status),
  size_bytes = c(file.info(raw_tar)$size, file.info(tcr_file)$size),
  md5 = unname(tools::md5sum(c(raw_tar, tcr_file)))
)

fwrite(
  download_manifest,
  file.path(results_dir, "11_download_manifest.csv")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(results_dir, "11_sessionInfo.txt")
)

message("Step 11 complete.")
message("Plate files: ", nrow(plate_inventory))
message("Cell columns: ", nrow(cell_metadata))
message("Unique patients: ", uniqueN(cell_metadata$patient_id, na.rm = TRUE))
message("Results: ", results_dir)

print(dataset_qc)
print(marker_availability)