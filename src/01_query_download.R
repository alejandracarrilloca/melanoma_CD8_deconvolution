#!/usr/bin/env Rscript

# ==============================================================================
# Download TCGA STAR-Counts RNA-seq data for the CD43 deconvolution project
# ==============================================================================
#
# This script:
#   1. Queries each TCGA project separately with TCGAbiolinks.
#   2. Downloads open-access STAR - Counts gene-expression files.
#   3. Prepares one SummarizedExperiment per cancer.
#   4. Extracts the GDC TPM assay without log transformation.
#   5. Saves sample metadata and a genes-by-samples TPM matrix.
#
# The TPM matrices produced here are intended for later CIBERSORTx preparation.
# Do not submit the RDS files directly to CIBERSORTx.
#
# Run from the project root with:
#   Rscript src/01_download_tcga_expression.R
#
# To run only selected projects, provide their TCGA abbreviations:
#   Rscript src/01_download_tcga_expression.R GBM PAAD STAD
#
# Required Bioconductor packages:
#   TCGAbiolinks, SummarizedExperiment
#
# One-time installation:
#   if (!requireNamespace("BiocManager", quietly = TRUE))
#       install.packages("BiocManager")
#   BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment"))
#
# Optional but recommended for faster compressed TSV writing:
#   install.packages("data.table")
# ==============================================================================

options(stringsAsFactors = FALSE)
options(timeout = max(3600, getOption("timeout")))

required_packages <- c("TCGAbiolinks", "SummarizedExperiment")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them with BiocManager before running this script."
  )
}

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
})

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------

project_config <- data.frame(
  cancer = c(
    "GBM", "LAML", "PAAD", "STAD", "TGCT",
    "THYM", "SKCM", "UVM", "LUAD", "LUSC"
  ),
  project = c(
    "TCGA-GBM", "TCGA-LAML", "TCGA-PAAD", "TCGA-STAD", "TCGA-TGCT",
    "TCGA-THYM", "TCGA-SKCM", "TCGA-UVM", "TCGA-LUAD", "TCGA-LUSC"
  ),
  cd43_direction = c(
    "overexpressed", "overexpressed", "overexpressed", "overexpressed",
    "overexpressed", "overexpressed", "overexpressed", "overexpressed",
    "underexpressed", "underexpressed"
  ),
  sample_type = c(
    "Primary Tumor",
    "Primary Blood Derived Cancer - Peripheral Blood",
    "Primary Tumor", "Primary Tumor", "Primary Tumor",
    "Primary Tumor", "Primary Tumor", "Primary Tumor",
    "Primary Tumor", "Primary Tumor"
  )
)

# If cancer abbreviations were provided as command-line arguments, run only those.
requested_cancers <- toupper(commandArgs(trailingOnly = TRUE))

if (length(requested_cancers) > 0L) {
  unknown_cancers <- setdiff(requested_cancers, project_config$cancer)

  if (length(unknown_cancers) > 0L) {
    stop(
      "Unknown cancer abbreviation(s): ",
      paste(unknown_cancers, collapse = ", "),
      "\nAllowed values: ",
      paste(project_config$cancer, collapse = ", ")
    )
  }

  project_config <- project_config[
    match(requested_cancers, project_config$cancer),
    ,
    drop = FALSE
  ]
}

# Reusable downloaded and processed expression data remain under data/.
# Script-specific summaries and reproducibility records are written under the
# matching numbered results directory.
data_output_root <- normalizePath(
  file.path(getwd(), "data", "tcga_star_counts"),
  mustWork = FALSE
)

results_dir <- normalizePath(
  file.path(getwd(), "results", "01"),
  mustWork = FALSE
)

raw_download_dir <- file.path(data_output_root, "gdc_download")
query_dir <- file.path(data_output_root, "queries")
se_dir <- file.path(data_output_root, "summarized_experiments")
tpm_dir <- file.path(data_output_root, "tpm_matrices")
metadata_dir <- file.path(data_output_root, "sample_metadata")

invisible(lapply(
  c(
    data_output_root, raw_download_dir, query_dir, se_dir,
    tpm_dir, metadata_dir, results_dir
  ),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

write.csv(
  project_config,
  file.path(results_dir, "01_project_configuration.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

timestamp_message <- function(...) {
  message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ...)
}

make_csv_safe <- function(x) {
  x <- as.data.frame(x, check.names = FALSE)

  collapse_list_element <- function(value) {
    if (is.null(value) || length(value) == 0L) {
      return(NA_character_)
    }

    flattened <- unlist(value, recursive = TRUE, use.names = FALSE)

    if (length(flattened) == 0L) {
      return(NA_character_)
    }

    paste(as.character(flattened), collapse = "|")
  }

  list_columns <- vapply(x, is.list, logical(1))

  x[list_columns] <- lapply(
    x[list_columns],
    function(column) {
      vapply(column, collapse_list_element, character(1))
    }
  )

  x
}

select_tpm_assay <- function(se) {
  available_assays <- SummarizedExperiment::assayNames(se)

  # Current GDC STAR-Counts data prepared by TCGAbiolinks normally uses
  # "tpm_unstrand". The fallback allows for minor naming changes.
  preferred_names <- c("tpm_unstrand", "tpm")
  exact_match <- preferred_names[preferred_names %in% available_assays]

  if (length(exact_match) > 0L) {
    return(exact_match[[1]])
  }

  partial_match <- grep("^tpm", available_assays, value = TRUE, ignore.case = TRUE)

  if (length(partial_match) == 1L) {
    return(partial_match)
  }

  stop(
    "Could not identify one TPM assay. Available assays: ",
    paste(available_assays, collapse = ", ")
  )
}

make_unique_gene_names <- function(row_data) {
  row_df <- as.data.frame(row_data)

  symbol_candidates <- c(
    "gene_name", "gene_symbol", "external_gene_name", "symbol"
  )
  symbol_column <- symbol_candidates[symbol_candidates %in% names(row_df)][1]

  id_candidates <- c("gene_id", "ensembl_gene_id")
  id_column <- id_candidates[id_candidates %in% names(row_df)][1]

  if (is.na(symbol_column)) {
    gene_symbols <- rownames(row_df)
  } else {
    gene_symbols <- as.character(row_df[[symbol_column]])
  }

  if (is.na(id_column)) {
    gene_ids <- rownames(row_df)
  } else {
    gene_ids <- as.character(row_df[[id_column]])
  }

  gene_ids <- sub("\\.[0-9]+$", "", gene_ids)
  missing_symbol <- is.na(gene_symbols) | gene_symbols == ""
  gene_symbols[missing_symbol] <- gene_ids[missing_symbol]

  # CIBERSORTx requires unique row identifiers. Duplicated symbols are retained
  # temporarily with unique suffixes. A later preparation script will collapse
  # them using a documented aggregation rule before CIBERSORTx submission.
  make.unique(gene_symbols, sep = "__duplicate_")
}

write_tpm_matrix <- function(tpm_matrix, row_data, output_file) {
  gene_names <- make_unique_gene_names(row_data)
  tpm_df <- data.frame(
    GeneSymbol = gene_names,
    as.data.frame(tpm_matrix, check.names = FALSE),
    check.names = FALSE
  )

  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(
      tpm_df,
      file = output_file,
      sep = "\t",
      quote = FALSE,
      na = "NA",
      compress = "gzip"
    )
  } else {
    connection <- gzfile(output_file, open = "wt")
    on.exit(close(connection), add = TRUE)
    write.table(
      tpm_df,
      file = connection,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE,
      na = "NA"
    )
  }
}

download_with_retries <- function(query, directory, attempts = 3L) {
  last_error <- NULL

  for (attempt in seq_len(attempts)) {
    timestamp_message("Download attempt ", attempt, " of ", attempts)

    result <- tryCatch(
      {
        GDCdownload(
          query = query,
          method = "api",
          directory = directory,
          files.per.chunk = 10
        )
        TRUE
      },
      error = function(e) {
        last_error <<- e
        timestamp_message("Download attempt failed: ", conditionMessage(e))
        FALSE
      }
    )

    if (isTRUE(result)) {
      return(invisible(TRUE))
    }
  }

  stop(
    "GDCdownload failed after ", attempts, " attempts. Last error: ",
    conditionMessage(last_error)
  )
}

process_project <- function(cancer, project, cd43_direction, sample_type) {
  timestamp_message("Starting ", project, " (", cancer, ")")

  query_file <- file.path(query_dir, paste0(cancer, "_query.rds"))
  query_results_file <- file.path(query_dir, paste0(cancer, "_query_results.csv"))
  se_file <- file.path(se_dir, paste0(cancer, "_STAR_Counts_SE.rds"))
  tpm_file <- file.path(tpm_dir, paste0(cancer, "_TPM_primary_samples.tsv.gz"))
  metadata_file <- file.path(
    metadata_dir,
    paste0(cancer, "_primary_sample_metadata.csv")
  )

  # A completed SummarizedExperiment is treated as the restart checkpoint.
  if (file.exists(se_file)) {
    timestamp_message("Prepared object already exists; loading ", se_file)
    se <- readRDS(se_file)
  } else {
    timestamp_message(
      "Querying ", project, " for sample type: ", sample_type
    )

    query <- GDCquery(
      project = project,
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts",
      sample.type = sample_type,
      access = "open"
    )

    query_results <- getResults(query)

    if (nrow(query_results) == 0L) {
      stop(
        "The query returned no files for ", project,
        " with sample type '", sample_type, "'."
      )
    }

    saveRDS(query, query_file)
    # Some GDC query fields are list columns. Flatten them only for the CSV
    # representation; the complete original query object is preserved in RDS.
    write.csv(
      make_csv_safe(query_results),
      query_results_file,
      row.names = FALSE
    )

    timestamp_message(
      "Query returned ", nrow(query_results), " files for ", project
    )

    download_with_retries(query, raw_download_dir)

    timestamp_message("Preparing SummarizedExperiment for ", project)
    se <- GDCprepare(
      query = query,
      directory = raw_download_dir,
      summarizedExperiment = TRUE
    )

    saveRDS(se, se_file, compress = FALSE)
    timestamp_message("Saved prepared object: ", se_file)
  }

  # Validate that the prepared object contains samples.
  if (ncol(se) == 0L) {
    stop("Prepared object contains no samples for ", project)
  }

  tpm_assay_name <- select_tpm_assay(se)
  tpm_matrix <- SummarizedExperiment::assay(se, tpm_assay_name)

  if (anyNA(tpm_matrix)) {
    warning(project, " TPM matrix contains missing values.")
  }

  if (any(tpm_matrix < 0, na.rm = TRUE)) {
    stop(project, " TPM matrix contains negative values.")
  }

  # Save metadata with explicit project annotations.
  sample_metadata <- as.data.frame(SummarizedExperiment::colData(se))
  sample_metadata$analysis_cancer <- cancer
  sample_metadata$analysis_project <- project
  sample_metadata$cd43_tumor_vs_normal_direction <- cd43_direction
  sample_metadata$requested_sample_type <- sample_type
  sample_metadata$prepared_sample_id <- colnames(se)

  # colData can also contain list-valued fields that base write.csv cannot
  # encode directly. Preserve the complete metadata in the SE RDS and flatten
  # list values only in the human-readable CSV export.
  write.csv(
    make_csv_safe(sample_metadata),
    metadata_file,
    row.names = FALSE
  )
  write_tpm_matrix(tpm_matrix, SummarizedExperiment::rowData(se), tpm_file)

  spn_candidates <- which(
    make_unique_gene_names(SummarizedExperiment::rowData(se)) == "SPN"
  )

  if (length(spn_candidates) != 1L) {
    warning(
      project, " contains ", length(spn_candidates),
      " rows with the exact gene symbol SPN after annotation."
    )
  }

  timestamp_message(
    "Completed ", project,
    ": ", nrow(tpm_matrix), " genes x ", ncol(tpm_matrix), " samples; ",
    "TPM assay = ", tpm_assay_name
  )

  data.frame(
    cancer = cancer,
    project = project,
    cd43_direction = cd43_direction,
    requested_sample_type = sample_type,
    genes = nrow(tpm_matrix),
    samples = ncol(tpm_matrix),
    tpm_assay = tpm_assay_name,
    spn_exact_rows = length(spn_candidates),
    se_file = se_file,
    tpm_file = tpm_file,
    metadata_file = metadata_file,
    status = "completed",
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
# Run all requested projects
# ------------------------------------------------------------------------------

run_summary <- vector("list", nrow(project_config))

for (i in seq_len(nrow(project_config))) {
  config_row <- project_config[i, , drop = FALSE]

  run_summary[[i]] <- tryCatch(
    process_project(
      cancer = config_row$cancer,
      project = config_row$project,
      cd43_direction = config_row$cd43_direction,
      sample_type = config_row$sample_type
    ),
    error = function(e) {
      timestamp_message(
        "ERROR in ", config_row$project, ": ", conditionMessage(e)
      )

      data.frame(
        cancer = config_row$cancer,
        project = config_row$project,
        cd43_direction = config_row$cd43_direction,
        requested_sample_type = config_row$sample_type,
        genes = NA_integer_,
        samples = NA_integer_,
        tpm_assay = NA_character_,
        spn_exact_rows = NA_integer_,
        se_file = NA_character_,
        tpm_file = NA_character_,
        metadata_file = NA_character_,
        status = paste0("failed: ", conditionMessage(e)),
        stringsAsFactors = FALSE
      )
    }
  )

  current_summary <- do.call(rbind, run_summary[seq_len(i)])
  write.csv(
    current_summary,
    file.path(results_dir, "01_download_run_summary.csv"),
    row.names = FALSE
  )
}

final_summary <- do.call(rbind, run_summary)

timestamp_message("All requested projects have been attempted.")
print(final_summary[, c("cancer", "project", "samples", "status")], row.names = FALSE)

if (any(final_summary$status != "completed")) {
  timestamp_message(
    "At least one project failed. Review results/01/",
    "01_download_run_summary.csv and rerun the failed project abbreviation(s)."
  )
}

session_file <- file.path(results_dir, "01_sessionInfo.txt")
writeLines(capture.output(sessionInfo()), session_file)