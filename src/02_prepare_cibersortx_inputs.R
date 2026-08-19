#!/usr/bin/env Rscript

# ==============================================================================
# Prepare TCGA TPM matrices and CD43 groups for CIBERSORTx
# ==============================================================================
#
# Input (created by 01_download_tcga_expression.R):
#   data/tcga_star_counts/tpm_matrices/<CANCER>_TPM_primary_samples.tsv.gz
#   data/tcga_star_counts/sample_metadata/<CANCER>_primary_sample_metadata.csv
#
# Output:
#   data/cibersortx_input/mixture_files/<CANCER>_CIBERSORTx_TPM.txt
#   data/cibersortx_input/cd43_groups/<CANCER>_CD43_groups.csv
#   results/02/02_<CANCER>_preparation_QC.csv
#   results/02/02_preparation_summary.csv
#   results/02/02_project_configuration.csv
#   results/02/02_sessionInfo.txt
#
# Run every available cancer:
#   Rscript src/02_prepare_cibersortx_inputs.R
#
# Run selected cancers:
#   Rscript src/02_prepare_cibersortx_inputs.R PAAD LUAD
#
# Required CRAN package:
#   data.table
# ==============================================================================

options(stringsAsFactors = FALSE)

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop(
    "Package 'data.table' is required. Install it with: ",
    "install.packages('data.table')"
  )
}

suppressPackageStartupMessages(library(data.table))

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------

project_config <- data.table(
  cancer = c(
    "GBM", "LAML", "PAAD", "STAD", "TGCT",
    "THYM", "SKCM", "UVM", "LUAD", "LUSC"
  ),
  project = c(
    "TCGA-GBM", "TCGA-LAML", "TCGA-PAAD", "TCGA-STAD", "TCGA-TGCT",
    "TCGA-THYM", "TCGA-SKCM", "TCGA-UVM", "TCGA-LUAD", "TCGA-LUSC"
  ),
  cd43_direction = c(
    rep("overexpressed", 8),
    rep("underexpressed", 2)
  )
)

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

  project_config <- project_config[match(requested_cancers, cancer)]
}

input_root <- normalizePath(
  file.path(getwd(), "data", "tcga_star_counts"),
  mustWork = FALSE
)

prepared_data_root <- normalizePath(
  file.path(getwd(), "data", "cibersortx_input"),
  mustWork = FALSE
)

results_dir <- normalizePath(
  file.path(getwd(), "results", "02"),
  mustWork = FALSE
)

tpm_input_dir <- file.path(input_root, "tpm_matrices")
metadata_input_dir <- file.path(input_root, "sample_metadata")
mixture_output_dir <- file.path(prepared_data_root, "mixture_files")
groups_output_dir <- file.path(prepared_data_root, "cd43_groups")

invisible(lapply(
  c(prepared_data_root, mixture_output_dir, groups_output_dir, results_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

fwrite(
  project_config,
  file.path(results_dir, "02_project_configuration.csv")
)

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

timestamp_message <- function(...) {
  message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ...)
}

coerce_expression_columns <- function(x, columns) {
  for (column in columns) {
    if (!is.numeric(x[[column]])) {
      converted <- suppressWarnings(as.numeric(x[[column]]))

      bad_values <- is.na(converted) & !is.na(x[[column]])
      if (any(bad_values)) {
        stop(
          "Column '", column, "' contains ", sum(bad_values),
          " values that cannot be converted to numeric."
        )
      }

      set(x, j = column, value = converted)
    }
  }

  x
}

prepare_one_cancer <- function(cancer, project, cd43_direction) {
  timestamp_message("Preparing ", cancer)

  tpm_input_file <- file.path(
    tpm_input_dir,
    paste0(cancer, "_TPM_primary_samples.tsv.gz")
  )
  metadata_input_file <- file.path(
    metadata_input_dir,
    paste0(cancer, "_primary_sample_metadata.csv")
  )

  mixture_output_file <- file.path(
    mixture_output_dir,
    paste0(cancer, "_CIBERSORTx_TPM.txt")
  )
  groups_output_file <- file.path(
    groups_output_dir,
    paste0(cancer, "_CD43_groups.csv")
  )
  qc_output_file <- file.path(
    results_dir,
    paste0("02_", cancer, "_preparation_QC.csv")
  )

  if (!file.exists(tpm_input_file)) {
    timestamp_message("Skipping ", cancer, ": TPM input is not available yet")

    return(data.table(
      cancer = cancer,
      project = project,
      cd43_direction = cd43_direction,
      input_genes = NA_integer_,
      output_genes = NA_integer_,
      samples = NA_integer_,
      cd43_low = NA_integer_,
      cd43_intermediate = NA_integer_,
      cd43_high = NA_integer_,
      q1_tpm = NA_real_,
      q3_tpm = NA_real_,
      status = "skipped: TPM input not available"
    ))
  }

  expression_data <- fread(
    tpm_input_file,
    check.names = FALSE,
    showProgress = TRUE
  )

  if (ncol(expression_data) < 2L) {
    stop(cancer, " input must contain one gene column and at least one sample.")
  }

  if (names(expression_data)[1] != "GeneSymbol") {
    stop(
      cancer, " first column must be named 'GeneSymbol'; found '",
      names(expression_data)[1], "'."
    )
  }

  input_gene_count <- nrow(expression_data)
  sample_columns <- setdiff(names(expression_data), "GeneSymbol")

  if (anyDuplicated(sample_columns)) {
    duplicated_samples <- unique(sample_columns[duplicated(sample_columns)])
    stop(
      cancer, " contains duplicated sample columns: ",
      paste(duplicated_samples, collapse = ", ")
    )
  }

  expression_data <- coerce_expression_columns(expression_data, sample_columns)

  expression_matrix <- as.matrix(expression_data[, ..sample_columns])

  if (anyNA(expression_matrix)) {
    stop(cancer, " expression matrix contains missing values.")
  }

  if (any(!is.finite(expression_matrix))) {
    stop(cancer, " expression matrix contains non-finite values.")
  }

  if (any(expression_matrix < 0)) {
    stop(cancer, " expression matrix contains negative TPM values.")
  }

  # Undo the temporary suffix introduced by script 01, then collapse Ensembl
  # rows that map to the same official gene symbol by summing their TPM values.
  # Summation preserves the total abundance assigned to the shared symbol.
  expression_data[, GeneSymbol := trimws(as.character(GeneSymbol))]
  expression_data[, GeneSymbol := sub(
    "__duplicate_[0-9]+$", "", GeneSymbol
  )]

  invalid_symbol <- (
    is.na(expression_data$GeneSymbol) |
      expression_data$GeneSymbol == "" |
      grepl("^ENSG[0-9]+", expression_data$GeneSymbol)
  )

  removed_invalid_symbols <- sum(invalid_symbol)
  expression_data <- expression_data[!invalid_symbol]

  genes_before_collapse <- nrow(expression_data)
  duplicated_symbol_rows <- sum(duplicated(expression_data$GeneSymbol))

  collapsed_data <- expression_data[
    ,
    lapply(.SD, sum),
    by = GeneSymbol,
    .SDcols = sample_columns
  ]

  collapsed_matrix <- as.matrix(collapsed_data[, ..sample_columns])
  all_zero <- rowSums(collapsed_matrix) == 0
  removed_all_zero_genes <- sum(all_zero)
  collapsed_data <- collapsed_data[!all_zero]

  if (anyDuplicated(collapsed_data$GeneSymbol)) {
    stop(cancer, " still contains duplicated symbols after collapsing.")
  }

  spn_row <- which(collapsed_data$GeneSymbol == "SPN")

  if (length(spn_row) != 1L) {
    stop(
      cancer, " must contain exactly one SPN row after gene cleanup; found ",
      length(spn_row), "."
    )
  }

  spn_tpm <- as.numeric(unlist(
    collapsed_data[spn_row, ..sample_columns],
    use.names = FALSE
  ))
  names(spn_tpm) <- sample_columns
  spn_log2 <- log2(spn_tpm + 1)

  q1_tpm <- unname(quantile(spn_tpm, probs = 0.25, type = 7))
  q3_tpm <- unname(quantile(spn_tpm, probs = 0.75, type = 7))
  q1_log2 <- unname(quantile(spn_log2, probs = 0.25, type = 7))
  q3_log2 <- unname(quantile(spn_log2, probs = 0.75, type = 7))

  if (q1_tpm == q3_tpm) {
    stop(
      cancer, " has identical Q1 and Q3 SPN TPM thresholds (", q1_tpm,
      "); quartile groups cannot be separated."
    )
  }

  cd43_group <- fifelse(
    spn_tpm <= q1_tpm,
    "CD43_low",
    fifelse(spn_tpm >= q3_tpm, "CD43_high", "intermediate")
  )

  group_data <- data.table(
    sample_id = sample_columns,
    cancer = cancer,
    project = project,
    cd43_tumor_vs_normal_direction = cd43_direction,
    SPN_TPM = spn_tpm,
    SPN_log2_TPM_plus_1 = spn_log2,
    CD43_group = cd43_group,
    Q1_TPM = q1_tpm,
    Q3_TPM = q3_tpm,
    Q1_log2_TPM_plus_1 = q1_log2,
    Q3_log2_TPM_plus_1 = q3_log2
  )

  # Append the TCGAbiolinks sample metadata without changing sample order.
  if (file.exists(metadata_input_file)) {
    sample_metadata <- fread(metadata_input_file, check.names = FALSE)

    if ("prepared_sample_id" %in% names(sample_metadata)) {
      metadata_index <- match(group_data$sample_id, sample_metadata$prepared_sample_id)

      if (anyNA(metadata_index)) {
        warning(
          cancer, ": metadata did not match ", sum(is.na(metadata_index)),
          " expression sample(s)."
        )
      }

      extra_metadata_columns <- setdiff(
        names(sample_metadata),
        c(
          "prepared_sample_id", "analysis_cancer", "analysis_project",
          "cd43_tumor_vs_normal_direction"
        )
      )

      matched_metadata <- sample_metadata[
        metadata_index,
        ..extra_metadata_columns
      ]
      group_data <- cbind(group_data, matched_metadata)
    } else {
      warning(
        cancer,
        ": metadata file does not contain 'prepared_sample_id'; not joined."
      )
    }
  } else {
    warning(cancer, ": sample metadata file is not available; groups saved alone.")
  }

  # CIBERSORTx requires a plain, tab-delimited, genes-by-samples mixture file.
  fwrite(
    collapsed_data,
    mixture_output_file,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )
  fwrite(group_data, groups_output_file)

  group_counts <- table(
    factor(
      cd43_group,
      levels = c("CD43_low", "intermediate", "CD43_high")
    )
  )

  qc_data <- data.table(
    cancer = cancer,
    project = project,
    cd43_direction = cd43_direction,
    input_genes = input_gene_count,
    genes_after_symbol_filter = genes_before_collapse,
    removed_invalid_symbols = removed_invalid_symbols,
    duplicated_symbol_rows_collapsed = duplicated_symbol_rows,
    removed_all_zero_genes = removed_all_zero_genes,
    output_genes = nrow(collapsed_data),
    samples = length(sample_columns),
    spn_min_tpm = min(spn_tpm),
    spn_median_tpm = median(spn_tpm),
    spn_max_tpm = max(spn_tpm),
    q1_tpm = q1_tpm,
    q3_tpm = q3_tpm,
    q1_log2_tpm_plus_1 = q1_log2,
    q3_log2_tpm_plus_1 = q3_log2,
    cd43_low = unname(group_counts["CD43_low"]),
    cd43_intermediate = unname(group_counts["intermediate"]),
    cd43_high = unname(group_counts["CD43_high"]),
    mixture_file = mixture_output_file,
    groups_file = groups_output_file,
    status = "completed"
  )

  fwrite(qc_data, qc_output_file)

  timestamp_message(
    "Completed ", cancer, ": ", nrow(collapsed_data), " genes x ",
    length(sample_columns), " samples; Q1 = ", signif(q1_tpm, 5),
    " TPM; Q3 = ", signif(q3_tpm, 5), " TPM"
  )

  qc_data
}

# ------------------------------------------------------------------------------
# Run the preparation
# ------------------------------------------------------------------------------

preparation_results <- vector("list", nrow(project_config))

for (i in seq_len(nrow(project_config))) {
  config_row <- project_config[i]

  preparation_results[[i]] <- tryCatch(
    prepare_one_cancer(
      cancer = config_row$cancer,
      project = config_row$project,
      cd43_direction = config_row$cd43_direction
    ),
    error = function(e) {
      timestamp_message("ERROR in ", config_row$cancer, ": ", conditionMessage(e))

      data.table(
        cancer = config_row$cancer,
        project = config_row$project,
        cd43_direction = config_row$cd43_direction,
        input_genes = NA_integer_,
        output_genes = NA_integer_,
        samples = NA_integer_,
        cd43_low = NA_integer_,
        cd43_intermediate = NA_integer_,
        cd43_high = NA_integer_,
        q1_tpm = NA_real_,
        q3_tpm = NA_real_,
        status = paste0("failed: ", conditionMessage(e))
      )
    }
  )

  current_summary <- rbindlist(
    preparation_results[seq_len(i)],
    fill = TRUE,
    use.names = TRUE
  )
  fwrite(
    current_summary,
    file.path(results_dir, "02_preparation_summary.csv")
  )
}

final_summary <- rbindlist(preparation_results, fill = TRUE, use.names = TRUE)

timestamp_message("CIBERSORTx preparation finished for all requested cancers.")
print(
  final_summary[
    ,
    .(
      cancer, samples, output_genes, cd43_low,
      cd43_intermediate, cd43_high, status
    )
  ]
)

writeLines(
  capture.output(sessionInfo()),
  file.path(results_dir, "02_sessionInfo.txt")
)