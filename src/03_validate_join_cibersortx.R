#!/usr/bin/env Rscript

# ==============================================================================
# Validate CIBERSORTx LM22 results and join them to CD43 expression groups
# ==============================================================================
#
# Inputs:
#   data/cibersortx_output/<CANCER>/CIBERSORTx_<CANCER>_Adjusted.txt
#   data/cibersortx_input/cd43_groups/<CANCER>_CD43_groups.csv
#
# Outputs:
#   data/analysis_ready/cancer_tables/<CANCER>_analysis_ready.tsv.gz
#   data/analysis_ready/all_cancers_analysis_ready.tsv.gz
#   results/03/03_<CANCER>_deconvolution_QC.csv
#   results/03/03_deconvolution_QC_summary.csv
#   results/03/03_project_configuration.csv
#   results/03/03_sessionInfo.txt
#   results/03/figures/03_CIBERSORTx_fit_quality_by_cancer.png
#
# Run all cancers:
#   Rscript src/03_validate_join_cibersortx.R
#
# Run selected cancers:
#   Rscript src/03_validate_join_cibersortx.R GBM PAAD
# ==============================================================================

options(stringsAsFactors = FALSE)

required_packages <- c("data.table", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required package(s): ", paste(missing_packages, collapse = ", "),
    ". Install with install.packages()."
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

project_config <- data.table(
  cancer = c(
    "GBM", "LAML", "PAAD", "STAD", "TGCT",
    "THYM", "SKCM", "UVM", "LUAD", "LUSC"
  ),
  project = c(
    "TCGA-GBM", "TCGA-LAML", "TCGA-PAAD", "TCGA-STAD", "TCGA-TGCT",
    "TCGA-THYM", "TCGA-SKCM", "TCGA-UVM", "TCGA-LUAD", "TCGA-LUSC"
  ),
  expected_samples = c(372L, 151L, 178L, 412L, 150L, 120L, 103L, 80L, 540L, 511L),
  cd43_direction = c(rep("overexpressed", 8), rep("underexpressed", 2))
)

requested_cancers <- toupper(commandArgs(trailingOnly = TRUE))

if (length(requested_cancers) > 0L) {
  unknown <- setdiff(requested_cancers, project_config$cancer)

  if (length(unknown) > 0L) {
    stop(
      "Unknown cancer abbreviation(s): ", paste(unknown, collapse = ", "),
      "\nAllowed values: ", paste(project_config$cancer, collapse = ", ")
    )
  }

  project_config <- project_config[match(requested_cancers, cancer)]
}

fractions_root <- normalizePath(
  file.path(getwd(), "data", "cibersortx_output"),
  mustWork = FALSE
)
groups_root <- normalizePath(
  file.path(getwd(), "data", "cibersortx_input", "cd43_groups"),
  mustWork = FALSE
)
analysis_data_root <- normalizePath(
  file.path(getwd(), "data", "analysis_ready"),
  mustWork = FALSE
)
results_dir <- normalizePath(
  file.path(getwd(), "results", "03"),
  mustWork = FALSE
)
cancer_output_dir <- file.path(analysis_data_root, "cancer_tables")
figure_output_dir <- file.path(results_dir, "figures")

invisible(lapply(
  c(analysis_data_root, cancer_output_dir, results_dir, figure_output_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

fwrite(
  project_config,
  file.path(results_dir, "03_project_configuration.csv")
)

timestamp_message <- function(...) {
  message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ...)
}

normalize_column_name <- function(x) {
  tolower(gsub("[^[:alnum:]]", "", x))
}

find_one_column <- function(column_names, normalized_target, description) {
  normalized_names <- normalize_column_name(column_names)
  matches <- which(normalized_names == normalized_target)

  if (length(matches) != 1L) {
    stop(
      "Expected exactly one ", description, " column; found ",
      length(matches), ". Available columns: ",
      paste(column_names, collapse = ", ")
    )
  }

  column_names[matches]
}

coerce_numeric_columns <- function(x, columns, context) {
  for (column in columns) {
    original <- x[[column]]
    converted <- suppressWarnings(as.numeric(original))
    invalid <- is.na(converted) & !is.na(original)

    if (any(invalid)) {
      stop(
        context, ": column '", column, "' contains ", sum(invalid),
        " nonnumeric value(s)."
      )
    }

    set(x, j = column, value = converted)
  }

  x
}

process_one_cancer <- function(
    cancer,
    project,
    expected_samples,
    cd43_direction,
    fraction_sum_tolerance = 1e-4,
    fraction_value_tolerance = 1e-8) {

  timestamp_message("Validating ", cancer)

  fractions_file <- file.path(
    fractions_root,
    cancer,
    paste0("CIBERSORTx_", cancer, "_Adjusted.txt")
  )
  groups_file <- file.path(groups_root, paste0(cancer, "_CD43_groups.csv"))
  joined_file <- file.path(
    cancer_output_dir,
    paste0(cancer, "_analysis_ready.tsv.gz")
  )
  qc_file <- file.path(
    results_dir,
    paste0("03_", cancer, "_deconvolution_QC.csv")
  )

  if (!file.exists(fractions_file)) {
    stop(cancer, ": fraction result is missing: ", fractions_file)
  }

  if (!file.exists(groups_file)) {
    stop(cancer, ": CD43 group file is missing: ", groups_file)
  }

  fractions <- fread(fractions_file, check.names = FALSE)
  groups <- fread(groups_file, check.names = FALSE)

  if (ncol(fractions) < 26L) {
    stop(
      cancer, ": expected a sample column, 22 fractions, and three QC columns; ",
      "found only ", ncol(fractions), " columns."
    )
  }

  if (nrow(fractions) != expected_samples) {
    stop(
      cancer, ": expected ", expected_samples,
      " samples but fraction file contains ", nrow(fractions), "."
    )
  }

  # The first CIBERSORTx column is normally named Mixture.
  first_column <- names(fractions)[1]
  setnames(fractions, first_column, "sample_id")
  fractions[, sample_id := as.character(sample_id)]

  p_value_column <- find_one_column(names(fractions), "pvalue", "P-value")
  correlation_column <- find_one_column(
    names(fractions), "correlation", "Correlation"
  )
  rmse_column <- find_one_column(names(fractions), "rmse", "RMSE")

  setnames(
    fractions,
    c(p_value_column, correlation_column, rmse_column),
    c("CIBERSORTx_P_value", "CIBERSORTx_Correlation", "CIBERSORTx_RMSE")
  )

  qc_columns <- c(
    "CIBERSORTx_P_value", "CIBERSORTx_Correlation", "CIBERSORTx_RMSE"
  )
  cell_columns <- setdiff(names(fractions), c("sample_id", qc_columns))

  if (length(cell_columns) != 22L) {
    stop(
      cancer, ": expected 22 LM22 fraction columns but found ",
      length(cell_columns), ". Columns classified as fractions: ",
      paste(cell_columns, collapse = ", ")
    )
  }

  fractions <- coerce_numeric_columns(
    fractions,
    c(cell_columns, qc_columns),
    paste0(cancer, " CIBERSORTx result")
  )

  if (anyDuplicated(fractions$sample_id)) {
    stop(cancer, ": duplicated sample identifiers in fraction results.")
  }

  fraction_matrix <- as.matrix(fractions[, ..cell_columns])

  if (anyNA(fraction_matrix) || any(!is.finite(fraction_matrix))) {
    stop(cancer, ": immune fractions contain missing or non-finite values.")
  }

  negative_values <- sum(fraction_matrix < -fraction_value_tolerance)
  above_one_values <- sum(fraction_matrix > 1 + fraction_value_tolerance)
  fraction_sums <- rowSums(fraction_matrix)
  bad_fraction_sums <- sum(abs(fraction_sums - 1) > fraction_sum_tolerance)

  if (negative_values > 0L) {
    stop(cancer, ": found ", negative_values, " fraction values below zero.")
  }

  if (above_one_values > 0L) {
    stop(cancer, ": found ", above_one_values, " fraction values above one.")
  }

  if (bad_fraction_sums > 0L) {
    stop(
      cancer, ": ", bad_fraction_sums,
      " sample(s) have relative fractions that do not sum to one within ",
      fraction_sum_tolerance, "."
    )
  }

  if (anyNA(fractions$CIBERSORTx_P_value) ||
      anyNA(fractions$CIBERSORTx_Correlation) ||
      anyNA(fractions$CIBERSORTx_RMSE)) {
    stop(cancer, ": one or more CIBERSORTx QC metrics are missing.")
  }

  if (!"sample_id" %in% names(groups)) {
    stop(cancer, ": CD43 group table does not contain 'sample_id'.")
  }

  groups[, sample_id := as.character(sample_id)]

  if (anyDuplicated(groups$sample_id)) {
    stop(cancer, ": duplicated sample identifiers in CD43 group table.")
  }

  missing_groups <- setdiff(fractions$sample_id, groups$sample_id)
  extra_groups <- setdiff(groups$sample_id, fractions$sample_id)

  if (length(missing_groups) > 0L || length(extra_groups) > 0L) {
    stop(
      cancer, ": sample identifiers do not match. Missing CD43 groups for ",
      length(missing_groups), " fraction sample(s); group table contains ",
      length(extra_groups), " unmatched sample(s)."
    )
  }

  groups <- groups[match(fractions$sample_id, sample_id)]

  # Prevent duplicated column names when appending the rich TCGA metadata.
  extra_group_columns <- setdiff(names(groups), "sample_id")
  duplicated_names <- intersect(extra_group_columns, names(fractions))

  if (length(duplicated_names) > 0L) {
    setnames(
      groups,
      duplicated_names,
      paste0("metadata_", duplicated_names)
    )
  }

  joined <- cbind(
    fractions,
    groups[, setdiff(names(groups), "sample_id"), with = FALSE]
  )

  joined[, LM22_fraction_sum := fraction_sums]
  joined[, CIBERSORTx_pass_P_0_05 := CIBERSORTx_P_value <= 0.05]
  joined[, analysis_include_extreme_quartiles := CD43_group %in% c(
    "CD43_low", "CD43_high"
  )]

  if (!all(joined$cancer == cancer)) {
    stop(cancer, ": joined metadata contains an inconsistent cancer label.")
  }

  fwrite(
    joined,
    joined_file,
    sep = "\t",
    quote = FALSE,
    na = "NA",
    compress = "gzip"
  )

  group_counts <- table(
    factor(
      joined$CD43_group,
      levels = c("CD43_low", "intermediate", "CD43_high")
    )
  )

  qc_summary <- data.table(
    cancer = cancer,
    project = project,
    cd43_direction = cd43_direction,
    samples = nrow(joined),
    lm22_cell_types = length(cell_columns),
    cd43_low = unname(group_counts["CD43_low"]),
    cd43_intermediate = unname(group_counts["intermediate"]),
    cd43_high = unname(group_counts["CD43_high"]),
    p_le_0_05_n = sum(joined$CIBERSORTx_pass_P_0_05),
    p_le_0_05_percent = 100 * mean(joined$CIBERSORTx_pass_P_0_05),
    p_value_min = min(joined$CIBERSORTx_P_value),
    p_value_median = median(joined$CIBERSORTx_P_value),
    p_value_max = max(joined$CIBERSORTx_P_value),
    correlation_min = min(joined$CIBERSORTx_Correlation),
    correlation_mean = mean(joined$CIBERSORTx_Correlation),
    correlation_median = median(joined$CIBERSORTx_Correlation),
    correlation_max = max(joined$CIBERSORTx_Correlation),
    rmse_min = min(joined$CIBERSORTx_RMSE),
    rmse_mean = mean(joined$CIBERSORTx_RMSE),
    rmse_median = median(joined$CIBERSORTx_RMSE),
    rmse_max = max(joined$CIBERSORTx_RMSE),
    fraction_sum_min = min(fraction_sums),
    fraction_sum_median = median(fraction_sums),
    fraction_sum_max = max(fraction_sums),
    negative_fraction_values = negative_values,
    above_one_fraction_values = above_one_values,
    bad_fraction_sum_samples = bad_fraction_sums,
    analysis_ready_file = joined_file,
    status = "completed"
  )

  fwrite(qc_summary, qc_file)

  timestamp_message(
    "Completed ", cancer, ": ", nrow(joined), " samples; ",
    sum(joined$CIBERSORTx_pass_P_0_05), " (",
    sprintf("%.1f", 100 * mean(joined$CIBERSORTx_pass_P_0_05)),
    "%) with P <= 0.05; median correlation = ",
    sprintf("%.3f", median(joined$CIBERSORTx_Correlation)),
    "; median RMSE = ", sprintf("%.3f", median(joined$CIBERSORTx_RMSE))
  )

  list(joined = joined, qc = qc_summary)
}

all_joined <- vector("list", nrow(project_config))
all_qc <- vector("list", nrow(project_config))

for (i in seq_len(nrow(project_config))) {
  config <- project_config[i]

  result <- tryCatch(
    process_one_cancer(
      cancer = config$cancer,
      project = config$project,
      expected_samples = config$expected_samples,
      cd43_direction = config$cd43_direction
    ),
    error = function(e) {
      timestamp_message("ERROR in ", config$cancer, ": ", conditionMessage(e))

      list(
        joined = NULL,
        qc = data.table(
          cancer = config$cancer,
          project = config$project,
          cd43_direction = config$cd43_direction,
          samples = NA_integer_,
          status = paste0("failed: ", conditionMessage(e))
        )
      )
    }
  )

  all_joined[[i]] <- result$joined
  all_qc[[i]] <- result$qc

  current_qc <- rbindlist(all_qc[seq_len(i)], fill = TRUE, use.names = TRUE)
  fwrite(
    current_qc,
    file.path(results_dir, "03_deconvolution_QC_summary.csv")
  )
}

qc_summary <- rbindlist(all_qc, fill = TRUE, use.names = TRUE)
completed_tables <- Filter(Negate(is.null), all_joined)

if (length(completed_tables) > 0L) {
  combined_data <- rbindlist(completed_tables, fill = TRUE, use.names = TRUE)
  combined_file <- file.path(
    analysis_data_root,
    "all_cancers_analysis_ready.tsv.gz"
  )

  fwrite(
    combined_data,
    combined_file,
    sep = "\t",
    quote = FALSE,
    na = "NA",
    compress = "gzip"
  )

  timestamp_message(
    "Combined analysis-ready table contains ", nrow(combined_data),
    " samples across ", uniqueN(combined_data$cancer), " cancer(s)."
  )

  # Match the aesthetic of the longitudinal CIBERSORTx QC figure used in the
  # previous project, replacing treatment visit with TCGA cancer type.
  fit_plot_data <- melt(
    combined_data[
      ,
      .(
        sample_id,
        cancer,
        Correlation = CIBERSORTx_Correlation,
        RMSE = CIBERSORTx_RMSE
      )
    ],
    id.vars = c("sample_id", "cancer"),
    measure.vars = c("Correlation", "RMSE"),
    variable.name = "metric",
    value.name = "value"
  )

  fit_plot_data[, cancer := factor(cancer, levels = project_config$cancer)]
  fit_plot_data[, metric := factor(
    metric,
    levels = c("Correlation", "RMSE")
  )]

  fit_plot <- ggplot(
    fit_plot_data,
    aes(x = cancer, y = value)
  ) +
    geom_boxplot(
      outlier.shape = NA,
      fill = "#DCE8EE",
      colour = "#315D6B",
      linewidth = 0.5,
      width = 0.65
    ) +
    geom_jitter(
      shape = 21,
      fill = "#76B7D5",
      colour = "black",
      stroke = 0.15,
      width = 0.14,
      alpha = 0.30,
      size = 1
    ) +
    facet_wrap(
      ~ metric,
      scales = "free_y",
      ncol = 1
    ) +
    labs(
      title = "CIBERSORTx fit quality across TCGA cancers",
      x = "Cancer type",
      y = NULL
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 8, colour = "grey30"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.ticks = element_blank(),
      axis.line = element_line(colour = "grey25", linewidth = 0.4),
      strip.text = element_text(
        face = "bold",
        size = 9,
        margin = margin(t = 4, b = 4)
      ),
      strip.background = element_blank(),
      panel.grid = element_blank(),
      panel.spacing = grid::unit(0.9, "lines"),
      plot.margin = margin(8, 8, 8, 8)
    )

  ggsave(
    file.path(figure_output_dir, "03_CIBERSORTx_fit_quality_by_cancer.png"),
    fit_plot,
    width = 9,
    height = 7,
    dpi = 600
  )

  timestamp_message(
    "Saved CIBERSORTx fit-quality figure in ", figure_output_dir
  )
}

print(
  qc_summary[
    ,
    .(
      cancer, samples, p_le_0_05_n, p_le_0_05_percent,
      correlation_median, rmse_median, status
    )
  ]
)

writeLines(
  capture.output(sessionInfo()),
  file.path(results_dir, "03_sessionInfo.txt")
)