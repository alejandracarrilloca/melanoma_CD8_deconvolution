#!/usr/bin/env Rscript

# ==============================================================================
# Test CD43-high versus CD43-low immune differences and continuous SPN
# associations across TCGA cancers
# ==============================================================================
#
# Inputs:
#   data/analysis_ready/cancer_tables/<CANCER>_analysis_ready.tsv.gz
#
# Outputs:
#   results/04/04_<CANCER>_CD43_high_vs_low.csv
#   results/04/04_<CANCER>_SPN_continuous_correlations.csv
#   results/04/04_all_cancers_CD43_high_vs_low.csv
#   results/04/04_all_cancers_SPN_correlations.csv
#   results/04/04_cross_cancer_consistency_summary.csv
#   results/04/04_cross_cancer_spearman_rho_matrix.csv
#   results/04/04_cross_cancer_rank_biserial_matrix.csv
#   results/04/04_analysis_configuration.csv
#   results/04/04_sessionInfo.txt
#
# Run from the project root with:
#   Rscript src/04_test_cd43_immune_associations.R
# ==============================================================================

options(stringsAsFactors = FALSE)

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Install data.table with install.packages('data.table')")
}

suppressPackageStartupMessages(library(data.table))

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------

cancers <- c(
  "GBM", "LAML", "PAAD", "STAD", "TGCT",
  "THYM", "SKCM", "UVM", "LUAD", "LUSC"
)

input_dir <- file.path(
  getwd(),
  "data",
  "analysis_ready",
  "cancer_tables"
)

results_dir <- file.path(
  getwd(),
  "results",
  "04"
)

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

analysis_configuration <- data.table(
  cancer = cancers,
  input_file = file.path(
    input_dir,
    paste0(cancers, "_analysis_ready.tsv.gz")
  ),
  categorical_comparison = "CD43_high_vs_CD43_low",
  continuous_predictor = "SPN_log2_TPM_plus_1",
  primary_analysis_set = "all_samples",
  sensitivity_analysis_set = "CIBERSORTx_P_value_le_0.05",
  multiple_testing_method = "Benjamini-Hochberg"
)

fwrite(
  analysis_configuration,
  file.path(results_dir, "04_analysis_configuration.csv")
)

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

rank_biserial <- function(high, low) {
  n_high <- length(high)
  n_low <- length(low)
  ranks <- rank(c(high, low), ties.method = "average")
  u_high <- sum(ranks[seq_len(n_high)]) - n_high * (n_high + 1) / 2
  2 * u_high / (n_high * n_low) - 1
}

analyze_dataset <- function(dat, cancer, analysis_set) {
  p_index <- match("CIBERSORTx_P_value", names(dat))

  if (is.na(p_index) || p_index != 24L) {
    stop(cancer, ": could not identify the expected 22 LM22 columns.")
  }

  cell_types <- names(dat)[2:23]

  extreme <- dat[CD43_group %in% c("CD43_low", "CD43_high")]
  high <- extreme[CD43_group == "CD43_high"]
  low <- extreme[CD43_group == "CD43_low"]

  categorical <- rbindlist(lapply(cell_types, function(cell) {
    x <- high[[cell]]
    y <- low[[cell]]
    test <- suppressWarnings(wilcox.test(x, y, exact = FALSE))

    data.table(
      cancer = cancer,
      analysis_set = analysis_set,
      cell_type = cell,
      n_high = length(x),
      n_low = length(y),
      median_high = median(x),
      median_low = median(y),
      median_difference_high_minus_low = median(x) - median(y),
      mean_high = mean(x),
      mean_low = mean(y),
      rank_biserial = rank_biserial(x, y),
      wilcoxon_p = test$p.value
    )
  }))

  categorical[, wilcoxon_FDR := p.adjust(wilcoxon_p, method = "BH")]
  categorical[, significant_FDR_0_05 := wilcoxon_FDR <= 0.05]
  categorical[, direction := fifelse(
    rank_biserial > 0,
    "higher_in_CD43_high",
    fifelse(
      rank_biserial < 0,
      "lower_in_CD43_high",
      "no_direction"
    )
  )]

  continuous <- rbindlist(lapply(cell_types, function(cell) {
    test <- suppressWarnings(cor.test(
      dat$SPN_log2_TPM_plus_1,
      dat[[cell]],
      method = "spearman",
      exact = FALSE
    ))

    data.table(
      cancer = cancer,
      analysis_set = analysis_set,
      cell_type = cell,
      n = nrow(dat),
      spearman_rho = unname(test$estimate),
      spearman_p = test$p.value
    )
  }))

  continuous[, spearman_FDR := p.adjust(spearman_p, method = "BH")]
  continuous[, significant_FDR_0_05 := spearman_FDR <= 0.05]
  continuous[, direction := fifelse(
    spearman_rho > 0,
    "positive",
    fifelse(
      spearman_rho < 0,
      "negative",
      "no_direction"
    )
  )]

  list(
    categorical = categorical,
    continuous = continuous
  )
}

# ------------------------------------------------------------------------------
# Run cancer-specific analyses
# ------------------------------------------------------------------------------

categorical_results <- list()
continuous_results <- list()

for (cancer in cancers) {
  message("Analyzing ", cancer)

  input_file <- file.path(
    input_dir,
    paste0(cancer, "_analysis_ready.tsv.gz")
  )

  if (!file.exists(input_file)) {
    stop("Missing input: ", input_file)
  }

  dat <- fread(input_file, check.names = FALSE)

  primary <- analyze_dataset(
    dat,
    cancer,
    "primary_all_samples"
  )

  sensitivity <- analyze_dataset(
    dat[CIBERSORTx_pass_P_0_05 == TRUE],
    cancer,
    "sensitivity_P_le_0_05"
  )

  categorical_results[[cancer]] <- rbindlist(list(
    primary$categorical,
    sensitivity$categorical
  ))

  continuous_results[[cancer]] <- rbindlist(list(
    primary$continuous,
    sensitivity$continuous
  ))

  fwrite(
    categorical_results[[cancer]],
    file.path(
      results_dir,
      paste0("04_", cancer, "_CD43_high_vs_low.csv")
    )
  )

  fwrite(
    continuous_results[[cancer]],
    file.path(
      results_dir,
      paste0("04_", cancer, "_SPN_continuous_correlations.csv")
    )
  )
}

# ------------------------------------------------------------------------------
# Combine cancers and summarize cross-cancer consistency
# ------------------------------------------------------------------------------

categorical_all <- rbindlist(categorical_results)
continuous_all <- rbindlist(continuous_results)

fwrite(
  categorical_all,
  file.path(results_dir, "04_all_cancers_CD43_high_vs_low.csv")
)

fwrite(
  continuous_all,
  file.path(results_dir, "04_all_cancers_SPN_correlations.csv")
)

primary_cat <- categorical_all[
  analysis_set == "primary_all_samples"
]

primary_cor <- continuous_all[
  analysis_set == "primary_all_samples"
]

consistency <- primary_cor[
  ,
  .(
    cancers_tested = .N,
    positive_cancers = sum(spearman_rho > 0),
    negative_cancers = sum(spearman_rho < 0),
    significant_positive_cancers = sum(
      spearman_rho > 0 & spearman_FDR <= 0.05
    ),
    significant_negative_cancers = sum(
      spearman_rho < 0 & spearman_FDR <= 0.05
    ),
    median_rho = median(spearman_rho),
    mean_rho = mean(spearman_rho)
  ),
  by = cell_type
][order(-abs(median_rho))]

fwrite(
  consistency,
  file.path(results_dir, "04_cross_cancer_consistency_summary.csv")
)

rho_matrix <- dcast(
  primary_cor,
  cell_type ~ cancer,
  value.var = "spearman_rho"
)

effect_matrix <- dcast(
  primary_cat,
  cell_type ~ cancer,
  value.var = "rank_biserial"
)

fwrite(
  rho_matrix,
  file.path(results_dir, "04_cross_cancer_spearman_rho_matrix.csv")
)

fwrite(
  effect_matrix,
  file.path(results_dir, "04_cross_cancer_rank_biserial_matrix.csv")
)

# ------------------------------------------------------------------------------
# Complete run
# ------------------------------------------------------------------------------

message("Analysis completed. Results saved in: ", results_dir)
print(consistency)

writeLines(
  capture.output(sessionInfo()),
  file.path(results_dir, "04_sessionInfo.txt")