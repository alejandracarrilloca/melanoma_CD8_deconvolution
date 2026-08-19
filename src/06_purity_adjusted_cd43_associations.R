#!/usr/bin/env Rscript

# ==============================================================================
# Description:  Joins published TCGA tumor-purity estimates to the analysis-
#               ready CD43/CIBERSORTx dataset and evaluates whether tumor purity
#               explains the observed CD43-associated immune composition. The
#               analysis is restricted to the eight solid-tumor cohorts; LAML
#               and THYM remain context-specific because conventional purity
#               and LM22 infiltration have different biological meanings there.
#
#               Four analyses are performed for every available purity method:
#               (1) CD43-high versus CD43-low purity comparison, (2) continuous
#               SPN-purity Spearman correlation, (3) purity-adjusted continuous
#               SPN-cell-fraction regression, and (4) purity-adjusted CD43-high
#               versus CD43-low regression. Models use arcsine-square-root-
#               transformed relative fractions. Unadjusted models are fit on
#               the identical complete-case samples used by adjusted models.
#
# Inputs:       data/analysis_ready/all_cancers_analysis_ready.tsv.gz
#               TCGAbiolinks::Tumor.purity
#
# Outputs:      data/purity_analysis/
#                 all_samples_with_purity.tsv.gz
#                 reference/TCGA_mastercalls.abs_tables_JSedit.fixed.txt
#
#               results/06/
#                 06_purity_matching_summary.csv
#                 06_purity_method_coverage_summary.csv
#                 06_primary_purity_method.csv
#                 06_CD43_group_purity_tests.csv
#                 06_SPN_purity_correlations.csv
#                 06_continuous_models_unadjusted_vs_purity_adjusted.csv
#                 06_categorical_models_unadjusted_vs_purity_adjusted.csv
#                 06_sessionInfo.txt
#                 figures/06A_primary_purity_by_CD43_group.png
#                 figures/06B_primary_purity_adjusted_continuous_effects.png
#                 figures/06C_primary_purity_adjusted_effect_heatmap.png
#
# Usage:        Rscript src/06_purity_adjusted_cd43_associations.R
#
# Dependencies: TCGAbiolinks, data.table, ggplot2, scales
# ==============================================================================

options(stringsAsFactors = FALSE)

required_packages <- c("TCGAbiolinks", "data.table", "ggplot2", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Install missing package(s): ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

solid_cancers <- c("GBM", "PAAD", "STAD", "TGCT", "SKCM", "UVM", "LUAD", "LUSC")
legacy_purity_methods <- c("CPE", "ABSOLUTE", "LUMP", "ESTIMATE", "IHC")
pancan_method <- "PANCAN_ABSOLUTE"
purity_methods <- c(pancan_method, legacy_purity_methods)

cd43_low_fill <- "#79e875"
cd43_low_edge <- "#3bb53d"
cd43_high_fill <- "#FC4E92"
cd43_high_edge <- "#C8326D"
negative_fill <- "#927CEB"
negative_edge <- "#6552BD"
neutral_fill <- "#C9CDD2"
neutral_edge <- "#7B8086"

input_file <- file.path(
  getwd(), "data", "analysis_ready",
  "all_cancers_analysis_ready.tsv.gz")

data_output_dir <- file.path(getwd(), "data", "purity_analysis")
results_dir <- file.path(getwd(), "results", "06")
figure_dir <- file.path(results_dir, "figures")
reference_dir <- file.path(data_output_dir, "reference")
absolute_file <- file.path(reference_dir, "TCGA_mastercalls.abs_tables_JSedit.fixed.txt")
absolute_url <- paste0(
  "https://api.gdc.cancer.gov/data/",
  "4f277128-f793-4354-a13d-30cc7fe9f6b5")

invisible(lapply(
  c(data_output_dir, reference_dir, results_dir, figure_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

if (!file.exists(input_file)) stop("Missing input: ", input_file)

normalize_tcga_sample <- function(x) {
  x <- gsub("[._]", "-", toupper(as.character(x)))
  position <- regexpr(
    "TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}[A-Z]",
    x,
    perl = TRUE
  )
  result <- rep(NA_character_, length(x))
  found <- position > 0L
  match_length <- attr(position, "match.length")
  result[found] <- substr(
    x[found],
    position[found],
    position[found] + match_length[found] - 1L
  )
  result
}

rank_biserial <- function(high, low) {
  n_high <- length(high)
  n_low <- length(low)
  ranks <- rank(c(high, low), ties.method = "average")
  u_high <- sum(ranks[seq_len(n_high)]) - n_high * (n_high + 1) / 2
  2 * u_high / (n_high * n_low) - 1
}

extract_model <- function(dat, outcome, predictor, purity_column, adjusted) {
  model_data <- data.table(
    y = asin(sqrt(pmin(pmax(dat[[outcome]], 0), 1))),
    predictor = dat[[predictor]],
    purity = dat[[purity_column]]
  )
  model_data <- model_data[complete.cases(model_data)]
  if (nrow(model_data) < 20L || sd(model_data$purity) == 0) {
    return(data.table(n = nrow(model_data), estimate = NA_real_, SE = NA_real_, p_value = NA_real_))
  }

  if (predictor == "SPN_log2_TPM_plus_1") {
    model_data[, predictor_z := as.numeric(scale(predictor))]
    model_data[, purity_z := as.numeric(scale(purity))]
    fit <- if (adjusted) {
      lm(y ~ predictor_z + purity_z, data = model_data)
    } else {
      lm(y ~ predictor_z, data = model_data)
    }
    coefficient_name <- "predictor_z"
  } else {
    model_data[, predictor := factor(
      predictor,
      levels = c("CD43_low", "CD43_high")
    )]
    model_data[, purity_z := as.numeric(scale(purity))]
    fit <- if (adjusted) {
      lm(y ~ predictor + purity_z, data = model_data)
    } else {
      lm(y ~ predictor, data = model_data)
    }
    coefficient_name <- "predictorCD43_high"
  }

  coefs <- summary(fit)$coefficients
  if (!coefficient_name %in% rownames(coefs)) {
    return(data.table(n = nrow(model_data), estimate = NA_real_, SE = NA_real_, p_value = NA_real_))
  }
  data.table(
    n = nrow(model_data),
    estimate = unname(coefs[coefficient_name, "Estimate"]),
    SE = unname(coefs[coefficient_name, "Std. Error"]),
    p_value = unname(coefs[coefficient_name, "Pr(>|t|)"])
  )
}

# Load the packaged purity table without relying on attachment to the search path.
purity_environment <- new.env(parent = globalenv())
data("Tumor.purity", package = "TCGAbiolinks", envir = purity_environment)
if (!exists("Tumor.purity", envir = purity_environment, inherits = FALSE)) {
  stop("TCGAbiolinks::Tumor.purity could not be loaded.")
}
purity <- as.data.table(get("Tumor.purity", envir = purity_environment))

required_purity_columns <- c("Sample.ID", "Cancer.type", legacy_purity_methods)
missing_purity_columns <- setdiff(required_purity_columns, names(purity))
if (length(missing_purity_columns) > 0L) {
  stop("Tumor.purity is missing: ", paste(missing_purity_columns, collapse = ", "))
}

# Some TCGAbiolinks/Bioconductor builds store the published purity columns as
# factors or character strings. Convert through character explicitly so that
# factor level codes are never mistaken for purity values. Non-numeric missing
# value markers are converted to NA.
parse_purity_value <- function(x) {
  raw <- trimws(as.character(x))
  missing_marker <- is.na(raw) | raw %in% c(
    "", "NA", "N/A", "--", "[Not Available]", "Not Available"
  )
  uses_percent <- grepl("%", raw, fixed = TRUE)
  raw <- gsub(",", ".", raw, fixed = TRUE)
  number_position <- regexpr(
    "[-+]?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)",
    raw,
    perl = TRUE
  )
  parsed <- rep(NA_real_, length(raw))
  has_number <- number_position > 0L & !missing_marker
  parsed[has_number] <- suppressWarnings(as.numeric(
    regmatches(raw, number_position)[has_number]
  ))
  parsed[uses_percent & !is.na(parsed)] <- parsed[uses_percent & !is.na(parsed)] / 100
  parsed
}

purity[, (legacy_purity_methods) := lapply(.SD, parse_purity_value),
       .SDcols = legacy_purity_methods]

numeric_counts <- vapply(
  purity[, ..legacy_purity_methods],
  function(x) sum(is.finite(x)),
  integer(1)
)
message(
  "Usable purity values after numeric conversion: ",
  paste(names(numeric_counts), numeric_counts, sep = "=", collapse = "; ")
)
if (sum(numeric_counts) == 0L) {
  stop("Purity conversion produced no numeric values; inspect TCGAbiolinks::Tumor.purity.")
}

purity[, sample_barcode := normalize_tcga_sample(Sample.ID)]
purity <- unique(purity, by = "sample_barcode")

# Download the open-access PanCanAtlas ABSOLUTE purity/ploidy mastercalls once
# and cache them for reproducible reruns. The UUID is the official GDC file ID.
if (!file.exists(absolute_file) || file.info(absolute_file)$size == 0L) {
  message("Downloading PanCanAtlas ABSOLUTE mastercalls...")
  download_ok <- tryCatch({
    suppressWarnings(download.file(
      absolute_url,
      absolute_file,
      mode = "wb",
      method = "libcurl",
      quiet = FALSE
    ))
    file.exists(absolute_file) && file.info(absolute_file)$size > 0L
  }, error = function(e) {
    message("Automatic download failed: ", conditionMessage(e))
    FALSE
  })
  if (!download_ok) {
    stop(
      "PanCanAtlas ABSOLUTE mastercalls could not be downloaded. Run:\n",
      "curl -L --retry 3 '", absolute_url, "' -o '", absolute_file, "'\n",
      "and rerun this script."
    )
  }
}

absolute_calls <- fread(absolute_file, check.names = FALSE)
normalized_absolute_names <- tolower(gsub("[^a-z0-9]+", "_", names(absolute_calls)))
sample_candidates <- which(normalized_absolute_names %in% c(
  "sample", "sample_id", "sampleid", "tumor_sample_barcode"
))
if (length(sample_candidates) == 0L) {
  sample_candidates <- grep("sample", normalized_absolute_names)
}
purity_candidates <- which(normalized_absolute_names == "purity")
if (length(sample_candidates) == 0L || length(purity_candidates) == 0L) {
  stop(
    "Could not identify sample and purity columns in ", absolute_file,
    ". Columns found: ", paste(names(absolute_calls), collapse = ", ")
  )
}
absolute_join <- absolute_calls[, .(
  sample_barcode = normalize_tcga_sample(get(names(absolute_calls)[sample_candidates[1L]])),
  PANCAN_ABSOLUTE = parse_purity_value(get(names(absolute_calls)[purity_candidates[1L]]))
)]
absolute_join <- absolute_join[
  !is.na(sample_barcode) & is.finite(PANCAN_ABSOLUTE) &
    PANCAN_ABSOLUTE >= 0 & PANCAN_ABSOLUTE <= 1
]
absolute_join <- absolute_join[, .(
  PANCAN_ABSOLUTE = median(PANCAN_ABSOLUTE)
), by = sample_barcode]
message("Usable PanCanAtlas ABSOLUTE values: ", nrow(absolute_join))

analysis_data <- fread(input_file, check.names = FALSE)
analysis_data[, sample_barcode := normalize_tcga_sample(sample_id)]
analysis_data <- analysis_data[cancer %in% solid_cancers]

p_index <- match("CIBERSORTx_P_value", names(analysis_data))
if (is.na(p_index) || p_index != 24L) {
  stop("Could not identify the expected 22 LM22 columns in the analysis-ready table.")
}
cell_types <- names(analysis_data)[2:23]

purity_join <- purity[, c("sample_barcode", "Cancer.type", legacy_purity_methods), with = FALSE]
setnames(purity_join, "Cancer.type", "purity_cancer")
joined <- merge(analysis_data, purity_join, by = "sample_barcode", all.x = TRUE, sort = FALSE)
joined <- merge(joined, absolute_join, by = "sample_barcode", all.x = TRUE, sort = FALSE)

fwrite(
  joined,
  file.path(data_output_dir, "all_samples_with_purity.tsv.gz"),
  sep = "\t",
  compress = "gzip"
)

matching_summary <- rbindlist(lapply(purity_methods, function(method) {
  joined[, .(
    total_samples = .N,
    matched_samples = sum(!is.na(get(method))),
    matched_percent = 100 * mean(!is.na(get(method)))
  ), by = cancer][, purity_method := method]
}))
setcolorder(matching_summary, c("cancer", "purity_method", "total_samples", "matched_samples", "matched_percent"))
fwrite(
  matching_summary,
  file.path(results_dir, "06_purity_matching_summary.csv")
)

# Select a primary method only when it has usable coverage within every cancer.
# Total pan-cancer matches alone are not sufficient because a method concentrated
# in a few cohorts would create a biased apparent primary analysis.
method_coverage <- matching_summary[, .(
  total_matched_samples = sum(matched_samples),
  cancers_with_any = sum(matched_samples > 0L),
  cancers_with_at_least_20 = sum(matched_samples >= 20L),
  cancers_with_at_least_50_percent = sum(matched_percent >= 50),
  minimum_matched_samples = min(matched_samples),
  minimum_matched_percent = min(matched_percent)
), by = purity_method]
fwrite(
  method_coverage,
  file.path(results_dir, "06_purity_method_coverage_summary.csv")
)

primary_candidates <- method_coverage[
  cancers_with_at_least_20 == length(solid_cancers)
][order(
  -as.integer(purity_method == pancan_method),
  -cancers_with_at_least_50_percent,
  -minimum_matched_percent,
  -total_matched_samples
)]
if (nrow(primary_candidates) == 0L) {
  stop(
    "No purity method has at least 20 matched samples in every solid-tumor cohort. ",
    "A pan-cancer primary purity analysis would be biased. Review ",
    file.path(results_dir, "06_purity_matching_summary.csv"), " and ",
    file.path(results_dir, "06_purity_method_coverage_summary.csv"), "."
  )
}
primary_purity_method <- primary_candidates$purity_method[1L]
message(
  "Primary purity method selected after per-cancer coverage checks: ",
  primary_purity_method,
  " (", primary_candidates$total_matched_samples[1L], " matched samples; minimum ",
  sprintf("%.1f", primary_candidates$minimum_matched_percent[1L]),
  "% within a cancer)"
)
fwrite(
  data.table(
    primary_purity_method = primary_purity_method,
    matched_samples = primary_candidates$total_matched_samples[1L],
    cancers_with_at_least_20 = primary_candidates$cancers_with_at_least_20[1L],
    minimum_matched_percent = primary_candidates$minimum_matched_percent[1L]
  ),
  file.path(results_dir, "06_primary_purity_method.csv")
)

# Purity differences between the extreme CD43 groups.
group_tests <- rbindlist(lapply(purity_methods, function(method) {
  rbindlist(lapply(solid_cancers, function(cancer_name) {
    dat <- joined[cancer == cancer_name & CD43_group %in% c("CD43_low", "CD43_high")]
    high <- dat[CD43_group == "CD43_high"][[method]]
    low <- dat[CD43_group == "CD43_low"][[method]]
    high <- high[is.finite(high)]
    low <- low[is.finite(low)]
    if (length(high) < 5L || length(low) < 5L) return(NULL)
    test <- suppressWarnings(wilcox.test(high, low, exact = FALSE))
    data.table(
      cancer = cancer_name,
      purity_method = method,
      n_high = length(high),
      n_low = length(low),
      median_high = median(high),
      median_low = median(low),
      median_difference_high_minus_low = median(high) - median(low),
      rank_biserial = rank_biserial(high, low),
      p_value = test$p.value
    )
  }), fill = TRUE)
}), fill = TRUE)
group_tests[, FDR := p.adjust(p_value, method = "BH"), by = purity_method]
fwrite(
  group_tests,
  file.path(results_dir, "06_CD43_group_purity_tests.csv")
)

# Continuous SPN-purity correlations.
purity_correlations <- rbindlist(lapply(purity_methods, function(method) {
  rbindlist(lapply(solid_cancers, function(cancer_name) {
    dat <- joined[cancer == cancer_name]
    keep <- is.finite(dat$SPN_log2_TPM_plus_1) & is.finite(dat[[method]])
    if (sum(keep) < 10L) return(NULL)
    test <- suppressWarnings(cor.test(
      dat$SPN_log2_TPM_plus_1[keep], dat[[method]][keep],
      method = "spearman", exact = FALSE
    ))
    data.table(
      cancer = cancer_name,
      purity_method = method,
      n = sum(keep),
      spearman_rho = unname(test$estimate),
      p_value = test$p.value
    )
  }), fill = TRUE)
}), fill = TRUE)
purity_correlations[, FDR := p.adjust(p_value, method = "BH"), by = purity_method]
fwrite(
  purity_correlations,
  file.path(results_dir, "06_SPN_purity_correlations.csv")
)

# Purity-adjusted continuous and categorical cell-fraction models.
fit_grid <- function(predictor, extreme_only) {
  results <- rbindlist(lapply(purity_methods, function(method) {
    rbindlist(lapply(solid_cancers, function(cancer_name) {
      dat <- joined[cancer == cancer_name]
      if (extreme_only) dat <- dat[CD43_group %in% c("CD43_low", "CD43_high")]
      rbindlist(lapply(cell_types, function(cell) {
        unadjusted <- extract_model(dat, cell, predictor, method, FALSE)
        adjusted <- extract_model(dat, cell, predictor, method, TRUE)
        data.table(
          cancer = cancer_name,
          purity_method = method,
          cell_type = cell,
          n = adjusted$n,
          unadjusted_estimate = unadjusted$estimate,
          unadjusted_SE = unadjusted$SE,
          unadjusted_p = unadjusted$p_value,
          adjusted_estimate = adjusted$estimate,
          adjusted_SE = adjusted$SE,
          adjusted_p = adjusted$p_value
        )
      }))
    }))
  }), fill = TRUE)

  results[, unadjusted_FDR := p.adjust(unadjusted_p, method = "BH"),
          by = .(cancer, purity_method)]
  results[, adjusted_FDR := p.adjust(adjusted_p, method = "BH"),
          by = .(cancer, purity_method)]
  results[, attenuation_percent := fifelse(
    is.na(unadjusted_estimate) | unadjusted_estimate == 0,
    NA_real_,
    100 * (1 - abs(adjusted_estimate) / abs(unadjusted_estimate))
  )]
  results[, direction_retained := sign(unadjusted_estimate) == sign(adjusted_estimate)]
  results
}

continuous_models <- fit_grid("SPN_log2_TPM_plus_1", FALSE)
categorical_models <- fit_grid("CD43_group", TRUE)
fwrite(
  continuous_models,
  file.path(
    results_dir,
    "06_continuous_models_unadjusted_vs_purity_adjusted.csv"
  )
)
fwrite(
  categorical_models,
  file.path(
    results_dir,
    "06_categorical_models_unadjusted_vs_purity_adjusted.csv"
  )
)

# Primary purity measure by CD43 group.
primary_plot_data <- joined[
  CD43_group %in% c("CD43_low", "CD43_high") &
    !is.na(get(primary_purity_method))
]
primary_plot_data[, purity_value := get(primary_purity_method)]
primary_plot_data[, cancer := factor(cancer, levels = solid_cancers)]
primary_plot_data[, CD43_group := factor(
  CD43_group,
  levels = c("CD43_low", "CD43_high"),
  labels = c("CD43 low", "CD43 high")
)]

p1 <- ggplot(
  primary_plot_data,
  aes(CD43_group, purity_value, fill = CD43_group, colour = CD43_group)
) +
  geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.38) +
  geom_jitter(shape = 21, width = 0.12, size = 1.1, stroke = 0.25,
              alpha = 0.42) +
  facet_wrap(~ cancer, ncol = 4) +
  scale_fill_manual(
    values = c(
      "CD43 low" = alpha(cd43_low_fill, 0.42),
      "CD43 high" = alpha(cd43_high_fill, 0.42)
    ),
    guide = "none"
  ) +
  scale_colour_manual(
    values = c("CD43 low" = cd43_low_edge, "CD43 high" = cd43_high_edge),
    guide = "none"
  ) +
  scale_y_continuous(labels = label_percent(accuracy = 1)) +
  labs(
    title = "Tumor purity in CD43-low and CD43-high tumors",
    subtitle = paste0(primary_purity_method, " purity across solid-tumor cohorts"),
    x = NULL, y = paste0(primary_purity_method, " tumor purity")
  ) +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 13),
        strip.background = element_blank(), strip.text = element_text(face = "bold"),
        axis.text = element_text(colour = "grey30"), axis.ticks = element_blank())
ggsave(file.path(figure_dir, "06A_primary_purity_by_CD43_group.png"), p1,
       width = 12, height = 7, dpi = 600, bg = "white")

# Direct unadjusted-versus-adjusted effect comparison for the primary method.
primary_effects <- continuous_models[
  purity_method == primary_purity_method & !is.na(adjusted_estimate)
]
if (nrow(primary_effects) == 0L) {
  stop("No adjusted models could be fitted for primary method ", primary_purity_method, ".")
}
primary_effects[, effect_direction := factor(
  fifelse(
    adjusted_estimate > 0, "Positive",
    fifelse(adjusted_estimate < 0, "Negative", "Zero")
  ),
  levels = c("Negative", "Zero", "Positive")
)]
primary_effects[, adjusted_significant := adjusted_FDR <= 0.05]

p2 <- ggplot(primary_effects, aes(unadjusted_estimate, adjusted_estimate)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey70", linewidth = 0.5) +
  geom_hline(yintercept = 0, colour = "grey90", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey90", linewidth = 0.3) +
  geom_point(
    aes(
      fill = effect_direction,
      colour = effect_direction,
      alpha = adjusted_significant
    ),
    shape = 21, size = 2.2, stroke = 0.4
  ) +
  facet_wrap(~ cancer, ncol = 4, scales = "free") +
  scale_fill_manual(
    values = c(
      "Negative" = alpha(negative_fill, 0.60),
      "Zero" = alpha(neutral_fill, 0.60),
      "Positive" = alpha(cd43_high_fill, 0.60)
    ),
    name = "Adjusted effect"
  ) +
  scale_colour_manual(
    values = c(
      "Negative" = negative_edge,
      "Zero" = neutral_edge,
      "Positive" = cd43_high_edge
    ),
    name = "Adjusted effect"
  ) +
  scale_alpha_manual(
    values = c(`TRUE` = 0.95, `FALSE` = 0.35),
    labels = c(`TRUE` = "Yes", `FALSE` = "No"),
    name = "Adjusted FDR <= 0.05"
  ) +
  labs(title = "Effect of purity adjustment on SPN-cell associations",
       x = "Unadjusted standardized SPN coefficient",
       y = paste0(primary_purity_method, "-adjusted standardized SPN coefficient")) +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 13),
        strip.background = element_blank(), strip.text = element_text(face = "bold"),
        axis.text = element_text(colour = "grey30"), axis.ticks = element_blank())
ggsave(file.path(figure_dir, "06B_primary_purity_adjusted_continuous_effects.png"), p2,
       width = 12, height = 8, dpi = 600, bg = "white")

# Contiguous adjusted-effect heatmap.
heat <- copy(primary_effects)
cell_order <- heat[, .(median_effect = median(adjusted_estimate)),
                   by = cell_type][order(median_effect)]$cell_type
heat[, cancer := factor(cancer, levels = solid_cancers)]
heat[, cell_type := factor(cell_type, levels = cell_order)]
effect_limit <- max(abs(heat$adjusted_estimate), na.rm = TRUE)
p3 <- ggplot(heat, aes(cancer, cell_type, fill = adjusted_estimate)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_point(data = heat[adjusted_FDR <= 0.05], shape = 8, size = 1.4,
             colour = "black") +
  scale_fill_gradient2(low = negative_fill, mid = "#F7F7F7", high = cd43_high_fill,
                       midpoint = 0, limits = c(-effect_limit, effect_limit),
                       name = paste0(primary_purity_method, "-adjusted\nSPN effect")) +
  labs(title = "CD43-immune associations after tumor-purity adjustment",
       subtitle = "Continuous SPN models using all tumors; * indicates BH-FDR <= 0.05",
       x = NULL, y = NULL) +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 13),
        axis.text.x = element_text(angle = 45, hjust = 1, colour = "grey30"),
        axis.text.y = element_text(size = 8, colour = "grey30"),
        axis.ticks = element_blank(), axis.line = element_blank())
ggsave(file.path(figure_dir, "06C_primary_purity_adjusted_effect_heatmap.png"), p3,
       width = 10, height = 9, dpi = 600, bg = "white")

message("Purity analysis completed: ", results_dir)
print(matching_summary)
writeLines(
  capture.output(sessionInfo()),
  file.path(results_dir, "06_sessionInfo.txt")
)
