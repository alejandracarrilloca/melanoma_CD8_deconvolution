#!/usr/bin/env Rscript

# Step 3: Test associations between bulk SPN/CD43 expression and complementary
# immune-infiltration measurements in TCGA-SKCM and TCGA-UVM.
#
# Analyses are performed separately by cancer. The script reports:
#   1. Spearman correlations using continuous log2(SPN TPM + 1)
#   2. CD43-high versus CD43-low Wilcoxon comparisons
#   3. Standardized linear models adjusted for tumor purity
#   4. CD8 models adjusted for PTPRC to distinguish CD8-specific association
#      from general leukocyte infiltration
#
# Input:
#   data/analysis_ready/skcm_uvm_immune_scores/
#     08_tcga_analysis_with_immune_scores.tsv.gz
#
# Outputs:
#   results/09/
#     09_spn_immune_spearman_results.csv
#     09_cd43_high_vs_low_immune_comparisons.csv
#     09_spn_adjusted_immune_models.csv
#     09_spn_immune_association_qc.csv
#     09_sessionInfo.txt
#     figures/09_<CANCER>_SPN_vs_<MEASURE>.png
#     figures/09_<CANCER>_SPN_immune_association_scatterplots.png
#
# Run from the project root with:
#   Rscript src/09_test_spn_immune_associations.R

options(stringsAsFactors = FALSE)

plot_style_version <- "gray_lavender_directionality_style_v3"

input_file <- file.path(
  "data",
  "analysis_ready",
  "skcm_uvm_immune_scores",
  "08_tcga_analysis_with_immune_scores.tsv.gz"
)
results_directory <- file.path("results", "09")
cancers <- c("SKCM", "UVM")

figure_directory <- file.path(results_directory, "figures")
invisible(lapply(
  c(results_directory, figure_directory),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop(
    "The ggplot2 package is required for the project-standard scatterplots. ",
    "Install it with install.packages('ggplot2') or the appropriate server package manager."
  )
}

read_tsv <- function(path) {
  if (!file.exists(path)) stop("Missing input file: ", path)
  read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA", "N/A", "Not Reported", "not reported", "[Not Available]")
  )
}

safe_scale <- function(x) {
  x <- as.numeric(x)
  if (sum(!is.na(x)) < 3L || sd(x, na.rm = TRUE) == 0) return(rep(NA_real_, length(x)))
  as.numeric(scale(x))
}

immune_measures <- c(
  "PTPRC_log2_TPM_plus_1",
  "T_cell_marker_score",
  "CD8_marker_score",
  "T cells CD8",
  "relative_total_T_cell_fraction",
  "relative_total_lymphocyte_fraction",
  "relative_CD8_within_T_cells",
  "tumor_purity_primary"
)

measure_labels <- c(
  PTPRC_log2_TPM_plus_1 = "PTPRC expression",
  T_cell_marker_score = "T-cell marker score",
  CD8_marker_score = "CD8 marker score",
  `T cells CD8` = "CIBERSORTx CD8 fraction",
  relative_total_T_cell_fraction = "Total T-cell fraction",
  relative_total_lymphocyte_fraction = "Total lymphocyte fraction",
  relative_CD8_within_T_cells = "CD8 within T cells",
  tumor_purity_primary = "Tumor purity"
)

required_columns <- c(
  "sample_id", "cancer", "CD43_group", "SPN_TPM",
  "tumor_purity_primary", immune_measures
)

cat("Reading:", input_file, "\n")
x <- read_tsv(input_file)

missing_columns <- setdiff(required_columns, names(x))
if (length(missing_columns) > 0L) {
  stop("Missing required column(s): ", paste(missing_columns, collapse = ", "))
}
if (anyDuplicated(x$sample_id)) stop("Duplicate sample_id values in input table")
if (!all(x$cancer %in% cancers)) stop("Unexpected cancer value in input table")

numeric_columns <- unique(c("SPN_TPM", immune_measures))
for (variable in numeric_columns) {
  original <- x[[variable]]
  converted <- suppressWarnings(as.numeric(as.character(original)))
  if (any(is.na(converted) & !is.na(original))) {
    stop(variable, " contains non-numeric values")
  }
  x[[variable]] <- converted
}

x$SPN_log2_TPM_plus_1_analysis <- log2(x$SPN_TPM + 1)

# -----------------------------------------------------------------------------
# 1. Spearman correlations with continuous SPN expression
# -----------------------------------------------------------------------------

spearman_results <- do.call(rbind, lapply(cancers, function(cancer) {
  y <- x[x$cancer == cancer, , drop = FALSE]

  do.call(rbind, lapply(immune_measures, function(outcome) {
    keep <- complete.cases(y[, c("SPN_log2_TPM_plus_1_analysis", outcome), drop = FALSE])
    n_complete <- sum(keep)

    if (n_complete < 3L ||
        length(unique(y$SPN_log2_TPM_plus_1_analysis[keep])) < 2L ||
        length(unique(y[[outcome]][keep])) < 2L) {
      rho <- NA_real_
      p_value <- NA_real_
      status <- "insufficient_variation"
    } else {
      test <- suppressWarnings(cor.test(
        y$SPN_log2_TPM_plus_1_analysis[keep],
        y[[outcome]][keep],
        method = "spearman",
        exact = FALSE
      ))
      rho <- unname(test$estimate)
      p_value <- test$p.value
      status <- "ok"
    }

    data.frame(
      cancer = cancer,
      outcome = outcome,
      outcome_label = unname(measure_labels[outcome]),
      n = n_complete,
      spearman_rho = rho,
      p_value = p_value,
      status = status,
      stringsAsFactors = FALSE
    )
  }))
}))

spearman_results$FDR <- ave(
  spearman_results$p_value,
  spearman_results$cancer,
  FUN = function(p) p.adjust(p, method = "BH")
)

# -----------------------------------------------------------------------------
# 2. CD43-high versus CD43-low comparisons
# -----------------------------------------------------------------------------

group_results <- do.call(rbind, lapply(cancers, function(cancer) {
  y <- x[x$cancer == cancer & x$CD43_group %in% c("CD43_high", "CD43_low"), , drop = FALSE]

  do.call(rbind, lapply(immune_measures, function(outcome) {
    high <- y[[outcome]][y$CD43_group == "CD43_high"]
    low <- y[[outcome]][y$CD43_group == "CD43_low"]
    high <- high[!is.na(high)]
    low <- low[!is.na(low)]

    if (length(high) < 2L || length(low) < 2L) {
      wilcoxon_p <- NA_real_
      rank_biserial <- NA_real_
      status <- "insufficient_data"
    } else {
      test <- suppressWarnings(wilcox.test(high, low, exact = FALSE))
      # For the two-sample R implementation, the reported W is the
      # Mann-Whitney U statistic for the first group (CD43_high).
      mann_whitney_u <- unname(test$statistic)
      rank_biserial <- 2 * mann_whitney_u / (length(high) * length(low)) - 1
      wilcoxon_p <- test$p.value
      status <- "ok"
    }

    data.frame(
      cancer = cancer,
      outcome = outcome,
      outcome_label = unname(measure_labels[outcome]),
      n_high = length(high),
      n_low = length(low),
      median_high = if (length(high) == 0L) NA_real_ else median(high),
      median_low = if (length(low) == 0L) NA_real_ else median(low),
      median_difference_high_minus_low = if (length(high) == 0L || length(low) == 0L) {
        NA_real_
      } else {
        median(high) - median(low)
      },
      rank_biserial = rank_biserial,
      wilcoxon_p = wilcoxon_p,
      status = status,
      stringsAsFactors = FALSE
    )
  }))
}))

group_results$FDR <- ave(
  group_results$wilcoxon_p,
  group_results$cancer,
  FUN = function(p) p.adjust(p, method = "BH")
)

# -----------------------------------------------------------------------------
# 3. Adjusted standardized linear models
# -----------------------------------------------------------------------------

model_definitions <- data.frame(
  outcome = c(
    "PTPRC_log2_TPM_plus_1",
    "T_cell_marker_score",
    "CD8_marker_score",
    "T cells CD8",
    "relative_total_T_cell_fraction",
    "relative_total_lymphocyte_fraction",
    "relative_CD8_within_T_cells",
    "CD8_marker_score",
    "T cells CD8",
    "relative_CD8_within_T_cells"
  ),
  covariate = c(
    rep("tumor_purity_primary", 7),
    rep("PTPRC_log2_TPM_plus_1", 3)
  ),
  model_question = c(
    rep("SPN association adjusted for tumor purity", 7),
    rep("CD8 association adjusted for general leukocyte expression", 3)
  ),
  stringsAsFactors = FALSE
)

adjusted_results <- do.call(rbind, lapply(cancers, function(cancer) {
  y <- x[x$cancer == cancer, , drop = FALSE]

  do.call(rbind, lapply(seq_len(nrow(model_definitions)), function(i) {
    outcome <- model_definitions$outcome[i]
    covariate <- model_definitions$covariate[i]
    variables <- c("SPN_log2_TPM_plus_1_analysis", outcome, covariate)
    keep <- complete.cases(y[, variables, drop = FALSE])
    d <- y[keep, variables, drop = FALSE]

    if (nrow(d) < 10L || any(vapply(d, function(z) length(unique(z)) < 2L, logical(1)))) {
      return(data.frame(
        cancer = cancer,
        outcome = outcome,
        outcome_label = unname(measure_labels[outcome]),
        covariate = covariate,
        model_question = model_definitions$model_question[i],
        n = nrow(d),
        standardized_SPN_beta = NA_real_,
        standard_error = NA_real_,
        t_value = NA_real_,
        p_value = NA_real_,
        adjusted_R_squared = NA_real_,
        predictor_covariate_correlation = NA_real_,
        status = "insufficient_data_or_variation",
        stringsAsFactors = FALSE
      ))
    }

    model_data <- data.frame(
      outcome_z = safe_scale(d[[outcome]]),
      SPN_z = safe_scale(d$SPN_log2_TPM_plus_1_analysis),
      covariate_z = safe_scale(d[[covariate]])
    )

    fit <- lm(outcome_z ~ SPN_z + covariate_z, data = model_data)
    fit_summary <- summary(fit)
    coefficient <- fit_summary$coefficients["SPN_z", ]

    data.frame(
      cancer = cancer,
      outcome = outcome,
      outcome_label = unname(measure_labels[outcome]),
      covariate = covariate,
      model_question = model_definitions$model_question[i],
      n = nrow(model_data),
      standardized_SPN_beta = unname(coefficient["Estimate"]),
      standard_error = unname(coefficient["Std. Error"]),
      t_value = unname(coefficient["t value"]),
      p_value = unname(coefficient["Pr(>|t|)"]),
      adjusted_R_squared = fit_summary$adj.r.squared,
      predictor_covariate_correlation = cor(
        model_data$SPN_z,
        model_data$covariate_z,
        method = "pearson"
      ),
      status = "ok",
      stringsAsFactors = FALSE
    )
  }))
}))

adjusted_results$FDR <- ave(
  adjusted_results$p_value,
  adjusted_results$cancer,
  FUN = function(p) p.adjust(p, method = "BH")
)

# -----------------------------------------------------------------------------
# 4. QC summaries and plots
# -----------------------------------------------------------------------------

qc <- do.call(rbind, lapply(cancers, function(cancer) {
  y <- x[x$cancer == cancer, , drop = FALSE]
  data.frame(
    cancer = cancer,
    n_samples = nrow(y),
    n_unique_samples = length(unique(y$sample_id)),
    n_cd43_low = sum(y$CD43_group == "CD43_low", na.rm = TRUE),
    n_intermediate = sum(y$CD43_group == "intermediate", na.rm = TRUE),
    n_cd43_high = sum(y$CD43_group == "CD43_high", na.rm = TRUE),
    n_complete_all_tested_measures = sum(complete.cases(y[, c(
      "SPN_log2_TPM_plus_1_analysis", immune_measures
    ), drop = FALSE])),
    stringsAsFactors = FALSE
  )
}))

safe_filename <- function(value) {
  value <- gsub("[^A-Za-z0-9]+", "_", value)
  value <- gsub("^_+|_+$", "", value)
  tolower(value)
}

plot_palette <- c(
  `All other tumors` = "#E6E6E6",
  `CD43-high tumors` = "#C599FF"
)

combined_panel_measures <- c(
  "PTPRC_log2_TPM_plus_1",
  "T_cell_marker_score",
  "CD8_marker_score",
  "T cells CD8",
  "relative_total_lymphocyte_fraction",
  "tumor_purity_primary"
)

for (cancer in cancers) {
  y <- x[x$cancer == cancer, , drop = FALSE]
  combined_panels <- list()

  for (outcome in immune_measures) {
    keep <- complete.cases(y[, c(
      "SPN_log2_TPM_plus_1_analysis", outcome, "CD43_group"
    ), drop = FALSE])
    plot_data <- y[keep, , drop = FALSE]

    plot_data$point_class <- factor(
      ifelse(
        plot_data$CD43_group == "CD43_high",
        "CD43-high tumors",
        "All other tumors"
      ),
      levels = c("All other tumors", "CD43-high tumors")
    )

    # Draw highlighted CD43-high samples last so they remain visible.
    plot_data <- plot_data[order(plot_data$point_class), , drop = FALSE]

    result_row <- spearman_results[
      spearman_results$cancer == cancer & spearman_results$outcome == outcome,
      , drop = FALSE
    ]

    statistics_label <- sprintf(
      "TCGA-%s, n = %d\nSpearman rho = %.2f, BH FDR = %.3g",
      cancer,
      nrow(plot_data),
      result_row$spearman_rho,
      result_row$FDR
    )

    figure <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = SPN_log2_TPM_plus_1_analysis,
        y = .data[[outcome]],
        fill = point_class
      )
    ) +
      ggplot2::geom_hline(
        yintercept = 0,
        color = "#BDBDBD",
        linetype = "dashed",
        linewidth = 0.65
      ) +
      ggplot2::geom_smooth(
        ggplot2::aes(group = 1, linetype = "Linear fit"),
        method = "lm",
        formula = y ~ x,
        se = TRUE,
        color = "#222222",
        fill = "#E6E6E6",
        linewidth = 0.85,
        alpha = 0.60,
        inherit.aes = TRUE
      ) +
      ggplot2::geom_point(
        shape = 21,
        size = 3.0,
        color = "#222222",
        stroke = 0.30,
        alpha = 0.84
      ) +
      ggplot2::scale_fill_manual(
        values = plot_palette,
        drop = FALSE,
        name = NULL
      ) +
      ggplot2::scale_linetype_manual(
        values = c(`Linear fit` = "dashed"),
        name = NULL
      ) +
      ggplot2::annotate(
        "text",
        x = -Inf,
        y = Inf,
        label = statistics_label,
        hjust = -0.04,
        vjust = 1.15,
        size = 3.7,
        fontface = "bold",
        color = "#222222",
        lineheight = 1.08
      ) +
      ggplot2::expand_limits(y = 0) +
      ggplot2::labs(
        x = "log2(SPN TPM + 1)",
        y = unname(measure_labels[outcome])
      ) +
      ggplot2::theme_classic(base_size = 11) +
      ggplot2::theme(
        aspect.ratio = 1,
        plot.background = ggplot2::element_rect(fill = "white", color = NA),
        panel.background = ggplot2::element_rect(fill = "white", color = NA),
        axis.title = ggplot2::element_text(
          size = 14,
          face = "bold",
          color = "#222222"
        ),
        axis.text = ggplot2::element_text(size = 10.5, color = "#222222"),
        axis.line = ggplot2::element_line(color = "#222222", linewidth = 0.9),
        axis.ticks = ggplot2::element_line(color = "#222222", linewidth = 0.75),
        axis.ticks.length = grid::unit(5, "pt"),
        panel.grid = ggplot2::element_blank(),
        legend.position = c(0.035, 0.035),
        legend.justification = c(0, 0),
        legend.direction = "vertical",
        legend.text = ggplot2::element_text(size = 9.2, color = "#222222"),
        legend.key = ggplot2::element_blank(),
        legend.background = ggplot2::element_rect(fill = "white", color = NA),
        legend.spacing.y = grid::unit(0, "pt"),
        plot.margin = ggplot2::margin(12, 14, 12, 12)
      ) +
      ggplot2::guides(
        fill = ggplot2::guide_legend(
          order = 1,
          override.aes = list(shape = 21, size = 2.7, alpha = 1)
        ),
        linetype = ggplot2::guide_legend(
          order = 2,
          override.aes = list(color = "#222222", linewidth = 0.8)
        )
      )

    output_path <- file.path(
      figure_directory,
      paste0(
        "09_", cancer, "_SPN_vs_", safe_filename(outcome), ".png"
      )
    )

    ggplot2::ggsave(
      filename = output_path,
      plot = figure,
      width = 7.2,
      height = 7.2,
      units = "in",
      dpi = 600,
      bg = "white"
    )

    if (outcome %in% combined_panel_measures) {
      panel_figure <- figure

      # Retain one internal legend, as in the reference figure, and avoid
      # repeating it in every panel of the combined 2 x 3 layout.
      if (outcome != combined_panel_measures[1]) {
        panel_figure <- panel_figure + ggplot2::theme(legend.position = "none")
      }

      combined_panels[[outcome]] <- panel_figure
    }
  }

  combined_output <- file.path(
    figure_directory,
    paste0("09_", cancer, "_SPN_immune_association_scatterplots.png")
  )

  grDevices::png(
    filename = combined_output,
    width = 16,
    height = 10.7,
    units = "in",
    res = 400,
    bg = "white"
  )
  grid::grid.newpage()
  panel_layout <- grid::grid.layout(nrow = 2, ncol = 3)
  grid::pushViewport(grid::viewport(layout = panel_layout))

  for (panel_index in seq_along(combined_panel_measures)) {
    panel_name <- combined_panel_measures[panel_index]
    panel_row <- ceiling(panel_index / 3)
    panel_column <- ((panel_index - 1) %% 3) + 1

    print(
      combined_panels[[panel_name]],
      vp = grid::viewport(
        layout.pos.row = panel_row,
        layout.pos.col = panel_column
      ),
      newpage = FALSE
    )
  }
  grid::popViewport()
  grDevices::dev.off()

  cat("Styled combined PNG:", combined_output, "\n")
}

write.csv(
  spearman_results,
  file.path(results_directory, "09_spn_immune_spearman_results.csv"),
  row.names = FALSE
)
write.csv(
  group_results,
  file.path(results_directory, "09_cd43_high_vs_low_immune_comparisons.csv"),
  row.names = FALSE
)
write.csv(
  adjusted_results,
  file.path(results_directory, "09_spn_adjusted_immune_models.csv"),
  row.names = FALSE
)
write.csv(
  qc,
  file.path(results_directory, "09_spn_immune_association_qc.csv"),
  row.names = FALSE
)
capture.output(
  sessionInfo(),
  file = file.path(results_directory, "09_sessionInfo.txt")
)

cat("\nStep 3 complete.\n")
cat("Plot style:", plot_style_version, "\n")
cat("\nQC summary:\n")
print(qc, row.names = FALSE)
cat("\nSpearman results:\n")
print(spearman_results, row.names = FALSE)
cat("\nCD43-high versus CD43-low results:\n")
print(group_results, row.names = FALSE)
cat("\nAdjusted SPN models:\n")
print(adjusted_results, row.names = FALSE)
cat("\nPositive coefficients indicate that higher SPN is associated with a higher outcome value.\n")
cat("For tumor purity, a negative association indicates higher SPN in less-pure tumors.\n")