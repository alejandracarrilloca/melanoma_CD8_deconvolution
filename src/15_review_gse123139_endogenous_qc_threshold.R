#!/usr/bin/env Rscript

# ==============================================================================
# Step 15: Review endogenous QC thresholds after ERCC removal
# ==============================================================================
#
# Step 14 proved that 92 plate-specific rows are ERCC spike-ins. This script
# uses the saved endogenous QC metrics to reassess the UMI cutoff before the
# final SingleCellExperiment is rebuilt. It does not modify or filter data.
#
# Input:
#   data/analysis_ready/scRNA/GSE123139/
#     14_GSE123139_endogenous_cell_qc.tsv.gz
#
# Outputs:
#   results/15/
#     15_endogenous_threshold_retention.csv
#     15_failure_reasons_by_threshold.csv
#     15_endogenous_qc_quantiles.csv
#     15_figure_manifest.csv
#     figures/15A_endogenous_library_complexity.png
#     figures/15B_endogenous_umi_distributions.png
#     figures/15C_tumor_retention_by_umi_threshold.png
#
# Run from the project root:
#   Rscript src/15_review_gse123139_endogenous_qc_threshold.R
# ==============================================================================

options(stringsAsFactors = FALSE)

required_packages <- c("data.table", "ggplot2", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

project_root <- normalizePath(getwd(), mustWork = TRUE)
analysis_dir <- file.path(
  project_root, "data", "analysis_ready", "scRNA", "GSE123139"
)
results_dir <- file.path(project_root, "results", "15")
figures_dir <- file.path(results_dir, "figures")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

qc_file <- file.path(
  analysis_dir, "14_GSE123139_endogenous_cell_qc.tsv.gz"
)
if (!file.exists(qc_file)) stop("Missing Step 14 QC table: ", qc_file)

qc <- fread(qc_file, check.names = FALSE)
required_columns <- c(
  "cell_id", "source_group", "plate_configuration",
  "endogenous_n_counts", "endogenous_n_genes",
  "endogenous_mitochondrial_percent", "endogenous_high_library_outlier"
)
missing_columns <- setdiff(required_columns, names(qc))
if (length(missing_columns) > 0L) {
  stop("Step 14 QC table is missing: ", paste(missing_columns, collapse = ", "))
}

minimum_genes <- 200L
maximum_mitochondrial_percent <- 20
umi_thresholds <- c(500L, 600L, 700L, 800L, 900L, 1000L, 1250L, 1500L)

retention <- rbindlist(lapply(umi_thresholds, function(threshold) {
  qc[, {
    pass <- endogenous_n_counts >= threshold &
      endogenous_n_genes >= minimum_genes &
      !is.na(endogenous_mitochondrial_percent) &
      endogenous_mitochondrial_percent <= maximum_mitochondrial_percent &
      !endogenous_high_library_outlier

    .(
      minimum_umi = threshold,
      minimum_genes = minimum_genes,
      maximum_mitochondrial_percent = maximum_mitochondrial_percent,
      n_input_wells = .N,
      n_pass = sum(pass),
      percent_pass = 100 * mean(pass)
    )
  }, by = source_group]
}))

failure_reasons <- rbindlist(lapply(umi_thresholds, function(threshold) {
  qc[, {
    low_umi <- endogenous_n_counts < threshold
    low_genes <- endogenous_n_genes < minimum_genes
    high_mito <- is.na(endogenous_mitochondrial_percent) |
      endogenous_mitochondrial_percent > maximum_mitochondrial_percent
    high_library <- endogenous_high_library_outlier
    pass <- !low_umi & !low_genes & !high_mito & !high_library

    .(
      minimum_umi = threshold,
      n_pass = sum(pass),
      n_fail = sum(!pass),
      fail_low_umi = sum(low_umi),
      fail_low_genes = sum(low_genes),
      fail_high_mitochondrial = sum(high_mito),
      fail_high_library_outlier = sum(high_library),
      fail_only_low_umi = sum(
        low_umi & !low_genes & !high_mito & !high_library
      ),
      fail_only_low_genes = sum(
        !low_umi & low_genes & !high_mito & !high_library
      )
    )
  }, by = source_group]
}))

quantile_probabilities <- c(0, 0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 1)
quantiles <- qc[, {
  count_quantiles <- quantile(
    endogenous_n_counts, probs = quantile_probabilities, na.rm = TRUE
  )
  gene_quantiles <- quantile(
    endogenous_n_genes, probs = quantile_probabilities, na.rm = TRUE
  )
  rbind(
    data.table(
      metric = "endogenous_n_counts",
      probability = quantile_probabilities,
      value = as.numeric(count_quantiles)
    ),
    data.table(
      metric = "endogenous_n_genes",
      probability = quantile_probabilities,
      value = as.numeric(gene_quantiles)
    )
  )
}, by = .(source_group, plate_configuration)]

fwrite(retention, file.path(results_dir, "15_endogenous_threshold_retention.csv"))
fwrite(
  failure_reasons,
  file.path(results_dir, "15_failure_reasons_by_threshold.csv")
)
fwrite(quantiles, file.path(results_dir, "15_endogenous_qc_quantiles.csv"))

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

set.seed(123139)
plot_cells <- if (nrow(qc) > 40000L) qc[sample(.N, 40000L)] else copy(qc)

plot_15a <- ggplot(
  plot_cells,
  aes(x = log10(endogenous_n_counts + 1), y = endogenous_n_genes)
) +
  geom_point(
    aes(fill = source_group),
    shape = 21, size = 1.4, stroke = 0.2, color = gray_outline, alpha = 0.45
  ) +
  geom_vline(
    xintercept = log10(500 + 1),
    color = gray_outline,
    linetype = "dotted",
    linewidth = 0.7
  ) +
  geom_vline(
    xintercept = log10(800 + 1),
    color = lavender_dark,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  geom_vline(
    xintercept = log10(1000 + 1),
    color = "#222222",
    linetype = "solid",
    linewidth = 0.7
  ) +
  geom_hline(yintercept = minimum_genes, color = "#222222", linetype = "dashed") +
  scale_fill_manual(values = c("Tumor" = lavender, "PBMC" = gray_fill)) +
  facet_wrap(~source_group) +
  labs(
    title = "Endogenous library complexity after ERCC removal",
    subtitle = "Vertical guides mark 500, 800, and 1,000 endogenous UMIs; 40,000 wells plotted at most",
    x = "log10(endogenous UMIs + 1)",
    y = "Endogenous detected genes",
    fill = "Sample source"
  ) +
  theme_project

plot_15b <- ggplot(
  qc[endogenous_n_counts > 0],
  aes(x = log10(endogenous_n_counts + 1), fill = source_group, color = source_group)
) +
  geom_density(alpha = 0.30, linewidth = 0.7) +
  geom_vline(
    xintercept = log10(500 + 1),
    color = gray_outline,
    linetype = "dotted",
    linewidth = 0.7
  ) +
  geom_vline(
    xintercept = log10(800 + 1),
    color = lavender_dark,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  geom_vline(
    xintercept = log10(1000 + 1),
    color = "#222222",
    linetype = "solid",
    linewidth = 0.7
  ) +
  scale_fill_manual(values = c("Tumor" = lavender, "PBMC" = gray_fill)) +
  scale_color_manual(values = c("Tumor" = lavender_dark, "PBMC" = gray_outline)) +
  facet_wrap(~plate_configuration, scales = "free_y") +
  labs(
    title = "Endogenous UMI distributions across plate configurations",
    subtitle = "The ERCC-containing plates are assessed using ERCC-excluded counts",
    x = "log10(endogenous UMIs + 1)",
    y = "Density",
    fill = "Sample source",
    color = "Sample source"
  ) +
  theme_project

plot_15c <- ggplot(
  retention[source_group == "Tumor"],
  aes(x = minimum_umi, y = n_pass)
) +
  geom_line(color = lavender_dark, linewidth = 1) +
  geom_point(shape = 21, fill = lavender, color = lavender_dark, size = 3) +
  scale_x_continuous(breaks = umi_thresholds) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Tumor-cell retention across endogenous UMI thresholds",
    subtitle = "All thresholds also require >=200 genes, <=20% mitochondrial UMIs, and no endogenous high-library outlier flag",
    x = "Minimum endogenous UMIs",
    y = "Retained tumor cells"
  ) +
  theme_project +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

figure_paths <- file.path(
  figures_dir,
  c(
    "15A_endogenous_library_complexity.png",
    "15B_endogenous_umi_distributions.png",
    "15C_tumor_retention_by_umi_threshold.png"
  )
)

ggsave(figure_paths[1L], plot_15a, width = 11, height = 6.5, dpi = 300, bg = "white")
ggsave(figure_paths[2L], plot_15b, width = 11, height = 6.5, dpi = 300, bg = "white")
ggsave(figure_paths[3L], plot_15c, width = 9, height = 6, dpi = 300, bg = "white")

fwrite(
  data.table(
    figure = c("15A", "15B", "15C"),
    file = basename(figure_paths),
    purpose = c(
      "Inspect corrected library complexity and candidate cutoffs",
      "Compare corrected endogenous UMI distributions by plate configuration",
      "Quantify tumor-cell retention sensitivity to the UMI cutoff"
    )
  ),
  file.path(results_dir, "15_figure_manifest.csv")
)

writeLines(capture.output(sessionInfo()), file.path(results_dir, "15_sessionInfo.txt"))

message("Step 15 threshold review complete.")
message("Tumor retention:")
print(retention[source_group == "Tumor"])
message("PBMC retention:")
print(retention[source_group == "PBMC"])