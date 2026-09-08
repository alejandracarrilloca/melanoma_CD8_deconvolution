#!/usr/bin/env Rscript

# ==============================================================================
# Step 14: Audit GSE123139 QC after removal of plate-specific features
# ==============================================================================
#
# Step 13 showed that 58 plates contain 55,765 genes while 146 plates contain
# 55,857 genes. The 92 additional features may be ERCC spike-ins. Because Step
# 12 calculated library size across every row in each plate, those features may
# have influenced the >=1,000-UMI filter and the plate-specific high-library
# outlier flag. This script quantifies that influence before normalization.
#
# This is an audit only: it does not alter the Step 13 SingleCellExperiment.
# It calculates endogenous/common-gene QC metrics for all 78,336 wells, compares
# the original and corrected retention decisions, and saves a revised QC table
# that can be used to rebuild the filtered matrix if necessary.
#
# Inputs:
#   data/scRNA/GSE123139/raw_counts/GSM*_*.txt.gz
#   data/analysis_ready/scRNA/GSE123139/
#     12_GSE123139_cell_qc.tsv.gz
#     13_GSE123139_common_gene_symbols.tsv.gz
#
# Outputs:
#   data/analysis_ready/scRNA/GSE123139/
#     14_GSE123139_endogenous_cell_qc.tsv.gz
#   results/14/
#     14_excluded_feature_inventory.csv
#     14_qc_decision_summary.csv
#     14_retention_comparison_by_source.csv
#     14_retention_comparison_by_plate_configuration.csv
#     14_endogenous_qc_threshold_sensitivity.csv
#     14_cells_with_changed_qc_decision.tsv.gz
#     14_dataset_qc.csv
#     14_figure_manifest.csv
#     14_sessionInfo.txt
#     figures/*.png
#
# Run from the project root:
#   Rscript src/14_audit_gse123139_endogenous_qc.R
# ==============================================================================

options(stringsAsFactors = FALSE)

required_packages <- c("data.table", "ggplot2")
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

# ------------------------------------------------------------------------------
# Paths and QC definition
# ------------------------------------------------------------------------------

project_root <- normalizePath(getwd(), mustWork = TRUE)
raw_dir <- file.path(project_root, "data", "scRNA", "GSE123139", "raw_counts")
analysis_dir <- file.path(
  project_root, "data", "analysis_ready", "scRNA", "GSE123139"
)
results_dir <- file.path(project_root, "results", "14")
figures_dir <- file.path(results_dir, "figures")

dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

cell_qc_file <- file.path(analysis_dir, "12_GSE123139_cell_qc.tsv.gz")
common_gene_file <- file.path(
  analysis_dir, "13_GSE123139_common_gene_symbols.tsv.gz"
)
output_qc_file <- file.path(
  analysis_dir, "14_GSE123139_endogenous_cell_qc.tsv.gz"
)

minimum_umi <- 1000L
minimum_genes <- 200L
maximum_mitochondrial_percent <- 20
high_library_nmads <- 4

for (path in c(cell_qc_file, common_gene_file)) {
  if (!file.exists(path)) stop("Missing required input: ", path)
}

count_files <- sort(list.files(
  raw_dir,
  pattern = "^GSM[0-9]+_.*\\.txt\\.gz$",
  full.names = TRUE
))

if (length(count_files) == 0L) {
  stop("No GSE123139 count files found in: ", raw_dir)
}

qc <- fread(cell_qc_file, check.names = FALSE)
common_genes <- fread(common_gene_file)[[1L]]
common_genes <- as.character(common_genes)

required_columns <- c(
  "cell_id", "count_file", "gsm_accession", "source_group", "subject_id",
  "specimen_id", "n_counts", "n_genes", "mitochondrial_counts",
  "mitochondrial_percent", "high_library_outlier"
)
missing_columns <- setdiff(required_columns, names(qc))

if (length(missing_columns) > 0L) {
  stop("Step 12 QC table is missing: ", paste(missing_columns, collapse = ", "))
}
if (anyDuplicated(qc$cell_id)) stop("Duplicated cell IDs occur in Step 12 QC.")
if (anyDuplicated(common_genes)) stop("Duplicated common genes are present.")

# ------------------------------------------------------------------------------
# Helpers
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

safe_percent <- function(numerator, denominator) {
  fifelse(denominator > 0, 100 * numerator / denominator, NA_real_)
}

# ------------------------------------------------------------------------------
# Count features excluded by the common-gene harmonization
# ------------------------------------------------------------------------------

correction_list <- vector("list", length(count_files))
excluded_inventory_list <- vector("list", length(count_files))
plate_feature_qc_list <- vector("list", length(count_files))

for (i in seq_along(count_files)) {
  path <- count_files[[i]]
  file_name <- basename(path)
  message("Auditing plate ", i, "/", length(count_files), ": ", file_name)

  cell_ids <- read_cell_header(path)
  count_table <- fread(
    path,
    skip = 1L,
    header = FALSE,
    col.names = c("gene_symbol", cell_ids),
    check.names = FALSE,
    showProgress = FALSE
  )

  if (ncol(count_table) != length(cell_ids) + 1L) {
    stop("Count-file width does not match header: ", file_name)
  }

  genes <- as.character(count_table[[1L]])
  excluded_rows <- which(!genes %chin% common_genes)
  excluded_genes <- genes[excluded_rows]

  if (length(excluded_rows) == 0L) {
    excluded_counts <- numeric(length(cell_ids))
    excluded_detected <- integer(length(cell_ids))
  } else {
    excluded_matrix <- as.matrix(
      count_table[excluded_rows, -1L, with = FALSE]
    )
    storage.mode(excluded_matrix) <- "numeric"
    excluded_counts <- as.numeric(colSums(excluded_matrix))
    excluded_detected <- as.integer(colSums(excluded_matrix > 0))

    excluded_inventory_list[[i]] <- data.table(
      count_file = file_name,
      gene_symbol = excluded_genes
    )
    rm(excluded_matrix)
  }

  correction_list[[i]] <- data.table(
    cell_id = cell_ids,
    count_file = file_name,
    excluded_feature_counts = excluded_counts,
    excluded_features_detected = excluded_detected
  )

  plate_feature_qc_list[[i]] <- data.table(
    count_file = file_name,
    n_input_features = length(genes),
    n_common_features = sum(genes %chin% common_genes),
    n_excluded_features = length(excluded_genes),
    excluded_features_all_ercc = if (length(excluded_genes) == 0L) {
      NA
    } else {
      all(grepl("^ERCC-", excluded_genes, ignore.case = TRUE))
    }
  )

  rm(count_table)
  if (i %% 10L == 0L) invisible(gc())
}

corrections <- rbindlist(correction_list)
plate_feature_qc <- rbindlist(plate_feature_qc_list)

if (nrow(corrections) != nrow(qc) || anyDuplicated(corrections$cell_id)) {
  stop("Excluded-feature correction table does not match the Step 12 cells.")
}

if (length(Filter(Negate(is.null), excluded_inventory_list)) > 0L) {
  excluded_inventory_long <- rbindlist(excluded_inventory_list, fill = TRUE)
  excluded_inventory <- excluded_inventory_long[, .(
    n_plates_present = uniqueN(count_file),
    ercc_name_pattern = grepl("^ERCC-", gene_symbol, ignore.case = TRUE)
  ), by = gene_symbol][order(gene_symbol)]
} else {
  excluded_inventory <- data.table(
    gene_symbol = character(),
    n_plates_present = integer(),
    ercc_name_pattern = logical()
  )
}

fwrite(
  excluded_inventory,
  file.path(results_dir, "14_excluded_feature_inventory.csv")
)
fwrite(
  plate_feature_qc,
  file.path(results_dir, "14_plate_feature_configuration.csv")
)

# ------------------------------------------------------------------------------
# Recalculate endogenous/common-gene QC for every well
# ------------------------------------------------------------------------------

qc <- merge(
  qc,
  corrections,
  by = c("cell_id", "count_file"),
  all.x = TRUE,
  all.y = FALSE,
  sort = FALSE
)

if (nrow(qc) != nrow(corrections) || anyNA(qc$excluded_feature_counts)) {
  stop("Failed to attach excluded-feature counts to every cell.")
}

qc <- merge(qc, plate_feature_qc, by = "count_file", all.x = TRUE, sort = FALSE)

qc[, endogenous_n_counts := n_counts - excluded_feature_counts]
qc[, endogenous_n_genes := n_genes - excluded_features_detected]

if (any(qc$endogenous_n_counts < 0) || any(qc$endogenous_n_genes < 0)) {
  stop("Corrected endogenous QC metrics became negative.")
}

qc[, excluded_count_fraction := fifelse(
  n_counts > 0, excluded_feature_counts / n_counts, NA_real_
)]
qc[, endogenous_mitochondrial_percent := safe_percent(
  mitochondrial_counts, endogenous_n_counts
)]
qc[, endogenous_log10_counts_plus_1 := log10(endogenous_n_counts + 1)]

qc[, endogenous_high_library_cutoff := {
  center <- median(endogenous_log10_counts_plus_1, na.rm = TRUE)
  spread <- mad(endogenous_log10_counts_plus_1, center = center, na.rm = TRUE)
  10^(center + high_library_nmads * spread) - 1
}, by = gsm_accession]

qc[, endogenous_high_library_outlier :=
  endogenous_n_counts > endogenous_high_library_cutoff
]

qc[, original_final_qc_pass :=
  n_counts >= minimum_umi &
  n_genes >= minimum_genes &
  !is.na(mitochondrial_percent) &
  mitochondrial_percent <= maximum_mitochondrial_percent &
  !high_library_outlier
]

qc[, endogenous_final_qc_pass :=
  endogenous_n_counts >= minimum_umi &
  endogenous_n_genes >= minimum_genes &
  !is.na(endogenous_mitochondrial_percent) &
  endogenous_mitochondrial_percent <= maximum_mitochondrial_percent &
  !endogenous_high_library_outlier
]

qc[, qc_decision := fcase(
  original_final_qc_pass & endogenous_final_qc_pass, "Pass in both",
  original_final_qc_pass & !endogenous_final_qc_pass, "Lost after correction",
  !original_final_qc_pass & endogenous_final_qc_pass, "Gained after correction",
  default = "Fail in both"
)]

qc[, plate_configuration := ifelse(
  n_excluded_features == 0L,
  paste0(n_input_features, " common features only"),
  paste0(n_input_features, " features; ", n_excluded_features, " excluded")
)]

fwrite(qc, output_qc_file, sep = "\t")
fwrite(
  qc[original_final_qc_pass != endogenous_final_qc_pass],
  file.path(results_dir, "14_cells_with_changed_qc_decision.tsv.gz"),
  sep = "\t"
)

# ------------------------------------------------------------------------------
# Summary tables
# ------------------------------------------------------------------------------

decision_summary <- qc[, .(
  n_cells = .N,
  percent_cells = 100 * .N / nrow(qc)
), by = qc_decision][order(qc_decision)]

retention_by_source <- qc[, .(
  n_cell_columns = .N,
  original_pass = sum(original_final_qc_pass),
  endogenous_pass = sum(endogenous_final_qc_pass),
  net_change = sum(endogenous_final_qc_pass) - sum(original_final_qc_pass),
  lost_after_correction = sum(
    original_final_qc_pass & !endogenous_final_qc_pass
  ),
  gained_after_correction = sum(
    !original_final_qc_pass & endogenous_final_qc_pass
  )
), by = source_group][order(source_group)]

retention_by_configuration <- qc[, .(
  n_cell_columns = .N,
  median_excluded_count_fraction = median(
    excluded_count_fraction, na.rm = TRUE
  ),
  maximum_excluded_count_fraction = max(
    excluded_count_fraction, na.rm = TRUE
  ),
  original_pass = sum(original_final_qc_pass),
  endogenous_pass = sum(endogenous_final_qc_pass),
  net_change = sum(endogenous_final_qc_pass) - sum(original_final_qc_pass),
  lost_after_correction = sum(
    original_final_qc_pass & !endogenous_final_qc_pass
  ),
  gained_after_correction = sum(
    !original_final_qc_pass & endogenous_final_qc_pass
  )
), by = .(source_group, plate_configuration)][
  order(source_group, plate_configuration)
]

threshold_grid <- CJ(
  minimum_umi = c(500L, 800L, 1000L, 1500L),
  minimum_genes = c(100L, 200L),
  maximum_mitochondrial_percent = c(10, 20, 30)
)

threshold_sensitivity <- rbindlist(lapply(
  seq_len(nrow(threshold_grid)),
  function(i) {
    threshold <- threshold_grid[i]
    qc[, .(
      n_cell_columns = .N,
      n_pass = sum(
        endogenous_n_counts >= threshold$minimum_umi &
        endogenous_n_genes >= threshold$minimum_genes &
        !is.na(endogenous_mitochondrial_percent) &
        endogenous_mitochondrial_percent <=
          threshold$maximum_mitochondrial_percent &
        !endogenous_high_library_outlier
      ),
      percent_pass = 100 * mean(
        endogenous_n_counts >= threshold$minimum_umi &
        endogenous_n_genes >= threshold$minimum_genes &
        !is.na(endogenous_mitochondrial_percent) &
        endogenous_mitochondrial_percent <=
          threshold$maximum_mitochondrial_percent &
        !endogenous_high_library_outlier
      ),
      minimum_umi = threshold$minimum_umi,
      minimum_genes = threshold$minimum_genes,
      maximum_mitochondrial_percent =
        threshold$maximum_mitochondrial_percent
    ), by = source_group]
  }
))

setcolorder(
  threshold_sensitivity,
  c(
    "source_group", "minimum_umi", "minimum_genes",
    "maximum_mitochondrial_percent", "n_cell_columns", "n_pass",
    "percent_pass"
  )
)

dataset_qc <- data.table(
  metric = c(
    "plate_files", "cell_columns", "common_features",
    "unique_excluded_features", "excluded_features_matching_ercc_pattern",
    "plates_without_excluded_features", "plates_with_excluded_features",
    "original_final_qc_pass", "endogenous_final_qc_pass",
    "lost_after_correction", "gained_after_correction", "net_retention_change",
    "minimum_endogenous_count_fraction_among_original_pass",
    "median_endogenous_count_fraction_among_original_pass",
    "maximum_excluded_count_fraction_all_cells"
  ),
  value = as.character(c(
    length(count_files), nrow(qc), length(common_genes),
    nrow(excluded_inventory), sum(excluded_inventory$ercc_name_pattern),
    sum(plate_feature_qc$n_excluded_features == 0L),
    sum(plate_feature_qc$n_excluded_features > 0L),
    sum(qc$original_final_qc_pass), sum(qc$endogenous_final_qc_pass),
    sum(qc$original_final_qc_pass & !qc$endogenous_final_qc_pass),
    sum(!qc$original_final_qc_pass & qc$endogenous_final_qc_pass),
    sum(qc$endogenous_final_qc_pass) - sum(qc$original_final_qc_pass),
    min(
      qc$endogenous_n_counts[qc$original_final_qc_pass] /
        qc$n_counts[qc$original_final_qc_pass],
      na.rm = TRUE
    ),
    median(
      qc$endogenous_n_counts[qc$original_final_qc_pass] /
        qc$n_counts[qc$original_final_qc_pass],
      na.rm = TRUE
    ),
    max(qc$excluded_count_fraction, na.rm = TRUE)
  ))
)

fwrite(decision_summary, file.path(results_dir, "14_qc_decision_summary.csv"))
fwrite(
  retention_by_source,
  file.path(results_dir, "14_retention_comparison_by_source.csv")
)
fwrite(
  retention_by_configuration,
  file.path(results_dir, "14_retention_comparison_by_plate_configuration.csv")
)
fwrite(
  threshold_sensitivity,
  file.path(results_dir, "14_endogenous_qc_threshold_sensitivity.csv")
)
fwrite(dataset_qc, file.path(results_dir, "14_dataset_qc.csv"))

# ------------------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------------------

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

plot_14a <- ggplot(
  plot_cells,
  aes(x = log10(n_counts + 1), y = log10(endogenous_n_counts + 1))
) +
  geom_abline(slope = 1, intercept = 0, color = "#222222", linetype = "dashed") +
  geom_point(
    aes(fill = n_excluded_features > 0L),
    shape = 21, size = 1.5, stroke = 0.25, color = gray_outline, alpha = 0.55
  ) +
  scale_fill_manual(
    values = c("FALSE" = gray_fill, "TRUE" = lavender),
    labels = c("FALSE" = "Common features only", "TRUE" = "Extra features removed")
  ) +
  facet_wrap(~source_group) +
  labs(
    title = "Original and common-gene library sizes",
    subtitle = "Points below the diagonal contained counts from plate-specific features",
    x = "log10(original UMIs + 1)",
    y = "log10(common-gene UMIs + 1)",
    fill = NULL
  ) +
  theme_project

plot_14b <- ggplot(
  qc[n_excluded_features > 0L & is.finite(excluded_count_fraction)],
  aes(x = 100 * excluded_count_fraction, fill = source_group, color = source_group)
) +
  geom_density(alpha = 0.30, linewidth = 0.7) +
  scale_fill_manual(values = c("Tumor" = lavender, "PBMC" = gray_fill)) +
  scale_color_manual(values = c("Tumor" = lavender_dark, "PBMC" = gray_outline)) +
  labs(
    title = "Contribution of plate-specific features to library size",
    subtitle = "Only plates containing features outside the common gene set are shown",
    x = "Counts from excluded features (%)",
    y = "Density",
    fill = "Sample source",
    color = "Sample source"
  ) +
  theme_project

decision_plot_data <- qc[, .N, by = .(source_group, qc_decision)]
plot_14c <- ggplot(
  decision_plot_data,
  aes(x = qc_decision, y = N, fill = source_group)
) +
  geom_col(position = "dodge", color = gray_outline, linewidth = 0.25) +
  scale_fill_manual(values = c("Tumor" = lavender, "PBMC" = gray_fill)) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Effect of endogenous-count correction on QC retention",
    subtitle = "Final thresholds: >=1,000 UMIs, >=200 genes, <=20% mitochondrial UMIs",
    x = NULL,
    y = "Number of wells",
    fill = "Sample source"
  ) +
  theme_project +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

figure_paths <- file.path(
  figures_dir,
  c(
    "14A_original_vs_endogenous_library_size.png",
    "14B_excluded_count_fraction.png",
    "14C_qc_decision_changes.png"
  )
)

ggsave(figure_paths[1L], plot_14a, width = 11, height = 6.5, dpi = 300, bg = "white")
ggsave(figure_paths[2L], plot_14b, width = 9, height = 6, dpi = 300, bg = "white")
ggsave(figure_paths[3L], plot_14c, width = 10, height = 6, dpi = 300, bg = "white")

fwrite(
  data.table(
    figure = c("14A", "14B", "14C"),
    file = basename(figure_paths),
    purpose = c(
      "Compare original and common-gene library sizes",
      "Quantify counts contributed by excluded features",
      "Show cell-level changes in final QC decisions"
    )
  ),
  file.path(results_dir, "14_figure_manifest.csv")
)

writeLines(capture.output(sessionInfo()), file.path(results_dir, "14_sessionInfo.txt"))

message("Step 14 audit complete.")
message("Excluded features: ", nrow(excluded_inventory))
message(
  "Excluded features matching ^ERCC-: ",
  sum(excluded_inventory$ercc_name_pattern), "/", nrow(excluded_inventory)
)
message("Original QC pass: ", sum(qc$original_final_qc_pass))
message("Endogenous QC pass: ", sum(qc$endogenous_final_qc_pass))
message("Decision summary:")
print(decision_summary)
message("Retention by source:")
print(retention_by_source)