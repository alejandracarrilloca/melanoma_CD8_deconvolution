#!/usr/bin/env Rscript

# ==============================================================================
# Step 10: Quantify the cellular source of SPN/CD43 in melanoma scRNA-seq
# ==============================================================================
#
# Dataset:
#   GSE115978 (human melanoma single-cell RNA-seq)
#
# Purpose:
#   1. Validate cell and gene identifiers across annotations, counts, and TPM.
#   2. Quantify SPN expression and detection across annotated cell types.
#   3. Identify malignant SPN-positive/PTPRC-negative cells.
#   4. Confirm that those cells retain melanoma-lineage marker expression.
#   5. Create patient-by-cell-type summaries for non-cell-level follow-up.
#
# Inputs:
#   data/scRNA/GSE115978/GSE115978_cell.annotations.csv.gz
#   data/scRNA/GSE115978/GSE115978_counts.csv.gz
#   data/scRNA/GSE115978/GSE115978_tpm.csv.gz
#
# Reusable output:
#   data/analysis_ready/scRNA/GSE115978/
#     10_GSE115978_selected_marker_cell_table.tsv.gz
#
# Results:
#   results/10/*.csv
#   results/10/*.txt
#   results/10/figures/*.png
#
# Notes:
#   - Raw-count detection is used for detected/not-detected classifications.
#   - TPM is used for expression magnitude and log2(TPM + 1) plots.
#   - Cell-level results are descriptive because cells from one patient are not
#     statistically independent. Patient-by-cell-type summaries are therefore
#     also saved for subsequent analyses.
#   - SPN-positive/PTPRC-negative malignant cells support tumor-cell SPN
#     expression, but do not alone prove it because dropout, ambient RNA, and
#     annotation uncertainty remain possible.
#
# Run from the project root:
#   Rscript src/10_gse115978_spn_cellular_source.R
# ==============================================================================

options(stringsAsFactors = FALSE)

required_packages <- c("data.table", "ggplot2", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install the missing package(s): ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

# ------------------------------------------------------------------------------
# File locations
# ------------------------------------------------------------------------------

input_dir <- file.path(getwd(), "data", "scRNA", "GSE115978")
results_dir <- file.path(getwd(), "results", "10")
figure_dir <- file.path(results_dir, "figures")
analysis_dir <- file.path(
  getwd(), "data", "analysis_ready", "scRNA", "GSE115978"
)

annotation_file <- file.path(
  input_dir, "GSE115978_cell.annotations.csv.gz"
)
count_file <- file.path(input_dir, "GSE115978_counts.csv.gz")
tpm_file <- file.path(input_dir, "GSE115978_tpm.csv.gz")

invisible(lapply(
  c(results_dir, figure_dir, analysis_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

required_files <- c(annotation_file, count_file, tpm_file)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
  stop(
    "Missing required GSE115978 input file(s): ",
    paste(missing_files, collapse = ", ")
  )
}

# ------------------------------------------------------------------------------
# Marker definitions
# ------------------------------------------------------------------------------

marker_sets <- list(
  source = c("SPN", "PTPRC"),
  t_cell = c("CD3D", "CD3E", "TRAC"),
  cd8 = c("CD8A", "CD8B"),
  cytotoxic = c("NKG7", "GNLY", "GZMB", "PRF1"),
  melanoma = c("MLANA", "PMEL", "SOX10", "MITF", "TYR"),
  myeloid = c("LST1", "CD68"),
  b_cell = c("MS4A1", "CD79A"),
  fibroblast = c("COL1A1", "COL1A2"),
  endothelial = c("PECAM1", "VWF")
)

target_genes <- unique(unlist(marker_sets, use.names = FALSE))

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

normalize_gene_field <- function(x) {
  x <- sub("^\\ufeff", "", x)
  x <- sub('^"', "", x)
  sub('"$', "", x)
}

read_matrix_header <- function(path) {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection))

  header_line <- readLines(connection, n = 1L, warn = FALSE)
  if (length(header_line) != 1L) {
    stop("Could not read a header from: ", path)
  }

  header <- strsplit(header_line, ",", fixed = TRUE)[[1L]]
  header <- gsub('^"|"$', "", header)

  if (length(header) < 2L) {
    stop("Matrix header has fewer than two columns: ", path)
  }

  header
}

read_selected_gene_rows <- function(path, genes, chunk_size = 500L) {
  message("Scanning selected marker rows: ", path)

  connection <- gzfile(path, open = "rt")
  on.exit(close(connection))

  header_line <- readLines(connection, n = 1L, warn = FALSE)
  if (length(header_line) != 1L) {
    stop("Could not read matrix header: ", path)
  }

  selected_lines <- character()

  repeat {
    block <- readLines(
      connection,
      n = chunk_size,
      warn = FALSE
    )

    if (length(block) == 0L) {
      break
    }

    first_field <- sub(",.*$", "", block)
    first_field <- normalize_gene_field(first_field)
    keep <- first_field %chin% genes

    if (any(keep)) {
      selected_lines <- c(selected_lines, block[keep])
    }
  }

  if (length(selected_lines) == 0L) {
    stop("None of the requested marker genes were found in: ", path)
  }

  selected <- fread(
    text = paste(c(header_line, selected_lines), collapse = "\n"),
    sep = ",",
    header = TRUE,
    check.names = FALSE,
    showProgress = FALSE
  )

  setnames(selected, 1L, "gene")
  selected[, gene := normalize_gene_field(as.character(gene))]

  duplicated_genes <- selected[duplicated(gene) | duplicated(gene, fromLast = TRUE),
                               unique(gene)]
  if (length(duplicated_genes) > 0L) {
    stop(
      "Duplicate requested gene rows found in ", basename(path), ": ",
      paste(duplicated_genes, collapse = ", ")
    )
  }

  selected
}

selected_rows_to_cells <- function(selected, genes, suffix) {
  matrix_cells <- names(selected)[-1L]
  expression_matrix <- as.matrix(selected[, -1L, with = FALSE])
  storage.mode(expression_matrix) <- "double"
  rownames(expression_matrix) <- selected$gene

  output <- data.table(cells = matrix_cells)

  for (gene in genes) {
    output_name <- paste0(gene, suffix)

    if (gene %in% rownames(expression_matrix)) {
      output[[output_name]] <- as.numeric(expression_matrix[gene, ])
    } else {
      output[[output_name]] <- NA_real_
    }
  }

  output
}

available_columns <- function(dat, genes, suffix) {
  candidates <- paste0(genes, suffix)
  candidates[
    candidates %in% names(dat) &
      vapply(candidates, function(column) {
        !all(is.na(dat[[column]]))
      }, logical(1))
  ]
}

mean_log2_tpm_score <- function(dat, genes) {
  columns <- available_columns(dat, genes, "_TPM")

  if (length(columns) == 0L) {
    return(rep(NA_real_, nrow(dat)))
  }

  rowMeans(log2(as.matrix(dat[, ..columns]) + 1), na.rm = TRUE)
}

any_count_detected <- function(dat, genes) {
  columns <- available_columns(dat, genes, "_count")

  if (length(columns) == 0L) {
    return(rep(NA, nrow(dat)))
  }

  rowSums(as.matrix(dat[, ..columns]) > 0, na.rm = TRUE) > 0
}

project_theme <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, colour = "grey30"),
      plot.caption = element_text(size = 8, colour = "grey40", hjust = 0),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "grey25"),
      axis.line = element_line(colour = "grey20", linewidth = 0.45),
      axis.ticks = element_line(colour = "grey30", linewidth = 0.35),
      panel.grid = element_blank(),
      legend.title = element_blank(),
      legend.position = "bottom",
      plot.margin = margin(8, 10, 8, 8)
    )
}

save_png <- function(plot, filename, width, height) {
  output_path <- file.path(figure_dir, filename)
  ggsave(
    filename = output_path,
    plot = plot,
    width = width,
    height = height,
    dpi = 600,
    bg = "white"
  )
  message("Saved: ", output_path)
  output_path
}

# ------------------------------------------------------------------------------
# Read annotations and selected expression rows
# ------------------------------------------------------------------------------

message("Reading annotations: ", annotation_file)
annotations <- fread(annotation_file, check.names = FALSE)

required_annotation_columns <- c(
  "cells", "samples", "cell.types", "treatment.group", "Cohort"
)
missing_annotation_columns <- setdiff(
  required_annotation_columns,
  names(annotations)
)

if (length(missing_annotation_columns) > 0L) {
  stop(
    "Annotation table is missing column(s): ",
    paste(missing_annotation_columns, collapse = ", ")
  )
}

if (anyDuplicated(annotations$cells)) {
  stop("Annotation cell identifiers are not unique.")
}

count_header <- read_matrix_header(count_file)
tpm_header <- read_matrix_header(tpm_file)
count_cells <- count_header[-1L]
tpm_cells <- tpm_header[-1L]

if (!identical(count_cells, tpm_cells)) {
  stop("Count and TPM matrix cell identifiers are not identical and ordered.")
}

if (!setequal(annotations$cells, tpm_cells)) {
  stop(
    "Annotation and expression-matrix cell identifiers do not contain the ",
    "same cells."
  )
}

count_selected <- read_selected_gene_rows(count_file, target_genes)
tpm_selected <- read_selected_gene_rows(tpm_file, target_genes)

if (!identical(names(count_selected)[-1L], count_cells)) {
  stop("Selected count rows do not preserve the count-matrix header.")
}

if (!identical(names(tpm_selected)[-1L], tpm_cells)) {
  stop("Selected TPM rows do not preserve the TPM-matrix header.")
}

count_cells_table <- selected_rows_to_cells(
  count_selected,
  target_genes,
  "_count"
)
tpm_cells_table <- selected_rows_to_cells(
  tpm_selected,
  target_genes,
  "_TPM"
)

cell_data <- merge(
  tpm_cells_table,
  count_cells_table,
  by = "cells",
  all = TRUE,
  sort = FALSE
)
cell_data <- merge(
  cell_data,
  annotations,
  by = "cells",
  all.x = TRUE,
  sort = FALSE
)

cell_data[, matrix_order__ := match(cells, tpm_cells)]
setorder(cell_data, matrix_order__)
cell_data[, matrix_order__ := NULL]

if (nrow(cell_data) != length(tpm_cells) || anyNA(cell_data$cell.types)) {
  stop("Cell-level table construction failed or produced missing annotations.")
}

# ------------------------------------------------------------------------------
# Harmonize labels and calculate marker evidence
# ------------------------------------------------------------------------------

cell_type_labels <- c(
  "Mal" = "Malignant",
  "T.CD8" = "CD8 T cells",
  "T.CD4" = "CD4 T cells",
  "T.cell" = "T cells, unspecified",
  "B.cell" = "B cells",
  "Macrophage" = "Macrophages",
  "CAF" = "Fibroblasts",
  "Endo." = "Endothelial cells",
  "NK" = "NK cells",
  "?" = "Unclassified"
)

cell_data[, cell_type_original := cell.types]
cell_data[, cell_type := unname(cell_type_labels[cell.types])]
cell_data[is.na(cell_type), cell_type := cell_type_original]

if (all(is.na(cell_data$SPN_count)) || all(is.na(cell_data$SPN_TPM))) {
  stop("SPN was not available in both selected count and TPM data.")
}

if (all(is.na(cell_data$PTPRC_count)) || all(is.na(cell_data$PTPRC_TPM))) {
  stop("PTPRC was not available in both selected count and TPM data.")
}

cell_data[, SPN_log2_TPM_plus_1 := log2(SPN_TPM + 1)]
cell_data[, PTPRC_log2_TPM_plus_1 := log2(PTPRC_TPM + 1)]
cell_data[, SPN_detected := SPN_count > 0]
cell_data[, PTPRC_detected := PTPRC_count > 0]

cell_data[, melanoma_marker_score := mean_log2_tpm_score(
  cell_data, marker_sets$melanoma
)]
cell_data[, t_cell_marker_score := mean_log2_tpm_score(
  cell_data, marker_sets$t_cell
)]
cell_data[, cd8_marker_score := mean_log2_tpm_score(
  cell_data, marker_sets$cd8
)]
cell_data[, melanoma_marker_detected := any_count_detected(
  cell_data, marker_sets$melanoma
)]

cell_data[, SPN_PTPRC_quadrant := fcase(
  SPN_detected & PTPRC_detected, "SPN+ PTPRC+",
  SPN_detected & !PTPRC_detected, "SPN+ PTPRC-",
  !SPN_detected & PTPRC_detected, "SPN- PTPRC+",
  default = "SPN- PTPRC-"
)]

cell_data[, malignant_SPN_PTPRC_negative :=
  cell_type == "Malignant" &
    SPN_detected &
    !PTPRC_detected]

cell_data[, source_highlight := fifelse(
  malignant_SPN_PTPRC_negative,
  "Malignant SPN+ / PTPRC-",
  "All other cells"
)]

# ------------------------------------------------------------------------------
# Marker availability and QC
# ------------------------------------------------------------------------------

marker_availability <- data.table(
  gene = target_genes,
  marker_set = vapply(target_genes, function(gene) {
    names(marker_sets)[vapply(marker_sets, function(set) gene %in% set,
                             logical(1))][1L]
  }, character(1)),
  present_in_counts = target_genes %in% count_selected$gene,
  present_in_tpm = target_genes %in% tpm_selected$gene
)

fwrite(
  marker_availability,
  file.path(results_dir, "10_marker_availability.csv")
)

cell_type_counts <- cell_data[, .(
  n_cells = .N,
  n_samples = uniqueN(samples),
  percent_of_dataset = 100 * .N / nrow(cell_data)
), by = cell_type][order(-n_cells)]

fwrite(
  cell_type_counts,
  file.path(results_dir, "10_cell_type_counts.csv")
)

qc_summary <- data.table(
  metric = c(
    "annotation_cells",
    "count_matrix_cells",
    "tpm_matrix_cells",
    "cell_ids_identical_between_count_and_tpm",
    "cell_sets_identical_between_annotations_and_expression",
    "unique_samples",
    "annotated_cell_types",
    "requested_marker_genes",
    "markers_present_in_counts",
    "markers_present_in_tpm",
    "SPN_detected_cells",
    "PTPRC_detected_cells",
    "malignant_cells",
    "malignant_SPN_positive_PTPRC_negative_cells"
  ),
  value = as.character(c(
    nrow(annotations),
    length(count_cells),
    length(tpm_cells),
    identical(count_cells, tpm_cells),
    setequal(annotations$cells, tpm_cells),
    uniqueN(cell_data$samples),
    uniqueN(cell_data$cell_type),
    length(target_genes),
    sum(marker_availability$present_in_counts),
    sum(marker_availability$present_in_tpm),
    sum(cell_data$SPN_detected, na.rm = TRUE),
    sum(cell_data$PTPRC_detected, na.rm = TRUE),
    sum(cell_data$cell_type == "Malignant"),
    sum(cell_data$malignant_SPN_PTPRC_negative, na.rm = TRUE)
  ))
)

fwrite(qc_summary, file.path(results_dir, "10_dataset_qc.csv"))

# ------------------------------------------------------------------------------
# Cellular-source summaries
# ------------------------------------------------------------------------------

source_summary <- cell_data[, .(
  n_cells = .N,
  n_samples = uniqueN(samples),
  SPN_detected_n = sum(SPN_detected, na.rm = TRUE),
  SPN_detected_percent = 100 * mean(SPN_detected, na.rm = TRUE),
  SPN_TPM_mean = mean(SPN_TPM, na.rm = TRUE),
  SPN_TPM_median = median(SPN_TPM, na.rm = TRUE),
  SPN_TPM_maximum = max(SPN_TPM, na.rm = TRUE),
  total_SPN_TPM = sum(SPN_TPM, na.rm = TRUE),
  PTPRC_detected_percent = 100 * mean(PTPRC_detected, na.rm = TRUE),
  SPN_positive_PTPRC_negative_n = sum(
    SPN_detected & !PTPRC_detected,
    na.rm = TRUE
  ),
  SPN_positive_PTPRC_negative_percent = 100 * mean(
    SPN_detected & !PTPRC_detected,
    na.rm = TRUE
  ),
  melanoma_marker_detected_percent = 100 * mean(
    melanoma_marker_detected,
    na.rm = TRUE
  ),
  melanoma_marker_score_median = median(
    melanoma_marker_score,
    na.rm = TRUE
  )
), by = cell_type]

source_summary[, SPN_TPM_contribution_percent :=
  100 * total_SPN_TPM / sum(total_SPN_TPM)]
setorder(source_summary, -SPN_TPM_mean)

fwrite(
  source_summary,
  file.path(results_dir, "10_spn_source_by_cell_type.csv")
)

quadrant_summary <- cell_data[, .N, by = .(
  cell_type,
  SPN_PTPRC_quadrant
)]
quadrant_summary[, percent_within_cell_type := 100 * N / sum(N),
                 by = cell_type]
setorder(quadrant_summary, cell_type, SPN_PTPRC_quadrant)

fwrite(
  quadrant_summary,
  file.path(results_dir, "10_spn_ptprc_quadrants_by_cell_type.csv")
)

sample_cell_type_summary <- cell_data[, .(
  n_cells = .N,
  SPN_detected_percent = 100 * mean(SPN_detected, na.rm = TRUE),
  SPN_TPM_mean = mean(SPN_TPM, na.rm = TRUE),
  SPN_TPM_median = median(SPN_TPM, na.rm = TRUE),
  PTPRC_detected_percent = 100 * mean(PTPRC_detected, na.rm = TRUE),
  PTPRC_TPM_mean = mean(PTPRC_TPM, na.rm = TRUE),
  melanoma_marker_score_mean = mean(melanoma_marker_score, na.rm = TRUE),
  cd8_marker_score_mean = mean(cd8_marker_score, na.rm = TRUE)
), by = .(
  samples,
  treatment.group,
  Cohort,
  cell_type
)]

fwrite(
  sample_cell_type_summary,
  file.path(results_dir, "10_spn_source_by_sample_and_cell_type.csv")
)

malignant_sample_evidence <- cell_data[cell_type == "Malignant", .(
  n_malignant_cells = .N,
  SPN_detected_n = sum(SPN_detected, na.rm = TRUE),
  SPN_detected_percent = 100 * mean(SPN_detected, na.rm = TRUE),
  SPN_positive_PTPRC_negative_n = sum(
    SPN_detected & !PTPRC_detected,
    na.rm = TRUE
  ),
  SPN_positive_PTPRC_negative_percent = 100 * mean(
    SPN_detected & !PTPRC_detected,
    na.rm = TRUE
  ),
  melanoma_marker_detected_percent = 100 * mean(
    melanoma_marker_detected,
    na.rm = TRUE
  ),
  SPN_TPM_mean = mean(SPN_TPM, na.rm = TRUE),
  melanoma_marker_score_mean = mean(melanoma_marker_score, na.rm = TRUE)
), by = .(
  samples,
  treatment.group,
  Cohort
)]

fwrite(
  malignant_sample_evidence,
  file.path(results_dir, "10_malignant_spn_evidence_by_sample.csv")
)

cell_output_file <- file.path(
  analysis_dir,
  "10_GSE115978_selected_marker_cell_table.tsv.gz"
)
fwrite(cell_data, cell_output_file, sep = "\t")

# ------------------------------------------------------------------------------
# Figure 10A: SPN expression by annotated cell type
# ------------------------------------------------------------------------------

cell_type_order <- source_summary[order(SPN_TPM_median)]$cell_type
plot_data <- copy(cell_data)
plot_data[, cell_type := factor(cell_type, levels = cell_type_order)]
plot_data[, malignant_status := fifelse(
  cell_type == "Malignant",
  "Malignant cells",
  "Other cell types"
)]

fill_values <- c(
  "Other cell types" = "#EEEEEE",
  "Malignant cells" = "#C9A3FF"
)
border_values <- c(
  "Other cell types" = "#8A8A8A",
  "Malignant cells" = "#7551A8"
)

expression_plot <- ggplot(
  plot_data,
  aes(x = cell_type, y = SPN_log2_TPM_plus_1, fill = malignant_status)
) +
  geom_violin(
    aes(colour = malignant_status),
    scale = "width",
    trim = TRUE,
    linewidth = 0.35
  ) +
  geom_boxplot(
    width = 0.16,
    outlier.shape = NA,
    fill = "white",
    colour = "grey20",
    linewidth = 0.35
  ) +
  scale_fill_manual(values = fill_values) +
  scale_colour_manual(values = border_values) +
  labs(
    title = "SPN expression across GSE115978 melanoma cell types",
    subtitle = paste0(
      "Cell-level distribution; malignant cells are highlighted in lavender"
    ),
    x = NULL,
    y = "log2(SPN TPM + 1)",
    caption = paste0(
      "Cell-level distributions are descriptive; cells from the same patient ",
      "are not independent."
    )
  ) +
  coord_flip() +
  project_theme()

expression_path <- save_png(
  expression_plot,
  "10A_GSE115978_SPN_expression_by_cell_type.png",
  width = 8.5,
  height = 6.8
)

# ------------------------------------------------------------------------------
# Figure 10B: SPN detection by annotated cell type
# ------------------------------------------------------------------------------

detection_plot_data <- copy(source_summary)
detection_plot_data[, cell_type := factor(
  cell_type,
  levels = source_summary[order(SPN_detected_percent)]$cell_type
)]
detection_plot_data[, malignant_status := fifelse(
  cell_type == "Malignant",
  "Malignant cells",
  "Other cell types"
)]

detection_plot <- ggplot(
  detection_plot_data,
  aes(x = cell_type, y = SPN_detected_percent, fill = malignant_status)
) +
  geom_col(
    aes(colour = malignant_status),
    width = 0.68,
    linewidth = 0.45
  ) +
  geom_text(
    aes(label = paste0(round(SPN_detected_percent, 1), "%")),
    hjust = -0.10,
    size = 3.1,
    colour = "grey20"
  ) +
  scale_fill_manual(values = fill_values) +
  scale_colour_manual(values = border_values) +
  scale_y_continuous(
    limits = c(0, max(detection_plot_data$SPN_detected_percent) * 1.16),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    title = "SPN detection frequency by cell type",
    subtitle = "Detection is defined as SPN raw count > 0",
    x = NULL,
    y = "SPN-detected cells (%)",
    caption = "Lavender identifies the annotated malignant compartment."
  ) +
  coord_flip() +
  project_theme()

detection_path <- save_png(
  detection_plot,
  "10B_GSE115978_SPN_detection_by_cell_type.png",
  width = 8.5,
  height = 6.8
)

# ------------------------------------------------------------------------------
# Figure 10C: SPN versus PTPRC source evidence
# ------------------------------------------------------------------------------

correlation_test <- suppressWarnings(cor.test(
  cell_data$SPN_log2_TPM_plus_1,
  cell_data$PTPRC_log2_TPM_plus_1,
  method = "spearman",
  exact = FALSE
))

scatter_plot_data <- copy(cell_data)
scatter_plot_data[, source_highlight := factor(
  source_highlight,
  levels = c("All other cells", "Malignant SPN+ / PTPRC-")
)]

source_scatter <- ggplot(
  scatter_plot_data,
  aes(x = SPN_log2_TPM_plus_1, y = PTPRC_log2_TPM_plus_1)
) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    colour = "grey20",
    fill = "grey85",
    linewidth = 0.7,
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = 0,
    colour = "grey70",
    linewidth = 0.4,
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = 0,
    colour = "grey70",
    linewidth = 0.4,
    linetype = "dashed"
  ) +
  geom_point(
    aes(fill = source_highlight, colour = source_highlight),
    shape = 21,
    size = 2.0,
    stroke = 0.45,
    alpha = 0.78
  ) +
  scale_fill_manual(values = c(
    "All other cells" = "#EEEEEE",
    "Malignant SPN+ / PTPRC-" = "#C9A3FF"
  )) +
  scale_colour_manual(values = c(
    "All other cells" = "#8A8A8A",
    "Malignant SPN+ / PTPRC-" = "#7551A8"
  )) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    hjust = 1.06,
    vjust = 1.35,
    fontface = "bold",
    size = 3.6,
    label = paste0(
      "GSE115978, n = ", nrow(cell_data),
      "\nSpearman rho = ", round(unname(correlation_test$estimate), 2),
      "\nMalignant SPN+/PTPRC- = ",
      sum(cell_data$malignant_SPN_PTPRC_negative)
    )
  ) +
  labs(
    title = "SPN and PTPRC expression across melanoma single cells",
    subtitle = paste0(
      "Lavender highlights malignant cells with SPN detected but no PTPRC ",
      "raw counts"
    ),
    x = "log2(SPN TPM + 1)",
    y = "log2(PTPRC TPM + 1)",
    caption = paste0(
      "The dashed line is a descriptive linear fit. PTPRC non-detection may ",
      "reflect biological absence or single-cell dropout."
    )
  ) +
  project_theme()

scatter_path <- save_png(
  source_scatter,
  "10C_GSE115978_SPN_vs_PTPRC_source_scatter.png",
  width = 7.4,
  height = 7.0
)

# ------------------------------------------------------------------------------
# Figure 10D: estimated contribution of each sampled cell compartment
# ------------------------------------------------------------------------------

contribution_plot_data <- copy(source_summary)
contribution_plot_data[, cell_type := factor(
  cell_type,
  levels = source_summary[order(SPN_TPM_contribution_percent)]$cell_type
)]
contribution_plot_data[, malignant_status := fifelse(
  cell_type == "Malignant",
  "Malignant cells",
  "Other cell types"
)]

contribution_plot <- ggplot(
  contribution_plot_data,
  aes(
    x = cell_type,
    y = SPN_TPM_contribution_percent,
    fill = malignant_status
  )
) +
  geom_col(
    aes(colour = malignant_status),
    width = 0.68,
    linewidth = 0.45
  ) +
  geom_text(
    aes(label = paste0(round(SPN_TPM_contribution_percent, 1), "%")),
    hjust = -0.10,
    size = 3.1,
    colour = "grey20"
  ) +
  scale_fill_manual(values = fill_values) +
  scale_colour_manual(values = border_values) +
  scale_y_continuous(
    limits = c(
      0,
      max(contribution_plot_data$SPN_TPM_contribution_percent) * 1.16
    ),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    title = "Sampled-cell contribution to total SPN TPM",
    subtitle = paste0(
      "Contribution depends on both per-cell expression and the number of ",
      "cells sampled"
    ),
    x = NULL,
    y = "Contribution to summed SPN TPM (%)",
    caption = paste0(
      "This is a dataset-composition summary, not an estimate of true tumor ",
      "tissue abundance."
    )
  ) +
  coord_flip() +
  project_theme()

contribution_path <- save_png(
  contribution_plot,
  "10D_GSE115978_sampled_cell_SPN_contribution.png",
  width = 8.5,
  height = 6.8
)

# ------------------------------------------------------------------------------
# Manifest, interpretation guardrails, and completion
# ------------------------------------------------------------------------------

figure_manifest <- data.table(
  figure = c("10A", "10B", "10C", "10D"),
  description = c(
    "SPN expression distribution by annotated cell type",
    "SPN raw-count detection frequency by annotated cell type",
    "SPN versus PTPRC scatter highlighting malignant SPN+/PTPRC- cells",
    "Sampled-cell contribution to summed SPN TPM"
  ),
  file = c(
    expression_path,
    detection_path,
    scatter_path,
    contribution_path
  )
)

fwrite(
  figure_manifest,
  file.path(results_dir, "10_figure_manifest.csv")
)

interpretation_notes <- c(
  "Primary question: Which annotated cell populations express SPN in GSE115978?",
  "Evidence supporting tumor-cell expression includes SPN detection in annotated malignant cells, absence of PTPRC counts in the same cells, and retained melanoma-lineage marker expression.",
  "This analysis does not prove that bulk TCGA SPN originates primarily from malignant cells.",
  "Single-cell dropout can create apparent PTPRC-negative cells, and ambient RNA can create low-level SPN detection.",
  "Cell numbers are sampling-dependent, so summed-TPM contribution is not equivalent to tissue composition.",
  "Patient-by-cell-type summaries should be prioritized over cell-level hypothesis tests in follow-up analyses.",
  "The next step is CD8-state refinement and signature construction using an independent CD8-rich melanoma dataset such as GSE123139."
)

writeLines(
  interpretation_notes,
  file.path(results_dir, "10_interpretation_notes.txt")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(results_dir, "10_sessionInfo.txt")
)

message("Step 10 complete.")
message("Cell-level reusable table: ", cell_output_file)
message("Results: ", results_dir)

cat("\nDataset QC:\n")
print(qc_summary)
cat("\nSPN source summary:\n")
print(source_summary)
cat("\nMalignant sample evidence:\n")
print(malignant_sample_evidence)
