
#!/usr/bin/env Rscript

# ==============================================================================
# Step 12: Calculate cell-level QC metrics for GSE123139
# ==============================================================================
#
# Purpose:
#   Stream through the 204 MARS-seq plate files, calculate cell-level quality
#   metrics, attach the Step 11 metadata, and evaluate candidate QC thresholds.
#
# Important:
#   This script does not remove cells and does not build the final expression
#   object. The candidate_qc_pass column is diagnostic and will be reviewed
#   before the filtered sparse matrix is created in the next step.
#
# Inputs:
#   data/scRNA/GSE123139/raw_counts/GSM*_*.txt.gz
#   data/analysis_ready/scRNA/GSE123139/
#     11_GSE123139_cell_metadata.csv.gz
#
# Outputs:
#   results/12/12_dataset_qc.csv
#   results/12/12_qc_summary_by_source.csv
#   results/12/12_qc_summary_by_plate.csv
#   results/12/12_gene_inventory_by_plate.csv
#   results/12/12_qc_threshold_sensitivity.csv
#   results/12/12_marker_detection_summary.csv
#   results/12/12_figure_manifest.csv
#   results/12/12_sessionInfo.txt
#   results/12/figures/12A_cell_qc_distributions.png
#   results/12/figures/12B_library_size_vs_detected_genes.png
#   results/12/figures/12C_candidate_qc_retention_by_plate.png
#   data/analysis_ready/scRNA/GSE123139/
#     12_GSE123139_cell_qc.tsv.gz
#
# Run from the project root:
#   Rscript src/12_gse123139_cell_qc.R
# ==============================================================================

options(stringsAsFactors = FALSE)

required_packages <- c("data.table", "ggplot2")
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
  library(ggplot2)
})

# ------------------------------------------------------------------------------
# Paths and configuration
# ------------------------------------------------------------------------------

project_root <- normalizePath(getwd(), mustWork = TRUE)

raw_counts_dir <- file.path(
  project_root,
  "data", "scRNA", "GSE123139", "raw_counts"
)

analysis_ready_dir <- file.path(
  project_root,
  "data", "analysis_ready", "scRNA", "GSE123139"
)

results_dir <- file.path(project_root, "results", "12")
figures_dir <- file.path(results_dir, "figures")

dir.create(analysis_ready_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

metadata_file <- file.path(
  analysis_ready_dir,
  "11_GSE123139_cell_metadata.csv.gz"
)

if (!file.exists(metadata_file)) {
  stop("Missing Step 11 cell metadata: ", metadata_file)
}

plate_files <- sort(list.files(
  raw_counts_dir,
  pattern = "^GSM[0-9]+_.+\\.txt(\\.gz)?$",
  full.names = TRUE,
  recursive = TRUE
))

if (length(plate_files) == 0L) {
  stop("No GSE123139 plate count files found in: ", raw_counts_dir)
}

cell_metadata <- fread(metadata_file)

# GEO's "patient id" field contains specimen-level suffixes in this study
# (for example, p2-4-LN-2IT). Preserve it as specimen_id and derive the shared
# subject prefix for patient-aware downstream analyses.
cell_metadata[, specimen_id := patient_id]
cell_metadata[, subject_id := sub(
  "^([pP][0-9]+).*",
  "\\1",
  specimen_id,
  perl = TRUE
)]

candidate_min_umi <- 800
candidate_min_genes <- 200
candidate_max_mito_percent <- 20
high_library_nmads <- 4

marker_genes <- c(
  "PTPRC", "CD3D", "CD3E", "TRAC", "CD4", "CD8A", "CD8B", "SPN",
  "PDCD1", "HAVCR2", "LAG3", "TIGIT", "TOX", "CXCL13",
  "TCF7", "IL7R", "CCR7", "NKG7", "GNLY", "PRF1", "GZMB",
  "FGFBP2", "MKI67", "TOP2A", "STMN1"
)

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

safe_percent <- function(numerator, denominator) {
  fifelse(denominator > 0, 100 * numerator / denominator, NA_real_)
}

marker_counts <- function(count_matrix, genes, marker) {
  rows <- which(genes == marker)

  if (length(rows) == 0L) {
    return(rep(NA_real_, ncol(count_matrix)))
  }

  if (length(rows) == 1L) {
    return(as.numeric(count_matrix[rows, ]))
  }

  colSums(count_matrix[rows, , drop = FALSE])
}

read_cell_header <- function(path) {
  connection <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }

  on.exit(close(connection))
  header_line <- readLines(connection, n = 1L, warn = FALSE)

  if (length(header_line) != 1L || !nzchar(header_line)) {
    stop("Could not read the cell header from: ", path)
  }

  strsplit(header_line, "\t", fixed = TRUE)[[1L]]
}

median_or_na <- function(x) {
  if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
}

# ------------------------------------------------------------------------------
# Stream through plates and calculate cell-level metrics
# ------------------------------------------------------------------------------

qc_list <- vector("list", length(plate_files))
gene_inventory_list <- vector("list", length(plate_files))
reference_genes <- NULL
all_genes <- character()

for (index in seq_along(plate_files)) {
  path <- plate_files[index]
  message(
    "Processing plate ", index, "/", length(plate_files),
    ": ", basename(path)
  )

  cell_ids <- read_cell_header(path)

  count_table <- fread(
    path,
    skip = 1L,
    header = FALSE,
    col.names = c("gene_symbol", cell_ids),
    check.names = FALSE
  )

  if (ncol(count_table) != length(cell_ids) + 1L) {
    stop(
      "Count-file width does not match its cell header: ",
      basename(path)
    )
  }

  genes <- as.character(count_table[[1L]])

  if (is.null(reference_genes)) {
    reference_genes <- genes
  }

  gene_inventory_list[[index]] <- data.table(
    count_file = basename(path),
    n_gene_rows = length(genes),
    n_unique_genes = uniqueN(genes),
    duplicated_gene_rows = sum(duplicated(genes)),
    identical_to_reference_order = identical(genes, reference_genes),
    identical_to_reference_gene_set = setequal(genes, reference_genes)
  )

  all_genes <- union(all_genes, genes)

  count_matrix <- as.matrix(count_table[, -1L, with = FALSE])
  storage.mode(count_matrix) <- "numeric"

  if (!identical(cell_ids, colnames(count_matrix))) {
    stop("Cell IDs were altered while reading: ", basename(path))
  }

  total_counts <- colSums(count_matrix)
  detected_genes <- colSums(count_matrix > 0)

  mitochondrial_rows <- grepl("^MT-", genes, ignore.case = FALSE)
  ribosomal_rows <- grepl("^RP[SL][0-9]", genes, ignore.case = FALSE)

  mitochondrial_counts <- if (any(mitochondrial_rows)) {
    colSums(count_matrix[mitochondrial_rows, , drop = FALSE])
  } else {
    rep(0, ncol(count_matrix))
  }

  ribosomal_counts <- if (any(ribosomal_rows)) {
    colSums(count_matrix[ribosomal_rows, , drop = FALSE])
  } else {
    rep(0, ncol(count_matrix))
  }

  plate_qc <- data.table(
    cell_id = cell_ids,
    n_counts = as.numeric(total_counts),
    n_genes = as.integer(detected_genes),
    mitochondrial_counts = as.numeric(mitochondrial_counts),
    mitochondrial_percent = safe_percent(
      mitochondrial_counts,
      total_counts
    ),
    ribosomal_counts = as.numeric(ribosomal_counts),
    ribosomal_percent = safe_percent(ribosomal_counts, total_counts)
  )

  for (marker in marker_genes) {
    plate_qc[, paste0(marker, "_count") := marker_counts(
      count_matrix,
      genes,
      marker
    )]
  }

  qc_list[[index]] <- plate_qc

  rm(count_table, count_matrix, plate_qc)
  if (index %% 10L == 0L) {
    invisible(gc())
  }
}

qc_metrics <- rbindlist(qc_list, use.names = TRUE, fill = TRUE)
gene_inventory <- rbindlist(gene_inventory_list, use.names = TRUE, fill = TRUE)

fwrite(
  gene_inventory,
  file.path(results_dir, "12_gene_inventory_by_plate.csv")
)

if (anyDuplicated(qc_metrics$cell_id)) {
  stop("Duplicated cell IDs were found in the calculated QC metrics.")
}

# Preserve the Step 11 metadata ordering.
qc <- merge(
  cell_metadata,
  qc_metrics,
  by = "cell_id",
  all.x = TRUE,
  all.y = FALSE,
  sort = FALSE
)

if (nrow(qc) != nrow(cell_metadata)) {
  stop("Cell count changed during the metadata/QC join.")
}

missing_qc <- sum(is.na(qc$n_counts))
if (missing_qc > 0L) {
  stop(missing_qc, " metadata cells did not receive QC metrics.")
}

# ------------------------------------------------------------------------------
# Diagnostic flags; no cells are removed in this script
# ------------------------------------------------------------------------------

qc[, log10_counts_plus_1 := log10(n_counts + 1)]
qc[, log10_genes_plus_1 := log10(n_genes + 1)]

qc[, high_library_cutoff := {
  center <- median(log10_counts_plus_1, na.rm = TRUE)
  spread <- mad(log10_counts_plus_1, center = center, na.rm = TRUE)
  10^(center + high_library_nmads * spread) - 1
}, by = gsm_accession]

qc[, high_library_outlier := n_counts > high_library_cutoff]
qc[, zero_library := n_counts == 0]

qc[, candidate_qc_pass :=
  n_counts >= candidate_min_umi &
  n_genes >= candidate_min_genes &
  !is.na(mitochondrial_percent) &
  mitochondrial_percent <= candidate_max_mito_percent &
  !high_library_outlier
]

qc[, source_group := fifelse(
  tolower(sample_source) == "tumor",
  "Tumor",
  fifelse(tolower(sample_source) == "pbmc", "PBMC", "Other/unknown")
)]

qc_output_file <- file.path(
  analysis_ready_dir,
  "12_GSE123139_cell_qc.tsv.gz"
)

fwrite(qc, qc_output_file, sep = "\t")

# ------------------------------------------------------------------------------
# Summary tables
# ------------------------------------------------------------------------------

dataset_qc <- data.table(
  metric = c(
    "plate_files",
    "cell_columns",
    "cells_with_qc_metrics",
    "zero_count_wells",
    "candidate_qc_pass",
    "candidate_qc_fail",
    "candidate_min_umi",
    "candidate_min_genes",
    "candidate_max_mito_percent",
    "high_library_nmads",
    "unique_specimens",
    "unique_subjects",
    "genes_in_union_across_plates",
    "plates_with_reference_gene_order",
    "plates_with_reference_gene_set",
    "mitochondrial_genes_detected",
    "ribosomal_genes_detected"
  ),
  value = as.character(c(
    length(plate_files),
    nrow(qc),
    sum(!is.na(qc$n_counts)),
    sum(qc$zero_library),
    sum(qc$candidate_qc_pass),
    sum(!qc$candidate_qc_pass),
    candidate_min_umi,
    candidate_min_genes,
    candidate_max_mito_percent,
    high_library_nmads,
    uniqueN(qc$specimen_id),
    uniqueN(qc$subject_id),
    length(all_genes),
    sum(gene_inventory$identical_to_reference_order),
    sum(gene_inventory$identical_to_reference_gene_set),
    sum(grepl("^MT-", all_genes)),
    sum(grepl("^RP[SL][0-9]", all_genes))
  ))
)

fwrite(dataset_qc, file.path(results_dir, "12_dataset_qc.csv"))

qc_summary_by_source <- qc[, .(
  n_cell_columns = .N,
  n_zero_count_wells = sum(zero_library),
  median_counts = median_or_na(n_counts),
  median_genes = median_or_na(n_genes),
  median_mitochondrial_percent = median_or_na(mitochondrial_percent),
  n_candidate_qc_pass = sum(candidate_qc_pass),
  candidate_qc_pass_percent = 100 * mean(candidate_qc_pass)
), by = source_group][order(source_group)]

fwrite(
  qc_summary_by_source,
  file.path(results_dir, "12_qc_summary_by_source.csv")
)

qc_summary_by_plate <- qc[, .(
  sample_source = unique(sample_source)[1L],
  subject_id = unique(subject_id)[1L],
  patient_id = unique(patient_id)[1L],
  facs_gate = unique(facs_gate)[1L],
  n_cell_columns = .N,
  median_counts = median_or_na(n_counts),
  median_genes = median_or_na(n_genes),
  median_mitochondrial_percent = median_or_na(mitochondrial_percent),
  n_candidate_qc_pass = sum(candidate_qc_pass),
  candidate_qc_pass_percent = 100 * mean(candidate_qc_pass)
), by = .(gsm_accession, amplification_batch, plate_id)]

setorder(qc_summary_by_plate, sample_source, patient_id, gsm_accession)

fwrite(
  qc_summary_by_plate,
  file.path(results_dir, "12_qc_summary_by_plate.csv")
)

threshold_grid <- CJ(
  minimum_umi = c(200, 500, 800, 1000, 1500),
  minimum_genes = c(100, 200),
  maximum_mitochondrial_percent = c(10, 20, 30)
)

threshold_sensitivity <- rbindlist(lapply(
  seq_len(nrow(threshold_grid)),
  function(index) {
    threshold <- threshold_grid[index]

    qc[, .(
      n_cell_columns = .N,
      n_pass = sum(
        n_counts >= threshold$minimum_umi &
        n_genes >= threshold$minimum_genes &
        !is.na(mitochondrial_percent) &
        mitochondrial_percent <=
          threshold$maximum_mitochondrial_percent &
        !high_library_outlier
      ),
      percent_pass = 100 * mean(
        n_counts >= threshold$minimum_umi &
        n_genes >= threshold$minimum_genes &
        !is.na(mitochondrial_percent) &
        mitochondrial_percent <=
          threshold$maximum_mitochondrial_percent &
        !high_library_outlier
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
    "maximum_mitochondrial_percent", "n_cell_columns",
    "n_pass", "percent_pass"
  )
)

fwrite(
  threshold_sensitivity,
  file.path(results_dir, "12_qc_threshold_sensitivity.csv")
)

marker_detection_summary <- rbindlist(lapply(marker_genes, function(marker) {
  column <- paste0(marker, "_count")

  qc[, .(
    gene = marker,
    n_cells_detected = sum(get(column) > 0, na.rm = TRUE),
    percent_cells_detected = 100 * mean(get(column) > 0, na.rm = TRUE),
    n_candidate_pass_detected = sum(
      candidate_qc_pass & get(column) > 0,
      na.rm = TRUE
    ),
    percent_candidate_pass_detected = 100 * mean(
      get(column)[candidate_qc_pass] > 0,
      na.rm = TRUE
    )
  ), by = source_group]
}))

fwrite(
  marker_detection_summary,
  file.path(results_dir, "12_marker_detection_summary.csv")
)

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

distribution_data <- rbindlist(list(
  qc[, .(
    source_group,
    metric = "log10(UMI counts + 1)",
    value = log10_counts_plus_1
  )],
  qc[, .(
    source_group,
    metric = "Detected genes",
    value = as.numeric(n_genes)
  )],
  qc[, .(
    source_group,
    metric = "Mitochondrial UMIs (%)",
    value = mitochondrial_percent
  )]
))

plot_12a <- ggplot(
  distribution_data[is.finite(value)],
  aes(x = value, fill = source_group, color = source_group)
) +
  geom_density(alpha = 0.30, linewidth = 0.7) +
  facet_wrap(~metric, scales = "free", ncol = 1) +
  scale_fill_manual(values = c(
    "Tumor" = lavender,
    "PBMC" = gray_fill,
    "Other/unknown" = "#D8D8D8"
  )) +
  scale_color_manual(values = c(
    "Tumor" = lavender_dark,
    "PBMC" = gray_outline,
    "Other/unknown" = "#777777"
  )) +
  labs(
    title = "GSE123139 cell-level quality distributions",
    subtitle = "All sorted wells are shown before final filtering",
    x = NULL,
    y = "Density",
    fill = "Sample source",
    color = "Sample source"
  ) +
  theme_project

figure_12a <- file.path(
  figures_dir,
  "12A_cell_qc_distributions.png"
)

ggsave(
  figure_12a,
  plot_12a,
  width = 9,
  height = 11,
  dpi = 300,
  bg = "white"
)

set.seed(123139)
plot_cells <- if (nrow(qc) > 40000L) {
  qc[sample(.N, 40000L)]
} else {
  copy(qc)
}

setorder(plot_cells, candidate_qc_pass)

plot_12b <- ggplot(
  plot_cells,
  aes(x = log10_counts_plus_1, y = n_genes)
) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "#222222",
    fill = gray_fill,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  geom_point(
    aes(fill = candidate_qc_pass),
    shape = 21,
    size = 1.7,
    stroke = 0.35,
    color = gray_outline,
    alpha = 0.65
  ) +
  scale_fill_manual(
    values = c("FALSE" = gray_fill, "TRUE" = lavender),
    labels = c("FALSE" = "Other wells", "TRUE" = "Candidate QC pass")
  ) +
  facet_wrap(~source_group) +
  labs(
    title = "Library complexity across GSE123139 cells",
    subtitle = paste0(
      "Candidate pass: ≥", candidate_min_umi, " UMIs, ≥",
      candidate_min_genes, " genes, ≤", candidate_max_mito_percent,
      "% mitochondrial UMIs; 40,000 cells plotted at most"
    ),
    x = "log10(UMI counts + 1)",
    y = "Detected genes",
    fill = NULL
  ) +
  theme_project

figure_12b <- file.path(
  figures_dir,
  "12B_library_size_vs_detected_genes.png"
)

ggsave(
  figure_12b,
  plot_12b,
  width = 12,
  height = 6,
  dpi = 300,
  bg = "white"
)

plate_plot_data <- copy(qc_summary_by_plate)
plate_plot_data[, plate_label := paste0(gsm_accession, " | ", patient_id)]
setorder(plate_plot_data, sample_source, candidate_qc_pass_percent)
plate_plot_data[, plate_label := factor(plate_label, levels = plate_label)]

plot_12c <- ggplot(
  plate_plot_data,
  aes(
    x = plate_label,
    y = candidate_qc_pass_percent,
    fill = sample_source
  )
) +
  geom_col(color = gray_outline, linewidth = 0.25) +
  coord_flip() +
  scale_fill_manual(values = c(
    "Tumor" = lavender,
    "PBMC" = gray_fill
  )) +
  labs(
    title = "Candidate QC retention across GSE123139 plates",
    subtitle = "Large plate-to-plate differences will be reviewed before filtering",
    x = NULL,
    y = "Candidate QC pass (%)",
    fill = "Sample source"
  ) +
  theme_project +
  theme(
    axis.text.y = element_text(size = 4.5),
    legend.position = "bottom"
  )

figure_12c <- file.path(
  figures_dir,
  "12C_candidate_qc_retention_by_plate.png"
)

ggsave(
  figure_12c,
  plot_12c,
  width = 10,
  height = 28,
  dpi = 300,
  bg = "white"
)

figure_manifest <- data.table(
  figure = c("12A", "12B", "12C"),
  file = c(figure_12a, figure_12b, figure_12c),
  description = c(
    "Cell-level QC distributions by sample source",
    "Library size versus detected genes with candidate-pass cells highlighted",
    "Candidate QC retention percentage across plates"
  )
)

fwrite(
  figure_manifest,
  file.path(results_dir, "12_figure_manifest.csv")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(results_dir, "12_sessionInfo.txt")
)

message("Step 12 complete.")
message("Cell QC table: ", qc_output_file)
message("Candidate QC pass: ", sum(qc$candidate_qc_pass), "/", nrow(qc))
message("Results: ", results_dir)

print(dataset_qc)
print(qc_summary_by_source)
