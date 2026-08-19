#!/usr/bin/env Rscript

# Step 2: Calculate complementary immune-infiltration measurements for
# TCGA-SKCM and TCGA-UVM.
#
# Important: CIBERSORTx absolute-mode values are intentionally excluded because
# diagnostic inspection showed that all cell scores were approximately zero
# (about 1e-17), making them technically invalid as abundance measurements.
#
# Input:
#   data/analysis_ready/skcm_uvm_tables/07_tcga_analysis_table.tsv.gz
#
# Outputs:
#   data/analysis_ready/skcm_uvm_immune_scores/
#     08_tcga_analysis_with_immune_scores.tsv.gz
#     08_tcga_immune_infiltration_scores.tsv.gz
#     08_<CANCER>_analysis_with_immune_scores.tsv.gz
#
#   results/08/
#     08_immune_infiltration_score_qc.csv
#     08_immune_infiltration_score_missingness.csv
#     08_immune_infiltration_score_summary.csv
#     08_sessionInfo.txt
#     figures/
#
# Run from the project root with:
#   Rscript src/08_build_immune_infiltration_scores.R

options(stringsAsFactors = FALSE)

input_table <- file.path(
  "data",
  "analysis_ready",
  "skcm_uvm_tables",
  "07_tcga_analysis_table.tsv.gz"
)
tpm_directory <- "data/tcga_star_counts/tpm_matrices"
data_output_directory <- file.path(
  "data",
  "analysis_ready",
  "skcm_uvm_immune_scores"
)
results_directory <- file.path("results", "08")
figure_directory <- file.path(results_directory, "figures")
cancers <- c("SKCM", "UVM")

invisible(lapply(
  c(data_output_directory, results_directory, figure_directory),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

required_markers <- c("SPN", "PTPRC", "CD3D", "CD3E", "TRAC", "CD8A", "CD8B")
t_cell_columns <- c(
  "T cells CD8",
  "T cells CD4 naive",
  "T cells CD4 memory resting",
  "T cells CD4 memory activated",
  "T cells follicular helper",
  "T cells regulatory (Tregs)",
  "T cells gamma delta"
)

lymphocyte_columns <- c(
  "B cells naive",
  "B cells memory",
  "Plasma cells",
  t_cell_columns,
  "NK cells resting",
  "NK cells activated"
)

read_tsv <- function(path) {
  if (!file.exists(path)) stop("Missing input file: ", path)
  read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA", "N/A", "Not Reported", "not reported", "[Not Available]")
  )
}

write_tsv_gz <- function(x, path) {
  connection <- gzfile(path, open = "wt")
  on.exit(close(connection), add = TRUE)
  write.table(x, connection, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

as_numeric_strict <- function(x, label) {
  output <- suppressWarnings(as.numeric(as.character(x)))
  newly_missing <- is.na(output) & !is.na(x)
  if (any(newly_missing)) {
    stop(label, " contains ", sum(newly_missing), " non-numeric value(s)")
  }
  output
}

cat("Reading Step 1 analysis table:", input_table, "\n")
analysis <- read_tsv(input_table)

required_analysis_columns <- unique(c(
  "sample_id", "patient", "cancer", "CD43_group", "SPN_TPM", "tumor_purity_primary",
  "CIBERSORTx_P_value", t_cell_columns, lymphocyte_columns
))
missing_analysis_columns <- setdiff(required_analysis_columns, names(analysis))
if (length(missing_analysis_columns) > 0) {
  stop(
    "Step 1 table is missing required column(s): ",
    paste(missing_analysis_columns, collapse = ", ")
  )
}

if (anyDuplicated(analysis$sample_id)) stop("Duplicate sample_id values in Step 1 table")
if (!all(analysis$cancer %in% cancers)) stop("Unexpected cancer value in Step 1 table")

for (column in unique(c("SPN_TPM", "tumor_purity_primary", "CIBERSORTx_P_value",
                        t_cell_columns, lymphocyte_columns))) {
  analysis[[column]] <- as_numeric_strict(analysis[[column]], column)
}

marker_tables <- list()
marker_qc <- list()

for (cancer in cancers) {
  tpm_file <- file.path(tpm_directory, paste0(cancer, "_TPM_primary_samples.tsv.gz"))
  cat("Reading TPM matrix:", tpm_file, "\n")
  tpm <- read_tsv(tpm_file)

  gene_column <- names(tpm)[1]
  genes <- as.character(tpm[[gene_column]])
  marker_counts <- table(factor(genes, levels = required_markers))

  if (any(marker_counts != 1L)) {
    stop(
      cancer, ": each required marker must occur exactly once. Counts: ",
      paste(names(marker_counts), marker_counts, sep = "=", collapse = ", ")
    )
  }

  sample_columns <- setdiff(names(tpm), gene_column)
  expected_samples <- analysis$sample_id[analysis$cancer == cancer]
  missing_in_tpm <- setdiff(expected_samples, sample_columns)
  extra_in_tpm <- setdiff(sample_columns, expected_samples)

  if (length(missing_in_tpm) > 0 || length(extra_in_tpm) > 0) {
    stop(
      cancer, ": TPM/sample mismatch. Missing in TPM: ", length(missing_in_tpm),
      "; extra in TPM: ", length(extra_in_tpm)
    )
  }

  marker_matrix <- as.matrix(tpm[match(required_markers, genes), expected_samples, drop = FALSE])
  storage.mode(marker_matrix) <- "numeric"
  rownames(marker_matrix) <- required_markers

  if (anyNA(marker_matrix)) stop(cancer, ": missing marker TPM values")
  if (any(marker_matrix < 0)) stop(cancer, ": negative marker TPM values")

  marker_log2 <- log2(marker_matrix + 1)
  marker_table <- data.frame(
    sample_id = expected_samples,
    cancer = cancer,
    SPN_TPM_from_matrix = marker_matrix["SPN", ],
    SPN_log2_TPM_plus_1_from_matrix = marker_log2["SPN", ],
    PTPRC_TPM = marker_matrix["PTPRC", ],
    PTPRC_log2_TPM_plus_1 = marker_log2["PTPRC", ],
    CD3D_log2_TPM_plus_1 = marker_log2["CD3D", ],
    CD3E_log2_TPM_plus_1 = marker_log2["CD3E", ],
    TRAC_log2_TPM_plus_1 = marker_log2["TRAC", ],
    CD8A_log2_TPM_plus_1 = marker_log2["CD8A", ],
    CD8B_log2_TPM_plus_1 = marker_log2["CD8B", ],
    T_cell_marker_score = colMeans(marker_log2[c("CD3D", "CD3E", "TRAC"), , drop = FALSE]),
    CD8_marker_score = colMeans(marker_log2[c("CD8A", "CD8B"), , drop = FALSE]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  marker_tables[[cancer]] <- marker_table
  marker_qc[[cancer]] <- data.frame(
    cancer = cancer,
    n_analysis_samples = length(expected_samples),
    n_tpm_samples = length(sample_columns),
    n_missing_analysis_samples_in_tpm = length(missing_in_tpm),
    n_extra_tpm_samples = length(extra_in_tpm),
    all_markers_present_once = all(marker_counts == 1L),
    stringsAsFactors = FALSE
  )
}

markers <- do.call(rbind, marker_tables)
rownames(markers) <- NULL

marker_order <- match(analysis$sample_id, markers$sample_id)
if (anyNA(marker_order)) stop("Internal error while joining TPM marker scores")
markers <- markers[marker_order, ]
if (!identical(analysis$cancer, markers$cancer)) stop("Cancer mismatch after marker join")

# Confirm that SPN extracted directly from the TPM matrix matches the value
# already carried in the Step 1 table.
spn_difference <- abs(analysis$SPN_TPM - markers$SPN_TPM_from_matrix)
if (any(spn_difference > 1e-8, na.rm = TRUE)) {
  stop("SPN TPM mismatch between Step 1 table and TPM matrix; maximum difference = ",
       max(spn_difference, na.rm = TRUE))
}

analysis$SPN_TPM_from_matrix <- markers$SPN_TPM_from_matrix
analysis$SPN_log2_TPM_plus_1_from_matrix <- markers$SPN_log2_TPM_plus_1_from_matrix
analysis$PTPRC_TPM <- markers$PTPRC_TPM
analysis$PTPRC_log2_TPM_plus_1 <- markers$PTPRC_log2_TPM_plus_1
analysis$CD3D_log2_TPM_plus_1 <- markers$CD3D_log2_TPM_plus_1
analysis$CD3E_log2_TPM_plus_1 <- markers$CD3E_log2_TPM_plus_1
analysis$TRAC_log2_TPM_plus_1 <- markers$TRAC_log2_TPM_plus_1
analysis$CD8A_log2_TPM_plus_1 <- markers$CD8A_log2_TPM_plus_1
analysis$CD8B_log2_TPM_plus_1 <- markers$CD8B_log2_TPM_plus_1
analysis$T_cell_marker_score <- markers$T_cell_marker_score
analysis$CD8_marker_score <- markers$CD8_marker_score

analysis$relative_total_T_cell_fraction <- rowSums(analysis[, t_cell_columns, drop = FALSE])
analysis$relative_total_lymphocyte_fraction <- rowSums(
  analysis[, lymphocyte_columns, drop = FALSE]
)
analysis$relative_CD8_within_T_cells <- ifelse(
  analysis$relative_total_T_cell_fraction > 0,
  analysis[["T cells CD8"]] / analysis$relative_total_T_cell_fraction,
  NA_real_
)
analysis$inverse_tumor_purity <- 1 - analysis$tumor_purity_primary

fraction_columns <- c(
  "T cells CD8", "relative_total_T_cell_fraction",
  "relative_total_lymphocyte_fraction", "relative_CD8_within_T_cells"
)
if (any(analysis[, fraction_columns[1:3], drop = FALSE] < -1e-8, na.rm = TRUE) ||
    any(analysis[, fraction_columns[1:3], drop = FALSE] > 1 + 1e-8, na.rm = TRUE)) {
  stop("One or more relative CIBERSORTx fractions are outside [0, 1]")
}
if (any(analysis$tumor_purity_primary < 0 | analysis$tumor_purity_primary > 1, na.rm = TRUE)) {
  stop("tumor_purity_primary contains values outside [0, 1]")
}

score_columns <- c(
  "SPN_TPM", "PTPRC_TPM", "PTPRC_log2_TPM_plus_1",
  "T_cell_marker_score", "CD8_marker_score",
  "T cells CD8", "relative_total_T_cell_fraction",
  "relative_total_lymphocyte_fraction", "relative_CD8_within_T_cells",
  "tumor_purity_primary", "inverse_tumor_purity"
)

score_table <- analysis[, c(
  "sample_id", "patient", "cancer", "CD43_group", "CIBERSORTx_P_value",
  score_columns
), drop = FALSE]

qc <- do.call(rbind, marker_qc)
qc$max_absolute_SPN_TPM_difference <- vapply(
  cancers,
  function(cancer) max(spn_difference[analysis$cancer == cancer], na.rm = TRUE),
  numeric(1)
)
qc$n_cibersortx_relative_pass_P_0_05 <- vapply(
  cancers,
  function(cancer) sum(
    analysis$CIBERSORTx_P_value[analysis$cancer == cancer] <= 0.05,
    na.rm = TRUE
  ),
  integer(1)
)
qc$n_with_primary_purity <- vapply(
  cancers,
  function(cancer) sum(!is.na(analysis$tumor_purity_primary[analysis$cancer == cancer])),
  integer(1)
)
qc$cibersortx_absolute_mode_status <- "excluded_all_scores_approximately_zero"
qc$estimate_status <- "not_calculated_package_unavailable"

missingness <- do.call(rbind, lapply(cancers, function(cancer) {
  y <- analysis[analysis$cancer == cancer, , drop = FALSE]
  do.call(rbind, lapply(score_columns, function(variable) {
    n_missing <- sum(is.na(y[[variable]]))
    data.frame(
      cancer = cancer,
      variable = variable,
      n_total = nrow(y),
      n_missing = n_missing,
      percent_missing = 100 * n_missing / nrow(y),
      stringsAsFactors = FALSE
    )
  }))
}))

summary_rows <- do.call(rbind, lapply(cancers, function(cancer) {
  y <- analysis[analysis$cancer == cancer, , drop = FALSE]
  do.call(rbind, lapply(score_columns, function(variable) {
    values <- y[[variable]]
    data.frame(
      cancer = cancer,
      variable = variable,
      n = sum(!is.na(values)),
      minimum = if (all(is.na(values))) NA_real_ else min(values, na.rm = TRUE),
      median = if (all(is.na(values))) NA_real_ else median(values, na.rm = TRUE),
      mean = if (all(is.na(values))) NA_real_ else mean(values, na.rm = TRUE),
      maximum = if (all(is.na(values))) NA_real_ else max(values, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
}))

combined_output <- file.path(
  data_output_directory,
  "08_tcga_analysis_with_immune_scores.tsv.gz"
)
score_output <- file.path(
  data_output_directory,
  "08_tcga_immune_infiltration_scores.tsv.gz"
)

write_tsv_gz(analysis, combined_output)
write_tsv_gz(score_table, score_output)

for (cancer in cancers) {
  write_tsv_gz(
    analysis[analysis$cancer == cancer, , drop = FALSE],
    file.path(
      data_output_directory,
      paste0("08_", cancer, "_analysis_with_immune_scores.tsv.gz")
    )
  )
}

write.csv(
  qc,
  file.path(results_directory, "08_immune_infiltration_score_qc.csv"),
  row.names = FALSE
)
write.csv(
  missingness,
  file.path(results_directory, "08_immune_infiltration_score_missingness.csv"),
  row.names = FALSE
)
write.csv(
  summary_rows,
  file.path(results_directory, "08_immune_infiltration_score_summary.csv"),
  row.names = FALSE
)
capture.output(
  sessionInfo(),
  file = file.path(results_directory, "08_sessionInfo.txt")
)

cat("\nStep 2 complete.\n")
cat("Combined table:", combined_output, "\n")
cat("Dimensions:", nrow(analysis), "rows x", ncol(analysis), "columns\n")
cat("\nQC summary:\n")
print(qc, row.names = FALSE)
cat("\nScore missingness:\n")
print(missingness, row.names = FALSE)
cat("\nCIBERSORTx absolute mode: excluded (all values approximately zero).\n")
cat("ESTIMATE: not calculated (package unavailable); not required for this step.\n")
