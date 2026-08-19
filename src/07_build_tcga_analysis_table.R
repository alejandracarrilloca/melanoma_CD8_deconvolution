#!/usr/bin/env Rscript
# Revised version: PANCAN_ABSOLUTE purity and cancer-specific clinical endpoints.

# Build the focused SKCM/UVM analysis table for the CD43-CD8 project.
#
# Input:
#   data/purity_analysis/all_samples_with_purity.tsv.gz
#
# Outputs:
#   data/analysis_ready/skcm_uvm_tables/07_tcga_analysis_table.tsv.gz
#   data/analysis_ready/skcm_uvm_tables/07_SKCM_analysis_table.tsv.gz
#   data/analysis_ready/skcm_uvm_tables/07_UVM_analysis_table.tsv.gz
#   results/07/07_tcga_analysis_table_qc.csv
#   results/07/07_tcga_analysis_table_missingness.csv
#   results/07/07_sessionInfo.txt
#
# Run from the project root with:
#   Rscript src/07_build_tcga_analysis_table.R

options(stringsAsFactors = FALSE)

input_file <- "data/purity_analysis/all_samples_with_purity.tsv.gz"
data_output_dir <- file.path(
  "data",
  "analysis_ready",
  "skcm_uvm_tables"
)
results_dir <- file.path("results", "07")

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

invisible(lapply(
  c(data_output_dir, results_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

message("Reading: ", input_file)

dat <- read.delim(
  input_file,
  check.names = FALSE,
  na.strings = c("", "NA", "N/A", "Not Reported", "not reported", "[Not Available]"),
  quote = "\"",
  comment.char = ""
)

required_columns <- c(
  "sample_id",
  "cancer",
  "SPN_TPM",
  "SPN_log2_TPM_plus_1",
  "CD43_group",
  "T cells CD8",
  "CIBERSORTx_P_value",
  "CIBERSORTx_Correlation",
  "CIBERSORTx_RMSE"
)

missing_required <- setdiff(required_columns, names(dat))

if (length(missing_required) > 0L) {
  stop(
    "Required columns are missing from the input: ",
    paste(missing_required, collapse = ", ")
  )
}

dat <- dat[dat$cancer %in% c("SKCM", "UVM"), , drop = FALSE]

if (nrow(dat) == 0L) {
  stop("No SKCM or UVM rows were found in the input table.")
}

immune_columns <- c(
  "B cells naive",
  "B cells memory",
  "Plasma cells",
  "T cells CD8",
  "T cells CD4 naive",
  "T cells CD4 memory resting",
  "T cells CD4 memory activated",
  "T cells follicular helper",
  "T cells regulatory (Tregs)",
  "T cells gamma delta",
  "NK cells resting",
  "NK cells activated",
  "Monocytes",
  "Macrophages M0",
  "Macrophages M1",
  "Macrophages M2",
  "Dendritic cells resting",
  "Dendritic cells activated",
  "Mast cells resting",
  "Mast cells activated",
  "Eosinophils",
  "Neutrophils"
)

candidate_columns <- c(
  # Identifiers and CD43 variables
  "sample_barcode",
  "sample_id",
  "patient",
  "sample",
  "barcode",
  "cancer",
  "project",
  "SPN_TPM",
  "SPN_log2_TPM_plus_1",
  "CD43_group",
  "Q1_TPM",
  "Q3_TPM",
  "Q1_log2_TPM_plus_1",
  "Q3_log2_TPM_plus_1",
  "cd43_tumor_vs_normal_direction",

  # LM22 fractions and deconvolution QC
  immune_columns,
  "LM22_fraction_sum",
  "CIBERSORTx_P_value",
  "CIBERSORTx_Correlation",
  "CIBERSORTx_RMSE",
  "CIBERSORTx_pass_P_0_05",
  "analysis_include_extreme_quartiles",

  # Sample context
  "shortLetterCode",
  "definition",
  "sample_type",
  "specimen_type",
  "tumor_descriptor",
  "requested_sample_type",
  "classification_of_tumor",
  "tissue_or_organ_of_origin",
  "site_of_resection_or_biopsy",
  "sites_of_involvement",
  "paper_ALL_PRIMARY_VS_METASTATIC",
  "paper_REGIONAL_VS_PRIMARY",
  "paper_CURATED_TCGA_SPECIMEN_SITE",
  "paper_CURATED_DISTANT_ANATOMIC_SITE",
  "paper_CURATED_TCGA_SPECIMEN_Distant",

  # Clinical variables
  "age_at_diagnosis",
  "age_at_index",
  "sex_at_birth",
  "race",
  "ethnicity",
  "ajcc_pathologic_stage",
  "ajcc_clinical_stage",
  "ajcc_pathologic_t",
  "ajcc_pathologic_n",
  "ajcc_pathologic_m",
  "ajcc_clinical_t",
  "ajcc_clinical_n",
  "ajcc_clinical_m",
  "metastasis_at_diagnosis",
  "melanoma_known_primary",
  "ulceration_indicator",
  "clark_level",
  "paper_CURATED_PATHOLOGIC_STAGE_AJCC7_AT_DIAGNOSIS_COMPLEX",
  "paper_CURATED_PATHOLOGIC_STAGE_AJCC7_AT_DIAGNOSIS_SIMPLE",
  "prior_malignancy",
  "prior_treatment",
  "year_of_diagnosis",

  # Survival candidates
  "vital_status",
  "days_to_death",
  "days_to_last_follow_up",
  "cause_of_death",
  "paper_CURATED_VITAL_STATUS",
  "paper_CURATED_DAYS_TO_DEATH_OR_LAST_FU",
  "paper_CURATED_TCGA_DAYS_TO_DEATH_OR_LAST_FU",
  "paper_CURATED_MELANOMA_SPECIFIC_VITAL_STATUS..0....ALIVE.OR.CENSORED...1....DEAD.OF.MELANOMA..",
  "paper_Death..Metastasis",

  # Purity estimates
  "purity_cancer",
  "CPE",
  "ABSOLUTE",
  "LUMP",
  "ESTIMATE",
  "IHC",
  "PANCAN_ABSOLUTE",
  "paper_PURITY..ABSOLUTE.",
  "paper_Purity",
  "paper_Purity..FACETS."
)

selected_columns <- candidate_columns[candidate_columns %in% names(dat)]
analysis <- dat[, selected_columns, drop = FALSE]

to_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

first_nonmissing <- function(...) {
  values <- list(...)
  out <- rep(NA_character_, length(values[[1L]]))

  for (value in values) {
    value <- as.character(value)
    use <- is.na(out) & !is.na(value) & nzchar(trimws(value))
    out[use] <- value[use]
  }

  out
}

get_column <- function(data, column, default = NA_character_) {
  if (column %in% names(data)) {
    data[[column]]
  } else {
    rep(default, nrow(data))
  }
}

# Standardized identifiers ---------------------------------------------------

analysis$patient_id <- first_nonmissing(
  get_column(analysis, "patient"),
  substr(as.character(analysis$sample_id), 1L, 12L)
)

analysis$sample_barcode_standard <- first_nonmissing(
  get_column(analysis, "sample_id"),
  get_column(analysis, "barcode"),
  get_column(analysis, "sample_barcode")
)

# Standardized clinical fields ----------------------------------------------

age_at_index <- to_numeric(get_column(analysis, "age_at_index"))
age_at_diagnosis_days <- to_numeric(get_column(analysis, "age_at_diagnosis"))

analysis$age_years <- age_at_index
use_diagnosis_age <- is.na(analysis$age_years) & !is.na(age_at_diagnosis_days)
analysis$age_years[use_diagnosis_age] <-
  age_at_diagnosis_days[use_diagnosis_age] / 365.25

analysis$sex <- tolower(trimws(as.character(get_column(analysis, "sex_at_birth"))))
analysis$sex[!analysis$sex %in% c("female", "male")] <- NA_character_

analysis$sample_context <- first_nonmissing(
  get_column(analysis, "paper_ALL_PRIMARY_VS_METASTATIC"),
  get_column(analysis, "paper_CURATED_TCGA_SPECIMEN_SITE"),
  get_column(analysis, "requested_sample_type"),
  get_column(analysis, "sample_type"),
  get_column(analysis, "definition")
)

analysis$stage_primary <- first_nonmissing(
  get_column(analysis, "paper_CURATED_PATHOLOGIC_STAGE_AJCC7_AT_DIAGNOSIS_SIMPLE"),
  get_column(analysis, "ajcc_pathologic_stage"),
  get_column(analysis, "ajcc_clinical_stage")
)

# Preliminary GDC overall-survival variables --------------------------------

vital_status <- tolower(trimws(as.character(get_column(analysis, "vital_status"))))
days_to_death <- to_numeric(get_column(analysis, "days_to_death"))
days_to_last_follow_up <- to_numeric(get_column(analysis, "days_to_last_follow_up"))

analysis$os_event_gdc <- ifelse(
  vital_status == "dead",
  1L,
  ifelse(vital_status == "alive", 0L, NA_integer_)
)

analysis$os_time_days_gdc <- ifelse(
  analysis$os_event_gdc == 1L,
  days_to_death,
  days_to_last_follow_up
)

analysis$os_time_years_gdc <- analysis$os_time_days_gdc / 365.25

# Retain a separate curated SKCM endpoint candidate without silently mixing it
# with the GDC endpoint. The endpoint source will be chosen during survival work.
analysis$skcm_curated_followup_days <- to_numeric(
  get_column(analysis, "paper_CURATED_TCGA_DAYS_TO_DEATH_OR_LAST_FU")
)

analysis$skcm_melanoma_specific_event <- to_numeric(
  get_column(
    analysis,
    "paper_CURATED_MELANOMA_SPECIFIC_VITAL_STATUS..0....ALIVE.OR.CENSORED...1....DEAD.OF.MELANOMA.."
  )
)

skcm_curated_vital_status <- tolower(
  trimws(as.character(get_column(analysis, "paper_CURATED_VITAL_STATUS")))
)

analysis$skcm_curated_os_event <- ifelse(
  skcm_curated_vital_status == "dead",
  1L,
  ifelse(skcm_curated_vital_status == "alive", 0L, NA_integer_)
)

analysis$skcm_curated_os_time_days <- to_numeric(
  get_column(analysis, "paper_CURATED_DAYS_TO_DEATH_OR_LAST_FU")
)

analysis$uvm_death_metastasis_category <- as.character(
  get_column(analysis, "paper_Death..Metastasis")
)

# This is a cross-sectional outcome category, not a time-to-metastasis endpoint.
# "Death, other" and "Death, unknown" are not counted as UM-metastasis events.
analysis$uvm_um_metastasis_or_death_event <- ifelse(
  analysis$uvm_death_metastasis_category %in%
    c("Alive, with UM metastasis", "Death, metastatic UM"),
  1L,
  ifelse(
    analysis$uvm_death_metastasis_category == "Alive, no UM metastasis",
    0L,
    NA_integer_
  )
)

# Standardized deconvolution and purity fields -------------------------------

analysis$cibersortx_pass <-
  to_numeric(analysis$CIBERSORTx_P_value) <= 0.05

analysis$extreme_quartile <- analysis$CD43_group %in% c("CD43_low", "CD43_high")

analysis$tumor_purity_cpe <- to_numeric(get_column(analysis, "CPE"))
analysis$tumor_purity_pancan_absolute <- to_numeric(
  get_column(analysis, "PANCAN_ABSOLUTE")
)
analysis$tumor_purity_primary <- analysis$tumor_purity_pancan_absolute
analysis$tumor_purity_primary_method <- ifelse(
  is.na(analysis$tumor_purity_primary),
  NA_character_,
  "PANCAN_ABSOLUTE"
)

# Place standardized variables at the beginning -----------------------------

standard_columns <- c(
  "sample_barcode_standard",
  "patient_id",
  "cancer",
  "project",
  "SPN_TPM",
  "SPN_log2_TPM_plus_1",
  "CD43_group",
  "extreme_quartile",
  "T cells CD8",
  "cibersortx_pass",
  "CIBERSORTx_P_value",
  "CIBERSORTx_Correlation",
  "CIBERSORTx_RMSE",
  "tumor_purity_primary",
  "tumor_purity_primary_method",
  "tumor_purity_pancan_absolute",
  "tumor_purity_cpe",
  "sample_context",
  "age_years",
  "sex",
  "stage_primary",
  "os_event_gdc",
  "os_time_days_gdc",
  "os_time_years_gdc",
  "skcm_curated_followup_days",
  "skcm_melanoma_specific_event",
  "skcm_curated_os_event",
  "skcm_curated_os_time_days",
  "uvm_death_metastasis_category",
  "uvm_um_metastasis_or_death_event"
)

standard_columns <- standard_columns[standard_columns %in% names(analysis)]
analysis <- analysis[, c(standard_columns, setdiff(names(analysis), standard_columns)), drop = FALSE]

# Validation -----------------------------------------------------------------

duplicate_sample_ids <- duplicated(analysis$sample_barcode_standard) |
  duplicated(analysis$sample_barcode_standard, fromLast = TRUE)

if (any(duplicate_sample_ids)) {
  duplicate_values <- unique(analysis$sample_barcode_standard[duplicate_sample_ids])
  stop(
    "Duplicate sample identifiers detected: ",
    paste(head(duplicate_values, 10L), collapse = ", ")
  )
}

expected_counts <- c(SKCM = 103L, UVM = 80L)
observed_counts <- table(analysis$cancer)

for (cancer in names(expected_counts)) {
  observed <- if (cancer %in% names(observed_counts)) observed_counts[[cancer]] else 0L
  if (observed != expected_counts[[cancer]]) {
    warning(
      cancer,
      ": expected ", expected_counts[[cancer]],
      " rows based on the previous analysis table, but found ", observed, "."
    )
  }
}

if (any(to_numeric(analysis$SPN_TPM) < 0, na.rm = TRUE)) {
  stop("Negative SPN TPM values were detected.")
}

if (any(analysis$tumor_purity_primary < 0 | analysis$tumor_purity_primary > 1, na.rm = TRUE)) {
  warning("Some PANCAN_ABSOLUTE tumor-purity values fall outside [0, 1].")
}

# Write outputs ---------------------------------------------------------------

write_tsv_gz <- function(data, path) {
  connection <- gzfile(path, open = "wt")
  on.exit(close(connection), add = TRUE)
  write.table(
    data,
    file = connection,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
}

combined_path <- file.path(
  data_output_dir,
  "07_tcga_analysis_table.tsv.gz"
)
write_tsv_gz(analysis, combined_path)

for (cancer in c("SKCM", "UVM")) {
  cancer_path <- file.path(
    data_output_dir,
    paste0("07_", cancer, "_analysis_table.tsv.gz")
  )
  write_tsv_gz(analysis[analysis$cancer == cancer, , drop = FALSE], cancer_path)
}

qc_summary <- do.call(
  rbind,
  lapply(c("SKCM", "UVM"), function(cancer) {
    x <- analysis[analysis$cancer == cancer, , drop = FALSE]

    data.frame(
      cancer = cancer,
      n_samples = nrow(x),
      n_subjects = length(unique(x$patient_id)),
      n_duplicate_sample_ids = sum(duplicated(x$sample_barcode_standard)),
      n_cd43_low = sum(x$CD43_group == "CD43_low", na.rm = TRUE),
      n_intermediate = sum(x$CD43_group == "intermediate", na.rm = TRUE),
      n_cd43_high = sum(x$CD43_group == "CD43_high", na.rm = TRUE),
      n_cibersortx_pass = sum(x$cibersortx_pass, na.rm = TRUE),
      n_with_primary_purity = sum(!is.na(x$tumor_purity_primary)),
      n_with_age = sum(!is.na(x$age_years)),
      n_with_sex = sum(!is.na(x$sex)),
      n_with_stage = sum(!is.na(x$stage_primary)),
      n_with_os_event = sum(!is.na(x$os_event_gdc)),
      n_with_os_time = sum(!is.na(x$os_time_days_gdc)),
      n_deaths_gdc = sum(x$os_event_gdc == 1L, na.rm = TRUE),
      n_with_skcm_curated_os = sum(
        !is.na(x$skcm_curated_os_event) & !is.na(x$skcm_curated_os_time_days)
      ),
      n_with_uvm_outcome_category = sum(!is.na(x$uvm_death_metastasis_category)),
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  qc_summary,
  file.path(results_dir, "07_tcga_analysis_table_qc.csv"),
  row.names = FALSE,
  na = "NA"
)

missingness_variables <- c(
  "SPN_TPM",
  "T cells CD8",
  "CIBERSORTx_P_value",
  "tumor_purity_primary",
  "age_years",
  "sex",
  "stage_primary",
  "sample_context",
  "os_event_gdc",
  "os_time_days_gdc",
  "skcm_curated_os_event",
  "skcm_curated_os_time_days",
  "skcm_melanoma_specific_event",
  "uvm_death_metastasis_category",
  "uvm_um_metastasis_or_death_event"
)

missingness <- do.call(
  rbind,
  lapply(c("SKCM", "UVM"), function(cancer) {
    x <- analysis[analysis$cancer == cancer, , drop = FALSE]

    do.call(
      rbind,
      lapply(missingness_variables, function(variable) {
        n_missing <- sum(is.na(x[[variable]]))
        data.frame(
          cancer = cancer,
          variable = variable,
          n_total = nrow(x),
          n_missing = n_missing,
          percent_missing = 100 * n_missing / nrow(x),
          stringsAsFactors = FALSE
        )
      })
    )
  })
)

write.csv(
  missingness,
  file.path(results_dir, "07_tcga_analysis_table_missingness.csv"),
  row.names = FALSE,
  na = "NA"
)

capture.output(
  sessionInfo(),
  file = file.path(results_dir, "07_sessionInfo.txt")
)

message("Step 1 complete.")
message("Combined table: ", combined_path)
message("Dimensions: ", nrow(analysis), " rows x ", ncol(analysis), " columns")
message("Cancer counts:")
print(table(analysis$cancer, useNA = "ifany"))
message("CD43 groups by cancer:")
print(table(analysis$cancer, analysis$CD43_group, useNA = "ifany"))
message("QC summary:")
print(qc_summary)
