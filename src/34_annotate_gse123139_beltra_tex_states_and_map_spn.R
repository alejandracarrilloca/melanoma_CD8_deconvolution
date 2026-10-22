#!/usr/bin/env Rscript

# ==============================================================================
# Step 34: Project the Beltra four-stage exhausted-CD8 hierarchy into GSE123139
# ==============================================================================
#
# Biological framework (Beltra et al., Immunity 2020):
#   Tex_prog1 <-> Tex_prog2 -> Tex_int -> Tex_term
#
# Important design decisions
#   - SPN is never used to define exhaustion eligibility, signature scores,
#     state calls, margins, or the developmental order.
#   - The original Ly108(SLAMF6)/CD69 protein gates are not copied literally
#     into sparse scRNA-seq. Instead, multi-gene transcriptional programs are
#     scored.
#   - Two independent, frozen signature systems are used:
#       1. Beltra-conserved mouse-to-human orthologue programs.
#       2. Human-melanoma-adapted programs describing the same biology.
#   - The primary label requires exact agreement between the two systems.
#   - Ordinary memory/non-Tex cells remain an external reference and are not
#     automatically called Tex_prog1 merely because they express TCF7/IL7R.
#   - Ambiguous cells remain ambiguous rather than being forced into a state.
#   - Patients, not cells, are the independent units for SPN inference.
#
# This is a cross-sectional projection of a published developmental framework.
# It does not prove that an individual human cell moved through these states.
# ==============================================================================

options(stringsAsFactors = FALSE)
options(error = function() {
  cat("\n--- STEP 34 TRACEBACK ---\n", file = stderr())
  traceback(20L)
  quit(save = "no", status = 1L, runLast = FALSE)
})

required_packages <- c(
  "data.table", "SingleCellExperiment", "SummarizedExperiment", "S4Vectors",
  "Matrix", "ggplot2"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

project_root <- normalizePath(getwd(), mustWork = TRUE)
analysis_dir <- file.path(
  project_root, "data", "analysis_ready", "scRNA", "GSE123139"
)
results_dir <- file.path(project_root, "results", "34")
figures_dir <- file.path(results_dir, "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

input_candidates <- file.path(analysis_dir, c(
  "33_GSE123139_spn_clonal_expansion_sce.rds",
  "32_GSE123139_spn_continuum_clone_consolidated_sce.rds",
  "31_GSE123139_cd8_continuum_clonal_states_sce.rds",
  "28_GSE123139_tumor_cd8_within_state_spn_biology_sce.rds"
))
input_sce_file <- input_candidates[file.exists(input_candidates)][1L]
if (is.na(input_sce_file)) {
  stop("No completed Step 28/31/32/33 SCE was found in ", analysis_dir)
}
output_sce_file <- file.path(
  analysis_dir, "34_GSE123139_beltra_tex_states_spn_sce.rds"
)

random_seed <- 123139L
bootstrap_iterations <- 5000L
minimum_signature_genes <- 4L
minimum_call_margin <- 0.15
minimum_patient_state_cells <- 10L
minimum_patient_model_cells <- 100L
minimum_detected_exhaustion_genes <- 1L

base_states <- c("CD8_memory_transitional", "CD8_dysfunctional")
tex_states <- c("Tex_prog1", "Tex_prog2", "Tex_int", "Tex_term")
all_state_levels <- c("Memory_nonTex", tex_states, "Ambiguous_Tex")
state_order <- setNames(0:4, c("Memory_nonTex", tex_states))
adjacent_pairs <- list(
  c("Tex_prog2", "Tex_prog1"),
  c("Tex_int", "Tex_prog2"),
  c("Tex_term", "Tex_int")
)

state_colors <- c(
  Memory_nonTex = "#4C92C3",
  Tex_prog1 = "#68B0A8",
  Tex_prog2 = "#D7B44A",
  Tex_int = "#E77C4E",
  Tex_term = "#7651B9",
  Ambiguous_Tex = "#A9A9A9"
)

# Frozen before examining SPN by state. These are conservative one-to-one
# human orthologues of the subset biology described by Beltra et al.
beltra_orthologue_signatures <- list(
  Tex_prog1 = c(
    "TCF7", "MYB", "IL7R", "SELL", "CXCR5", "ICOS", "CD28",
    "SLAMF6", "IRF7", "OAS1", "STAT1"
  ),
  Tex_prog2 = c(
    "TCF7", "SLAMF6", "CXCR5", "CD28", "MKI67", "CCNB1", "CCNB2",
    "CDC20", "TOP2A", "ANXA2", "ITGB7", "ALCAM"
  ),
  Tex_int = c(
    "CX3CR1", "TBX21", "ZEB2", "ID2", "PRDM1", "GZMA", "GZMB",
    "PRF1", "KLRG1", "NKG7", "CTSW"
  ),
  Tex_term = c(
    "PDCD1", "LAG3", "TIGIT", "CD244", "ENTPD1", "CD101", "CD38",
    "TOX", "EOMES", "HAVCR2"
  )
)

# Human-melanoma-adapted signatures describing the same four biological states.
# These emphasize genes reproducibly used for human tumor CD8 annotation and
# genes already audited in this project's SPN-independent Step 27 programs.
human_melanoma_signatures <- list(
  Tex_prog1 = c(
    "TCF7", "SLAMF6", "IL7R", "CCR7", "SELL", "LEF1", "MAL", "LTB",
    "CXCR5", "CD28", "ICOS"
  ),
  Tex_prog2 = c(
    "TCF7", "SLAMF6", "MKI67", "TOP2A", "STMN1", "TYMS", "UBE2C",
    "TUBA1B", "HMGB2", "ANXA2", "ALCAM"
  ),
  Tex_int = c(
    "TBX21", "CX3CR1", "ZEB2", "ID2", "GZMK", "GZMA", "CCL5",
    "PRF1", "GZMB", "NKG7", "CTSW", "FGFBP2"
  ),
  Tex_term = c(
    "TOX", "PDCD1", "HAVCR2", "LAG3", "TIGIT", "ENTPD1", "CXCL13",
    "LAYN", "EOMES", "PRDM1", "CD101"
  )
)

memory_nontex_signature <- c(
  "TCF7", "LEF1", "CCR7", "IL7R", "SELL", "MAL", "LTB", "NOSIP"
)
exhaustion_entry_signature <- c(
  "TOX", "PDCD1", "TIGIT", "LAG3", "HAVCR2", "ENTPD1", "CXCL13",
  "LAYN", "BATF"
)

stopifnot(!"SPN" %in% unique(c(
  unlist(beltra_orthologue_signatures),
  unlist(human_melanoma_signatures), memory_nontex_signature,
  exhaustion_entry_signature
)))

as_flag <- function(x) !is.na(x) & as.logical(x)

zscore_safe <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

safe_wilcox_p <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3L || all(x == 0)) return(NA_real_)
  tryCatch(
    suppressWarnings(stats::wilcox.test(x, mu = 0, exact = FALSE)$p.value),
    error = function(e) NA_real_
  )
}

bootstrap_median_ci <- function(x, seed_offset = 0L) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(lower = NA_real_, upper = NA_real_))
  if (length(x) == 1L) return(c(lower = x, upper = x))
  set.seed(random_seed + as.integer(seed_offset))
  boot <- replicate(
    bootstrap_iterations,
    stats::median(sample(x, length(x), replace = TRUE)
  ))
  setNames(
    as.numeric(stats::quantile(boot, c(0.025, 0.975), names = FALSE)),
    c("lower", "upper")
  )
}

summarize_effect <- function(x, comparison, effect_scale, seed_offset = 0L) {
  x <- x[is.finite(x)]
  ci <- bootstrap_median_ci(x, seed_offset)
  data.table(
    comparison = comparison,
    effect_scale = effect_scale,
    n_patients = as.integer(length(x)),
    median_effect = if (length(x)) as.numeric(stats::median(x)) else NA_real_,
    bootstrap_ci_lower = unname(ci["lower"]),
    bootstrap_ci_upper = unname(ci["upper"]),
    percent_positive = if (length(x)) 100 * mean(x > 0) else NA_real_,
    wilcoxon_signed_rank_p = safe_wilcox_p(x)
  )
}

theme_project <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.18)),
      plot.subtitle = element_text(color = "#5A5A5A"),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )
}

save_png <- function(filename, plot, width = 10, height = 8) {
  ggsave(
    file.path(figures_dir, filename), plot = plot, width = width,
    height = height, units = "in", dpi = 180, bg = "white"
  )
}

message("Loading completed CD8 object: ", basename(input_sce_file))
sce <- readRDS(input_sce_file)
required_assays <- c("counts", "logcounts")
missing_assays <- setdiff(
  required_assays, SummarizedExperiment::assayNames(sce)
)
if (length(missing_assays)) {
  stop("Required assays are missing: ", paste(missing_assays, collapse = ", "))
}

required_metadata <- c(
  "subject_id", "step25_retained_cd8", "step25_base_state",
  "endogenous_n_counts", "step26_spn_detected", "step26_spn_logcounts",
  "tcr_productive_primary", "tcr_clonotype_key"
)
missing_metadata <- setdiff(
  required_metadata, colnames(SummarizedExperiment::colData(sce))
)
if (length(missing_metadata)) {
  stop("Required metadata are missing: ", paste(missing_metadata, collapse = ", "))
}

cell <- as.data.table(SummarizedExperiment::colData(sce))
cell[, cell_id := colnames(sce)]
cell[, subject_id := as.character(subject_id)]
cell[, step25_base_state := as.character(step25_base_state)]
cell[, retained_cd8 := as_flag(step25_retained_cd8)]
cell[, spn_detected := as.integer(as_flag(step26_spn_detected))]
cell[, spn_logcounts := as.numeric(step26_spn_logcounts)]
cell[, log1p_endogenous_umis := log1p(as.numeric(endogenous_n_counts))]
cell[, productive_tcr := as_flag(tcr_productive_primary)]
cell[, clonotype_key := as.character(tcr_clonotype_key)]

retained_index <- which(
  cell$retained_cd8 & cell$step25_base_state %chin% base_states
)
if (!length(retained_index)) stop("No retained Step 25 CD8 cells were found.")
analysis_dt <- copy(cell[retained_index])
logcounts <- SummarizedExperiment::assay(sce, "logcounts")

all_signature_lists <- list(
  Beltra_orthologue = beltra_orthologue_signatures,
  Human_melanoma = human_melanoma_signatures
)

signature_presence <- rbindlist(lapply(names(all_signature_lists), function(source) {
  signatures <- all_signature_lists[[source]]
  rbindlist(lapply(names(signatures), function(state) {
    requested <- unique(signatures[[state]])
    present <- intersect(requested, rownames(sce))
    data.table(
      signature_source = source,
      tex_state = state,
      requested_genes = as.integer(length(requested)),
      present_genes = as.integer(length(present)),
      genes_present = paste(present, collapse = ";"),
      genes_absent = paste(setdiff(requested, present), collapse = ";")
    )
  }))
}))
fwrite(signature_presence,
       file.path(results_dir, "34_signature_gene_presence.csv"))
if (any(signature_presence$present_genes < minimum_signature_genes)) {
  stop("At least one Tex signature has fewer than ", minimum_signature_genes,
       " genes present. Inspect 34_signature_gene_presence.csv.")
}

all_scoring_genes <- unique(c(
  unlist(beltra_orthologue_signatures),
  unlist(human_melanoma_signatures),
  memory_nontex_signature, exhaustion_entry_signature
))
all_scoring_genes <- intersect(all_scoring_genes, rownames(sce))
expression_matrix <- as.matrix(
  logcounts[all_scoring_genes, retained_index, drop = FALSE]
)

# Standardize each gene across retained CD8 cells, clip extreme values, then
# average genes within each predefined program.
gene_means <- rowMeans(expression_matrix)
gene_sds <- apply(expression_matrix, 1L, stats::sd)
gene_sds[!is.finite(gene_sds) | gene_sds == 0] <- 1
expression_z <- sweep(expression_matrix, 1L, gene_means, "-")
expression_z <- sweep(expression_z, 1L, gene_sds, "/")
expression_z[expression_z > 5] <- 5
expression_z[expression_z < -5] <- -5

score_signature <- function(genes) {
  present <- intersect(genes, rownames(expression_z))
  if (length(present) < minimum_signature_genes) {
    return(rep(NA_real_, ncol(expression_z)))
  }
  colMeans(expression_z[present, , drop = FALSE])
}

analysis_dt[, step34_memory_nontex_score := score_signature(memory_nontex_signature)]
analysis_dt[, step34_exhaustion_entry_score :=
              score_signature(exhaustion_entry_signature)]

for (source in names(all_signature_lists)) {
  source_prefix <- if (source == "Beltra_orthologue") "beltra" else "human"
  for (state in tex_states) {
    column_name <- paste0("step34_", source_prefix, "_score_", state)
    analysis_dt[, (column_name) := score_signature(
      all_signature_lists[[source]][[state]]
    )]
  }
}

# Patient-centering prevents patients with globally shifted program scores from
# dominating state calls.
score_columns <- grep("^step34_(beltra|human)_score_", names(analysis_dt),
                      value = TRUE)
for (column in c(score_columns, "step34_memory_nontex_score",
                 "step34_exhaustion_entry_score")) {
  centered <- paste0(column, "_patient_z")
  analysis_dt[, (centered) := zscore_safe(get(column)), by = subject_id]
}

# Exhaustion-lineage eligibility. A cell is eligible if it belongs to the
# inherited dysfunctional base state, or if it has above-patient-average
# exhaustion-entry expression plus at least one detected exhaustion gene.
entry_present <- intersect(exhaustion_entry_signature, rownames(sce))
entry_detected_count <- Matrix::colSums(
  SummarizedExperiment::assay(sce, "counts")[
    entry_present, retained_index, drop = FALSE
  ] > 0
)
analysis_dt[, step34_detected_exhaustion_genes := as.integer(entry_detected_count)]
analysis_dt[, step34_tex_eligible :=
  step25_base_state == "CD8_dysfunctional" |
    (step34_exhaustion_entry_score_patient_z > 0 &
       step34_detected_exhaustion_genes >= minimum_detected_exhaustion_genes)
]

call_from_scores <- function(dt, source_prefix) {
  columns <- paste0(
    "step34_", source_prefix, "_score_", tex_states, "_patient_z"
  )
  score_matrix <- as.matrix(dt[, ..columns])
  winner_index <- max.col(score_matrix, ties.method = "first")
  sorted_scores <- t(apply(score_matrix, 1L, sort, decreasing = TRUE))
  margin <- sorted_scores[, 1L] - sorted_scores[, 2L]
  call <- tex_states[winner_index]
  call[!dt$step34_tex_eligible] <- "Memory_nonTex"
  call[dt$step34_tex_eligible & (!is.finite(margin) |
                                  margin < minimum_call_margin)] <-
    "Ambiguous_Tex"
  list(call = call, margin = margin)
}

beltra_call <- call_from_scores(analysis_dt, "beltra")
human_call <- call_from_scores(analysis_dt, "human")
analysis_dt[, `:=`(
  step34_beltra_state = beltra_call$call,
  step34_beltra_margin = beltra_call$margin,
  step34_human_state = human_call$call,
  step34_human_margin = human_call$margin
)]

analysis_dt[, step34_consensus_state := fifelse(
  step34_beltra_state == "Memory_nonTex" &
    step34_human_state == "Memory_nonTex",
  "Memory_nonTex",
  fifelse(
    step34_beltra_state %chin% tex_states &
      step34_human_state == step34_beltra_state,
    step34_beltra_state,
    "Ambiguous_Tex"
  )
)]
analysis_dt[, step34_consensus_confidence := pmin(
  step34_beltra_margin, step34_human_margin
)]
analysis_dt[, step34_consensus_order :=
  as.numeric(state_order[step34_consensus_state])]

annotation_concordance <- analysis_dt[, .(
  n_cells = as.integer(.N),
  percent_retained_cd8 = 100 * .N / nrow(analysis_dt),
  n_patients = uniqueN(subject_id)
), by = .(step34_beltra_state, step34_human_state,
          step34_consensus_state)][order(-n_cells)]

state_summary <- analysis_dt[, .(
  n_cells = as.integer(.N),
  percent_retained_cd8 = 100 * .N / nrow(analysis_dt),
  n_patients = uniqueN(subject_id),
  n_patients_with_at_least_10_cells = sum(table(subject_id) >= 10L),
  median_call_confidence = as.numeric(median(
    step34_consensus_confidence, na.rm = TRUE
  )),
  spn_detected_cells = as.integer(sum(spn_detected == 1L)),
  spn_detected_percent = 100 * mean(spn_detected == 1L),
  mean_spn_logcounts = mean(spn_logcounts),
  median_endogenous_umis = as.numeric(median(endogenous_n_counts))
), by = step34_consensus_state][
  order(match(step34_consensus_state, all_state_levels))
]

# Biological validation: expected programs should peak in their corresponding
# state. SPN is deliberately absent from these summaries.
validation_genes <- unique(c(
  "TCF7", "SLAMF6", "CD69", "MKI67", "TOP2A", "TBX21", "CX3CR1",
  "GZMB", "PRF1", "TOX", "EOMES", "PDCD1", "HAVCR2", "ENTPD1",
  "CXCL13"
))
validation_genes <- intersect(validation_genes, rownames(sce))
validation_expression <- as.matrix(
  logcounts[validation_genes, retained_index, drop = FALSE]
)
validation_long <- rbindlist(lapply(validation_genes, function(gene) {
  data.table(
    subject_id = analysis_dt$subject_id,
    step34_consensus_state = analysis_dt$step34_consensus_state,
    marker = gene,
    expression = as.numeric(validation_expression[gene, ])
  )
}))
validation_patient <- validation_long[
  step34_consensus_state %chin% c("Memory_nonTex", tex_states),
  .(n_cells = as.integer(.N), mean_expression = mean(expression)),
  by = .(subject_id, step34_consensus_state, marker)
][n_cells >= minimum_patient_state_cells]
validation_summary <- validation_patient[, .(
  n_patients = as.integer(.N),
  median_patient_mean_expression = as.numeric(median(mean_expression)),
  q25_patient_mean_expression = as.numeric(quantile(mean_expression, 0.25)),
  q75_patient_mean_expression = as.numeric(quantile(mean_expression, 0.75))
), by = .(step34_consensus_state, marker)]

# Patient-balanced SPN summaries.
patient_state_spn <- analysis_dt[
  step34_consensus_state %chin% c("Memory_nonTex", tex_states), .(
    n_cells = as.integer(.N),
    spn_detected_percent = 100 * mean(spn_detected == 1L),
    mean_spn_logcounts = mean(spn_logcounts),
    median_endogenous_umis = as.numeric(median(endogenous_n_counts))
  ), by = .(subject_id, step34_consensus_state)
][n_cells >= minimum_patient_state_cells]

patient_balanced_spn <- patient_state_spn[, .(
  n_patients = as.integer(.N),
  median_patient_cells = as.numeric(median(n_cells)),
  median_patient_spn_detected_percent = as.numeric(
    median(spn_detected_percent)
  ),
  q25_patient_spn_detected_percent = as.numeric(
    quantile(spn_detected_percent, 0.25)
  ),
  q75_patient_spn_detected_percent = as.numeric(
    quantile(spn_detected_percent, 0.75)
  ),
  median_patient_mean_spn_logcounts = as.numeric(
    median(mean_spn_logcounts)
  )
), by = step34_consensus_state][
  order(match(step34_consensus_state, c("Memory_nonTex", tex_states)))
]

# Adjacent state contrasts use paired patients only.
adjacent_contrast_rows <- list()
for (i in seq_along(adjacent_pairs)) {
  pair <- adjacent_pairs[[i]]
  higher <- pair[1L]
  lower <- pair[2L]
  temp <- patient_state_spn[
    step34_consensus_state %chin% c(higher, lower)
  ]
  detection_wide <- dcast(
    temp, subject_id ~ step34_consensus_state,
    value.var = "spn_detected_percent"
  )
  expression_wide <- dcast(
    temp, subject_id ~ step34_consensus_state,
    value.var = "mean_spn_logcounts"
  )
  if (all(c(higher, lower) %in% names(detection_wide))) {
    adjacent_contrast_rows[[length(adjacent_contrast_rows) + 1L]] <-
      summarize_effect(
        detection_wide[[higher]] - detection_wide[[lower]],
        paste0(higher, " minus ", lower, ": SPN detection"),
        "SPN-detection percentage points", 100L + i
      )
  }
  if (all(c(higher, lower) %in% names(expression_wide))) {
    adjacent_contrast_rows[[length(adjacent_contrast_rows) + 1L]] <-
      summarize_effect(
        expression_wide[[higher]] - expression_wide[[lower]],
        paste0(higher, " minus ", lower, ": SPN expression"),
        "Mean SPN logcounts", 200L + i
      )
  }
}
adjacent_contrasts <- rbindlist(adjacent_contrast_rows, fill = TRUE)
adjacent_contrasts[, fdr := p.adjust(wilcoxon_signed_rank_p, method = "BH")]

# Patient-specific ordinal trend models across the four Tex states. Memory is
# excluded: this directly asks whether SPN changes within the Tex hierarchy.
tex_model_dt <- analysis_dt[
  step34_consensus_state %chin% tex_states &
    is.finite(step34_consensus_order) & is.finite(log1p_endogenous_umis)
]
fit_patient_trend <- function(patient_dt) {
  if (nrow(patient_dt) < minimum_patient_model_cells ||
      uniqueN(patient_dt$step34_consensus_state) < 3L ||
      uniqueN(patient_dt$spn_detected) < 2L) {
    return(data.table(
      n_cells = as.integer(nrow(patient_dt)), n_states = as.integer(
        uniqueN(patient_dt$step34_consensus_state)
      ), detection_log_odds = NA_real_, expression_coefficient = NA_real_
    ))
  }
  patient_dt[, depth_z__ := zscore_safe(log1p_endogenous_umis)]
  patient_dt[, order_z__ := zscore_safe(step34_consensus_order)]
  detection_fit <- tryCatch(
    suppressWarnings(glm(
      spn_detected ~ order_z__ + depth_z__, data = patient_dt,
      family = quasibinomial()
    )), error = function(e) NULL
  )
  expression_fit <- tryCatch(
    lm(spn_logcounts ~ order_z__ + depth_z__, data = patient_dt),
    error = function(e) NULL
  )
  data.table(
    n_cells = as.integer(nrow(patient_dt)),
    n_states = as.integer(uniqueN(patient_dt$step34_consensus_state)),
    detection_log_odds = if (is.null(detection_fit)) NA_real_ else
      as.numeric(coef(detection_fit)["order_z__"]),
    expression_coefficient = if (is.null(expression_fit)) NA_real_ else
      as.numeric(coef(expression_fit)["order_z__"])
  )
}

patient_trend_models <- tex_model_dt[, fit_patient_trend(copy(.SD)),
                                     by = subject_id]
trend_summary <- rbindlist(list(
  summarize_effect(
    patient_trend_models$detection_log_odds,
    "Ordered Beltra Tex state versus SPN detection; depth adjusted",
    "SPN-detection log-odds per patient-SD state order", 301L
  ),
  summarize_effect(
    patient_trend_models$expression_coefficient,
    "Ordered Beltra Tex state versus SPN expression; depth adjusted",
    "SPN logcounts per patient-SD state order", 302L
  )
))
trend_summary[, fdr := p.adjust(wilcoxon_signed_rank_p, method = "BH")]

# Productive-clonotype state sharing provides connectivity support, not temporal
# proof. Clones are always patient-specific.
clone_cells <- analysis_dt[
  productive_tcr == TRUE & !is.na(clonotype_key) & nzchar(clonotype_key) &
    step34_consensus_state %chin% tex_states
]
clone_sharing <- if (nrow(clone_cells)) {
  clone_cells[, .(
    clone_cells = as.integer(.N),
    n_tex_states = as.integer(uniqueN(step34_consensus_state)),
    tex_states_occupied = paste(
      tex_states[tex_states %chin% unique(step34_consensus_state)],
      collapse = ";"
    ),
    minimum_state_order = min(step34_consensus_order),
    maximum_state_order = max(step34_consensus_order),
    state_order_span = max(step34_consensus_order) - min(step34_consensus_order)
  ), by = .(subject_id, clonotype_key)]
} else {
  data.table(
    subject_id = character(), clonotype_key = character(),
    clone_cells = integer(), n_tex_states = integer(),
    tex_states_occupied = character(), minimum_state_order = numeric(),
    maximum_state_order = numeric(), state_order_span = numeric()
  )
}
clone_sharing_summary <- data.table(
  productive_tex_clonotypes = as.integer(nrow(clone_sharing)),
  productive_tex_cells = as.integer(nrow(clone_cells)),
  clones_spanning_at_least_two_tex_states = as.integer(
    sum(clone_sharing$n_tex_states >= 2L)
  ),
  percent_clones_spanning_at_least_two_tex_states = if (nrow(clone_sharing))
    100 * mean(clone_sharing$n_tex_states >= 2L) else NA_real_,
  clones_spanning_prog1_to_term = as.integer(sum(
    grepl("Tex_prog1", clone_sharing$tex_states_occupied, fixed = TRUE) &
      grepl("Tex_term", clone_sharing$tex_states_occupied, fixed = TRUE)
  ))
)

# Write Step 34 annotations back to the complete SCE.
metadata_columns <- grep("^step34_", names(analysis_dt), value = TRUE)
for (column in metadata_columns) {
  output <- rep(NA, ncol(sce))
  if (is.logical(analysis_dt[[column]])) output <- rep(NA, ncol(sce))
  output[retained_index] <- analysis_dt[[column]]
  SummarizedExperiment::colData(sce)[[column]] <- output
}

annotation_accounting_validation <-
  nrow(analysis_dt) == length(retained_index) &&
  all(analysis_dt$step34_consensus_state %chin% all_state_levels) &&
  !anyNA(analysis_dt$step34_consensus_state)
spn_independent_annotation_validation <-
  !any(grepl("SPN", unique(c(
    unlist(beltra_orthologue_signatures),
    unlist(human_melanoma_signatures),
    memory_nontex_signature, exhaustion_entry_signature
  )), ignore.case = FALSE))

dataset_qc <- data.table(
  metric = c(
    "retained_cd8_cells", "subjects", "tex_eligible_cells",
    "memory_nontex_consensus_cells", "tex_prog1_consensus_cells",
    "tex_prog2_consensus_cells", "tex_int_consensus_cells",
    "tex_term_consensus_cells", "ambiguous_tex_cells",
    "exact_annotation_agreement_percent", "patients_with_trend_models",
    "productive_tex_clonotypes", "annotation_accounting_validation",
    "spn_independent_annotation_validation", "signature_exported",
    "output_sce_size_mb"
  ),
  value = c(
    nrow(analysis_dt), uniqueN(analysis_dt$subject_id),
    sum(analysis_dt$step34_tex_eligible),
    sum(analysis_dt$step34_consensus_state == "Memory_nonTex"),
    sum(analysis_dt$step34_consensus_state == "Tex_prog1"),
    sum(analysis_dt$step34_consensus_state == "Tex_prog2"),
    sum(analysis_dt$step34_consensus_state == "Tex_int"),
    sum(analysis_dt$step34_consensus_state == "Tex_term"),
    sum(analysis_dt$step34_consensus_state == "Ambiguous_Tex"),
    100 * mean(analysis_dt$step34_beltra_state ==
                 analysis_dt$step34_human_state),
    sum(is.finite(patient_trend_models$detection_log_odds)),
    nrow(clone_sharing), as.integer(annotation_accounting_validation),
    as.integer(spn_independent_annotation_validation), 0L, NA_real_
  )
)

message("Writing Step 34 tables")
fwrite(dataset_qc, file.path(results_dir, "34_dataset_qc.csv"))
fwrite(annotation_concordance,
       file.path(results_dir, "34_annotation_concordance.csv"))
fwrite(state_summary, file.path(results_dir, "34_consensus_state_summary.csv"))
fwrite(validation_summary,
       file.path(results_dir, "34_marker_validation_summary.csv"))
fwrite(patient_state_spn,
       file.path(results_dir, "34_patient_state_spn.csv"))
fwrite(patient_balanced_spn,
       file.path(results_dir, "34_patient_balanced_spn_summary.csv"))
fwrite(adjacent_contrasts,
       file.path(results_dir, "34_adjacent_state_spn_contrasts.csv"))
fwrite(patient_trend_models,
       file.path(results_dir, "34_patient_ordinal_trend_models.csv"))
fwrite(trend_summary,
       file.path(results_dir, "34_ordinal_trend_summary.csv"))
fwrite(clone_sharing,
       file.path(results_dir, "34_clonotype_tex_state_sharing.tsv.gz"),
       sep = "\t")
fwrite(clone_sharing_summary,
       file.path(results_dir, "34_clonotype_tex_state_sharing_summary.csv"))

message("Rendering Step 34 figures")
if ("UMAP_step21" %in% SingleCellExperiment::reducedDimNames(sce)) {
  umap <- SingleCellExperiment::reducedDim(sce, "UMAP_step21")
  plot_dt <- data.table(
    UMAP1 = umap[retained_index, 1L], UMAP2 = umap[retained_index, 2L],
    state = factor(analysis_dt$step34_consensus_state,
                   levels = all_state_levels)
  )
  set.seed(random_seed)
  plot_dt <- plot_dt[sample(.N)]
  p_umap <- ggplot(plot_dt, aes(UMAP1, UMAP2, color = state)) +
    geom_point(size = 0.55, alpha = 0.72) +
    scale_color_manual(values = state_colors, drop = FALSE) +
    labs(
      title = "Beltra exhaustion-state projection in tumor CD8 cells",
      subtitle = "Consensus of conserved-orthologue and human-melanoma signatures; SPN excluded",
      color = "Consensus state"
    ) + theme_project()
  save_png("34A_beltra_consensus_states_umap.png", p_umap, 12, 8)
}

plot_spn <- copy(patient_state_spn)
plot_spn[, step34_consensus_state := factor(
  step34_consensus_state, levels = c("Memory_nonTex", tex_states)
)]
p_spn <- ggplot(
  plot_spn,
  aes(step34_consensus_state, spn_detected_percent,
      group = subject_id)
) +
  geom_line(color = "#B8B8B8", alpha = 0.45) +
  geom_point(aes(color = step34_consensus_state), size = 2, alpha = 0.80) +
  stat_summary(fun = median, geom = "line", aes(group = 1),
               color = "#222222", linewidth = 1.2) +
  stat_summary(fun = median, geom = "point", color = "#222222", size = 3) +
  scale_color_manual(values = state_colors, guide = "none") +
  labs(
    title = "Patient-balanced SPN detection across Beltra CD8 states",
    subtitle = "Lines connect states observed with at least 10 cells in the same patient",
    x = NULL, y = "SPN-detected cells (%)"
  ) + theme_project() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_png("34B_patient_spn_across_beltra_states.png", p_spn, 11, 8)

heatmap_dt <- copy(validation_summary)
heatmap_dt[, step34_consensus_state := factor(
  step34_consensus_state, levels = c("Memory_nonTex", tex_states)
)]
heatmap_dt[, marker_z := zscore_safe(median_patient_mean_expression),
           by = marker]
p_heatmap <- ggplot(
  heatmap_dt,
  aes(step34_consensus_state, marker, fill = marker_z)
) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "#2C7BB6", mid = "white", high = "#7B3294", midpoint = 0,
    name = "Patient-balanced\nmarker z-score"
  ) +
  labs(
    title = "Biological validation of projected Beltra states",
    subtitle = "Expected marker ordering; SPN is not included",
    x = NULL, y = NULL
  ) + theme_project() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_png("34C_beltra_state_marker_validation.png", p_heatmap, 11, 8)

message("Saving the Step 34 annotated object")
saveRDS(sce, output_sce_file, compress = FALSE)
output_size_mb <- file.info(output_sce_file)$size / 1024^2
dataset_qc[metric == "output_sce_size_mb", value := output_size_mb]
fwrite(dataset_qc, file.path(results_dir, "34_dataset_qc.csv"))

if (!annotation_accounting_validation) {
  stop("Step 34 annotation accounting validation failed.")
}
if (!spn_independent_annotation_validation) {
  stop("SPN entered a Step 34 annotation signature.")
}

message("Step 34 complete: Beltra Tex hierarchy projected and SPN mapped")
