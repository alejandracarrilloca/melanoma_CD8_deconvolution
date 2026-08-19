# CD43-associated immune infiltration across TCGA cancers

This repository investigates the relationship between **SPN/CD43 expression**, immune-cell composition, and CD8 T-cell infiltration in cancer. The current focused hypothesis is that melanoma tumors with higher bulk CD43 expression contain more CD8 T cells and that these CD8 populations may acquire dysfunctional or exhausted states associated with adverse clinical outcomes.

Because CD43 is strongly expressed by leukocytes, the workflow also tests whether the bulk SPN signal reflects immune infiltration rather than tumor-cell expression. TCGA bulk RNA-seq analyses are therefore combined with tumor-purity adjustment and, in the next phase, melanoma single-cell RNA-seq.

## Analysis overview

The pipeline is organized so that every script represents one numbered analysis step:

```text
src/00_*.R  →  results/00/
src/01_*.R  →  results/01/
...
src/09_*.R  →  results/09/
src/10_*.R  →  results/10/
```

Reusable intermediate datasets are stored under `data/`, while statistical summaries, quality-control tables, figures, and session information are stored under `results/<step>/`. Figures are saved as PNG files under `results/<step>/figures/`.

## Cancer cohorts

The bulk analysis includes ten TCGA cancer cohorts:

- GBM: Glioblastoma multiforme
- LAML: Acute myeloid leukemia
- LUAD: Lung adenocarcinoma
- LUSC: Lung squamous cell carcinoma
- PAAD: Pancreatic adenocarcinoma
- SKCM: Skin cutaneous melanoma
- STAD: Stomach adenocarcinoma
- TGCT: Testicular germ cell tumors
- THYM: Thymoma
- UVM: Uveal melanoma

The focused melanoma analyses in Steps 07–10 use TCGA-SKCM and TCGA-UVM.


## Requirements

- R 4.3 or later
- Bioconductor 3.18 for R 4.3
- Access to the NCI Genomic Data Commons
- CIBERSORTx access for the external LM22 deconvolution step
- Sufficient memory and storage for TCGA and single-cell matrices

All commands below must be executed from the repository root.

## Step 00 — Create the project R environment

**Script:** `src/00_setup_project_r_environment.R`

Creates a project-local `renv` environment and installs the CRAN and Bioconductor packages required for the bulk, survival, and single-cell analyses. The R 4.3 configuration uses a compatible CRAN snapshot and Bioconductor 3.18.

```bash
Rscript src/00_setup_project_r_environment.R
```

Main outputs:

```text
renv.lock
renv/
.Rprofile
results/00/00_required_r_packages.csv
results/00/00_installed_r_package_versions.csv
results/00/00_r_environment_summary.txt
results/00/00_sessionInfo.txt
```

To restore the environment later:

```bash
Rscript -e 'renv::restore(prompt = FALSE)'
```

## Step 01 — Download and prepare TCGA RNA-seq data

**Script:** `src/01_query_download.R`

Queries TCGA STAR-counts data through the GDC, selects the requested primary tumor samples, extracts untransformed TPM expression values, and saves sample metadata. Serialized `SummarizedExperiment` objects are retained as checkpoints to avoid repeated downloads.

```bash
Rscript src/01_query_download.R
```

Main reusable outputs:

```text
data/tcga_star_counts/queries/
data/tcga_star_counts/summarized_experiments/
data/tcga_star_counts/tpm_matrices/
data/tcga_star_counts/sample_metadata/
```

Run records are written to:

```text
results/01/
```

## Step 02 — Prepare CIBERSORTx mixtures and CD43 groups

**Script:** `src/02_prepare_cibersortx_inputs.R`

Prepares one genes-by-samples TPM mixture matrix per cancer for CIBERSORTx. Duplicate gene symbols are handled using a documented aggregation rule. Samples are stratified within each cancer using SPN expression quartiles:

- `CD43_low`: SPN TPM at or below the first quartile
- `CD43_high`: SPN TPM at or above the third quartile
- `intermediate`: samples between the two quartiles

```bash
Rscript src/02_prepare_cibersortx_inputs.R
```

Main outputs:

```text
data/cibersortx_input/mixture_files/
data/cibersortx_input/cd43_groups/
results/02/
```

## External CIBERSORTx step

The mixture files from Step 02 are analyzed outside R using CIBERSORTx with the LM22 signature matrix. The returned result files must be organized as:

```text
data/cibersortx_output/<CANCER>/CIBERSORTx_<CANCER>_Adjusted.txt
```

These files are external inputs and are not recreated by Scripts 00–10.

## Step 03 — Validate and join CIBERSORTx results

**Script:** `src/03_validate_join_cibersortx.R`

Validates the 22 LM22 immune-cell fractions, CIBERSORTx fit statistics, sample identifiers, and fraction sums. It then joins the relative immune fractions to the SPN expression groups and TCGA metadata.

```bash
Rscript src/03_validate_join_cibersortx.R
```

Main reusable outputs:

```text
data/analysis_ready/all_cancers_analysis_ready.tsv.gz
data/analysis_ready/cancer_tables/<CANCER>_analysis_ready.tsv.gz
```

QC outputs and the PNG fit-quality figure are written to:

```text
results/03/
results/03/figures/
```

## Step 04 — Test CD43–immune-cell associations

**Script:** `src/04_cd43_inmune_cell_compositions.R`

Tests the association between CD43 and each of the 22 relative LM22 immune-cell fractions using two complementary approaches:

1. Wilcoxon rank-sum tests comparing CD43-high and CD43-low tumors.
2. Spearman correlations using continuous `log2(SPN TPM + 1)`.

Rank-biserial correlations quantify the high-versus-low effect sizes. Benjamini–Hochberg false-discovery-rate correction is applied within each analysis.

```bash
Rscript src/04_cd43_inmune_cell_compositions.R
```

Outputs include per-cancer results, combined cross-cancer tables, effect matrices, and a cross-cancer consistency summary:

```text
results/04/
```

## Step 05 — Visualize cross-cancer immune associations

**Script:** `src/05_plot_cd43_immune_associations.R`

Creates cross-cancer heatmaps and per-cancer butterfly profiles summarizing the direction and magnitude of CD43-associated immune-cell differences. All figures are exported as PNG files.

```bash
Rscript src/05_plot_cd43_immune_associations.R
```

Main outputs:

```text
results/05/05_composition_summary.csv
results/05/05_figure_manifest.csv
results/05/figures/
results/05/figures/per_cancer/
```

## Step 06 — Adjust CD43 associations for tumor purity

**Script:** `src/06_purity_adjusted_cd43_associations.R`

Joins published TCGA tumor-purity estimates to the analysis-ready dataset and tests whether purity explains the observed CD43–immune associations. The purity analysis is restricted to the eight solid-tumor cohorts because conventional tumor purity has a different interpretation in LAML and THYM.

The script evaluates:

1. Tumor-purity differences between CD43-high and CD43-low samples.
2. Continuous correlations between SPN expression and purity.
3. SPN–immune-fraction models adjusted for purity.
4. CD43-high versus CD43-low models adjusted for purity.

Relative immune fractions are transformed using the arcsine square-root transformation. Adjusted and unadjusted models use the same complete-case samples.

```bash
Rscript src/06_purity_adjusted_cd43_associations.R
```

Main reusable output:

```text
data/purity_analysis/all_samples_with_purity.tsv.gz
```

Statistical results and PNG figures are written to:

```text
results/06/
results/06/figures/
```

## Step 07 — Build the focused SKCM/UVM analysis table

**Script:** `src/07_build_tcga_analysis_table.R`

Creates a compact, validated analysis table for TCGA-SKCM and TCGA-UVM. It combines SPN expression, CD43 groups, LM22 fractions, tumor purity, sample context, demographics, tumor stage, and available clinical outcome variables. The script also checks sample uniqueness and reports variable coverage and missingness.

```bash
Rscript src/07_build_tcga_analysis_table.R
```

Main reusable outputs:

```text
data/analysis_ready/skcm_uvm_tables/07_tcga_analysis_table.tsv.gz
data/analysis_ready/skcm_uvm_tables/07_SKCM_analysis_table.tsv.gz
data/analysis_ready/skcm_uvm_tables/07_UVM_analysis_table.tsv.gz
```

QC outputs are written to:

```text
results/07/
```

## Step 08 — Calculate complementary immune-infiltration measurements

**Script:** `src/08_calculate_immune_infiltration_scores.R`

Adds expression-based and composition-based immune measurements for SKCM and UVM:

- PTPRC/CD45 expression
- T-cell marker score: mean log2 expression of `CD3D`, `CD3E`, and `TRAC`
- CD8 marker score: mean log2 expression of `CD8A` and `CD8B`
- CIBERSORTx relative CD8 fraction
- Total relative T-cell fraction
- Total relative lymphocyte fraction
- CD8 fraction within the inferred T-cell compartment
- Tumor purity and inverse tumor purity

CIBERSORTx absolute-mode scores are intentionally excluded because diagnostic inspection showed that they were approximately zero and therefore technically invalid as abundance measurements.

```bash
Rscript src/08_calculate_immune_infiltration_scores.R
```

Main reusable outputs:

```text
data/analysis_ready/skcm_uvm_immune_scores/08_tcga_analysis_with_immune_scores.tsv.gz
data/analysis_ready/skcm_uvm_immune_scores/08_tcga_immune_infiltration_scores.tsv.gz
data/analysis_ready/skcm_uvm_immune_scores/08_<CANCER>_analysis_with_immune_scores.tsv.gz
```

QC, missingness, and summary tables are written to:

```text
results/08/
```

## Step 09 — Test SPN associations with immune infiltration

**Script:** `src/09_test_spn_immune_infiltration_associations.R`

Performs the focused melanoma tests separately in SKCM and UVM:

1. Spearman correlations between continuous SPN expression and immune measurements.
2. CD43-high versus CD43-low Wilcoxon comparisons.
3. Standardized linear models adjusted for tumor purity.
4. CD8 models adjusted for PTPRC expression to test whether the CD8 association extends beyond general leukocyte infiltration.

```bash
Rscript src/09_test_spn_immune_infiltration_associations.R
```

Main outputs:

```text
results/09/09_spn_immune_spearman_results.csv
results/09/09_cd43_high_vs_low_immune_comparisons.csv
results/09/09_spn_adjusted_immune_models.csv
results/09/09_spn_immune_association_qc.csv
results/09/figures/
```

The scatterplots use the project-standard gray and lavender style and are saved only as PNG files.

## Reproducibility notes

- Every analysis step records `sessionInfo()` in its corresponding results directory.
- Benjamini–Hochberg correction is used where multiple immune-cell populations are tested.
- Relative LM22 fractions are compositional and should not be interpreted as absolute cell counts.
- TCGA data and other large matrices should be downloaded or regenerated rather than committed directly to Git.
- External CIBERSORTx outputs, accessions, parameter settings, and logs should be retained for provenance.

