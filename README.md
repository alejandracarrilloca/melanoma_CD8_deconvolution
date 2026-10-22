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
