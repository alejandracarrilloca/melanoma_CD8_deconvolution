#!/usr/bin/env Rscript

# ==============================================================================
# Create and lock the project-local R package environment
# ==============================================================================
#
# This script initializes renv in the current repository and installs the
# packages required for:
#   - TCGA download and SummarizedExperiment processing
#   - CIBERSORTx input preparation and downstream statistics
#   - publication-quality plotting
#   - SKCM/UVM survival analysis
#   - melanoma scRNA-seq processing, annotation and pseudobulk analysis
#   - pathway and gene-set analysis
#
# Run once from the repository root:
#   Rscript src/00_setup_project_r_environment.R
#
# Recreate the same environment later with:
#   Rscript -e 'renv::restore(prompt = FALSE)'
#
# Outputs:
#   renv.lock
#   renv/
#   .Rprofile
#   results/00/00_required_r_packages.csv
#   results/00/00_installed_r_package_versions.csv
#   results/00/00_r_environment_summary.txt
#   results/00/00_sessionInfo.txt
#
# The renv package library is project-local. Commit renv.lock, .Rprofile and
# renv/activate.R to version control, but do not commit renv/library/.
# ==============================================================================

options(stringsAsFactors = FALSE)
options(timeout = max(3600, getOption("timeout")))

# Prevent renv from rewriting the source repository into a platform-specific
# Posit Package Manager binary URL. The server identifies itself as CentOS 9,
# for which the historical binary endpoint is not available.
Sys.setenv(RENV_CONFIG_PPM_ENABLED = "FALSE")

# R 4.3.3 is paired with Bioconductor 3.18. Use a dated CRAN snapshot from
# the R 4.3 era so that a future rerun does not resolve packages requiring
# newer R headers or a newer Linux toolchain.
cran_snapshot_date <- "2024-06-28"
cran_repository <- paste0(
  "https://packagemanager.posit.co/cran/",
  cran_snapshot_date
)
options(repos = c(CRAN = cran_repository))

project_root <- normalizePath(getwd(), mustWork = TRUE)
results_directory <- file.path(project_root, "results", "00")
bootstrap_library <- file.path(project_root, ".renv-bootstrap-library")
project_cache <- file.path(project_root, ".renv-cache")

# The default renv cache is placed under the user's home directory. Large
# Bioconductor experiment packages exceeded the writable home-cache capacity
# on this HPC and curl stopped with error code 23. Keep the cache beside the
# project on /export/space3 instead.
Sys.setenv(RENV_PATHS_CACHE = project_cache)

invisible(lapply(
  c(
    file.path(project_root, "src"),
    file.path(project_root, "data"),
    file.path(project_root, "results"),
    results_directory,
    bootstrap_library,
    project_cache
  ),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

# The shared Bioconda R library is read-only. Add a small writable bootstrap
# library inside the repository before attempting to install renv itself.
.libPaths(unique(c(bootstrap_library, .libPaths())))

# Keep the temporary bootstrap library out of version control. renv will create
# and manage its final project library separately under renv/library/.
gitignore_file <- file.path(project_root, ".gitignore")
gitignore_entries <- c(
  ".renv-bootstrap-library/",
  ".renv-cache/"
)
gitignore_lines <- if (file.exists(gitignore_file)) {
  readLines(gitignore_file, warn = FALSE)
} else {
  character()
}

missing_gitignore_entries <- setdiff(gitignore_entries, gitignore_lines)

if (length(missing_gitignore_entries) > 0L) {
  writeLines(
    c(gitignore_lines, missing_gitignore_entries),
    gitignore_file
  )
}

timestamp_message <- function(...) {
  message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ...)
}

if (getRversion() < "4.3.0") {
  stop(
    "R >= 4.3.0 is required. Current version: ",
    as.character(getRversion())
  )
}

# Install sequentially. Parallel renv installation on this shared filesystem
# caused packages to collide in renv/staging/1, producing truncated archives,
# disappearing working directories and dependencies that could not see Rcpp.
options(Ncpus = 1L)
Sys.setenv(MAKEFLAGS = "-j1")

# ------------------------------------------------------------------------------
# Direct project dependencies
# ------------------------------------------------------------------------------

cran_packages <- c(
  # Core data handling and reproducible analysis
  "renv",
  "data.table",
  "dplyr",
  "tidyr",
  "tibble",
  "readr",
  "stringr",
  "purrr",
  "forcats",
  "broom",

  # Plotting
  "ggplot2",
  "scales",
  "ggrepel",
  "patchwork",
  "cowplot",

  # Survival analysis
  "survival",

  # Single-cell analysis
  "Seurat",
  "SeuratObject",
  "harmony",
  "future",
  "future.apply",
  "hdf5r",

  # Gene-set resources
  "msigdbr"
)

bioconductor_packages <- c(
  # TCGA and expression containers
  "TCGAbiolinks",
  "SummarizedExperiment",
  "SingleCellExperiment",
  "Biobase",

  # GEO and annotation resources
  "GEOquery",
  "AnnotationDbi",
  "org.Hs.eg.db",
  "biomaRt",

  # Single-cell QC, normalization and annotation. scater is intentionally not
  # required: on this HPC it pulls ggrastr -> Cairo/ragg/textshaping, which need
  # unavailable system development headers. Seurat, scuttle and scran provide
  # the QC, normalization and visualization functionality used by this project.
  "scran",
  "scuttle",
  "SingleR",
  "celldex",
  "DropletUtils",
  "glmGamPoi",
  "MAST",
  "UCell",

  # Differential-expression and pseudobulk analysis
  "edgeR",
  "limma",
  "DESeq2",

  # Gene-set and pathway analysis
  "GSVA",
  "fgsea",
  "clusterProfiler",

  # Heatmaps and scalable computation
  "ComplexHeatmap",
  "BiocParallel",
  "BiocNeighbors",
  "BiocSingular"
)

package_manifest <- rbind(
  data.frame(
    package = cran_packages,
    source = "CRAN",
    stringsAsFactors = FALSE
  ),
  data.frame(
    package = bioconductor_packages,
    source = "Bioconductor",
    stringsAsFactors = FALSE
  )
)

package_manifest <- package_manifest[!duplicated(package_manifest$package), ]

write.csv(
  package_manifest,
  file.path(results_directory, "00_required_r_packages.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# Bootstrap and configure renv
# ------------------------------------------------------------------------------

if (!requireNamespace("renv", quietly = TRUE)) {
  timestamp_message(
    "Installing renv into writable bootstrap library: ",
    bootstrap_library
  )
  install.packages(
    "renv",
    repos = cran_repository,
    lib = bootstrap_library
  )
}

if (!requireNamespace("renv", quietly = TRUE)) {
  stop("renv could not be installed.")
}

renv::consent(provided = TRUE)

activate_file <- file.path(project_root, "renv", "activate.R")
lock_file <- file.path(project_root, "renv.lock")

if (!file.exists(activate_file)) {
  timestamp_message("Initializing a bare project-local renv environment")
  renv::init(
    project = project_root,
    bare = TRUE,
    restart = FALSE
  )
} else {
  timestamp_message("Using the existing renv project infrastructure")
}

# Do not call renv::activate() or renv::load() here. When this setup script was
# launched non-interactively with Rscript --vanilla, activation spent many
# minutes rescanning the project without reaching installation. The project
# argument supplied to renv operations is sufficient; explicitly put its
# library first for requireNamespace() and validation in this process.
project_library <- renv::paths$library(project = project_root)
dir.create(project_library, recursive = TRUE, showWarnings = FALSE)
.libPaths(unique(c(
  project_library,
  setdiff(.libPaths(), bootstrap_library)
)))

# Assert the compatible snapshot for this installation instead of inheriting a
# repository recorded by a partially completed setup attempt.
options(repos = c(CRAN = cran_repository))

if (normalizePath(.libPaths()[1], mustWork = FALSE) !=
    normalizePath(project_library, mustWork = FALSE)) {
  stop(
    "The active installation library is not the project renv library.\n",
    "Expected: ", project_library, "\n",
    "Observed: ", .libPaths()[1]
  )
}

# R 4.3 corresponds to Bioconductor 3.18.
renv::settings$bioconductor.version(
  "3.18",
  project = project_root
)
renv::settings$ppm.enabled(
  FALSE,
  project = project_root
)

# This setup script builds the environment from the manifest below and writes a
# fresh lockfile only after validation succeeds. To reproduce an already
# completed environment, run renv::restore() separately as documented above.

# ------------------------------------------------------------------------------
# Install all declared packages into the project library
# ------------------------------------------------------------------------------

# Pin the two packages most sensitive to the R/toolchain versions on this HPC.
# fs 2.1.0 failed to compile against the server's system headers, while the
# current Seurat release pulled a dependency stack intended for newer R builds.
cran_specs <- c(
  "fs@1.6.4",
  "SeuratObject@5.0.2",
  "Seurat@5.1.0",
  setdiff(cran_packages, c("Seurat", "SeuratObject"))
)
cran_specs <- unique(cran_specs)
bioconductor_specs <- paste0("bioc::", bioconductor_packages)

timestamp_message(
  "Installing or validating ",
  length(cran_specs),
  " direct CRAN dependencies"
)

renv::install(
  packages = cran_specs,
  project = project_root,
  prompt = FALSE
)

timestamp_message(
  "Installing or validating ",
  length(bioconductor_specs),
  " direct Bioconductor dependencies sequentially"
)

for (package_spec in bioconductor_specs) {
  timestamp_message("Bioconductor dependency: ", package_spec)
  renv::install(
    packages = package_spec,
    project = project_root,
    prompt = FALSE
  )
}

# ------------------------------------------------------------------------------
# Validate the completed environment
# ------------------------------------------------------------------------------

all_required_packages <- package_manifest$package

package_available <- vapply(
  all_required_packages,
  requireNamespace,
  logical(1),
  quietly = TRUE
)

version_table <- package_manifest
version_table$installed <- package_available
version_table$version <- vapply(
  version_table$package,
  function(package_name) {
    if (!requireNamespace(package_name, quietly = TRUE)) {
      return(NA_character_)
    }
    as.character(utils::packageVersion(package_name))
  },
  character(1)
)

version_table$library_path <- vapply(
  version_table$package,
  function(package_name) {
    if (!requireNamespace(package_name, quietly = TRUE)) {
      return(NA_character_)
    }
    dirname(dirname(find.package(package_name)))
  },
  character(1)
)

write.csv(
  version_table,
  file.path(results_directory, "00_installed_r_package_versions.csv"),
  row.names = FALSE,
  na = "NA"
)

missing_after_install <- version_table$package[!version_table$installed]

if (length(missing_after_install) > 0L) {
  stop(
    "The following required package(s) remain unavailable: ",
    paste(missing_after_install, collapse = ", "),
    "\nInspect the installation output for missing system dependencies."
  )
}

timestamp_message("Creating the reproducible package lockfile")
renv::snapshot(
  project = project_root,
  type = "all",
  prompt = FALSE
)

bioconductor_version <- tryCatch(
  as.character(renv::bioconductor.version(project = project_root)),
  error = function(e) NA_character_
)

environment_summary <- c(
  paste0("Project root: ", project_root),
  paste0("R version: ", R.version.string),
  paste0("Platform: ", R.version$platform),
  paste0("Bioconductor version: ", bioconductor_version),
  paste0("CRAN snapshot: ", cran_snapshot_date),
  paste0("CRAN repository: ", cran_repository),
  paste0("Bootstrap library: ", bootstrap_library),
  paste0("Project renv cache: ", project_cache),
  paste0("Project library: ", renv::paths$library(project = project_root)),
  paste0("Direct CRAN packages: ", length(cran_packages)),
  paste0("Direct Bioconductor packages: ", length(bioconductor_packages)),
  paste0("All required packages available: ", all(package_available)),
  paste0("Lockfile: ", lock_file)
)

writeLines(
  environment_summary,
  file.path(results_directory, "00_r_environment_summary.txt")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(results_directory, "00_sessionInfo.txt")
)

timestamp_message("Project-local R environment setup is complete")
cat(paste(environment_summary, collapse = "\n"), "\n")
cat(
  "\nTo reproduce this environment on another system, run:\n",
  "Rscript -e 'renv::restore(prompt = FALSE)'\n",
  sep = ""
)