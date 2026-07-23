# cellNexus: Quality control, annotation, aggregation and analytical layers for the Human Cell Atlas data

This repository contains all analysis code and workflows accompanying the cellNexus article. The article describes a harmonised resource of over 44 million human single-cell RNA-seq cells drawn from [CELLxGENE Census](https://chanzuckerberg.github.io/cellxgene-census/), processed through standardised quality control, consistent normalisation, and unified abundance layers (raw counts, counts per million, normalised expression, and pseudobulk) to enable reliable cross-dataset analyses. Each subdirectory below corresponds to a distinct analytical component of the paper and contains the scripts or notebooks needed to reproduce the associated figures and results.

---

## Summary

The cellNexus resource aggregates data at scale so that metadata, QC flags, and expression representations are comparable across studies and tissues. This repository documents exactly how those representations were built and evaluated:

- **Metadata pipeline** — how raw CELLxGENE Census data were ingested, quality-controlled, and unified.
- **Quality control** — tissue-level QC landscapes and density-based diagnostics used in the paper.
- **Cell typing** — the immune-graph consensus mapping methodology and results.
- **Doublet-aware DE** — differential expression accounting for doublet contamination.
- **Sex prediction** — a classifier that predicts donor sex from expression data.
- **Ethnicity prediction** — harmonisation and benchmarking of ethnicity labels.
- **Cell–cell communication** — ligand–receptor interaction analyses across tissues.

---

## Repository layout

| Path | Purpose in the article |
|------|------------------------|
| `metadata/` | **Reference pipeline — not intended to be re-run locally.** Stepwise R scripts (`step1`–`step9`) that built the cellNexus resource from scratch: downloading CELLxGENE Census data from AWS, extracting and harmonising metadata with **cellxgenedp**, splitting large samples, running **HPCell**-based QC at scale, unifying metadata, and preparing single-cell and pseudobulk caches (local and cloud). These scripts required large-scale high-performance computing (HPC) with substantial parallelisation and storage resources; they are provided for transparency and documentation, not for direct reproduction. |
| `quality_control/` | Scripts and a Quarto report (`tissue_landscape.Qmd`) that produce the tissue-level QC landscape figures in the paper, including technology-by-tissue bubble plots, cell/sample count scatterplots, library-size density plots, age distributions, and diagnostics for samples expression distribution. |
| `cell_typing/` | Quarto notebooks documenting the immune-graph consensus cell-typing methodology (`cell_typing_method.qmd`), summary statistics of the mapping results (`cell_typing_results.qmd`), a memory B-cell case study showing how consensus annotation improves pseudobulk analysis (`B_memory_case_study.qmd`), and an alternative immune graph built from HCAO (`HCAO_tree.qmd`). |
| `doublet_DE/` | Quarto notebook and helper scripts for doublet-aware differential expression analysis, demonstrating the impact of doublet contamination on DE results. |
| `sex_prediction/` | R Markdown analysis that trains and evaluates a classifier for predicting donor sex from single-cell expression data. |
| `ethnicity_prediction/` | R/Quarto analysis (`Ethnicity_CellNexus.qmd`) and Python classifier benchmark notebook (`Classifier_benchmark.ipynb`) for harmonising ethnicity labels and benchmarking prediction models; fixed processed inputs are available from Zenodo. |
| `communication/` | R Markdown analysis (`communication_analyses.Rmd`) of cell–cell communication via ligand–receptor interactions across tissues using the cellNexus data. |
| `annotate_harmony_clusters.R` | Shared utility for annotating Harmony-derived clusters used across analyses. |
| `tissue_color_utils.R`, `tissue_*.csv` | Shared colour palettes and tissue label → macrocategory mappings used consistently across all figures. |

---

## Software

**R** (4.5 recommended) with [BiocManager](https://bioconductor.org/) for Bioconductor packages.

Install the development package as needed, for example:
```r
# Install BiocManager if not already installed
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("cellNexus")
```
OR 
```r
remotes::install_github("MangiolaLaboratory/cellNexus")
```

---

## Data

Source observational data are drawn from the [**CELLxGENE Census**](https://chanzuckerberg.github.io/cellxgene-census/) (see manuscript for the Census snapshot / policy). Processed layers and metadata consumed by **cellNexus** are accessed via the **cellNexus** R package and paths documented in the paper’s Data availability statement—replace local paths in this repo with those resources.

---

## Licence

Code in this repository is licensed under the **GNU General Public License v3.0** — see [`LICENSE`](LICENSE).

---

## Contact

Documentation for the **cellNexus** R and Python API is on the pkgdown site at [**cellnexus.org**](https://cellnexus.org/). For bug reports and package development, use issues or discussions on the main [**cellNexus** GitHub repository](https://github.com/MangiolaLaboratory/cellNexus) where appropriate.
