# cellNexus — companion code

This repository contains analysis code and workflows accompanying the cellNexus article: a harmonised integration of large-scale human single-cell data with standardised quality control, consistent normalisation, and unified abundance layers (single-cell counts, counts-per-million, normalised expression, and pseudobulk) to support cross-dataset analyses and downstream integration.

---

## Summary

cellNexus aggregates over 40 million human cells, processed with harmonised pipelines so that metadata, QC flags, and expression representations are comparable across studies. This repository holds scripts and notebooks used for figures and supplementary analyses: metadata preparation from Census releases, tissue-level QC landscapes, cell-type consensus mapping, doublet-aware differential expression, sex and ethnicity prediction experiments, and cell–cell communication analyses.

---

## Repository layout

| Path | Description |
|------|-------------|
| `metadata/` | Stepwise R scripts (numbered `step1`–`step10`) for CELLxGENE metadata extraction by cellxgenedp, CELLxGENE-Census API data download, large-sample handling, HPCell-based quality control processing, unified metadata, single-cell and pseudobulk local and cloud cache preparation, and API tests against the **cellNexus** R package. |
| `quality_control/` | Quarto report for tissue landscape QC figures. |
| `cell_typing/` | Quarto notebooks for cell-typing methodology and results (consensus / immune atlas mapping). |
| `doublet_DE/` | Doublet-aware differential expression and supporting R helpers/themes. |
| `sex_prediction/` | Sex prediction analysis. |
| `ethnicity_prediction/` | Ethnicity prediction. |
| `communication/` | Cell–cell communication analyses. |
| `annotate_harmony_clusters.R` | Harmony cluster annotation utilities. |
| `tissue_color_utils.R`, `tissue_*.csv` | Shared colour maps and tissue label → macrocategory mappings. |

Paths inside many scripts point to institutional HPC or personal scratch directories. **You must edit those paths** (and optionally Census version pins such as `2024-07-01`) to match your environment before re-running.

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
