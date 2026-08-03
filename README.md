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

### `metadata/` — Metadata harmonisation pipeline

Stepwise R scripts that built the cellNexus resource from scratch: downloading CELLxGENE Census data from AWS, extracting and harmonising metadata with **cellxgenedp**, splitting large samples, running **HPCell**-based QC at scale, unifying metadata, and preparing single-cell and pseudobulk caches. Steps 2–9 are orchestrated with **{targets}** and dispatched across SLURM workers via **{crew.cluster}**; they require a SLURM-based HPC cluster with substantial parallelisation and storage resources.

| File | Description |
|------|-------------|
| `step1_downloads_census_datasets_from_aws.R` | Identifies the stable CELLxGENE Census snapshot (2024-07-01), generates AWS S3 download commands, and runs them in parallel to fetch all dataset `.h5ad` files. |
| `step2_cellxgene_to_metadata.R` | Extracts and harmonises cell- and sample-level metadata from Census datasets using the **{cellxgenedp}** API; writes `sample_metadata.parquet`, `cell_ids_for_metadata.parquet`, `age_days.parquet`, and `tissue_grouped.parquet`. Dispatched across SLURM workers. |
| `step3_split-large_samples.R` | Identifies samples exceeding size thresholds and splits their cell IDs so that downstream per-sample h5ad files remain manageable. Uses the Census remote API. |
| `step4_split_census_anndata_base_on_sample_id.R` | Splits each per-dataset AnnData file into per-sample `.h5ad` files using the observation join IDs from step 2. Parallelised across SLURM workers via **{crew.cluster}**. |
| `step5_identify_census_sample_counts_distribution.R` | Computes diagnostic metrics on per-sample raw count matrices (value range, integrality, min/median ratio, positive mode) to classify the appropriate count transformation for each sample. Dispatched across SLURM workers via **{crew.cluster}**. Requires a SLURM-based HPC cluster. |
| `step6_execute_hpcell_on_census_and_defining_data_tranformation.R` | Runs the **HPCell** QC pipeline on all per-sample `.h5ad` files and applies a decision-tree count-transformation classifier based on the metrics from step 5. Dispatched across a tiered hierarchy of SLURM workers. Requires a SLURM-based HPC cluster. |
| `step7_prepare_local_cache_splitting_du_dataset_and_cell_type.R` | Builds the cellNexus local single-cell cache from HPCell QC output: assigns stable hash-based file IDs, merges all metadata via DuckDB SQL, and writes four expression layers (raw counts, CPM, per-gene rank, SCTransform) to the cache directory. Dispatched across SLURM workers via **{crew.cluster}**. Requires a SLURM-based HPC cluster. |
| `step7_supp_dataset_cell_map.R` | Supplementary to step 7: builds a stable integer-indexed mapping between cell IDs and their cellNexus single-cell file IDs; writes a compressed parquet dictionary used in downstream cache preparation. Dispatched across SLURM workers via **{crew.cluster}**. |
| `step8_prepare_pseudobulk_local_cache.R` | Aggregates QC-processed single-cell data into pseudobulk `SummarizedExperiment` objects stratified by dataset, sample, and consensus cell type; writes chunked `.h5ad` files to the cellNexus local pseudobulk cache. Dispatched across a tiered hierarchy of SLURM workers. |
| `step9_unify_and_update_sce_metadata.R` | Cleans and finalises the cellNexus metadata: renames columns, executes DuckDB SQL harmonisation passes, and writes the unified metadata parquet and HDF5 files ready for shipping to the cellNexus package. |
| `step10_query_cellNexus_data_local_and_cloud.R` | Validation script: exercises the public cellNexus R API against both the cloud cache and a local cache to confirm all expression layers, metadata, and QC flags are accessible. Run after steps 7–9 outputs have been uploaded. |
| `census_accepted_assays_2024-07-01.csv` | Reference list of assay types accepted from the July 2024 Census snapshot; used as a filter in the metadata-extraction steps. |

---

### `quality_control/` — Tissue-level QC figures

Scripts and a Quarto report producing the tissue-level QC landscape figures in the paper.

| File | Description |
|------|-------------|
| `tissue_landscape.Qmd` | Generates Figure 3A. Main QC Quarto report: assembles all tissue-level QC landscape figures, including technology-by-tissue bubble plots, cell/sample count scatterplots, library-size density distributions, age distributions, and alluvial composition plots. |
| `Fig2.R` | Assembles Figure 2 panel from individual QC plot components. |
| `Fig3.R` | Assembles Figure 3 panel. |
| `Supplementary Figure 1.R` | Generates Supplementary Figure 1. |
| `2A_alluvial_plots.R` | Alluvial plots showing cell-type composition shifts across datasets. |
| `2B_age_distribution_by_sex_ethnicity_per_tissue.R` | Age distribution broken down by sex and ethnicity for each tissue. |
| `2C_scatterplot_cell_sample_count_of_dataset.R` | Scatter plot of per-dataset cell and sample counts. |
| `2D_lib_size_density.R` | Library-size density plots comparing raw and normalised data. |
| `CAQ_age_analysis_functions.R` | Helper functions for age-stratified CuratedAtlasQuery analyses, shared with `doublet_DE/`. |
| `clean_metadata.R` | Metadata cleaning utilities; calls `get_metadata()` from the cellNexus package. |
| `plot_custom_theme.R` | Custom ggplot2 theme applied consistently across all QC figures. |

---

### `cell_typing/` — Immune-graph consensus cell typing

Quarto notebooks documenting the cell-typing methodology and results. Input files are available from [Zenodo (record 19909778)](https://zenodo.org/records/19909778).

| File | Description |
|------|-------------|
| `cell_typing_method.qmd` | Immune-graph construction and consensus cell-typing methodology: graph structure, reference datasets, and diagnostic characterisation. |
| `cell_typing_method.html` | Pre-rendered HTML report for `cell_typing_method.qmd`. |
| `cell_typing_results.qmd` | Summary statistics of the consensus cell-typing results; produces `caq_celltype_level_map.csv` and `cell_annotation_new.parquet`. |
| `cell_typing_results.html` | Pre-rendered HTML report for `cell_typing_results.qmd`. |
| `B_memory_case_study.qmd` | Memory B-cell case study showing how cellNexus consensus annotation improves pseudobulk DE analysis. |
| `B_memory_case_study.html` | Pre-rendered HTML report for `B_memory_case_study.qmd`. |
| `HCAO_tree.qmd` | Builds and annotates an alternative immune graph from the Human Cell Atlas Ontology (HCAO) and runs cell typing with it. |
| `HCAO_tree.html` | Pre-rendered HTML report for `HCAO_tree.qmd`. |

---

### `doublet_DE/` — Doublet-aware differential expression

| File | Description |
|------|-------------|
| `cellNexus_doublets_DE.qmd` | Quarto notebook demonstrating the impact of doublet contamination on differential expression results and the benefit of doublet-aware DE. |
| `CAQ_age_analysis_functions.R` | Helper functions for age-stratified analyses (shared with `quality_control/`). |
| `plot_custom_theme.R` | Custom ggplot2 theme (shared with `quality_control/`). |

---

### `sex_prediction/` — Donor sex prediction

| File | Description |
|------|-------------|
| `sex_prediction.Rmd` | Trains and evaluates a classifier for predicting donor sex from single-cell expression data. Processed pseudobulk input (`pseudobulk_se.h5ad`) available from [Zenodo (record 17060800)](https://zenodo.org/records/17060800). |

---

### `ethnicity_prediction/` — Ethnicity label harmonisation and benchmarking

Fixed processed inputs for the recommended entry point are available from [Zenodo (record 19917564)](https://zenodo.org/records/19917564).

| File | Description |
|------|-------------|
| `Ethnicity_CellNexus.qmd` | Main R/Quarto analysis: loads the processed pseudobulk object, harmonises ethnicity labels, defines high-confidence reference and query sets, exports `.h5ad` files for the Python benchmarking step, and generates Figure 4 A,B and C, Supplementary Figures S3, S8, S9, S10 and S11. |
| `Classifier_benchmark.ipynb` | Python notebook that benchmarks multiple classifiers under dataset-held-out and tissue-held-out cross-validation and generates final query-label predictions. |
| `etnicity_pipeline_targets.R` | Optional upstream preprocessing pipeline: regenerates the processed pseudobulk object from source using a **{targets}** pipeline. Requires a SLURM-based HPC cluster (job scheduling via `crew.cluster`). |
| `utils.R` | Helper functions used by `Ethnicity_CellNexus.qmd`. |
| `Ethnicity_CellNexus.pdf` | Ethnicity analysis report. |

---

### `communication/` — Cell–cell communication

| File | Description |
|------|-------------|
| `communication_analyses.Rmd` | R Markdown analysis of ligand–receptor interactions across tissues using the cellNexus data. |
| `_targets.R` | **{targets}** pipeline definition that orchestrates the communication analysis steps. |
| `logcounts_riboRegressed_bySample.h5` | Pre-processed log-counts expression matrix (ribosomal-gene regressed, per-sample) used as input to the communication analysis. |

---

### Root-level shared utilities

| File | Description |
|------|-------------|
| `annotate_harmony_clusters.R` | Shared utility for annotating Harmony-derived clusters; used across multiple analyses. |
| `tissue_color_utils.R` | Shared colour palettes mapping tissues to consistent hex colours across all figures. |
| `tissue_groups_to_labels.csv` | Mapping from tissue group identifiers to human-readable tissue labels. |
| `tissue_labels_to_macrocategories.csv` | Mapping from fine-grained tissue labels to broad macrocategory groupings used in figures. |

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

Source observational data are drawn from the [**CELLxGENE Census**](https://chanzuckerberg.github.io/cellxgene-census/) (see manuscript for the Census snapshot / policy). Processed layers and metadata consumed by **cellNexus** are accessed via the **cellNexus** R package and paths documented in the paper’s Data availability statement.

---

## Licence

Code in this repository is licensed under the **GNU General Public License v3.0** — see [`LICENSE`](LICENSE).

---

## Contact

Documentation for the **cellNexus** R and Python API is on the pkgdown site at [**cellnexus.org**](https://cellnexus.org/). For bug reports and package development, use issues or discussions on the main [**cellNexus** GitHub repository](https://github.com/MangiolaLaboratory/cellNexus) where appropriate.
