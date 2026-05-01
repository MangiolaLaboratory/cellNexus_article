# Cell typing workflow

This directory contains CellNexus workflows for cell typing: including the work
flow for cell typing with the heuristic immune graph, the workflow for case study
of memory B cells, and the workflow for cell typing with HCAO immune graph.

## Inputs and Downloads

The following input files can be downloaded from [https://zenodo.org/records/19909778](https://zenodo.org/records/19909778)
and unzipped into this folder as the input file for workflows.

- `cell_annotation.parquet`
- `cell_type_harmonisation_non_immune.csv`
- `data_driven_consensus.parquet`
- `graphs/immune_map_azimuth.csv`
- `graphs/immune_map_blueprint.csv`
- `graphs/immune_map_cellxgene.csv`
- `graphs/immune_map_monaco.csv`
- `graphs/immune_tree.csv`

## Files

- `cell_typing_method.qmd`: score ethnicity labels, selects marker genes, defines reference/query sets, exports `.h5ad` files, runs Seurat `LabelTransfer`, and generates figures.
- `cell_typing_method.html`: score ethnicity labels, selects marker genes, defines reference/query sets, exports `.h5ad` files, runs Seurat `LabelTransfer`, and generates figures.
- `cell_typing_results.qmd`: benchmarks classifiers under dataset-held-out and tissue-held-out cross-validation, evaluates ROC and per-ethnicity performance, and predicts query labels with the final model.
- `cell_typing_results.html`: benchmarks classifiers under dataset-held-out and tissue-held-out cross-validation, evaluates ROC and per-ethnicity performance, and predicts query labels with the final model.
- `B_memory_case_study.qmd`:
- `B_memory_case_study.html`:
- `HCAO_tree.qmd`:
- `HCAO_tree.html`:

## Requirements

- R >= 4.3
- Quarto

## Install R packages

```r
install.packages(c(
  "tidyverse", "tidybulk", "Seurat", "caret", "ggpubr", "yardstick",
  "nnet", "ggalluvial", "ggsci", "patchwork", "igraph", "scico",
  "cowplot", "magick", "ggplotify", "ggrepel", "pheatmap", "scales",
  "uwot", "formatR"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "SummarizedExperiment", "SingleCellExperiment", "edgeR",
  "DelayedMatrixStats", "HDF5Array", "scater", "standR", "bluster",
  "BiocNeighbors", "zellkonverter", "org.Hs.eg.db",
  "ComplexHeatmap", "circlize", "cellNexus"
))

```

## Run order

1. Render `cell_typing_method.qmd` to create report and:
  - PH

2. Render `cell_typing_results.qmd` to create report and:
  - PH

3. Render `B_memory_case_study.qmd` to create report.

4. Render `HCAO_tree.qmd` to create report and:
  - HCAO_immune_graph_vertices.csv
