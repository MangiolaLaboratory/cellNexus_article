# Ethnicity prediction workflow

This directory contains the CellNexus ethnicity-label harmonisation, benchmarking, and prediction workflow. It combines an R/Quarto analysis for marker selection and data export with a Python notebook for classifier benchmarking and final prediction.

## Inputs

A `data_10Mar`.

## Files

- `Ethnicity_CellNexus.qmd`: score ethnicity labels, selects marker genes, defines reference/query sets, exports `.h5ad` files, runs Seurat `LabelTransfer`, and generates figures.
- `classifier.ipynb`: benchmarks classifiers under dataset-held-out and tissue-held-out cross-validation, evaluates ROC and per-ethnicity performance, and predicts query labels with the final model.

## Requirements

- R >= 4.3
- Python >= 3.10
- Quarto
- Jupyter Notebook or JupyterLab

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

## Install Python packages

```
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip

pip install \
  jupyter notebook ipykernel \
  scanpy anndata numpy pandas scipy scikit-learn \
  matplotlib seaborn torch torchviz
```

## Run order

1. Render Ethnicity_CellNexus.qmd to create:
  - data/sce_relabel.h5ad
  - data/sce_query.h5ad

2. Run classifier.ipynb to create:
  - overall_performance.csv
  - pereth_performance.csv
  - roc_points.csv
  - auc_by_class.csv
  - histgbdt_Log1pTMM_S_PCA_query_predictions.csv

3. Re-render Ethnicity_CellNexus.qmd to update figures and summaries using the benchmark outputs.

## Current benchmark
Models compared:
  - Logistic regression
  - Linear SVM
  - Random forest
  - HistGradientBoosting
  - Multilayer perceptron
  - Seurat LabelTransfer

Cross-validation settings:
  - dataset-held-out
  - tissue-held-out
Final model:
  - HistGBDT__Log1pTMM_S_PCA

## Notes
The Quarto document and notebook are coupled: the notebook reads .h5ad files exported by the Quarto workflow, and the Quarto report reads benchmark CSVs written by the notebook.
This workflow expects project-specific input data and helper code from the full repository layout.
