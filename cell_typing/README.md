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

- `cell_typing_method.qmd`: diagnostic and characteristic of the immune graph and the references
- `cell_typing_method.html`: associated html report for `cell_typing_method.qmd`
- `cell_typing_results.qmd`: summary statistics of consensus cell typing
- `cell_typing_results.html`: associated html report for `cell_typing_results.qmd`
- `B_memory_case_study.qmd`: showcasing how cellNexus consensus annotation improves pseudobulk analysis
- `B_memory_case_study.html`: associated html report for `B_memory_case_study.qmd`
- `HCAO_tree.qmd`: build and annotate alternative immune graph from HCAO and run cell typing with it
- `HCAO_tree.html`: associated html report for `HCAO_tree.qmd`

## Requirements

- R >= 4.3
- Quarto

## Install R packages

The following packages can be installed from CRAN:

```r

install.packages(c(
  "tidyverse", "duckdb", "arrow", "patchwork","igraph", "DT","tidyHeatmap",
  "tidygraph","ggraph","scico","ggalluvial","ggrepel","tidybulk","igraph","rdflib"
))
```

The following packages can be installed from Bioconductor:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "BiocParallel", "ComplexHeatmap",
  "SummarizedExperiment", "SingleCellExperiment","scuttle",
  "edgeR", "limma", "variancePartition"
  "org.Hs.eg.db", "cellNexus","vissE"
))
```

The following package can be installed from GitHub:

```r
devtools::install_github("stemangiola/CuratedAtlasQueryR")
```

## Run order

1. Render `cell_typing_method.qmd` to create report.

2. Render `cell_typing_results.qmd` to create report and:
  - `caq_celltype_level_map.csv`
  - `cell_annotation_new.parquet`
  
3. Render `B_memory_case_study.qmd` to create report.

4. Render `HCAO_tree.qmd` to create report and:
  - `HCAO_immune_graph_vertices.csv`
  - `HCAO_immune_graph_edges.csv`
