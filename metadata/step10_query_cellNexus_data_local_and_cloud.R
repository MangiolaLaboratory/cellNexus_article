# This script validates the final cellNexus resource by exercising the public
# R API against both the cloud cache and a local cache. It checks metadata
# retrieval, quality-cell filtering, all four expression layers (counts, cpm,
# rank, sct), pseudobulk access, and QC-flag distributions. It should be run
# after the outputs of steps 7–9 have been uploaded to the Nectar cloud store.
library(dplyr)
library(cellNexus)
library(stringr)
library(zellkonverter)

cache <- tempdir()

x <- get_metadata()
x <- x |>
  keep_quality_cells()
x <- x |> filter(
  cell_type_unified_ensemble |>
    str_like("%bone%")
)

# Test all sce assay
sce <- x |>
  get_single_cell_experiment(assays = c("counts", "cpm", "rank", "sct"))
sce

pseudobulk <- x |>
  get_pseudobulk()
pseudobulk

# Check the number of cells per dataset
x |> dplyr::count(dataset_id)

# Check the number of cells per sample_id
x |>
  dplyr::count(sample_id) |>
  dplyr::count(n > 10)

# Check whether numeric cell_id strategy is implemented
x |>
  select(cell_id) |>
  arrange(cell_id)

# Check QC metrics
get_metadata() |>
  filter(
    cell_type_unified_ensemble |>
      str_like("%bone%")
  ) |>
  dplyr::count(empty_droplet, alive, scDblFinder.class)
