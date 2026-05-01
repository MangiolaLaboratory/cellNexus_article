# This scripts test cellNexus API with new data generated in STEP_7_unify_and_update_sce_metadata.R
# Before running this script, it assumes that files generated in step6-8 are uploaded to Nectar cloud.
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
