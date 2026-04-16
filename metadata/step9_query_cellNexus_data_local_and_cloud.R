# This scripts test cellNexus API with new data generated in ~/git_control/cellNexus/dev/STEP_7_unify_cell_metadata.R
# Before running this script, it assumes that files generated in step6-8 are uploaded to Nectar cloud.
library(dplyr)
library(cellNexus)
library(stringr)
library(zellkonverter)

cache = "~/scratch/cache_temp"

x = get_metadata()
x = x |>
  keep_quality_cells()
x = x |> dplyr::filter(
    self_reported_ethnicity == "African" &
      assay |> stringr::str_like("%10x%") & 
      tissue == "lung parenchyma" &
      cell_type |> stringr::str_like("%CD4%")
  )


# One anndata
anndata = readH5AD("~/scratch/cellNexus/cellxgene/01-07-2024/counts/9e62207287ebeaa020d3e92d17b01f8e___1.h5ad", reader="R",use_hdf5 = T)
anndata

# Test SCE
sce = x |> 
  get_single_cell_experiment()

# TEST CPM
cpm = x |> 
  get_single_cell_experiment(assays = "cpm")

# TEST RANK
rank = x |> 
  get_single_cell_experiment(assays = "rank")

# TEST SCT
sct = x |> get_single_cell_experiment(assays = "sct")

pseudobulk = x |> 
  get_pseudobulk() 

# Check the number of cells per dataset
x |> dplyr::count(dataset_id)

# Check the number of cells per sample_id
x |> dplyr::count(sample_id) |> dplyr::count(n>10)

# Check whether cell_id strategt is implemented
x |> select(cell_id) |> arrange(cell_id)

# Check QC metrics
x |> dplyr::count(empty_droplet, alive, scDblFinder.class)

# Check empty droplet ratio
x |> dplyr::count(empty_droplet) |>
  collect() |> mutate(n_cells = sum(n), pct = n / n_cells * 100)

# Check alive ratio
x |> filter(!empty_droplet) |> dplyr::count(alive) |> 
  collect() |> mutate(n_cells = sum(n), pct = n / n_cells * 100)

# Check doublet ratio
x |> filter(!empty_droplet, alive) |> dplyr::count(scDblFinder.class) |> 
  collect() |> mutate(n_cells = sum(n), pct = n / n_cells * 100)

