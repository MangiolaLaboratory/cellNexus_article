# This script runs the HPCell quality-control pipeline on all per-sample h5ad
# files produced in step 4, and determines the appropriate count transformation
# for each sample (raw counts, CPM, or SCTransform) based on a decision tree
# applied to the count-distribution metrics computed in step 5 supplementary.
# The pipeline is orchestrated with {targets} and dispatched across a tiered
# hierarchy of SLURM workers. It requires HPC infrastructure and is not intended to be run locally.

library(dplyr)
library(tibble)
library(glue)
library(purrr)
library(stringr)
library(HPCell)
library(arrow)
library(targets)
library(crew)
library(crew.cluster)
library(duckdb)
directory = "split_h5ad_based_on_sample_id/2024-07-01/" # MODIFY HERE: directory containing per-sample h5ad files
downloaded_samples_tbl <- read_parquet("metadata_cellxgenedp_Apr_2024/census_samples_to_download_groups_MODIFIED.parquet") # MODIFY HERE: input samples metadata parquet
downloaded_samples_tbl <- downloaded_samples_tbl |>
  dplyr::rename(cell_number = list_length) |>
  mutate(cell_number = cell_number |> as.integer(),
         file_name = glue("{directory}{sample_2}.h5ad") |> as.character())

sample_tbl = downloaded_samples_tbl |>
  filter(!dataset_id %in% c("99950e99-2758-41d2-b2c9-643edcdf6d82", "9fcb0b73-c734-40a5-be9c-ace7eea401c9" )) |>
  left_join(
    cellxgenedp::datasets() |>
      dplyr::select(dataset_id, x_approximate_distribution) |>
      distinct(), by = "dataset_id", copy = TRUE) |>
  mutate(cell_number = cell_number |> as.integer(),
         file_name = glue("{directory}{sample_2}.h5ad") |> as.character()) |> 
  
  left_join(
    cellNexus::get_metadata() |> # MODIFY HERE: cellNexus local cache directory
      cellNexus::join_census_table() |>
      distinct(sample_id, assay) ,
    by = c("sample_2" = "sample_id"),
    copy = T
  ) |>
  # Propositional set up expressed genes threshold for panel technologies 500/20K = x/462
  mutate(feature_thresh = ifelse(assay == "BD Rhapsody Targeted mRNA", 11, 200))
  
sample_summary_df = tar_read(sample_summary_df, store = "calculate_census_raw_counts_target_store/_targets/") |> # MODIFY HERE: targets store from step5_identify_census_sample_counts_distribution.qmd 
  bind_rows() |>
  mutate(max_gt_20 = ifelse(max_val > 20, TRUE, FALSE))

# Impute distribution decision tree
impute_x_approximate_distribution <- function(df,
                                              counts_gap_threshold,
                                              pos_mode_threshold) {
  df |>
    dplyr::mutate(
      inferred_distribution = dplyr::case_when(
        
        # 0) When counts gap between 0 and next min value >= threshold
        !has_negative & !max_gt_20 & !all_integer & !has_floating &
          (counts_gap_min_mean >= counts_gap_threshold) & (positive_mode > pos_mode_threshold) ~ "double_log1p",
        
        # 1) Small counts gap
        !has_negative & !max_gt_20 & !all_integer & !has_floating &
          !(
            (counts_gap_min_mean >= counts_gap_threshold) &
              (positive_mode > pos_mode_threshold)
            
          ) ~ "log1p",
        
        # 2) No negatives, has large values
        !has_negative & max_gt_20 & !all_integer & !has_floating ~ "raw_limit_max_to_10",
        
        # 3) Large values, integer counts
        !has_negative & max_gt_20 & all_integer & !has_floating ~ "raw_limit_max_to_10",
        
        # 4) Has negatives, compressed range
        has_negative & !max_gt_20 & !all_integer & !has_floating ~ "raw_limit_max_to_10",
        
        # 5) Has negatives and large values
        has_negative & max_gt_20 & !all_integer & !has_floating ~ "raw_limit_max_to_10",
        
        # fallback
        TRUE ~ NA_character_
      )
    )
}

sample_summary_df = sample_summary_df |> impute_x_approximate_distribution(0.25, 1) |> 
  mutate(count_upper_bound = case_when(
    # 0) When counts gap between 0 and next min value >= 0.25, double log. Max value before exp is 10.
    inferred_distribution == "double_log1p" ~ 10,
    
    # 1) make 10 as max before exp
    inferred_distribution == "log1p" ~ 10,
    
    # 2,3,5) transform algo picks up negative value. should always scale max to 10
    # 4) Has negatives, no large values, no integer, no floating. Counts peak at 10
    inferred_distribution == "raw_limit_max_to_10" ~ 10
     
  )) |>
  # Inverse distribution
  mutate(method_to_apply = case_when(inferred_distribution == "double_log1p" ~ "safe_expm1",
                                     inferred_distribution == "log1p" ~ "expm1",
                                     inferred_distribution == "raw_limit_max_to_10" ~ "identity_with_max_limit"))

sample_tbl = sample_tbl |> left_join(sample_summary_df |>
                                       select(sample_2,
                                              method_to_apply,
                                              dataset_id,
                                              count_upper_bound),
                                     by = c("sample_2", "dataset_id")) |> 
  
  # This should be addressed in step 5. NA method_to_apply has max counts =0, fall in to case 4)
  mutate(method_to_apply = if_else(is.na(method_to_apply), "identity_with_max_limit", method_to_apply),
         count_upper_bound = if_else(is.na(count_upper_bound), 10, count_upper_bound),
         )

sample_tbl = sample_tbl |>
  
  select(file_name, cell_number, dataset_id, sample_2, method_to_apply, assay, count_upper_bound, feature_thresh)

sample_tbl |> write_parquet("updated_transform_sample_tbl_2024_Jul.parquet")
sample_tbl <- read_parquet("updated_transform_sample_tbl_2024_Jul.parquet")  # MODIFY HERE: output path for sliced_sample_tbl RDS

# Enable sample_names.rds to store sample names for the input
sample_names <-
  sample_tbl |> 
  pull(file_name) |> 
  set_names(sample_tbl |> pull(sample_2))
functions = sample_tbl |> pull(method_to_apply)
feature_thresh = sample_tbl |> pull(feature_thresh)
count_upper_bound = sample_tbl |> pull(count_upper_bound)


my_store = "2024-07-01/process_updated_samples_transform_hpcell_target_store_v1" # MODIFY HERE: HPCell targets store (used throughout this script)

new_elastic <- function(name, mem_gb, time_min, workers, crashes_max, cpus_per_task = 1, backup = NULL) {
  crew_controller_slurm(
    name = name,
    workers = workers,
    crashes_max = crashes_max,
    seconds_idle = 30,
    options_cluster = crew_options_slurm(
      memory_gigabytes_required = mem_gb,
      cpus_per_task = cpus_per_task,
      time_minutes = time_min
    ),
    backup = backup
  )
}
elastic_160 <- new_elastic("elastic_160", 160, 60 * 24, workers = 8,  crashes_max = 2)
elastic_120  <- new_elastic("elastic_120",  120,  60 * 8,  workers = 16, crashes_max = 1, cpus_per_task = 1, backup = elastic_160)
elastic_80  <- new_elastic("elastic_80",   80,  60 * 8,  workers = 24, crashes_max = 1, cpus_per_task = 1, backup = elastic_120)
elastic_40  <- new_elastic("elastic_40",   40,  60 * 4,  workers = 32, crashes_max = 1, cpus_per_task = 1, backup = elastic_80)
elastic_20  <- new_elastic("elastic_20",   20,  60 * 4,  workers = 48, crashes_max = 1, cpus_per_task = 1, backup = elastic_40)
elastic_10   <- new_elastic("elastic_10",   10, 60 * 4,  workers = 150, crashes_max = 2, cpus_per_task = 1, backup = elastic_20)

elastic_5_minimal   <- new_elastic("elastic_5_minimal",     5, 60 * 4,  workers = 300, crashes_max = 2, cpus_per_task = 1, backup = elastic_10)

# Group for targets (small → large)
controllers <- crew_controller_group(
  elastic_10, elastic_20, elastic_40, elastic_80, elastic_120, elastic_160, elastic_5_minimal
)

job::job({
  
  library(HPCell)
  
  sample_names |>
    initialise_hpc(
      store = my_store,
      gene_nomenclature = "ensembl",
      data_container_type = "anndata",
      computing_resources = list(
        elastic_5_minimal, elastic_10, elastic_20, elastic_40, elastic_80, elastic_120, elastic_160
      ),
      default_controller = "elastic_5_minimal", 
      verbosity = "summary",
      update = "never", 
      #update = "thorough", 
      error = "continue",
      garbage_collection = 100, 
      workspace_on_error = TRUE
      
    ) |> 
    transform_assay(fx = functions, target_output = "sce_transformed", scale_max = count_upper_bound) |>

    # # Remove empty outliers based on RNA count threshold per cell
    remove_empty_threshold(target_input = "sce_transformed", RNA_feature_threshold = feature_thresh) |>

    # Annotation
    annotate_cell_type(target_input = "sce_transformed", azimuth_reference = "pbmcref") |>

    # Cell type harmonisation
    celltype_consensus_constructor(target_input = "sce_transformed",
                                   target_output = "cell_type_concensus_tbl") |>

    # Alive identification
    remove_dead_scuttle(target_input = "sce_transformed", target_annotation = "cell_type_concensus_tbl",
                        group_by = "cell_type_unified_ensemble") |>

    # Doublets identification
    remove_doublets_scDblFinder(target_input = "sce_transformed") |>
    
    # SCT
    normalise_abundance_seurat_SCT(target_input = "sce_transformed", factors_to_regress = c(
      "subsets_Mito_percent",
      "subsets_Ribo_percent")) |>
    
    # Pseudobulk
    calculate_pseudobulk(target_input = "sce_transformed",
                         group_by = "cell_type_unified_ensemble") |>

    # # metacell
    # cluster_metacell(target_input = "sce_transformed",  group_by = "cell_type_unified_ensemble") |>

    # # Cell Chat
    # ligand_receptor_cellchat(target_input = "sce_transformed",
    #                          group_by = "cell_type_unified_ensemble") |>
    
    print()
  
  
})

# Sample metadata
library(arrow)
library(dplyr)
library(duckdb)
library(targets)

# Write annotation light
# MODIFY HERE: base cell_metadata parquet path inside the SQL string below
con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
cell_metadata <- 
  tbl(
    con,
    sql("SELECT * FROM read_parquet('metadata_cellxgenedp_Apr_2024/cell_metadata.parquet')")
  ) |>
  mutate(cell_ = paste0(cell_, "___", dataset_id)) |> 
  select(cell_, observation_joinid, contains("cell_type"), dataset_id,  self_reported_ethnicity, tissue, donor_id,  sample_id, is_primary_data, assay)

empty_droplet = 
  tar_read(empty_tbl, store = "2024-07-01/process_updated_samples_transform_hpcell_target_store_v1") |> # MODIFY HERE: HPCell targets store (must match my_store above)
  bind_rows() |>
  dplyr::rename(cell_ = .cell)

alive_cells = 
  tar_read(alive_tbl, store = "2024-07-01/process_updated_samples_transform_hpcell_target_store_v1") |> # MODIFY HERE: HPCell targets store (must match my_store above)
  bind_rows() |>
  select(-any_of(c("cell_type_unified_ensemble",  "observation_originalid"))) |>
  dplyr::rename(cell_ = .cell)

doublet_cells =
  tar_read(doublet_tbl, store ="2024-07-01/process_updated_samples_transform_hpcell_target_store_v1") |> # MODIFY HERE: HPCell targets store (must match my_store above)
  bind_rows() |>
  dplyr::rename(cell_ = .cell)

metacell = 
  tar_read(metacell_tbl, store = "2024-07-01/process_updated_samples_transform_hpcell_target_store_v1") |> # MODIFY HERE: HPCell targets store (must match my_store above)
  bind_rows() |> 
  dplyr::rename(cell_ = cell) |> 
  dplyr::rename_with(
    ~ stringr::str_replace(.x, "^gamma", "metacell_"),
    starts_with("gamma")
  )

# Save cell type concensus tbl from HPCell output to disk
cell_type_concensus_tbl = tar_read(cell_type_concensus_tbl, store = "2024-07-01/process_updated_samples_transform_hpcell_target_store_v1") |> # MODIFY HERE: HPCell targets store (must match my_store above)
  bind_rows() |> 
  dplyr::rename(cell_ = .cell)

cell_type_concensus_tbl = cell_type_concensus_tbl |> mutate(cell_type_unified_ensemble = 
                                    ifelse(is.na(cell_type_unified_ensemble),
                                           "Unknown",
                                           cell_type_unified_ensemble)) 

# This command needs a big memory machine
cell_metadata_joined = cell_metadata |> 
  left_join(empty_droplet, copy=TRUE) |>  
  left_join(cell_type_concensus_tbl, copy=TRUE) |>
  left_join(alive_cells, copy=TRUE) |> 
  left_join(doublet_cells, copy=TRUE)
# |>
#   left_join(metacell, copy=TRUE)

cell_metadata_joined |> filter(is.na(blueprint_first_labels_fine))

cell_metadata_joined2 = cell_metadata_joined |>
  mutate(
    cell_type_unified_ensemble    = coalesce(cell_type_unified_ensemble,    "Unknown"),
    data_driven_ensemble          = coalesce(data_driven_ensemble,          "Unknown"),
    blueprint_first_labels_fine   = coalesce(blueprint_first_labels_fine,   "Other"),
    monaco_first_labels_fine      = coalesce(monaco_first_labels_fine,      "Other"),
    azimuth_predicted_celltype_l2 = coalesce(azimuth_predicted_celltype_l2, "Other"),
    azimuth                       = coalesce(azimuth,                       "Other"),
    blueprint                     = coalesce(blueprint,                     "Other"),
    monaco                        = coalesce(monaco,                        "Other")
  )

output_path <- "cell_annotation_2024_Jul.parquet" # MODIFY HERE: output cell annotation parquet (used as input to step7)
final_sql <- dbplyr::sql_render(cell_metadata_joined2)
DBI::dbExecute(con, sprintf(
  "COPY (%s) TO '%s' (FORMAT PARQUET, COMPRESSION 'zstd')",
  final_sql, output_path
))

# Cellchat output
ligand_receptor_tbl = tar_read(ligand_receptor_tbl, store = "2024-07-01/process_updated_samples_transform_hpcell_target_store_v1") |> bind_rows() # MODIFY HERE: HPCell targets store (must match my_store above)
# save
con <- dbConnect(duckdb::duckdb(), dbdir = "cellNexus_lr_signaling_pathway_strength.duckdb") # MODIFY HERE: output DuckDB file for ligand-receptor results
duckdb::dbWriteTable(con, "lr_pathway_table", ligand_receptor_tbl, overwrite = TRUE)
dbDisconnect(con)


# Helper function to save parquet read by duckdb to parquet on disk
# write_parquet_to_parquet = function(data_tbl, output_parquet, compression = "gzip") {
#   
#   # Establish connection to DuckDB in-memory database
#   con_write <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
#   
#   # Register `data_tbl` within the DuckDB connection (this doesn't load it into memory)
#   duckdb::duckdb_register(con_write, "data_tbl_view", data_tbl)
#   
#   # Use DuckDB's COPY command to write `data_tbl` directly to Parquet with compression
#   copy_query <- paste0("
#   COPY data_tbl_view TO '", output_parquet, "' (FORMAT PARQUET, COMPRESSION '", compression, "');
#   ")
#   
#   # Execute the COPY command
#   dbExecute(con_write, copy_query)
#   
#   # Unregister the temporary view
#   duckdb::duckdb_unregister(con_write, "data_tbl_view")
#   
#   # Disconnect from the database
#   dbDisconnect(con_write, shutdown = TRUE)
# }
