library(targets)
library(tidyverse)
# Get summary table in parallel
summary_store = "sct_failed_samples_raw_counts_summary_target_store"
tar_script({
  library(dplyr)
  library(SummarizedExperiment)
  library(zellkonverter)
  library(crew)
  library(crew.cluster)
  library(duckdb)
  
  
  computing_resources = crew_controller_slurm(
    name = "elastic",
    workers = 300,
    tasks_max = 20,
    seconds_idle = 30,
    crashes_error = 10,
    options_cluster = crew_options_slurm(
      #memory_gigabytes_required = c(30, 45, 45, 90, 200, 1000, 1500), 
      memory_gigabytes_required = c(90, 130, 170, 200, 1000, 1500), 
      cpus_per_task = c(2, 2, 2, 2, 2, 2), 
      time_minutes = c(30*4, 30*4, 30*4, 60*4, 60*24, 60*24),
      verbose = T
    ))
  
  tar_option_set(
    memory = "transient",
    garbage_collection = TRUE,
    storage = "worker",
    retrieval = "worker",
    format = "qs",
    cue = tar_cue(mode = "never"),
    error = "continue",
    controller = computing_resources
  )
  
  get_sample_summary_stats <- function(files) {
    sce = readH5AD(files, reader = "R", use_hdf5 = T)
    
    if (ncol(sce) == 0) return(NULL)
    
    assay_name = sce@assays |> names() |> magrittr::extract(1)
    counts_mat = sce |> assay(assay_name) |> as.matrix()
    counts_vec = as.numeric(counts_mat)
    sample_id = basename(files)
    
    # Perform checks
    min_val = min(counts_vec, na.rm = TRUE)
    max_val = max(counts_vec, na.rm = TRUE)
    has_negative = min_val < 0
    max_gt_10 = max_val > 10
    
    # Check if all integers
    tol = 1e-5
    all_integer = all(counts_vec == floor(counts_vec), na.rm = TRUE)
    
    #should've excluded integer, but I will fix the output dataframe
    #has_floating = all(abs(counts_vec - round(counts_vec)) < tol, na.rm = TRUE)
    has_floating = !all_integer && all(abs(counts_vec - round(counts_vec)) < tol, na.rm = TRUE)
    
    tbl = tibble::tibble(
      sample_id = sample_id,
      min_val = min_val,
      max_val = max_val,
      has_negative = has_negative,
      max_gt_10 = max_gt_10,
      all_integer = all_integer,
      has_floating = has_floating,
      n_cells = ncol(counts_mat),
      n_genes = nrow(counts_mat)
    )
    
    
    tbl
  }
  
  list(
    tar_target(
      files,
      arrow::read_parquet("metadata/failed_sct_samples_to_rerun_Jul_2024.parquet") |> pull(file_name),
      deployment = "main"
    ),
    
    # Get raw counts matrix with sample_id
    tar_target(
      sample_summary_df,
      get_sample_summary_stats(files),
      pattern = map(files),
      iteration = "list"
    )
  )
}, ask = FALSE,  script = glue("{summary_store}/_targets.R"))


job::job({
  
  tar_make(
    # callr_function = NULL,
    reporter = "summary",
    script = glue("{summary_store}/_targets.R"),
    store = glue("{summary_store}/_targets")
  )
  
})