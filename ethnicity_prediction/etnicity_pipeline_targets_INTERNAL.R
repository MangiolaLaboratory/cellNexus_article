library(targets)

setwd("~/Documents/GitHub/article_cellNexus/ethnicity_prediction")



tar_script({
  
  setwd("~/Documents/GitHub/article_cellNexus/ethnicity_prediction")

  library(tidyverse)
  library(targets)
  library(tarchetypes)
  library(glue)
  library(qs)
  library(crew)
  library(crew.cluster)
  
  # Set file path -----
  hdf5_path = "pseudobulk_sample_is_immune"

  tar_option_set(
    
    
    memory = "transient", 
    garbage_collection = 100, 
    storage = "worker", 
    retrieval = "worker", 
    error = "continue", 
    
    #cue = tar_cue(mode = "never"), 
    
    workspace_on_error = TRUE,
    format = "qs",
    
    debug = "estimates_chunk",
    
    controller = crew_controller_group(
      
      crew_controller_slurm(
        name = "elastic",
        workers = 500,
        tasks_max = 20,
        seconds_idle = 30,
        crashes_max = 7,
        options_cluster = crew_options_slurm(
          script_lines = '#SBATCH -A saigencir003',
          memory_gigabytes_required = c(5, 10, 20, 40, 80, 160), 
          cpus_per_task = 2, 
          time_minutes = c(60*4, 60*4, 60*4, 60*4, 60*24, 60*24),
          verbose = T
        )
      ),
      crew_controller_slurm(
        name = "elastic_big",
        workers = 150,
        tasks_max = 20,
        seconds_idle = 30,
        crashes_max = 5,
        options_cluster = crew_options_slurm(
          script_lines = '#SBATCH -A saigencir003',
          memory_gigabytes_required = c(80, 160), 
          cpus_per_task = 2, 
          time_minutes = c(60*24, 60*24),
          verbose = T
        )
      ),
      crew_controller_slurm(
        name = "elastic_big_30_cores",
        workers = 150,
        tasks_max = 20,
        seconds_idle = 30,
        crashes_max = 5,
        options_cluster = crew_options_slurm(
          script_lines = '#SBATCH -A saigencir003',
          memory_gigabytes_required = c(160), 
          cpus_per_task = 30, 
          verbose = T
        )
      )
    )
    
    
  )
  
  
  #-----------------------#
  # Functions
  #-----------------------#  
  
  #' Remove Unwanted Effects from a brmsfit Model
  #'
  #' This function calculates posterior residuals from a \code{brmsfit} model and combines them with 
  #' factor-specific fitted values (potentially excluding random effects or other parts of the model), 
  #' thereby producing adjusted outcomes that highlight the contribution of a specified factor or subset 
  #' of model terms.
  #'
  #' @param fit A \code{brmsfit} object, resulting from a model fitted by \code{\link[brms]{brm}}.
  #' @param newdata A data frame or list containing new data. Passed to \code{\link[brms]{fitted}} 
  #'   to obtain factor-specific fitted values at specified covariate levels.
  #' @param robust A logical value indicating whether to use robust (median-based) summaries rather 
  #'   than means. Defaults to \code{FALSE}.
  #' @param correct_by_offset A logical value indicating whether to divide the residuals by 
  #'   \code{exp(offset)} (from \code{fit$data$offset}). Defaults to \code{TRUE}.
  #' @param re_formula A formula specifying which random effects (if any) to include when generating 
  #'   fitted values. Defaults to \code{~0}, which removes random effects and thus isolates the 
  #'   contribution of fixed effects in the new data.
  #'
  #' @return A \code{tibble} containing posterior summaries of:
  #'   \itemize{
  #'     \item Adjusted outcomes (prefix: \code{adjusted___}): The combined values of the specified 
  #'     factor's fitted counts and the residuals.
  #'     \item Residuals (prefix: \code{residuals___}): The model's posterior residuals, possibly 
  #'     normalised by the offset.
  #'     \item Fitted values for the factor (prefix: \code{fitted___}): The model's fitted values based 
  #'     on the \code{re_formula} and provided \code{newdata}.
  #'   }
  #'
  #' @details
  #' The function proceeds as follows:
  #' \enumerate{
  #'   \item Extracts posterior residuals (via \code{\link[brms]{residuals}}).
  #'   \item (Optionally) divides these residuals by the exponential of the offset, if \code{correct_by_offset = TRUE}.
  #'   \item Obtains new fitted values from the model (via \code{\link[brms]{fitted}}), usually excluding random effects 
  #'         by specifying \code{re_formula = ~0}.
  #'   \item Adds these residuals to the factor-specific fitted values to obtain adjusted outcomes 
  #'         that highlight the contribution of the factor of interest.
  #'   \item Summarises all these draws (residuals, fitted values, adjusted outcomes) and returns them 
  #'         in a single \code{tibble}.
  #' }
  #'
  #' This method is particularly useful for examining how a factor or other subset of the model 
  #' affects the outcome when other model components (e.g., random intercepts) are removed. 
  #' It can assist in visualising or quantifying the partial contribution of certain terms.
  #'
  #' @examples
  #' \dontrun{
  #' # Suppose 'fit' is a brmsfit model object predicting a 'counts' outcome
  #' # We create a new data frame 'some_data' for which we want partial predictions
  #' adjusted_results <- remove_unwanted_effect(
  #'   fit,
  #'   newdata = some_data,
  #'   robust = TRUE,
  #'   correct_by_offset = TRUE,
  #'   re_formula = ~0
  #' )
  #' }
  #'
  #' @importFrom magrittr %>%
  #' @importFrom dplyr bind_cols
  #' @importFrom tibble as_tibble
  #' @importFrom brms posterior_summary fitted residuals
  #'
  #' @export
  remove_unwanted_effect = function(fit, newdata, robust = FALSE, correct_by_offset = T, re_formula = ~0){
    
    # Calculate residuals: observed counts minus fitted values, normalised by exp(offset)
    # This places residuals on a consistent scale, making them addable to adjusted predictions later.
    fitted_residuals =   fit |> residuals(robust = robust, summary = FALSE) 
    
    # Correct by offset
    if(correct_by_offset)
      fitted_residuals = fitted_residuals |>
        sweep(2, fit$data$offset |> exp(), FUN = "/")
    
    # Extract fitted values for the specified factor only, removing random effects by setting re_formula = ~0
    # 'resp = factor' focuses on the selected response variable (factor)
    fitted_values_ethnicity <- fitted(fit, newdata = newdata, re_formula = re_formula, summary = FALSE, offset=0)
    
    # Adjusted counts are obtained by adding the factor-specific fitted values and the normalised residuals
    adjusted_counts = fitted_values_ethnicity + fitted_residuals
    
    # Summarise residuals into a tibble, prefixed to denote their source
    fitted_residuals_tbl = 
      fitted_residuals |> 
      posterior_summary(robust = robust) |> 
      as_tibble()
    fitted_residuals_tbl |> colnames() = paste0("residuals___", fitted_residuals_tbl |> colnames())
    
    # Summarise the factor-only fitted values into a tibble, prefixed accordingly
    fitted_values_ethnicity_tbl = 
      fitted_values_ethnicity |> 
      posterior_summary(robust = robust) |> 
      as_tibble()
    fitted_values_ethnicity_tbl |> colnames() = paste0("fitted___", fitted_values_ethnicity_tbl |> colnames())
    
    # Summarise the adjusted counts (factor + residuals) into a tibble, prefixed for clarity
    adjusted_counts_tbl = 
      adjusted_counts |> 
      posterior_summary(robust = robust) |> 
      as_tibble()
    adjusted_counts_tbl |> colnames() = paste0("adjusted___", adjusted_counts_tbl |> colnames())
    
    # Combine all three resulting tables into one tibble
    adjusted_counts_tbl |> 
      bind_cols(fitted_residuals_tbl) |> 
      bind_cols(fitted_values_ethnicity_tbl)
  }
  
  remove_unwanted_effect_new = function(fit, newdata, robust = FALSE, correct_by_offset = T, re_formula = ~0){
    
    # Calculate residuals: observed counts minus fitted values, normalised by exp(offset)
    # This places residuals on a consistent scale, making them addable to adjusted predictions later.
    fitted_residuals =   fit |> predictive_error(robust = robust, summary = FALSE, offset = 0) 
    
    # Correct by offset
    if(correct_by_offset)
      fitted_residuals = fitted_residuals |>
        sweep(2, fit$data$offset |> exp(), FUN = "/")
    
    # Extract fitted values for the specified factor only, removing random effects by setting re_formula = ~0
    # 'resp = factor' focuses on the selected response variable (factor)
    fitted_values_ethnicity <- posterior_epred(fit, newdata = newdata, re_formula = re_formula,  offset=0)
    
    # Adjusted counts are obtained by adding the factor-specific fitted values and the normalised residuals
    adjusted_counts = fitted_values_ethnicity + fitted_residuals
    
    # Summarise residuals into a tibble, prefixed to denote their source
    fitted_residuals_tbl = 
      fitted_residuals |> 
      posterior_summary(robust = robust) |> 
      as_tibble()
    fitted_residuals_tbl |> colnames() = paste0("residuals___", fitted_residuals_tbl |> colnames())
    
    # Summarise the factor-only fitted values into a tibble, prefixed accordingly
    fitted_values_ethnicity_tbl = 
      fitted_values_ethnicity |> 
      posterior_summary(robust = robust) |> 
      as_tibble()
    fitted_values_ethnicity_tbl |> colnames() = paste0("fitted___", fitted_values_ethnicity_tbl |> colnames())
    
    # Summarise the adjusted counts (factor + residuals) into a tibble, prefixed for clarity
    adjusted_counts_tbl = 
      adjusted_counts |> 
      posterior_summary(robust = robust) |> 
      as_tibble()
    adjusted_counts_tbl |> colnames() = paste0("adjusted___", adjusted_counts_tbl |> colnames())
    
    # Combine all three resulting tables into one tibble
    adjusted_counts_tbl |> 
      bind_cols(fitted_residuals_tbl) |> 
      bind_cols(fitted_values_ethnicity_tbl)
  }
  
  get_adjusted_matrix = function(summary_df, column_adjusted){
    
    column_adjusted = enquo(column_adjusted)
    
    m = 
      summary_df |>
      unnest(!!column_adjusted) |> 
      # dplyr::filter(analysis == "observed_proportion") |> 
      select(.feature, adjusted___Estimate, sample_id) |> 
      pivot_wider(names_from = sample_id, values_from = adjusted___Estimate) |> 
      tidybulk:::as_matrix(rownames = ".feature") |> 
      as("sparseMatrix")  |> 
      Matrix::Matrix(sparse = T)
    
    # Cap infinite
    max_rm_infinite = 
      m |> 
      _[!m |> is.infinite()] |> 
      quantile(0.999)
    
    m |> 
      _[m > max_rm_infinite] = 
      max_rm_infinite
    
    m |> 
      _[m < 0] = 
      0
    
    return(m)
  }
  
  #-----------------------#
  # Pipeline
  #-----------------------#
  
  list(
    
    
    # tar_target(
    #   result_directory,
    #   "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/de_ethnicity_pseudobulk_sample"
    # ),
    # tar_target(
    #   glmGamPoi_overdispersions,
    #   {
    #     glmGamPoi_overdispersions  = readRDS("/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/glmGamPoi_all_samples_no_subsampling_cellNexus_1_0_6.rds")$overdispersions
    #     glmGamPoi_overdispersions[glmGamPoi_overdispersions>1e5] = max(glmGamPoi_overdispersions[glmGamPoi_overdispersions<1e5])
    #     glmGamPoi_overdispersions
    #   }, 
    #   deployment = "main"
    #   
    # ),
    
    # This target loads and processes the pseudobulk sample data. It imports a HDF5 SummarizedExperiment, 
    # applies filters to retain shared genes, immune cells, and samples marked for analysis, integrates age metadata,
    # filters for common genes and samples with an appropriate number of detected genes, computes the mean library size, 
    # selects a reference sample, and performs normalisation and scaling.
   
     tar_target(
      # pseudobulk_sample ------
      pseudobulk_sample,
      {
        message('TAR: pseudobulk_sample START')
        system(glue("~/bin/rclone copy box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/pseudobulk_sample_is_immune {getwd()}/"))
        
        se = 
          loadHDF5SummarizedExperiment(hdf5_path) |> 
          filter(is_gene_shared) |> 
          
          #---------------------------------#
          # Edit or add more filters here for analyses
          #---------------------------------#
          filter(is_immune & do_analyse) 
        

        # Age
        se = 
          se |> 
          filter(age_days > 365) |> 
          mutate(age_years = age_days / 365) |> 
          mutate(age_bin = dplyr::case_when(
            age_years < 3 ~ "Infancy",
            age_years < 12 ~ "Childhood",
            age_years < 20 ~ "Adolescence",
            age_years < 40 ~ "Young Adulthood",
            age_years < 50 ~ "Middle Age",
            age_years < 60 ~ "Senior_50",
            age_years < 70 ~ "Senior_60",
            age_years >= 70 ~ "Senior_70",
            TRUE ~ NA_character_
          )) |> 
          mutate(age_decade = ceiling(age_years/10) |> as.integer() |> as.character())  
          
        # Ethnicity
        # THIS CHANGE WAS SUGGESTED BY THE REVIEWER
        se = 
          se |> 
          mutate(ethnicity_groups = if_else(ethnicity_groups=="Japanese", "East Asian", ethnicity_groups)) |> 
          mutate(ethnicity_groups = if_else(ethnicity_groups=="Hispanic/Latin American", "Hispanic", ethnicity_groups))
        
        # CHEN'S PIPELINE START
          
        # Filter samples that have enough genes > 0 but not too many
        samples_with_right_number_of_detected_genes =
          (se |> assay() > 0) |>
          colSums() |>
          divide_by(nrow(se)) |>
          dplyr::between(0.1, 1)
        
        se = se[, samples_with_right_number_of_detected_genes]
        
        cli::cli_alert_info("\nCalculating reference sample for scaling gene counts...\n")
        
        # Compute mean library size
        mean_library_size <- se |>
          assay("counts") |>
          _[nrow(se) |> seq_len() |> sample(size = 2000), ] |>
          colSums() |>
          mean()
        
        # Optional: retrieve the sample name (column name in the SummarizedExperiment)
        reference_index <-
          se |>
          assay("counts") |>
          colSums() |>
          {
            \(x) abs(x - mean_library_size)
          }() |>                          # Calculate absolute difference from the mean
          which.min()                     # Identify the smallest difference
        reference_sample <- colnames(se)[reference_index]
        
        reference_path = glue::glue(
          "{targets::tar_config_get('store')}/reference_sample.rds"
        )
        saveRDS(reference_sample, file = reference_path)
        cli::cli_alert_info("\nReference sample saved: {file.exists(reference_path)}")
        
        design =
          se |>
          
          # Discretise the age for the following operation
          # mutate(is_old_individual = age_days > 50 * 365) |>
          
          # This is to resolve some confounders to preserve the genes.
          # In this case we care about data variability, not the actual meaning of the variables
          tidybulk::resolve_complete_confounders_of_non_interest(ethnicity_groups, dataset_id) |> 
          colData() |>
          droplevels() |>
          model.matrix( ~ ethnicity_groups___altered, data = _)
        
        h <- hat(design)
        # MinSampleSize <- 1 / max(h)
        # 
        # counts_matrix <- assay(se, "counts")
        # lib.size <- colSums(counts_matrix)
        # CPM <- edgeR::cpm(counts_matrix, lib.size = lib.size)
        # 
        # quantile_cpm = CPM %>%
        #   apply(
        #     1, quantile, 
        #     probs =  1 - MinSampleSize %>%
        #       ceiling() %>%
        #       {. / ncol(se)}
        #   )
        # 
        # mini_cpm_threshold = quantile_cpm %>% quantile(0.5) %>% unname()
        # cli::cli_alert_info("\nMini_cpm_threshold = {mini_cpm_threshold}")
        # 
        se =
          se |>
          keep_abundant(
            design = design,
            minimum_count_per_million = 70
          ) |>
          
          # Get scaling factor
          scale_abundance(method = "TMMwsp", reference_sample = reference_sample) |>

        # CHEN'S PIPELINE END
          
        
          # Drop sex unknown as causes problem during fit
          mutate(
            sex = if_else(sex |> is.na(), "unknown", sex),
            ethnicity_groups = if_else(ethnicity_groups |> is.na(), "Other/Unknown", ethnicity_groups)
          ) |> 
          filter(sex != "unknown") |> 
          filter(!age_bin |> is.na()) |> 
          
          # Eliminate complete confounders
          tidybulk:::resolve_complete_confounders_of_non_interest(assay_groups, dataset_id, disease_groups) |> 
          
          # sibrary size factor is the reciproque of the multiplier (correction factor)
          mutate(offset = log(1/multiplier)) |> 
          
          # Set intercept
          mutate(
            ethnicity_groups = fct_relevel(ethnicity_groups, "European"),
            assay_groups___altered = fct_relevel(assay_groups___altered, "10x Genomics 3"),
            disease_groups___altered = fct_relevel(disease_groups___altered, "Normal"),
            age_bin = fct_relevel(age_bin, "Adolescence")
          ) 
        
        # # Add dispersion
        # rowData(se)  = 
        #   rowData(se) |> 
        #   as_tibble(rownames = ".feature") |> 
        #   left_join(glmGamPoi_overdispersions |> enframe(name = ".feature", value = "dispersion")) |> 
        #   data.frame(row.names = ".feature") |> DataFrame()
        
        message('TAR: pseudobulk_sample COMPLETE')
        se
        
      }, 
      packages = c("tidybulk", "HDF5Array", "tidySummarizedExperiment", "magrittr", "tibble", "forcats"),
      resources = tar_resources(crew = tar_resources_crew("elastic_big")),
      memory = "persistent", 
      error = "stop"
    ),
    
    # pseudobulk_sample_id ------
    # This target extracts unique sample ids from the pseudobulk sample  
    tar_target(
      pseudobulk_sample_id,
      pseudobulk_sample |> colnames(),
      packages = c( "tidySummarizedExperiment", "targets", "purrr", "dplyr"),
      resources = tar_resources(crew = tar_resources_crew("elastic"))
    ),
    
    # feature_df ------
    # This target extracts unique features from the pseudobulk sample and groups them into 
    # chunks for parallel processing.
    tar_target(
      feature_df, 
      pseudobulk_sample |> 
        distinct(.feature)|>
        # testing genes that ran for long time
        # filter(.feature %in% c(
        #   "ENSG00000175274", "ENSG00000213221", "ENSG00000101405", "ENSG00000003756",
        #   "ENSG00000236859", "ENSG00000104825", "ENSG00000203497", "ENSG00000143110",
        #   "ENSG00000077454", "ENSG00000104231"
        # )) |>
        # head(1) %>% 
        group_by(.feature) |> 
        tar_group(), 
      iteration = "group",
      packages = c( "tidySummarizedExperiment", "targets", "purrr", "dplyr"),
      resources = tar_resources(crew = tar_resources_crew("elastic"))
    ),
    
    # se_df -----
    # This target creates a list-column of SummarizedExperiment objects,
    # with each object corresponding to a distinct feature.
    tar_target(
      se_df, 
      feature_df |> mutate(se = map(.feature, ~ 
                                      pseudobulk_sample[.x, , drop=FALSE]
      ))  , 
      pattern = map(feature_df),
      packages = c( "brms", "glue"),
      resources = tar_resources(crew = tar_resources_crew("elastic"))
    ),
    
    # estimates_chunk ------
    # This target fits Bayesian models on chunks of the data. It processes each feature's data, handles missing values,
    # defines the model specification with priors, and runs the Bayesian inference using the brm function.
    tar_target(
      estimates_chunk, 
      
      se_df |> mutate(brms_fit = map(se, ~ {
        
        data = 
          .x |>
          as_tibble() |> 
          mutate(counts = counts |> as.integer()) |> 
          droplevels()
        
        # Drop NA counts. Not sure why they are there. E.g.:
        # $ .feature             <chr> "ENSG00000134419"
        # $ .sample              <chr> "3bfa31867cc1c823e0cb2f1ff24df26b___1"
        # $ counts               <int> NA
        # $ gene_presence        <int> 25
        # $ counts_scaled        <dbl> 38316.51
        # $ sample_id            <chr> "3bfa31867cc1c823e0cb2f1ff24df26b"
        # $ is_immune            <int> 1
        # $ do_analyse           <lgl> TRUE
        # $ donor_id             <chr> "one_Ten"
        # $ title                <chr> "Individual Single-Cell RNA-seq PBMC Data from Schulte-Schrepping et al."
        # $ dataset_id           <chr> "5e717147-0f75-4de1-8bd2-6fda01b8d75f"
        # $ collection_id        <chr> "b9fc3d70-5a72-4479-a046-c2cc1ab19efc"
        # $ age_days             <int> 10585
        # $ sex                  <chr> "male"
        # $ ethnicity_groups     <fct> European
        # $ tissue_groups        <chr> "blood"
        # $ assay_groups         <fct> 10x Genomics 3
        # $ disease_groups       <fct> Normal
        # $ age_bin <fct> Young Adulthood
        # $ TMM                  <dbl> 2.538842
        # $ multiplier           <dbl> 1.394764e-05
        # $ offset               <dbl> -11.1802
        # $ is_gene_shared       <lgl> TRUE
        # $ .abundant            <lgl> TRUE
        # $ dispersion           <dbl> 0.4960158
        n_NAs = data |> filter(counts |> is.na()) |> nrow()
        if(n_NAs > 0){
          warning(glue("You have {n_NAs} NAs in counts. They have been filtered out"))
          data = 
            data |> 
            filter(!counts |> is.na()) |> 
            droplevels()
        }
        
        # Manually revise data colnames to suit brms bug
        colnames(data) = colnames(data) |> stringr::str_replace_all("_+", "_")
        
        # # Check if dispersion estimation has failed
        # if(data |> pull(dispersion) |> unique() |> is.na()){
        #   warning("The dispersion calculation has failed. 1 is given as default prior.")
        #   data = data |> mutate(dispersion = 1)
        # }
        
        # Define the model formula
        formula <- bf(
          
          # Formula for counts
          counts ~ 1 + offset(offset) + age_bin*sex + disease_groups_altered + ethnicity_groups + assay_groups_altered + 
            (1 | dataset_id_altered) + 
            (1 + age_bin*sex + ethnicity_groups | tissue_groups),
          
          # Formula for dispersion
          shape ~ 1 + disease_groups_altered + assay_groups_altered + ethnicity_groups + (1 | tissue_groups)  # Model 'shape' as a function of scaled 'disp'
          
          # Using the externally, eBayes inferred overdispersion
          # shape ~ 1 + offset(log(1/dispersion))
        )
        
        # prior = c(
        #   prior(normal(i, 5), class = Intercept),
        #   prior(normal(0, 2), class = Intercept, dpar = shape),
        #   prior(normal(0, 5), class = b),
        #   prior(normal(0, 2), class = b, dpar = shape)
        # ) |> 
        #   substitute(env = list(i = mean(log1p(data$counts / exp(data$offset))))) |> 
        #   eval()
        
        # prior = prior(normal(-0.0002056948, 0.07690437))
        
        # HPC pipeline: param V1:
        # prior = c(
        #   prior(student_t(4.45496, 0.008599254, 1.143344), class = "b"),
        #   prior(student_t(18.16242, 0.07952513, 0.9926044), class = "b", dpar = "shape"),
        #   prior(normal(5.441626, 2.25460683), class = "Intercept"),
        #   prior(normal(0.1459487, 0.8347875), class = "Intercept", dpar = "shape"),
        #   prior(student_t(4.7655009	, 0.887529, 0.8684176), class = "sd", lb = 0),
        #   prior(student_t(53.08894, 0.9080073, 0.2870678), class = "sd", dpar = "shape", lb = 0),
        #   prior(beta(0.5541155, 9.337894), class = "zi")
        # ) 
        
        # HPC pipeline: param V2:
        prior = c(
          prior(student_t(6.153327, 0.06161134, 0.9263627), class = "b"),
          prior(student_t(40.51669, 0.07603337, 0.8252114), class = "b", dpar = "shape"),
          prior(normal(6.057503, 2.438534), class = "Intercept"),
          prior(normal(0.4260793, 1.470536), class = "Intercept", dpar = "shape"),
          prior(student_t(52.19541	, 0.5703259, 0.4147664), class = "sd", lb = 0),
          prior(normal(0.8670409, 0.1779553), class = "sd", dpar = "shape", lb = 0),
          prior(beta(0.5381488, 10.3577433), class = "zi", lb = 0, ub = 1)
        )
        
        chains = 2
        
        Kc <- 39
        Kc_shape <- 28
        M_1 <- 1; N_1 <- 105
        M_2 <- 19; N_2 <- 26
        M_3 <- 1; N_3 <- 26
        
        inits <- lapply(1:chains, function(i) {
          list(
            # Fixed effects for count part
            b = 0.06161134 + 0.9263627 * rt(Kc, 6.153327),
            Intercept = rnorm(1, 6.057503, 2.438534),
            
            # Fixed effects for shape submodel
            b_shape = 0.07603337 + 0.8252114 * rt(Kc_shape, 40.51669),
            Intercept_shape = rnorm(1, 0.4260793, 1.470536),
            
            # Zero-inflation probability
            zi = rbeta(1, 0.5381488, 10.3577433),
            
            # Group-level standard deviations and effects
            sd_1 = abs(0.5703259 + 0.4147664 * rt(M_1, 52.19541)),      # count
            z_1 = replicate(M_1, rnorm(N_1, mean = 0 , sd = 0.08547970), simplify = FALSE),
            
            sd_2 = abs(0.5703259 + 0.4147664 * rt(M_2, 52.19541)),      # zi
            z_2 = matrix(rnorm(M_2 * N_2, mean = 0 , sd = 0.08547970), nrow = M_2, ncol = N_2),
            L_2 = diag(M_2),                                            # no correlation (identity)
            
            sd_3 = abs(rnorm(M_3, 0.8670409, 0.1779553)),               # shape
            z_3 = replicate(M_3, rnorm(N_3, mean = 0 , sd = 0.08547970), simplify = FALSE)
          )
        })
        
        # chains = 2
        # inits <- list(Intercept = mean(log1p(data$counts / exp(data$offset))))
        # inits <- replicate(chains, inits, simplify = FALSE)
        
        
        brm(
          formula = formula,
          data = data,
          family = zero_inflated_negbinomial(),
          prior = prior,
          chains = chains,
          cores = pmax(as.numeric(parallelly::availableCores()), 2), #, threads = 2,
          warmup = 300, 
          refresh = 10,
          backend = "cmdstanr", 
          #sparse = TRUE,
          #save_model = glue("{external_directory}~/temp.rds"),
          #algorithm = "pathfinder",
          init = inits,
          iter = 400  # Increase iterations for better convergence
        )
        
      })) |> 
        
        # Drop data because it is withn the brms object
        select(-se), 
      pattern = map(se_df),
      packages = c( "brms", "glue", "stringr", "dplyr", "purrr", "SummarizedExperiment", "tidySummarizedExperiment"),
      resources = tar_resources(crew = tar_resources_crew("elastic")),
      cue = tar_cue(mode = "never")
      
    ),
    
    ## summary ----- 
    # This target summarises the fitted Bayesian models by performing hypothesis tests for ethnicity contrasts 
    # and extracting convergence diagnostics (Rhat) for the ethnicity parameters.
    tar_target(
      summary,
      estimates_chunk |>
        mutate(summary = map(brms_fit, ~ .x |> hypothesis(
          c(
            "Europeans" = "(ethnicity_groupsAfrican
    + ethnicity_groupsEastAsian
    + ethnicity_groupsHispanicDLatinAmerican
    + ethnicity_groupsSouthAsian
    + `ethnicity_groupsJapanese`) / 5 = 0",
            "EastAsian" = "(
       ethnicity_groupsAfrican
     + ethnicity_groupsHispanicDLatinAmerican
     + ethnicity_groupsSouthAsian
     + `ethnicity_groupsJapanese`
     - 5 * ethnicity_groupsEastAsian
     ) / 5 = 0",
            "SouthAsian" = "(
       ethnicity_groupsAfrican
     + ethnicity_groupsHispanicDLatinAmerican
     + ethnicity_groupsEastAsian
     + `ethnicity_groupsJapanese`
     - 5 * ethnicity_groupsSouthAsian
     ) / 5 = 0",
            "African" = "(
       ethnicity_groupsEastAsian
     + ethnicity_groupsHispanicDLatinAmerican
     + ethnicity_groupsSouthAsian
     + `ethnicity_groupsJapanese`
     - 5 * ethnicity_groupsAfrican
     ) / 5 = 0",
            "HispanicDLatinAmerican" = "(
       ethnicity_groupsAfrican
     + ethnicity_groupsEastAsian
     + ethnicity_groupsSouthAsian
     + `ethnicity_groupsJapanese`
     - 5 * ethnicity_groupsHispanicDLatinAmerican
     ) / 5 = 0",
            
            "Japanese" = "(
       ethnicity_groupsAfrican
     + ethnicity_groupsHispanicDLatinAmerican
     + ethnicity_groupsSouthAsian
     + ethnicity_groupsEastAsian
     - 5 * `ethnicity_groupsJapanese`
     ) / 5 = 0"
          ),
          #        c(
          #          "African" = "(ethnicity_groupsEuropean
          #  + ethnicity_groupsEastAsian
          #  + ethnicity_groupsHispanicDLatinAmerican
          #  + ethnicity_groupsSouthAsian
          #  + `ethnicity_groupsJapanese`) / 5 = 0",
          #          
          #          "Europeans" = "(
          #   ethnicity_groupsEastAsian
          # + ethnicity_groupsHispanicDLatinAmerican
          # + ethnicity_groupsSouthAsian
          # + `ethnicity_groupsJapanese`
          # - 5 * ethnicity_groupsEuropean
          # ) / 5 = 0",
          #          
          #          "EastAsian" = "(
          #   ethnicity_groupsEuropean
          # + ethnicity_groupsHispanicDLatinAmerican
          # + ethnicity_groupsSouthAsian
          # + `ethnicity_groupsJapanese`
          # - 5 * ethnicity_groupsEastAsian
          # ) / 5 = 0",
          #          
          #          "SouthAsian" = "(
          #   ethnicity_groupsEuropean
          # + ethnicity_groupsHispanicDLatinAmerican
          # + ethnicity_groupsEastAsian
          # + `ethnicity_groupsJapanese`
          # - 5 * ethnicity_groupsSouthAsian
          # ) / 5 = 0",
          #          
          #          "HispanicDLatinAmerican" = "(
          #   ethnicity_groupsEuropean
          # + ethnicity_groupsEastAsian
          # + ethnicity_groupsSouthAsian
          # + `ethnicity_groupsJapanese`
          # - 5 * ethnicity_groupsHispanicDLatinAmerican
          # ) / 5 = 0",
          #          
          #          "Japanese" = "(
          #   ethnicity_groupsEuropean
          # + ethnicity_groupsHispanicDLatinAmerican
          # + ethnicity_groupsSouthAsian
          # + ethnicity_groupsEastAsian
          # - 5 * `ethnicity_groupsJapanese`
          # ) / 5 = 0"
          #        ),
          
          # Median instead and mad of mean and sd
          robust=TRUE,
          alpha = 0.1
        )
        )) |>
        
        mutate(
          
          summary_tissue = map(
            brms_fit,  function(x) {
              
              params = x$fit %>% summary() |> _[[1]] |> rownames()
              params = params[grepl("^r_tissue_groups\\[.*?,Intercept\\]$", params)] %>% sub("^r_", "", .) %>% paste0("`", . , "`")
              tissue_names <- sub("`tissue_groups\\[(.*),Intercept\\]`", "\\1", params)
              
              equations <- sapply(seq_along(params), function(i) {
                this_tissue <- tissue_names[i]
                this_param <- params[i]
                other_params <- params[-i]
                avg_expr <- paste0("(", paste(other_params, collapse = " + "), ")/", length(other_params))
                eq <- paste0(this_param, " - ", avg_expr, " = 0")
                eq
              })
              names(equations) <- tissue_names
              
              return(
                x |> hypothesis(equations, class = "r", robust=TRUE, alpha = 0.1)
              )
              
            }
          )
          
        ) %>% 
        
        mutate(Rhat_ethnicity = map_dbl(brms_fit, 
                                        ~ summary(.x)$fixed |> 
                                          as_tibble(rownames = "par") |> 
                                          filter(par |> str_detect("ethnicity")) |> 
                                          pull(Rhat) |>
                                          max()
        )) |> 
        
        mutate(Rhat_tissue = map_dbl(brms_fit, 
                                     ~ summary(.x)$random$tissue_groups |> 
                                       as_tibble() |> 
                                       pull(Rhat) |>
                                       max()
        )) %>%  
        
        select(-brms_fit),
      
      pattern = map(estimates_chunk),
      packages = c( "brms", "glue", "dplyr", "purrr", "rstan", "magrittr", "stringr"),
      resources = tar_resources(crew = tar_resources_crew("elastic"))
    ),
    
    ## effect_removed -----
    # This target generates adjusted model estimates by removing unwanted effects from the fitted Bayesian models,
    # thereby isolating the effects of interest. Here, nuisance covariates are set to NA and removed from the predictions.
    # This target produces adjusted estimates from the Bayesian models, removing unwanted effects while retaining 
    # the tissue group random effect, thus preserving variability associated with tissue-specific factors.
    tar_target(
      effect_removed, 
      estimates_chunk |> 
        mutate(brms_fit_adjusted_ethnicity = map(brms_fit, ~ .x |> remove_unwanted_effect(
          newdata = .x$data |> mutate(assay_groups_altered=NA, sex = NA, age_bin = NA, disease_groups_altered = NA, dataset_id_altered = NA), # age_bin*sex + disease_groups + ethnicity_groups + assay_groups
          robust = TRUE, 
          re_formula = ~ 0
        ))) |> 
        mutate(brms_fit_adjusted_ethnicity_new = map(brms_fit, ~ .x |> remove_unwanted_effect_new(
          newdata = .x$data |> mutate(assay_groups_altered=NA, sex = NA, age_bin = NA, disease_groups_altered = NA, dataset_id_altered = NA), # age_bin*sex + disease_groups + ethnicity_groups + assay_groups
          robust = FALSE, correct_by_offset = FALSE,
          re_formula = ~ 0
        ))) |> 
        mutate(brms_fit_adjusted_tissue = map(brms_fit, ~ .x |> remove_unwanted_effect(
          newdata = .x$data |> mutate(assay_groups=NA, sex = NA, age_bin = NA, disease_groups = NA, ethnicity_groups = NA), # age_bin*sex + disease_groups + ethnicity_groups + assay_groups
          robust = TRUE, 
          re_formula = ~ (1 | tissue_groups)
        ))) |> 
        mutate(brms_fit_adjusted_tissue_new = map(brms_fit, ~ .x |> remove_unwanted_effect_new(
          newdata = .x$data |> mutate(assay_groups_altered=NA, ethnicity_groups = NA, sex = NA, age_bin = NA, disease_groups_altered = NA, dataset_id_altered = NA), # age_bin*sex + disease_groups + ethnicity_groups + assay_groups
          robust = FALSE, correct_by_offset = FALSE,
          re_formula = ~ (1 | tissue_groups)
        ))) |> 
        mutate(brms_fit_adjusted_ethnicity_estimate = map(brms_fit_adjusted_ethnicity, ~ {
          
          df = .x |> as_tibble()
          if (nrow(df) == length(pseudobulk_sample_id)){
            return(
              df |>
                select(adjusted___Estimate) |> #, adjusted___Q2.5, adjusted___Q97.5) |> 
                mutate(sample_id = pseudobulk_sample_id)
            )
          }else{
            return(NULL)
          }
          
        })) |> 
        mutate(brms_fit_adjusted_ethnicity_new_estimate = map(brms_fit_adjusted_ethnicity_new, ~ {
          
          df = .x |> as_tibble()
          if (nrow(df) == length(pseudobulk_sample_id)){
            return(
              df |>
                select(adjusted___Estimate) |> #, adjusted___Q2.5, adjusted___Q97.5) |> 
                mutate(sample_id = pseudobulk_sample_id)
            )
          }else{
            return(NULL)
          }
          
        })) |> 
        mutate(brms_fit_adjusted_tissue_estimate = map(brms_fit_adjusted_tissue, ~ {
          
          df = .x |> as_tibble()
          if (nrow(df) == length(pseudobulk_sample_id)){
            return(
              df |>
                select(adjusted___Estimate) |> #, adjusted___Q2.5, adjusted___Q97.5) |> 
                mutate(sample_id = pseudobulk_sample_id)
            )
          }else{
            return(NULL)
          }
          
        })) |>
        mutate(brms_fit_adjusted_tissue_new_estimate = map(brms_fit_adjusted_tissue_new, ~ {
          
          df = .x |> as_tibble()
          if (nrow(df) == length(pseudobulk_sample_id)){
            return(
              df |>
                select(adjusted___Estimate) |> #, adjusted___Q2.5, adjusted___Q97.5) |> 
                mutate(sample_id = pseudobulk_sample_id)
            )
          }else{
            return(NULL)
          }
          
        })) |>
        select(-brms_fit),
      
      pattern = map(estimates_chunk),
      packages = c( "brms", "glue", "dplyr", "purrr", "rstan"),
      resources = tar_resources(crew = tar_resources_crew("elastic"))
    ),
    
    # adjusted_matrix -----
    tar_target(
      adjusted_assay_ethnicity,
      get_adjusted_matrix(effect_removed, brms_fit_adjusted_ethnicity_estimate),
      packages = c( "brms", "glue", "dplyr", "purrr", "rstan", "magrittr", "stringr", "tidySummarizedExperiment") ,
      resources = tar_resources(
        crew = tar_resources_crew("elastic_big_30_cores")
      )
    ),
    
    tar_target(
      adjusted_assay_ethnicity_new,
      get_adjusted_matrix(effect_removed, brms_fit_adjusted_ethnicity_new_estimate),
      packages = c( "brms", "glue", "dplyr", "purrr", "rstan", "magrittr", "stringr", "tidySummarizedExperiment") ,
      resources = tar_resources(
        crew = tar_resources_crew("elastic_big_30_cores")
      )
    ),
    
    tar_target(
      adjusted_assay_tissue,
      get_adjusted_matrix(effect_removed, brms_fit_adjusted_tissue_estimate),
      packages = c( "brms", "glue", "dplyr", "purrr", "rstan", "magrittr", "stringr", "tidySummarizedExperiment") ,
      resources = tar_resources(
        crew = tar_resources_crew("elastic_big_30_cores")
      )
    ),
    
    tar_target(
      adjusted_assay_tissue_new,
      get_adjusted_matrix(effect_removed, brms_fit_adjusted_tissue_new_estimate),
      packages = c( "brms", "glue", "dplyr", "purrr", "rstan", "magrittr", "stringr", "tidySummarizedExperiment") ,
      resources = tar_resources(
        crew = tar_resources_crew("elastic_big_30_cores")
      )
    )
    
  ) # end ) of all target list
  
  
}, ask = FALSE)
