library(targets)

setwd("ethnicity_prediction")


job::job({
  library(targets)
  
  tar_script({
    
    library(tidyverse)
    library(targets)
    library(tarchetypes)
    library(glue)
    library(qs)
    library(crew)
    library(crew.cluster)
    
    # Helper (optional) to avoid repetition
    new_elastic <- function(name, mem_gb, time_min, workers, crashes_max, cpus_per_task = 2, backup = NULL) {
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
    
    # Small → large, with fallbacks to the next size up
    elastic_160 <- new_elastic("elastic_160", 160, 60 * 24, workers = 8,  crashes_max = 2)
    elastic_80  <- new_elastic("elastic_80",   80,  60 * 4,  workers = 16, crashes_max = 1, cpus_per_task = 8, backup = elastic_160)
    elastic_40  <- new_elastic("elastic_40",   40,  60 * 4,  workers = 24, crashes_max = 1, cpus_per_task = 8, backup = elastic_80)
    elastic_20  <- new_elastic("elastic_20",   20,  60 * 4,  workers = 32, crashes_max = 1, cpus_per_task = 8, backup = elastic_40)
    elastic_10  <- new_elastic("elastic_10",   10,  60 * 4,  workers = 48, crashes_max = 1, cpus_per_task = 8, backup = elastic_20)
    elastic_5   <- new_elastic("elastic_5",     5, 60 * 4,  workers = 150, crashes_max = 6, cpus_per_task = 8, backup = elastic_10)
    
    elastic_5_minimal   <- new_elastic("elastic_5_minimal",     5, 60 * 4,  workers = 300, crashes_max = 6, cpus_per_task = 2, backup = elastic_10)
    
    
    # Group for targets (small → large)
    controllers <- crew_controller_group(
      elastic_5, elastic_10, elastic_20, elastic_40, elastic_80, elastic_160, elastic_5_minimal
    )
    
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
      
      controller = controllers, 
      resources = tar_resources(
        crew = tar_resources_crew(controller = "elastic_5_minimal")
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
    remove_unwanted_effect_new = function(fit,
                                          newdata,
                                          robust = FALSE,
                                          correct_by_offset = T,
                                          re_formula = ~ 0) {
      # Calculate residuals: observed counts minus fitted values, normalised by exp(offset)
      # This places residuals on a consistent scale, making them addable to adjusted predictions later.
      fitted_residuals =   fit |> predictive_error(robust = robust,
                                                   summary = FALSE,
                                                   offset = 0)
      
      # Correct by offset
      if (correct_by_offset)
        fitted_residuals = fitted_residuals |>
          sweep(2, fit$data$offset |> exp(), FUN = "/")
      
      # Extract fitted values for the specified factor only, removing random effects by setting re_formula = ~0
      # 'resp = factor' focuses on the selected response variable (factor)
      fitted_values_ethnicity <- posterior_epred(fit,
                                                 newdata = newdata,
                                                 re_formula = re_formula,
                                                 offset = 0)
      
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
    
    remove_unwanted_effect = function(fit,
                                      newdata,
                                      robust = FALSE,
                                      correct_by_offset = T,
                                      re_formula = ~ 0) {
      # Calculate residuals: observed counts minus fitted values, normalised by exp(offset)
      # This places residuals on a consistent scale, making them addable to adjusted predictions later.
      fitted_residuals =   fit |> residuals(robust = robust, summary = FALSE)
      
      # Correct by offset
      if (correct_by_offset)
        fitted_residuals = fitted_residuals |>
          sweep(2, fit$data$offset |> exp(), FUN = "/")
      
      # Extract fitted values for the specified factor only, removing random effects by setting re_formula = ~0
      # 'resp = factor' focuses on the selected response variable (factor)
      fitted_values_ethnicity <- fitted(
        fit,
        newdata = newdata,
        re_formula = re_formula,
        summary = FALSE,
        offset = 0
      )
      
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
          hdf5_path = glue("{targets::tar_config_get('store')}/pseudobulk_sample_is_immune")
          
          
          system(glue("~/bin/rclone copy box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/pseudobulk_sample_is_immune {hdf5_path}/"))
          
          
          se = 
            loadHDF5SummarizedExperiment(hdf5_path) |> 
            filter(is_gene_shared) |> 
            
            #---------------------------------#
            # Edit or add more filters here for analyses
            #---------------------------------#
            filter(is_immune & do_analyse) |>
            
            
            #---------------------------------#
            # THIS IS A TEST FILTERING TO SEE IF RERUNNING THE MODEL
            # ON HIGH CONFIDENCE DATA GIVES BETTER PREDICTION
            # PERFORMANCE
            #---------------------------------#
            left_join(read_csv(low_consensus_filtering_Ning), by = "sample_id") |> 
            dplyr::rename(filter_out = Filter) |>
            filter(!filter_out)  
            
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
        packages = c("tidybulk", "HDF5Array", "tidySummarizedExperiment", "magrittr", "tibble", "forcats", "glue"),
        memory = "persistent", 
        error = "stop", 
        resources = tar_resources(crew = tar_resources_crew("elastic_80"))
      ),
      
      tar_target(
        low_consensus_filtering_Ning, 
        "~/labHead/cellNexus_article/ethnicity_prediction/sample_filtering_based_on_ethnicity_clearness.csv", 
        format = "file"
        ),
      
      # pseudobulk_sample_id ------
      # This target extracts unique sample ids from the pseudobulk sample  
      tar_target(
        pseudobulk_sample_id,
        pseudobulk_sample |> colnames(),
        packages = c( "tidySummarizedExperiment", "targets", "purrr", "dplyr")
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
        packages = c( "tidySummarizedExperiment", "targets", "purrr", "dplyr")
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
        packages = c( "brms", "glue")
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
            prior(student_t(3, i, 1.5), class = Intercept),
            prior(student_t(3, 0, 1), class = Intercept, dpar = shape),
            prior(student_t(3, 0, 5), class = b),
            prior(student_t(3, 0, 2), class = b, dpar = shape)
            # prior(beta(0.5381488, 10.3577433), class = "zi", lb = 0, ub = 1) # addition zi from V2
          ) |>
            substitute(env = list(i = mean(log1p(
              data$counts / exp(data$offset)
            )))) |>
            eval()
          
          chains = 2
          
          # dynamically extract param from stan data
          # code used from brm
          bterms <- brmsterms(
            formula = brms:::validate_formula(
              formula,
              data = data,
              family = zero_inflated_negbinomial(),
              autocor = NULL,
              sparse = NULL,
              cov_ranef = NULL
            )
          )
          bframe <- brms:::brmsframe(bterms, data)
          sdata <- brms:::.standata(
            bframe,
            data = data,
            prior = prior,
            data2 = NULL,
            stanvars = NULL,
            threads = NULL
          )
          
          Kc <- sdata$Kc
          Kc_shape <- sdata$Kc_shape
          M_1 <- sdata$M_1
          N_1 <- sdata$N_1
          M_2 <- sdata$M_2
          N_2 <- sdata$N_2
          M_3 <- sdata$M_3
          N_3 <- sdata$N_3
          
          inits <- lapply(1:chains, function(i) {
            list(
              #### revert v0 prior
              # Fixed effects for count part
              b = rnorm(Kc, 0, 5),
              # dynamically set mu for intercept
              Intercept = rnorm(1, mean(log1p(
                data$counts / exp(data$offset)
              )), 1.5),
              
              # Fixed effects for shape submodel
              b_shape = rnorm(Kc_shape, 0, 2),
              Intercept_shape = rnorm(1, 0, 1)
            
            )
          })
          
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
    + ethnicity_groupsHispanic
    + ethnicity_groupsSouthAsian) / 4 = 0",
              "EastAsian" = "(
       ethnicity_groupsAfrican
     + ethnicity_groupsHispanic
     + ethnicity_groupsSouthAsian
     - 4 * ethnicity_groupsEastAsian
     ) / 4 = 0",
              "SouthAsian" = "(
       ethnicity_groupsAfrican
     + ethnicity_groupsHispanic
     + ethnicity_groupsEastAsian
     - 4 * ethnicity_groupsSouthAsian
     ) / 4 = 0",
              "African" = "(
       ethnicity_groupsEastAsian
     + ethnicity_groupsHispanic
     + ethnicity_groupsSouthAsian
     - 4 * ethnicity_groupsAfrican
     ) / 4 = 0",
              "HispanicDLatinAmerican" = "(
       ethnicity_groupsAfrican
     + ethnicity_groupsEastAsian
     + ethnicity_groupsSouthAsian
     - 4 * ethnicity_groupsHispanic
     ) / 4 = 0"           
      ),
     
            
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
        packages = c( "brms", "glue", "dplyr", "purrr", "rstan", "magrittr", "stringr")
      ),
      
      ## effect_removed -----
      # This target generates adjusted model estimates by removing unwanted effects from the fitted Bayesian models,
      # thereby isolating the effects of interest. Here, nuisance covariates are set to NA and removed from the predictions.
      # This target produces adjusted estimates from the Bayesian models, removing unwanted effects while retaining 
      # the tissue group random effect, thus preserving variability associated with tissue-specific factors.
      tar_target(
        effect_removed, 
        estimates_chunk |> 
          mutate(brms_fit_adjusted_ethnicity_new = map(brms_fit, ~ .x |> remove_unwanted_effect_new(
            newdata = .x$data |> mutate(assay_groups_altered=NA, sex = NA, age_bin = NA, disease_groups_altered = NA, dataset_id_altered = NA), # age_bin*sex + disease_groups + ethnicity_groups + assay_groups
            robust = FALSE, correct_by_offset = FALSE,
            re_formula = ~ 0
          ))) |> 
          # mutate(brms_fit_adjusted_tissue_new = map(brms_fit, ~ .x |> remove_unwanted_effect_new(
          #   newdata = .x$data |> mutate(assay_groups_altered=NA, ethnicity_groups = NA, sex = NA, age_bin = NA, disease_groups_altered = NA, dataset_id_altered = NA), # age_bin*sex + disease_groups + ethnicity_groups + assay_groups
          #   robust = FALSE, correct_by_offset = FALSE,
          #   re_formula = ~ (1 | tissue_groups)
          # ))) |> 
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
          
          
          # mutate(brms_fit_adjusted_tissue_new_estimate = map(brms_fit_adjusted_tissue_new, ~ {
          #   
          #   df = .x |> as_tibble()
          #   if (nrow(df) == length(pseudobulk_sample_id)){
          #     return(
          #       df |>
          #         select(adjusted___Estimate) |> #, adjusted___Q2.5, adjusted___Q97.5) |> 
          #         mutate(sample_id = pseudobulk_sample_id)
          #     )
          #   }else{
          #     return(NULL)
          #   }
          #   
          # })) |>
          
          mutate(brms_fit_adjusted_ethnicity = map(brms_fit, ~ .x |> remove_unwanted_effect(
            newdata = .x$data |> mutate(assay_groups_altered=NA, sex = NA, age_bin = NA, disease_groups_altered = NA, dataset_id_altered = NA), # age_bin*sex + disease_groups + ethnicity_groups + assay_groups
            robust = FALSE, correct_by_offset = FALSE,
            re_formula = ~ 0
          ))) |> 
          # mutate(brms_fit_adjusted_tissue_new = map(brms_fit, ~ .x |> remove_unwanted_effect_new(
          #   newdata = .x$data |> mutate(assay_groups_altered=NA, ethnicity_groups = NA, sex = NA, age_bin = NA, disease_groups_altered = NA, dataset_id_altered = NA), # age_bin*sex + disease_groups + ethnicity_groups + assay_groups
          #   robust = FALSE, correct_by_offset = FALSE,
          #   re_formula = ~ (1 | tissue_groups)
          # ))) |> 
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
          
          select(-brms_fit),
        
        pattern = map(estimates_chunk),
        packages = c( "brms", "glue", "dplyr", "purrr", "rstan"), 
        resources = tar_resources(crew = tar_resources_crew("elastic_5_minimal"))
      ),
      
      tar_target(
        adjusted_assay_ethnicity_new,
        get_adjusted_matrix(effect_removed, brms_fit_adjusted_ethnicity_new_estimate),
        packages = c( "brms", "glue", "dplyr", "purrr", "rstan", "magrittr", "stringr", "tidySummarizedExperiment") 
      ),
      tar_target(
        adjusted_assay_ethnicity,
        get_adjusted_matrix(effect_removed, brms_fit_adjusted_ethnicity_estimate),
        packages = c( "brms", "glue", "dplyr", "purrr", "rstan", "magrittr", "stringr", "tidySummarizedExperiment") 
      )
      # ,
      # 
      # tar_target(
      #   adjusted_assay_tissue_new,
      #   get_adjusted_matrix(effect_removed, brms_fit_adjusted_tissue_new_estimate),
      #   packages = c( "brms", "glue", "dplyr", "purrr", "rstan", "magrittr", "stringr", "tidySummarizedExperiment") 
      # )
      
    )
    
    
  }, ask = FALSE, script = "~/labHead/cellNexus_article/ethnicity_prediction/_targets_fileted_by_Ning.R")
  
  
  tar_make(
    script = "~/labHead/cellNexus_article/ethnicity_prediction/_targets_fileted_by_Ning.R", 
    store = "~/labHead/cellNexus_article/ethnicity_prediction/_targets_fileted_by_Ning", 
    reporter = "verbose_positives")
})

tar_workspace(
  effect_removed_145922f1fc9c6b0f, 
  script = "~/labHead/cellNexus_article/ethnicity_prediction/_targets_fileted_by_Ning.R", 
  store = "~/labHead/cellNexus_article/ethnicity_prediction/_targets_fileted_by_Ning"
  )


tar_poll(
  store = "~/labHead/cellNexus_article/ethnicity_prediction/_targets_fileted_by_Ning"
)

tar_progress_summary(
  store = "~/labHead/cellNexus_article/ethnicity_prediction/_targets_fileted_by_Ning"
)

tar_meta(store = "~/labHead/cellNexus_article/ethnicity_prediction/_targets_fileted_by_Ning") |> 
  filter(name |> str_detect("effect_removed_")) |> 
  filter(bytes>0) |> 
  arrange(desc(time)) |> 
  slice(1) |> 
  pull(name) 

tar_read_raw("effect_removed_b9d568e7b8d9603b", store = "~/labHead/cellNexus_article/ethnicity_prediction/_targets_fileted_by_Ning")

library(tidyverse)
library(purrr)
library(targets)
library(SummarizedExperiment)

pseudobulk_sample = tar_read(
  pseudobulk_sample, 
  store = "~/labHead/cellNexus_article/ethnicity_prediction/_targets_fileted_by_Ning"
)


effect_removed = 
  
  # FOR INCOMPLETE PIPELINE
  tar_read(effect_removed,  store = "~/labHead/cellNexus_article/ethnicity_prediction/_targets_fileted_by_Ning") |> 
  
  mutate(max_nrow = map_int(brms_fit_adjusted_ethnicity, nrow) |> max()) |> 
  filter(map_int(brms_fit_adjusted_ethnicity, nrow) == max_nrow ) |> 
  
  
  mutate(brms_fit_adjusted = map(brms_fit_adjusted_ethnicity, ~ .x |> 
                                   select(adjusted___Estimate) |> #, adjusted___Q2.5, adjusted___Q97.5) |> 
                                   mutate(sample_id = colnames(pseudobulk_sample)), 
                                 .progress = TRUE
  )) |>
  select(.feature, brms_fit_adjusted) |> 
  unnest(brms_fit_adjusted) |>
  
  pivot_wider(names_from = sample_id, values_from = adjusted___Estimate) |> 
  tidybulk:::as_matrix(rownames = ".feature") |> 
  as("sparseMatrix")  |> 
  Matrix::Matrix(sparse = T)

effect_removed |> saveRDS("effect_removed_fileted_by_Ning.rds")
effect_removed = readRDS("effect_removed_fileted_by_Ning.rds")


pseudobulk_sample = pseudobulk_sample[rownames(effect_removed),, drop=FALSE ] 
assay(pseudobulk_sample, "counts_adjusted_ethnicity") = effect_removed
pseudobulk_sample = pseudobulk_sample[((assay(pseudobulk_sample, "gene_presence") > 0) |> rowSums() > (ncol(pseudobulk_sample) * 0.8)),,drop=FALSE ]

library(tidyomics)
library(magrittr)
summaries = 
    tar_read(
    summary,
    store =  "~/labHead/cellNexus_article/ethnicity_prediction/_targets_fileted_by_Ning"
  ) |>
  
  mutate(summary = map(summary, ~ .x %$% hypothesis |> as_tibble())) |>
  select(-summary_tissue) |> 
  unnest(summary) |>
  filter(Rhat_ethnicity |> dplyr::between(0.90, 1.1)) |> 
  select(.feature, Hypothesis, Estimate, Est.Error, CI.Lower, CI.Upper, Star) |> 
  pivot_wider(names_from = Hypothesis, values_from = -c(.feature, Hypothesis)) |> 
  # filter(Star == "*") |>
  filter(.feature %in% rownames(pseudobulk_sample))
  # |> 
  # mutate(closest_to_zero = pmin(abs(CI.Lower), abs(CI.Upper))) |>
  # add_count(.feature) |> 
  # filter(n < 5) |> 
  # with_groups(Hypothesis, ~ .x |> arrange(desc(closest_to_zero)) |> dplyr::slice(1:50))

pseudobulk_sample |> HDF5Array::saveHDF5SummarizedExperiment("pseudobulk_sample_fileted_by_Ning")
pseudobulk_sample = HDF5Array::loadHDF5SummarizedExperiment("pseudobulk_sample_fileted_by_Ning")

pseudobulk_sample_for_PCA = pseudobulk_sample


pseudobulk_sample_for_PCA = 
  pseudobulk_sample_for_PCA |> 
  select(-contains("PC"), -contains("tSNE"), -contains("UMAP")) |> 
  _[summaries |> pull(.feature) |> unique(), , drop=FALSE] 
  
 rd = rowData(pseudobulk_sample_for_PCA) |> as_tibble(rownames = ".feature") |> left_join(summaries) |> as("DFrame")
  rownames(rd) = rd$.feature
  rd = rd[,-1]
 rowData(pseudobulk_sample_for_PCA) = rd
 
 pseudobulk_sample_for_PCA = 
   pseudobulk_sample_for_PCA |> 
   filter(!ethnicity_groups %in% c("Other/Unknown", "Native American & Pacific Islander")) |> 
  mutate(counts_adjusted_ethnicity = pmax(counts_adjusted_ethnicity, 0)) |> 
  tidybulk::reduce_dimensions(method = "PCA", .abundance = counts_adjusted_ethnicity, .dims = 20 ) |> 
  tidybulk::reduce_dimensions(method = "tSNE", .abundance = counts_adjusted_ethnicity, initial_dims = 10, .dims = 2)  |> 
  tidybulk::reduce_dimensions(method = "UMAP", .abundance = counts_adjusted_ethnicity, pca = 10, .dims = 2, calculate_for_pca_dimensions = NULL)  

# pseudobulk_sample_for_PCA  = pseudobulk_sample_for_PCA |> filter(PC1 < 60)
# pseudobulk_sample_for_PCA  = pseudobulk_sample_for_PCA |> filter(PC4 > -20)

pseudobulk_sample_for_PCA |> HDF5Array::saveHDF5SummarizedExperiment("pseudobulk_sample_for_PCA_fileted_by_Ning")

pseudobulk_sample_for_PCA |> 
  as("SingleCellExperiment") |> 
  # anndataR::write_h5ad("~/labHead/cellNexus_article/ethnicity_prediction/pseudobulk_sample_for_PCA_adjusted_ethnicity.h5ad")
  zellkonverter::writeH5AD(file = "~/labHead/cellNexus_article/ethnicity_prediction/pseudobulk_sample_for_PCA_adjusted_ethnicity_fileted_by_Ning.h5ad", compression = "gzip")

system("~/bin/rclone copy ~/labHead/cellNexus_article/ethnicity_prediction/pseudobulk_sample_for_PCA_fileted_by_Ning box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/removal_unwanted_effects/pseudobulk_sample_for_PCA_fileted_by_Ning")
