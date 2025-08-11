
# Load required libraries, installing if not present
required_pkgs <- c(
  "GenomicRanges",
  "org.Hs.eg.db",
  "EnsDb.Hsapiens.v86",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",  # Pre-built gene coordinate table
  "GenomicFeatures",
  "AnnotationDbi",
  "tidymodels",
  "randomForest",
  "vip",
  "dplyr",
  "tibble",
  "airway",
  "zellkonverter",
  "scuttle",
  "edgeR",
  "parallel",     # For parallel processing
  "doParallel",   # For parallel backend
  "foreach",      # For parallel loops
  "future",       # Alternative parallel backend
  "future.apply",  # Future-based parallel apply
  "tidySummarizedExperiment",
  "tidyprint"
)

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing package:", pkg, "\n")
    if (pkg %in% c(
      "EnsDb.Hsapiens.v86",
      "org.Hs.eg.db",
      "GenomicRanges",
      "TxDb.Hsapiens.UCSC.hg38.knownGene",
      "GenomicFeatures",
      "AnnotationDbi",
      "airway",
      "scuttle",
      "edgeR"
    )) {
      # Bioconductor packages
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cran.r-project.org/")
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      # CRAN packages
      install.packages(pkg, repos = "https://cran.r-project.org/")
    }
  }
  library(pkg, character.only = TRUE)
}

# Set up parallel processing for tidymodels
cat("Setting up parallel processing...\n")
n_cores <- parallel::detectCores()
cat("Detected", n_cores, "CPU cores\n")

# Use 75% of available cores to avoid overwhelming the system
n_workers <- max(1, floor(n_cores * 0.75))
cat("Using", n_workers, "workers for parallel processing\n")

# Primary: Use future backend (more reliable for tidymodels)
library(future)
library(future.apply)
plan(multisession, workers = n_workers)
cat("Using future backend with", n_workers, "workers\n")

# Backup: Also set up doParallel for compatibility
library(doParallel)
cl <- makePSOCKcluster(n_workers)
registerDoParallel(cl)

# Verify parallel setup
cat("Future plan:", class(plan())[1], "\n")
cat("doParallel backend:", getDoParName(), "\n")
cat("doParallel workers:", getDoParWorkers(), "\n")

# Test parallel processing
cat("Testing parallel processing...\n")
test_result <- future_lapply(1:4, function(i) {
  Sys.sleep(1)  # Simulate work
  paste("Worker", i, "completed")
})
cat("Parallel test result:", paste(unlist(test_result), collapse = ", "), "\n")

# Alternative: Use future backend (uncomment if doParallel doesn't work)
# library(future)
# plan(multisession, workers = n_workers)
# cat("Using future backend with", n_workers, "workers\n")


#' Add chromosome coordinates using pre-built annotation tables (FASTER)
#'
#' This function uses pre-built Bioconductor annotation tables for much faster
#' coordinate lookup compared to querying EnsDb databases.
#'
#' @param se A SummarizedExperiment object with gene symbols or Ensembl IDs as rownames
#' @param genome_version Character, genome version to use ("hg38" or "hg19", default: "hg38")
#' @param verbose Logical, whether to print summary information (default: TRUE)
#'
#' @return A SummarizedExperiment object with genomic coordinates added to rowRanges
#'
#' @examples
#' se <- add_chromosome_coordinates_fast(se)
#' se <- add_chromosome_coordinates_fast(se, genome_version = "hg19")
add_chromosome_coordinates_fast <- function(se, genome_version = "hg38", verbose = TRUE) {
  
  # Validate input
  if (!is(se, "SummarizedExperiment")) {
    stop("Input must be a SummarizedExperiment object")
  }
  
  if (is.null(rownames(se))) {
    stop("SummarizedExperiment object must have rownames (gene symbols or Ensembl IDs)")
  }
  
  # Load appropriate TxDb package
  if (genome_version == "hg38") {
    if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
      stop("TxDb.Hsapiens.UCSC.hg38.knownGene package is required for hg38")
    }
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  } else if (genome_version == "hg19") {
    if (!requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene", quietly = TRUE)) {
      stop("TxDb.Hsapiens.UCSC.hg19.knownGene package is required for hg19")
    }
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  } else {
    stop("genome_version must be 'hg38' or 'hg19'")
  }
  
  # Get gene identifiers from SE object
  gene_ids <- rownames(se)
  
  if (verbose) {
    cat("Processing", length(gene_ids), "genes using pre-built", genome_version, "annotations...\n")
    cat("Sample gene IDs:", paste(head(gene_ids, 5), collapse = ", "), "\n")
  }
  
  # Check if we have Ensembl IDs or gene symbols
  is_ensembl <- any(grepl("^ENSG", gene_ids))
  
  # Extract gene coordinates from TxDb
  genes_gr <- genes(txdb)
  
  if (verbose) {
    cat("Loaded", length(genes_gr), "genes from", genome_version, "database\n")
  }
  
  if (is_ensembl) {
    # Clean Ensembl IDs (remove version numbers and suffixes)
    clean_ensembl_ids <- sub("_.*", "", sub("\\..*", "", gene_ids))
    
    if (verbose) {
      cat("Detected Ensembl gene IDs\n")
      cat("Sample cleaned IDs:", paste(head(clean_ensembl_ids, 5), collapse = ", "), "\n")
    }
    
    # Convert Entrez IDs to Ensembl IDs for matching
    entrez_to_ensembl <- mapIds(org.Hs.eg.db, keys = names(genes_gr), 
                                column = "ENSEMBL", keytype = "ENTREZID", 
                                multiVals = "first")
    
    # Add Ensembl IDs to the genes GRanges
    genes_gr$ensembl_id <- entrez_to_ensembl
    
      # Match user's Ensembl IDs with database
  if (verbose) {
    cat("Matching gene IDs to coordinates...\n")
  }
  matched_indices <- match(clean_ensembl_ids, genes_gr$ensembl_id)
    
  } else {
    # Gene symbols - convert Entrez IDs to symbols
    if (verbose) {
      cat("Detected gene symbols\n")
    }
    
    entrez_to_symbol <- mapIds(org.Hs.eg.db, keys = names(genes_gr), 
                               column = "SYMBOL", keytype = "ENTREZID", 
                               multiVals = "first")
    
    # Add gene symbols to the genes GRanges
    genes_gr$gene_symbol <- entrez_to_symbol
    
      # Match user's gene symbols with database
  if (verbose) {
    cat("Matching gene symbols to coordinates...\n")
  }
  matched_indices <- match(gene_ids, genes_gr$gene_symbol)
  }
  
  # Count matches
  genes_with_coords <- sum(!is.na(matched_indices))
  genes_without_coords <- sum(is.na(matched_indices))
  
  if (verbose) {
    cat("Matched", genes_with_coords, "out of", length(gene_ids), "genes\n")
    if (genes_without_coords > 0) {
      cat("Warning:", genes_without_coords, "genes could not be matched to coordinates\n")
    }
  }
  
  # Create GRanges object for all genes
  if (genes_with_coords > 0) {
    # Create initial GRanges with placeholder values
    feature_ranges <- GRanges(
      seqnames = rep("unknown", length(gene_ids)),
      ranges = IRanges(start = rep(1, length(gene_ids)), 
                       end = rep(1, length(gene_ids))),
      strand = rep("*", length(gene_ids)),
      gene_id = rep(NA_character_, length(gene_ids)),
      gene_name = rep(NA_character_, length(gene_ids))
    )
    
    # Set names to original gene IDs
    names(feature_ranges) <- gene_ids
    
    # Fill in coordinates for matched genes
    valid_matches <- !is.na(matched_indices)
    
    if (any(valid_matches)) {
      matched_coords <- genes_gr[matched_indices[valid_matches]]
      matched_gene_indices <- which(valid_matches)
      
      if (verbose) {
        cat("Assigning coordinates to", length(matched_gene_indices), "genes...\n")
      }
      
      for (i in seq_along(matched_gene_indices)) {
        gene_idx <- matched_gene_indices[i]
        coord_data <- matched_coords[i]
        
        feature_ranges[gene_idx] <- GRanges(
          seqnames = seqnames(coord_data),
          ranges = ranges(coord_data),
          strand = strand(coord_data),
          gene_id = names(coord_data),  # Entrez ID
          gene_name = if(is_ensembl) coord_data$ensembl_id else coord_data$gene_symbol
        )
      }
    }
    
    # Add the genomic coordinates to rowRanges
    rowRanges(se) <- feature_ranges
    
    # Display summary
    if (verbose) {
      cat("Added chromosome coordinates for", genes_with_coords, "out of", 
          length(gene_ids), "genes\n")
      
      valid_seqnames <- seqnames(rowRanges(se))[!is.na(seqnames(rowRanges(se)))]
      if (length(valid_seqnames) > 0) {
        cat("Chromosomes present:", paste(unique(as.character(valid_seqnames)), collapse = ", "), "\n")
      }
    }
    
  } else {
    if (verbose) {
      cat("Warning: No genes could be matched to coordinates. Check gene identifiers.\n")
    }
    return(se)  # Return original SE object unchanged
  }
  
  return(se)
}

# Read SE file
cat("Reading H5AD file...\n")
se <- zellkonverter::readH5AD("~/Downloads/pseudobulk_se.h5ad", reader = "R", use_hdf5 = TRUE )

# Add chromosome coordinates using the FAST function with pre-built tables
cat("Adding chromosome coordinates...\n")
se <- add_chromosome_coordinates_fast(se, genome_version = "hg38")

# Alternative: Use the original EnsDb function (slower but more flexible)
# se <- add_chromosome_coordinates(se)

# This operation takes to much memory we need to do this in chunks
# use map and print progress based on the chunk size
# Split column wise into chunks of 1000
cols <- colnames(se)
col_chunks <- split(cols, ceiling(seq_along(cols)/1000))

cat("Processing", length(col_chunks), "chunks for pseudobulk aggregation...\n")

# Map over the chunks with progress tracking
res <- map(col_chunks, function(chunk) {
  se_chunk <- se[,chunk]
  se_chunk <- scuttle::aggregateAcrossCells(  
    se_chunk,
    ids = colData(se_chunk)$sample_id,
    store.number = TRUE,
    BPPARAM = BiocParallel::MulticoreParam(workers = parallel::detectCores())
  )
  se_chunk
})

cat("Combining chunks...\n")
se <- do.call(cbind, res)

se |> saveRDS("~/Documents/GitHub/article_cellNexus/pseudobulk_se_X_Y_subset_sample_id.rds")

se = readRDS("~/Documents/GitHub/article_cellNexus/pseudobulk_se_X_Y_subset_sample_id.rds")

# rm(res)

# Filter out samples with 0 total counts
cat("Filtering samples with zero counts...\n")
se <- se[, colSums(assay(se, "counts")) > 0]

# Optional convenience assay: logCPM
cat("Calculating logCPM values...\n")
assay(se, "logCPM") <- edgeR::cpm(assay(se, "counts"), log = TRUE, prior.count = 1)


# Calculate rank matrix using singscore from bioconductor
library(singscore)
cat("Calculating gene ranks using singscore...\n")
assay(se, "singscore") <- assay(se, "logCPM") |> singscore::rankGenes()


# Filter just the sex chromosomes
cat("Filtering to sex chromosomes only...\n")
se <- se[seqnames(rowRanges(se)) %in% c("chrX", "chrY"), ]
cat("Retained", nrow(se), "sex chromosome genes\n")


# Load tidymodels libraries
library(tidymodels)
library(randomForest)
library(vip)
library(dplyr)

# Set memory limits for parallel processing
options(future.globals.maxSize = 2 * 1024^3)  # 2GB limit
options(future.plan = "multisession", workers = 2)  # Reduce workers to 2

# Build a random forest that is able to predict the sex of the samples
# This should be robust to some simple covariates such as batch and batch effect correction and replicate
# This should be robust to label swap in the training data

#' Build a robust sex prediction model using random forest
#'
#' This function creates a tidymodels workflow for sex prediction using
#' gene expression data from sex chromosomes. The model is designed to be
#' robust to batch effects and potential label swaps in training data.
#'
#' @param se A SummarizedExperiment object with sex chromosome genes
#' @param sex_column Character, name of column containing sex labels (default: "sex")
#' @param batch_columns Character vector, names of batch/covariate columns to include
#' @param train_prop Numeric, proportion of data to use for training (default: 0.8)
#' @param trees Numeric, number of trees in random forest (default: 1000)
#' @param min_n Numeric, minimum observations in terminal nodes (default: 5)
#' @param mtry Numeric, number of predictors to sample at each split (default: NULL, will tune)
#' @param seed Numeric, random seed for reproducibility (default: 123)
#'
#' @return A list containing the trained workflow, performance metrics, and data splits
build_sex_prediction_model <- function(se,
                                       my_assay = "logCPM",
                                       sex_column = "sex",
                                       batch_columns = c("batch", "replicate"),
                                       train_prop = 0.8,
                                       trees = 1000,
                                       min_n = 5,
                                       mtry = NULL,
                                       seed = 123) {
  
  set.seed(seed)
  
  # Extract expression data and metadata
  expr_data <- t(assay(se, my_assay)) # transpose so samples are rows
  metadata <- colData(se) %>% as.data.frame()
  
  # Combine expression and metadata
  model_data <- cbind(metadata, expr_data) %>%
    as_tibble()
  
  # Check if sex column exists
  if (!sex_column %in% colnames(model_data)) {
    stop(paste("Sex column", sex_column, "not found in metadata"))
  }
  
  # Convert sex to factor for classification and remove unknown/NA values
  model_data[[sex_column]] <- as.factor(model_data[[sex_column]])
  
  # Remove samples with unknown, NA, or ambiguous sex labels
  valid_sex_labels <- c("female", "male", "Female", "Male", "F", "M", "f", "m")
  model_data <- model_data %>%
    filter(!is.na(!!sym(sex_column)) & 
           as.character(!!sym(sex_column)) %in% valid_sex_labels)
  
  # Standardize sex labels to lowercase
  model_data[[sex_column]] <- tolower(as.character(model_data[[sex_column]]))
  model_data[[sex_column]] <- case_when(
    model_data[[sex_column]] %in% c("female", "f") ~ "female",
    model_data[[sex_column]] %in% c("male", "m") ~ "male",
    TRUE ~ NA_character_
  )
  model_data <- model_data %>% filter(!is.na(!!sym(sex_column)))
  model_data[[sex_column]] <- as.factor(model_data[[sex_column]])
  
  if (nrow(model_data) == 0) {
    stop("No samples with valid sex labels found. Check your sex column values.")
  }
  
  # Select features: sex chromosomes genes + batch covariates
  feature_cols <- c(sex_column, 
                    intersect(batch_columns, colnames(model_data)),
                    colnames(expr_data))
  
  model_data <- model_data %>%
    select(all_of(feature_cols)) %>%
    # Remove any columns with all NA or zero variance (handle factors properly)
    select(where(~ {
      !all(is.na(.x)) && 
      if(is.numeric(.x)) {
        var(.x, na.rm = TRUE) > 0
      } else if(is.factor(.x) || is.character(.x)) {
        length(unique(.x[!is.na(.x)])) > 1  # More than one unique value
      } else {
        TRUE
      }
    }))
  
  cat("Using", ncol(model_data) - 1, "features for prediction\n")
  cat("Target variable levels:", paste(levels(model_data[[sex_column]]), collapse = ", "), "\n")
  
  # Create data split
  data_split <- initial_split(model_data, prop = train_prop, strata = !!sym(sex_column))
  train_data <- training(data_split)
  test_data <- testing(data_split)
  
  cat("Training samples:", nrow(train_data), "| Test samples:", nrow(test_data), "\n")
  
  # Check class distribution
  class_counts <- table(train_data[[sex_column]])
  cat("Class distribution in training:", paste(names(class_counts), "=", class_counts, collapse = ", "), "\n")
  
  # Check for minimum viable dataset size and class balance
  min_class_size <- min(class_counts)
  if (nrow(train_data) < 10) {
    warning("Very small training dataset (", nrow(train_data), " samples). Results may be unreliable.")
  }
  if (min_class_size < 3) {
    warning("Very few samples in smallest class (", min_class_size, "). Consider collecting more data.")
  }
  
  # Create recipe for preprocessing
  # This helps with batch effect robustness
  recipe_spec <- recipe(as.formula(paste(sex_column, "~ .")), data = train_data) %>%
    # Handle any remaining missing values first
    step_impute_median(all_numeric_predictors()) %>%
    # Remove zero variance predictors BEFORE normalization
    step_zv(all_predictors()) %>%
    # Center and scale all numeric predictors (gene expression)
    step_normalize(all_numeric_predictors()) %>%
    # Remove highly correlated predictors to reduce overfitting
    step_corr(all_numeric_predictors(), threshold = 0.95) %>%
    # Optional: PCA for dimensionality reduction (adjust for small sample size)
    step_pca(all_numeric_predictors(), 
             num_comp = min(10, floor(nrow(train_data) * 0.5), ncol(expr_data) - 1),
             threshold = 0.95)
  
  # Define random forest model with tuning
  if (is.null(mtry)) {
    rf_spec <- rand_forest(
      trees = trees,
      min_n = min_n,
      mtry = tune()
    ) %>%
      set_engine("randomForest", 
                 importance = TRUE,
                 do.trace = FALSE,
                 keep.forest = TRUE) %>%
      set_mode("classification")
  } else {
    rf_spec <- rand_forest(
      trees = trees,
      min_n = min_n,
      mtry = mtry
    ) %>%
      set_engine("randomForest", 
                 importance = TRUE,
                 do.trace = FALSE,
                 keep.forest = TRUE) %>%
      set_mode("classification")
  }
  
  # Create workflow
  rf_workflow <- workflow() %>%
    add_recipe(recipe_spec) %>%
    add_model(rf_spec)
  
  # Cross-validation for robustness assessment (adjust folds for small datasets)
  n_folds <- min(5, min_class_size, floor(nrow(train_data) / 2))
  if (n_folds < 2) {
    warning("Dataset too small for cross-validation. Using simple validation.")
    n_folds <- 2
  }
  
  cv_folds <- vfold_cv(train_data, v = n_folds, strata = !!sym(sex_column))
  cat("Using", n_folds, "fold cross-validation\n")
  
  # Tune hyperparameters if mtry is not specified
  if (is.null(mtry)) {
    cat("Tuning hyperparameters...\n")
    
    # Create tuning grid (adjust for small sample size)
    max_mtry <- min(10, ncol(expr_data), floor(sqrt(ncol(expr_data))))
    rf_grid <- grid_regular(
      mtry(range = c(1, max_mtry)),
      levels = min(3, max_mtry)  # Fewer levels for small datasets
    )
    
    cat("Testing", nrow(rf_grid), "hyperparameter combinations with", n_folds, "fold CV...\n")
    cat("Using future plan:", class(plan())[1], "with", nbrOfWorkers(), "workers\n")
    cat("Using doParallel with", getDoParWorkers(), "workers\n")
    
    # Tune the model with parallel processing
    rf_tune_results <- tune_grid(
      rf_workflow,
      resamples = cv_folds,
      grid = rf_grid,
      metrics = metric_set(accuracy, roc_auc, sensitivity, specificity),
      control = control_grid(
        verbose = FALSE,
        allow_par = TRUE,
        parallel_over = "everything"
      )
    )
    
    # Select best parameters
    best_params <- select_best(rf_tune_results, metric = "roc_auc")
    rf_workflow <- finalize_workflow(rf_workflow, best_params)
    
    cat("Best mtry:", best_params$mtry, "\n")
  }
  
  # Fit final model
  cat("Training final model with", trees, "trees...\n")
  
  # Fit the model
  final_fit <- fit(rf_workflow, train_data)
  
  cat("Model training completed!\n")
  
  # Evaluate on test set
  cat("Evaluating model on test set...\n")
  test_predictions <- predict(final_fit, test_data, type = "prob") %>%
    bind_cols(predict(final_fit, test_data)) %>%
    bind_cols(test_data %>% select(!!sym(sex_column)))
  
  # Calculate performance metrics
  test_metrics <- test_predictions %>%
    metrics(truth = !!sym(sex_column), estimate = .pred_class) %>%
    bind_rows(
      test_predictions %>%
        roc_auc(truth = !!sym(sex_column), paste0(".pred_", levels(model_data[[sex_column]])[1]))
    )
  
  # Cross-validation performance for robustness assessment
  cat("Performing cross-validation for robustness assessment...\n")
  cv_metrics <- fit_resamples(
    rf_workflow,
    resamples = cv_folds,
    metrics = metric_set(accuracy, roc_auc, sensitivity, specificity),
    control = control_resamples(
      verbose = FALSE,
      allow_par = FALSE,  # Disable parallel processing to avoid memory issues
      parallel_over = "everything"
    )
  ) %>%
    collect_metrics()
  
  # Print results
  cat("\n=== Model Performance ===\n")
  cat("Cross-validation results (mean ± sd):\n")
  print(cv_metrics)
  
  cat("\nTest set performance:\n")
  print(test_metrics)
  
  # Print final summary
  cat("\n", strrep("=", 60), "\n")
  cat("MODEL TRAINING COMPLETED\n")
  cat(strrep("=", 60), "\n")
  cat("Training samples:", nrow(train_data), "\n")
  cat("Test samples:", nrow(test_data), "\n")
  cat("Features used:", ncol(model_data) - 1, "\n")
  cat("Cross-validation folds:", n_folds, "\n")
  
  # Return comprehensive results
  list(
    workflow = final_fit,
    recipe = recipe_spec,
    test_metrics = test_metrics,
    cv_metrics = cv_metrics,
    test_predictions = test_predictions,
    data_split = data_split,
    feature_importance = NULL  # Will be added if model is fitted
  )
}

# Function to create tissue-specific validation sets
create_tissue_validation_split <- function(se, 
                                          sex_column = "sex",
                                          tissue_column = "tissue_groups",
                                          male_tissues = c("prostate"),
                                          female_tissues = c("breast", "female reproductive system", "ovary"),
                                          train_prop = 0.8,
                                          seed = 123) {
  
  set.seed(seed)
  
  # Get metadata
  metadata <- colData(se) %>% as.data.frame()
  
  # Check if required columns exist
  if (!sex_column %in% colnames(metadata)) {
    stop(paste("Sex column", sex_column, "not found in metadata"))
  }
  if (!tissue_column %in% colnames(metadata)) {
    stop(paste("Tissue column", tissue_column, "not found in metadata"))
  }
  
  # Filter for valid sex labels (exclude unknown)
  valid_sex_labels <- c("female", "male", "Female", "Male", "F", "M", "f", "m")
  metadata <- metadata %>%
    filter(!is.na(!!sym(sex_column)) & 
           as.character(!!sym(sex_column)) %in% valid_sex_labels &
           as.character(!!sym(sex_column)) != "unknown")
  
  # Filter out embryonic/fetal samples (age_days < 365) and samples with fetal/embryonic titles
  metadata <- metadata %>%
    filter(
      # Keep samples with age_days >= 365
      (!is.na(age_days) & age_days >= 365) |
      # Keep samples with missing age_days but exclude fetal/embryonic titles
      (is.na(age_days) & 
       !grepl("fetal|embryo|placenta|cord blood", tolower(title), ignore.case = TRUE))
    )
  
  cat("Samples after filtering embryonic/fetal (age_days >= 365):", nrow(metadata), "\n")
  cat("Samples with missing age_days (kept):", sum(is.na(metadata$age_days)), "\n")
  
  # Standardize sex labels
  metadata[[sex_column]] <- tolower(as.character(metadata[[sex_column]]))
  metadata[[sex_column]] <- case_when(
    metadata[[sex_column]] %in% c("female", "f") ~ "female",
    metadata[[sex_column]] %in% c("male", "m") ~ "male",
    TRUE ~ NA_character_
  )
  metadata <- metadata %>% filter(!is.na(!!sym(sex_column)))
  
  # Create tissue categories
  metadata$tissue_category <- case_when(
    metadata[[tissue_column]] %in% male_tissues ~ "male_specific",
    metadata[[tissue_column]] %in% female_tissues ~ "female_specific",
    TRUE ~ "other"
  )
  
  # Print tissue distribution
  cat("=== Tissue Distribution ===\n")
  tissue_dist <- table(metadata[[tissue_column]], metadata[[sex_column]])
  print(tissue_dist)
  
  cat("\n=== Tissue Categories ===\n")
  cat_dist <- table(metadata$tissue_category, metadata[[sex_column]])
  print(cat_dist)
  
  # Create validation sets
  # Use male-specific tissues for male validation
  male_validation <- metadata %>%
    filter(!!sym(sex_column) == "male" & tissue_category == "male_specific")
  
  # Use female-specific tissues for female validation  
  female_validation <- metadata %>%
    filter(!!sym(sex_column) == "female" & tissue_category == "female_specific")
  
  # Use other tissues for training
  training_samples <- metadata %>%
    filter(tissue_category == "other")
  
  # Ensure we have enough samples in each set
  if (nrow(male_validation) == 0) {
    warning("No male-specific tissue samples found. Using random male samples for validation.")
    male_samples <- metadata %>% filter(!!sym(sex_column) == "male")
    if (nrow(male_samples) > 0) {
      male_validation <- male_samples %>% slice_sample(n = min(50, nrow(male_samples)))
      training_samples <- metadata %>% filter(!sample_id %in% male_validation$sample_id)
    }
  }
  
  if (nrow(female_validation) == 0) {
    warning("No female-specific tissue samples found. Using random female samples for validation.")
    female_samples <- metadata %>% filter(!!sym(sex_column) == "female")
    if (nrow(female_samples) > 0) {
      female_validation <- female_samples %>% slice_sample(n = min(50, nrow(female_samples)))
      training_samples <- metadata %>% filter(!sample_id %in% female_validation$sample_id)
    }
  }
  
  # Combine validation sets
  validation_samples <- bind_rows(male_validation, female_validation)
  
  cat("\n=== Sample Split ===\n")
  cat("Training samples:", nrow(training_samples), "\n")
  cat("Validation samples:", nrow(validation_samples), "\n")
  cat("  - Male-specific tissue validation:", nrow(male_validation), "\n")
  cat("  - Female-specific tissue validation:", nrow(female_validation), "\n")
  
  # Create sample IDs for splitting
  train_ids <- training_samples$sample_id
  validation_ids <- validation_samples$sample_id
  
  return(list(
    train_ids = train_ids,
    validation_ids = validation_ids,
    training_metadata = training_samples,
    validation_metadata = validation_samples
  ))
}

# Train the sex prediction model
# Note: Adjust the sex_column and batch_columns parameters based on your actual metadata
cat("\n", strrep("=", 60), "\n")
cat("BUILDING SEX PREDICTION MODEL\n")
cat(strrep("=", 60), "\n")

# Create tissue-specific validation split
cat("Creating tissue-specific validation split...\n")
cat("Note: Only prostate tissue available for male reproductive validation\n")
tissue_split <- create_tissue_validation_split(
  se = se,
  sex_column = "sex",
  tissue_column = "tissue_groups",
  male_tissues = c("prostate"),  # Only male reproductive tissue available
  female_tissues = c("breast", "female reproductive system", "ovary"),
  seed = 123
)

# Filter SE objects for training and validation
se_train <- se[, colnames(se) %in% tissue_split$train_ids]
se_validation <- se[, colnames(se) %in% tissue_split$validation_ids]

# Filter validation set to exclude unknown sex labels
colData(se_validation)$sex <- colData(se_validation)$sex |> droplevels()  

cat("Training SE dimensions:", dim(se_train), "\n")
cat("Validation SE dimensions:", dim(se_validation), "\n")

# If training set is too large, sample a subset to avoid memory issues
if (ncol(se_train) > 10000) {
  cat("Training set is large, sampling 10000 samples to avoid memory issues...\n")
  set.seed(123)
  sample_cols <- sample(colnames(se_train), 10000)
  se_train <- se_train[, sample_cols]
  cat("Sampled training SE dimensions:", dim(se_train), "\n")
}

# Clean up memory
gc()

# Save session state
save.image("~/Documents/GitHub/article_cellNexus/sex_prediction_session_state.RData")

# Build model on training data
cat("Building sex prediction model on training data...\n")
sex_model_results <- build_sex_prediction_model(
  se = se_train,
  my_assay = "singscore",
  sex_column = "sex",
  batch_columns = c("assay", "tissue_groups", "age_days"),  # Added age_days
  train_prop = 0.8,
  trees = 500,  # Reduced from 1000 to save memory
  seed = 123
)


#=== Model Performance ===
#Cross-validation results (mean ± sd):
# A tibble: 4 × 6
#  .metric     .estimator  mean     n  std_err .config             
#  <chr>       <chr>      <dbl> <int>    <dbl> <chr>               
#1 accuracy    binary     0.985     5 0.000523 Preprocessor1_Model1
#2 roc_auc     binary     0.994     5 0.000810 Preprocessor1_Model1
#3 sensitivity binary     0.984     5 0.00258  Preprocessor1_Model1
#4 specificity binary     0.985     5 0.00166  Preprocessor1_Model1

#Test set performance:
# A tibble: 3 × 3
#  .metric  .estimator .estimate
#  <chr>    <chr>          <dbl>
#1 accuracy binary         0.991
#2 kap      binary         0.982
#3 roc_auc  binary         0.996

# ============================================================ 
#MODEL TRAINING COMPLETED
#============================================================ 
#Training samples: 7999 
#Test samples: 2001 
#Features used: 370 
#Cross-validation folds: 5

# save model
saveRDS(sex_model_results, "~/Documents/GitHub/article_cellNexus/sex_prediction_model.rds")


# Validate on tissue-specific validation set
cat("\n=== Tissue-Specific Validation ===\n")
cat("Validating on sex-specific tissues...\n")

# Get predictions on validation set
validation_predictions <- predict(sex_model_results$workflow, 
                                 new_data = t(assay(se_validation, "singscore")) %>%
                                   as.data.frame() %>%
                                   bind_cols(colData(se_validation) %>% as.data.frame()),
                                 type = "prob") %>%
  bind_cols(predict(sex_model_results$workflow, 
                   new_data = t(assay(se_validation, "singscore")) %>%
                     as.data.frame() %>%
                     bind_cols(colData(se_validation) %>% as.data.frame()))) %>%
  bind_cols(colData(se_validation) %>% as.data.frame() %>% select(sex, tissue_groups))


cat("Validation predictions after filtering unknown sex:", nrow(validation_predictions), "\n")
cat("Sex levels in validation predictions:", paste(levels(validation_predictions$sex), collapse = ", "), "\n")
  
  # Calculate validation metrics
  validation_metrics <- validation_predictions %>%
    metrics(truth = sex, estimate = .pred_class) %>%
    bind_rows(
      validation_predictions %>%
        roc_auc(truth = sex, .pred_female)
    )
  
  cat("Tissue-specific validation performance:\n")
  print(validation_metrics)
  
  # Show performance by tissue type
  cat("\nPerformance by tissue type:\n")
  tissue_performance <- validation_predictions %>%
    group_by(tissue_groups) %>%
    summarise(
      n_samples = n(),
      accuracy = sum(sex == .pred_class) / n(),
      .groups = "drop"
    ) %>%
    arrange(desc(n_samples))
  
    print(tissue_performance)
  
  # Extract the trained model for future predictions
trained_model <- sex_model_results$workflow

# Display feature importance (if available)
if ("randomForest" %in% class(extract_fit_engine(trained_model))) {
  cat("\n=== Top 10 Most Important Features ===\n")
  importance_scores <- extract_fit_engine(trained_model)$importance
  if (!is.null(importance_scores)) {
    top_features <- importance_scores[order(importance_scores[,"MeanDecreaseGini"], decreasing = TRUE)[1:10], ]
    print(top_features)
  }
}

# Plot loadings of top_features PCA
# Plot feature importance (if available)
cat("Top 10 most important features:\n")
print(top_features)

# Plot ROC curve on training data
library(ggplot2)
library(yardstick)

# Get training predictions by applying model to training data
# First, let's get the training data from the original SE object
train_data <- t(assay(se_train, "singscore")) %>%
  as.data.frame() %>%
  bind_cols(colData(se_train) %>% as.data.frame() %>% select(sex, assay, tissue_groups, age_days))

# Get predictions on training data
train_predictions <- predict(sex_model_results$workflow, 
                           new_data = train_data,
                           type = "prob") %>%
  bind_cols(sex = train_data$sex) %>%
  # Drop unused factor levels
  mutate(sex = droplevels(sex))

# Compute ROC curve data using training predictions
roc_data <- roc_curve(train_predictions, truth = sex, .pred_female)

# Plot ROC curve
ggplot(roc_data, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(color = "blue", size = 1.2) +
  geom_abline(linetype = "dashed", color = "gray") +
  labs(
    title = "ROC Curve on Training Data",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  theme_minimal()

# Save model for future use (optional)
# saveRDS(sex_model_results, "sex_prediction_model.rds")

cat("\n", strrep("=", 60), "\n")
cat("ANALYSIS COMPLETED SUCCESSFULLY\n")
cat(strrep("=", 60), "\n")
cat("Final model saved and ready for predictions\n")
# now select the sex == unknown and predict the sex
se_unknown <- se[, colData(se)$sex == "unknown" | is.na(colData(se)$sex)]

# Check if age_days is available in the unknown samples
cat("Unknown samples with age_days:", sum(!is.na(colData(se_unknown)$age_days)), "\n")
cat("Unknown samples without age_days:", sum(is.na(colData(se_unknown)$age_days)), "\n")

# Filter unknown samples to match training data (age_days >= 365 or missing age_days but not fetal/embryonic)
se_unknown_filtered <- se_unknown[, 
  (!is.na(colData(se_unknown)$age_days) & colData(se_unknown)$age_days >= 365) |
  (is.na(colData(se_unknown)$age_days) & 
   !grepl("fetal|embryo|placenta|cord blood", tolower(colData(se_unknown)$title), ignore.case = TRUE))
]

cat("Unknown samples after filtering embryonic/fetal (age_days >= 365):", ncol(se_unknown_filtered), "\n")
cat("Samples excluded (embryonic/fetal or missing age):", ncol(se_unknown) - ncol(se_unknown_filtered), "\n")

# Use the filtered dataset for predictions
se_unknown <- se_unknown_filtered

se_unknown_predictions <- predict(sex_model_results$workflow, 
                                  new_data = t(assay(se_unknown, "singscore")) %>%
                                    as.data.frame() %>%
                                    bind_cols(colData(se_unknown) %>% as.data.frame()),
                                  type = "prob")

# Get class predictions (not just probabilities)
se_unknown_class_predictions <- predict(sex_model_results$workflow, 
                                       new_data = t(assay(se_unknown, "singscore")) %>%
                                         as.data.frame() %>%
                                         bind_cols(colData(se_unknown) %>% as.data.frame()))

# add the predictions to the colData
colData(se_unknown)$sex_prediction <- se_unknown_class_predictions$.pred_class

# count tissue_groups and sex_prediction
colData(se_unknown) |>
 as_tibble() |> 
 dplyr::count(tissue, tissue_groups, sex_prediction, age_days) |> 
 print(n=99)

# save the tissue_prediction_counts
write.csv(tissue_prediction_counts, "~/Documents/GitHub/article_cellNexus/tissue_prediction_counts.csv", row.names = FALSE)

# save the se_unknown
saveRDS(se_unknown, "~/Documents/GitHub/article_cellNexus/se_unknown.rds")

