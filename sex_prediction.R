
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
  "BiocParallel",
  "parallel",     # For parallel processing
  "doParallel",   # For parallel backend
  "foreach",      # For parallel loops
  "future",       # Alternative parallel backend
  "future.apply",  # Future-based parallel apply
  "tidySummarizedExperiment",
  "tidyprint",
  "stringr",
  "purrr",
  "singscore",
  "ggplot2",
  "yardstick",
  "patchwork",
  "rlang"
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
#' @param se A `SummarizedExperiment` with gene symbols or Ensembl IDs as `rownames(se)`
#' @param genome_version Character, genome version to use ("hg38" or "hg19", default: "hg38")
#' @param verbose Logical, whether to print summary information (default: TRUE)
#'
#' @return A `SummarizedExperiment` with genomic coordinates added to `rowRanges(se)`
#'
#' @examples
#' \dontrun{
#' se <- add_chromosome_coordinates_fast(se)
#' se <- add_chromosome_coordinates_fast(se, genome_version = "hg19")
#' }
#'
#' @importFrom GenomicFeatures genes
#' @importFrom GenomicRanges GRanges seqnames ranges strand
#' @importFrom IRanges IRanges
#' @importFrom AnnotationDbi mapIds
#' @importFrom SummarizedExperiment rowRanges rowRanges<- 
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

#' Create tissue-specific validation sets for sex prediction
#'
#' This function splits a SummarizedExperiment object into tissue-specific validation sets
#' based on sex and tissue type. It filters out samples with ambiguous or missing sex labels,
#' removes embryonic/fetal samples (unless age_days is missing and not in fetal/placental tissues),
#' and categorizes tissues as male-specific, female-specific, or other.
#'
#' @param se A SummarizedExperiment object containing sample metadata.
#' @param sex_column The name of the column in colData(se) indicating sex (default: "sex").
#' @param tissue_column The name of the column in colData(se) indicating tissue type (default: "tissue_groups").
#' @param male_tissues Character vector of tissue names considered male-specific (default: c("prostate")).
#' @param female_tissues Character vector of tissue names considered female-specific (default: c("breast", "female reproductive system", "ovary")).
#' @param train_prop Proportion of data to use for training (default: 0.8).
#' @param seed Random seed for reproducibility (default: 123).
#'
#' @return A list containing data frames for male-specific and female-specific validation sets,
#'         as well as the filtered metadata with tissue categories.
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr filter case_when %>% arrange group_by summarise
#' @importFrom tibble as_tibble
#' @importFrom stats setNames
#' @importFrom methods as
#' @export
create_tissue_validation_split <- function(se, 
                                          sex_column = "sex",
                                          tissue_column = "tissue_groups",
                                          male_tissues = c("prostate"),
                                          female_tissues = c("breast", "female reproductive system", "ovary"),
                                          train_prop = 0.8,
                                          seed = 123) {
  set.seed(seed)

  metadata <- colData(se) %>% as.data.frame()

  if (!sex_column %in% colnames(metadata)) {
    stop(paste("Sex column", sex_column, "not found in metadata"))
  }
  if (!tissue_column %in% colnames(metadata)) {
    stop(paste("Tissue column", tissue_column, "not found in metadata"))
  }

  valid_sex_labels <- c("female", "male", "Female", "Male", "F", "M", "f", "m")
  metadata <- metadata %>%
    filter(!is.na(!!sym(sex_column)) & 
             as.character(!!sym(sex_column)) %in% valid_sex_labels &
             as.character(!!sym(sex_column)) != "unknown")

  metadata <- metadata %>%
    filter((!is.na(age_days) & age_days >= 365) |
             (is.na(age_days) & 
                !grepl("fetal|embryo|placenta|cord blood", tolower(title), ignore.case = TRUE)))

  cat("Samples after filtering embryonic/fetal (age_days >= 365):", nrow(metadata), "\n")
  cat("Samples with missing age_days (kept):", sum(is.na(metadata$age_days)), "\n")

  metadata[[sex_column]] <- tolower(as.character(metadata[[sex_column]]))
  metadata[[sex_column]] <- case_when(
    metadata[[sex_column]] %in% c("female", "f") ~ "female",
    metadata[[sex_column]] %in% c("male", "m") ~ "male",
    TRUE ~ NA_character_
  )
  metadata <- metadata %>% filter(!is.na(!!sym(sex_column)))

  metadata$tissue_category <- case_when(
    metadata[[tissue_column]] %in% male_tissues ~ "male_specific",
    metadata[[tissue_column]] %in% female_tissues ~ "female_specific",
    TRUE ~ "other"
  )

  cat("=== Tissue Distribution ===\n")
  tissue_dist <- table(metadata[[tissue_column]], metadata[[sex_column]])
  print(tissue_dist)

  cat("\n=== Tissue Categories ===\n")
  cat_dist <- table(metadata$tissue_category, metadata[[sex_column]])
  print(cat_dist)

  male_validation <- metadata %>%
    filter(!!sym(sex_column) == "male" & tissue_category == "male_specific")
  female_validation <- metadata %>%
    filter(!!sym(sex_column) == "female" & tissue_category == "female_specific")
  training_samples <- metadata %>%
    filter(tissue_category == "other")

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

  validation_samples <- bind_rows(male_validation, female_validation)

  cat("\n=== Sample Split ===\n")
  cat("Training samples:", nrow(training_samples), "\n")
  cat("Validation samples:", nrow(validation_samples), "\n")
  cat("  - Male-specific tissue validation:", nrow(male_validation), "\n")
  cat("  - Female-specific tissue validation:", nrow(female_validation), "\n")

  train_ids <- training_samples$sample_id
  validation_ids <- validation_samples$sample_id

  return(list(
    train_ids = train_ids,
    validation_ids = validation_ids,
    training_metadata = training_samples,
    validation_metadata = validation_samples
  ))
}

#' Preprocess a SummarizedExperiment for sex prediction and create splits
#'
#' Performs coordinate annotation, pseudobulk aggregation in chunks, logCPM
#' calculation, singscore ranking, sex-chromosome filtering, and produces
#' training, test, and validation `SummarizedExperiment` objects.
#'
#' @param se A `SummarizedExperiment` (cells/samples as columns)
#' @param genome_version Genome version for coordinate mapping (default: "hg38")
#' @param my_assay Name of assay to generate/use for modeling (default: "singscore")
#' @param chunk_size Integer chunk size for pseudobulk aggregation (default: 1000)
#' @param sex_column Column name in `colData(se)` with sex labels (default: "sex")
#' @param tissue_column Column name in `colData(se)` with tissue labels (default: "tissue_groups")
#' @param train_prop Proportion of labeled data used for training (default: 0.8)
#' @param seed Random seed (default: 123)
#'
#' @return A list with `se_train`, `se_test`, `se_validation`, and `se_processed`
#'
#' @examples
#' \dontrun{
#' splits <- preprocess_for_sex_prediction(se)
#' se_train <- splits$se_train
#' se_test  <- splits$se_test
#' se_validation <- splits$se_validation
#' }
#'
#' @importFrom scuttle aggregateAcrossCells
#' @importFrom BiocParallel MulticoreParam
#' @importFrom SummarizedExperiment colData assay assay<-
#' @importFrom edgeR cpm
#' @importFrom singscore rankGenes
#' @importFrom GenomicRanges rowRanges seqnames
#' @importFrom dplyr filter select mutate arrange bind_rows
#' @importFrom rsample initial_split training testing
#' @importFrom tibble as_tibble
#' @importFrom purrr map
preprocess_for_sex_prediction <- function(
  se,
  genome_version = "hg38",
  my_assay = "singscore",
  chunk_size = 1000,
  sex_column = "sex",
  tissue_column = "tissue_groups",
  train_prop = 0.8,
  seed = 123
) {
  set.seed(seed)

  # Add chromosome coordinates
  cat("Adding chromosome coordinates...\n")
  se <- add_chromosome_coordinates_fast(se, genome_version = genome_version)

  # Aggregate across cells in chunks to avoid high memory use
  cat("Preparing pseudobulk in chunks of", chunk_size, "columns...\n")
  cols <- colnames(se)
  col_chunks <- split(cols, ceiling(seq_along(cols)/chunk_size))
  cat("Processing", length(col_chunks), "chunks for pseudobulk aggregation...\n")

  res <- purrr::map(col_chunks, function(chunk) {
    se_chunk <- se[, chunk]
    scuttle::aggregateAcrossCells(
      se_chunk,
      ids = colData(se_chunk)$sample_id,
      store.number = TRUE,
      BPPARAM = BiocParallel::MulticoreParam(workers = parallel::detectCores())
    )
  }, .progress = TRUE)

  cat("Combining chunks...\n")
  se_pb <- do.call(cbind, res)

  # Persist intermediate (optional)
  try({
    saveRDS(se_pb, "~/Documents/GitHub/article_cellNexus/pseudobulk_se_X_Y_subset_sample_id.rds")
  }, silent = TRUE)

  # Filter out samples with zero counts
  cat("Filtering samples with zero counts...\n")
  se_pb <- se_pb[, colSums(assay(se_pb, "counts")) > 0]

  # Create logCPM assay
  cat("Calculating logCPM values...\n")
  assay(se_pb, "logCPM") <- edgeR::cpm(assay(se_pb, "counts"), log = TRUE, prior.count = 1)

  # Rank genes using singscore
  library(singscore)
  cat("Calculating gene ranks using singscore...\n")
  assay(se_pb, "singscore") <- assay(se_pb, "logCPM") |> singscore::rankGenes()

  # Keep sex chromosomes only
  cat("Filtering to sex chromosomes only...\n")
  se_pb <- se_pb[seqnames(rowRanges(se_pb)) %in% c("chrX", "chrY"), ]
  cat("Retained", nrow(se_pb), "sex chromosome genes\n")

  # Build tissue-specific validation split on full processed SE
  cat("Creating tissue-specific validation split...\n")
  split_info <- create_tissue_validation_split(
    se = se_pb,
    sex_column = sex_column,
    tissue_column = tissue_column,
    seed = seed
  )

  se_validation <- se_pb[, colnames(se_pb) %in% split_info$validation_ids]
  colData(se_validation)[[sex_column]] <- droplevels(colData(se_validation)[[sex_column]])

  # Labeled data for train/test (exclude unknown)
  valid_sex_labels <- c("female", "male", "Female", "Male", "F", "M", "f", "m")
  labeled_mask <- !is.na(colData(se_pb)[[sex_column]]) &
    as.character(colData(se_pb)[[sex_column]]) %in% valid_sex_labels &
    !(colnames(se_pb) %in% split_info$validation_ids)

  se_labeled <- se_pb[, labeled_mask]

  # Create initial train/test split based on column metadata
  meta_labeled <- as.data.frame(colData(se_labeled))
  meta_labeled$.__sample_id__ <- colnames(se_labeled)

  rs <- rsample::initial_split(
    meta_labeled,
    prop = train_prop,
    strata = !!rlang::sym(sex_column)
  )
  train_ids <- rsample::training(rs)$.__sample_id__
  test_ids  <- rsample::testing(rs)$.__sample_id__

  se_train <- se_labeled[, colnames(se_labeled) %in% train_ids]
  se_test  <- se_labeled[, colnames(se_labeled) %in% test_ids]

  # Optional downsampling to cap size
  if (ncol(se_train) > 10000) {
    cat("Training set is large, sampling 10000 samples to avoid memory issues...\n")
    set.seed(seed)
    sample_cols <- sample(colnames(se_train), 10000)
    se_train <- se_train[, sample_cols]
    cat("Sampled training SE dimensions:", dim(se_train), "\n")
  }

  cat("Training SE dimensions:", dim(se_train), "\n")
  cat("Test SE dimensions:", dim(se_test), "\n")
  cat("Validation SE dimensions:", dim(se_validation), "\n")

  invisible(list(
    se_train = se_train,
    se_test = se_test,
    se_validation = se_validation,
    se_processed = se_pb
  ))
}

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
#' Train a robust sex prediction model using tidymodels/randomForest
#'
#' Builds a preprocessing recipe and random forest classifier on sex-chromosome
#' features and optional batch covariates, with optional tuning of `mtry` and
#' cross-validation. Can accept an external test set; otherwise, performs an
#' internal split.
#'
#' @param se A `SummarizedExperiment` for training
#' @param se_test Optional `SummarizedExperiment` for testing (default: NULL)
#' @param my_assay Name of assay to use (default: "logCPM")
#' @param sex_column Outcome column name in metadata (default: "sex")
#' @param batch_columns Character vector of additional metadata predictors
#' @param train_prop Train proportion if internal split is used (default: 0.8)
#' @param trees Number of trees (default: 1000)
#' @param min_n Minimum observations per terminal node (default: 5)
#' @param mtry Optional mtry value; if NULL, tuned via grid (default: NULL)
#' @param seed Random seed (default: 123)
#'
#' @return A list with trained `workflow`, `recipe`, metrics and predictor metadata
#'
#' @examples
#' \dontrun{
#' fit <- build_sex_prediction_model(se_train, se_test = se_test, my_assay = "singscore")
#' }
#'
#' @importFrom dplyr select filter mutate
#' @importFrom tibble as_tibble
#' @importFrom yardstick metrics roc_auc accuracy sensitivity specificity
#' @importFrom rsample initial_split training testing vfold_cv
#' @importFrom recipes recipe step_impute_median step_unknown step_novel step_zv step_normalize step_corr all_numeric_predictors all_nominal_predictors all_predictors
#' @importFrom parsnip rand_forest set_engine set_mode fit
#' @importFrom workflows workflow add_recipe add_model
#' @importFrom tune tune tune_grid grid_regular select_best finalize_workflow control_grid
#' @importFrom tune fit_resamples
#' @importFrom dplyr bind_cols
build_sex_prediction_model <- function(se, se_test = NULL, my_assay = "singscore", sex_column = "sex", 
                                       batch_columns = c("assay", "tissue_groups", "age_days"), 
                                       train_prop = 0.8, trees = 500, seed = 123) {
  
  set.seed(seed)
  
  # Extract expression data from specified assay
  expr_data <- t(assay(se, my_assay))
  
  # Get metadata
  metadata <- colData(se) %>% as.data.frame()
  
  # Filter out "unknown" sex labels and standardize to female/male
  valid_sex <- metadata[[sex_column]] %in% c("female", "male")
  if (sum(valid_sex) == 0) {
    stop("No valid sex labels (female/male) found in the dataset")
  }
  
  # Subset data to valid sex labels
  expr_data <- expr_data[valid_sex, ]
  metadata <- metadata[valid_sex, ]
  
  # Standardize sex labels
  metadata[[sex_column]] <- factor(metadata[[sex_column]], levels = c("female", "male"))
  
  # Combine expression data with batch columns
  train_data <- expr_data %>%
    as.data.frame() %>%
    bind_cols(metadata %>% select(all_of(c(sex_column, batch_columns))))
  
  # Check for required columns
  missing_cols <- setdiff(c(sex_column, batch_columns), colnames(train_data))
  if (length(missing_cols) > 0) {
    stop("The required columns ", paste(missing_cols, collapse = ", "), " are missing.")
  }
  
  # Remove rows with any NA values in key columns
  complete_cases <- complete.cases(train_data %>% select(all_of(c(sex_column, batch_columns))))
  if (sum(complete_cases) == 0) {
    stop("No complete cases found after removing NA values")
  }
  
  train_data <- train_data[complete_cases, ]
  
  # Check for infinite values and replace with NA
  numeric_cols <- sapply(train_data, is.numeric)
  for (col in names(numeric_cols)[numeric_cols]) {
    train_data[[col]][is.infinite(train_data[[col]])] <- NA
  }
  
  # Remove rows with any remaining NA values
  complete_cases <- complete.cases(train_data)
  if (sum(complete_cases) == 0) {
    stop("No complete cases found after handling infinite values")
  }
  
  train_data <- train_data[complete_cases, ]
  
  # Final check: ensure no "unknown" or NA values in sex column
  train_data <- train_data %>%
    filter(!is.na(!!sym(sex_column)) & !!sym(sex_column) != "unknown")
  
  if (nrow(train_data) == 0) {
    stop("No data remaining after filtering out unknown/NA sex labels")
  }
  
  # Convert sex to factor with only female/male levels
  train_data[[sex_column]] <- factor(train_data[[sex_column]], levels = c("female", "male"))
  
  cat("Final class distribution:", paste(names(table(train_data[[sex_column]])), "=", table(train_data[[sex_column]]), collapse = ", "), "\n")
  
  # Identify feature columns (exclude target and batch columns)
  feature_cols <- setdiff(colnames(train_data), c(sex_column, batch_columns))
  
  # Remove zero-variance features
  feature_data <- train_data[, feature_cols, drop = FALSE]
  zero_var_features <- sapply(feature_data, function(x) {
    if (is.numeric(x)) {
      var(x, na.rm = TRUE) == 0
    } else {
      length(unique(x)) <= 1
    }
  })
  
  if (sum(zero_var_features) > 0) {
    cat("Removing", sum(zero_var_features), "zero-variance features\n")
    feature_cols <- feature_cols[!zero_var_features]
  }
  
  cat("Using", length(feature_cols), "features for prediction\n")
  
  # Create recipe
  recipe_formula <- as.formula(paste(sex_column, "~", paste(c(feature_cols, batch_columns), collapse = " + ")))
  
  recipe <- recipe(recipe_formula, data = train_data) %>%
    step_impute_median(all_numeric_predictors()) %>%
    step_zv(all_predictors()) %>%
    step_normalize(all_numeric_predictors()) %>%
    step_corr(all_numeric_predictors(), threshold = 0.95)
  
  # Create random forest model
  rf_model <- rand_forest(
    trees = trees,
    mtry = tune()
  ) %>%
    set_engine("randomForest") %>%
    set_mode("classification")
  
  # Create workflow
  workflow <- workflow() %>%
    add_recipe(recipe) %>%
    add_model(rf_model)
  
  # Determine cross-validation strategy
  min_class_size <- min(table(train_data[[sex_column]]))
  n_folds <- min(5, min_class_size)  # Use fewer folds for small datasets
  
  if (n_folds < 2) {
    warning("Very few samples in smallest class (", min_class_size, "). Consider collecting more data.")
    n_folds <- 2
  }
  
  if (nrow(train_data) < 20) {
    warning("Dataset too small for cross-validation. Using simple validation.")
    # Use simple train/test split instead
    set.seed(seed)
    split_indices <- sample(1:nrow(train_data), size = floor(0.8 * nrow(train_data)))
    train_split <- train_data[split_indices, ]
    test_split <- train_data[-split_indices, ]
    
    # Fit model on training data
    fitted_workflow <- workflow %>%
      fit(data = train_split)
    
    # Make predictions on test data
    predictions <- predict(fitted_workflow, new_data = test_split, type = "prob")
    predictions$.pred_class <- predict(fitted_workflow, new_data = test_split)$.pred_class
    predictions$truth <- test_split[[sex_column]]
    
    # Calculate metrics
    metrics <- predictions %>%
      metrics(truth = truth, estimate = .pred_class) %>%
      bind_rows(
        roc_auc(predictions, truth = truth, .pred_female)
      )
    
    return(list(
      workflow = fitted_workflow,
      metrics = metrics,
      predictions = predictions,
      train_data = train_split,
      test_data = test_split
    ))
  }
  
  # Create cross-validation folds
  cv_folds <- vfold_cv(train_data, v = n_folds, strata = sex_column)
  
  cat("Training samples:", nrow(train_data), "| Test samples:", ifelse(is.null(se_test), "None", ncol(se_test)), "\n")
  cat("Class distribution in training:", paste(names(table(train_data[[sex_column]])), "=", table(train_data[[sex_column]]), collapse = ", "), "\n")
  cat("Using", n_folds, "fold cross-validation\n")
  
  # Define tuning grid
  mtry_values <- seq(1, min(20, length(feature_cols)), by = 1)
  if (length(mtry_values) == 0) mtry_values <- 1
  
  tuning_grid <- grid_regular(
    mtry(range = c(1, min(20, length(feature_cols)))),
    levels = min(3, length(mtry_values))
  )
  
  cat("Tuning hyperparameters...\n")
  cat("Testing", nrow(tuning_grid), "hyperparameter combinations with", n_folds, "fold CV...\n")
  
  # Configure parallel processing safely
  n_workers <- parallel::detectCores() - 1  # Use fewer workers to avoid conflicts
  cat("Using", n_workers, "workers for parallel processing\n")
  
  # Set up parallel processing with reduced memory usage
  plan(multisession, workers = n_workers)
  
  # Tune the model with controlled parallel processing
  rf_tune_results <- workflow %>%
    tune_grid(
      resamples = cv_folds,
      grid = tuning_grid,
      metrics = metric_set(accuracy, roc_auc, sensitivity, specificity),
      control = control_resamples(
        allow_par = TRUE, 
        verbose = TRUE,
        parallel_over = "everything"  # Only parallelize over resamples, not everything
      )
    )
  
  # Check if tuning was successful
  if (nrow(rf_tune_results) == 0 || all(is.na(rf_tune_results$.metrics))) {
    stop("All models failed during tuning. Check data quality and try with fewer features.")
  }
  
  # Select best model
  best_model <- select_best(rf_tune_results, metric = "roc_auc")
  
  # Finalize workflow with best parameters
  final_workflow <- workflow %>%
    finalize_workflow(best_model)
  
  # Fit final model on full training data
  final_fit <- final_workflow %>%
    fit(data = train_data)
  
  # Make predictions on training data
  train_predictions <- predict(final_fit, new_data = train_data, type = "prob")
  train_predictions$.pred_class <- predict(final_fit, new_data = train_data)$.pred_class
  train_predictions$truth <- train_data[[sex_column]]
  
  # Calculate training metrics
  train_metrics <- train_predictions %>%
    metrics(truth = truth, estimate = .pred_class) %>%
    bind_rows(
      roc_auc(train_predictions, truth = truth, .pred_female)
    )
  
  # If test data is provided, evaluate on it
  test_metrics <- NULL
  test_predictions <- NULL
  
  if (!is.null(se_test)) {
    # Prepare test data
    test_expr_data <- t(assay(se_test, my_assay))
    test_metadata <- colData(se_test) %>% as.data.frame()
    
    # Filter test data
    test_valid_sex <- test_metadata[[sex_column]] %in% c("female", "male")
    if (sum(test_valid_sex) > 0) {
      test_expr_data <- test_expr_data[test_valid_sex, ]
      test_metadata <- test_metadata[test_valid_sex, ]
      test_metadata[[sex_column]] <- factor(test_metadata[[sex_column]], levels = c("female", "male"))
      
      test_data <- test_expr_data %>%
        as.data.frame() %>%
        bind_cols(test_metadata %>% select(all_of(c(sex_column, batch_columns))))
      
      # Apply same preprocessing as training data
      test_complete_cases <- complete.cases(test_data %>% select(all_of(c(sex_column, batch_columns))))
      test_data <- test_data[test_complete_cases, ]
      
      # Handle infinite values
      for (col in names(numeric_cols)[numeric_cols]) {
        if (col %in% colnames(test_data)) {
          test_data[[col]][is.infinite(test_data[[col]])] <- NA
        }
      }
      
      test_complete_cases <- complete.cases(test_data)
      test_data <- test_data[test_complete_cases, ]
      
      if (nrow(test_data) > 0) {
        test_predictions <- predict(final_fit, new_data = test_data, type = "prob")
        test_predictions$.pred_class <- predict(final_fit, new_data = test_data)$.pred_class
        test_predictions$truth <- test_data[[sex_column]]
        
        test_metrics <- test_predictions %>%
          metrics(truth = truth, estimate = .pred_class) %>%
          bind_rows(
            roc_auc(test_predictions, truth = truth, .pred_female)
          )
      }
    }
  }
  
  return(list(
    workflow = final_fit,
    train_metrics = train_metrics,
    test_metrics = test_metrics,
    train_predictions = train_predictions,
    test_predictions = test_predictions,
    tune_results = rf_tune_results,
    best_params = best_model
  ))
}

# Predict function: add sex_predicted to colData of a SummarizedExperiment
#' Predict sex for all samples and write result to colData(se)$sex_predicted
#'
#' Takes a trained model (as returned by `build_sex_prediction_model`) and a
#' `SummarizedExperiment`, constructs the predictor matrix from metadata and the
#' specified assay, reconciles missing predictors, and writes a `sex_predicted`
#' column into `colData(se)`.
#'
#' @param se A `SummarizedExperiment` to annotate with predictions
#' @param trained_model Object returned by `build_sex_prediction_model`
#' @param my_assay Assay name to use for prediction (default: "singscore")
#' @param sex_column Outcome column name (used for factor level alignment)
#'
#' @return The same `SummarizedExperiment` with `colData(se)$sex_predicted`
#'
#' @examples
#' \dontrun{
#' se_pred <- predict_sex_for_se(se_validation, sex_model_fitted, my_assay = "singscore")
#' }
#'
#' @importFrom SummarizedExperiment colData assay
#' @importFrom tibble as_tibble
#' @importFrom dplyr select
predict_sex_for_se <- function(se,
                               trained_model,
                               my_assay = "singscore",
                               sex_column = "sex") {
  stopifnot("workflow" %in% class(trained_model))
  wf <- trained_model

  # Build new_data from SE
  new_data <- cbind(
    colData(se) %>% as.data.frame(),
    t(assay(se, my_assay))
  ) %>% as_tibble()

  # Get the feature names from the workflow's recipe
  recipe <- extract_preprocessor(wf)
  feature_names <- recipe$var_info$variable[recipe$var_info$role == "predictor"]
  
  # Ensure required predictors exist; create placeholders when missing
  missing_cols <- setdiff(feature_names, colnames(new_data))
  if (length(missing_cols) > 0) {
    cat("Warning: Missing columns in prediction data:", length(missing_cols), "columns\n")
    cat("First few missing columns:", head(missing_cols), "\n")
    
    # Add missing columns with NA values
    for (col in missing_cols) {
      new_data[[col]] <- NA_real_
    }
  }

  # Keep only predictors used by the model
  new_data <- new_data %>% dplyr::select(dplyr::all_of(feature_names))

  # Generate class predictions
  class_pred <- predict(wf, new_data = new_data, type = "class")
  colData(se)$sex_predicted <- class_pred$.pred_class
  se
}


# Load tidymodels libraries
library(tidymodels)
library(randomForest)
library(vip)
library(dplyr)

# Set memory limits for parallel processing
options(future.globals.maxSize = 4 * 1024^3)  # 4GB limit - more conservative
options(future.plan = "multisession", workers = 4)  # Use 4 workers instead of 10

 

## === Part 1: Preprocessing to obtain train / test / validation ===
cat("Reading H5AD file...\n")
se_raw <- zellkonverter::readH5AD("~/Downloads/pseudobulk_se.h5ad", reader = "R", use_hdf5 = TRUE)

splits <- preprocess_for_sex_prediction(
  se = se_raw,
  genome_version = "hg38",
  my_assay = "singscore",
  chunk_size = 1000,
  sex_column = "sex",
  tissue_column = "tissue_groups",
  train_prop = 0.8,
  seed = 123
)

# Save the splits
saveRDS(splits, "~/Documents/GitHub/article_cellNexus/sex_prediction_splits.rds")

splits <- readRDS("~/Documents/GitHub/article_cellNexus/sex_prediction_splits.rds")

se_train <- splits$se_train
se_test <- splits$se_test
se_validation <- splits$se_validation


## === Part 2: Train model and keep trained workflow for prediction ===
cat("Building sex prediction model on training data...\n")
sex_model_results <- build_sex_prediction_model(
  se = se_train,
  se_test = se_test,
  my_assay = "singscore",
  sex_column = "sex",
  batch_columns = c("assay", "tissue_groups", "age_days"),
  train_prop = 0.8,
  trees = 500,
  seed = 123
)

# save model
saveRDS(sex_model_results, "~/Documents/GitHub/article_cellNexus/sex_prediction_model.rds")

# Print training metrics
cat("\n=== Training Performance ===\n")
print(sex_model_results$train_metrics)

if (!is.null(sex_model_results$test_metrics)) {
  cat("\n=== Test Performance ===\n")
  print(sex_model_results$test_metrics)
}


# === Model Performance Summary ===
# The following results summarize the model's performance on both cross-validation and test sets.
#
# Cross-validation results (mean ± sd):
#   - accuracy:    0.971 ± 0.00178
#   - roc_auc:     0.982 ± 0.000965
#   - sensitivity: 0.984 ± 0.00152
#   - specificity: 0.955 ± 0.00302
#
# Test set performance:
#   - accuracy: 0.970
#   - kappa:    0.939
#   - roc_auc:  0.977
#
# Model training details:
#   - Training samples: 7899
#   - Test samples:     1975
#   - Features used:    371
#   - Cross-validation folds: 5


## === Part 3: Prediction function and validation application ===
cat("\n=== Tissue-Specific Validation ===\n")
cat("Validating on sex-specific tissues...\n")

# Add sex_predicted to validation SE using the helper function
se_validation_pred <- predict_sex_for_se(
  se = se_validation,
  trained_model = sex_model_results$workflow,
  my_assay = "singscore",
  sex_column = "sex"
)

# Build evaluation tibble
validation_predictions <- colData(se_validation_pred) %>%
  as.data.frame() %>%
  as_tibble() %>%
  mutate(.pred_class = factor(sex_predicted, levels = levels(sex)))


cat("Validation predictions after filtering unknown sex:", nrow(validation_predictions), "\n")
cat("Sex levels in validation predictions:", paste(levels(validation_predictions$sex), collapse = ", "), "\n")
  
  # Calculate validation metrics (classification only here)
  validation_metrics <- validation_predictions %>%
    metrics(truth = sex, estimate = .pred_class)
  
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

# === ROC Curve and Feature Importance ===
cat("\n=== ROC Curve and Feature Importance ===\n")

# Use training predictions from model results
train_predictions <- sex_model_results$train_predictions

# Compute ROC curve data using training predictions
roc_data <- roc_curve(train_predictions, truth = truth, .pred_female)

plot_roc = ggplot(roc_data, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(color = "#377EB8", size = 1.2) +
  geom_abline(linetype = "dashed", color = "gray") +
  labs(
    title = "ROC Curve on Training Data",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  theme_minimal()

# Plot feature importance
library(ggplot2)
library(org.Hs.eg.db)
library(AnnotationDbi)

 # Extract feature importance from the model
importance_scores <- extract_fit_engine(sex_model_results$workflow)$importance
  # Convert to data frame for plotting
  importance_df <- data.frame(
    gene = rownames(importance_scores),
    importance = importance_scores[, "MeanDecreaseGini"]
  ) %>%
    arrange(desc(importance)) %>%
    head(20)  # Top 20 genes

  plot_importance <- 
    importance_df |>
    # Clear ensemble _X prefix using tidyverse
    mutate(gene = str_remove(gene, "_X$")) |>
    # add gene symbol using annotation
    mutate(symbol = mapIds(org.Hs.eg.db,
                          keys = gene,
                          keytype = "ENSEMBL",
                          column = "SYMBOL",
                          multiVals = "first"
    ), symbol = if_else(is.na(symbol), gene, symbol)
    ) |>

    ggplot(aes(x = reorder(symbol, importance), y = importance)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      coord_flip() +
      labs(
        title = "Top 20 Most Important Genes for Sex Prediction",
        x = "Gene",
        y = "Importance (Mean Decrease Gini)"
      ) +
      theme_minimal() +
     theme(axis.text.y = element_text(size = 8))


detach("package:org.Hs.eg.db", unload = TRUE)

# Now predict for the unknown samples
# remotes::install_github("MangiolaLaboratory/cellNexus", force = TRUE)
# library(cellNexus)

# get the unknown samples
#se_unknown_samples <- 
#get_metadata() |>
#filter(sex == "unknown" | is.na(sex)) |>
#get_pseudobulk()

# Save model for future use (optional)
# saveRDS(sex_model_fitted, "sex_prediction_model.rds")

cat("\n", strrep("=", 60), "\n")
cat("ANALYSIS COMPLETED SUCCESSFULLY\n")
cat(strrep("=", 60), "\n")
cat("Final model saved and ready for predictions\n")
# now select the sex == unknown and predict the sex using the helper function
se_all <- splits$se_processed
se_unknown <- se_all[, colData(se_all)$sex == "unknown" | is.na(colData(se_all)$sex)]

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


# Add predictions directly to SE (use filtered version)
se_unknown_filtered <- predict_sex_for_se(
  se = se_unknown_filtered,
  trained_model = sex_model_results$workflow,
  my_assay = "singscore",
  sex_column = "sex"
)

# count tissue_groups and sex_prediction
colData(se_unknown_filtered) |>
 as_tibble() |> 
  dplyr::count(tissue_groups, sex_predicted, age_days) |> 
 print(n=99)

# Plot bar chart of predicted sex for unknown samples across tissues
# Use red/blue color scheme and reorder tissues within ggplot call
library(forcats)

# Plot bar chart of predicted sex for unknown samples across tissues
plot_unknown_prediction <- 
  colData(se_unknown_filtered) |> 
  as_tibble() |>
  dplyr::count(tissue_groups, sex_predicted) |>
  arrange(desc(n)) |>
  ggplot(aes(x = reorder(tissue_groups, n, decreasing = TRUE), y = n, fill = sex_predicted)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("female" = "#E41A1C", "male" = "#377EB8")) +
    labs(
      title = "Prediction of Unknown Samples Across Tissues",
      x = "Tissue",
      y = "Count",
      fill = "Predicted Sex"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # save the unknown counts data
  colData(se_unknown_filtered) |> 
    as.data.frame() |>
    as_tibble() |>
    count(tissue_groups, sex_predicted) |>
    arrange(desc(n)) |>
    write.csv("~/Documents/GitHub/article_cellNexus/unknown_prediction_counts.csv", row.names = FALSE)
} else {
  cat("Warning: sex_predicted column not found in se_unknown_filtered. Skipping unknown prediction plot.\n")
  plot_unknown_prediction <- NULL
}

# Use patchwork to combine the plots
library(patchwork)

p = plot_roc + plot_importance + plot_unknown_prediction 

# legend grouped to bottom
p = p + plot_layout(guides = "collect")

ggsave(
  "~/Documents/GitHub/article_cellNexus/sex_prediction_plots.pdf", 
  p, units = "mm", 
  width = 90, height = 45, scale = 3
)
