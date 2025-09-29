#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(scater)
  library(scran)
  library(igraph)
  library(HDF5Array)
})

safe_require <- function(pkg) {
  suppressWarnings(suppressPackageStartupMessages(requireNamespace(pkg, quietly = TRUE)))
}

message_time <- function(...) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(sprintf("[%s] %s", ts, paste0(..., collapse = " ")))
}

detect_species <- function(gene_ids) {
  subset_ids <- head(gene_ids, 1000)
  if (any(grepl("^ENSG", subset_ids))) return("human")
  if (any(grepl("^ENSMUSG", subset_ids))) return("mouse")
  # Fallback heuristic on gene symbols
  human_like <- sum(grepl("^[A-Z0-9]+$", subset_ids))
  mouse_like <- sum(grepl("^[A-Z][a-z0-9]+$", subset_ids))
  if (human_like >= mouse_like) "human" else "mouse"
}

add_qc_metrics <- function(sce, species = "auto") {
  if (species == "auto") species <- detect_species(rownames(sce))
  
  # Define gene patterns by species
  if (species == "human") {
    ribo_pattern <- "^RP[SL]"
    mito_pattern <- "^MT-"
  } else if (species == "mouse") {
    ribo_pattern <- "^Rp[sl]"
    mito_pattern <- "^mt-"
  } else {
    # Generic patterns for other species
    ribo_pattern <- "^[Rr][Pp][SsLl]"
    mito_pattern <- "^[Mm][Tt]-?"
  }
  
  message_time(sprintf("Identifying ribosomal genes using pattern: %s", ribo_pattern))
  message_time(sprintf("Identifying mitochondrial genes using pattern: %s", mito_pattern))
  
  # Find ribosomal genes
  ribo_genes <- grep(ribo_pattern, rownames(sce), value = TRUE)
  
  if (length(ribo_genes) == 0) {
    # Try alternative patterns
    alt_patterns <- c("^RPS", "^RPL", "^Rps", "^Rpl")
    for (pattern in alt_patterns) {
      ribo_genes <- grep(pattern, rownames(sce), value = TRUE)
      if (length(ribo_genes) > 0) {
        message_time(sprintf("Found %d ribosomal genes using pattern: %s", length(ribo_genes), pattern))
        break
      }
    }
  } else {
    message_time(sprintf("Found %d ribosomal genes", length(ribo_genes)))
  }
  
  # Find mitochondrial genes
  mito_genes <- grep(mito_pattern, rownames(sce), value = TRUE)
  
  if (length(mito_genes) == 0) {
    # Try alternative patterns
    alt_mito_patterns <- c("^MT-", "^mt-", "^Mt-")
    for (pattern in alt_mito_patterns) {
      mito_genes <- grep(pattern, rownames(sce), value = TRUE)
      if (length(mito_genes) > 0) {
        message_time(sprintf("Found %d mitochondrial genes using pattern: %s", length(mito_genes), pattern))
        break
      }
    }
  } else {
    message_time(sprintf("Found %d mitochondrial genes", length(mito_genes)))
  }
  
  # Calculate QC metrics
  subsets_list <- list()
  if (length(ribo_genes) > 0) {
    subsets_list$Ribo <- ribo_genes
  }
  if (length(mito_genes) > 0) {
    subsets_list$Mito <- mito_genes
  }
  
  if (length(subsets_list) > 0) {
    qc_metrics <- scuttle::perCellQCMetrics(sce, subsets = subsets_list)
    
    if (length(ribo_genes) > 0) {
      sce$percent_ribo <- qc_metrics$subsets_Ribo_percent
      message_time(sprintf("Added ribosomal percentage (mean: %.2f%%)", mean(sce$percent_ribo)))
    } else {
      message_time("Warning: No ribosomal genes found. Setting percent_ribo to 0.")
      sce$percent_ribo <- 0
    }
    
    if (length(mito_genes) > 0) {
      sce$percent_mito <- qc_metrics$subsets_Mito_percent
      message_time(sprintf("Added mitochondrial percentage (mean: %.2f%%)", mean(sce$percent_mito)))
    } else {
      message_time("Warning: No mitochondrial genes found. Setting percent_mito to 0.")
      sce$percent_mito <- 0
    }
  } else {
    message_time("Warning: No ribosomal or mitochondrial genes found.")
    sce$percent_ribo <- 0
    sce$percent_mito <- 0
  }
  
  return(sce)
}

regress_qc_effects <- function(sce, assay_name = "logcounts") {
  # Check which QC metrics are available
  has_ribo <- "percent_ribo" %in% names(SummarizedExperiment::colData(sce))
  has_mito <- "percent_mito" %in% names(SummarizedExperiment::colData(sce))
  
  if (!has_ribo && !has_mito) {
    stop("Neither ribosomal nor mitochondrial percentages found. Run add_qc_metrics() first.")
  }
  
  covariates <- c()
  if (has_ribo) covariates <- c(covariates, "percent_ribo")
  if (has_mito) covariates <- c(covariates, "percent_mito")
  
  message_time(sprintf("Regressing out QC effects (%s) from expression data", paste(covariates, collapse = ", ")))
  
  # Get expression matrix and covariates
  expr_mat <- SummarizedExperiment::assay(sce, assay_name)
  covariate_data <- SummarizedExperiment::colData(sce)[, covariates, drop = FALSE]
  
  # Regress QC effects gene by gene
  if (safe_require("BiocParallel") && !is.null(bp)) {
    message_time("Using parallel processing for QC regression")
    corrected_expr <- BiocParallel::bplapply(seq_len(nrow(expr_mat)), function(i) {
      gene_expr <- expr_mat[i, ]
      if (stats::var(gene_expr) > 0) {
        # Build regression formula
        if (length(covariates) == 1) {
          covar_vals <- covariate_data[, 1]
          if (stats::var(covar_vals) > 0) {
            fit <- stats::lm(gene_expr ~ covar_vals)
            residuals(fit) + mean(gene_expr)  # Add back mean to preserve scale
          } else {
            gene_expr  # No correction if no variance in covariate
          }
        } else {
          # Multiple covariates
          df_temp <- data.frame(gene_expr = gene_expr, covariate_data)
          if (any(apply(covariate_data, 2, function(x) stats::var(x) > 0))) {
            fit <- stats::lm(gene_expr ~ ., data = df_temp)
            residuals(fit) + mean(gene_expr)
          } else {
            gene_expr
          }
        }
      } else {
        gene_expr  # No correction if no variance in gene expression
      }
    }, BPPARAM = bp)
    corrected_expr <- do.call(rbind, corrected_expr)
  } else {
    message_time("Using sequential processing for QC regression")
    corrected_expr <- t(apply(expr_mat, 1, function(gene_expr) {
      if (stats::var(gene_expr) > 0) {
        # Build regression formula
        if (length(covariates) == 1) {
          covar_vals <- covariate_data[, 1]
          if (stats::var(covar_vals) > 0) {
            fit <- stats::lm(gene_expr ~ covar_vals)
            residuals(fit) + mean(gene_expr)  # Add back mean to preserve scale
          } else {
            gene_expr  # No correction if no variance in covariate
          }
        } else {
          # Multiple covariates
          df_temp <- data.frame(gene_expr = gene_expr, covariate_data)
          if (any(apply(covariate_data, 2, function(x) stats::var(x) > 0))) {
            fit <- stats::lm(gene_expr ~ ., data = df_temp)
            residuals(fit) + mean(gene_expr)
          } else {
            gene_expr
          }
        }
      } else {
        gene_expr  # No correction if no variance in gene expression
      }
    }))
  }
  
  rownames(corrected_expr) <- rownames(expr_mat)
  colnames(corrected_expr) <- colnames(expr_mat)
  
  # Add corrected expression as new assay
  SummarizedExperiment::assay(sce, "logcounts_qc_regressed") <- corrected_expr
  message_time(sprintf("Added 'logcounts_qc_regressed' assay with %s effects removed", paste(covariates, collapse = " and ")))
  
  return(sce)
}

get_matrix <- function(sce, assay_name = c("logcounts", "counts")) {
  assay_name <- assay_name[assay_name %in% SummarizedExperiment::assayNames(sce)][1]
  if (is.na(assay_name)) stop("Neither 'logcounts' nor 'counts' assays are available.")
  SummarizedExperiment::assay(sce, assay_name)
}

collapse_ensembl_to_symbol <- function(sce, species) {
  ids <- rownames(sce)
  if (!any(grepl("^ENS", ids))) return(list(expr = get_matrix(sce), genes = ids))

  if (species == "human") {
    if (!safe_require("AnnotationDbi") || !safe_require("org.Hs.eg.db")) {
      message_time("org.Hs.eg.db not installed; proceeding with Ensembl IDs (matching to ref may be poor).")
      return(list(expr = get_matrix(sce), genes = ids))
    }
    map <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = ids, columns = c("SYMBOL"), keytype = "ENSEMBL")
  } else {
    if (!safe_require("AnnotationDbi") || !safe_require("org.Mm.eg.db")) {
      message_time("org.Mm.eg.db not installed; proceeding with Ensembl IDs (matching to ref may be poor).")
      return(list(expr = get_matrix(sce), genes = ids))
    }
    map <- AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys = ids, columns = c("SYMBOL"), keytype = "ENSEMBL")
  }

  map <- map[!is.na(map$SYMBOL) & map$ENSEMBL %in% ids, , drop = FALSE]
  if (nrow(map) == 0) {
    message_time("No SYMBOL mappings found; proceeding with Ensembl IDs.")
    return(list(expr = get_matrix(sce), genes = ids))
  }

  expr <- get_matrix(sce)
  # Choose one Ensembl per SYMBOL based on highest average expression
  by_symbol <- split(map$ENSEMBL, map$SYMBOL)
  keep_ensembl <- vapply(by_symbol, function(ens_vec) {
    if (length(ens_vec) == 1L) return(ens_vec)
    cand <- ens_vec[ens_vec %in% rownames(expr)]
    if (length(cand) == 0L) return(ens_vec[1L])
    means <- rowMeans(expr[cand, , drop = FALSE])
    cand[which.max(means)]
  }, character(1))

  picked <- unique(keep_ensembl)
  expr2 <- expr[picked, , drop = FALSE]
  # Set rownames to corresponding SYMBOL
  sym_lookup <- setNames(names(keep_ensembl), keep_ensembl)
  new_genes <- unname(sym_lookup[rownames(expr2)])
  rownames(expr2) <- new_genes
  # Keep unique symbols only
  dup <- duplicated(rownames(expr2))
  if (any(dup)) expr2 <- expr2[!dup, , drop = FALSE]
  list(expr = expr2, genes = rownames(expr2))
}


find_top_markers <- function(sce, groups, n_top = 10L, lfc = 0.5, BPPARAM = NULL, assay_name = "logcounts") {
  n_clusters <- length(unique(groups))
  message_time(sprintf("Finding markers per cluster (direction = 'up') for %d clusters", n_clusters))
  if (!is.null(BPPARAM)) {
    message_time(sprintf("Using %d cores for parallel marker detection", BiocParallel::bpnworkers(BPPARAM)))
    res <- scran::findMarkers(sce, groups = groups, direction = "up", lfc = lfc, BPPARAM = BPPARAM, assay.type = assay_name)
  } else {
    message_time("Running marker detection on single core")
    res <- scran::findMarkers(sce, groups = groups, direction = "up", lfc = lfc, assay.type = assay_name)
  }
  message_time("Processing marker results")
  cluster_ids <- names(res)
  out_list <- lapply(seq_along(cluster_ids), function(i) {
    cl <- cluster_ids[i]
    if (i %% 5 == 0 || i == length(cluster_ids)) {
      message_time(sprintf("Processing cluster %s (%d/%d)", cl, i, length(cluster_ids)))
    }
    df <- as.data.frame(res[[cl]])
    df$gene <- rownames(res[[cl]])
    # Order robustly: FDR then AUC (if present) then summary.logFC (desc)
    ord <- order(df$FDR, if ("AUC" %in% names(df)) -df$AUC else df$FDR, -df$summary.logFC)
    df2 <- df[ord, , drop = FALSE]
    head_df <- head(df2, min(n_top, nrow(df2)))
    data.frame(
      cluster = cl,
      gene = head_df$gene,
      FDR = head_df$FDR,
      summary.logFC = head_df$summary.logFC,
      AUC = if ("AUC" %in% names(head_df)) head_df$AUC else NA_real_,
      row.names = NULL
    )
  })
  message_time("Marker detection completed")
  do.call(rbind, out_list)
}

singleR_labels <- function(sce, species, BPPARAM = NULL) {
  if (!safe_require("SingleR") || !safe_require("celldex")) {
    message_time("SingleR/celldex not installed; skipping SingleR labeling.")
    return(NULL)
  }
  converted <- collapse_ensembl_to_symbol(sce, species)
  expr2 <- converted$expr
  # Choose reference
  ref <- NULL
  ref_labels_col <- NULL
  if (species == "human") {
    if (safe_require("celldex")) {
      # Prefer immune-focused reference if available
      if ("MonacoImmuneData" %in% getNamespaceExports("celldex")) {
        ref <- celldex::MonacoImmuneData()
        ref_labels_col <- if ("label.fine" %in% names(SummarizedExperiment::colData(ref))) "label.fine" else "label.main"
      } else if ("BlueprintEncodeData" %in% getNamespaceExports("celldex")) {
        ref <- celldex::BlueprintEncodeData()
        ref_labels_col <- if ("label.fine" %in% names(SummarizedExperiment::colData(ref))) "label.fine" else "label.main"
      } else if ("HumanPrimaryCellAtlasData" %in% getNamespaceExports("celldex")) {
        ref <- celldex::HumanPrimaryCellAtlasData()
        ref_labels_col <- if ("label.fine" %in% names(SummarizedExperiment::colData(ref))) "label.fine" else "label.main"
      }
    }
  } else {
    if (safe_require("celldex") && "MouseRNAseqData" %in% getNamespaceExports("celldex")) {
      ref <- celldex::MouseRNAseqData()
      ref_labels_col <- if ("label.fine" %in% names(SummarizedExperiment::colData(ref))) "label.fine" else "label.main"
    }
  }
  if (is.null(ref)) {
    message_time("No suitable celldex reference found; skipping SingleR labeling.")
    return(NULL)
  }
  message_time("Running SingleR for per-cell labels")
  if (!is.null(BPPARAM)) {
    pred <- SingleR::SingleR(test = expr2, ref = ref, labels = SummarizedExperiment::colData(ref)[[ref_labels_col]], BPPARAM = BPPARAM)
  } else {
    pred <- SingleR::SingleR(test = expr2, ref = ref, labels = SummarizedExperiment::colData(ref)[[ref_labels_col]])
  }
  labels <- pred$labels
  data.frame(cell = colnames(sce), label = labels, stringsAsFactors = FALSE)
}

consensus_labels_by_cluster <- function(cell_labels_df, clusters) {
  if (is.null(cell_labels_df)) return(NULL)
  labs <- split(cell_labels_df$label, clusters)
  cl <- names(labs)
  m <- lapply(labs, function(v) {
    tb <- sort(table(v), decreasing = TRUE)
    top <- names(tb)[1]
    prop <- as.numeric(tb[1]) / sum(tb)
    data.frame(label = top, proportion = prop, stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, m)
  res$cluster <- cl
  res[, c("cluster", "label", "proportion")]
}

marker_sets_immune <- function() {
  list(
    T_cell = c("CD3D","CD3E","CD2","TRAC","IL7R"),
    CD4_T = c("CD4","IL7R","CCR7","LTB"),
    CD8_T = c("CD8A","CD8B","GZMK","GZMB","PRF1","NKG7"),
    NK = c("NKG7","GNLY","KLRD1","KLRK1","PRF1"),
    B_cell = c("MS4A1","CD79A","CD79B","CD74","BANK1"),
    Plasma = c("SDC1","MZB1","XBP1","JCHAIN","PRDM1"),
    Monocyte = c("LYZ","S100A8","S100A9","VCAN","FCN1","CTSS"),
    Macrophage = c("C1QA","C1QB","C1QC","MRC1","CD68","MARCO","MSR1"),
    Dendritic = c("FCER1A","CST3","LILRA4","GZMB"),
    Neutrophil = c("S100A8","S100A9","CXCL8","FCGR3B"),
    Erythrocyte = c("HBB","HBA1","HBA2","ALAS2"),
    Platelet = c("PPBP","PF4","GP9","ITGA2B")
  )
}

heuristic_labels_from_markers <- function(sce, clusters, species, BPPARAM = NULL) {
  message_time("Computing heuristic labels from marker sets")
  # Ensure gene symbols for marker matching
  conv <- collapse_ensembl_to_symbol(sce, species)
  expr2 <- conv$expr
  genes <- rownames(expr2)

  sets <- marker_sets_immune()
  cl_levels <- levels(factor(clusters))
  message_time(sprintf("Scoring %d clusters against %d marker sets", length(cl_levels), length(sets)))
  
  # Parallelize cluster scoring if BPPARAM provided
  if (!is.null(BPPARAM) && safe_require("BiocParallel")) {
    message_time(sprintf("Using %d cores for heuristic labeling", BiocParallel::bpnworkers(BPPARAM)))
    scores <- BiocParallel::bplapply(cl_levels, function(cl) {
      idx <- which(clusters == cl)
      mean_expr <- rowMeans(expr2[, idx, drop = FALSE])
      set_scores <- vapply(names(sets), function(nm) {
        g <- intersect(sets[[nm]], genes)
        if (length(g) == 0L) return(NA_real_)
        mean(mean_expr[g])
      }, numeric(1))
      data.frame(cluster = cl, label = names(set_scores), score = set_scores, row.names = NULL)
    }, BPPARAM = BPPARAM)
  } else {
    scores <- lapply(cl_levels, function(cl) {
      idx <- which(clusters == cl)
      mean_expr <- rowMeans(expr2[, idx, drop = FALSE])
      set_scores <- vapply(names(sets), function(nm) {
        g <- intersect(sets[[nm]], genes)
        if (length(g) == 0L) return(NA_real_)
        mean(mean_expr[g])
      }, numeric(1))
      data.frame(cluster = cl, label = names(set_scores), score = set_scores, row.names = NULL)
    })
  }
  
  score_df <- do.call(rbind, scores)
  # Pick top label per cluster
  top_df <- do.call(rbind, lapply(split(score_df, score_df$cluster), function(df) df[which.max(df$score), , drop = FALSE]))
  message_time("Heuristic labeling completed")
  list(scores = score_df, top = top_df)
}

## Auto-load the HDF5-based SummarizedExperiment into macro_muscle when sourced
macro_muscle <- HDF5Array::loadHDF5SummarizedExperiment("macro_muscle_communication_5_tissues_hdf5")

## Step-by-step execution below

# Step-by-step execution
if (!safe_require("scuttle")) stop("Package 'scuttle' is required for normalization; please install it.")

# Configure parallelism
bp <- NULL
if (safe_require("BiocParallel")) {
  if (.Platform$OS.type == "windows") {
    bp <- BiocParallel::SnowParam(workers = 12, type = "SOCK")
  } else {
    bp <- BiocParallel::MulticoreParam(workers = 12)
  }
}

message_time("Preparing assays and detecting species")
species_detected <- detect_species(rownames(macro_muscle))
message_time("Species:", species_detected)

# Add QC metrics (ribosomal and mitochondrial) and regress out their effects
# macro_muscle <- add_qc_metrics(macro_muscle, species_detected)
# macro_muscle <- regress_qc_effects(macro_muscle, "logcounts")

# Use existing harmony clusters
message_time(sprintf("Using existing clusters: %d clusters found", length(unique(macro_muscle$cluster_harmony))))

# macro_muscle |> HDF5Array::quickResaveHDF5SummarizedExperiment()

# Find markers using QC-corrected expression
markers_df <- find_top_markers(macro_muscle, groups = macro_muscle$cluster_harmony, n_top = 15L, BPPARAM = bp, assay_name = "logcounts")
S4Vectors::metadata(macro_muscle)$cluster_markers <- markers_df

cell_labels_df <- singleR_labels(macro_muscle, species_detected, BPPARAM = bp)
cl_cons <- NULL
if (!is.null(cell_labels_df)) {
  macro_muscle$cell_label_singleR <- cell_labels_df$label[match(colnames(macro_muscle), cell_labels_df$cell)]
  cl_cons <- consensus_labels_by_cluster(cell_labels_df, macro_muscle$cluster_harmony)
  if (!is.null(cl_cons) && nrow(cl_cons) > 0) {
    label_map <- stats::setNames(cl_cons$label, cl_cons$cluster)
    macro_muscle$cluster_harmony_label_singleR <- label_map[as.character(macro_muscle$cluster_harmony)]
  }
}

heur <- heuristic_labels_from_markers(macro_muscle, macro_muscle$cluster_harmony, species_detected, BPPARAM = bp)
if (!is.null(heur)) {
  S4Vectors::metadata(macro_muscle)$cluster_marker_scores <- heur$scores
  if (!is.null(heur$top) && nrow(heur$top) > 0) {
    heur_map <- stats::setNames(heur$top$label, heur$top$cluster)
    macro_muscle$cluster_harmony_label_markers <- heur_map[as.character(macro_muscle$cluster_harmony)]
  }
}

# Package all results into a single list and save
results <- list(
  macro_muscle = macro_muscle,
  markers_df = markers_df,
  cell_labels_df = cell_labels_df,
  cl_cons = cl_cons,
  heur = heur
)
saveRDS(results, "macro_muscle_communication_5_tissues_hdf5_annotated.rds")

#results = readRDS("macro_muscle_communication_5_tissues_hdf5_annotated.rds")

# Analysis and Visualization
message_time("Analyzing cluster identities and creating visualizations")

# Summary of cluster identities
cluster_summary <- data.frame(
  cluster = levels(factor(results$macro_muscle$cluster_harmony)),
  n_cells = as.numeric(table(results$macro_muscle$cluster_harmony)),
  stringsAsFactors = FALSE
)

# Add marker-based labels
  cluster_summary$marker_label <- results$heur$top$label[match(cluster_summary$cluster, results$heur$top$cluster)]
  cluster_summary$marker_score <- round(results$heur$top$score[match(cluster_summary$cluster, results$heur$top$cluster)], 3)


# Add SingleR consensus labels if available
if (!is.null(results$cl_cons)) {
  cluster_summary$singleR_label <- results$cl_cons$label[match(cluster_summary$cluster, results$cl_cons$cluster)]
  cluster_summary$singleR_proportion <- round(results$cl_cons$proportion[match(cluster_summary$cluster, results$cl_cons$cluster)], 3)
}

# Get top 3 markers per cluster for interpretation (convert to gene symbols)
top_markers_per_cluster <- split(results$markers_df, results$markers_df$cluster)

# Get gene symbol mapping if using Ensembl IDs
species_detected <- detect_species(rownames(results$macro_muscle))
gene_mapping <- NULL
if (any(grepl("^ENS", results$markers_df$gene))) {
  if (species_detected == "human" && safe_require("AnnotationDbi") && safe_require("org.Hs.eg.db")) {
    tryCatch({
      gene_mapping <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, 
                                           keys = unique(results$markers_df$gene), 
                                           columns = c("SYMBOL"), 
                                           keytype = "ENSEMBL")
      gene_mapping <- gene_mapping[!is.na(gene_mapping$SYMBOL), ]
    }, error = function(e) NULL)
  } else if (species_detected == "mouse" && safe_require("AnnotationDbi") && safe_require("org.Mm.eg.db")) {
    tryCatch({
      gene_mapping <- AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, 
                                           keys = unique(results$markers_df$gene), 
                                           columns = c("SYMBOL"), 
                                           keytype = "ENSEMBL")
      gene_mapping <- gene_mapping[!is.na(gene_mapping$SYMBOL), ]
    }, error = function(e) NULL)
  }
}

cluster_summary$top_markers <- vapply(cluster_summary$cluster, function(cl) {
  if (cl %in% names(top_markers_per_cluster)) {
    top3_genes <- head(top_markers_per_cluster[[cl]]$gene, 3)
    # Convert to symbols if mapping available
    if (!is.null(gene_mapping)) {
      top3_symbols <- gene_mapping$SYMBOL[match(top3_genes, gene_mapping$ENSEMBL)]
      top3_symbols[is.na(top3_symbols)] <- top3_genes[is.na(top3_symbols)]  # Keep original if no mapping
      paste(top3_symbols, collapse = ", ")
    } else {
      paste(top3_genes, collapse = ", ")
    }
  } else {
    "No markers"
  }
}, character(1))

# Literature-based interpretation
interpret_cluster <- function(marker_label, singleR_label, top_markers) {
  # Priority: marker-based > SingleR > top markers
  if (!is.na(marker_label) && marker_label != "") {
    return(switch(marker_label,
      "Macrophage" = "Tissue-resident macrophages (C1Q+ phenotype)",
      "Monocyte" = "Inflammatory monocytes/monocyte-derived macrophages",
      "T_cell" = "T lymphocytes (CD3+)",
      "CD4_T" = "CD4+ T helper cells",
      "CD8_T" = "CD8+ T cytotoxic cells",
      "NK" = "Natural killer cells",
      "B_cell" = "B lymphocytes",
      "Plasma" = "Plasma cells/plasmablasts", 
      "Dendritic" = "Dendritic cells",
      "Neutrophil" = "Neutrophils",
      "Erythrocyte" = "Erythrocytes/red blood cells",
      "Platelet" = "Platelets/megakaryocytes",
      paste("Unknown", marker_label, "signature")
    ))
  } else if (!is.na(singleR_label) && singleR_label != "") {
    return(paste("SingleR:", singleR_label))
  } else {
    return(paste("Unknown (top markers:", strtrim(top_markers, 50), ")"))
  }
}

cluster_summary$literature_identity <- mapply(
  interpret_cluster, 
  if ("marker_label" %in% names(cluster_summary)) cluster_summary$marker_label else NA_character_,
  if ("singleR_label" %in% names(cluster_summary)) cluster_summary$singleR_label else NA_character_,
  cluster_summary$top_markers,
  USE.NAMES = FALSE
)

# Print summary
message_time("Cluster Identity Summary:")
print(cluster_summary)

# Create visualizations
message_time("Creating visualizations")

# 1. UMAP plots
if ("UMAP_HARMONY" %in% SingleCellExperiment::reducedDimNames(results$macro_muscle)) {
  
  # Plot by cluster number
  p1 <- scater::plotReducedDim(
    results$macro_muscle, 
    dimred = "UMAP_HARMONY", 
    colour_by = "cluster_harmony",
    point_size = 0.3,
    point_alpha = 0.7,
    rasterise = TRUE
  ) + 
  ggplot2::labs(title = "Clusters (Harmony UMAP)", color = "Cluster") +
  ggplot2::theme_minimal()

# List the gene symbol markers per cluster
for (cl in names(top_markers_per_cluster)) {
  marker_genes <- head(top_markers_per_cluster[[cl]]$gene, 30)  # Top 10 markers
  # Convert to symbols if mapping available
  if (!is.null(gene_mapping)) {
    marker_symbols <- gene_mapping$SYMBOL[match(marker_genes, gene_mapping$ENSEMBL)]
    marker_symbols[is.na(marker_symbols)] <- marker_genes[is.na(marker_symbols)]  # Keep original if no mapping
    message_time(sprintf("Cluster %s: %s", cl, paste(marker_symbols, collapse = ", ")))
  } else {
    message_time(sprintf("Cluster %s: %s", cl, paste(marker_genes, collapse = ", ")))
  }
}

# For cluster 4 plot a faceted UMAP with these markers
# PECAM1, VWF, KDR/FLT1, CLDN5 alongside SPARCL1.
p2 <- scater::plotReducedDim(
  results$macro_muscle, 
  dimred = "UMAP_HARMONY", 
  colour_by = "cluster_harmony_label_markers",
  point_size = 0.3,
  point_alpha = 0.7,
  rasterise = TRUE
) + 
ggplot2::labs(title = "Marker-based Labels (Harmony UMAP)", color = "Cell Type") +
ggplot2::theme_minimal() +
ggplot2::theme(legend.position = "bottom")

#' Plot gene expression markers on reduced dimensions
#'
#' @param sce SingleCellExperiment object
#' @param markers Character vector of gene symbols to plot
#' @param dimred Name of reduced dimension to use (default: "UMAP_HARMONY")
#' @param ncol Number of columns for plot arrangement (default: 3)
#' @param point_size Size of points (default: 0.2)
#' @param point_alpha Transparency of points (default: 0.6)
#' @param color_scale Color scale function (default: ggplot2::scale_color_viridis_c)
#' @param rasterise Whether to rasterise points (default: TRUE)
#' @param verbose Whether to print available markers (default: TRUE)
#'
#' @return Combined ggplot object with marker expression plots
plot_marker_expression <- function(sce, 
                                   markers, 
                                   dimred = "UMAP_HARMONY",
                                   ncol = 3,
                                   point_size = 0.2,
                                   point_alpha = 0.6,
                                   color_scale = ggplot2::scale_color_viridis_c,
                                   rasterise = TRUE,
                                   verbose = TRUE) {
  
  # Load required libraries
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("org.Hs.eg.db package is required")
  }
  
  # Convert gene symbols to Ensembl IDs if needed
  ensembl_ids <- tryCatch({
    AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, 
                         keys = markers, 
                         column = "ENSEMBL", 
                         keytype = "SYMBOL", 
                         multiVals = "first")
  }, error = function(e) {
    if (verbose) cat("Warning: Could not map gene symbols to Ensembl IDs\n")
    return(NULL)
  })
  
  # Check which markers are available in the data
  available_markers <- intersect(c(markers, ensembl_ids), rownames(sce))
  
  if (verbose) {
    cat("Available markers:", paste(available_markers, collapse = ", "), "\n")
    if (length(available_markers) == 0) {
      cat("No markers found in the dataset. Available genes start with:", 
          paste(head(rownames(sce), 5), collapse = ", "), "\n")
    }
  }
  
  # Create expression plots for available markers
  if (length(available_markers) > 0) {
    marker_plots <- list()
    
    for (marker in available_markers) {
      # Get original symbol name for title if this is an Ensembl ID
      marker_title <- marker
      if (marker %in% ensembl_ids) {
        symbol_match <- names(ensembl_ids)[ensembl_ids == marker & !is.na(ensembl_ids)]
        if (length(symbol_match) > 0) {
          marker_title <- paste0(symbol_match[1], " (", marker, ")")
        }
      }
      
      marker_plots[[marker]] <- scater::plotReducedDim(
        sce,
        dimred = dimred,
        colour_by = marker,
        point_size = point_size,
        point_alpha = point_alpha,
        rasterise = rasterise
      ) +
      color_scale(name = "Expression") +
      ggplot2::labs(title = marker_title) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        legend.position = "bottom",
        plot.title = ggplot2::element_text(size = 10, hjust = 0.5)
      )
    }
    
    # Combine all marker plots
    combined_plot <- patchwork::wrap_plots(marker_plots, ncol = ncol)
    
    return(combined_plot)
    
  } else {
    if (verbose) cat("No markers found in the dataset\n")
    return(NULL)
  }
}

# Plot fibroblast markers using the function 
fibroblast_marker_plot <- plot_marker_expression(
  sce = results$macro_muscle,
  #markers = c("PECAM1", "VWF", "KDR", "FLT1", "CLDN5", "SPARCL1"),  # endothelial
  #markers = c("RGS5","PDGFRB","CSPG4","MCAM","ABCC9","KCNJ8","NOTCH3"),  # pericytes
  markers = c("COL1A1","COL1A2","COL3A1","DCN","LUM","PDGFRA","THY1","FN1","POSTN"),  # fibroblasts
  dimred = "UMAP_HARMONY",
  ncol = 3,
  verbose = TRUE
)

if (!is.null(fibroblast_marker_plot)) {
  print(fibroblast_marker_plot)
}



  # Plot by marker-based labels if available
  if ("cluster_harmony_label_markers" %in% names(SummarizedExperiment::colData(results$macro_muscle))) {
    p2 <- scater::plotReducedDim(
      results$macro_muscle, 
      dimred = "UMAP_HARMONY", 
      colour_by = "cluster_harmony_label_markers",
      point_size = 0.3,
      point_alpha = 0.7,
      rasterise = TRUE
    ) + 
    ggplot2::labs(title = "Marker-based Labels (Harmony UMAP)", color = "Cell Type") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
  }
  
  # Plot by SingleR labels if available
  if ("cluster_harmony_label_singleR" %in% names(SummarizedExperiment::colData(results$macro_muscle))) {
    p3 <- scater::plotReducedDim(
      results$macro_muscle, 
      dimred = "UMAP_HARMONY", 
      colour_by = "cluster_harmony_label_singleR",
      point_size = 0.3,
      point_alpha = 0.7,
      rasterise = TRUE
    ) + 
    ggplot2::labs(title = "SingleR Labels (Harmony UMAP)", color = "Cell Type") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
  }
}

# 2. Cluster composition barplot
if (safe_require("ggplot2")) {
  p4 <- ggplot2::ggplot(cluster_summary, ggplot2::aes(x = reorder(cluster, -n_cells), y = n_cells)) +
    ggplot2::geom_col(fill = "steelblue", alpha = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = n_cells), vjust = -0.5, size = 3) +
    ggplot2::labs(x = "Cluster", y = "Number of Cells", title = "Cluster Sizes") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

# 3. Marker expression heatmap for top markers
if (safe_require("ComplexHeatmap") && safe_require("circlize")) {
  # Get top 5 markers per cluster
  top5_markers <- do.call(c, lapply(split(results$markers_df, results$markers_df$cluster), function(x) head(x$gene, 5)))
  top5_markers <- unique(top5_markers)
  
  if (length(top5_markers) > 0 && all(top5_markers %in% rownames(results$macro_muscle))) {
    # Calculate cluster means for heatmap
    cluster_means <- sapply(levels(factor(results$macro_muscle$cluster_harmony)), function(cl) {
      cells <- which(results$macro_muscle$cluster_harmony == cl)
      rowMeans(SummarizedExperiment::assay(results$macro_muscle, "logcounts")[top5_markers, cells, drop = FALSE])
    })
    
    # Create heatmap
    col_fun <- circlize::colorRamp2(c(0, max(cluster_means)/2, max(cluster_means)), c("blue", "white", "red"))
    
    h1 <- ComplexHeatmap::Heatmap(
      cluster_means,
      name = "Mean Expression",
      col = col_fun,
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      show_row_names = TRUE,
      show_column_names = TRUE,
      row_names_gp = grid::gpar(fontsize = 8),
      column_names_gp = grid::gpar(fontsize = 10),
      heatmap_legend_param = list(title = "Log Expression")
    )
  }
}

message_time("Analysis complete. Objects created:")
message_time("- results: saved to macro_muscle_communication_5_tissues_hdf5_annotated.rds")
message_time("- cluster_summary: detailed cluster annotations")
message_time("- p1, p2, p3: UMAP plots") 
message_time("- p4: cluster size barplot")
message_time("- h1: marker expression heatmap")
message_time("- results$macro_muscle: annotated SCE with new colData columns")