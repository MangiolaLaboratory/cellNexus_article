# Tissue Color Utilities
# Common color scale for all tissue coloring across the project

library(RColorBrewer)
library(dplyr)


#' Get comprehensive tissue color mapping
#' 
#' This function provides a consistent color mapping for all tissues across the project.
#' It handles tissue name variations and provides both full names and abbreviations.
#' 
#' @param tissues Character vector of tissue names (optional)
#' @param return_abbrev Whether to return abbreviated tissue names as well
#' @return Named character vector of colors, or list with colors and abbreviations
get_tissue_colors <- function(tissues = NULL, return_abbrev = FALSE) {
  
  # Comprehensive tissue name mapping (full names to abbreviations)
  tissue_abbrev_mapping <- c(
    "cardiovascular system" = "cardiovascular",
    "female reproductive system" = "female repro", 
    "integumentary system (skin)" = "skin",
    "cerebral lobes and cortical areas" = "cerebral cortex",
    "brainstem and cerebellar structures" = "brainstem/cerebellum",
    "limbic and basal systems" = "limbic/basal",
    "muscular system (skeletal muscles)" = "skeletal muscle",
    "respiratory system" = "respiratory",
    "renal system" = "renal",
    "digestive system (general)" = "digestive",
    "large intestine" = "large intestine",
    "small intestine" = "small intestine", 
    "nasal, oral, and pharyngeal regions" = "nasal/oral",
    "digestive tract junctions and connections" = "digestive junctions",
    "peritoneal and abdominal cavity structures" = "peritoneal",
    "miscellaneous glands" = "glands",
    "adipose tissue" = "adipose",
    "connective tissue" = "connective tissue",
    "endothelial" = "endothelial",
    "epithelial" = "epithelial",
    "lymphatic system" = "lymphatic",
    "endocrine system" = "endocrine",
    "vasculature" = "vasculature",
    "sensory-related structures" = "sensory",
    "trachea" = "trachea",
    "breast" = "breast",
    "prostate" = "prostate",
    "ovary" = "ovary",
    "thymus" = "thymus",
    "spleen" = "spleen",
    "liver" = "liver",
    "stomach" = "stomach",
    "oesophagus" = "oesophagus",
    "gallbladder" = "gallbladder",
    "pancreas" = "pancreas",
    "bone marrow" = "bone marrow",
    "blood" = "blood",
    "population" = "population"
  )
  
  # Get all unique tissues from the mapping
  all_tissues <- names(tissue_abbrev_mapping)
  
  # If specific tissues provided, filter to those
  if (!is.null(tissues)) {
    # Handle both full names and abbreviations
    matched_tissues <- c()
    for (tissue in tissues) {
      if (tissue %in% all_tissues) {
        matched_tissues <- c(matched_tissues, tissue)
      } else if (tissue %in% tissue_abbrev_mapping) {
        # Find full name from abbreviation
        full_name <- names(tissue_abbrev_mapping)[tissue_abbrev_mapping == tissue]
        matched_tissues <- c(matched_tissues, full_name)
      }
    }
    all_tissues <- unique(matched_tissues)
  }
  
  # Create comprehensive color palette
  # Combine multiple RColorBrewer palettes for maximum distinctiveness
  tissue_colors <- c(
    brewer.pal(8, "Dark2"),    # 8 colors
    brewer.pal(8, "Set2"),     # 8 colors  
    brewer.pal(8, "Set1"),     # 8 colors
    brewer.pal(8, "Accent"),   # 8 colors
    brewer.pal(8, "Paired"),   # 8 colors
    brewer.pal(8, "Set3"),     # 8 colors
    brewer.pal(8, "Pastel1"),  # 8 colors
    brewer.pal(8, "Pastel2")   # 8 colors
  )
  
  # Remove any duplicates and ensure we have enough colors
  tissue_colors <- unique(tissue_colors)
  
  # Create color mapping for tissues
  n_tissues <- length(all_tissues)
  if (n_tissues > length(tissue_colors)) {
    # If we need more colors, extend the palette
    tissue_colors <- rep(tissue_colors, length.out = n_tissues)
  }
  
  # Create named vector
  tissue_color_map <- setNames(tissue_colors[1:n_tissues], all_tissues)
  
  if (return_abbrev) {
    return(list(
      colors = tissue_color_map,
      abbreviations = tissue_abbrev_mapping[names(tissue_color_map)]
    ))
  } else {
    return(tissue_color_map)
  }
}

#' Get tissue color scale for ggplot
#' 
#' @param tissues Character vector of tissue names (optional)
#' @param type Either "color" or "fill" for scale type
#' @param drop Whether to drop unused levels
#' @return ggplot2 scale function
get_tissue_scale <- function(tissues = NULL, type = "color", drop = FALSE) {
  color_map <- get_tissue_colors(tissues)
  
  if (type == "color") {
    return(ggplot2::scale_color_manual(values = color_map, drop = drop))
  } else if (type == "fill") {
    return(ggplot2::scale_fill_manual(values = color_map, drop = drop))
  } else {
    stop("type must be 'color' or 'fill'")
  }
}

#' Get tissue abbreviation for display
#' 
#' @param tissue_name Full tissue name (can be a vector)
#' @return Abbreviated tissue name (same length as input)
get_tissue_abbrev <- function(tissue_name) {
  tissue_abbrev_mapping <- c(
    "cardiovascular system" = "cardiovascular",
    "female reproductive system" = "female repro", 
    "integumentary system (skin)" = "skin",
    "cerebral lobes and cortical areas" = "cerebral cortex",
    "brainstem and cerebellar structures" = "brainstem/cerebellum",
    "limbic and basal systems" = "limbic/basal",
    "muscular system (skeletal muscles)" = "skeletal muscle",
    "respiratory system" = "respiratory",
    "renal system" = "renal",
    "digestive system (general)" = "digestive",
    "large intestine" = "large intestine",
    "small intestine" = "small intestine", 
    "nasal, oral, and pharyngeal regions" = "nasal/oral",
    "digestive tract junctions and connections" = "digestive junctions",
    "peritoneal and abdominal cavity structures" = "peritoneal",
    "miscellaneous glands" = "glands",
    "adipose tissue" = "adipose",
    "connective tissue" = "connective tissue",
    "endothelial" = "endothelial",
    "epithelial" = "epithelial",
    "lymphatic system" = "lymphatic",
    "endocrine system" = "endocrine",
    "vasculature" = "vasculature",
    "sensory-related structures" = "sensory",
    "trachea" = "trachea",
    "breast" = "breast",
    "prostate" = "prostate",
    "ovary" = "ovary",
    "thymus" = "thymus",
    "spleen" = "spleen",
    "liver" = "liver",
    "stomach" = "stomach",
    "oesophagus" = "oesophagus",
    "gallbladder" = "gallbladder",
    "pancreas" = "pancreas",
    "bone marrow" = "bone marrow",
    "blood" = "blood",
    "population" = "population"
  )
  
  # Handle both single values and vectors
  result <- tissue_name
  for (i in seq_along(tissue_name)) {
    if (tissue_name[i] %in% names(tissue_abbrev_mapping)) {
      result[i] <- tissue_abbrev_mapping[tissue_name[i]]
    }
  }
  
  return(result)
}

# Example usage:
# tissue_colors <- get_tissue_colors()
# tissue_scale_color <- get_tissue_scale(type = "color")
# tissue_scale_fill <- get_tissue_scale(type = "fill")
