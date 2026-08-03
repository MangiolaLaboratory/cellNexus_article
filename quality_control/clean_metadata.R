library(cellNexus)
library(dplyr)
library(duckdb)
library(dbplyr)
library(tibble)
source("plot_custom_theme.R")
source("CAQ_age_analysis_functions.R")

# Manipulate metadata
cell_metadata = get_metadata()
census_metadata <- cellNexus:::get_census_metadata("2024-07-01")
con <- dbplyr::remote_con(cell_metadata)
duckdb::duckdb_register_arrow(con, "census_metadata", census_metadata)
cell_metadata <- cell_metadata |>
  dplyr::left_join(tbl(con, "census_metadata")) |>
  # This dataset is removed due to high rates of low quality cells
  filter(assay != "ScaleBio single cell RNA sequencing")

nfeatures_df <- cellNexus:::get_cellxgene_metadata("dataset") |>
  dplyr::select(dplyr::where(~ !is.list(.x)))

cell_metadata <- cell_metadata |>
  dplyr::left_join(nfeatures_df,
                   by = "dataset_id",
                   copy = TRUE) 

age_groups_tbl <- cell_metadata |>
  distinct(sample_id, age_days, sex) |>
  mutate(sex = ifelse(is.na(sex), "unknown", sex)) |>
  as_tibble() |>
  mutate(
    age_groups = coarse_age_bin(age_days, sex),
    age_groups_fine = age_bin(age_days, sex)
  )

cell_metadata <- cell_metadata |>
  mutate(
    tissue_groups = ifelse(tissue %in% c("nose skin", "scalp"),
      "integumentary system (skin)",
      tissue_groups
    ),
    age_years = age_days / 365.25,
    sex = ifelse(is.na(sex), "unknown", sex),
    scaled_nCount_RNA = log(nCount_RNA)
  ) |>
  left_join(ethnicity_grouped, copy = T) |>
  left_join(assay_data_grouped, copy = T) |>
  left_join(disease_data_grouped, copy = T) |>
  left_join(disease_data_grouped_coarse, copy = T) |>
  left_join(age_groups_tbl, copy = T) |>
  left_join(shorten_technology_label, copy = T)

tissue_group_conversion_tbl <- tibble::tibble(
  tissue_groups = c(
    "small intestine",
    "spleen",
    "large intestine",
    "sensory-related structures",
    "adipose tissue",
    "vasculature",
    "male reproductive system (other)",
    "umbilical cord blood",
    "stomach",
    "gastrointestinal accessory organs",
    "nasal, oral, and pharyngeal regions",
    "connective tissue",
    "muscular system (skeletal muscles)",
    "breast",
    "prostate",
    "digestive tract junctions and connections",
    "respiratory system",
    "lymphatic system",
    "epithelium and mucosal tissues",
    "oesophagus",
    "trachea",
    "ovary",
    "peritoneal and abdominal cavity structures",
    "blood",
    "brainstem and cerebellar structures",
    "integumentary system (skin)",
    "miscellaneous glands",
    "thymus",
    "limbic and basal systems",
    "endocrine system",
    "renal system",
    "general brain and major structures",
    "cerebral lobes and cortical areas",
    "bone marrow",
    "female reproductive system",
    "cardiovascular system",
    "digestive system (general)"
  ),
  tissue_groups_short = c(
    "Small Intestine",
    "Spleen",
    "Large Intestine",
    "Sensory",
    "Adipose",
    "Vasculature",
    "Male Repro.",
    "Cord Blood",
    "Stomach",
    "GI Accessory",
    "Nasal/Oral/Pharyngeal",
    "Connective",
    "Muscle",
    "Breast",
    "Prostate",
    "GI Junctions",
    "Respiratory",
    "Lymphatic",
    "Epithelium/Mucosa",
    "Oesophagus",
    "Trachea",
    "Ovary",
    "Peritoneal/Abdominal",
    "Blood",
    "Brainstem/Cerebellum",
    "Skin",
    "Misc. Glands",
    "Thymus",
    "Limbic/Basal",
    "Endocrine",
    "Renal",
    "Brain (General)",
    "Cerebral Lobes",
    "Bone Marrow",
    "Female Repro.",
    "Cardiovascular",
    "Digestive"
  )
)

cell_metadata <- cell_metadata |>
  left_join(tissue_group_conversion_tbl, by = "tissue_groups", copy = TRUE) |>
  mutate(tissue_groups = tissue_groups_short) |>
  select(-tissue_groups_short)
