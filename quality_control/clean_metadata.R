library(cellNexus)
library(dplyr)
library(duckdb)
library(dbplyr)
library(tibble)
library(here)
library(tidyr)
source(here("quality_control/plot_custom_theme.R"))
source(here("quality_control/CAQ_age_analysis_functions.R"))
source(here("quality_control/utils.R"))

# Manipulate metadata
cell_metadata = get_metadata()
census_metadata <- cellNexus:::get_census_metadata("2024-07-01")
con <- dbplyr::remote_con(cell_metadata)
duckdb::duckdb_register_arrow(con, "census_metadata", census_metadata)
cell_metadata <- cell_metadata |>
  dplyr::left_join(tbl(con, "census_metadata")) |>
  # This dataset is removed due to high rates of low quality cells
  filter(assay != "ScaleBio single cell RNA sequencing")

age_groups_tbl <- cell_metadata |>
  distinct(sample_id, age_days, sex) |>
  mutate(sex = ifelse(is.na(sex), "unknown", sex)) |>
  as_tibble() |>
  mutate(
    age_groups = coarse_age_bin(age_days, sex),
    age_groups_fine = age_bin(age_days, sex)
  )

tissue_group_lookup <- get_tissue_grouped()
tissue_group_default <- tissue_group_lookup |>
  filter(is.na(sex)) |>
  select(-sex) |>
  rename(tissue_groups_default = tissue_groups)
tissue_group_by_sex <- tissue_group_lookup |>
  filter(!is.na(sex)) |>
  rename(tissue_groups_by_sex = tissue_groups)

cell_metadata <- cell_metadata |>
  left_join(tissue_group_default, by = "tissue", copy = TRUE) |>
  left_join(tissue_group_by_sex, by = c("tissue", "sex"), copy = TRUE) |>
  mutate(
    tissue_groups = coalesce(tissue_groups_by_sex, tissue_groups_default),
    tissue_groups = ifelse(tissue %in% c("nose skin", "scalp"),
      "integumentary system (skin)",
      tissue_groups
    ),
    age_years = age_days / 365.25,
    sex = ifelse(is.na(sex), "unknown", sex),
    logNorm_nCount_RNA = log(nCount_RNA)
  ) |>
  select(-tissue_groups_default, -tissue_groups_by_sex) |>
  left_join(ethnicity_grouped, copy = T) |>
  left_join(assay_data_grouped, copy = T) |>
  left_join(disease_data_grouped, copy = T) |>
  left_join(disease_data_grouped_coarse, copy = T) |>
  left_join(age_groups_tbl, copy = T) |>
  left_join(shorten_technology_label, copy = T)

cell_metadata <- cell_metadata |>
  left_join(tissue_group_conversion_tbl, by = "tissue_groups", copy = TRUE) |>
  rename(tissue_groups_shorten = tissue_groups_short)
