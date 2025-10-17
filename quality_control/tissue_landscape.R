library(cellNexus)
library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)

# Import theme_multipanel
source("https://gist.githubusercontent.com/stemangiola/fc67b08101df7d550683a5100106561c/raw/a0853a1a4e8a46baf33bad6268b09001d49faf51/ggplot_theme_multipanel")

my_metadata =
    get_metadata() |>
    left_join(readr::read_csv("tissue_groups_to_labels.csv"), copy = TRUE) |>
    filter(!is.na(tissue_label))


empty_droplet_proportion = 
    my_metadata |> 
    count(empty_droplet, tissue_label, sample_id, dataset_id) |>
    mutate(empty_proportion = n/sum(n)) |> 
    summarise(mean_empty_proportion = mean(empty_proportion, na.rm = TRUE), .by=tissue_label)

dead_proportion = 
    my_metadata |> 
    count(alive, tissue_label, sample_id, dataset_id) |> 
    mutate(alive = if_else(is.na(alive), FALSE, alive)) |> 
    mutate(dead_proportion = n/sum(n)) |> 
    summarise(mean_dead_proportion = mean(dead_proportion, na.rm = TRUE), .by=tissue_label)

doublet_proportion = 
    my_metadata |> 
    count(scDblFinder.class, tissue_label, sample_id, dataset_id) |> 
    mutate(doublet_proportion = n/sum(n)) |> 
    summarise(mean_doublet_proportion = mean(doublet_proportion, na.rm = TRUE), .by=tissue_label)


feature_count = 
    my_metadata |> 
    summarise(mean_feature_count = mean(feature_count, na.rm = TRUE), .by=c(tissue_label, sample_id)) |> 
    summarise(mean_feature_count = mean(mean_feature_count, na.rm = TRUE), .by=tissue_label)

donor_size = 
    my_metadata |> 
    distinct(tissue_label, donor_id) |> 
    count(tissue_label) |> 
    rename(donor_count = n) 

dataset_size = 
    my_metadata |> 
    distinct(tissue_label, dataset_id) |> 
    count(tissue_label) |> 
    rename(dataset_size = n) 

cell_count_mean = 
    my_metadata |> 
    count(tissue_label, sample_id, dataset_id) |> 
    summarise(mean_cell_count = mean(n, na.rm = TRUE), .by=tissue_label)

cell_count = 
    my_metadata |> 
    count(tissue_label)

missing_annotation_ethnicity_proportion = 
    my_metadata |> 
    mutate(is_missing_annotation_ethnicity = self_reported_ethnicity == "unknown") |> 
    distinct(donor_id, tissue_label, is_missing_annotation_ethnicity) |> 
    count(tissue_label, is_missing_annotation_ethnicity) |> 
    pivot_wider(names_from = is_missing_annotation_ethnicity, values_from = n) |> 
    mutate(proportion_missing_annotation_ethnicity = `TRUE`/sum(`TRUE`, `FALSE`)) |> 
    mutate(proportion_missing_annotation_ethnicity = if_else(is.na(proportion_missing_annotation_ethnicity), 0, proportion_missing_annotation_ethnicity)) |> 
    select(tissue_label, proportion_missing_annotation_ethnicity)

missing_annotation_sex_proportion = 
    my_metadata |> 
  filter(age_days >= 365) |> 
  filter(!grepl("fetal|embryo|placenta|cord blood", tolower(title), ignore.case = TRUE)) |> 
  mutate(is_missing_annotation_sex = sex == "unknown" | is.na(sex)) |> 
  distinct(donor_id, tissue_label, is_missing_annotation_sex) |> 
  count(tissue_label, is_missing_annotation_sex) |> 
  pivot_wider(names_from = is_missing_annotation_sex, values_from = n) |> 
  mutate(proportion_missing_annotation_sex = `TRUE`/sum(`TRUE`, `FALSE`)) |> 
  mutate(proportion_missing_annotation_sex = if_else(is.na(proportion_missing_annotation_sex), 0, proportion_missing_annotation_sex)) |> 
    select(tissue_label, proportion_missing_annotation_sex) 


qc_table = 
    empty_droplet_proportion |> 
    left_join(dead_proportion) |> 
    left_join(doublet_proportion) |> 
    left_join(feature_count) |> 
    left_join(donor_size) |> 
    left_join(cell_count_mean) |> 
    left_join(cell_count) |> 
    left_join(missing_annotation_ethnicity_proportion) |> 
    left_join(missing_annotation_sex_proportion) |>
    left_join(dataset_size) |>
    filter(!is.na(tissue_label))

# Reshape data to long format for visualization

# Which columns have NA?
 qc_table_collected = qc_table |> collect()
na_columns = colnames(qc_table_collected)[colSums(is.na(qc_table_collected)) > 0]
print(na_columns)

# source("tissue_color_utils.R")



qc_table_long = 
  qc_table |> 
  select(tissue_label, mean_empty_proportion, mean_dead_proportion, mean_doublet_proportion, mean_feature_count, 
         donor_count, mean_cell_count, n, 
         proportion_missing_annotation_ethnicity, proportion_missing_annotation_sex, dataset_size) |> 
  
  # Drop NA tissue_label
  as_tibble() |> 
  pivot_longer(cols = -tissue_label, names_to = "metric", values_to = "value") |>
 # mutate(tissue_label = get_tissue_abbrev(tissue_label)) |>

  # abbreviate metrics names smartly, e.g. keeping empty, alive, doublet, feature_count, donor_count, cell_count, missing_annotation_ethnicity, missing_annotation_sex
  mutate(metric_label = case_when(
    metric == "mean_empty_proportion" ~ "Empty",
    metric == "mean_dead_proportion" ~ "Alive",
    metric == "mean_doublet_proportion" ~ "Doublet",
    metric == "mean_feature_count" ~ "Features",
    metric == "donor_count" ~ "Donors",
    metric == "mean_cell_count" ~ "Cells mean",
    metric == "proportion_missing_annotation_ethnicity" ~ "Ethnicity",
    metric == "proportion_missing_annotation_sex" ~ "Sex",
    metric == "n" ~ "Cells",
    metric == "dataset_size" ~ "Datasets",
  )) |>
      mutate(metric_label = metric_label |> fct_relevel("Cells mean", "Donors", "Datasets", "Features", "Alive", "Empty", "Doublet", "Ethnicity", "Sex")) 

# mutate value to 0 if na
qc_table_long = qc_table_long |> mutate(value = if_else(is.na(value), 0, value))

# Create bar plots with metrics as columns and tissues as rows
plot_barplots = qc_table_long |> 

    filter(metric_label != "Cells") |>
    left_join(readr::read_csv("tissue_labels_to_macrocategories.csv"), copy = TRUE) |>
    mutate(macrocategory = if_else(is.na(macrocategory), "Other", macrocategory)) |>
  ggplot(aes(x = value, y = tissue_label)) +
  geom_col() +
  facet_grid(macrocategory ~ metric_label, scales = "free", space = "free_y") +
  theme_multipanel +
  theme(
    strip.text = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  ) +
  labs(x = "Value") 

# n should be a bubble plot that you patchwork with the other plots
plot_n = 
  qc_table_long |> 
  left_join(readr::read_csv("tissue_labels_to_macrocategories.csv"), copy = TRUE) |>
  mutate(macrocategory = if_else(is.na(macrocategory), "Other", macrocategory)) |>
  filter(metric == "n") |> 
  ggplot(aes(x = 1, y = tissue_label)) +
  geom_point(aes(size = value)) +
  facet_grid(macrocategory ~ metric_label, scales = "free", space = "free_y") +
  theme_multipanel +
  theme(
    strip.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  ) +
  labs(y = "Tissue") 

# patchwork the plots
library(patchwork)
plot_n + plot_barplots + plot_annotation(tag_levels = 'A') + plot_layout(guides = 'collect', nrow = 1, widths = c(1, 8)) &
  theme(plot.margin = margin(0, 0, 0, 0, "pt"), legend.position = "bottom")

# save within directory
ggsave("quality_control/tissue_landscape.pdf", width = 14, height = 10, plot = last_plot())
