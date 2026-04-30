# Disease - organ - technology alluvial plot
library(ggalluvial)
library(stringr)
library(ggfittext)
library(dplyr)
library(purrr)
library(targets)
library(ggplot2)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(hrbrthemes)
library(viridis)
library(packcircles)
library(ggrepel)
library(patchwork)
library(cellNexus)
source("plot_custom_theme.R")
source("clean_metadata.R")

tissue_cols <- c(
  "Blood" = "#91AD9E",
  "Respiratory" = "#979771",
  "Cerebral Lobes" = "#B0B1B6",
  "Breast" = "#C2CEDC",
  "Skin" = "#D4BAAD",
  "Renal" = "#E1CCB1"
)

cell_metadata <- cell_metadata |> filter(assay_groups != "ScaleBio")

### Panel A: Alluvial plot of tissue, disease, assay groups
# Do it separately, tissue groups level

# Extract top tissues donoted by donors
top_tissues_in_donors <- cell_metadata |>
  group_by(tissue_groups) |>
  summarise(donor_count = n_distinct(donor_id), .groups = "drop") |>
  collect() |>
  arrange(desc(donor_count)) |>
  head() |>
  pull(tissue_groups)

# Disease groups has little donors
diseases_to_exclude <- cell_metadata |>
  group_by(disease_groups_coarse) |>
  summarise(n = n_distinct(donor_id)) |>
  arrange(n) |>
  collect() |>
  head(3) |>
  pull(disease_groups_coarse)

# Get tissue donor counts
tissue_disease_pairs <- cell_metadata |>
  filter(tissue_groups %in% top_tissues_in_donors) |>
  # Exclude disease groups has little donors
  filter(!disease_groups_coarse %in% diseases_to_exclude) |>
  group_by(disease_groups_coarse, tissue_groups, assay_groups) |>
  summarise(donor_count = n_distinct(donor_id), .groups = "drop") |>
  collect() |>
  mutate(
    disease_groups_coarse = str_wrap(disease_groups_coarse, width = 20),
    tissue_groups = str_wrap(tissue_groups, width = 20)
    # assay_groups = str_wrap(assay_groups, width = 20)
  ) |>
  mutate(
    disease_groups_coarse = as.factor(disease_groups_coarse),
    tissue_groups = as.factor(tissue_groups),
    assay_groups = as.factor(assay_groups)
  )

# Extract order for three levels
tissue_order <- tissue_disease_pairs |>
  group_by(tissue_groups) |>
  summarise(total = sum(donor_count)) |>
  arrange(desc(total)) |>
  pull(tissue_groups)

tech_order <- tissue_disease_pairs |>
  group_by(assay_groups) |>
  summarise(total = sum(donor_count)) |>
  arrange(desc(total)) |>
  pull(assay_groups)

disease_order <- tissue_disease_pairs |>
  group_by(disease_groups_coarse) |>
  summarise(total = sum(donor_count)) |>
  arrange(desc(total)) |>
  pull(disease_groups_coarse)

# Refactor
tissue_disease_pairs <- tissue_disease_pairs %>%
  mutate(
    tissue_groups = factor(tissue_groups, levels = tissue_order),
    assay_groups = factor(assay_groups, levels = tech_order),
    disease_groups_coarse = factor(disease_groups_coarse, levels = disease_order)
  ) |>
  arrange(tissue_groups)

# Assign different text size to alluvial node
tissue_disease_pairs <- tissue_disease_pairs |> mutate(
  tissue_text_size = case_when(
    tissue_groups == "Blood" ~ 4,
    tissue_groups == "Respiratory" ~ 3,
    tissue_groups == "Cerebral Lobes" ~ 2.5,
    tissue_groups == "Breast" ~ 2,
    tissue_groups == "Skin" ~ 2,
    tissue_groups == "Renal" ~ 1.5
  ),
  assay_text_size = case_when(
    assay_groups == "10x Genomics 3" ~ 2.5,
    assay_groups == "10x Genomics 5" ~ 2.5,
    assay_groups == "Plate based Technologies" ~ 1.5,
    assay_groups == "Other Technologies" ~ 1.5,
    assay_groups == "Smart seq" ~ 1.5,
    assay_groups == "Microwell Technologies" ~ 1,
    assay_groups == "Drop based Technologies" ~ 1
  ),
  disease_text_size = case_when(
    disease_groups_coarse == "Normal" ~ 3,
    disease_groups_coarse == "COVID-19" ~ 3,
    disease_groups_coarse == "Glioblastoma" ~ 2,
    disease_groups_coarse == "Lung AdenoCA" ~ 2,
    disease_groups_coarse == "Cancer" ~ 2,
    disease_groups_coarse == "Lupus" ~ 2,
    disease_groups_coarse == "Respiratory" ~ 1.5,
    disease_groups_coarse == "Metabolic/Other" ~ 1.5,
    disease_groups_coarse == "Infectious/Immune" ~ 1
  )
)


tissue_size <- tissue_disease_pairs |>
  distinct(tissue_groups, tissue_text_size) |>
  deframe()
disease_size <- tissue_disease_pairs |>
  distinct(disease_groups_coarse, disease_text_size) |>
  deframe()
assay_size <- tissue_disease_pairs |>
  distinct(assay_groups, assay_text_size) |>
  deframe()
tissue_disease_size <- c(tissue_size, disease_size)
tissue_assay_size <- c(tissue_size, assay_size)


# Plotting
tissue_disease_alluvial <- tissue_disease_pairs |>
  ggplot(aes(
    axis1 = tissue_groups, axis2 = disease_groups_coarse,
    y = donor_count
  )) +
  geom_alluvium(aes(fill = tissue_groups), width = 0.5, show.legend = F) +
  geom_stratum(aes(fill = tissue_groups), width = 0.5, color = "black", stat = "stratum", na.rm = TRUE, show.legend = F) +
  geom_text(
    stat = "stratum",
    aes(
      label = after_stat(stratum),
      size  = after_stat(stratum)
    ),
    show.legend = F
  ) +
  scale_size_manual(values = tissue_disease_size) + # use custom sizes
  labs(x = NULL, y = NULL, title = NULL) +
  scale_x_discrete(limits = c("Tissue Groups", "Disease Groups"), expand = c(0.05, 0.05)) +
  scale_y_continuous(breaks = seq(0, 5471, by = 500), expand = c(0, 0)) +
  coord_cartesian(clip = "off") + # allow text fits outside of the node
  theme_multipanel +
  scale_fill_manual(values = tissue_cols, na.value = NA) +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "grey95", color = NA)
  )
# guides(fill = "none", color = "none")

tissue_assay_alluvial <- tissue_disease_pairs |>
  ggplot(aes(
    axis1 = assay_groups, axis2 = tissue_groups,
    y = donor_count
  )) +
  geom_alluvium(aes(fill = tissue_groups), width = 0.5, show.legend = F) +
  geom_stratum(aes(fill = tissue_groups), width = 0.5, color = "black", na.rm = TRUE, show.legend = F) +
  geom_text(
    stat = "stratum",
    aes(
      label = after_stat(stratum),
      size  = after_stat(stratum)
    ),
    show.legend = F
  ) +
  scale_size_manual(values = tissue_assay_size) + # use custom sizes
  labs(x = "", y = expression("Donor Count"), "", title = "") +
  scale_x_discrete(limits = c("Technologies", "Tissue Groups"), expand = c(0.05, 0.05)) +
  scale_y_continuous(breaks = seq(0, 5471, by = 500), expand = c(0, 0)) +
  coord_cartesian(clip = "off") + # allow text fits outside of the node
  theme_multipanel +
  scale_fill_manual(values = tissue_cols, na.value = NA) +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "grey95", color = NA)
  )
# guides(fill = "none", color = "none")
