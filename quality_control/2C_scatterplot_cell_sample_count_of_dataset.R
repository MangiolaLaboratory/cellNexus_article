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

tech_cols <- c(
  "10x Genomics 3" = "#a6cee3",
  "Multi-technologies" = "#1f78b4",
  "10x Genomics 5" = "#b2df8a",
  "Other Technologies" = "#33a02c",
  "Smart seq" = "#fb9a99",
  "Microwell Technologies" = "#e31a1c",
  "TruDrop" = "#fdbf6f",
  "Drop based Technologies" = "#ff7f00",
  "Plate based Technologies" = "#cab2d6"
)

# Step 1: Count alive cells per donor
donor_summary <- cell_metadata %>%
  filter(alive == TRUE) %>%
  group_by(dataset_id, donor_id) %>%
  summarise(alive_cells = n(), .groups = "drop")

# Step 2: Aggregate to dataset-level
dataset_summary <- donor_summary %>%
  group_by(dataset_id) %>%
  summarise(
    num_donor = n_distinct(donor_id),
    avg_alive_cells = mean(alive_cells),
    .groups = "drop"
  ) |>
  arrange(avg_alive_cells)

# Make dataset a factor ordered by avg_alive_cells
dataset_summary <- dataset_summary %>%
  as_tibble() %>%
  mutate(dataset_id = factor(dataset_id, levels = dataset_id[order(avg_alive_cells)]))

# Label the extreme case
label_data <- dataset_summary |>
  filter(num_donor == max(num_donor) |
    # avg_alive_cells == max(avg_alive_cells) |
    num_donor == min(num_donor) & avg_alive_cells == min(avg_alive_cells) | dataset_id == "39ed7d98-676d-4b8d-9d0a-0f3b60914ead") |>
  left_join(cell_metadata |> distinct(dataset_id, citation, explorer_url, title), by = "dataset_id", copy = T) |>
  mutate(title = stringr::str_wrap(title, width = 25)) |>
  mutate(author = case_when(
    explorer_url == "https://cellxgene.cziscience.com/e/d567b692-c374-4628-a508-8008f6778f22.cxg/" ~ "Kanemaru et al",
    explorer_url == "https://cellxgene.cziscience.com/e/3faad104-2ab8-4434-816d-474d8d2641db.cxg/" ~ "Yazar et al",
    explorer_url == "https://cellxgene.cziscience.com/e/2c820d53-cbd7-4e0a-be7a-a0ad1989a98f.cxg/" ~ "Scale Biosciences, Inc",
    # explorer_url == "https://cellxgene.cziscience.com/e/f7c1c579-2dc0-47e2-ba19-8165c5a0e353.cxg/" ~ "Cao et al",
    explorer_url == "https://cellxgene.cziscience.com/e/39ed7d98-676d-4b8d-9d0a-0f3b60914ead.cxg/" ~ "Smith, R. S. et al"
  ))

# Dataset and sequence tech
dataset_tech_df <- cell_metadata |>
  distinct(dataset_id, assay_groups) |>
  collect()

label_ids <- unique(label_data$dataset_id)

# Plot with dataset on x-axis instead of num_donor
alive_per_donor_scatter <- dataset_summary |>
  left_join(label_data |> select(dataset_id, author), copy = T) |>
  left_join(dataset_tech_df, copy = T) |>
  group_by(dataset_id) |>
  mutate(
    assay_groups = ifelse(n_distinct(assay_groups) > 1, "Multi-technologies", assay_groups)
  ) |>
  ungroup() |>
  distinct() |>
  mutate(author = ifelse(author |> is.na(), "", author)) |>
  ggplot(aes(x = num_donor, y = avg_alive_cells, group = author, color = assay_groups)) +
  geom_point(size = 0.1) +
  geom_point( # overlay solid black for labelled datasets
    data = \(d) d |> filter(dataset_id %in% label_ids),
    aes(x = num_donor, y = avg_alive_cells),
    inherit.aes = FALSE,
    color = "black", size = 1
  ) +
  scale_y_log10() +
  scale_x_log10() +
  labs(
    x = expression("Donor Count"),
    y = expression("Average Number of Alive Cells per Donor"),
  ) +
  geom_text_repel(
    data = label_data,
    aes(label = author),
    color = "black",
    size = 2,
    box.padding = 0.8,
    segment.color = "grey50",
    max.overlaps = Inf
    # nudge_y = -0.2
  ) +
  scale_color_manual(
    values = tech_cols,
    guide = guide_legend(nrow = 3)
  ) +
  theme_multipanel +
  theme(
    # Adjust legend size
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 5),
    legend.key.height = unit(0.05, "cm"),
    legend.key.width = unit(0.05, "cm")
  )

alive_per_donor_scatter
