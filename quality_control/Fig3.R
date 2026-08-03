# Standard: 183mm width x max 230mm height
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
library(ggalluvial)
library(stringr)
library(ggfittext)
library(zellkonverter)
library(grid)
library(duckdb)

setwd("quality_control/")
source("plot_custom_theme.R")
source("CAQ_age_analysis_functions.R")
source("clean_metadata.R")

### Figure 3B
# Hard to track samples density ---------------------------------------------------
samples <- c("734b3e20708fa099e9af65400d3d30ba", "3d459a68fe217c54e4180fe60b7a756e", "dc34c7f382dd3e246b1d39ee019e6fe3")

get_raw_counts_tbl_per_sample <- function(sample) {
  raw_counts <- zellkonverter::readH5AD(paste0(sample, ".h5ad"),
    reader = "R",
    use_hdf5 = TRUE
  ) |>
    assay("X") |>
    colSums() |>
    as.numeric()
  df <- tibble(
    counts = raw_counts,
    sample_id = sample
  )
  df |> mutate(method = "raw")
}

get_transformed_counts_tbl_per_sample <- function(sample) {
  transformed_counts <- cell_metadata |>
    filter(sample_id == sample) |>
    get_single_cell_experiment() |>
    counts() |>
    colSums() |>
    as.numeric()

  df <- tibble(counts = transformed_counts, sample_id = sample) |>
    mutate(method = "transformed")
  df
}

raw_counts_tbl <- map(samples, ~ get_raw_counts_tbl_per_sample(.x)) |>
  bind_rows()

transformed_counts_tbl <- map(samples, ~ get_transformed_counts_tbl_per_sample(.x)) |>
  bind_rows() |>
  left_join(cell_metadata |>
              distinct(sample_id, inverse_transform),
            by = "sample_id",
            copy = T)


p1 <- raw_counts_tbl |> ggplot(aes(x = counts, group = sample_id)) +
  geom_density(alpha = 0.4) +
  labs(
    x = "Counts",
    y = "Density"
  ) +
  theme_multipanel +
  scale_fill_manual(values = friendly_cols) +
  scale_color_manual(values = friendly_cols) +
  theme( # legend.position = "none",
    axis.text.x = element_text(size = 5, hjust = 1)
  )

p2 <- transformed_counts_tbl |> ggplot(aes(x = counts, group = sample_id, linetype = inverse_transform)) +
  geom_density(alpha = 0.4) +
  labs(
    x = "Transformed Counts",
    y = "Density",
    linetype = "Method"
  ) +
  theme_multipanel +
  # scale_fill_manual(values = friendly_cols) +
  # scale_color_manual(values = friendly_cols) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 5, hjust = 1),
    axis.title.x = element_text(hjust = 0.5)
  )

p1 + p2

### Figure 3C: Cell quality percentage by tech
cell_quality_cols <- c(
  "Empty Droplet" = "#4C72B0",
  "Dead"     = "#DD8452",
  "Doublet"       = "#55A868"
)

# Combine bad quality cells percentagte by tech
per_sample_rescaled <- cell_metadata |>
  group_by(sample_id, assay_shorten) |>
  summarise(
    total_n   = n(),
    empty_n   = sum(empty_droplet, na.rm = TRUE),
    dead_n    = sum(!empty_droplet & !alive, na.rm = TRUE),
    doublet_n = sum(!empty_droplet & scDblFinder.class == "doublet", na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    non_empty_n  = total_n - empty_n,
    frac_empty   = empty_n / total_n,
    frac_dead    = dead_n / total_n,
    frac_doublet = doublet_n / total_n
  )

tech_summary <- per_sample_rescaled |>
  group_by(assay_shorten) |>
  summarise(
    mean_frac_empty   = mean(frac_empty),
    mean_frac_dead    = mean(frac_dead),
    mean_frac_doublet = mean(frac_doublet),
    n_samples         = n(),
    .groups = "drop"
  ) |>
  mutate(
    mean_frac_good = 1 - mean_frac_empty - mean_frac_dead - mean_frac_doublet
  ) |>
  collect()

tech_long <- tech_summary |>
  select(assay_shorten, mean_frac_good, mean_frac_doublet, mean_frac_dead, mean_frac_empty) |>
  pivot_longer(-assay_shorten, names_to = "category", values_to = "frac") |>
  mutate(
    pct = frac * 100,
    category = recode(category,
                      mean_frac_good    = "Alive singlet",
                      mean_frac_doublet = "Doublet",
                      mean_frac_dead    = "Dead cell",
                      mean_frac_empty   = "Empty droplet"
    ),
    category = factor(category, levels = c("Alive singlet", "Doublet", "Dead cell", "Empty droplet"))
  )

assay_order <- tech_long |>
  filter(category == "Alive singlet") |>
  arrange(pct) |>
  pull(assay_shorten)

tech_long <- tech_long |>
  mutate(assay_shorten = factor(assay_shorten, levels = assay_order)) |>
  filter(category != "Alive singlet")

ggplot(tech_long, aes(x = assay_shorten, y = pct, fill = category)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c(
    "Doublet"       = "#F2C14E",
    "Dead cell"     = "#F26B5B",
    "Empty droplet" = "#495867"
  )) +
  labs(
    #title = "True composition of all droplets, by assay",
    #subtitle = "Per-sample rescaled fractions averaged across samples",
    x = NULL, y = "Mean percentage of all droplets across samples", fill = NULL
  ) +
  theme_multipanel +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "top"
  )

ggsave(here("quality_control", "Figure3C.pdf"), width = 120, height = 70, plot = last_plot(), units = "mm")


# Figure 3D: Doublet percentage plot
source("clean_metadata.R")
df <- cell_metadata |>
  filter(!is.na(cell_type_unified_ensemble)) |>
  group_by(sample_id, cell_type_unified_ensemble) |>
  summarise(
    total_cells = n(),
    doublet_cells = sum(scDblFinder.class == "doublet"),
    .groups = "drop"
  ) |>
  mutate(doublet_pct = 100 * doublet_cells / total_cells) |>
  group_by(cell_type_unified_ensemble) |>
  summarise(
    grand_avg_doublet_pct = mean(doublet_pct, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  ) |>
  filter(grand_avg_doublet_pct > 3)

doublet_pct_p <- df |>
  ggplot(aes(
    x = reorder(cell_type_unified_ensemble, -grand_avg_doublet_pct),
    y = grand_avg_doublet_pct
  )) +
  geom_bar(position = "stack", stat = "identity") +
  labs(
    x = "Cell Type",
    y = "Avg Doublet % across Samples"
  ) +
  theme_multipanel +
  theme(axis.text.x = element_text(size = 5, angle = 30, hjust = 1))

doublet_pct_p
