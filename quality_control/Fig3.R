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

source("plot_custom_theme.R")
source("CAQ_age_analysis_functions.R")
source("clean_metadata.R")

### Figure 3B
# Hard to track samples density ---------------------------------------------------
file_path = "~/"
samples = c("734b3e20708fa099e9af65400d3d30ba","3d459a68fe217c54e4180fe60b7a756e", "dc34c7f382dd3e246b1d39ee019e6fe3")

get_raw_counts_tbl_per_sample <- function(sample) {
  raw_counts = zellkonverter::readH5AD(paste0(file_path, sample,".h5ad"),
                                       reader = "R",
                                       use_hdf5 = TRUE) |> assay("X") |> colSums() |> as.numeric()
  df = tibble(counts = raw_counts,
              sample_id = sample)
  df |> mutate(method = "raw")
}

get_transformed_counts_tbl_per_sample <- function(sample) {
  transformed_counts = cell_metadata |> filter(sample_id == sample) |>
    get_single_cell_experiment() |> counts() |> colSums() |> as.numeric()
  
  df = tibble(counts = transformed_counts, sample_id = sample) |> mutate(method = "transformed")
  df
}

sliced_sample_tbl = readRDS("data/sliced_sample_tbl.rds")
raw_counts_tbl <- map(samples , ~get_raw_counts_tbl_per_sample(.x)) |> bind_rows()

transformed_counts_tbl = map(samples, ~get_transformed_counts_tbl_per_sample(.x)) |> bind_rows() |> 
  left_join(sliced_sample_tbl |> select(sample_2, method_to_apply), by = c("sample_id" = "sample_2"))


p1 <- raw_counts_tbl |> ggplot(aes(x = counts, group = sample_id)) +
  geom_density(alpha = 0.4) +
  labs(
    x = "Counts",
    y = "Density"
  ) +
  theme_multipanel +
  scale_fill_manual(values = friendly_cols) +
  scale_color_manual(values = friendly_cols) +
  theme(#legend.position = "none",
    axis.text.x = element_text(size = 5, hjust = 1)) 

p2 <- transformed_counts_tbl |> ggplot(aes(x = counts, group = sample_id, linetype = method_to_apply)) +
  geom_density(alpha = 0.4) +
  labs(
    x = "Transformed Counts",
    y = "Density",
    linetype = "Method"
  ) +
  theme_multipanel +
  #scale_fill_manual(values = friendly_cols) +
  #scale_color_manual(values = friendly_cols) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 5, hjust = 1),
        axis.title.x = element_text(hjust = 0.5))

p1 + p2

### Figure 3C: Cell quality percentage by tech
cell_quality_cols <- c("Dead" = "#C32626",
                       "Doublet" = "#E1B529",
                       "Empty Droplet" = "#0D6EAE")

# Combine bad quality cells percentagte by tech
qc_per_assay_metadata <- cell_metadata |>
  
  # ScaleBio has lots of Doublet
  filter(assay_shorten != "ScaleBio") |>
  
  group_by(assay_shorten) |>
  summarise(
    total_cells = n(),
    empty_cells = sum(ifelse(empty_droplet == "TRUE", 1, 0)),
    alive_cells = sum(ifelse(alive == TRUE, 1, 0)),
    dead_cells = sum(ifelse(empty_droplet == "FALSE" & alive |> is.na(), 1, 0)), # To calculate the true value, need to exclude empty droplet cells first
    doublet = sum(ifelse(empty_droplet == "FALSE" & scDblFinder.class == "doublet", 1, 0)) # To calculate the true value, need to exclude empty droplet cells first
  ) |>
  mutate(
    `Empty Droplet` = round(empty_cells / total_cells * 100, 2),
    `Alive Cells` = round(alive_cells / total_cells * 100, 2),
    `Dead` = round(dead_cells / total_cells * 100, 2),
    `Doublet` = round(doublet / total_cells * 100, 2)
  ) |> 
  collect() |> 
  #mutate(assay_shorten = stringr::str_wrap(assay_shorten, width = 10)) |>
  tidyr::pivot_longer(cols = c("Empty Droplet","Alive Cells","Dead", "Doublet"), 
                      names_to = "Group", values_to = "Percentage")

qc_plot <- qc_per_assay_metadata |>
  filter(Group != "Alive Cells") |> 
  ggplot(aes(x = reorder(assay_shorten, -Percentage), 
             y = Percentage,
             fill = Group)) +
  geom_bar(position="stack", stat="identity") +
  theme_multipanel +
  scale_fill_manual(values = cell_quality_cols) + 
  xlab("Technology") +
  ylab("Cell Percentage") + 
  theme(axis.text.x = element_text(size = 5, angle = 30, hjust = 1))

qc_plot

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
  ggplot(aes(x = reorder(cell_type_unified_ensemble, -grand_avg_doublet_pct), 
             y = grand_avg_doublet_pct)) + 
  geom_bar(position="stack", stat="identity"
           ) +
  labs(
    x = "Cell Type",
    y = "Avg Doublet % across Samples"
  ) +
  theme_multipanel +
  theme(axis.text.x = element_text(size = 5, angle = 30, hjust = 1))

doublet_pct_p
