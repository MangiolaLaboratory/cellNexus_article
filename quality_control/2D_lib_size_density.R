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
library(ggridges)

set.seed(123)
source("plot_custom_theme.R")
source("clean_metadata.R")
# 1. Density plot of library size color by technology ---------
# Total library size per sample 
density_df = cell_metadata |>
  group_by(assay_shorten, sample_id) |>
  summarise(
    total_library_size = sum(nCount_RNA, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop") |> 
  collect()

# Extract density peak
label_df <- density_df %>%
  mutate(log_total_library_size = log10(total_library_size)) %>%
  group_by(assay_shorten) %>%
  summarise(peak = list({
    d <- density(log_total_library_size, na.rm = TRUE)
    tibble(x = d$x, y = d$y) |> slice_max(y, n = 1)
  })) %>%
  tidyr::unnest(peak) %>% 
  mutate(mode_original_x = 10^x) 

p_library_size <- density_df |> 
  left_join(label_df) |>
  mutate(y_norm = log10(y)) |> 
  ggplot(aes(x = total_library_size, y = reorder(assay_shorten, total_library_size, decreasing=T), fill = y_norm)) + 
  geom_density_ridges_gradient(scale = 1.5, linewidth=0.1) + 
  scale_x_log10() +
  scale_fill_viridis(name = expression("log"[10]*"(Density Peak)"), option="D",
                     guide = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  labs(x="Samples Total Library Size", y=NULL) + 
  theme_multipanel +
  theme(
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 4),
    legend.key.height = unit(0.2, "cm"),
    legend.key.width = unit(0.3, "cm")
  ) 
  #theme(legend.position = "none") +
  #guides(fill = "none", color = "none")


p_library_size


