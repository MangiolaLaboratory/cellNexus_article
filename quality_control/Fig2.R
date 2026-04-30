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
library(ggsci)

setwd("quality_control/")
source("plot_custom_theme.R")
source("clean_metadata.R")

### Panel A: Alluvial plot of tissue, disease, assay groups
source("2A_alluvial_plots.R")

### Panel B: Age distribution bar
source("2B_age_distribution_by_sex_ethnicity_per_tissue.R")

### Panel C: scatter plot
source("2C_scatterplot_cell_sample_count_of_dataset.R")

### Panel D: Sample library size density
source("2D_lib_size_density.R")

# save alluvial plot separately
# save heatmap separately
tissue_disease_alluvial |> ggsave(
  filename = "/Users/shen.m/Documents/GitHub/cellNexus-preprint/scripts/plots/tissue_disease_alluvial.pdf",
  width = 60 * 1.5, height = 100 * 1.5, units = "mm"
)

tissue_assay_alluvial |> ggsave(
  filename = "/Users/shen.m/Documents/GitHub/cellNexus-preprint/scripts/plots/tissue_assay_alluvial.pdf",
  width = 60 * 1.5, height = 100 * 1.5, units = "mm"
)

# Patchwork to assemble
alluvial <- (tissue_assay_alluvial + tissue_disease_alluvial) +
  plot_layout(widths = c(1, 1)) & theme(legend.position = "none")
left_bottom <- (alive_per_donor_scatter + p_library_size) + plot_layout(
  guides = "collect",
  widths = c(2, 2)
)
plot_left <- alluvial / plot_spacer() / left_bottom +
  plot_layout(heights = c(5, 0.01, 2))


plot_right <- age_sex_ethnicity_tissue_bar
# + plot_layout(guides = 'auto')
# &
#   theme(legend.position = "bottom")


# Plotting
first_plot <-
  (plot_left | plot_right) +
    plot_layout(guides = "auto", width = c(57, 40)) &
    theme(
      plot.margin = margin(0, 0, 0, 0, "pt"),
      panel.spacing = unit(0, "cm"),
      legend.position = "bottom",
      legend.key.size = unit(0.2, "cm"),
      # legend.spacing.x = unit(0.1, "cm"),    # horizontal gap between keys
      # legend.spacing.y = unit(0.1, "cm"),    # vertical gap (if multi-row legend)
      legend.text = element_text(margin = margin(r = 1, l = 1)), # tighten text padding
      legend.margin = margin(t = 0, b = 0, l = -1, r = 2) # trim margin around the whole legend box
    )


ggsave(
  plot = first_plot,
  filename = "Figure2_dataset_content.pdf",
  width = 183,
  height = 230,
  units = "mm",
  dpi = 600
)
