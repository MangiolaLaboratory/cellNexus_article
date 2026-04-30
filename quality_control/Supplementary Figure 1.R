library(cellNexus)
library(dplyr)
library(ggplot2)

setwd("quality_control/")
source("clean_metadata.R")

method_colors <- c(
  "exp" = "#1f77b4",
  "expm1" = "#ff7f0e",
  "identity" = "#2ca02c"
)

# Transformation type and technologies relationship
x <- cell_metadata |>
  distinct(donor_id, sample_id, assay_shorten, assay_groups, inverse_transform) |>
  filter(!is.na(assay_shorten)) |>
  collect()

method_proportion_plot <- x |>
  ggplot(aes(x = assay_shorten, fill = inverse_transform)) +
  geom_bar(position = "fill") +
  scale_y_continuous(
    labels = function(x) x * 100, # convert proportions to 0–100
    limits = c(0, 1), # ensure full range 0–100%
    breaks = seq(0, 1, 0.25)
  ) +
  scale_fill_manual(values = method_colors) +
  labs(
    title = "Proportion of samples by transformation method per technology",
    x = "Technology",
    y = "Sample proportion (%)",
    fill = "Method"
  ) +
  theme_multipanel +
  theme(
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 4),
    legend.key.height = unit(0.2, "cm"),
    legend.key.width = unit(0.3, "cm"),
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)
  )

# ggsave(
#   plot = method_proportion_plot,
#   filename = "method_proportion_plot.pdf",
#   width = 80 * 1.5,
#   height = 60 * 1.5,
#   units = "mm",
#   dpi = 600
# )
