library(ggplot2)
library(ggsci)
# Theme
friendly_cols <- dittoSeq::dittoColors()

npg_cols <- c(pal_npg("nrc")(7),"gray40","gray")
my_colors <- setNames(npg_cols, 
  c("European",
    "East Asian",
    "Japanese",
    "South Asian",
    "African",
    "Hispanic/Latin American",
    "Native American & Pacific Islander",
    "LowConfidenceLabel",
    "Other/Unknown")
)

theme_multipanel =
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(size=0.1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    strip.background = element_blank(),
    axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0), size = 7),
    axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0), size = 7),
    panel.spacing.x=unit(0.1, "lines"),
    axis.text.x = element_text(size=6),
    axis.text.y = element_text(size=6),
    strip.text.x = element_text(size = 7),
    strip.text.y = element_text(size = 7),
    
    # legend
    legend.key.size = unit(5, 'mm'),
    legend.title = element_text(size=7),
    legend.text = element_text(size=6),
    
    # Avoid text clipping for facets. Currently not merged remotes::install_github("tidyverse/ggplot2#4223")
    strip.clip = "off",
    
    # Title
    plot.title = element_text(size=7),
    
    axis.line.x = element_line(size=0.2),
    axis.line.y = element_line(size=0.2),
    axis.ticks.x = element_line(size=0.2),
    axis.ticks.y = element_line(size=0.2)
  )

# Useful if we have a lot of decimals
# scale_x_continuous( labels = dropLeadingZero)
dropLeadingZero <- function(l){  stringr::str_replace(l, '0(?=.)', '') }

# Patchwork
# +  plot_layout(guides = 'collect' ) + plot_annotation(tag_levels = c('A')) & theme( plot.margin = margin(0, 0, 0, 0, "pt"),  legend.key.size = unit(0.2, 'cm'))

custom_theme <-
  list(
    scale_fill_manual(values = friendly_cols),
    scale_color_manual(values = friendly_cols),
    theme_bw() +
      theme(
        panel.border = element_blank(),
        axis.line = element_line(),
        panel.grid.major = element_line(size = 0.2),
        panel.grid.minor = element_line(size = 0.1),
        text = element_text(size = 12),
        legend.position = "bottom",
        strip.background = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
        axis.title.y = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)
      )
  )
