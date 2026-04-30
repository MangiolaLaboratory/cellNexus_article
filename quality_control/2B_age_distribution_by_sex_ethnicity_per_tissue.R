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
source("CAQ_age_analysis_functions.R")
source("clean_metadata.R")

donor_metadata <- cell_metadata |>
  distinct(donor_id, self_reported_ethnicity, age_days, sex, tissue_groups) |>
  as_tibble()

age_groups_vector <- coarse_age_bin(age_days = donor_metadata |> pull(age_days), sex = donor_metadata |> pull(sex))
donor_metadata <- donor_metadata |>
  mutate(age_groups = age_groups_vector) |>
  left_join(ethnicity_grouped)


donor_count_df <- donor_metadata |>
  dplyr::count(ethnicity_groups, age_groups, sex, tissue_groups) |>
  mutate(age_groups = ifelse(is.na(age_groups), "Unknown", age_groups)) |>
  as_tibble() |>
  mutate(
    tissue_groups = stringr::str_replace_all(tissue_groups, "/", "/\n"),
    tissue_groups = stringr::str_wrap(tissue_groups, width = 10)
  )

# Reorder tissue_groups by total donor count (n)
tissue_order <- donor_count_df |>
  group_by(tissue_groups) |>
  summarise(total_n = sum(n), .groups = "drop") |>
  arrange(desc(total_n)) |>
  pull(tissue_groups)

# calculate top tissues
top_tissues <- donor_count_df |>
  group_by(tissue_groups) |>
  summarise(total_n = sum(n), .groups = "drop") |>
  arrange(desc(total_n)) %>%
  # Only filter Top Tissue groups for better visualisation purpose
  slice_head(n = 25) %>%
  pull(tissue_groups)


df <- donor_count_df |>
  filter(age_groups != "Unknown") |>
  filter(tissue_groups %in% top_tissues) |>
  mutate(tissue_groups = factor(tissue_groups, levels = tissue_order)) |>
  mutate(n = ifelse(sex == "female", -n, n)) |> # Negative for Female to diverge
  dplyr::rename(donor_count = n)

# # double log
# df$signed_double_log10_value <- sign(df$donor_count) * log10(1 + log10(1 + abs(df$donor_count)))
# df$signed_log10_value <- sign(df$donor_count) * log10(abs(df$donor_count))
df$age_groups <- factor(df$age_groups, levels = c(">60", "30-60", "<30"))
df <- df |> mutate(Sign = sign(donor_count))

npg_cols <- c(pal_npg("nrc")(7), "gray40", "gray")
ethnicity_colors <- setNames(
  npg_cols,
  c(
    "European",
    "East Asian",
    "Japanese",
    "South Asian",
    "African",
    "Hispanic/Latin American",
    "Native American & Pacific Islander",
    "LowConfidenceLabel",
    "Other/Unknown"
  )
)

age_sex_ethnicity_tissue_bar <- df |>
  ggplot(aes(x = age_groups, y = log(abs(donor_count)) * Sign, fill = ethnicity_groups)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  facet_grid(tissue_groups ~ ., switch = "y") +
  scale_fill_manual(
    values = ethnicity_colors,
    guide = guide_legend(nrow = 4)
  ) +
  scale_color_manual(values = ethnicity_colors) +
  labs(
    title = "Signed Log-transformed Donor Counts by Age Group (Female/Male)",
    x = "",
    y = expression(log(abs(donor[count])) %.% Sign),
    fill = "Ethnicity"
  ) +
  theme_multipanel +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 6, face = "bold"),
    strip.background = element_blank(),

    # Adjust legend size
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 5),
    legend.key.height = unit(0.05, "cm"),
    legend.key.width = unit(0.05, "cm")
  )


age_sex_ethnicity_tissue_bar
# ggsave(plot = age_sex_ethnicity_tissue_bar,
#        filename = "~/projects/cellNexus/preprint/plots/age_distribution_by_sex_ethnicity_per_tissue.png",
#        width = 20,
#        height = 40,
#        units = "in")
