library(pacman)
p_load(glue, ggplot2,ggalluvial, glue, qs, tidyverse)

# Session/options
g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)

# Create output directory
dir.create("./outputs/fg_presence_by_region", showWarnings = FALSE)

# Define parameters
regions <- c("High_Boreal", "Mid_Boreal", "Low_Subarctic", 
             "High_Subarctic", "Low_Arctic_north")
gcms <- c("CanESM2", "CCSM4", "INM-CM4")
yrs <- c(2011, 2031, 2091)
path <- "./tables/fg_abun/2025"
files <- list.files('./tables/fg_abun/2025')

library(tidyverse)
library(glue)
library(qs)

regions <- c("High Boreal", "Mid-Boreal", "Low Subarctic", 
             "High Subarctic", "Low Arctic north")
gcms <- c("CanESM2", "CCSM4", "INM-CM4")
yrs <- c(2011, 2031, 2091)

# Read pixel file (once)
pix_file <- "./inputs/pix_df_regionLU_1km.qs"
pix_sample <- qread(pix_file)

# Function to read and join each fg file with pixel metadata
read_fg_with_meta <- function(gcm, yr) {
  fg_file <- glue("./tables/fg_abun/2025/fg_{gcm}_{yr}.qs")
  if (!file.exists(fg_file)) return(NULL)
  
  fg_df <- qread(fg_file)
  df <- left_join(pix_sample, fg_df, by = "pixelID") %>%
    mutate(gcm = gcm, yr = yr, region = ECO3_NAM_1)
  return(df)
}

# Read and join all combinations
fg_all <- crossing(gcm = gcms, yr = yrs) %>%
  mutate(data = map2(gcm, yr, read_fg_with_meta)) %>%
  filter(!map_lgl(data, is.null)) %>%
  pull(data) %>%
  bind_rows()

# Keep only relevant columns
fg_long <- fg_all %>%
  pivot_longer(cols = starts_with(c("INS", "OMN", "GRA")),
               names_to = "FG", values_to = "presence") %>%
  mutate(presence = if_else(presence > 0, 1, 0))  # binarize


# Create transition label for 2011 → 2031
transitions_31 <- fg_long %>%
  filter(yr %in% c(2011, 2031)) %>%
  select(pixelID, region, gcm, yr, FG, presence) %>%
  pivot_wider(names_from = yr, values_from = presence, names_prefix = "yr_") %>%
  mutate(transition = case_when(
    yr_2011 == 0 & yr_2031 == 1 ~ "appear",
    yr_2011 == 1 & yr_2031 == 0 ~ "disappear",
    yr_2011 == 1 & yr_2031 == 1 ~ "persist",
    TRUE ~ "absent"
  ))

transitions_91 <- fg_long %>%
  filter(yr %in% c(2011, 2091)) %>%
  select(pixelID, region, gcm, yr, FG, presence) %>%
  pivot_wider(names_from = yr, values_from = presence, names_prefix = "yr_") %>%
  mutate(transition = case_when(
    yr_2011 == 0 & yr_2091 == 1 ~ "appear",
    yr_2011 == 1 & yr_2091 == 0 ~ "disappear",
    yr_2011 == 1 & yr_2091 == 1 ~ "persist",
    TRUE ~ "absent"
  ))
# Step 1: Prepare summary_df with fg_rank for plotting
# Define transition styles
transition_styles <- c(
  "appear" = "dotdash",
  "disappear" = "dotted",
  "persist" = "solid"
)

# Function to prepare the data
make_arrow_df <- function(df, year) {
  df %>%
    pivot_longer(cols = c("appear", "disappear", "persist"),
                 names_to = "transition", values_to = "n") %>%
    filter(n > 0) %>%
    mutate(
      FG = factor(FG, levels = fg_levels),
      fg_rank = as.numeric(FG),  # y-axis
      x_start = 2011,
      x_end = year,
      region = case_when(
        region == "Mid-Boreal"       ~ "Mid_Boreal",
        region == "High Boreal"      ~ "High_Boreal",
        region == "Low Subarctic"    ~ "Low_Subarctic",
        region == "High Subarctic"   ~ "High_Subarctic",
        region == "Low Arctic north" ~ "Low_Arctic_north",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(region)) %>%
    mutate(region = factor(region, levels = region_order))
}

# Create combined dataset
df_31 <- make_arrow_df(summary_31, 2031)
df_91 <- make_arrow_df(summary_91, 2091)
arrow_data <- bind_rows(df_31, df_91)

# Custom FG color palette
fg_palette_named <- c(
  "INS-Sal" = "#1B9E77", "OMN-For" = "#D95F02", "INS-Gle" = "#7570B3",
  "INS-Exc" = "#E7298A", "INS-For" = "#66A61E", "OMN-Sca" = "#E6AB02",
  "OMN-Gle" = "#A6761D", "GRA-Gle" = "#A6CEE3", "INS-Scr" = "#FB9A99",
  "OMN-Exc" = "#B2DF8A"
)

# Prepare data
arrow_data2 <- arrow_data %>%
  filter(transition %in% c("appear", "disappear")) %>%
  mutate(
    transition = factor(transition, levels = c("appear", "disappear")),
    FG = fct_rev(FG)  # So highest rank is top of plot
  )

# Plot with legend for transition line types
fg_plot <- ggplot(arrow_data2, aes(x = x_start, xend = x_end, y = FG, yend = FG, color = FG, linetype = transition)) +
  geom_segment(arrow = arrow(length = unit(0.15, "inches")),
               linewidth = 0.6) +
  scale_color_manual(values = fg_palette_named) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  scale_x_continuous(breaks = c(2011, 2031, 2091), limits = c(2005, 2095)) +
  facet_grid(gcm ~ region) +
  labs(
    title = "Functional Group Transitions by Region and GCM",
    x = "Year",
    y = "Functional Group",
    color = "FG",
    linetype = "Transition"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 10),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(1, "lines"),
    legend.position = "right"
  )

out_dir <- './figs/persistance'
ggsave(glue("{out_dir}/FG_persistance_byGCM_region.png"),
       fg_plot, width = 12, height = 8, dpi = 300)

