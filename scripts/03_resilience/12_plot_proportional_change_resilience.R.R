# ================================
# Cross-/Within-scale Δ & SE plots
# ================================

# Libraries ---------------------------------------------------------------
library(pacman)
p_load(glue, ggplot2, qs, sf, tidyverse, cowplot)

# Session/options
g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)

# Paths -------------------------------------------------------------------
dir_path <- './outputs/Redundancy/2025'
out_dir  <- './figs/Redundancy/2025/delta_change'
tbl_dir  <- './tables/redundancy_delta_2025'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tbl_dir,  recursive = TRUE, showWarnings = FALSE)

# Data --------------------------------------------------------------------
# All redundancy/diversity metrics (long)
# Read files --------------------------------------------------------------
dir_path <- ('./outputs/Redundancy/2025')
# Read all redundancy metrics 
data <- qs::qread(glue('{dir_path}/allLogredMetrics_allgcms_allyrs.qs'))
# Read ecoregions file and drop geometry column 
pix_region <- qs::qread('./inputs/pix_df_regionLU_1km.qs')

pix_region<-pix_region %>% select(-geometry)


data <- left_join(data, pix_region, by =  'pixelID')

# Quick sanity (optional)
# print(data %>% count(yr))
# print(data %>% summarise(n=n(), yrs=n_distinct(yr), gcms=n_distinct(gcm), metrics=n_distinct(metric)))

# --- Helper NA-safe stats -----------------------------------------------
mean_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
sd_na   <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) NA_real_ else sd(x)
}

# 1) Means per yr × gcm × region × metric --------------------------------
mean_values <- data %>%
  group_by(yr, gcm, ECO3_NAM_1, metric) %>%
  summarise(
    mean_value = mean_na(value),
    sd_value   = sd_na(value),
    n_pix      = sum(is.finite(value)),
    se_value   = ifelse(n_pix > 1, sd_value / sqrt(n_pix), NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    yr         = as.integer(yr),
    mean_value = ifelse(is.nan(mean_value), NA_real_, mean_value),
    sd_value   = ifelse(is.nan(sd_value),   NA_real_, sd_value),
    se_value   = ifelse(is.nan(se_value),   NA_real_, se_value)
  )

# 2) Pivot wide by year (still per GCM) -----------------------------------
pivot_data <- mean_values %>%
  select(gcm, ECO3_NAM_1, metric, yr, mean_value, se_value) %>%
  tidyr::pivot_wider(
    id_cols    = c(gcm, ECO3_NAM_1, metric),
    names_from = yr,
    values_from = c(mean_value, se_value),
    names_glue = "{.value}_{yr}"
  )

# 3) Absolute change (Δ) per GCM ------------------------------------------
calc_delta <- function(df, yr1, yr2) {
  m1  <- df[[paste0("mean_value_", yr1)]]
  m2  <- df[[paste0("mean_value_", yr2)]]
  se1 <- df[[paste0("se_value_", yr1)]]
  se2 <- df[[paste0("se_value_", yr2)]]
  tibble(
    ECO3_NAM_1 = df$ECO3_NAM_1,
    gcm        = df$gcm,
    metric     = df$metric,
    yr_pair    = glue("{yr1} - {yr2}"),
    delta_mean = ifelse(is.finite(m1) & is.finite(m2), m2 - m1, NA_real_),
    # optional within-GCM SE for the difference (pixel sampling), not used in plots:
    delta_se_withinGCM = ifelse(is.finite(se1) & is.finite(se2), sqrt(se1^2 + se2^2), NA_real_)
  )
}

delta_gcm <- bind_rows(
  calc_delta(pivot_data, "2011", "2031"),
  calc_delta(pivot_data, "2011", "2091")
)

# 4) Summarise across GCMs (model uncertainty) ----------------------------
delta_summary <- delta_gcm %>%
  group_by(ECO3_NAM_1, metric, yr_pair) %>%
  summarise(
    n_gcm        = dplyr::n(),                              # should be 3 everywhere
    mean_delta   = mean(delta_mean, na.rm = TRUE),
    sd_acrossGCM = ifelse(n_gcm > 1, sd(delta_mean, na.rm = TRUE), NA_real_),
    se_acrossGCM = ifelse(n_gcm > 1, sd_acrossGCM / sqrt(n_gcm),   NA_real_),
    .groups = "drop"
  ) %>%
  rename(delta_mean = mean_delta)

# Save tidy table
write_csv(delta_summary, glue("{tbl_dir}/delta_absolute_change_SE_acrossGCM.csv"))

# 5) Plot prep -------------------------------------------------------------
region_order <- c("Mid-Boreal", "High Boreal", "Low Subarctic", "High Subarctic", "Low Arctic north")

# a) Δ cross-scale vs Δ within-scale redundancy (scatter)
red_data <- delta_summary %>%
  filter(metric %in% c("within_scale_red", "cross_scale_red"),
         !is.na(ECO3_NAM_1),
         ECO3_NAM_1 %in% region_order) %>%
  select(ECO3_NAM_1, metric, yr_pair, delta_mean, se_acrossGCM) %>%
  pivot_wider(
    names_from  = metric,
    values_from = c(delta_mean, se_acrossGCM)
  ) %>%
  mutate(ECO3_NAM_1 = factor(ECO3_NAM_1, levels = region_order))

# b) Δ cross-scale diversity (points with SE)
div_data <- delta_summary %>%
  filter(metric == "cross_scale_div",
         !is.na(ECO3_NAM_1),
         ECO3_NAM_1 %in% region_order) %>%
  mutate(ECO3_NAM_1 = factor(ECO3_NAM_1, levels = region_order))

# Theme
custom_theme <- theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# 6) Figures ---------------------------------------------------------------

# (a) Redundancy scatter
plot_delta_changes <- ggplot(
  red_data,
  aes(x = delta_mean_cross_scale_red,
      y = delta_mean_within_scale_red,
      color = yr_pair)
) +
  geom_point(size = 3, alpha = 0.9) +
  geom_errorbarh(aes(
    xmin = delta_mean_cross_scale_red - se_acrossGCM_cross_scale_red,
    xmax = delta_mean_cross_scale_red + se_acrossGCM_cross_scale_red),
    height = 0.05, alpha = 0.6, na.rm = TRUE) +
  geom_errorbar(aes(
    ymin = delta_mean_within_scale_red - se_acrossGCM_within_scale_red,
    ymax = delta_mean_within_scale_red + se_acrossGCM_within_scale_red),
    width = 0.05, alpha = 0.6, na.rm = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  scale_color_manual(NULL, values = c("2011 - 2031" = "deeppink1", "2011 - 2091" = "blueviolet")) +
  labs(
    x = "Δ Cross-scale redundancy",
    y = "Δ Within-scale redundancy",
    title = "Absolute change (Δ): Cross- vs Within-scale redundancy"
  ) +
  facet_wrap(~ ECO3_NAM_1, nrow = 2) +
  custom_theme

# (b) Cross-scale diversity by region
plot_cross_scale_div <- ggplot(
  div_data,
  aes(x = ECO3_NAM_1, y = delta_mean, color = yr_pair)
) +
  geom_point(position = position_dodge(width = 0.5), size = 3, alpha = 0.9) +
  geom_errorbar(aes(ymin = delta_mean - se_acrossGCM,
                    ymax = delta_mean + se_acrossGCM),
                width = 0.2,
                position = position_dodge(width = 0.5),
                alpha = 0.6, na.rm = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_color_manual(NULL, values = c("2011 - 2031" = "deeppink1", "2011 - 2091" = "blueviolet")) +
  labs(
    x = "Region",
    y = "Δ Cross-scale diversity",
    title = "Absolute change (Δ) in cross-scale diversity by region"
  ) +
  custom_theme +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1))

# 7) Save & combine --------------------------------------------------------
ggsave(glue("{out_dir}/delta_change_redundancy_scatter.png"),
       plot_delta_changes, width = 10, height = 8, dpi = 300)

ggsave(glue("{out_dir}/delta_change_cross_scale_div.png"),
       plot_cross_scale_div, width = 10, height = 8, dpi = 300)

aligned <- align_plots(plot_delta_changes, plot_cross_scale_div, align = "v", axis = "l")
combined_figure <- plot_grid(aligned[[1]], aligned[[2]],
                             labels = c("a)", "b)"),
                             ncol = 1, rel_heights = c(1.2, 1))
ggsave(glue("{out_dir}/combined_delta_redundancy_diversity.png"),
       combined_figure, width = 12, height = 10, dpi = 300)

print(combined_figure)


# =========================
# Proportional change plots
# =========================

# Settings for a safe ratio
min_denom <- 1e-3   # minimum |2011| mean to accept proportional change
eps       <- 1e-8   # tiny stabilizer
out_dir_pc <- './figs/Redundancy/2025/prp_change'
dir.create(out_dir_pc, recursive = TRUE, showWarnings = FALSE)

# Helper: proportional change per GCM (guarding baselines)
calc_prop <- function(df, yr1, yr2, min_denom = 1e-3, eps = 1e-8) {
  m1  <- df[[paste0("mean_value_", yr1)]]
  m2  <- df[[paste0("mean_value_", yr2)]]
  se1 <- df[[paste0("se_value_", yr1)]]
  se2 <- df[[paste0("se_value_", yr2)]]
  delta <- m2 - m1
  use_prop <- is.finite(m1) & abs(m1) >= min_denom
  tibble(
    ECO3_NAM_1 = df$ECO3_NAM_1,
    gcm        = df$gcm,
    metric     = df$metric,
    yr_pair    = glue("{yr1} - {yr2}"),
    prp_change = ifelse(use_prop, delta/(abs(m1)+eps), NA_real_),
    # optional: within-GCM SE for the ratio (not used in plots)
    prp_change_se_withinGCM = ifelse(use_prop, sqrt(se1^2 + se2^2)/(abs(m1)+eps), NA_real_),
    base_2011  = m1
  )
}

# Per-GCM proportional change
prop_gcm <- bind_rows(
  calc_prop(pivot_data, "2011", "2031", min_denom, eps),
  calc_prop(pivot_data, "2011", "2091", min_denom, eps)
)

# Summarise across GCMs (model spread → SE across GCMs)
prop_summary <- prop_gcm %>%
  group_by(ECO3_NAM_1, metric, yr_pair) %>%
  summarise(
    n_gcm_valid  = sum(is.finite(prp_change)),
    sd_acrossGCM = ifelse(n_gcm_valid > 1, sd(prp_change, na.rm = TRUE), NA_real_),
    mean_prop    = ifelse(n_gcm_valid > 0, mean(prp_change, na.rm = TRUE), NA_real_),
    se_acrossGCM = ifelse(n_gcm_valid > 1, sd_acrossGCM / sqrt(n_gcm_valid), NA_real_),
    .groups = "drop"
  ) %>%
  rename(prp_change = mean_prop)

# Save tidy table
write_csv(prop_summary, glue("{tbl_dir}/proportional_change_SE_acrossGCM_thresholded.csv"))

# =========================
# Plot proportional change
# =========================

# Fixed region order (and drop any NA region rows)
region_order <- c("Mid-Boreal", "High Boreal", "Low Subarctic", "High Subarctic", "Low Arctic north")

prop_summary_clean <- prop_summary %>%
  filter(!is.na(ECO3_NAM_1), ECO3_NAM_1 %in% region_order) %>%
  mutate(ECO3_NAM_1 = factor(ECO3_NAM_1, levels = region_order))

# A) Cross- vs Within-scale redundancy (proportional change)
red_prop <- prop_summary_clean %>%
  filter(metric %in% c("within_scale_red", "cross_scale_red")) %>%
  select(ECO3_NAM_1, yr_pair, metric, prp_change, se_acrossGCM) %>%
  pivot_wider(
    names_from  = metric,
    values_from = c(prp_change, se_acrossGCM)
  )

plot_prp_redundancy <- ggplot(
  red_prop,
  aes(x = prp_change_cross_scale_red,
      y = prp_change_within_scale_red,
      color = yr_pair)
) +
  geom_point(size = 3, alpha = 0.9, na.rm = TRUE) +
  geom_errorbarh(
    aes(xmin = prp_change_cross_scale_red - se_acrossGCM_cross_scale_red,
        xmax = prp_change_cross_scale_red + se_acrossGCM_cross_scale_red),
    height = 0.05, alpha = 0.6, na.rm = TRUE
  ) +
  geom_errorbar(
    aes(ymin = prp_change_within_scale_red - se_acrossGCM_within_scale_red,
        ymax = prp_change_within_scale_red + se_acrossGCM_within_scale_red),
    width = 0.05, alpha = 0.6, na.rm = TRUE
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  scale_color_manual(NULL, values = c("2011 - 2031" = "deeppink1",
                                      "2011 - 2091" = "turquoise3")) +
  labs(
    x = "Proportional change — Cross-scale redundancy",
    y = "Proportional change — Within-scale redundancy",
    title = "Proportional change: Cross- vs Within-scale redundancy"
  ) +
  facet_wrap(~ ECO3_NAM_1, nrow = 2) +
  theme_bw() +
  theme(
    legend.position = 'right',
    legend.key.width = unit(2, 'line'),
    plot.title = element_text(size = 16, face = 'bold', hjust = 0.5, vjust = 0.7),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 11),
    strip.text = element_text(size = 14, margin = margin(3, 0, 3, 0))
  ) 

# B) Cross-scale diversity (proportional change) by region
div_prop <- prop_summary_clean %>%
  filter(metric == "cross_scale_div")

plot_prp_div <- ggplot(
  div_prop,
  aes(x = ECO3_NAM_1, y = prp_change, color = yr_pair)
) +
  geom_point(position = position_dodge(0.5), size = 3, alpha = 0.9, na.rm = TRUE) +
  geom_errorbar(aes(ymin = prp_change - se_acrossGCM,
                    ymax = prp_change + se_acrossGCM),
                width = 0.2,
                position = position_dodge(0.5),
                alpha = 0.6, na.rm = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_color_manual(NULL, values = c("2011 - 2031" = "deeppink1",
                                      "2011 - 2091" = "turquoise3")) +
  labs(
    x = "Region",
    y = "Proportional change — Cross-scale diversity",
    title = "Proportional change in cross-scale diversity by region"
  ) +
  theme_bw() +
  theme(
    legend.position = 'right',
    legend.key.width = unit(2, 'line'),
    plot.title = element_text(size = 16, face = 'bold', hjust = 0.5, vjust = 0.7),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 11),
    strip.text = element_text(size = 14, margin = margin(3, 0, 3, 0))
  ) 

# Save
out_dir_pc <- './figs/Redundancy/2025/prp_change'
if (!dir.exists(out_dir_pc)) dir.create(out_dir_pc, recursive = TRUE)

ggsave(file.path(out_dir_pc, "proportional_change_redundancy_scatter.png"),
       plot_prp_redundancy, width = 10, height = 8, dpi = 300)
ggsave(file.path(out_dir_pc, "proportional_change_cross_scale_div.png"),
       plot_prp_div, width = 10, height = 8, dpi = 300)

# Combine (optional)
aligned_pc <- cowplot::align_plots(plot_prp_redundancy, plot_prp_div, align = "v", axis = "l")
combined_pc <- cowplot::plot_grid(aligned_pc[[1]], aligned_pc[[2]],
                                  labels = c("a)", "b)"),
                                  ncol = 1, rel_heights = c(1.2, 1))
ggsave(file.path(out_dir_pc, "combined_proportional_redundancy_diversity.png"),
       combined_pc, width = 12, height = 10, dpi = 300)

print(combined_pc)
