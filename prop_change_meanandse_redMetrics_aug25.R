# Read files --------------------------------------------------------------
dir_path <- ('./outputs/Redundancy/2025')
# Read all redundancy metrics 
data <- qs::qread(glue('{dir_path}/allLogredMetrics_allgcms_allyrs.qs'))
# Read ecoregions file and drop geometry column 
pix_region <- qs::qread('./inputs/pix_df_regionLU_1km.qs')

pix_region<-pix_region %>% select(-geometry)


data <- left_join(data, pix_region, by =  'pixelID')

xsr <- data %>% filter(metric == 'cross_scale_red')
wsr <- data %>% filter(metric == 'within_scale_red')
xsd <- data %>% filter(metric == 'cross_scale_div')


# 1) Keep a single row per pixelID×gcm×metric×yr (no mean; take first if duplicates)
make_wide_with_changes <- function(df, metric_name) {
  # Baseline coords/regions from 2011 (kept in output)
  base_2011 <- df %>%
    filter(yr == 2011) %>%
    distinct(pixelID, gcm, lat, long, ECO3_NAM_1)
  
  wide <- df %>%
    filter(yr %in% c(2011, 2031, 2091)) %>%
    arrange(pixelID, gcm, yr) %>%
    distinct(pixelID, gcm, yr, .keep_all = TRUE) %>%  # ensure 1 row per key; no averaging
    select(pixelID, gcm, yr, value) %>%
    pivot_wider(names_from = yr, values_from = value, names_prefix = "yr") %>%
    left_join(base_2011, by = c("pixelID","gcm")) %>%
    mutate(
      prop_2011_2031 = if_else(!is.na(yr2011) & yr2011 != 0,
                               (yr2031 - yr2011) / yr2011, NA_real_),
      prop_2011_2091 = if_else(!is.na(yr2011) & yr2011 != 0,
                               (yr2091 - yr2011) / yr2011, NA_real_),
      metric = metric_name
    ) %>%
    relocate(metric, .before = pixelID) %>%
    # keep year values + changes + context columns
    select(metric, pixelID, gcm, lat, long, ECO3_NAM_1,
           yr2011, yr2031, yr2091, prop_2011_2031, prop_2011_2091)
  
  wide
}

# Build per-metric wide tables
xsr_wide <- make_wide_with_changes(xsr, "cross_scale_red")
wsr_wide <- make_wide_with_changes(wsr, "within_scale_red")
xsd_wide <- make_wide_with_changes(xsd, "cross_scale_div")

# Combine all metrics into a single df (keeps all year values + changes)
all_metrics_wide <- dplyr::bind_rows(xsr_wide, wsr_wide, xsd_wide) %>% 
  filter(!is.na(ECO3_NAM_1))

# Inputs assumed already built:
# xsr_wide (cross_scale_red), wsr_wide (within_scale_red), xsd_wide (cross_scale_div)
# each has: metric, pixelID, gcm, ECO3_NAM_1, prop_2011_2031, prop_2011_2091

# 0) Bind all metrics and go long on timeframes
melted <- melt(
  all_metrics_wide,
  id.vars = c("metric", "ECO3_NAM_1", "gcm"),  # keys to keep
  measure.vars = c("prop_2011_2031", "prop_2011_2091"),
  variable.name = "timeframe",
  value.name = "prop_change"
)

# ----- 2) Mean per GCM x Region (averaging pixels within region for that GCM) -----
#         This makes each GCM contribute equally when we compute SE across GCMs.
gcm_region <- melted %>%
  group_by(metric, ECO3_NAM_1, timeframe, gcm) %>%
  summarise(gcm_mean = mean(prop_change, na.rm = TRUE), .groups = "drop")

# ----- 3) Mean ± SE across GCMs (what you asked for) -----
region_stats <- gcm_region %>%
  group_by(metric, ECO3_NAM_1, timeframe) %>%
  summarise(
    n_gcm = sum(!is.na(gcm_mean)),
    mean  = ifelse(n_gcm > 0, mean(gcm_mean, na.rm = TRUE), NA_real_),
    sd    = ifelse(n_gcm > 1, sd(gcm_mean,   na.rm = TRUE), NA_real_),
    se    = ifelse(n_gcm > 1, sd / sqrt(n_gcm), NA_real_),
    .groups = "drop"
  )


# Pretty labels for timeframe
region_stats_labeled <- region_stats %>%
  mutate(
    yr_pair = recode(as.character(timeframe),
                     "prop_2011_2031" = "2011 - 2031",
                     "prop_2011_2091" = "2011 - 2091")
  )

# ---- A) Cross-scale vs Within-scale redundancy (mean ± SE across GCMs) ----
red_prop <- region_stats_labeled %>%
  filter(metric %in% c("cross_scale_red", "within_scale_red")) %>%
  select(ECO3_NAM_1, yr_pair, metric, mean, se) %>%
  pivot_wider(
    names_from  = metric,
    values_from = c(mean, se)
  )
# Optional: order regions (edit to your liking)
# region_order <- c("Mid-Boreal","High Boreal","Low Subarctic","High Subarctic","Low Arctic north")
# red_prop <- red_prop %>% mutate(ECO3_NAM_1 = fct_relevel(ECO3_NAM_1, region_order))

plot_prp_redundancy <- ggplot(
  red_prop,
  aes(x = mean_cross_scale_red,
      y = mean_within_scale_red,
      color = yr_pair)
) +
  geom_point(size = 3, alpha = 0.9, na.rm = TRUE) +
  geom_errorbarh(
    aes(xmin = mean_cross_scale_red - se_cross_scale_red,
        xmax = mean_cross_scale_red + se_cross_scale_red),
    height = 0.05, alpha = 0.7, na.rm = TRUE
  ) +
  geom_errorbar(
    aes(ymin = mean_within_scale_red - se_within_scale_red,
        ymax = mean_within_scale_red + se_within_scale_red),
    width = 0.05, alpha = 0.7, na.rm = TRUE
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

# ---- B) Cross-scale diversity (mean ± SE across GCMs) ----
csd_df <- region_stats_labeled %>%
  filter(metric == "cross_scale_div")

plot_prp_div <- ggplot(
  csd_df,
  aes(x = ECO3_NAM_1, y = mean, color = yr_pair)
) +
  geom_point(size = 3, position = position_dodge(width = 0.5), na.rm = TRUE) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.2, position = position_dodge(width = 0.5), alpha = 0.7, na.rm = TRUE
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_color_manual(NULL, values = c("2011 - 2031" = "deeppink1",
                                      "2011 - 2091" = "turquoise3")) +
  labs(
    x = "Region", y = "Proportional change — Cross-scale diversity",
    title = "Proportional change in cross-scale diversity by region"
  ) +
  theme_bw() +
  theme(
    legend.position = 'right',
    legend.key.width = unit(2, 'line'),
    plot.title = element_text(size = 16, face = 'bold', hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(size = 12, angle = 25, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

# Save directory
out_dir_pc <- './figs/Redundancy/2025/prp_change'
if (!dir.exists(out_dir_pc)) dir.create(out_dir_pc, recursive = TRUE)

# Save individual panels
ggsave(file.path(out_dir_pc, "proportional_change_redundancy_scatter.png"),
       plot_prp_redundancy, width = 10, height = 8, dpi = 300, bg = "white")
ggsave(file.path(out_dir_pc, "proportional_change_cross_scale_div.png"),
       plot_prp_div, width = 10, height = 8, dpi = 300, bg = "white")


library(cowplot)
library(ggplot2)

# Ensure output dir exists
if (!dir.exists(out_dir_pc)) dir.create(out_dir_pc, recursive = TRUE)

# --- Option 1: keep separate legends on each panel (simplest) ---
aligned_pc <- cowplot::align_plots(
  plot_prp_redundancy,
  plot_prp_div,
  align = "v", axis = "l"
)

combined_pc <- cowplot::plot_grid(
  aligned_pc[[1]], aligned_pc[[2]],
  labels = c("a)", "b)"),
  ncol = 1, rel_heights = c(1.2, 1), label_size = 16, label_fontface = "bold"
)

ggsave(file.path(out_dir_pc, "combined_proportional_redundancy_diversity.png"),
       combined_pc, width = 12, height = 10, dpi = 300, bg = "white")

print(combined_pc)

