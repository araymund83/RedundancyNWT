richness <- qs::qread('./tables/richnessBin1km/richness_allyrs_allgcms.qs')

# Step 1: Average richness across GCMs per pixel/year
rich_mean <- richness %>%
  group_by(x, y, yr) %>%
  summarise(mean_rich = mean(total, na.rm = TRUE), .groups = "drop")

# Step 2: Pivot wider to have columns for each year
rich_wide <- rich_mean %>%
  tidyr::pivot_wider(names_from = yr, values_from = mean_rich)

# Step 3: Compute deltas
rich_wide <- rich_wide %>%
  mutate(delta_2011_2031 = `2031` - `2011`,
         delta_2011_2091 = `2091` - `2011`,
         Longitude = x, Latitude = y)

color_pal <- c('#543005','#8c510a','#bf812d','#dfc27d','#f6e8c3','#c7eae5',
               '#80cdc1','#35978f', '#01665e', '#003c30')

# Step 4: Map function for ∆ richness
plot_delta_richness <- function(df, delta_col, years, out_file){
  lim <- max(abs(df[[delta_col]]), na.rm = TRUE)  # symmetric limits
  gg <- ggplot() +
    geom_tile(data = df, aes(x = Longitude, y = Latitude, fill = .data[[delta_col]])) +
    scale_fill_gradientn(colours = color_pal,
                         limits = c(-lim, lim),
                         na.value = "grey80",
                         name = "\u0394 Richness") +
    geom_sf(data = limt, fill = NA, col = '#36454F') +
    geom_sf(data = ecrg, fill = NA, col = '#36454F') +
    ggtitle("Change in Species Richness",
            subtitle = glue("2011–{years}")) +
    theme_bw() +
    theme(legend.position = 'bottom',
          legend.key.width = unit(2, 'line'),
          plot.title = element_text(size = 16, face = 'bold', hjust = 0),
          plot.subtitle = element_text(size = 14)) +
    labs(x = "Longitude", y = "Latitude") +
    coord_sf()
  
  ggsave(gg, filename = out_file, width = 4, height = 7, dpi = 300)
}

# Step 5: Save the short- and long-term maps
out_dir <- "./figs/richness/2025"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

plot_delta_richness(rich_wide, "delta_2011_2031", 2031,
                    glue("{out_dir}/DeltaRichness_2011_2031.png"))

plot_delta_richness(rich_wide, "delta_2011_2091", 2091,
                    glue("{out_dir}/DeltaRichness_2011_2091.png"))

# Existing panels
p_mean2031 <- glue("{out_dir}/MeanTBI_Bin_2031.png")
p_mean2091 <- glue("{out_dir}/MeanTBI_Bin_2091.png")
p_sig2031  <- glue("{out_dir}/SigCount_Bin_2031.png")
p_sig2091  <- glue("{out_dir}/SigCount_Bin_2091.png")

# New Δ richness panels you generated earlier
p_delta2031 <- glue("{out_dir}/DeltaRichness_2011_2031.png")
p_delta2091 <- glue("{out_dir}/DeltaRichness_2011_2091.png")

stopifnot(
  file.exists(p_mean2031), file.exists(p_mean2091),
  file.exists(p_sig2031),  file.exists(p_sig2091),
  file.exists(p_delta2031), file.exists(p_delta2091)
)

# Read
img_mean31 <- image_read(p_mean2031)
img_mean91 <- image_read(p_mean2091)
img_sig31  <- image_read(p_sig2031)
img_sig91  <- image_read(p_sig2091)
img_dlt31  <- image_read(p_delta2031)
img_dlt91  <- image_read(p_delta2091)

# Make sizes consistent (same width for all)
w <- 2200
img_mean31 <- image_scale(img_mean31, w)
img_mean91 <- image_scale(img_mean91, w)
img_sig31  <- image_scale(img_sig31,  w)
img_sig91  <- image_scale(img_sig91,  w)
img_dlt31  <- image_scale(img_dlt31,  w)
img_dlt91  <- image_scale(img_dlt91,  w)

# Build rows:
row_a <- image_append(c(img_mean31, img_mean91))  # a) Mean TBI 11–31 | 11–91
row_b <- image_append(c(img_sig31,  img_sig91))   # b) #Sig GCMs 11–31 | 11–91
row_c <- image_append(c(img_dlt31,  img_dlt91))   # c) Δ Richness 11–31 | 11–91

# One label per row (top-left of each row)
row_a <- image_annotate(row_a, "a)", size = 80, weight = 700,
                        gravity = "northwest", location = "+36+24", color = "black")
row_b <- image_annotate(row_b, "b)", size = 80, weight = 700,
                        gravity = "northwest", location = "+36+24", color = "black")
row_c <- image_annotate(row_c, "c)", size = 80, weight = 700,
                        gravity = "northwest", location = "+36+24", color = "black")

# Stack rows (3x2 final)
final <- image_append(c(row_a, row_b, row_c), stack = TRUE)

# Optional border
final <- image_border(final, color = "white", geometry = "10x10")

# Save
out_file <- glue("{out_dir}/Composite_a_MeanTBI_b_SigCount_c_DeltaRichness_3x2.png")
image_write(final, out_file)
message(glue("Saved composite: {out_file}"))
