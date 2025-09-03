data <- qs::qread(glue('./outputs/Redundancy/2025/allLogredMetrics_allgcms_allyrs.qs'))
pix_region <- qs::qread('./inputs/pix_df_regionLU_1km.qs')

pix_region<-pix_region %>% select(-geometry)


data <- left_join(data, pix_region, by =  'pixelID')

# Load map information ----------------------------------------------

targetCRS <- paste("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95",
                   "+x_0=0 +y_0=0 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")

limt <- sf::st_read('inputs/NT1_BCR6/NT1_BCR6_poly.shp') 
ecrg <- sf::st_read('inputs/ecoregions/NWT_ecoregions_5Regions.shp')
# Extract by mask for the ecoregions ---------------------------------------
plot(st_geometry(ecrg))
limt <- sf::st_transform(x = limt, crs = targetCRS)
ecrg <- sf::st_transform(x = ecrg, crs = targetCRS)
plot(st_geometry(ecrg))
plot(st_geometry(limt))

# new colors for maps
#color_pal <- c('#fff7fb','#ece2f0','#d0d1e6','#a6bddb','#67a9cf','#3690c0','#02818a','#016450')
color_pal <- c('#f7fcfd','#e0ecf4','#bfd3e6','#9ebcda','#8c96c6','#8c6bb1','#88419d','#810f7c','#4d004b')
#write the name of the metric
#vars <- c('cross_scale_red', 'within_scale_red', 'cross_scale_div')
vars <- unique(data$metric)
summary_stats <- data %>% 
  filter(metric %in% c("cross_scale_red", "within_scale_red", "cross_scale_div")) %>%
  group_by(yr, lat, long, metric) %>% 
  summarise(mean_value = mean(value, na.rm = TRUE), 
            sd_value = sd(value, na.rm = TRUE), 
            se_value = sd(value, na.rm = TRUE) / sqrt(n()), .groups = "drop")

summary_2011 <- summary_stats %>% 
  filter(yr == 2011)


make_facet_map_2011 <- function(summary_2011, limt, ecrg) {
  
  # Define better facet labels
  var_titles <- c(
    'cross_scale_red' = 'Cross-scale redundancy',
    'within_scale_red' = 'Within-scale redundancy',
    'cross_scale_div' = 'Cross-scale diversity'
  )
  facet_labels <- as_labeller(var_titles)
  
  # Global max for color scale
  max_val <- max(summary_2011$mean_value, na.rm = TRUE)
  
  # Main plot
  p <- ggplot(summary_2011) +
    geom_tile(aes(x = long, y = lat, fill = mean_value)) +
    scale_fill_gradientn(
      colours = color_pal,
      limits = c(0, max_val),
      na.value = 'grey80',
      name = "Mean value"
    ) +
    geom_sf(data = limt, fill = NA, color = '#36454F', linewidth = 0.3) +
    geom_sf(data = ecrg, fill = NA, color = '#36454F', linewidth = 0.3) +
    facet_wrap(~metric, labeller = facet_labels, nrow = 1) +
    ggtitle("Redundancy Metrics – 2011") +
    labs(x = "Longitude", y = "Latitude") +
    theme_bw() +
    theme(
      legend.position = 'bottom',
      legend.key.width = unit(2, 'line'),
      plot.title = element_text(size = 16, face = 'bold'),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.text = element_text(size = 11),
      legend.title = element_text(size = 12, face = 'bold'),
      strip.text = element_text(size = 13)
    ) +
    coord_sf(crs = st_crs(ecrg), expand = FALSE)
  
  # Save
  out_dir <- './figs/Redundancy/2011'
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  ggsave(filename = glue("{out_dir}/facet_metrics_mean_2011.png"),
         plot = p, width = 15, height = 5)
  
  return(p)
}

summary_2011 <- summary_stats %>% filter(yr == 2011)
facet_plot <- make_facet_map_2011(summary_2011, limt, ecrg)


make_and_save_map <- function(df, metric_name, title_label, filename_suffix) {
  df_metric <- df %>% filter(metric == metric_name, !is.na(mean_value))
  max_val <- max(df_metric$mean_value, na.rm = TRUE)
  
  p <- ggplot(df_metric) +
    geom_tile(aes(x = long, y = lat, fill = mean_value)) +
    scale_fill_gradientn(
      colours = color_pal,
      limits = c(0, max_val),
      name = "Mean value",
      na.value = "grey80"
    ) +
    geom_sf(data = limt, fill = NA, color = "#36454F", linewidth = 0.3) +
    geom_sf(data = ecrg, fill = NA, color = "#36454F", linewidth = 0.3) +
    ggtitle(title_label) +
    labs(x = "Longitude", y = "Latitude") +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.key.width = unit(2, "line"),
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 11, face = "bold")
    ) +
    coord_sf(crs = st_crs(ecrg), expand = FALSE)
  
  filename <- glue("{out_dir}/map_{filename_suffix}.png")
  ggsave(filename = filename, plot = p, width = 4, height = 6, dpi = 300)
  
  return(filename)
}

# Create and save individual maps
file1 <- make_and_save_map(summary_2011, "cross_scale_red", var_titles[["cross_scale_red"]], "a_cross_scale_red")
file2 <- make_and_save_map(summary_2011, "within_scale_red", var_titles[["within_scale_red"]], "b_within_scale_red")
file3 <- make_and_save_map(summary_2011, "cross_scale_div", var_titles[["cross_scale_div"]], "c_cross_scale_div")


library(magick)

# Read and annotate
img1 <- image_read(file1) %>% image_annotate("a)", size = 50, gravity = "northwest", location = "+30+30", color = "black", weight = 700)
img2 <- image_read(file2) %>% image_annotate("b)", size = 50, gravity = "northwest", location = "+30+30", color = "black", weight = 700)
img3 <- image_read(file3) %>% image_annotate("c)", size = 50, gravity = "northwest", location = "+30+30", color = "black", weight = 700)

# Stack vertically (or use image_append for side-by-side)
final_combined <- image_append(c(img1, img2, img3), stack = FALSE)

# Save final image
final_path <- glue("{out_dir}/combined_redundancy_maps_2011_labeled.png")
image_write(final_combined, path = final_path, format = "png", density = "300x300")
