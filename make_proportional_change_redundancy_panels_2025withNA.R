# Load libraries ----------------------------------------------------------
library(pacman)
p_load(glue, ggplot2, sf, tidyverse)

# Clear environment and reset options
g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)

# Read files --------------------------------------------------------------
dir_path <- ('./outputs/Redundancy/2025')
# Read all redundancy metrics 
fls <- list.files(dir_path, pattern = 'allLong',full.names = TRUE)

# Filter only the *_allyrs.qs files
allyrs_fls <- fls[str_detect(fls, "_allyrs\\.qs$")]
#read the files
all_long_metrics_df <- map_dfr(allyrs_fls, qs::qread)
#save the file for future
data<- qs::qread(glue('{dir_path}/allLogredMetrics_allgcms_allyrs.qs'))

# Read ecoregions file and drop geometry column 
pix_region <- qs::qread('./inputs/pix_df_regionLU_1km.qs')

pix_region <- pix_region %>% select(-geometry)


data <- left_join(data, pix_region, by = 'pixelID')

# Create a vector of the unique metrics 
metrics <- unique(data$metric)
yrs <- unique(data$yr)

# Load map information ----------------------------------------------
targetCRS <- paste("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95",
                   "+x_0=0 +y_0=0 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")

limt <- sf::st_read('inputs/NT1_BCR6/NT1_BCR6_poly.shp') 
ecrg <- sf::st_read('inputs/ecoregions/NWT_ecoregions_5Regions.shp')

# Extract by mask for the ecoregions 
limt <- sf::st_transform(x = limt, crs = targetCRS)
ecrg <- sf::st_transform(x = ecrg, crs = targetCRS)
plot(st_geometry(ecrg))
plot(st_geometry(limt))

# Calculate deltas data ----------------------------------------------------------
calculate_deltas <- function(summary_stats, metric, yr_pairs) {

  delta_list <- lapply(yr_pairs, function(yr_pair) {
    
    yr1 <- yr_pair[1]
    yr2 <- yr_pair[2]
   # browser()
    
    delta_data <- summary_stats %>%
      filter(metric == metric, yr %in% c(yr1, yr2)) %>%
      pivot_wider(names_from = yr, values_from = c(mean_value, sd_value, se_value))
    
    delta_data <- delta_data %>%
      mutate(delta_mean = get(paste0("mean_value_", yr2)) - get(paste0("mean_value_", yr1)),
             prp_change= delta_mean / get(paste0("mean_value_", yr1)), 
             yr_pair = glue('{yr1} - {yr2}'))
    
    return(delta_data)
  })
  
  bind_rows(delta_list)
}

# Define year pairs
yr_pairs <- list(c(yrs[1], yrs[2]), c(yrs[1], yrs[3]))


summary_stats <- data %>%
  filter(metric %in% c("cross_scale_red", "within_scale_red", "cross_scale_div")) %>%
  group_by(yr, lat, long, metric) %>%
  summarise(mean_value = mean(value, na.rm = TRUE),
            sd_value = sd(value, na.rm = TRUE),
            se_value = sd(value, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")

# Calculate deltas for each metric
deltas_list <- lapply(metrics, function(metric) {
  calculate_deltas(summary_stats, metric, yr_pairs)
})

# Combine deltas into a single dataframe
deltas_df <- bind_rows(deltas_list)
qs::qsave(deltas_df, './tables/redundancy_delta_2025/delta_df_prp_change.qs')
deltas_df <- qs::qread('./tables/redundancy_delta_2025/delta_df_prp_change.qs')

deltas_df <- deltas_df %>%
  filter(metric == 'within_scale_red')
# Calculate the maximum absolute value across all metrics
min_max_abs_value <- deltas_df %>%
  filter(is.finite(prp_change)) %>%
  group_by(metric) %>% 
  summarise(min_abs_prp_change = min(prp_change, na.rm = TRUE),
            max_abs_prp_change = max(prp_change, na.rm = TRUE))
   
  


color_pal <- c("#543005", "#8c510a", "#bf812d", 
               "white",
               "#35978f", "#01665e", "#003c30")

color_pal_se <- c("#f7f4f9", "#bfd3e6", '#88419d','#810f7c','#4d004b')

deltas_se <- deltas_df %>%
  filter(metric == 'within_scale_red' & yr_pair %in% c("2011 - 2031", "2011 - 2091")) %>%
  mutate(
    se_value = case_when(
      yr_pair == "2011 - 2031" ~ se_value_2031,
      yr_pair == "2011 - 2091" ~ se_value_2091,
      TRUE ~ NA_real_
    )
  )

# Define the function to create the plot
# Define the function to create the plot
make_map_se <- function(data, variable, var_titles, limt, ecrg) {
  # Generate title and subtitle
  
  metric <- as.character(unique(data$metric)[1])
  title <- var_titles[[metric]]

  # Get the min and max proportional changes for the metric
  
  max_value <- max(data[[variable]], na.rm = TRUE)

  # Create the plot
  p <- ggplot(data = data) +
    geom_tile(aes(x = long, y = lat, fill = !!sym(variable))) +
    scale_fill_gradientn(colors = color_pal_se,
                         limits =  c(0, max_value),
                         na.value = 'grey80',
                         guide = guide_colorbar(barwidth = 10, barheight = 1)) +
    geom_sf(data = limt, fill = NA, col = '#36454F') +
    geom_sf(data = ecrg, fill = NA, col = '#36454F') +
    ggtitle(label = glue('SE of proportional change in {title}')) +
    theme_bw() +
    theme(legend.position = 'bottom', 
          legend.key.width = unit(2, 'line'),
          plot.title = element_text(size = 16, face = 'bold', hjust = 0, vjust = 0.7),
          plot.subtitle = element_text(size = 14, hjust = 0, vjust = 0.5),
          axis.title = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 12, face = 'bold'),
          strip.text = element_text(size = 14)) +
    labs(x = 'Longitude', y = 'Latitude', fill = 'SE of proportional change') +
    facet_wrap(~yr_pair, labeller = labeller(yr_pair = as.character))
  
  # Save the plot object
  out_dir <- './figs/Redundancy/delta/2025'
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  delta_filename <- glue('{out_dir}/{metric}_{variable}_SE.png')
  ggsave(plot = p, filename = delta_filename, units = 'in', width = 6, height = 6, dpi = 300)
  
  return(p)
}

# Example of creating a map for the delta_mean of cross_scale_red for multiple year pairs
var_titles <- list(
  'cross_scale_red' = 'Cross-scale redundancy',
  'within_scale_red' = 'Within-scale redundancy',
  'cross_scale_div' = 'Cross-scale diversity'
)
metric <- metrics[2]

delta_data <- deltas_df %>%
  filter(metric == 'within_scale_red' )
                                  
make_map_se(deltas_se, variable = "se_value", var_titles = var_titles,
         limt = limt, ecrg = ecrg)


#put maps together

# Load magick
library(magick)

# Define image paths
#path1 <- './figs/Redundancy/delta/2025/cross_scale_red_prp_change.png'
#path2 <- './figs/Redundancy/delta/2025/within_scale_red_prp_change.png'
#path3 <- './figs/Redundancy/delta/2025/cross_scale_div_prp_change.png'


path1 <- './figs/Redundancy/delta/2025/cross_scale_red_se_value_SE.png'
path2 <- './figs/Redundancy/delta/2025/within_scale_red_se_value_SE.png'
path3 <- './figs/Redundancy/delta/2025/cross_scale_div_se_value_SE.png'

# Read images
img1 <- image_read(path1)
img2 <- image_read(path2)
img3 <- image_read(path3)

# Annotate with a), b), c)
img1_annot <- image_annotate(img1, "a)", size = 40, gravity = "northwest", 
                             location = "+10+10", color = "black", weight = 700)
img2_annot <- image_annotate(img2, "b)", size = 40, gravity = "northwest", 
                             location = "+10+10", color = "black", weight = 700)
img3_annot <- image_annotate(img3, "c)", size = 40, gravity = "northwest", 
                             location = "+10+10", color = "black", weight = 700)

# Combine images side by side
combined <- image_append(c(img1_annot, img2_annot, img3_annot))

# Save final combined image
# Set 300 dpi
# Save final image
output_path <- './figs/Redundancy/delta/2025/redundancy_SE_combined_labeled_300dpi2.png'
image_write(combined, path = output_path, density = "300x300")


make_map <- function(data, variable, var_titles, limt, ecrg, min_max_abs_values) {
  # Generate title and subtitle
  
  metric <- as.character(unique(data$metric)[1])
  title <- var_titles[[metric]]
  
  # Get the min and max proportional changes for the metric
  min_max_value <- min_max_abs_values %>% filter(metric == !!metric)
  min_abs_value <- min_max_value %>% pull(min_abs_prp_change)
  max_abs_value <- min_max_value %>% pull(max_abs_prp_change)
  
  
  # Create the plot
  p <- ggplot(data = data) +
    geom_tile(aes(x = long, y = lat, fill = !!sym(variable))) +
    scale_fill_gradientn(colors = color_pal,
                         values = scales::rescale(c(-max_abs_value, 0, max_abs_value)),
                         limits = c(-max_abs_value, max_abs_value),
                         na.value = 'grey80',
                         guide = guide_colorbar(barwidth = 10, barheight = 1)) +
    geom_sf(data = limt, fill = NA, col = '#36454F') +
    geom_sf(data = ecrg, fill = NA, col = '#36454F') +
    ggtitle(label = glue('Proportional change in {title}')) +
    theme_bw() +
    theme(legend.position = 'bottom', 
          legend.key.width = unit(2, 'line'),
          plot.title = element_text(size = 16, face = 'bold', hjust = 0, vjust = 0.7),
          plot.subtitle = element_text(size = 14, hjust = 0, vjust = 0.5),
          axis.title = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 12, face = 'bold'),
          strip.text = element_text(size = 14)) +
    labs(x = 'Longitude', y = 'Latitude', fill = 'Mean proportional change') +
    facet_wrap(~yr_pair, labeller = labeller(yr_pair = as.character))
  
  # Save the plot object
  out_dir <- './figs/Redundancy/delta/2025'
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  delta_filename <- glue('{out_dir}/{metric}_{variable}.png')
  ggsave(plot = p, filename = delta_filename, units = 'in', width = 6, height = 6, dpi = 300)
  
  return(p)
}