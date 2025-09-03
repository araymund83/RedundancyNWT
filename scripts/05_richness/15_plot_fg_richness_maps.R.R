library(pacman)
p_load(glue, ggplot2, glue, qs, tidyverse)

# Session/options
g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)

# Load data ---------------------------------------------------------------
data <- qs::qread('./tables/richnessFG/richnessFG_allyrs_allgcm.qs')
pixel_yx <- qs::qread('./inputs/xyPixelID.qs')
yrs <- c('2011', '2031','2091')
gcms <- c('CanESM2', 'CCSM4', 'INM-CM4')
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

#Define FG columns present in fg_all
fg_cols <- c("INS-Sal","OMN-For","INS-Gle","INS-Exc","INS-For",
             "OMN-Sca","OMN-Gle","GRA-Gle","INS-Scr","OMN-Exc")
data <- data %>% select (- c(geometry, lat, long, ECO4_NAM_1))

pixel_yx <- pixel_yx %>% rename(long = x, lat = y)

make_ggFG_rich <- function(data, gcms, yrs,
                           fg_cols,
                           save = FALSE) {
  
  
  # Palette
  color_pal <- c('#f7fcfd','#bfd3e6','#8c96c6',
                 '#88419d','#810f7c','#4d004b')
  
  # Compute FG richness once
  data_rich <- data %>% 
    mutate(fg_rich = rowSums(across(all_of(fg_cols), ~ as.numeric(.x) > 0), na.rm = TRUE)) %>% 
    filter(yr %in% yrs)
  
  data_rich_xy <- data_rich %>%
    left_join(pixel_yx, by = "pixelID")
  
  # Build one plot per GCM (returns a named list)
  ggplot_list <- lapply(gcms, function(gcm_i) {
    
    var_data <- data_rich_xy |>
      filter(gcm == gcm_i) |>
      mutate(yr = factor(yr, levels = sort(unique(yr))))
    
  
    max_rich <- max(var_data$fg_rich, na.rm = TRUE)

    ggplot() +
      geom_tile(data = var_data, aes(x = long, y = lat, fill = fg_rich)) +
      scale_fill_gradientn(colours = color_pal, limits = c(0, 10),
                           oob = scales::squish, name = "FG richness") +
      geom_sf(data = limt, fill = NA, col = '#36454F') +
      geom_sf(data = ecrg, fill = NA, col = '#36454F') +
      ggtitle(label = "Functional-Group Richness", subtitle = gcm_i) +
      theme_bw() +
      theme(
        legend.position = 'bottom',
        legend.key.width = unit(2, 'line'),
        plot.title = element_text(size = 16, face = 'bold', hjust = 0, vjust = 0.7),
        plot.subtitle = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12, face = 'bold'),
        strip.text = element_text(size = 14)
      ) +
      labs(x = 'Latitude', y = 'Longitude') +
      coord_sf(lims_method = 'geometry_bbox') +
      facet_wrap(~ yr)
  })
  
  names(ggplot_list) <- gcms
  out_dir = "./figs/FG_richness/2025"
  
  if (isTRUE(save)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    invisible(lapply(names(ggplot_list), function(nm) {
      ggsave(filename = glue("{out_dir}/FG_richness_{nm}.png"),
             plot = ggplot_list[[nm]], units = 'in', width = 10, height = 7, dpi = 300)
    }))
  }
  
  return(ggplot_list)
}

# Call the function -------------------------------------------------------

plots <- make_ggFG_rich(data = data, yrs= yrs, fg_cols= fg_cols, gcms = gcms, save = TRUE)

# Directory with the saved PNGs
out_dir <- "./figs/FG_richness/2025"

# Read the images
img_can <- image_read(file.path(out_dir, "FG_richness_CanESM2.png"))
img_ccs <- image_read(file.path(out_dir, "FG_richness_CCSM4.png"))
img_inm <- image_read(file.path(out_dir, "FG_richness_INM-CM4.png"))

# Combine them horizontally
combo <- image_append(c(img_inm, img_ccs, img_can), stack = TRUE)

# Save the combined figure
out_file <- file.path(out_dir, "FG_richness_ALLGCMs.png")
image_write(combo, out_file)

# Show in RStudio viewer / plot window
print(combo)
