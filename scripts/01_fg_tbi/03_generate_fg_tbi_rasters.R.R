# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(dplyr,glue, purrr, qs, RColorBrewer, stringr, sf, terra, tidyverse)

g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)

# Load data ---------------------------------------------------------------
#Load lat and long for each pixel
pixel_yx <- qs::qread('./inputs/xyPixelID.qs')
#load files for extracting TBI
path <- './outputs/tbiFGBin_1km/2025'
files <- list.files(path, pattern = '72spp', full.names = TRUE)
yrs <- c('2011', '2031','2091')
yr1 <- yrs[1]
yr2 <- yrs[3]
gcms <- c('CanESM2', 'CCSM4', 'INM-CM4')
year_range <- glue::glue("{yr1}-{yr2}")
#gcm <- gcms[1]

#load shapefiles and CRS for maps. 
targetCRS <- paste("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95",
                   "+x_0=0 +y_0=0 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
limt <- sf::st_read('inputs/NT1_BCR6/NT1_BCR6_poly.shp') 
ecrg <- sf::st_read('inputs/ecoregions/NWT_ecoregions_5Regions.shp')

# Extract by mask for the ecoregions 
plot(st_geometry(ecrg))
limt <- sf::st_transform(x = limt, crs = targetCRS)
ecrg <- sf::st_transform(x = ecrg, crs = targetCRS)
plot(st_geometry(ecrg))
plot(st_geometry(limt))

##read the TBI outputs 
read_TBI <- function(gcms, yr2, files){
  map(gcms, ~ {gcm <- .
    message(crayon::green('Loading TBI results for GCM:', gcm))
    fls <- grep(gcm, files, value = TRUE)
    fl1 <- grep(yr2, fls, value = TRUE)
    df1 <- qs::qread(fl1)
    return(df1)
  })
}

# ---- Load TBI results ----
data <- read_TBI(gcms, yr2, files)
names(data) <- gcms

#function to extract the TBI, signficance and if the change was + or - element for each gcm
# ---- Combine all GCMs into one tibble ----
combined_results <- imap_dfr(data, function(x, gcm) {
  tibble(
    TBI = x$TBI,
    p.TBI = x$p.TBI,
    Change = str_trim(x$BCD.mat$Change),  # <-- trim whitespace
    gcm = gcm,
    yrs = year_range
  )
}) %>%
  mutate(
    change2 = case_when(
      Change == "+" ~ "positive",
      Change == "-" ~ "negative",
      Change == "0" ~ "no_change",
      TRUE ~ NA_character_
    ),
    change3 = case_when(
      Change == "-" ~ TBI * -1,
      Change == "+" ~ TBI,
      Change == "0" ~ 0,
      TRUE ~ NA_real_
    )
  )

# ---- Save combined results ----
output_dir <- './tables/tbiFGBin/2025'
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
qs::qsave(combined_results, file = glue('{output_dir}/tbiFG_2025_allgcms_{yr1}_{yr2}.qs'))

tbi_df <- qs::qread(glue('{output_dir}/tbiFG_2025_allgcms_{yr1}_{yr2}.qs'))



make_TBI_FGBin_Ras <- function(tbi_df, yr1, yr2, var) {
  message("Joining spatial coordinates...")
  
  n_pixels <- nrow(pixel_yx)
  pixel_ids <- glue::glue("S_{1:n_pixels}")
  
  rasDF <- tbi_df %>%
    mutate(
      pixelID = rep(pixel_ids, times = length(unique(gcm))),
      yr = glue("{yr1}-{yr2}"),
      pValue2 = ifelse(p.TBI < 0.05, "significant", "no_sig"),
      tbi = TBI,
      p_val = p.TBI
    ) %>%
    left_join(pixel_yx, by = "pixelID") %>%
    mutate(
      Longitude = x,
      Latitude = y,
      gcm = gsub("INM.CM4", "INM-CM4", gcm)
    )
  
  # Save output
  out <- glue('./tables/tbiFGBin/2025')
  if (!dir.exists(out)) dir.create(out, recursive = TRUE)
  qs::qsave(rasDF, glue('{out}/TBIpvalFGBin_allgcms_{yr1}_{yr2}.qs'))
  
  # Plot
  message(glue("Creating {var} map for {yr1}-{yr2}..."))
  
  if (var == 'pValue2') {
    color_2vals <- c('#b3b3b3', '#4d004b')
    ggVar <- ggplot() +
      geom_tile(data = rasDF, aes(x = Longitude, y = Latitude, fill = pValue2)) +
      scale_fill_manual(values = color_2vals,
                        labels = c('No significant', 'Significant')) +
      geom_sf(data = limt, fill = NA, col = '#36454F') +
      geom_sf(data = ecrg, fill = NA, col = '#36454F') +
      ggtitle('pValue TBI for Functional Groups', subtitle = glue('{yr1}-{yr2}')) +
      theme_bw() +
      theme(legend.position = 'bottom',
            legend.key.width = unit(2, 'line'),
            plot.title = element_text(size = 16, face = 'bold'),
            plot.subtitle = element_text(size = 14),
            axis.title = element_text(size = 14),
            axis.text = element_text(size = 12),
            legend.text = element_text(size = 11),
            legend.title = element_text(size = 12, face = 'bold'),
            strip.text = element_text(size = 14)) +
      labs(x = 'Longitude', y = 'Latitude', fill = 'pValue') +
      coord_sf() +
      facet_wrap(~ factor(gcm, c('INM-CM4', 'CCSM4', 'CanESM2')))
    
    out_plot <- glue('./figs/TBI_FG_Ras_72spp/2025')
    if (!dir.exists(out_plot)) dir.create(out_plot, recursive = TRUE)
    ggsave(plot = ggVar, filename = glue('{out_plot}/{var}_FG_{yr1}_{yr2}_2.png'),
           units = 'in', width = 10, height = 7, dpi = 300)
  } else if (var == 'tbi') {
    color_pal <- c('#f7fcfd','#bfd3e6','#8c96c6','#88419d','#810f7c','#4d004b')
    ggTBI <- ggplot() +
      geom_tile(data = rasDF, aes(x = Longitude, y = Latitude, fill = tbi)) +
      scale_fill_gradientn(colours = color_pal, limits = c(0, 1)) +
      geom_sf(data = limt, fill = NA, col = '#36454F') +
      geom_sf(data = ecrg, fill = NA, col = '#36454F') +
      ggtitle('TBI index for Functional Groups', subtitle = glue('{yr1}-{yr2}')) +
      theme_bw() +
      theme(legend.position = 'bottom',
            legend.key.width = unit(2, 'line'),
            plot.title = element_text(size = 16, face = 'bold'),
            plot.subtitle = element_text(size = 14),
            axis.title = element_text(size = 14),
            axis.text = element_text(size = 12),
            legend.text = element_text(size = 11),
            legend.title = element_text(size = 12, face = 'bold'),
            strip.text = element_text(size = 14)) +
      labs(x = 'Longitude', y = 'Latitude', fill = 'TBI') +
      coord_sf(lims_method = 'geometry_bbox') +
      facet_wrap(~ factor(gcm, c('INM-CM4', 'CCSM4', 'CanESM2')))
    
    out_plot <- glue('./figs/TBI_FG_Ras_72spp/2025')
    if (!dir.exists(out_plot)) dir.create(out_plot, recursive = TRUE)
    ggsave(plot = ggTBI, filename = glue('{out_plot}/{var}_FG_{yr1}_{yr2}.png'),
           units = 'in', width = 10, height = 7, dpi = 300)
  }
}

# ---- Run plotting function ----
make_TBI_FGBin_Ras(
  tbi_df = tbi_df,
  yr1 = yr1,
  yr2 = yr2,
  var = 'pValue2'  ##var options are 'tbi' and 'pValue2'
)



 
 

 

 
                   