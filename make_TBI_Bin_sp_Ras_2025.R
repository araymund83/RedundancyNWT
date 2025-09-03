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
path <- './outputs/tbi_1km_Bin'
files <- list.files(path, pattern = 'tbiBin', full.names = TRUE)
yrs <- c('2011', '2031','2091')
yr1 <- yrs[1]
yr2 <- yrs[3]
gcms <- c('CanESM2', 'CCSM4', 'INM-CM4')

fls <- grep(yr2, files, value = TRUE)

#load shapefiles and CRS for maps. 
targetCRS <- paste("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95",
                   "+x_0=0 +y_0=0 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")


limt <- sf::st_read('inputs/NT1_BCR6/NT1_BCR6_poly.shp') 
ecrg <- sf::st_read('inputs/ecoregions/NWT_ecoregions_dissolvec.shp')

# Extract by mask for the ecoregions 
plot(st_geometry(ecrg))
limt <- sf::st_transform(x = limt, crs = targetCRS)
ecrg <- sf::st_transform(x = ecrg, crs = targetCRS)
plot(st_geometry(ecrg))
plot(st_geometry(limt))

##read the TBI outputs 
read_TBI <- function(gcms, yr2, files){
  map(gcms, ~ {
    gcm <- .
    message(crayon::green('Loading TBI results for GCM:', gcm))
    fls <- grep(gcm, files, value = TRUE)
    fl1 <- grep(yr2, fls, value = TRUE)
    df1 <- qs::qread(fl1)
    return(df1)
  })
}

# Apply the function ------------------------------------------------------
data <- read_TBI(gcms, yr2, fls)

data_2031 <- read_TBI(gcms, '2031', files)

# For 2091
data_2091 <- read_TBI(gcms, '2091', files)

make_TBI_Bin_Ras <- function(data, gcms, yr, var){
  # year label for filenames
  yr_label <- yr
  
  # Extract TBI values and p-values
  p_TBI <- lapply(data, function(x) x[["p.TBI"]])
  valsTBI <- lapply(data, function(x) x[["TBI"]])
  names(p_TBI) <- gcms
  names(valsTBI) <- gcms
  
  # Convert to data frames
  p_TBIDF <- data.frame(do.call(cbind, p_TBI)) %>%
    mutate(pixelID = glue('S_{1:nrow(.)}'))
  tbiDF <- data.frame(do.call(cbind, valsTBI)) %>%
    mutate(pixelID = glue('S_{1:nrow(.)}'))
  
  # Add coordinates
  p_TBIDF <- left_join(p_TBIDF, pixel_yx, by = 'pixelID')
  tbiDF <- left_join(tbiDF, pixel_yx, by = 'pixelID')
  
  # Pivot longer
  tbiDF <- tbiDF %>% pivot_longer(cols = 1:3, names_to = 'gcm', values_to = 'tbi')
  p_TBIDF <- p_TBIDF %>% pivot_longer(cols = 1:3, names_to = 'gcm', values_to = 'p_val')
  
  # Join both
  rasDF <- left_join(tbiDF, p_TBIDF, by = c("pixelID", "x", "y", "gcm"))
  sig_cutOff <- 0.05
  rasDF <- rasDF %>%
    mutate(yr = yr,
           gcm = gsub('INM.CM4', 'INM-CM4', gcm),
           pValue2 = case_when(p_val >= sig_cutOff ~ 'no_sig',
                               p_val < sig_cutOff ~ 'significant')) %>%
    rename(Longitude = x, Latitude = y)
  
  # Save raw table
  out <- glue('./tables/tbi1kmBin/2025')
  if (!dir.exists(out)) dir.create(out, recursive = TRUE)
  qs::qsave(rasDF, glue('{out}/TBIpvalBin_allgcms_{yr_label}.qs'))
  
  # MEAN TBI ---------------------------------------------------------------
  if (var == 'mean_tbi') {
    color_pal <- c('#f7fcfd','#bfd3e6','#8c96c6',
                   '#88419d','#810f7c','#4d004b')
    
    tbi_mean <- rasDF %>%
      group_by(pixelID, Longitude, Latitude) %>%
      summarise(mean_tbi = mean(tbi, na.rm = TRUE), .groups = "drop")
    
    ggMeanTBI <- ggplot() +
      geom_tile(data = tbi_mean, aes(x = Longitude, y = Latitude, fill = mean_tbi)) +
      scale_fill_gradientn(colours = color_pal, limits = c(0, 1), na.value = 'grey80') +
      geom_sf(data = limt, fill = NA, col = '#36454F') +
      geom_sf(data = ecrg, fill = NA, col = '#36454F') +
      ggtitle('Mean Species TBI', subtitle = glue('2011-{yr_label}')) +
      theme_bw() +
      theme(legend.position = 'bottom',
            legend.key.width = unit(2, 'line'),
            plot.title = element_text(size = 16, face = 'bold', hjust = 0, vjust = 0.7),
            plot.subtitle = element_text(size = 14)) +
      labs(x = 'Longitude', y = 'Latitude', fill = 'Mean TBI') +
      coord_sf()
    
    out_fig <- glue('./figs/TBI_bin_Ras_72spp/2025')
    if (!dir.exists(out_fig)) dir.create(out_fig, recursive = TRUE)
    ggsave(ggMeanTBI, filename = glue('{out_fig}/MeanTBI_Bin_{yr_label}.png'),
           width = 5, height = 7, dpi = 300)
  }
  
  # P-VALUE2 COUNT (MAP) ---------------------------------------------------
  if (var == 'pValue2_count') {
    sig_count <- rasDF %>%
      group_by(pixelID, Longitude, Latitude) %>%
      summarise(n_significant = sum(pValue2 == "significant"), .groups = "drop")
    
    color_pal <- c("0" = "white",
                   "1" = "#bfd3e6",
                   "2" = "#88419d",
                   "3" = "#4d004b")
    
    ggSigCount <- ggplot() +
      geom_tile(data = sig_count, aes(x = Longitude, y = Latitude, fill = factor(n_significant))) +
      scale_fill_manual(values = color_pal,
                        breaks = c("0","1","2","3"),
                        labels = c("0","1","2","3")) +
      geom_sf(data = limt, fill = NA, col = '#36454F') +
      geom_sf(data = ecrg, fill = NA, col = '#36454F') +
      ggtitle('Number of GCMs with significant TBI\n (p<0.05)',
              subtitle = glue('2011-{yr_label}')) +
      theme_bw() +
      theme(legend.position = 'bottom',
            legend.key.width = unit(2, 'line'),
            plot.title = element_text(size = 14, face = 'bold', hjust = 0, vjust = 0.7),
            plot.subtitle = element_text(size = 14)) +
      labs(x = 'Longitude', y = 'Latitude', fill = 'Numbero of GCMs') +
      coord_sf()
    
    out_fig <- glue('./figs/TBI_bin_Ras_72spp/2025')
    if (!dir.exists(out_fig)) dir.create(out_fig, recursive = TRUE)
    ggsave(ggSigCount, filename = glue('{out_fig}/SigCount_Bin_{yr_label}.png'),
           width = 5, height = 7, dpi = 300)
  }
}

# Load your lists of data first
data_2031 <- read_TBI(gcms, '2031', files)
data_2091 <- read_TBI(gcms, '2091', files)

# Mean TBI maps
make_TBI_Bin_Ras(data = data_2031, gcms = gcms, yr = "2031", var = "mean_tbi")
make_TBI_Bin_Ras(data = data_2091, gcms = gcms, yr = "2091", var = "mean_tbi")

# P-value significance count maps
make_TBI_Bin_Ras(data = data_2031, gcms = gcms, yr = "2031", var = "pValue2_count")
make_TBI_Bin_Ras(data = data_2091, gcms = gcms, yr = "2091", var = "pValue2_count")

# Combine 4 panels into one (2x2) with labels ---------------------------------
library(magick)
library(glue)

out_dir <- "./figs/TBI_bin_Ras_72spp/2025"

# Edit these if your filenames differ
p_mean2031 <- glue("{out_dir}/MeanTBI_Bin_2031.png")
p_mean2091 <- glue("{out_dir}/MeanTBI_Bin_2091.png")
p_sig2031  <- glue("{out_dir}/SigCount_Bin_2031.png")
p_sig2091  <- glue("{out_dir}/SigCount_Bin_2091.png")

stopifnot(file.exists(p_mean2031), file.exists(p_mean2091),
          file.exists(p_sig2031),  file.exists(p_sig2091))

# Read images
img_mean31 <- image_read(p_mean2031)
img_mean91 <- image_read(p_mean2091)
img_sig31  <- image_read(p_sig2031)
img_sig91  <- image_read(p_sig2091)
# Make sizes consistent (choose a target width that fits your output needs)
# Scale to same width
w <- 2200
img_mean31 <- image_scale(img_mean31, w)
img_mean91 <- image_scale(img_mean91, w)
img_sig31  <- image_scale(img_sig31,  w)
img_sig91  <- image_scale(img_sig91,  w)

# Build rows first
row_top <- image_append(c(img_mean31, img_mean91))   # Mean TBI
row_bot <- image_append(c(img_sig31,  img_sig91))    # Sig GCMs

# Add ONE label per row (on the left panel)
row_top <- image_annotate(row_top, "a)", size = 80, weight = 700,
                          gravity = "northwest", location = "+36+24", color = "black")
row_bot <- image_annotate(row_bot, "b)", size = 80, weight = 700,
                          gravity = "northwest", location = "+36+24", color = "black")

# Stack rows
final <- image_append(c(row_top, row_bot), stack = TRUE)
final <- image_border(final, color = "white", geometry = "10x10")

out_file <- glue("{out_dir}/Composite_MeanTBI_SigCount_2x2.png")
image_write(final, out_file)
message(glue("Saved composite: {out_file}"))
