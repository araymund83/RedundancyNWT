# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(dplyr,glue, purrr, qs, stringr, sf, terra, tidyverse)

g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)

# Load data ---------------------------------------------------------------
#Load lat and long for each pixel
pixel_yx <- qs::qread('./inputs/xyPixelID.qs')
path <- glue('./tables/tbi1kmBin')
files <- list.files(path, pattern = 'qs', full.names = TRUE)
yrs <- c('2011', '2031','2091')
yr1 <- yrs[1]
yr2 <- yrs[2]
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

##read the BCD outputs 
read_BCD <- function(gcms, yr2, files){
  lapply(gcms, function(gcm){
    message(crayon::green('Loading bcd results for GCM:', gcm))
    fls <- grep('bcd', files, value = TRUE)
    fls <- grep(gcm, files, value = TRUE)
    fl1 <- grep(yr2, fls, value = TRUE)
    df1 <- qs::qread(fl1)
    return(df1)
  })
}

# Apply the function ------------------------------------------------------
data <- read_BCD(gcms, yr2, files)
names(data) <- gcms


tbiDF <- qs::qread( glue('./tables/tbi1kmBin/TBIpvalBin_allgcms_{yr1}_{yr2}.qs'))
tbiDF <- tbiDF %>% mutate(gcm = recode(gcm, 'INM.CM4' = 'INM-CM4'))
# Create a function for left join
left_join_gcm <- function(df, gcm_df) {
  df %>%
    filter(gcm  %in% unique(gcm_df$gc)) %>%
    left_join(gcm_df, by = c("pixelID", "gcm" = "gc"))
}

# List of GCM names
gcm_names <- names(data)

# Use lapply to perform left join for each GCM
joined_list <- map(gcm_names, ~left_join_gcm(tbiDF, data[[.x]]))


# Combine the list of data frames into a single dataframe
tbiBCDDF <- bind_rows(joined_list)
qs::qsave(tbiBCDDF, glue('./tables/tbi1kmBin/tbiBCD_allgcms_{yr1}_{yr2}.qs'))


data <- (qs::qread('./tables/tbiFGBin/2025/tbixy_change_allgcms_allyrs.qs'))
# Combine year range
year_range <- glue("{yr1}-{yr2}")

# Filter the data correctly
data_filtered <- data %>%
  filter(as.character(yrs) == year_range)
data_filtered <- data_filtered %>%
  mutate(pValue2 = factor(pValue2, levels = c("no_sig", "significant")))



#color_pal <- c('#f7fcfd','#bfd3e6','#8c96c6',
  #             '#88419d','#810f7c','#4d004b')
color_pal <- c('#543005','#8c510a','#bf812d','#dfc27d','#f6e8c3','#c7eae5','#80cdc1','#35978f', '#01665e', '#003c30')

ggTBI <-  ggplot() +
  geom_tile(data =data_filtered, aes(x = Longitude, y = Latitude, fill = change3)) +
  scale_fill_gradientn(colours = color_pal) +
  #scale_color_brewer(type= 'seq', palette = 'PuBu') +
  geom_sf(data = limt, fill = NA, col = '#36454F') +
  geom_sf(data = ecrg,fill = NA,col = '#36454F') +
  ggtitle(label = glue('TBI Functional Group Change'),
          subtitle = glue('{yr1}-{yr2}')) +
  theme_bw() +
  theme(legend.position = 'bottom',legend.key.width = unit(2, 'line'),
        plot.title = element_text(size = 16, face = 'bold', hjust = 0, vjust = 0.7),
        plot.subtitle = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12, face = 'bold'),
        strip.text = element_text(size = 14)) +
  labs(x = 'Longitude', y = 'Latitude', fill = 'TBI') +
  coord_sf(lims_method = 'geometry_bbox') +
  facet_wrap(~ gcm)

out <-  glue('./figs/TBI_bin_Ras_72spp/2025')
ifelse(!dir.exists(out), dir.create(out, recursive = TRUE),
       print('Folder already exists'))
ggsave(plot = ggTBI, filename = glue('{out}/TBIChange_FG_{yr1}_{yr2}.png'), 
       units = 'in', width = 10, height = 7, dpi = 300) 


ggPval <- ggplot() +
  geom_tile(data = filter(data_filtered, pValue2 == "no_sig"),
            aes(x = Longitude, y = Latitude, fill = pValue2)) +
  geom_tile(data = filter(data_filtered, pValue2 == "significant"),
            aes(x = Longitude, y = Latitude, fill = pValue2)) +
  scale_fill_manual(
    values = c("no_sig" = "grey80", "significant" = "#570E55"),
    name = "Number of GCMs",
    labels = c("no_sig" = "Not significant", "significant" = "Significant"),
    drop = FALSE
  ) +
  geom_sf(data = limt, fill = NA, col = '#36454F') +
  geom_sf(data = ecrg, fill = NA, col = '#36454F') +
  ggtitle(label = "p-value of TBI for Functional Groups",
          subtitle = glue("{yr1}-{yr2}")) +
  theme_bw() +
  theme(
    legend.position = 'bottom',
    legend.key.width = unit(2, 'line'),
    plot.title = element_text(size = 16, face = 'bold', hjust = 0),
    plot.subtitle = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 12, face = 'bold'),
    strip.text = element_text(size = 14)
  ) +
  labs(x = 'Longitude', y = 'Latitude', fill = 'Significance') +
  coord_sf(lims_method = 'geometry_bbox') +
  facet_wrap(~ gcm)



out <- glue('./figs/TBI_bin_Ras_72spp/2025')
if (!dir.exists(out)) dir.create(out, recursive = TRUE)

ggsave(plot = ggPval, 
       filename = glue('{out}/pValue2_FG_{yr1}_{yr2}.png'), 
       units = 'in', width = 10, height = 7, dpi = 300)

#put together tbi change and pvalue figures 
combined_plot <- ggTBI / ggPval +
  plot_annotation(tag_levels = 'a')  # labels: a) and b)

