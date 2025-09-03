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
yr3 <- yrs[3]
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

data_dir <- './tables/tbiFGBin/2025'
data31 <- qs::qread(glue("{data_dir}/TBIpvalFGBin_allgcms_{yr1}_{yr2}.qs"))
data91 <- qs::qread(glue("{data_dir}/TBIpvalFGBin_allgcms_{yr1}_{yr3}.qs"))

combined_tbiFG <- bind_rows(data31, data91)

#create a column for change 

combined_tbiFG <- combined_tbiFG %>%
  mutate(
    Change = str_trim(Change),  # Ensure consistent formatting
    change3 = case_when(
      Change == "-" ~ TBI * -1,
      Change == "+" ~ TBI,
      Change == "0"  ~ 0,
      TRUE ~ NA_real_
    )
  )
#save the file 
qs::qsave(combined_tbiFG, './tables/tbiFGBin/2025/tbixy_change_allgcms_allyrs.qs')

tbiBCDDF <- combined_tbiFG %>% 
  select(pixelID, Longitude, Latitude, gcm, tbi, p_val, yr, change3)



mean_tbiChange <- tbiBCDDF %>%
  group_by(Latitude, Longitude, yr) %>%
  summarise(mean_tbiChange = mean(change3, na.rm = TRUE))

# Create a new column for facet labels
mean_tbiChange <- mean_tbiChange %>%
  mutate(facet_label = case_when(
    yr == "2031" ~ "2011-2031",
    yr == "2091" ~ "2011-2091"
  ))


color_pal <- c('#543005','#8c510a','#bf812d','#dfc27d','white','#c7eae5','#80cdc1','#35978f', '#01665e', '#003c30')

ggTBI <-  ggplot() +
  geom_tile(data =mean_tbiChange, aes(x = Longitude, y = Latitude, fill = mean_tbiChange)) +
  scale_fill_gradientn(colours = color_pal) +
  #scale_color_brewer(type= 'seq', palette = 'PuBu') +
  geom_sf(data = limt, fill = NA, col = '#36454F') +
  geom_sf(data = ecrg,fill = NA,col = '#36454F') +
  ggtitle(label = glue('Mean TBI species change')) +
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
  labs(x = 'Longitude', y = 'Latitude', fill = 'mean TBI') +
  coord_sf(lims_method = 'geometry_bbox') +
  facet_wrap(~ yr)


out <-  glue('./figs/TBI_FG_changemean/2025')
ifelse(!dir.exists(out), dir.create(out, recursive = TRUE),
       print('Folder already exists'))
ggsave(plot = ggTBI, filename = glue('{out}/TBI_mean_Change_Bin_{yr2}_{yr3}.png'), 
       units = 'in', width = 10, height = 7, dpi = 300) 


tbi_sig_count <- combined_tbiFG %>%
  mutate(
    significant = ifelse(p.TBI < 0.05, 1, 0),  # mark 1 if significant
    Change = str_trim(Change)
  ) %>%
  group_by(Latitude, Longitude, yr) %>%
  summarise(significant_count = sum(significant, na.rm = TRUE)) %>%
  ungroup()

color_pal <- c(
  "0" = "white",  # light lavender (0)
  "1" = "#BFD3E6",  # light pink-purple (1)
  "2" = "#88419D",  # medium purple (2)
  "3" = "#4D004B"   # dark purple (3)
)
ggSig <- ggplot() +
  geom_tile(data = tbi_sig_count, aes(x = Longitude, y = Latitude, fill = as.factor(significant_count))) +
  scale_fill_manual(values = color_pal, name = "Significant GCMs",
                    labels = c("0 GCMs", "1 GCM", "2 GCMs", "3 GCMs")) +
  geom_sf(data = limt, fill = NA, color = '#36454F') +
  geom_sf(data = ecrg, fill = NA, color = '#36454F') +
  ggtitle("Number of GCMs with significant TBI (p < 0.05)") +
  facet_wrap(~ yr, labeller = labeller(yr = c("2031" = "2011–2031", "2091" = "2011–2091"))) +
  theme_bw() +
  theme(
    legend.position = 'bottom',
    plot.title = element_text(size = 16, face = 'bold'),
    strip.text = element_text(size = 14),
    legend.key.width = unit(2, 'line'),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 12, face = 'bold')
  ) +
  labs(x = "Longitude", y = "Latitude", fill = "GCM Agreement") +
  coord_sf()


out <- './figs/TBI_FG_changemean/2025'
if (!dir.exists(out)) dir.create(out, recursive = TRUE)

ggsave(filename = glue("{out}/TBI_Significance_Count_Bin_{yr2}_{yr3}.png"),
       plot = ggSig, width = 10, height = 7, dpi = 300)
