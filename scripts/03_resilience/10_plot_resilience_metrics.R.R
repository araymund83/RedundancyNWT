# Load libraries ----------------------------------------------------------
library(pacman)

p_load(dplyr, glue, qs, sf, tidyverse)

# Clean environment  ------------------------------------------------------
g <- gc(reset = TRUE)
rm(list = ls())

# Load data ---------------------------------------------------------------
out_dir <- './outputs/Redundancy/2025'
xy_table  <- qs::qread('./inputs/xyLU.qs')
gcms <- c('CanESM2', 'CCSM4', 'INM-CM4')
gcm <- gcms[1]
yrs <- c('2011','2031','2091')
yr <- yrs[1]
files <- list.files(out_dir, pattern = 'scale',full.names = TRUE)
fls <- grep(gcm, files, value = TRUE)
fls <- grep(yr, fls, value = TRUE)
red_tbls <- lapply(fls, qread)

# Create a function to perform the transformation and join
transform_and_join <- function(file_list) {
   cross_scale_red <- qs::qread(file_list) %>%
    rownames_to_column('pixelID') %>%
    mutate(gcm = gcm, yr = yr)
  
  left_join(xy_table, cross_scale_red, by = 'pixelID')
}
#use map to apply the function to each file and combine the results 
combined_data_list <- map(fls, transform_and_join)

# Define the merging function
merge_tables <- function(df1, df2) {
  left_join(df1, df2, by = c("pixelID", "gcm", "yr", "lat", "long"))
}

# # Apply the function to each data frame in the list
red_metrics_df <- purrr::reduce(combined_data_list, merge_tables)
 
 #save dataframe 
qs::qsave(red_metrics_df, glue('./outputs/Redundancy/2025/all_redMetrics_{gcm}_{yr}.qs'))

 #use pivot longer to convert to long format
red_metricsDF <- red_metrics_df %>% 
   pivot_longer(cols = starts_with('cross_scale') | starts_with('within_scale'),
                names_to = 'metric',
                values_to = 'value')
 
#save dataframe 
qsave(red_metricsDF, glue('./outputs/Redundancy/2025/allLong_redMetrics_{gcm}_{yr}.qs'))

#this code joins all the tables to one big data frame. 
 CanESM2011 <- qs::qread('./outputs/Redundancy/2025/allLong_redMetrics_CanESM2_2011.qs')
 CanESM2031 <- qs::qread('./outputs/Redundancy/allLong_redMetrics_CanESM2_2031.qs')
 CanESM2091 <- qs::qread('./outputs/Redundancy/allLong_redMetrics_CanESM2_2091.qs')
 all_metrics_df <- bind_rows(CanESM2011, CanESM2031, CanESM2091) %>% 
 na.omit()
# 
 qs::qsave(all_metrics_df, glue('./outputs/Redundancy/allMetrics_{gcm}_allyrs.qs'))
 data <- qs::qread(glue('./outputs/Redundancy/allMetrics_{gcm}_allyrs.qs'))
 
 
 read_and_merge_gcm_metrics <- function(gcm, years = c("2011", "2031", "2091"), 
                                        input_dir = "./outputs/Redundancy/2025", 
                                        save_output = TRUE) {
   # Construct full paths to files
   file_paths <- glue("{input_dir}/allLong_redMetrics_{gcm}_{years}.qs")
   
   # Read files
   message(glue("Reading files for {gcm}..."))
   data_list <- map(file_paths, qs::qread)
   
   # Combine and remove NA
   all_metrics_df <- bind_rows(data_list) %>%
     na.omit()
   
   # Save output
   if (save_output) {
     output_path <- glue("{input_dir}/allMetrics_{gcm}_allyrs.qs")
     qs::qsave(all_metrics_df, output_path)
     message(glue("Saved combined file to: {output_path}"))
   }
   
   return(all_metrics_df)
 } 
 
 # Read and save all metrics for INM-CM4
 data_CanESM2 <- read_and_merge_gcm_metrics(gcm = "CanESM2")
 
 
 
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

# Function to create maps of redundancy metrics ---------------------------
make_ggVar <- function(data, gcm, vars){
  ggplot_list <- lapply(vars, function(fill_var){
    
     #Determine the title based on 'fill_var'
    title <-  switch(fill_var,
             'cross_scale_red' = 'Cross-scale redundancy',
             'within_scale_red' = 'Within-scale redundancy',
             'cross_scale_div' = 'Cross-scale diversity',
      ) 
  
    #define fill variable dinamically based on 'var' argument
    var_data <- data %>% filter(metric %in% fill_var)

  ggplot() +
  geom_tile(data = var_data, aes(x = long, y = lat, fill = value))+
  scale_fill_gradientn(colours = color_pal) +
  geom_sf(data = limt, fill = NA, col = '#36454F') +
  geom_sf(data = ecrg,fill = NA,col = '#36454F') +
  ggtitle(label = glue('{title}'),
          subtitle = gcm) +
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
  labs(x = 'Longitude', y = 'Latitude', fill = fill_var) +
  coord_sf(lims_method = 'geometry_bbox') +
    facet_wrap(~ factoryr)
  })
  return(ggplot_list)
}

# Call the function  ------------------------------------------------------
result_plots <- make_ggVar(data = data, gcm = gcm, vars = vars)

#Save each ggplot using lapply
out_dir <-  glue('./figs/redundancy/2025')
ifelse(!dir.exists(out_dir), dir.create(out_dir, recursive = TRUE),
       print('Folder already exists'))
#save each plot of the list
lapply(seq_along(result_plots),function(i){
  plot <-result_plots[[i]]
  var_name <- vars[[i]]
  filename <- glue('{out_dir}/{var_name}_{gcm}_{yr}_2.png')
  ggsave(plot = plot, filename = filename, 
       units = 'in', width = 10, height = 7, dpi = 300)
})
  

