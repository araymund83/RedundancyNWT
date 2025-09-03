# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(crayon, data.table, dplyr, glue, progress, qs, stringr, tidyverse, vegan)

g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)
# Load data ---------------------------------------------------------------
path <- './tables/fg_abun/2025'
files <- list.files(path, pattern = 'fg', full.names = TRUE)
yrs <- c('2011', '2031','2091')
yr1 <- yrs[1]
yr2 <- yrs[2]
yr3 <- yrs[3]
gcms <- c('CanESM2', 'CCSM4', 'INM-CM4')
gcm <- gcms[2]

#read regions pixels sampled
dir_path <- './tables/sampling/5Regions'
regions_fls <- list.files(dir_path, pattern = '10percent')
#read each file using lapply
file_data <-lapply(regions_fls, function(file){
  qs::qread(file.path(dir_path, file))
})
#Function to read and process files
read_and_process_files <- function(gcm, yr1, yr2, yr3, path, pix_sample) {
  
  fg_df_yr1 <- qread(glue("{path}/fg_{gcm}_{yr1}.qs"))
  fg_df_yr2 <- qread(glue("{path}/fg_{gcm}_{yr2}.qs"))
  fg_df_yr3 <- qread(glue("{path}/fg_{gcm}_{yr3}.qs"))
  
  fg_xy_yr1 <- left_join( pix_sample,fg_df_yr1, by = "pixelID") %>% 
    mutate(yr = yr1, gcm = gcm)
  fg_xy_yr2 <- left_join(pix_sample,fg_df_yr2,  by = "pixelID") %>% 
    mutate(yr = yr2, gcm = gcm)
  fg_xy_yr3 <- left_join(pix_sample, fg_df_yr3, by = "pixelID") %>% 
    mutate(yr = yr3, gcm = gcm)
  
  return(list(fg_xy_yr1, fg_xy_yr2, fg_xy_yr3))
}

run_nmds <- function(region, pix_sample, gcm, yr1, yr2, yr3, path) {
  
  excluded_files <- character(0)  # Initialize an empty vector to store excluded file names
  
  if (length(pix_sample) == 0) {
    cat('Skipping empty pix_sample for region:', region, '\n')
    excluded_files <- c(excluded_files, region)  # Add the excluded region name to the vector
    return()
  }
  
  
  
  region_name <- strsplit(region, "_10percent_sampling\\.qs")[[1]][1]
  cat('Processing region:', region_name, '\n')
  
  
  # Read data for yr1, yr2 and yr3
  fg_xy_yrs <- read_and_process_files(gcm, yr1, yr2, yr3, path, pix_sample)
  
  # Combine the data frames for yr1, yr2 and yr3
  fg_df_yrs <- do.call(rbind, lapply(fg_xy_yrs, function(x) x))
  # Filter rows where sum is not equal to 0
  fg_df_yrs <- fg_df_yrs %>%
    filter(rowSums(select(., 7:16)) != 0) %>%
    mutate(total_sum = rowSums(select(., 7:16)))
  
  # Check if fg_df_yrs has fewer than 5 rows
  if (nrow(fg_df_yrs) < 5) {
    cat('Skipping region', region_name, 'due to insufficient data\n')
    return()
  }
  
  #create fg matrix 
  fg_df_sp <- fg_df_yrs %>% select(c(7:16)) 
  
  #create sites matrix
  fg_df_sites <- fg_df_yrs %>% select(c(1:3,17,18))
  
  # Replace hyphens with underscores in column names for the fg matrix 
  colnames(fg_df_sp) <- gsub("-", "_", colnames(fg_df_sp))
  
  #Bootstrapping and testing for differences between the groups
  fit <- adonis2(fg_df_sp ~ yr, data = fg_df_sites, permutations = 999, 
                 method = "bray")
  
  #Check assumption of homogeneity of multivariate dispersion
  cat('Performing permanova for:', region_name, '\n')
  distances_data <- vegdist(fg_df_sp)
  homo_test <-anova(betadisper(distances_data, fg_df_sites$yr))
  
  #Perform NMDS and measure the time
  nmds_time <- system.time({
    nmds_out <- metaMDS(fg_df_sp, distance = 'bray', k = 2)
  })
  
  # Print the time taken
  cat("Processing time:", nmds_time[[3]], "seconds\n")
  
  #save nmds results
  qs::qsave(nmds_out, glue('./outputs/nmds_results_2025/nmds_{region_name}_{gcm}_10percent_allyrs.qs'))
  qs::qsave(fg_df_sites, glue('./outputs/nmds_results_2025/sites_{region_name}_{gcm}_10_percent_allyrs.qs'))
  qs::qsave(fit, glue('./outputs/bootstrapping_2025/{region_name}_{gcm}_10_percent_allyrs.qs'))
  qs::qsave(homo_test, glue('./outputs/permanova_2025/{region_name}_{gcm}_10_percent_allyrs.qs'))
  
  # Save the excluded file names to a file
  if (length(excluded_files) > 0) {
    writeLines(excluded_files, file = "excluded_files.txt", append = TRUE)
  }
  
}


file_data_subset <- file_data[[5]]

# Use lapply to run NMDS for each region
nmds_results <- lapply(seq_along(file_data_subset), function(i) {
  region <- regions_fls[i]
  pix_sample <- file_data[[i]]
  
  run_nmds(region, pix_sample, gcm, yrs[1], yrs[2], yrs[3], path)
})
