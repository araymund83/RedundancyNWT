## Loading R packages
library(pacman)
p_load(data.table, dplyr, ggplot2,  glue,  tidyverse, future, future.apply)

# Clean environment  ------------------------------------------------------
g <- gc(reset = TRUE)
rm(list = ls())

# Source functions --------------------------------------------------------
source("R/functions_redundancymetrics.R")

# Load data ---------------------------------------------------------------
path <- './inputs'

traitsFG <- qs::qread('./inputs/traitsMassFG_2025.qs')
traitsFG <- traitsFG %>% select(FG) 
data <- qs::qread('./outputs/listCommMatricesBin_1191_1km.qs')
yrs <- c('2011', '2031','2091')
yr <- yrs[3]
gcms <- c('CanESM2', 'CCSM4', 'INM-CM4')
gcm <- gcms[3]

splu <- traits %>% select(species, Habitat) %>% 
  mutate(spID = as.character(1:nrow(traits)))

# #convert character traits to categorical traits
# traits <- traits %>% mutate_at(c(7:14), ~scale(.) %>% as.vector) %>% # normalize numerical traits 
#   mutate(Masslog10 = log10(Mass)) %>% 
#   mutate(across(c(2:6,16), as.factor)) %>% 
#   column_to_rownames('species')
# traits2 <- traits %>% select(2:6,7:13,16)
# write.csv(traits2, glue('{path}/birdTraitsStandard_2023.csv'))
# 
# traitsFG<- traits %>%
#            mutate(FG = glue('{Feeding}-{Forage_Breed1}')) %>% 
#            select(FG, Masslog10) 
traitsFG <- qs::qread('./inputs/traitsMassFG_2025.qs')
traitsFG <- traitsFG %>% select(FG)

#read aggregations list 
aggs <- qs::qread(glue('./outputs/Aggregations/aggs_{gcm}_{yr}.qs'))

get_cross_scale_redundancy <- function(aggregations_list, functional_groups_df) {
  future_sapply(seq_along(aggregations_list), function(i){
    agg_list <- aggregations_list[[i]]
    if (!is.data.frame(agg_list) || nrow(agg_list) == 0) return(NA)
    
    FG_match <- agg_list %>%
      group_by(aggs) %>%
      mutate(FG = functional_groups_df$FG[match(species, rownames(functional_groups_df))])
    
    phi_values <- table(FG_match$FG, FG_match$aggs)
    phi_sum_per_FG <- rowSums(phi_values)
    mean(phi_sum_per_FG)
  }, future.seed = TRUE)
}

# Call the function 
plan(multisession, workers = 10)
result_cross_scale_redundancy <- sapply(aggs, get_cross_scale_redundancy,
                                               functional_groups_df = traitsFG,
                                               USE.NAMES = TRUE)
future:::ClusterRegistry('stop')

combined_df <- as.data.frame(result_cross_scale_redundancy)
result_cross_scale_redundancy2<- do.call(rbind, lapply(result_cross_scale_redundancy, data.frame))

colnames(result_cross_scale_redundancy2) <- 'cross_scale_red'
out_dir <- './outputs/Redundancy/2025'

if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}

qs::qsave(result_cross_scale_redundancy, glue('{out_dir}/cross_scale_red_{gcm}_{yr}.qs'))


################################################################################
#WHITHIN SCALE REDUNDANCY
################################################################################

get_within_scale_redundancy <- function(aggregations_list, functional_groups_df) {
  ## TODO: bring the sapply to the Global Env.
  within_scale_redundancy_values <- sapply(seq_along(aggregations_list), function(i) {
    ## TODO: define this function in a separate script sourced at the start.
    # Print progress
    cat("Processing pixel", i, "of", length(aggregations_list), "\n")
    
    # Extract the species in the current aggregation
    agg_list <- aggregations_list[[i]]
    # Check if agg_list is NULL
   
    if (!is.data.frame(agg_list) || nrow(agg_list) == 0) {
      warning("agg_list is empty in pixel ", i, ". Skipping.")
      return(NA)
    }
   
      # Match species to functional groups 
      agg_list <- agg_list %>%
        mutate(FG = functional_groups_df$FG[match(species, rownames(functional_groups_df))])
      #agg_list$FG <- functional_groups_df[match(agg_list$sp, rownames(functional_groups_df)), , drop = FALSE]
      
      # calculate the number of species per functional group within an aggregation
      
      agg_list <- agg_list %>%
        group_by(FG, aggs) %>%
        summarise(num_sp = n_distinct(species), .groups = 'drop') %>%
        group_by(FG) %>%
        summarise(avg_num_sp = mean(num_sp)) 
      
      within_scale_red <- mean(agg_list$avg_num_sp)
      return(within_scale_red)      
    })
    return(within_scale_redundancy_values)   
 }

result_within_scale_redundancy <- get_within_scale_redundancy(aggregations_list = aggs, 
                                                              functional_groups_df = traitsFG)  
names(result_within_scale_redundancy) <- names(aggs)
result_within_scale_redundancy<- do.call(rbind, lapply(result_within_scale_redundancy, data.frame))
colnames(result_within_scale_redundancy) <- 'within_scale_red'
result_within_scale_redundancy <- data.frame(result_within_scale_redundancy)

qs::qsave(result_within_scale_redundancy, glue('{out_dir}/within_scale_red_{gcm}_{yr}.qs'))

# Function to calculate Shannon diversity index for a single data frame
plan(multisession, workers = 8)
result_cross_scale_diversity2 <- future_sapply(aggs, get_cross_scale_diversity2, 
                                               functional_groups_df = traitsFG,
                                               USE.NAMES = TRUE) # already with names
future:::ClusterRegistry("stop")

# Combine into a single data frame
combined_df <- as.data.frame(result_cross_scale_diversity2)
colnames(combined_df) <- 'cross_scale_div'

out_dir <- './outputs/Redundancy'

if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}

qs::qsave(combined_df, glue('{out_dir}/cross_scale_diversity_{gcm}_{yr}.qs'))

