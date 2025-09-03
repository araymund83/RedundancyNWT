
# Function to calculate cross-scale redundancy for each pixel.
get_cross_scale_redundancy <- function(aggregations_list, functional_groups_df) {
  
  cross_scale_red_vals <- sapply(seq_along(aggregations_list), function(i){
    # Print progress
    cat("Processing pixel", i, "of", length(aggregations_list), "\n")
    
    # Extract the species in the current aggregation
    agg_list <- aggregations_list[[i]]
    
    if (!is.data.frame(agg_list) || nrow(agg_list) == 0) {
      warning("agg_list is empty in pixel ", i, ". Skipping.")
      return(NA)
    }
    
    # Match species to functional groups 
    FG_match <- agg_list %>% group_by(aggs) %>% 
      mutate(FG = functional_groups_df$FG[match(species, rownames(functional_groups_df))])
    #FG_match <- functional_groups_df$FG[match(agg_list$sp, rownames(functional_groups_df))]
    # Calculate the number of aggregations for each functional group
    phi_values <- table(FG_match$FG, FG_match$aggs)
    
    #get the count for each aggregation per functional group
    phi_sum_per_FG <- rowSums(phi_values)
    
    # Calculate the mean of count_by_aggs
    cross_scale_red <- mean(phi_sum_per_FG)
    
    
    # Check if phi_values is empty
    if (length(phi_values) == 0) {
      return(NA)
    }
    
    return(cross_scale_red)
  })
  
  return(cross_scale_red_vals)
}

# Call the function 
result_cross_scale_redundancy <- get_cross_scale_redundancy(aggregations_list = aggs, 
                                                            functional_groups_df = traitsFG)
names(result_cross_scale_redundancy) <- names(aggs)

result_cross_scale_redundancy<- do.call(rbind, lapply(result_cross_scale_redundancy, data.frame))
colnames(result_cross_scale_redundancy) <- 'cross_scale_red'
out_dir <- './outputs/Redundancy/2025'

if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}

qs::qsave(result_cross_scale_redundancy, glue('{out_dir}/cross_scale_red_{gcm}_{yr}.qs'))
#result_cross_scale_redundacy <- qs::qread(glue('{out_dir}/cross_scale_red_{gcm}_{yr}.qs'))

###############################################################################
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
#combined_df <- as.data.frame(result_within_scale_redundancy)
colnames(result_within_scale_redundancy) <- 'within_scale_red'


qs::qsave(result_within_scale_redundancy, glue('{out_dir}/within_scale_red_{gcm}_{yr}.qs'))

###############################################################################
#CROSS-SCALE DIVERSITY
################################################################################
get_cross_scale_diversity <- function(i, functional_groups_df){

  # Extract the species in the current aggregation
  agg_list <- i
  
  if (!is.data.frame(agg_list) || nrow(agg_list) == 0) {
    warning("agg_list is empty in pixel ", i, ". Skipping.")
    return(NA)
  }
  
  # Match species to functional groups 
  agg_list <- agg_list %>%
    mutate(FG = functional_groups_df$FG[match(species, rownames(functional_groups_df))])
  
  agg_list_summary <- agg_list %>%
    group_by(aggs) %>%
    summarize(total_species = n_distinct(species))
  
  diversity <- agg_list %>%
    group_by(aggs, FG) %>%
    summarize(noSpp = n_distinct(species))  %>% 
    # Merge with total species summary
    left_join(agg_list_summary, by = "aggs") %>%
    # Calculate relative species richness of each functional group 
    mutate(proportion = noSpp / total_species) %>% 
    # Calculate shanon index 
    mutate(shanon = -(proportion * log(proportion)))
  
  # Aggregate to get the mean result for each aggregation
  mean_result <- diversity %>%
    group_by(aggs) %>%
    summarize(mean_result = mean(shanon, na.rm = TRUE))
  
  # Calculate the mean diversity value across aggs
  shanon <- mean(mean_result$mean_result, na.rm = TRUE)
  
  return(shanon)
}

#Combine into a single data frame
result_cross_scale_diversity <- get_cross_scale_diversity(i = aggs,
                                                          functional_groups_df = traitsFG)

names(result_cross_scale_diversity) <- names(aggs)
# Combine into a single data frame
combined_df <- as.data.frame(result_cross_scale_diversity)
colnames(combined_df) <- 'cross_scale_div'

out_dir <- './outputs/Redundancy/2025'

if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}

qs::qsave(combined_df, glue('{out_dir}/cross_scale_diversity_{gcm}_{yr}.qs'))

