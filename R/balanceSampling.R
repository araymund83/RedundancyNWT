# #https://dickbrus.github.io/SpatialSamplingwithR/BalancedSpreaded.html
 
# regions <- unique(pix_xy$ECO4_NAM_1)
# pix_sd <- pix_xy %>% filter(ECO4_NUM_1 == 1)
# n <- 1000
# pi <- sampling::inclusionprobabilities(pix_sd$ECO4_NUM_1, n)
# X <- cbind(pix_sd$lat, pix_sd$long)
# set.seed(314)
# units <- BalancedSampling::lpm1(pi, X)
# sample <- pix_sd[units, ]

sample_pix_by_region <- function(pix_xy, sample){
 #Filter out rows with NA values in the ECO4_NAM_1 column
  pix_xy <- pix_xy[!is.na(pix_xy$ECO3_NAM_1),] 
  regions <- unique(pix_xy$ECO3_NAM_1)
   
  
  sampled_pixels <- lapply(regions, function(region){
      pix_sd <- pix_xy %>% filter(ECO3_NAM_1 == region)
      
      if (nrow(pix_sd) >= sample){
        pi <- sampling::inclusionprobabilities(pix_sd$ECO3_NUM_1, as.numeric(sample))
        X <- cbind(pix_sd$lat, pix_sd$long)
        set.seed(314)
        units <- BalancedSampling::lpm1(pi, X)
        pix_sample<- pix_sd[units,]
        qs::qsave(pix_sample, glue('./tables/sampling/{region}_pix_sampling.qs'))
      } else {
        message(glue('Region {region} has fewer than {sample} pixels'))
        return(NA)
      }
    })
  return(sampled_pixels)
}


# Use the function --------------------------------------------------------
#sampled_pixels <- sample_pix_by_region(pix_xy = pix_xy, sample = 1000)

################################################################################ 
#########function to sample 10% of the pixels within each ecoregion#############
################################################################################

sample_percent_by_region <- function(pix_xy, percent_to_sample = 0.1){
  #Filter out rows with NA values in the ECO4_NAM_1 column
  
  pix_xy <- pix_xy[!is.na(pix_xy$ECO3_NAM_1),] 
  regions<- unique(pix_xy$ECO3_NAM_1)
           
  sampled_pix <- lapply( regions, function(region){

    pix_sd <- pix_xy %>% filter(ECO3_NAM_1 == region)
    pix_sd$pix <- as.numeric(sub('S_',"", pix_sd$pixelID))
    total_pix<- nrow(pix_sd)
    n_samples_per_region <- round(total_pix * percent_to_sample)
    
    if(total_pix >= n_samples_per_region && is.numeric(n_samples_per_region)){
      pi <- sampling::inclusionprobabilities(pix_sd$pix, as.numeric(n_samples_per_region))
      X <- cbind(pix_sd$lat , pix_sd$long)
      set.seed(314)
      units <- BalancedSampling::lpm1(pi, X)
      sampled_data <- pix_sd[units,]
      qs::qsave(sampled_data, glue('./tables/sampling/{region}_10percent_sampling.qs'))
      return(sampled_data)
    } else {
      message(glue('Region {region} has fewer than {sample} pixels'))
      return(NA)
    }
    
  })
  sampled_pix <- do.call(rbind, sampled_pix)
  return(sampled_pix)
}
