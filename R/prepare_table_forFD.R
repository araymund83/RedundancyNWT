speciesFD_df <- function(df, gcm, yr){
 # message(crayon::green('Preparing site-species table for:' gcm))
 #browser()
  df_yr <- dplyr::select(df, starts_with(yr))
  df_yr <- df_yr %>%  mutate(pixelID = row_number())
  df_yr <- df_yr %>%   mutate(pixelID = paste0('S_', pixelID))
  colnames(df_yr) <- gsub(yr,'', colnames(df_yr))
  colnames(df_yr) <- gsub('_',' ', colnames(df_yr))
  df_yr <- as_tibble(df_yr)
  df_yr  <-  df_yr  %>% remove_rownames()
  df_yr <-  column_to_rownames( df_yr , var = 'pixelID')
  df_yr <- as.matrix( df_yr )
  qs::qsave(x = df_yr, file = glue('./tables/comm_10km/{gcm}_{yr}_10km_ssmatrix.qs'))
  return(df_yr)
  
}
