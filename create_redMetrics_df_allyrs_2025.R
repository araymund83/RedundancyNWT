dir_path <- ('./outputs/Redundancy/2025')

gcms <- c('CanESM2', 'CCSM4', 'INM-CM4')
gcm <- gcms[3]
yrs <- c('2011','2031','2091')
yr <- yrs[3]
# Read all redundancy metrics 
all_fls <- list.files(dir_path, pattern = 'CanESM2', full.names = TRUE)

csr <- qs::qread(glue('{dir_path}/cross_scale_red_{gcm}_{yr}.qs'))
wir <- qs::qread(glue('{dir_path}/within_scale_red_{gcm}_{yr}.qs'))
csd <- qs::qread(glue('{dir_path}/cross_scale_diversity_{gcm}_{yr}.qs'))


# Step 1: Convert all to data frames and assign pixelID by row number
csr_df <- csr %>%
  rownames_to_column("pixelID") %>%
  rename(cross_scale_red = 2)

wir_df <- wir %>%
  rownames_to_column("pixelID") %>%
  rename(within_scale_red = 2)

# csd has no matching rownames, so assign pixelIDs based on csr order
# csd_df <- csd %>%
#   as.data.frame() %>%
#   mutate(pixelID = csr_df$pixelID) %>%   # <- use pixelID from csr
#   rename(cross_scale_div = 1)

 csd_df <- csd %>%
   rownames_to_column("pixelID") %>%
   rename(cross_scale_div = 2)

#csd_df <- csd %>%
 # rename(pixelID = pixel_ID)

# Step 2: Combine all
combined_df <- csr_df %>%
  left_join(wir_df, by = "pixelID") %>%
  left_join(csd_df, by = "pixelID") %>%
  mutate(gcm = gcm, yr = as.integer(yr))

qs::qsave(combined_df, glue('{dir_path}/allMetrics_{gcm}_{yr}.qs'))

INM.CM4_2011 <- qs::qread(glue('{dir_path}/allMetrics_{gcm}_2011.qs'))
INM.CM4_2031 <- qs::qread(glue('{dir_path}/allMetrics_{gcm}_2031.qs'))
INM.CM4_2091 <- qs::qread(glue('{dir_path}/allMetrics_{gcm}_2091.qs'))


# Step 2: Combine into a single wide-format data frame
combined_df <- bind_rows(INM.CM4_2011, INM.CM4_2031, INM.CM4_2091)
qs::qsave(combined_df, glue('{dir_path}/allMetrics_{gcm}_allyrs.qs'))

# Step 3: Pivot to long format
# Assumes these columns exist: pixelID, gcm, yr, cross_scale_red, within_scale_red, cross_scale_div
long_df <- combined_df %>%
  pivot_longer(
    cols = c(cross_scale_red, within_scale_red, cross_scale_div),
    names_to = "metric",
    values_to = "value"
  )

qs::qsave(long_df, glue('{dir_path}/allLong_redMetrics_{gcm}_allyrs.qs'))
