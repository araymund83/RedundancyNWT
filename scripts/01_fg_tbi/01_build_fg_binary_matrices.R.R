# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(adespatial,crayon,dplyr,furrr, glue, future.apply, purrr, stringr,
               tictoc,tidyverse)

g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)
# Load data ---------------------------------------------------------------
path <- './tables/pres_absMat/'
files <- list.files(path, pattern = 'presAbs', full.names = TRUE)
yrs <- c('2011', '2031','2091')
yr1 <- yrs[1]
yr2 <- yrs[1]
gcms <- c('CanESM2', 'CCSM4', 'INM-CM4')
gcm <- gcms[3]


data <- qs::qread('./outputs/listBinCommMatrices_1191_1km.qs')
traitsFG <- qs::qread('./inputs/traitsMassFG_2025.qs')
traitsFG <- traitsFG %>% select(FG) %>% 
  rownames_to_column('species')

# Get Temporal Beta diversity Index (TBI) ---------------------------------
#extract data frames for yr1 (before) and yr2 (after)
df1 <- map(1:length(data), function(i) data[[i]][[1]])
names(df1) <- glue('{gcms}_{yr1}')
df2 <- map(1:length(data), function(i) data[[i]][[2]])
names(df2) <- glue('{gcms}_{yr2}')

data1<- df1[[3]] ## !!!change to select the gcm 1: CanESM2, 2:CCSM4, 3:INM-CM4

# Function to add a unique identifier to duplicate names
make_unique <- function(names) {
  count <- table(names)
  ave(names, names, FUN = function(x) ifelse(count[x] > 1, paste0(x, seq_along(x)), x))
}

# Replace column names in df1 with values from FG column in df2
test <- data1 %>%
  rename_with(~make_unique(traitsFG$FG[match(., traitsFG$species)]), everything())


# Function to add total column for a specific FG value
add_total_column <- function(data, fg) {
  df <- data[, grep(paste0("^", fg), colnames(data), value = TRUE), drop = FALSE]
  df$total <- rowSums(df)
  return(df)
}

# Get unique FG values from column names
FG_values <- unique(gsub("\\d", "", colnames(test)))

# Use lapply to apply the function to each FG value
grouped_dfs <- lapply(FG_values, function(fg) add_total_column(test, fg))

# Access the data frame for a specific FG value (e.g., "INS-Sal")
#df_ins_sal <- grouped_dfs[["INS-Sal"]]

# Extract the 'total' column from each data frame in grouped_dfs
total_columns <- lapply(grouped_dfs, function(df) df$total)

# Create a data frame with total columns
total_df <- as.data.frame(total_columns)

# Rename the columns of the total_df based on the original FG values
colnames(total_df) <- FG_values

total_df <- total_df %>% mutate(pixelID = glue('S_{1:nrow(total_df)}'))

out_dir <- './tables/fg_abun/2025'

if (!dir.exists(out_dir)) {dir.create(out_dir, recursive = TRUE)
} else {print('Folder already exists')}

qs::qsave(total_df, glue('{out_dir}/fg_{gcm}_{yr2}.qs'))

bin_df <- total_df %>% mutate(across(-pixelID, ~ifelse(. >=1, 1,0)))

out_dir <- './tables/fg_bin/2025'

if (!dir.exists(out_dir)) {dir.create(out_dir, recursive = TRUE)
} else {print('Folder already exists')}

qs::qsave(bin_df, glue('{out_dir}/fg_{gcm}_{yr2}.qs'))
