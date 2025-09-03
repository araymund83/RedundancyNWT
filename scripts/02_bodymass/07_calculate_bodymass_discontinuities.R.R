# Load libraries ----------------------------------------------------------
library(pacman)
p_load(dplyr, future, future.apply, ggplot2, glue, qs, tidyverse)
# Clean environment -------------------------------------------------------
g <- gc(reset = TRUE)
rm(list = ls())

# Source functions -------------------------------------------------------
source('./R/functions_discontinuities.R')

# Load data  --------------------------------------------------------------
xy_table <- qs::qread('./inputs/xyPixelID.qs')
colnames(xy_table) <- c('long', 'lat', 'pixelID')
Data <- qs::qread('./tables/logMassMat/allLogMassXY_recoded.qs')
Data[Data == 0] <- NA
Data <- Data %>% filter(yr == 1, gcm == 1)
Data <- Data %>% mutate(pixelID = glue('S_{1:nrow(Data)}'))

yrs <- c('2011', '2031','2091')
yr <- yrs[1]

gcms <- c('CanESM2', 'CCSM4', 'INM-CM4')
gcm <- gcms[2]

plan(multisession, workers = 16)
#Perform parallel computation
# discont_list <- future_lapply(1:nrow(Data), function(i){
#  Discontinuity_Barichievy(Data = Data[i,], start.column = 6,
#                        end.column = 77, IS_NA = TRUE, Resolution = 4000,
#                        Sample.N = 1000, gcm = gcm, yr = yr)
# })
  
discont_list <-  Discontinuity_Barichievy(
    Data = Data, start.column = 6,
    end.column = 77, IS_NA = TRUE, Resolution = 4000,
    Sample.N = 1000, gcm = gcm, yr = yr, parallel = TRUE, 
    nworkers = 16)


discont_list <- do.call(rbind, discont_list)
rownames(discont_list) <- NULL

# discont_list <- Discontinuity_Barichievy(Data, start.column = 6,
#                            end.column = 77, IS_NA = TRUE, Resolution = 4000,
#                            Sample.N = 1000, gcm = gcm, yr = yr, parallel = TRUE, 
#                            nworkers = 2)

# Filter out NULL entries in the output (if any)
non_null_entries <- sapply(discont_list, function(x) !is.null(x))

system.time(discont_list<- Discontinuity_Barichievy(Data = Data, start.column = 6,
                              end.column = 77, IS_NA = TRUE, Resolution = 4000,
                              Sample.N = 1000, gcm = gcm, yr = yr))
# Convert the result to a list
discont_list <- flatten(discont_list)
qs::qsave(discont_list, glue('./outputs/Discontinuities/Disc_{gcm}_{yr}.qs'))
test<- qs::qread( './outputs/Discontinuities/Disc_INM-CM4_2031.qs')



## create aggregations for each GCM and year 
yrs <- c('2011', '2031', '2091')
yr <- yrs[3]
gcms <- c('CanESM2', 'CCSM4', 'INM-CM4')
gcm <- gcms[1]
out_dir <- './outputs/Discontinuities'
files <- list.files(out_dir, pattern =  '.qs', full.names = TRUE)
file <- grep(gcm, files, value = TRUE)
file<- grep(yr, file, value = TRUE)

file_list <- lapply(file, qs::qread)


data <-qs::qread(glue('{out_dir}/Disc_{gcm}_{yr}.qs'))
#data <- purrr::flatten(data)  # flatten unlist the nested list 

pwrtab <- read.csv('./R/GRI.constant.power.table.csv')

# Group discontinuity results into body mass aggregations

result <- system.time({
  aggs_list <- lapply(data, function(data_element) {
    Aggregations(data, Tab = pwrtab, SetThresh = TRUE, constantThresh = 0.5, 
                 gcm = gcm, yr = yr)
  })
})
system.time(aggs <- Aggregations(data, Tab = pwrtab, SetThresh = TRUE, 
                         constantThresh = 0.5, gcm = gcm, yr = yr))


qs::qsave(aggs, glue( './outputs/Aggregations/aggs_{gcm}_{yr}.qs'))


