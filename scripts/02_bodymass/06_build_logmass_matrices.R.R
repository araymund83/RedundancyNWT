# Load libraries ----------------------------
require(pacman)
pacman::p_load(dplyr,glue,qs,tidyverse)

g <- gc(reset = TRUE)
rm(list = ls())

# Load data -----------------------------------------
data <- qs::qread('./outputs/listBinCommMatrices_1191_1km.qs')
yrs <- c('2011', '2031','2091')
yr1 <- yrs[1]
yr2 <- yrs[3]
yr <- yr2
gcms <- c('CanESM2', 'CCSM4', 'INM-CM4')
log10Mass <-qs::qread('./inputs/traitsMassFG_2023.qs') 
log10Mass<- log10Mass %>% select(Masslog10)

##extract data frames for yr1 (before) and yr2 (after)
df1 <- map(1:length(data), function(i) data[[i]][[1]])
names(df1) <- glue('{gcms}_{yr1}')
df2 <- map(1:length(data), function(i) data[[i]][[2]])
names(df2) <- glue('{gcms}_{yr2}')

species <- colnames(df1[[1]])

get_bodyMass<- function(species, mass, gcm, df , yr){
  bodyMassDF <- sapply(species, function(sp) {
    message(crayon::green('Making bodyMass matrix for:', sp))
  
    mass <- filter(log10Mass, rownames(log10Mass) == sp)
    df_sp <- df %>% select(all_of(sp))
    df_sp[df_sp[[sp]]== 1 , ]<- mass
    df_sp[df_sp[[sp]]!= mass , ]<- 0
    
    # Check if any NA values exist in the matrix
    if (any(is.na(df_sp))) {
      message(crayon::red('WARNING: NA values detected in the matrix for', sp))
    }

    return(df_sp)
  })

bodyMassDF <- bodyMassDF %>% as.data.frame() %>%  mutate(yr = yr, 
                                                         gcm = gcm,
                                                         pixelID = glue('S_{row_number()}'))
colnames(bodyMassDF)[1:72] <- species
out <- './tables/logMassMat'

 ifelse(!dir.exists(out), dir.create(out), print ('Directory already exists'))
 qs::qsave (bodyMassDF, glue('{out}/logMassMatrix_{gcm}_{yr}.qs'))
 return(bodyMassDF)
}

result <- get_bodyMass(species, log10Mass, gcm = gcms[[3]], df2[[3]], yr = yr)



                   