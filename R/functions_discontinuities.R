#=============================================================================
## modification of Chris Barichievy's Neutral Null & Discontinuity Detector

##### Discontinuity_Barichievy()
# Function based on Chris Barichievy's "discontinuity detector" code from his
# Barichievy et al. (2018) manuscript. In this function, wide form data (e.g., 
# bbs_mass file) with body mass data over years and different locations are 
# processed. Then, the processed data is fed into "neutral null" and 
# "discontinuity detector" code that are wrapped together.
# 
# Arguments;
#   Data         = a data.frame-like object in wide form (see bbs_mass file)  
#                  with log-transformed body mass data by year, by location 
#                  (i.e., BBS route)
#   
#   start.column = column in Data that species body masses begin
#   end.column   = number of columns in Data
#   IS_NA        = a logical value indicating the form of Data; that is, if 
#                  a species does NOT occur at a route in a given year, is 
#                  this indicated by a "0" or an "NA"? If "NA",  
#                  IS_NA = TRUE. See bbs_mass file for example.
#   Resolution   = See code in Barichievy et al. (2018)
#   Sample.N     = See code in Barichievy et al. (2018)
Discontinuity_Barichievy <- function(...,
                                     parallel = FALSE, nworkers = 1L) {
  dots <- list(...)
  missing <- c("Data", "gcm", "yr")
  missing <- missing[!missing %in% names(dots)]
  if (length(missing)) {
    stop(paste("Please provide", paste(missing, collapse = ", ")))
  }
  
  # Convert Data into a data.frame
  if(!is(dots$Data, 'data.frame')){
    dots$Data <- as.data.frame(dots$Data)
  } 
  
  # Empty list
  # output <- vector('list' , length = nrow(Data))
  
  #generate loop for the data
  args <- append(list(X = 1:nrow(dots$Data), FUN = DD), dots)
  
  if (parallel) {
    plan(multisession, workers = nworkers)
    on.exit(future:::ClusterRegistry("stop"), add = TRUE)
  
    args$future.seed <- TRUE
    output <- do.call(future_lapply, args)
  } else {
    output <- do.call(lapply, args)
  }
  
  # Assign names only to non-NULL entries
  names(output) <- paste(dots$Data$pixelID)
  
  return(output)
}


DD <- function(i, Data, start.column = 3, 
               end.column = ncol(Data), 
               IS_NA = FALSE,
               Resolution = 4000, 
               Sample.N = 1000, 
               gcm,
               yr) {
  # Transpose 
  df <- as.data.frame(t(Data[i, start.column:end.column]))
  colnames(df) <- 'log.mass'
  
  if(IS_NA) {
    # Remove NAs
    df <- na.omit(df)
    df$species <- rownames(df)
    df <- df[order(df$log.mass),, drop=FALSE]
  } else {
    # temp <- sort(df[-which(df == 0), ], decreasing = FALSE)
    # df <- as.matrix(data.frame(scientific.name = names(temp),
    #                            log.mass = temp))
    # 
    stop()
  }
  
  rdata <- data.frame(log.mass = as.numeric(df$log.mass))
  
  ##Add a check to skip rows with insufficient data 
  if (nrow(rdata) < 2) {
    cat('Skipping row', i, 'due to insufficient data.\n')
    output <- NA
    return(output)
  }
  ##Step 2: Generate the neutral null----
  #define the parameters
  Dmax <- max(rdata, na.rm = TRUE)
  Dmin <- min(rdata, na.rm = TRUE)
  ds <- (Dmax - Dmin) / Resolution
  MaxK <- (Dmax - Dmin) / 2
  MinK <- ds * 2
  
  #define h's to analyze
  ks <- seq(MinK, MaxK, by = 1 / Resolution)
  
  #generate matrix
  bws <- matrix(data = NA, nrow = length(ks), ncol = 1)
  
  for(j in seq_along(ks)) {
    # Add a check to skip non-positive bandwidth values
    if (ks[j] <= 0) {
      cat("Warning: Skipping non-positive bandwidth value.\n")
      next
    }
    
    #calculate KS density estimate
    KSdens <- density(rdata$log.mass, bw = ks[j], "gaussian", adjust = 1)
    
    #Test if the ksdensity is unimodal
    TF <- which(diff(sign(diff(KSdens$y))) == 2) + 1
    if (length(TF) == 0) { 
      bws[j] = 1 
    } else { 
      bws[j] = 0
    }
  }
  r = min(which(bws == 1))
  hnull = ks[r]
  
  # Check if hnull is a finite and positive value before using it
  if (is.finite(hnull) && hnull > 0) {
    NNull <- density(rdata$log.mass, bw = hnull, "gaussian", adjust = 1)
    
    N <- nrow(rdata)
    null.samples <- matrix(data = 0, ncol = Sample.N, nrow = N)
    
    for (k in 1:Sample.N) {
      rand.N <- sample(NNull$x, N, replace = TRUE, prob = NNull$y)
      null.samples[, k] <- sort(rand.N, decreasing = FALSE)
    }
    
    gaps.log10.data <- diff(rdata$log.mass)
    gaps.null.samples <- diff(null.samples, decreasing = FALSE)
    gap.percentile <- matrix(data = 0, nrow = length(gaps.log10.data), ncol = 1)
    
    for (l in seq_along(gaps.log10.data)) {
      gap.percentile[l] <- ecdf(gaps.null.samples[l, ])(gaps.log10.data[l])
    }
    #browser()
    pixelID <- Data$pixelID[i]
    lat <- Data$lat[i]
    long <- Data$long[i]
    
    Bootstrap.gaps <- rbind(gap.percentile, 0)
    output <- data.frame(species = df$species,
                              log.mass = rdata$log.mass,
                              d.value = Bootstrap.gaps,
                              pixelID = rep(pixelID, each = length(df$species)),
                              gcm = rep(gcm, each = length(df$species)),
                              yr = rep(yr, each = length(df$species)),
                              lat = rep(lat, each = length(df$species)),
                              long = rep(long, each = length(df$species)))
    
    row.names(output) <- NULL
    
    # Progress report...
    cat('Completed row', i, '\n')
  } else {
    # Skip this row if hnull is not finite or not positive
    output <- NA
    cat('Skipping row', i, 'due to non-finite or non-positive hnull.\n')
  }
  return(output)
}

# debug(Discontinuity_Barichievy)
# undebug(Discontinuity_Barichievy)


#=============================================================================
## Caleb's functions for processing Discontinuity_Barichievy output
#=============================================================================

##### GRI.breaks()
# Helper function used internally by the "Aggregations" function. 
GRI.breaks <- function(Data, 
                       columnName = 'd.value',
                       Tab = pwrtab,
                       SetThreshold = FALSE,
                       constant.threshold = 0.8) {
  
  # Convert to data.frame
  Data1 <- as.data.frame(Data)
  Tab1 <- as.data.frame(Tab)
  
  #check column names 
  colnames(Data1)

  if (SetThreshold) {
    threshold <- constant.threshold
  } else {
    # Find species richness of data
    rich <- signif(nrow(Data1), 1)
    
    # Get the threshold
    threshold <- Tab1[which(Tab1$richness == rich), "threshold"]
  }

  # Check if the column "d.value" exists in Data1
  if (columnName %in% colnames(Data1)) {
    # Continue with your code
    breaks <- data.frame(rbind(Data1[which(Data1[, "d.value"] >= threshold), ],
                               Data1[nrow(Data1), ]))
  } else {
    stop('Column ', columnName, ' not found in Data1.')
  }
  # browser()
  return(breaks)
}
# debug(GRI.breaks)
# undebug(GRI.breaks)

##### Aggregations()
# Function using the "constant power table" to identify significant breaks in 
# discontinuity analysis output from the "Discontinuity_Barichievy" function.
# This function is used internally by the "DA_loop" function in the 
# "SpatialRegimes_Tracking" script.
# 
# Arguments:
#   GRI.output = output from the "Discontinuity_Barichievy" function. 
#   Tab        = constant power table. See pwrtab file.
# 
# Returns:
#   A data.frame similar to GRI.output with the addition of a numeric "agg" 
#   column indicating the body mass aggregation membership for each species 
#   observed in a location (e.g., BBS route) in a given year.
Aggregations <- function(GRI.output, 
                         Tab = pwrtab,
                         SetThresh = FALSE,
                         constantThresh = 0.8,
                         yr, gcm) {
  
  #Find threshold breaks only for data frames with the 'd.value' column
  df <- lapply(GRI.output, function(x) {## for each pixel...
    if (is.data.frame(x) && !is.null(x) && nrow(x) > 0 && "d.value" %in% colnames(x)) {
      return(GRI.breaks(x, SetThreshold = SetThresh, 
                        constant.threshold = constantThresh))
    } else {
      return(NA)  # Return NA for invalid or empty data frames
    }
  })
  #Remove NULL entries from the list
  find.thresholds <- df[sapply(df, function(x) !is.null(x))]

  # For each threshold break attribute an aggregation ID
  ## this is unnecessary -- aggs column is never used.
  # find.thresholds <-  lapply(df, function(x) {
  #   if(nrow(x) > 0){
  #     return(cbind(x, aggs = 1:nrow(x)))
  # } else {
  #     return(x)
  #   }
  # })
    

  # Create an aggregation column
  GRI.output <- lapply(GRI.output, function(x) {
   # browser()
    if (is.list(x) && !is.null(x) && is.data.frame(x) && nrow(x)>0) {
      #print("Entering if condition")
      return(cbind(x, aggs = rep(0, nrow(x))))
     
    } else {
    #  print("Entering else condition")
      return(x)  # Return x as is if nrow(x) <= 0 or x is NA
    }
  })
    
  # Assign species to aggregations
  for (i in 1:length(find.thresholds)) { ## for each pixel
    temp.thresh <- find.thresholds[[i]]
    # Check if GRI.output[[i]] is a data frame and has the "log.mass" column
    if (is.data.frame(GRI.output[[i]]) && "log.mass" %in% colnames(GRI.output[[i]])) {
     
      for (j in 1:nrow(temp.thresh)) { ## for each aggregation threshold
        ## spp w/ bMass > threshold are placed in the NEXT aggregation (first aggregation is 0 here)
        ## so the discontinuity break is on the LAST species of each aggregation
       # browser()
        indices <- which(GRI.output[[i]]$log.mass > temp.thresh[j, "log.mass"])
           # Check if indices are not empty
        if (length(indices) > 0) {
          GRI.output[[i]][indices, "aggs"] <- j
        }
      }
      
      ## now make aggregation IDs from 1 to N.
      GRI.output[[i]]$aggs <- GRI.output[[i]]$aggs + 1
      
    } else {
      warning("Invalid structure or missing 'log.mass' column in GRI.output[[", i, "]]. Skipping.")
    }
  }

  # out_dir <- './outputs/Aggregations/discont_breaks'
  # ifelse(!dir.exists(out_dir), dir.create(out_dir, recursive = TRUE),
  #        print('Directory already exists'))
  
  # Iterate over find.thresholds
  # for (i in seq_along(find.thresholds)) {
  #   # Generate a unique filename for each iteration
  #   filename <- glue::glue('{out_dir}/discont_breaks_{i}_{gcm}_{yr}.qs')
  #   
  #   # Save the data frame at the current index in find.thresholds
  #   qs::qsave(find.thresholds[[i]], file = filename)
  #}
   return(GRI.output)
}
# debug(Aggregations)

##### Aggregations_ggplot()
# A helper function for creating unique idenitifiers for body mass 
# aggregations for each factor combination (e.g., year, route, aggregation).
#
# Arguments:
#   Aggs          = output from "Aggregations" function 
#   paste.columns = character vector indicating which columns of the 
#                   data.frame to paste together to generate the unique 
#                   identifier.
# 
# Returns:
#   A data frame ready for use in the "SpReg_DA_Plot" function in the 
#   "functions_MakePlots" script.
Aggregations_ggplot <- function(Aggs, paste.columns) {
  
  Aggs$ggplot_aggs <- do.call(paste, Aggs[ , paste.columns])
  
  return(Aggs)
}
# debug(Aggregations_ggplot)
# undebug(Aggregations_ggplo