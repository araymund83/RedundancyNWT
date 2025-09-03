# Load libraries ----------------------------------------------------------
library(pacman)
p_load(dplyr,furrr, ggnewscale, ggplot2, ghibli, glue, qs, scales, tidyverse)

# Clean environment -------------------------------------------------------
g <- gc(reset = TRUE)
rm(list = ls())

## Load data
yrs <- c('2011', '2031', '2091')
yr <- yrs[]
gcms <- c('CanESM2', 'CCSM4', 'INM-CM4')
gcm <- gcms[2]
out_dir <- './tables/sigAggs'

traitsFG <- qs::qread('./inputs/traitsMassFG_2025.qs')

traitsFG <- traitsFG %>% select(FG) %>% 
  rownames_to_column('species')

# Define FG categories
fg_levels <- c(
  "INS-Sal", "OMN-For", "INS-Gle", "INS-Exc", "INS-For",
  "OMN-Sca", "OMN-Gle", "GRA-Gle", "INS-Scr", "OMN-Exc"
)



agg_colors <- c('#A65041', '#E7CDC2', '#80A0C7', '#394165', '#A890A8','#B1934A', 
                '#DCA258', '#100F14', '#8B9DAF', '#EEDA9D', '#90A878')

# Define the color palette (10 colors for 10 FGs)
fg_palette <- c(
  "#1B9E77",  # INS-Sal
  "#D95F02",  # OMN-For
  "#7570B3",  # INS-Gle
  "#E7298A",  # INS-Exc
  "#66A61E",  # INS-For
  "#E6AB02",  # OMN-Sca
  "#A6761D",  # OMN-Gle
  "#A6CEE3",  # GRA-Gle
  "#FB9A99",  # INS-Scr
  "#B2DF8A"   # OMN-Exc
)

# Create named palette
fg_palette_named <- setNames(fg_palette, fg_levels)

file<- list.files(out_dir, pattern =  'aggs', full.names = TRUE)

make_plot_discontinuities <- function(data, gcm, yr) {
  # Filter the file for the correct GCM and Year
  file <- grep(gcm, file, value = TRUE)
  file <- grep(yr, file, value = TRUE)
  data <- lapply(file, qs::qread)
  data <- flatten(data)
  
  aggPlot <- lapply(seq_along(data), function(i) {
    df <- data[[i]]
    
    # Check if the column names in df and traitsFG match the required names
    if("species" %in% colnames(df)) {
      # If 'species' exists in df, rename the column in traitsFG to 'species'
      if("sp" %in% colnames(traitsFG)) {
        traitsFG <- traitsFG %>% rename(species = sp)
      }
    } else if("sp" %in% colnames(df)) {
      # If 'sp' exists in df, rename the column in traitsFG to 'sp'
      if("species" %in% colnames(traitsFG)) {
        traitsFG <- traitsFG %>% rename(sp = species)
      }
    }
    
    
    # Join with FG info
    df <- left_join(df, traitsFG, by = "species")
    
    if (any(is.na(df$FG))) {
      warning("Some species did not match functional groups.")
    }
    
    df$spRank <- order(df$log.mass)
    
    ## Get aggregation thresholds/breaks (last species on each aggregation)
    aggThresholds <- sapply(split(df, df$aggs), function(ddf) {
      ddf[nrow(ddf), "log.mass"]
    }, simplify = TRUE)
    
    # Create a data frame with aggThresholds and aggs to create boundaries for each aggregation
    thresholds_df <- data.frame(aggThresholds = aggThresholds, aggs = names(aggThresholds))
    
    # Create the rectangles: xmin, xmax, ymin, ymax for each aggregation
    rectangles_df <- data.frame(
      xmin = c(-Inf, aggThresholds[-length(aggThresholds)]),  # Left boundary of each agg
      xmax = aggThresholds,  # Right boundary of each agg
      ymin = min(df$spRank),  # Bottom of the plot (fixed)
      ymax = max(df$spRank),  # Top of the plot (fixed)
      fill_color = agg_colors[1:length(aggThresholds)]  # Color for each aggregation (from agg_colors vector)
    )
    
    # Create the base plot with points
    p <- ggplot() +
      # Plot points with color by Functional Group using df
      geom_point(data = df, aes(x = log.mass, y = spRank, color = FG), size = 5) +  # Points by FG color
      
      ggtitle(glue('Discontinuities breaks for {gcm}-{yr}'),
              subtitle = glue('{names(data)[i]}')) +
      theme_bw() +
      theme(legend.position = 'bottom',
            plot.title = element_text(size = 16, face = 'bold'),
            plot.subtitle = element_text(size = 14),
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 13, face = 'bold')) +
      labs(x = 'log10 body mass', y = 'Ranked body mass', 
           color = 'Functional Groups') +  # Set legend title for points to "Functional Groups"
      
      # Add colored rectangles for each aggregation using geom_rect()
      geom_rect(data = rectangles_df, 
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_color),
                alpha = 0.2, color = NA) +  # Rectangles for each aggregation with no border color
      
      scale_fill_identity() +  # Use the colors from agg_colors directly for fill
      scale_color_manual(values = fg_palette_named)  # Apply FG colors for points (FG)
    
    # Now add the dashed lines with a new color scale for aggregation boundaries
    p2 <- p + 
      ggnewscale::new_scale_color() +  # New color scale for the dashed lines
      geom_vline(data = thresholds_df, aes(xintercept = aggThresholds, color = factor(aggs)),
                 linetype = "dashed", size = 1) +  # Dashed lines for boundaries
      scale_color_manual(values = agg_colors) +  # Apply aggregation boundary colors (lines)
      labs(color = "Aggregations")  # Set legend title for dashed lines to "Aggregations"
    
    # Save the plot
    out <- glue('./figs/aggregations2025_2/tbiSig/{gcm}/{yr}')
    if (!dir.exists(out)) dir.create(out, recursive = TRUE)
    
    filename <- glue('{out}/aggs_{names(data)[i]}_{gcm}_{yr}.png')
    ggsave(filename, plot = p2, width = 12, height = 10, units = "in", dpi = 300)
    
    # Return the plot explicitly
    return(p2)
  })
  
  # Return all plots
  return(aggPlot)
}

lapply(aggPlot, print)
plots <- make_plot_discontinuities(data = file, gcm, yr)
