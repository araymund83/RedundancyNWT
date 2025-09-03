# Load libraries ----------------------------------------------------------
library(pacman)
p_load(data.table, dplyr, ggplot2, ggrepel, glue, qs, tidyverse, vegan)

path <- './outputs/nmds_results_2025'
files <- list.files(path, pattern = '10percent')
gcms <- c('CanESM2', 'CCSM4', 'INM-CM4')
gcm <- gcms[1]
files <- grep(gcm, files,value = TRUE)


make_nmds_plot <- function(files, path){
  make_plot <- function(file, path){

    #extract region name from file name
    region <- str_extract(file, "(?<=nmds_)[^_]+_[^_]+")

    # Apply conditional adjustments based on specific cases
    if (region == "Low_Arctic") {
      region <- "Low_Arctic_north"
    } else if (grepl("^Mid-Boreal", region)) {
      region <- "Mid-Boreal"
    }
  
    #read sites file
    sites <- qs::qread(glue('{path}/sites_{region}_{gcm}_10_percent_allyrs.qs'))
    
    #extract resulsts from NMDS
    nmds_out <- qread(glue('{path}/{file}'))
    NMDSoutpoints <- data.table(nmds_out$points)
    NMDSoutpoints$yr <- as.factor(sites$yr)
    NMDSoutpoints$pixelID <- sites$pixelID
    
    # Extract species scores
    species_scores <- as.data.frame(scores(nmds_out, "species"))
    species_scores$species <- rownames(species_scores)
    
    #Extract stress value 
    stress <-round(nmds_out$stress, 3)
    
    
    custom_theme <- theme_bw() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
 
    # Create plot
    p <- ggplot() + 
      geom_point(data = NMDSoutpoints, aes(x = MDS1, y = MDS2, #shape = yr, 
                                           colour = yr), size = 3, alpha = 0.4) +
      stat_density2d(data = NMDSoutpoints, aes(x = MDS1, y = MDS2, fill = yr), 
                     geom = 'polygon', alpha = 0.3, h = 0.5) +
      scale_colour_manual(name = NULL, values = c("2011" = "#FFB00D", "2031" = "#0072B2", "2091" = "#AA3377")) +
      scale_fill_manual(values =c("2011" = "#FFB00D", "2031" = "#0072B2", "2091" = "#AA3377")) +
      #geom_text(data = species_scores, aes(x = NMDS1, y = NMDS2, label = species), 
       #         alpha = 0.9) +  
      coord_fixed(ratio = 1) +
      theme(axis.title = element_text(size = 10, face = "bold", colour = "#4D4D4D"), 
            panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "#4D4D4D"), 
            axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
            legend.text = element_text(size = 9, colour = "#4D4D4D"),
            plot.subtitle = element_text(face = 'bold')) +
      guides(shape = guide_legend(title = NULL), colour = guide_legend(title = NULL), 
             fill = 'none') + 
      custom_theme +
      geom_text(aes(x = -1.7, y = 1.6, label = glue('Stress:{stress}')), 
                hjust = 0, vjust = 1, size = 4, color = 'grey20', fontface = 'bold',
                nudge_y = 0.05, nudge_x = 0.05) +
      ggtitle(glue("NMDS for {region}"),
              subtitle = glue('{gcm}')) + 
      geom_segment(data = species_scores, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2 - 0.1),
                   arrow = arrow(length = unit(0.1, "inches")), color = "grey25") +
      geom_text_repel(data = species_scores, aes(x = NMDS1, y = NMDS2, label = species), 
                       box.padding = 0.5, point.padding = 0.5)

   
     out <-  glue('./figs/nmds_sampling/10percent/5Regions/2025')
    ifelse(!dir.exists(out), dir.create(out, recursive = TRUE),
           print('Folder already exists'))
    # Save plot to PNG
    ggsave(plot = p, filename = glue('{out}/nmds_{region}_{gcm}_10percent_allyrs_2.png'), 
           units = 'in', width = 5, height = 5, dpi = 300) 
  }
  
  # Apply the create_plot function to each file in the list
  lapply(files, make_plot, path = path)
}

make_nmds_plot(files, path)


# permanova_tables --------------------------------------------------------

library(qs)
library(dplyr)
library(tidyr)
library(stringr)
library(glue)

# Set the directory containing the PERMANOVA results
path <- "outputs/permanova_2025"
files <- list.files(path, pattern = "_10_percent_allyrs\\.qs$", full.names = TRUE)
# Filter only CanESM2 files
fls <- files[grepl("CanESM2", files)]


# Function to extract statistics from each file
extract_permanova <- function(file) {
  result <- qread(file)
  result_df <- as.data.frame(result)
  
  # Remove extension and suffix
  filename <- basename(file) %>%
    str_remove("\\.qs$") %>%
    str_remove("_10_percent_allyrs$")
  
  # Split parts by "_"
  parts <- str_split(filename, "_")[[1]]
  
  # GCM is always the last element
  gcm <- parts[length(parts)]
  
  # Region is everything before the GCM
  region_name <- paste(parts[1:(length(parts) - 1)], collapse = "_")
  
  cat("Processing region:", region_name, "GCM:", gcm, "\n")
  
  result_clean <- result_df %>%
    mutate(
      Term = rownames(result_df),
      GCM = gcm,
      Region = region_name,
      df = Df,
      SS = round(`Sum Sq`, 2),
      MS = round(`Mean Sq`, 4),
      F = round(`F value`, 2),
      `p-value` = ifelse(is.na(`Pr(>F)`), NA,
                         ifelse(`Pr(>F)` < 0.01, "<0.01", round(`Pr(>F)`, 3)))
    ) %>%
    select(GCM, Region, Term, df, SS, MS, F, `p-value`)
  
  return(result_clean)
}

# Process all files
permanova_table <- bind_rows(lapply(fls, extract_permanova))
# View and save
print(permanova_table)
write.csv(permanova_table, "Table4_PERMANOVA_summary.csv", row.names = FALSE)
