rm(list=ls())
setwd(dir="/Users/daseinl/Coding/Predator-Prey-Interaction-Networks")
getwd()

if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("permute", quietly = TRUE)) install.packages("permute")
if (!requireNamespace("lattice", quietly = TRUE)) install.packages("lattice")
if (!requireNamespace("vegan", quietly = TRUE)) install.packages("vegan")
if (!requireNamespace("statnet.common", quietly = TRUE)) install.packages("statnet.common")
if (!requireNamespace("sna", quietly = TRUE)) install.packages("sna")
if (!requireNamespace("bipartite", quietly = TRUE)) install.packages("bipartite")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("ggpubr", quietly = TRUE)) install.packages("ggpubr")

library(dplyr)
library(tidyr)
library(tidyverse) # version 2.0.0
library(readxl)
library(permute) 
library(lattice)
library(vegan)
library(statnet.common)
library(sna)
library(bipartite) # version 2.18
library(ggplot2) # version 3.4.2l
library(reshape2)
library(ggpubr)


### Start define function
# The code to filter arthropods and ambiguous species by Ariane
filter_by_species <- function(pre_data){
  arthropoda_and_species_only <- data %>%
    dplyr::filter(phylum == "Arthropoda") %>% # only keep arthropods
    dplyr::filter(order != "Mesostigmata"&  # exclude arachnids which are mites and therefore not dietary items
                    order != "Sarcoptiformes"&
                    order != "Trombidiformes"&
                    species != "Pseudacanthotermes militaris") %>%
    filter(!is.na(species)) %>%  # only keep OTUs identified to genus or species level
    filter(!is.na(genus)) %>%
    distinct()
  
  # Get column named "OTU" and start only with G
  columns_to_keep <- grep("^G[0-9]+|^OTU$", names(arthropoda_and_species_only), value = TRUE)
  pro_data <- arthropoda_and_species_only[, columns_to_keep]
  return(pro_data)
}

# The RRA code to eliminate bias from small data by Ariane
filter_by_RAA <- function(data_pre_RAA){
  num_col <- ncol(data_pre_RAA)
  data_pro_RAA <- data_pre_RAA %>%
    dplyr::select(1:num_col) %>%
    column_to_rownames(var = "OTU") %>%
    mutate(across(1:num_col - 1, ~ ifelse(. < 0.003 * sum(.), 0, .))) %>% # apply RRA filtering threshold (here 0.3% of total reads)
    filter(rowSums(select(., 1:num_col - 1)) != 0) %>%  # filter out OTUs that do not occur following the RRA filtering
    rownames_to_column(var = "OTU")
  return(data_pro_RAA)
}

filter_by_site <- function(pre_data, meta_data, site_name){
  meta_data <- filter(meta_data, Site.Type == site_name)
  columns_to_select <- intersect(names(pre_data), meta_data$Sample_number)
  pro_data <- pre_data[, c("OTU", columns_to_select)]
  return(pro_data)
}

# Combine column using the meta data to species
combine_data <- function(bird_meta, data_RAA, data) {
  # Combine column using the meta data to species
  bird_species_mapping <- setNames(bird_meta$BirdLife.Name.CORRECT, bird_meta$Specimen_number)
  data_merged <- data.frame(OTU = data_RAA$OTU)
  unique_species <- unique(bird_species_mapping)
  
  for (species in unique_species) {
    specimen_numbers <- names(bird_species_mapping[bird_species_mapping == species]) # bird code
    columns_to_merge <- names(data_RAA) %in% specimen_numbers
    if (sum(columns_to_merge) > 0) { # get sum value of each species
      average_data <- rowSums(data_RAA[, columns_to_merge, drop = FALSE], na.rm = TRUE)
      data_merged[[species]] <- average_data
    }
  }
  
  # Map insect species
  insect_species_mapping <- setNames(data$species, data$OTU)
  data_merged$OTU <- insect_species_mapping[as.character(data_merged$OTU)]
  colnames(data_merged)[colnames(data_merged) == "OTU"] <- "Insect_species"
  
  # Summarize data by insect species
  data_merged <- data_merged %>% 
    group_by(Insect_species) %>%
    summarise(across(everything(), ~ sum(.x, na.rm = TRUE)))
  
  return(data_merged)
}

# produce a table containing bird forest coverage rate
merge_bird_forest_coverage <- function(meta){
  forest_coverage <- meta %>% 
    group_by(BirdLife.Name.CORRECT) %>% 
    summarise(average.Distance.to.Forest.km = mean(Distance.to.Forest.km))
  return(forest_coverage)
}

# produce a table containing insect role
merge_insect_role <- function(meta){
  role <- meta[, c("Species", "Role")] %>%
    filter(if_all(everything(), ~ !is.na(.) & . != "Unknown" & . != ""))
  return(role)
}


### Bipartite network visualization
# function for coloring the forest bird green
colour_sel <- function(plant_spp, first_col){
  col_names <- colnames(plant_spp)
  colours <- if_else(col_names %in% first_col,"#238E23", "darkgoldenrod2", 'grey')
  return(colours)
}

save_forest_category_bipartite <- function(forest_coverage, web, name){
  first_col <- forest_coverage %>% filter(average.Distance.to.Forest.km < 0)
  first_col<- unique(first_col$BirdLife.Name.CORRECT)
  
  time_stamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  jpeg(paste0("bipartite_plot_", name, "_", time_stamp, ".jpg"), width = 2000, height = 1000, units = 'px', res = 300) 
  plotweb(web, col.high = colour_sel(plant_spp = web, first_col= first_col),
          text.rot=90, 
          col.low= 'grey',
          bor.col.interaction = NA,
          #low.lab.dis = 0.01,
          #high.lab.dis = 0.3,
          #low.spacing = 0.02,
          #high.spacing = 0.8,
          labsize = 0.5,
          low.lablength = 0,
          # high.lablength = 0,
          col.interaction = colour_sel(plant_spp = web, first_col= first_col),
          y.lim=c(-0.6,2.5),
          method = 'normal',
          #add = TRUE
  )
  dev.off()
}

# produce the none-zero number of web
count_nonzero <- function(X, A, B) {
  count_A <- 0; count_B <- 0; count_shared <- 0
  for (i in 1:nrow(X)) {
    row_data <- X[i, ]
    has_nonzero_A <- any(row_data[A] != 0)
    has_nonzero_B <- any(row_data[B] != 0)
    if (has_nonzero_A && !has_nonzero_B) {count_A <- count_A + 1}
    else if (!has_nonzero_A && has_nonzero_B) {count_B <- count_B + 1}
    else if (has_nonzero_A && has_nonzero_B) {count_shared <- count_shared + 1}
  }
  return(list(count_A, count_shared, count_B))
}

# Visualize the bar chart for different insect
save_bar_chart_insect_role <- function(df_plot, name){
  time_stamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  pp <- ggbarplot(df_plot, x="Insect Role", y="proportion", color="black", fill="Bird Type") +
    scale_fill_manual(values = c("#238E23", "grey40", "darkgoldenrod2")) +
    coord_flip() +
    theme_minimal() +
    theme(axis.title.y = element_blank(),
          legend.title = element_blank(),
          legend.position = "bottom")
  ggsave(paste0("proportion_", name, "_", time_stamp, ".jpg"), pp)
}



### Start data processing
# load data
data <- read.csv("./009MRes Denian Li Project/Data/Ghana_fwh2_consensus_RAW.csv")
bird_meta <- read.csv("./009MRes Denian Li Project/Data/MistNetting_wrangled_fullcoor.csv")
insect_meta <- read.csv("./009MRes Denian Li Project/Data/Pest_annotation_Ghana_060223_cleaned.csv")
bird_meta_update <- read_excel("./009MRes Denian Li Project/Data/sample_BF_dist_cover_allcounts.xlsx", sheet="sample_BF_dist_cover_allcounts")

# Filter the data to keep only arthropods exclude arachnids which are mites and therefore not dietary items
data_filtered <- filter_by_species(data)
data_RAA <- filter_by_RAA(data_filtered) # RAA filter to eliminate bias from small data by threshold of 0.3%

# site_type <- c("Forest Interior", "Forest Edge", "Agriculture")
data_forest <- filter_by_site(data_RAA, bird_meta_update, "Forest Interior")
data_edge <- filter_by_site(data_RAA, bird_meta_update, "Forest Edge")
data_agriculture <- filter_by_site(data_RAA, bird_meta_update, "Agriculture")

# Merge bird forest coverage metadata and insect role metadata
bird_forest_coverage <- merge_bird_forest_coverage(bird_meta)
bird_forest_coverage <- bird_forest_coverage[order(bird_forest_coverage$average.Distance.to.Forest.km), ]
insect_role <- merge_insect_role(insect_meta)

# Combine metadata with the filtered data
data_in_site_type <- list(
  "Total" = combine_data(bird_meta, data_RAA, data), 
  "Forest Interior" = combine_data(bird_meta, data_forest, data), 
  "Forest Edge" = combine_data(bird_meta, data_edge, data), 
  "Agriculture" = combine_data(bird_meta, data_agriculture, data)
  )

# Main loop for data in different site type
for(name in names(data_in_site_type)){
  # Visualize bipartite plot
  data_merged <- data_in_site_type[[name]]
  valid_columns <- bird_forest_coverage$BirdLife.Name.CORRECT %in% colnames(data_merged)
  data_merged_reordered <- data_merged[, bird_forest_coverage$BirdLife.Name.CORRECT[valid_columns]]
  data_merged_reordered <- cbind(Insect_species=data_merged[["Insect_species"]], data_merged_reordered)
  bird_insect_web <- data_merged_reordered %>% column_to_rownames(var = 'Insect_species')
  
  save_forest_category_bipartite(bird_forest_coverage, bird_insect_web, name)
  
  # Visualize proportion bar chart
  forest_birds <- bird_forest_coverage %>%
    filter(average.Distance.to.Forest.km <= 0) %>%
    pull(BirdLife.Name.CORRECT) %>%
    intersect(colnames(bird_insect_web))
  
  non_forest_birds <- bird_forest_coverage %>%
    filter(average.Distance.to.Forest.km > 0) %>%
    pull(BirdLife.Name.CORRECT) %>%
    intersect(colnames(bird_insect_web))
  
  role_groups <- split(insect_role, insect_role$Role)
  bird_insect_web_by_role <- list()
  for(role in names(role_groups)) {
    species_in_role <- role_groups[[role]]$Species
    bird_insect_web_by_role[[role]] <- bird_insect_web %>% 
      filter(row.names(.) %in% species_in_role)
  }
  
  static_data <- cbind(
    c("Pest", "Non-pest", "Natural enemy", "Vector of animal disease"),
    rbind(count_nonzero(bird_insect_web_by_role$Pest, forest_birds, non_forest_birds),
          count_nonzero(bird_insect_web_by_role$`Non-pest`, forest_birds, non_forest_birds),
          count_nonzero(bird_insect_web_by_role$`Natural enemy`, forest_birds, non_forest_birds),
          count_nonzero(bird_insect_web_by_role$`Vector of animal disease`, forest_birds, non_forest_birds)))
  
  colnames(static_data) <- c("Insect Role", "Forest birds only", "Shared", "Non-forest birds only")
  static_df <- as.data.frame(static_data)
  static_df[, -1] <- lapply(static_df[, -1], as.numeric)
  df_percentages <- static_df %>%
    mutate(across(-`Insect Role`, ~ . / rowSums(static_df[, -1]))) %>%
    mutate(`Insect Role` = factor(`Insect Role`, levels = rev(`Insect Role`)))
  df_plot <- melt(df_percentages, id.vars = "Insect Role", variable.name = "Bird Type",
                  value.name = "proportion")
  
  save_bar_chart_insect_role(df_plot, name)
}







