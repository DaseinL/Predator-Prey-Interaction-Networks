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
  bird_species_mapping <- setNames(bird_meta$BirdLife.Name.CORRECT, bird_meta$Specimen.Number)
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

# Visualize bird community capture patterns across site types
statistic_richness <- function(matrix_list, bird_coverage) {
  results <- list()
  
  # classify birds
  classified_birds <- bird_coverage %>%
    mutate(Bird_Type = ifelse(average.Distance.to.Forest.km < 0, "Forest bird", "Non-Forest bird"))
  
  for (i in seq_along(matrix_list)) {
    matrix <- matrix_list[[i]]
    
    # count non-empty columns
    non_empty_cols <- colSums(matrix != 0) > 0
    non_empty_col_names <- colnames(matrix)[non_empty_cols]
    
    # filter non-empty birds
    non_empty_birds <- classified_birds %>% filter(BirdLife.Name.CORRECT %in% non_empty_col_names)
    
    # count forest and non-forest birds
    forest_birds <- non_empty_birds %>% filter(Bird_Type == "Forest bird") %>% nrow()
    non_forest_birds <- non_empty_birds %>% filter(Bird_Type == "Non-Forest bird") %>% nrow()
    
    results <- c(results, forest_birds, non_forest_birds)
  }
  
  return(results)
}

# produce a table containing bird forest coverage rate
merge_bird_forest_coverage <- function(meta){
  forest_coverage <- meta %>% 
    group_by(BirdLife.Name.CORRECT) %>% 
    summarise(average.Distance.to.Forest.km = mean(Bird.Distance.to.Forest.km))
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
colour_bird <- function(plant_spp, first_col){
  col_names <- colnames(plant_spp)
  colours <- if_else(col_names %in% first_col,"#238E23", "darkgoldenrod2")
  return(colours)
}

# function for coloring the pest red
colour_pest <- function(plant_spp, first_col){
  row_names <- rownames(plant_spp)
  colours <- if_else(row_names %in% first_col,"red", "grey")
  return(colours)
}

save_forest_category_bipartite <- function(forest_coverage, insect_rol, web, name){
  first_col <- forest_coverage %>% filter(average.Distance.to.Forest.km < 0)
  first_col <- unique(first_col$BirdLife.Name.CORRECT)
  pest_col <- insect_rol %>% filter(Role == "Pest")
  pest_col <- unique(pest_col$Species)
  
  time_stamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  jpeg(paste0("bipartite_plot_", name, "_", time_stamp, ".jpg"), width = 2000, height = 1000, units = 'px', res = 300) 
  plotweb(web, col.high = colour_bird(plant_spp = web, first_col= first_col),
          col.low = colour_pest(plant_spp = web, first_col= pest_col),
          text.rot=90, 
          bor.col.interaction = NA,
          #low.lab.dis = 0.01,
          #high.lab.dis = 0.3,
          #low.spacing = 0.02,
          #high.spacing = 0.8,
          labsize = 0.5,
          low.lablength = 0,
          # high.lablength = 0,
          col.interaction = colour_bird(plant_spp = web, first_col= first_col),
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
  ggsave(paste0("proportion_", name, "_", time_stamp, ".jpg"), pp,width = 6.67, height = 2.5, units = "in")
}

# the function to convert data to edges and nodes which can be deal with Geghi
convert_data_for_gephi <- function(interaction_data, bird_info_data, output_nodes_file, output_edges_file) {
  interaction_data <- as.data.frame(interaction_data)
  rownames(interaction_data) <- interaction_data[, 1]
  interaction_data <- interaction_data[, -1]
  edges <- list()
  
  for (insect in rownames(interaction_data)) {
    for (bird in colnames(interaction_data)) {
      interaction <- interaction_data[insect, bird]
      if (interaction != 0) {
        edges <- append(edges, list(data.frame(Source = insect, Target = bird, Weight = interaction)))
      }
    }
  }
  
  edges_df <- bind_rows(edges)
  insects_list <- rownames(interaction_data)
  birds_list <- colnames(interaction_data)
  
  # Add meta data to nodes
  bird_info_data <- bird_info_data %>%
    select(BirdLife.Name.CORRECT, average.Distance.to.Forest.km) %>%
    rename(Id = BirdLife.Name.CORRECT, averageDistanceToForest = average.Distance.to.Forest.km)
  birds_info <- data.frame(Id = birds_list, Type = 'Bird', stringsAsFactors = FALSE)
  birds_info <- left_join(birds_info, bird_info_data, by = "Id")
  birds_info <- birds_info %>%
    mutate(IsForestBird = ifelse(is.na(averageDistanceToForest) | averageDistanceToForest > 0, 'Non-Forest', 'Forest'))
  
  insects_info <- data.frame(Id = insects_list, Type = 'Insect', stringsAsFactors = FALSE)
  
  nodes_info <- bind_rows(insects_info, birds_info)
  nodes_df <- as.data.frame(nodes_info)
  
  write_csv(nodes_df, output_nodes_file)
  write_csv(edges_df, output_edges_file)
  print(paste("Nodes and edges have been saved to", output_nodes_file, "and", output_edges_file))
}


first_col_to_rowname <- function(net){
  net <- as.data.frame(net)
  rownames(net) <- net[, 1]
  net <- net[, -1]
  return(net)
}

# function to draw an extinction curve
extinction_curve <- function(matrix_network, insect_roles, extinction_sequence = NULL) {
  if (is.null(extinction_sequence)) {
    extinction_sequence <- sample(colnames(matrix_network))
    repeat_count <- 1000
  }else{
    repeat_count <- 1
  }
  remaining_pests_list <- rep(0, length(extinction_sequence) + 1)
  
  # 计算初始剩余害虫数量并保存在列表中
  for (t in 1:repeat_count){
    i <- 1
    matrix_network_copy <- matrix_network
    matrix_network_copy <- matrix_network_copy[rowSums(matrix_network_copy) != 0, , drop = FALSE]
    active_insects <- rownames(matrix_network_copy)
    initial_pests <- sum(active_insects %in% insect_roles$Species[insect_roles$Role == "Pest"])
    remaining_pests_list[i] <- remaining_pests_list[i] + initial_pests
    
    for (bird in extinction_sequence) {
      i <- i + 1
      matrix_network_copy <- matrix_network_copy[, colnames(matrix_network_copy) != bird, drop = FALSE] # extinction
      matrix_network_copy <- matrix_network_copy[rowSums(matrix_network_copy) != 0, , drop = FALSE] # drop pest out of control
      active_insects <- rownames(matrix_network_copy)
      
      remaining_pests <- sum(active_insects %in% insect_roles$Species[insect_roles$Role == "Pest"])
      remaining_pests_list[i] <- remaining_pests_list[i] + remaining_pests
    }
  }
  
  remaining_pests_list <- remaining_pests_list / repeat_count
  return(remaining_pests_list)
}

# function to visualize the extinction curve of pest control
save_extinction_curve <- function(pest_list, name, type){
  time_stamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  result_df <- data.frame(
    Step = 0:(length(pest_list) - 1),
    Remaining_Pests = unlist(pest_list)
  )
  p <- ggplot(result_df, aes(x = Step, y = Remaining_Pests)) +
    geom_line() +
    geom_point() +
    labs(title = paste0("Extinction Curve(Pest Control) in ", name, " Site"), x = "Extinction Step", y = "Pests under control") +
    theme_minimal()
  ggsave(paste0("extinction_curve_", type, "_", name, "_", time_stamp, ".jpg"), p)
  
}



### Start data processing
# load data
data <- read.csv("./009MRes Denian Li Project/Complete Diet Data/fwh_consensus.csv")
bird_meta <- read.csv("./009MRes Denian Li Project/Complete Diet Data/Ghana_Mistnetting_wrangled_14052024.csv")
insect_meta <- read.csv("./009MRes Denian Li Project/Data/Pest_annotation_Ghana_060223_cleaned.csv")
bird_meta_update <- read_excel("./009MRes Denian Li Project/Data/sample_BF_dist_cover_allcounts.xlsx", sheet="sample_BF_dist_cover_allcounts")

# Filter the data to keep only arthropods exclude arachnids which are mites and therefore not dietary items
data_filtered <- filter_by_species(data)
data_RAA <- filter_by_RAA(data_filtered) # RAA filter to eliminate bias from small data by threshold of 0.3%

# site_type <- c("Forest Interior", "Forest Edge", "Agriculture")
data_RAA_forest <- filter_by_site(data_RAA, bird_meta_update, "Forest Interior")
data_RAA_edge <- filter_by_site(data_RAA, bird_meta_update, "Forest Edge")
data_RAA_agriculture <- filter_by_site(data_RAA, bird_meta_update, "Agriculture")

# Merge bird forest coverage metadata and insect role metadata
bird_forest_coverage <- merge_bird_forest_coverage(bird_meta) %>%
  mutate(Bird.Type = ifelse(average.Distance.to.Forest.km <= 0, "Forest bird", "Non-Forest bird"))
bird_forest_coverage <- bird_forest_coverage[order(bird_forest_coverage$average.Distance.to.Forest.km), ]
insect_role <- merge_insect_role(insect_meta)

# Combine metadata with the filtered data
data_in_site_type <- list(
  "Total" = combine_data(bird_meta, data_RAA, data), 
  "Forest Interior" = combine_data(bird_meta, data_RAA_forest, data), 
  "Forest Edge" = combine_data(bird_meta, data_RAA_edge, data), 
  "Agriculture" = combine_data(bird_meta, data_RAA_agriculture, data)
  )

# Main loop for data in different site type
for(name in names(data_in_site_type)){
  data_merged <- data_in_site_type[[name]]
  valid_columns <- bird_forest_coverage$BirdLife.Name.CORRECT %in% colnames(data_merged)
  data_merged_reordered <- data_merged[, bird_forest_coverage$BirdLife.Name.CORRECT[valid_columns]]
  data_merged_reordered <- cbind(Insect_species=data_merged[["Insect_species"]], data_merged_reordered)
  bird_insect_web <- data_merged_reordered %>% column_to_rownames(var = 'Insect_species')
  
  # Visualize bipartite plot
  save_forest_category_bipartite(bird_forest_coverage, insect_role, bird_insect_web, name)
  
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
  
  # Save the node and edge data for Gephi
  nodes_o <- paste0("nodes_", gsub(" ", "_", name), ".csv")
  edges_o <- paste0("edges_", gsub(" ", "_", name), ".csv")
  convert_data_for_gephi(data_merged, bird_forest_coverage, nodes_o, edges_o)
  
  # Visualize the extinction curve
  result <- extinction_curve(bird_insect_web, insect_role)
  save_extinction_curve(result, name, "Random")
}


stop()













##############################################################################################################################
results <- statistic_richness(data_in_site_type, bird_forest_coverage)

data_stat <- data.frame(
  Site = rep(c("Total", "Forest Interior", "Forest Edge", "Agriculture"), each = 2),
  Bird = rep(c("Forest bird", "Non-Forest bird"), sites = 4),
  Metrics = unlist(results),
  SE = c(0, 0, 0, 0, 0, 0, 0, 0)
)

p <- ggplot(data_stat, aes(x = Site, y = Metrics, color = Bird, group = Bird)) +
  geom_point(size = 4, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = Metrics - SE, ymax = Metrics + SE), 
                width = 0.2, 
                position = position_dodge(0.5), 
                linewidth = 1) +
  geom_vline(xintercept = 1.5, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "gray") + 
  geom_vline(xintercept = 3.5, linetype = "dashed", color = "gray") + 
  theme_minimal() +
  labs(x = "", y = "Bird species richness") +
  scale_color_manual(values = c("Forest bird" = "#238E23", "Non-Forest bird" = "darkgoldenrod2")) +
  theme(
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    plot.background = element_blank()
  )

time_stamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
name = "species_richness"
ggsave(paste0("bird_static_", name, "_", time_stamp, ".jpg"), p)


data_individual_in_site_type <- list(
  "Total" = data_RAA, 
  "Forest Interior" = data_RAA_forest, 
  "Forest Edge" = data_RAA_edge, 
  "Agriculture" = data_RAA_agriculture
)

# 统计丰度并计算平均值和标准误
statistic_abundance <- function(matrix_list, bird_coverage) {
  results_mean <- list()
  results_se <- list()
  
  # 分类鸟类
  classified_birds <- bird_coverage %>%
    mutate(Bird_Type = ifelse(Distance.to.Forest.km < 0, "Forest bird", "Non-Forest bird"))
  
  for (i in seq_along(matrix_list)) {
    print(i)
    matrix <- matrix_list[[i]]
    non_empty_cols <- colSums(matrix != 0) > 0
    non_empty_col_names <- colnames(matrix)[non_empty_cols]
    
    # 筛选非空鸟类
    non_empty_birds <- classified_birds %>% filter(Sample_number %in% non_empty_col_names)
    
    forest_birds <- non_empty_birds %>% filter(Bird_Type == "Forest bird")
    non_forest_birds <- non_empty_birds %>% filter(Bird_Type == "Non-Forest bird")
    
    # 统计每种鸟类的丰度
    forest_abundance <- forest_birds %>%
      group_by(BirdLife.Name.CORRECT) %>%
      summarise(Abundance = n())
    
    non_forest_abundance <- non_forest_birds %>%
      group_by(BirdLife.Name.CORRECT) %>%
      summarise(Abundance = n())
    
    # 计算平均值和标准误
    forest_mean <- mean(forest_abundance$Abundance)
    if (nrow(forest_abundance) < 2) {forest_se <- 0}
    else {
      forest_se <- sd(forest_abundance$Abundance) / sqrt(nrow(forest_abundance))
    }
    
    non_forest_mean <- mean(non_forest_abundance$Abundance)
    if (nrow(non_forest_abundance) < 2) {non_forest_se <- 0}
    else {
      non_forest_se <- sd(non_forest_abundance$Abundance) / sqrt(nrow(non_forest_abundance))
    }
    
    results_mean <- c(results_mean, list(c(forest_mean, non_forest_mean)))
    results_se <- c(results_se, list(c(forest_se, non_forest_se)))
  }
  
  return(list(mean = results_mean, se = results_se))
}

results <- statistic_abundance(data_individual_in_site_type, bird_meta_update)

# 打印结果
print(results$mean)
print(results$se)

data_stat <- data.frame(
  Site = rep(c("Total", "Forest Interior", "Forest Edge", "Agriculture"), each = 2),
  Bird = rep(c("Forest bird", "Non-Forest bird"), sites = 4),
  Metrics = unlist(results$mean),
  SE = unlist(results$se)
)

p <- ggplot(data_stat, aes(x = Site, y = Metrics, color = Bird, group = Bird)) +
  geom_point(size = 4, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = Metrics - SE, ymax = Metrics + SE), 
                width = 0.2, 
                position = position_dodge(0.5), 
                linewidth = 1) +
  geom_vline(xintercept = 1.5, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "gray") + 
  geom_vline(xintercept = 3.5, linetype = "dashed", color = "gray") + 
  theme_minimal() +
  labs(x = "", y = "Bird abundance") +
  scale_color_manual(values = c("Forest bird" = "#238E23", "Non-Forest bird" = "darkgoldenrod2")) +
  theme(
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    plot.background = element_blank()
  )

time_stamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
name = "abundance"
ggsave(paste0("bird_static_", name, "_", time_stamp, ".jpg"), p)

##############################################################################################################################














##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################


















web <- data_in_site_type[["Total"]] %>%
  column_to_rownames(var = colnames(data_in_site_type[["Total"]])[1])

metrics <- c('weighted nestedness','links per species','linkage density', 'generality', 'H2' )

for(metrc in metrics){
  nestedness_result <- networklevel(web, index = metrc, weighted = FALSE)
  print(nestedness_result)
}


birds <- unique(bird_forest_coverage$Bird.Type)
metrics <- c('weighted nestedness','links per species','linkage density', 'generality', 'H2')

obs <- unlist(networklevel(web, index=c(metrics[[1]]), weighted = F))

nulls <- r2dexternal(100, empty(web), abun.higher = NULL, abun.lower = NULL) # replace 10 with 1,000, this is a lower number for the script to run quicker!

null <- unlist(sapply(nulls, networklevel, index=c(metrics[[1]])))

praw <- sum(null>obs) / length(null)
p_values <- rep(NA, 2)
p_values[1] <- ifelse(praw > 0.5, 1-praw, praw)


means <- mean(null)
est_HL[counter] <- means
est_LL[counter] <- NA
store_obs[counter] <- obs

means <- mean(null[1])
obs_1 <- obs[1]

est_HL[counter] <- means[1]
est_LL[counter] <- means[2]
store_obs[counter] <- obs[1]


# no_values <- length(insects)*length(metrics)
# 
# p_values <- rep(NA, no_values)
# est_HL <- rep(NA, no_values)
# est_LL <- rep(NA, no_values)
# 
# store_time   <- rep(NA, no_values)
# store_insect <- rep(NA, no_values)
# #store_site   <- rep(NA, no_values)
# store_obs <- rep(NA, no_values)
# store_metric <- rep(NA, no_values)
# counter <-1
# 
for(metric in metrics){
  for(bird in birds){
    
    
    LD1 <- as.data.frame(t(LD1))
    
    if(dim(LD1)[2] > 1) {
      
      obs <- unlist(networklevel(LD1, index=c(metric), weighted = F))
      
      nulls <- r2dexternal(10, empty(LD1), abun.higher = NULL, abun.lower = NULL) # replace 10 with 1,000, this is a lower number for the script to run quicker!
      
      null <- unlist(sapply(nulls, networklevel, index=c(metric)))
      
      praw <- sum(null>obs) / length(null)
      p_values[counter] <- ifelse(praw > 0.5, 1-praw, praw)  
      
      
      means <- mean(null)
      est_HL[counter] <- means
      est_LL[counter] <- NA
      store_obs[counter] <- obs
      
      means <- mean(null[1])
      obs_1 <- obs[1]
      
      est_HL[counter] <- means[1]
      est_LL[counter] <- means[2]
      store_obs[counter] <- obs[1] 
      #}
      
    }
    
    store_time[counter] <- time
    store_insect[counter] <- insect
    #store_site[counter] <- site
    store_metric[counter] <- metric
    
    
    cat(insect,"  ",time,"  ",p_values[counter],"  ",
        est_HL[counter],"  ",
        est_LL[counter],"  ",
        metric, '  ', "\n")
    counter <- counter +1
  }
}

# 创建灭绝曲线
# 将第一列作为行名
web <- data_in_site_type[[2]]

# 移除第一列
web <- web[, -1]

# 确保数据为数值型矩阵
web <- as.matrix(web)
extinction_results <- second.extinct(web, participant = "higher", method = "random")

# 绘制灭绝曲线
plot(extinction_results, main = "Bird-Insect Bipartite Network Extinction Curve",
     xlab = "Proportion of Birds Removed", ylab = "Proportion of Insects Remaining")





