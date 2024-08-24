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
if (!requireNamespace("iNEXT", quietly = TRUE)) install.packages("iNEXT")

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
library(igraph)
library(iNEXT)


### Start define function
# The code to filter arthropods and ambiguous species by Ariane
filter_by_species <- function(pre_data){
  pre_data$species <- gsub("sp.*", "sp.", pre_data$species)
  arthropoda_and_species_only <- pre_data %>%
    dplyr::filter(phylum == "Arthropoda") %>% # only keep arthropods
    dplyr::filter(order != "Mesostigmata"&  # exclude arachnids which are mites and therefore not dietary items
                    order != "Sarcoptiformes"&
                    order != "Trombidiformes"&
                    species != "Pseudacanthotermes militaris") %>%
    filter(!is.na(species)) %>%  # only keep OTUs identified to genus or species levelzhe
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

filter_by_site <- function(pre_data, meta_data, site_names) {
  filtered_meta_data <- filter(meta_data, Site.Type %in% site_names)
  columns_to_select <- intersect(names(pre_data), filtered_meta_data$Specimen.Number)
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
  
  data_merged$Insect_species <- gsub("sp.*", "sp.", data_merged$Insect_species)
  data_merged <- data_merged %>%
    group_by(Insect_species) %>%
    summarise(across(everything(), sum, na.rm = TRUE))
  
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
  # calculate the forest dependency of the first row
  first_rows <- meta %>%
    group_by(BirdLife.Name.CORRECT) %>%
    slice(1) %>%
    mutate(ForestDependency = Bird.Forest + 0.5 * Bird.Edge - Bird.Agriculture)
  
  # calculate average distance
  forest_coverage <- meta %>%
    group_by(BirdLife.Name.CORRECT) %>%
    summarise(average.Distance.to.Forest.km = mean(Bird.Distance.to.Forest.km))
  
  # get canopy
  canopy_data <- first_rows %>%
    select(BirdLife.Name.CORRECT, Bird.Forest.Cover.200m)
  
  # get abundance
  abundance_data <- first_rows %>%
    select(BirdLife.Name.CORRECT, Bird.Rarity.Score)
  
  # merge
  result <- first_rows %>%
    select(BirdLife.Name.CORRECT, ForestDependency) %>%
    left_join(forest_coverage, by = "BirdLife.Name.CORRECT") %>%
    left_join(canopy_data, by = "BirdLife.Name.CORRECT") %>%
    left_join(abundance_data, by = "BirdLife.Name.CORRECT")
  
  return(result)
}

# Calculate Robustness of pest annotation
calculate_robustness <- function(df, role_col = "Role", target_cols = c("Is the genus a plant pest?", "Known pest?", "Agricultural pest in SSA?", "Pest in Ghana?"), robustness_col = "Robustness") {
  for (i in 1:nrow(df)) {
    # Update Role column based on specific conditions
    if (df[i, role_col] == "Disease vector") {
      df[i, role_col] <- "Vector of animal disease"
    } else if (df[i, role_col] == "Forestry pest") {
      df[i, role_col] <- "Pest"
    } else if (df[i, role_col] == "Natual enemy") {
      df[i, role_col] <- "Natural enemy"
    } else if (df[i, role_col] == "Disease Vector") {
      df[i, role_col] <- "Vector of animal disease"
    } else if (df[i, role_col] == "Unknown") {
      df[i, role_col] <- "Non-pest"
    } 
    
    # Skip the row if Robustness column already has a value
    # if (!is.na(df[i, robustness_col])) next
    
    # Update Robustness column if Role is "Pest"
    # if (df[i, role_col] == "Pest") {
    #   yes_count <- sum(df[i, target_cols] == "yes", na.rm = TRUE)
    #   df[i, robustness_col] <- ifelse(yes_count > 0, yes_count, 1)
    # }
  }
  
  df_extracted <- df %>%
    select(Species = species.corrected, Role, Robustness)
  
  return(df_extracted)
}

# produce a table containing insect role
merge_insect_role <- function(meta) {
  # discard bold annotation
  meta$Species <- gsub("sp.*", "sp.", meta$Species)
  role <- meta %>%
    filter(!is.na(Species) & Species != "" & !is.na(Role) & Role != "Unknown" & Role != "") %>%
    select(Species, Role, Robustness)
  # discard copy species in order
  role_priority <- c("Pest", "Natural enemy", "Vector of animal disease", "Vector of plant disease", "Non-pest")
  role <- role %>%
    mutate(Role = factor(Role, levels = role_priority)) %>%
    group_by(Species) %>%
    summarise(Role = first(sort(Role, na.last = TRUE)),
              Robustness = first(Robustness)) %>%
    ungroup()
  role$Role <- as.character(role$Role)
  return(role)
}

update_insect_data <- function(df1, df2) {
  role_priority <- c("Pest", "Natural enemy", "Vector of animal disease", "Vector of plant disease", "Non-pest")
  combined_df <- bind_rows(df1, df2)
  processed_df <- combined_df %>%
    mutate(Role = factor(Role, levels = role_priority)) %>%
    arrange(Species, Role) %>%
    group_by(Species) %>%
    slice_min(order_by = Role, with_ties = FALSE) %>%
    ungroup()
  processed_df$Role <- as.character(processed_df$Role)
  return(processed_df)
}

# 扩展的 web_feature 函数
web_feature <- function(bird_insect_web, bird_forest_coverage, insect_role_all){
  # 将鸟名和昆虫名从bird_insect_web中提取出来
  bird_names <- colnames(bird_insect_web)
  insect_names <- rownames(bird_insect_web)
  
  # 匹配并统计森林鸟类和非森林鸟类
  bird_forest_coverage_filtered <- bird_forest_coverage %>%
    filter(BirdLife.Name.CORRECT %in% bird_names)
  
  forest_bird_names <- bird_forest_coverage_filtered %>% 
    filter(Bird.Type == "Forest bird") %>% 
    pull(BirdLife.Name.CORRECT)
  
  non_forest_bird_names <- bird_forest_coverage_filtered %>% 
    filter(Bird.Type == "Non-Forest bird") %>% 
    pull(BirdLife.Name.CORRECT)
  
  # 匹配并统计害虫和非害虫
  insect_role_filtered <- insect_role_all %>%
    filter(Species %in% insect_names)
  
  pest_insect_names <- insect_role_filtered %>%
    filter(Role == "Pest") %>%
    pull(Species)
  
  non_pest_insect_names <- insect_role_filtered %>%
    filter(Role != "Pest") %>%
    pull(Species)
  
  # 计算交互强度（和，不是度）
  forest_pest_interaction_strength <- sum(bird_insect_web[rownames(bird_insect_web) %in% pest_insect_names, colnames(bird_insect_web) %in% forest_bird_names])
  forest_non_pest_interaction_strength <- sum(bird_insect_web[rownames(bird_insect_web) %in% non_pest_insect_names, colnames(bird_insect_web) %in% forest_bird_names])
  
  non_forest_pest_interaction_strength <- sum(bird_insect_web[rownames(bird_insect_web) %in% pest_insect_names, colnames(bird_insect_web) %in% non_forest_bird_names])
  non_forest_non_pest_interaction_strength <- sum(bird_insect_web[rownames(bird_insect_web) %in% non_pest_insect_names, colnames(bird_insect_web) %in% non_forest_bird_names])
  
  # 统计数量
  forest_bird_count <- length(forest_bird_names)
  non_forest_bird_count <- length(non_forest_bird_names)
  pest_count <- length(pest_insect_names)
  total_insects_count <- length(insect_names)
  
  # 打印结果
  cat("Forest bird count:", forest_bird_count, "\n")
  cat("Non-Forest bird count:", non_forest_bird_count, "\n")
  cat("Pest count:", pest_count, "\n")
  cat("Total insect count:", total_insects_count, "\n")
  cat("Forest bird - Pest interaction strength:", forest_pest_interaction_strength, "\n")
  cat("Forest bird - Non-Pest interaction strength:", forest_non_pest_interaction_strength, "\n")
  cat("Non-Forest bird - Pest interaction strength:", non_forest_pest_interaction_strength, "\n")
  cat("Non-Forest bird - Non-Pest interaction strength:", non_forest_non_pest_interaction_strength, "\n")
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

  title(main = paste("Bipartite Network Plot for", name))
  # legend("right", legend = c("Forest Birds", "Non-Forest Birds", "Pests", "Other Insects"),
  #        col = c("#238E23", "darkgoldenrod2", "red", "grey"), pch = 15, cex = 0.8)
  
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
  ggsave(paste0("proportion_", name, "_", time_stamp, ".jpg"), pp,width = 6.67, height = 1.5, units = "in")
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

calculate_pest_percentage <- function(bird_species, web, insect_role) {
  bird_interactions <- web[, bird_species]
  non_zero_indices <- which(bird_interactions > 0)
  bird_interactions <- bird_interactions[non_zero_indices]
  insect_species_eaten <- rownames(web)[non_zero_indices]
  insect_roles <- insect_role[insect_role$Species %in% insect_species_eaten, ]
  pest_species_indices <- match(insect_roles$Species[insect_roles$Role == "Pest"], insect_species_eaten)
  pest_species_indices <- pest_species_indices[!is.na(pest_species_indices)]
  
  if (length(pest_species_indices) == 0) {
    return(list(
      pest_interaction_percentage = 0,
      pest_species_percentage = 0
    ))
  }
  
  pest_interactions <- sum(bird_interactions[pest_species_indices])
  total_interactions <- sum(bird_interactions)
  pest_species_count <- length(pest_species_indices)
  total_species_count <- length(unique(insect_roles$Species))

  pest_interaction_percentage <- (pest_interactions / total_interactions) * 100
  pest_species_percentage <- (pest_species_count / total_species_count) * 100
  
  return(c(
    pest_interaction_percentage = pest_interaction_percentage,
    pest_species_percentage = pest_species_percentage
  ))
}



calculate_bird_feature <- function(site_name, bird_insect_web, bird_insect_web_weighted, insect_role){
  # calculate pest proportion in individual bird species
  result_list <- list()
  for (bird_name in colnames(bird_insect_web)) {
    result_list[[bird_name]] <- calculate_pest_percentage(bird_name, bird_insect_web, insect_role) # species result 1
  }
  bird_species_pest_proportion <- as.data.frame(do.call(rbind, result_list)) 
  bird_species_professionalism <- specieslevel(bird_insect_web, level = "higher", index = "d") # species result 2
  
  result_modularity <- calculate_modularity(bird_insect_web) # network result 1
  result_nested <- networklevel(bird_insect_web, index = "NODF") # network result 2

  # bipartite projects to unipartite
  bird_insect_matrix <- as.matrix(bird_insect_web)
  network_specialization <- H2fun(bird_insect_web_weighted) # network result 3
  g <- graph_from_biadjacency_matrix(bird_insect_matrix)
  bird_projection <- bipartite_projection(g)
  bird_network <- bird_projection$proj2
  network_diameter <- diameter(bird_network, directed = FALSE, unconnected = FALSE) 
  total_birds <- length(V(bird_network))
  diet_overlap <- (network_diameter - 1) / (total_birds - 2) # network result 4
  diet_overlap <- networklevel(bird_insect_web, index = "niche overlap")[["niche.overlap.HL"]]
  
  bird_nodes <- V(bird_network)
  total_birds <- length(bird_nodes)
  
  # 使用 closeness 函数计算每个鸟类节点的紧密度
  closeness_scores <- closeness(bird_network, mode = "all", normalized = TRUE)
  
  # 将结果转换为数据框
  bird_closeness <- data.frame(
    closeness = closeness_scores,
    row.names = names(bird_nodes)
  )
  
  merged_df <- merge(bird_species_pest_proportion, bird_species_professionalism, by = "row.names")
  rownames(merged_df) <- merged_df$Row.names
  merged_df <- merged_df[, -1]
  merged_df <- merge(merged_df, bird_closeness, by = "row.names")
  rownames(merged_df) <- merged_df$Row.names
  merged_df <- merged_df[, -1]
  
  
  return(list(
    bird_metrics = merged_df,
    network_metrics = c(
      network_name = site_name,
      modularity = result_modularity,
      nestedness = result_nested[["NODF"]],
      specialization = network_specialization[["H2"]],
      diameter = diet_overlap
    )
  ))
}




calculate_modularity <- function(bird_insect_web) {
  # 将 bird_insect_web 转换为 igraph 对象，不使用权重
  g <- graph_from_incidence_matrix(bird_insect_web, weighted = NULL)
  
  # 使用适合二部图的社区检测算法：cluster_walktrap
  community <- cluster_walktrap(g)
  
  # 计算模块度
  modularity_value <- modularity(community)
  
  return(modularity_value)
  
  # 如果想获取每个节点的社区
  # return(membership(community))
}


check_correlation <- function(y) {
  # 检查输入是否为向量且非列表类型
  if (!is.vector(y) || is.list(y)) {
    stop("输入必须是一个数值向量。")
  }
  
  # 数据准备和清理：移除缺失值
  if (any(is.na(y))) {
    print("数据中存在缺失值，移除这些缺失值。")
    y <- na.omit(y)
  }
  
  # 创建向量 x 为 1 到 y 的长度
  x <- 1:length(y)
  
  # 初步探索分析：绘制散点图
  ggplot(data = data.frame(x, y), aes(x = x, y = y)) +
    geom_point() +
    labs(title = "Scatter plot of x and y", x = "Index (x)", y = "Value (y)") +
    theme_minimal()
  
  # 计算相关系数
  correlation <- cor(x, y, method = "pearson")
  print(paste("Pearson correlation coefficient:", correlation))
  
  # 统计显著性检验
  cor_test <- cor.test(x, y, method = "pearson")
  print(cor_test)
  
  # 结果解读
  if (cor_test$p.value < 0.10) {
    print("相关性显著")
  } else {
    print("相关性不显著")
  }
}

extract_beak_length <- function(abundance, bird_meta, bird_name_col, beak_length_col) {
  # 检查输入是否为向量
  if (!is.vector(abundance) || is.list(abundance)) {
    stop("abundance输入必须是一个数值向量。")
  }
  
  # 使用sapply提取每个种名的喙长度数据
  beak_length_nares <- sapply(unique(abundance), function(species) {
    # 在 bird_meta 中找到匹配的第一行
    index <- match(species, bird_meta[[bird_name_col]])
    # 如果存在匹配的种名，提取数据
    if (!is.na(index)) {
      return(bird_meta[[beak_length_col]][index])
    } else {
      return(NA)
    }
  })
  
  # 移除名称属性
  beak_length_nares <- unname(beak_length_nares)
  
  return(beak_length_nares)
}

probabilistic_selection <- function(original_sequence, prob_left = 0.8) {
  new_sequence <- vector()
  select_individual <- function(sequence, prob_left) {
    while (length(sequence) > 1) {
      n <- length(sequence)
      mid <- floor(n / 2)
      if (runif(1) < prob_left) {
        sequence <- sequence[1:mid]
      } else {
        sequence <- sequence[(mid + 1):n]
      }
    }
    return(sequence)
  }
  
  while (length(original_sequence) > 0) {
    selected <- select_individual(original_sequence, prob_left)
    new_sequence <- c(new_sequence, selected)
    original_sequence <- original_sequence[original_sequence != selected]
  }
  return(new_sequence)
}

# function to draw an extinction curve
extinction_curve <- function(matrix_network, insect_roles, repeat_count = 1000, extinction_sequence = NULL) {
  matrix_network <- matrix_network[rowSums(matrix_network) != 0, colSums(matrix_network) != 0, drop = FALSE]
  random_flag <- FALSE
  if (is.null(extinction_sequence)) {
    extinction_sequence <- sample(colnames(matrix_network))
    random_flag <- TRUE
  }
  extinction_sequence <- extinction_sequence[extinction_sequence %in% colnames(matrix_network)]
  remaining_pests_matrix <- matrix(0, nrow = repeat_count, ncol = length(extinction_sequence) + 1)
  
  for (t in 1:repeat_count){
    if (t %% (repeat_count / 10) == 0) {
      print(paste0(t/repeat_count*100, "%"))
    }
    i <- 1
    if(random_flag){
      extinction_sequence <- sample(extinction_sequence)
    }
    else{
      extinction_sequence <- probabilistic_selection(extinction_sequence, 0.8)
    }
    matrix_network_copy <- matrix_network
    matrix_network_copy <- matrix_network_copy[rowSums(matrix_network_copy) != 0, , drop = FALSE]
    active_insects <- rownames(matrix_network_copy)
    initial_pests <- sum(active_insects %in% insect_roles$Species[insect_roles$Role == "Pest"])
    
    remaining_pests_matrix[t, i] <- initial_pests
    
    # Start extinction for one sequence
    for (bird in extinction_sequence) {
      i <- i + 1
      matrix_network_copy <- matrix_network_copy[, colnames(matrix_network_copy) != bird, drop = FALSE] # extinction
      matrix_network_copy <- matrix_network_copy[rowSums(matrix_network_copy) != 0, , drop = FALSE] # drop pest out of control
      active_insects <- rownames(matrix_network_copy)
      remaining_pests <- sum(active_insects %in% insect_roles$Species[insect_roles$Role == "Pest"])
      remaining_pests_matrix[t, i] <- remaining_pests
    }
  }
  mean_remaining_pests <- colMeans(remaining_pests_matrix)
  sd_remaining_pests <- apply(remaining_pests_matrix, 2, sd)
  ci <- 1.96 * (sd_remaining_pests / sqrt(repeat_count))
  return(list(mean = mean_remaining_pests, error = ci))
}

plot_extinction_curve <- function(result) {
  # 提取均值和误差
  mean_values <- result$mean
  error_values <- result$error
  
  # 创建数据框用于绘图
  data <- data.frame(
    Step = 0:(length(mean_values) - 1),
    Mean = mean_values,
    Error = error_values
  )
  
  # 使用 ggplot2 进行可视化
  p <- ggplot(data, aes(x = Step, y = Mean)) +
    geom_line(color = "blue") +
    geom_point(color = "blue") +
    geom_errorbar(aes(ymin = Mean - Error, ymax = Mean + Error), width = 0.2, color = "red") +
    labs(title = "Extinction Curve",
         x = "Extinction Step",
         y = "Remaining Pests") +
    theme_minimal()
  
  print(p)
}



save_extinction_curve_by_metric <- function(extinction_data, plot_name) {
  time_stamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  
  normalize <- function(x) {
    return ((x - min(x)) / (max(x) - min(x)))
  }
  
  process_data <- function(data, name) {
    mean_data <- data$mean
    error_data <- data$error
    min_mean <- min(mean_data)
    max_mean <- max(mean_data)
    scale_factor <- max_mean - min_mean
    
    normalized_mean <- (mean_data - min_mean) / scale_factor
    normalized_error <- error_data / scale_factor
    
    df <- data.frame(
      Step = 0:(length(mean_data) - 1),
      Mean_Remaining_Pests = normalized_mean,
      Error_Remaining_Pests = normalized_error,
      Name = name
    )
    df$Step <- normalize(df$Step)
    return(df)
  }
  
  combined_df <- do.call(rbind, lapply(names(extinction_data), function(name) {
    process_data(extinction_data[[name]], name)
  }))
  
  colors <- c("Agriculture" = "darkgoldenrod2", 
              "Forest Interior" = "#238E23", 
              "Forest Edge" = "#91A724", 
              "Total" = "black")
  
  p <- ggplot(combined_df, aes(x = Step, y = Mean_Remaining_Pests, color = Name, fill = Name)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = Mean_Remaining_Pests - Error_Remaining_Pests, ymax = Mean_Remaining_Pests + Error_Remaining_Pests), alpha = 0, color = NA) +
    annotate("segment", x = 0, xend = 1, y = 1, yend = 0, color = "grey", linetype = "dashed") +
    labs(title = paste("Resilience of natural pest control under", plot_name, "order of bird"),
         x = "Proportion of predator species lost",
         y = "Proportion of predator-pest interactions intact") +
    theme_minimal() +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    theme(
      legend.title = element_blank(),
      legend.position = "right",
      legend.box.margin = margin(0, 10, 0, 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = "white"),
      axis.line = element_line(color = "black")
    )
  
  # 保存图形
  ggsave(paste0("extinction_curves_", plot_name, "_", time_stamp, ".jpg"), p,
         width = 7.5, height = 5, units = "in")
}

calculate_auc <- function(y){
  x <- seq(1, length(y))
  auc <- sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
  max_auc <- (max(x) - min(x)) *  max(y)  # 使用实际的 y 最大可能值
  return(auc / max_auc)
}

plot_auc_half_life <- function(data, metric_name) {
  # 筛选数据
  data_metric <- data %>% filter(metrics == metric_name)
  
  # 区分均值和半衰期数据
  avg_data <- data_metric %>% filter(grepl("avg", variable)) %>% mutate(value = value * 100)
  half_data <- data_metric %>% filter(grepl("half", variable))
  
  # 自定义颜色
  custom_colors <- c(
    "agriculture" = "darkgoldenrod2", 
    "forest" = "#238E23", 
    "edge" = "#91A724", 
    "total" = "black"
  )
  
  # 根据变量名的前缀映射颜色
  avg_data$color_group <- sub("_.*", "", avg_data$variable)
  half_data$color_group <- sub("_.*", "", half_data$variable)
  
  # 创建图表
  p <- ggplot() +
    geom_bar(data = avg_data, aes(x = variable, y = value, fill = color_group), stat = "identity", position = position_dodge(), alpha = 0.7) +
    geom_line(data = half_data, aes(x = variable, y = value, group = 1, color = "Half-life (%)"), size = 1) +
    geom_point(data = half_data, aes(x = variable, y = value, color = "Half-life (%)")) +
    scale_y_continuous(name = "Area Under Curve (%)", sec.axis = sec_axis(~ ., name = "Half-life (%)")) +
    scale_fill_manual(values = custom_colors) +
    scale_color_manual(values = c("Half-life (%)" = "blue")) +
    labs(title = paste("Metrics:", metric_name), x = "Variable", y = "Value") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      panel.background = element_rect(fill = "white", color = NA),  # 设置面板背景为白色
      plot.background = element_rect(fill = "white", color = NA),   # 设置绘图背景为白色
      panel.grid = element_blank(),  # 取消网格线
      axis.line = element_line(color = "black")  # 添加轴线
    )
  
  return(p)
}

# Plot all metrics of species
plot_species_metric <- function(species_result){
  convert_species_result_to_df <- function(species_result) {
    data_list <- list()
    for (site_type in names(species_result)) {
      # 提取当前 site_type 对应的数据框
      df <- species_result[[site_type]]
      
      # 创建一个新的数据框，包括 species 列和 site_type 列
      df_new <- data.frame(
        species = rownames(df),
        site_type = site_type,
        df
      )
      
      # 将处理后的数据框添加到列表中
      data_list[[site_type]] <- df_new
    }
    
    # 将所有小数据框合并成一个大的数据框
    final_df <- do.call(rbind, data_list)
    
    return(final_df)
  }
  
  # 使用该函数将 species_result 转换为数据框
  final_df <- convert_species_result_to_df(species_result)
  
  # 确保所有列是数值型
  final_df <- final_df %>%
    mutate(across(c(pest_interaction_percentage, pest_species_percentage, d, closeness), as.numeric))
  
  final_df[["pest_species_percentage"]] <- final_df[["pest_species_percentage"]] / 100
  final_df[["pest_interaction_percentage"]] <- final_df[["pest_interaction_percentage"]] / 100
  
  # 将需要绘图的四大类提取出来，并将数据转换为长格式
  plot_data <- final_df %>%
    select(species, site_type, pest_interaction_percentage, pest_species_percentage, d, closeness) %>%
    melt(id.vars = c("species", "site_type"), 
         measure.vars = c("pest_interaction_percentage", "pest_species_percentage", "d", "closeness"),
         variable.name = "feature", 
         value.name = "value")
  
  # 将 site_type 转换为因子，并设置因子级别顺序
  plot_data$site_type <- factor(plot_data$site_type, levels = c("Forest Interior", "Forest Edge", "Agriculture"))
  
  # 绘制箱线图，并按照 site_type 排序
  feature_levels <- levels(plot_data$feature)
  vline_positions <- seq(1.5, length(feature_levels) - 0.5, by = 1)
  
  p <- ggplot(plot_data, aes(x = feature, y = value, fill = site_type)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, position = position_dodge(width = 0.8)) + # 设置箱线图宽度和位置
    scale_fill_manual(values = c("Forest Interior" = "#238E23", "Forest Edge" = "#91A724", "Agriculture" = "darkgoldenrod2")) + # 设置填充颜色
    labs(title = "Standardized Box Plot of Bird Species Features by Site Type",
         x = "Feature",
         y = "Value of Species Metrics",
       fill = "Site Type") +  # 修改y轴标签名称
    scale_x_discrete(labels = c("Pest Dietary (PD)", "Pest Dietary (PD)", "Species Specialisation (d')", "Diet Overlap (DO)")) +  # 修改x轴标签顺序和名称
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + # 设置y轴范围，并将x轴移动到y=0
    theme_minimal() +
    theme(
      legend.position = "right",                  # 将图例放在右侧
      panel.grid = element_blank(),               # 去掉网格线
      axis.line.y = element_line(color = "black"),  # 增加y轴线
      axis.ticks.y = element_line(color = "black"), # 增加y轴刻度
      axis.line.x = element_line(color = "black", size = 0.5), # 增加x轴线
      panel.border = element_blank()              # 去掉边框
    ) +
    geom_vline(xintercept = vline_positions, linetype = "dashed", color = "grey", size = 0.7) # 添加每个大类之间的虚线
  
  
  print(p)
}



# Plot the all metrics of network
plot_network_metric <- function(basic_network_analysis) {
  colors <- c(
    "Forest Interior" = "#238E23",    # 绿色
    "Forest Edge" = "#91A724",        # 浅绿黄色
    "Agriculture" = "darkgoldenrod2"  # 暗金黄色
  )
  
  df <- basic_network_analysis %>%
    gather(key = "Metric", value = "Value", -network_name) %>%
    mutate(LandUseType = case_when(
      network_name == "Forest Interior" ~ "Forest Interior",
      network_name == "Forest Edge" ~ "Forest Edge",
      network_name == "Agriculture" ~ "Agriculture"
    ))
  
  # 将 Metric 设置为因子，并指定显示顺序
  df$Metric <- factor(df$Metric, levels = c("diameter", "modularity", "specialization", "nestedness", "auc", "half_life"))
  
  # 将 LandUseType 设置为因子，并指定顺序
  df$LandUseType <- factor(df$LandUseType, levels = c("Forest Interior", "Forest Edge", "Agriculture"))
  
  df$Value <- as.numeric(df$Value)
  
  # 数据标准化：按每个指标在不同用地类型内的值进行标准化
  df <- df %>%
    group_by(Metric) %>%
    mutate(StandardizedValue = (Value - mean(Value)) / sd(Value)) %>%
    ungroup()
  
  # 绘制条形图
  p <- ggplot(df, aes(x = Metric, y = StandardizedValue, fill = LandUseType)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black", width = 0.7) +  # 设置条形图宽度和位置
    geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 0.7) +  # 添加 y=0 的辅助线
    scale_fill_manual(values = colors) +  # 使用指定的颜色
    theme_minimal() +
    labs(title = "Changes in Metrics Across Different Land Use Types",
         x = "Metrics",
         y = "Standardized Value",
         fill = "Site Type") +
    scale_x_discrete(labels = c(
      "diameter" = expression(DO), 
      "modularity" = expression(Q), 
      "specialization" = expression(H[2] * "'"), 
      "nestedness" = expression(N), 
      "auc" = expression(AUC), 
      "half_life" = expression(t["1/2"])
    )) +
    theme(
      axis.line = element_line(color = "black"),  # 添加x和y轴线
      panel.grid = element_blank(),  # 移除网格线
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.ticks = element_line(color = "black"),  # 添加轴刻度
      axis.text = element_text(size = 10)
    ) +
    geom_vline(xintercept = seq(1.5, 5.5, by = 1), linetype = "dashed", color = "grey", size = 0.7)  # 添加每个大类之间的分割虚线
  
  # 显示图表
  print(p)
}

plot_network_metric(basic_network_analysis) 

# End of function defination

##############################################################################################################################
### Start data processing
# load data
data <- read.csv("./009MRes Denian Li Project/Complete Diet Data/fwh_consensus.csv")
bird_meta <- read.csv("./009MRes Denian Li Project/Complete Diet Data/Ghana_Mistnetting_wrangled_19062024.csv")
insect_meta <- read.csv("./009MRes Denian Li Project/Data/Pest_annotation_Ghana_060223_cleaned.csv")
insect_supply <- read_excel("./009MRes Denian Li Project/Pest annotation 22072024.xlsx", sheet="pest annotation")

# Filter the data to keep only arthropods exclude arachnids which are mites and therefore not dietary items
data_filtered <- filter_by_species(data)
data_RAA <- filter_by_RAA(data_filtered) # RAA filter to eliminate bias from small data by threshold of 0.3%

# site_type <- c("Forest Interior", "Forest Edge", "Agriculture")
data_RAA_forest <- filter_by_site(data_RAA, bird_meta, "Forest Interior")
data_RAA_edge <- filter_by_site(data_RAA, bird_meta, "Forest Edge")
data_RAA_agriculture <- filter_by_site(data_RAA, bird_meta, "Agriculture")
data_RAA_forest_total <- filter_by_site(data_RAA, bird_meta, c("Forest Interior", "Forest Edge"))

# Merge bird forest coverage metadata and insect role metadata
bird_forest_coverage <- merge_bird_forest_coverage(bird_meta) %>%
  mutate(Bird.Type = ifelse(average.Distance.to.Forest.km <= 0, "Forest bird", "Non-Forest bird"))
bird_forest_coverage <- bird_forest_coverage[order(bird_forest_coverage$average.Distance.to.Forest.km), ]
insect_role_all <- update_insect_data(merge_insect_role(insect_meta), calculate_robustness(insect_supply))
insect_role_confirmed <- insect_role_all %>%
  mutate(Role = ifelse(!is.na(Robustness) & Robustness != 4, "Non-pest", Role))

# Combine metadata with the filtered data
data_in_site_type <- list(
  # "Total" = combine_data(bird_meta, data_RAA, data),
  "Forest Interior" = combine_data(bird_meta, data_RAA_forest, data),
  "Forest Edge" = combine_data(bird_meta, data_RAA_edge, data),
  "Agriculture" = combine_data(bird_meta, data_RAA_agriculture, data)
  )
# data_in_site_type <- list(
#   "Forest" = combine_data(bird_meta, data_RAA_forest_total, data), 
#   "Agriculture" = combine_data(bird_meta, data_RAA_agriculture, data)
# )
results_list <- list()


# Main loop for data in different site type
is_all <- TRUE
# for(insect_role in list(insect_role_all, insect_role_confirmed)){
for(insect_role in list(insect_role_all)){
  results_pest_control_df <- data.frame(matrix(ncol = 9, nrow = 0), stringsAsFactors = FALSE) # make empty data frame to store stability pest control result
  names(results_pest_control_df) <- c("metrics", "total_avg", "total_half", "forest_avg", "forest_half", "edge_avg", "edge_half", "agriculture_avg", "agriculture_half")
  species_result <- list()
  basic_network_analysis <- data.frame(
    network_name = character(), 
    modularity = numeric(),
    nestedness = numeric(),
    specialization = numeric(),
    diameter = numeric(),
    stringsAsFactors = FALSE
  )
  
  for(name in names(data_in_site_type)){
    print(name)
    data_merged <- data_in_site_type[[name]]
    valid_columns <- bird_forest_coverage$BirdLife.Name.CORRECT %in% colnames(data_merged)
    data_merged_reordered <- data_merged[, bird_forest_coverage$BirdLife.Name.CORRECT[valid_columns]]
    data_merged_reordered <- cbind(Insect_species=data_merged[["Insect_species"]], data_merged_reordered)
    bird_insect_web_weight <- data_merged_reordered %>% column_to_rownames(var = 'Insect_species')
    bird_insect_web <- bird_insect_web_weight
    bird_insect_web[bird_insect_web != 0] <- 1
    
    # 将数据框转换为长格式
    data_long <- melt(bird_insect_web, varnames = c("Insect", "Bird"), value.name = "value")
    
    # 去除值为0的行
    data_long <- data_long[data_long$value > 0, ]
    
    # 重命名列名，符合所需格式
    colnames(data_long) <- c("to", "from", "value")
    
    # 查看结果
    print(data_long)
    
    stop()
    
    # # calculate the network feature
    # basic_analysis <- calculate_bird_feature(name, bird_insect_web, bird_insect_web_weight, insect_role)
    # species_result[[name]] <- basic_analysis$bird_metrics # used for draw
    # network_feature <- basic_analysis$network_metrics
    # basic_network_analysis <- rbind(basic_network_analysis, as.data.frame(t(network_feature), stringsAsFactors = FALSE))
    # 
    # # # statistics pest distribution
    # # pest_species <- insect_role$Species[insect_role$Role == "Pest"]
    # # pest_distribution[[length(pest_distribution) + 1]] <- sum(pest_species %in% rownames(bird_insect_web))
    # 
    # # Visualize bipartite plot
    # save_forest_category_bipartite(bird_forest_coverage, insect_role, bird_insect_web, name)
    # 
    # # Visualize proportion bar chart
    # forest_birds <- bird_forest_coverage %>%
    #   filter(average.Distance.to.Forest.km <= 0) %>%
    #   pull(BirdLife.Name.CORRECT) %>%
    #   intersect(colnames(bird_insect_web))
    # 
    # non_forest_birds <- bird_forest_coverage %>%
    #   filter(average.Distance.to.Forest.km > 0) %>%
    #   pull(BirdLife.Name.CORRECT) %>%
    #   intersect(colnames(bird_insect_web))
    # 
    # role_groups <- split(insect_role, insect_role$Role)
    # bird_insect_web_by_role <- list()
    # for(role in names(role_groups)) {
    #   species_in_role <- role_groups[[role]]$Species
    #   bird_insect_web_by_role[[role]] <- bird_insect_web %>%
    #     filter(row.names(.) %in% species_in_role)
    # }
    # 
    # 
    # # Plot the proportion of only pest and non-pest
    # static_data <- cbind(
    #   c("Pest", "Non-pest"),
    #   rbind(count_nonzero(bird_insect_web_by_role$Pest, forest_birds, non_forest_birds),
    #         count_nonzero(bird_insect_web_by_role$`Non-pest`, forest_birds, non_forest_birds)))
    # colnames(static_data) <- c("Insect Role", "Forest birds only", "Shared", "Non-forest birds only")
    # static_df <- as.data.frame(static_data)
    # static_df[, -1] <- lapply(static_df[, -1], as.numeric)
    # df_percentages <- static_df %>%
    #   mutate(across(-`Insect Role`, ~ . / rowSums(static_df[, -1]))) %>%
    #   mutate(`Insect Role` = factor(`Insect Role`, levels = rev(`Insect Role`)))
    # df_plot <- melt(df_percentages, id.vars = "Insect Role", variable.name = "Bird Type",
    #                 value.name = "proportion")
    # save_bar_chart_insect_role(df_plot, name)
    # 
    # # # Save the node and edge data for Gephi
    # # nodes_o <- paste0("nodes_", gsub(" ", "_", name), ".csv")
    # # edges_o <- paste0("edges_", gsub(" ", "_", name), ".csv")
    # # convert_data_for_gephi(data_merged, bird_forest_coverage, nodes_o, edges_o)
    # 
    # # Visualize the extinction curve in different site type
    # repeat_time <- 10000
    # results_list[["random"]][[name]] <- extinction_curve(bird_insect_web, insect_role, repeat_time, )
    # abundance <- bird_forest_coverage %>% arrange(Bird.Rarity.Score) %>% pull(BirdLife.Name.CORRECT)
    # results_list[["abundance"]][[name]] <- extinction_curve(bird_insect_web, insect_role, repeat_time, abundance)
    # forest_dependence_distance <- bird_forest_coverage %>% arrange(average.Distance.to.Forest.km) %>% pull(BirdLife.Name.CORRECT)
    # results_list[["distance_to_forest"]][[name]] <- extinction_curve(bird_insect_web, insect_role, repeat_time, forest_dependence_distance)
    # forest_dependence_canopy <- bird_forest_coverage %>% arrange(desc(Bird.Forest.Cover.200m)) %>% pull(BirdLife.Name.CORRECT)
    # results_list[["forest_canopy"]][[name]] <- extinction_curve(bird_insect_web, insect_role, repeat_time, forest_dependence_canopy)
  }

  # for(metric in names(results_list)){
  #   stability_pest_cotrol_result <- list(c(metric))
  #   save_extinction_curve_by_metric(results_list[[metric]], gsub("_", " ", metric))
  #   for(sites in names(results_list[[metric]])){
  #     normalize <- function(x) {
  #       return ((x - min(x)) / (max(x) - min(x)))
  #     }
  #     normalized_sequence <- normalize(results_list[[metric]][[sites]]$mean)
  #     stability_pest_cotrol_sum <- calculate_auc(normalized_sequence)
  #     stability_pest_cotrol_result[[length(stability_pest_cotrol_result) + 1]] <- stability_pest_cotrol_sum
  #     half_life_proportion <- which(normalized_sequence < 0.5)[1] / length(normalized_sequence) # half-life stability
  #     stability_pest_cotrol_result[[length(stability_pest_cotrol_result) + 1]] <- half_life_proportion*100
  #   }
  #   # results_pest_control_df <- rbind(results_pest_control_df, setNames(data.frame(t(unlist(stability_pest_cotrol_result)), stringsAsFactors = FALSE), c("metrics", "total_avg", "total_half", "forest_avg", "forest_half", "edge_avg", "edge_half", "agriculture_avg", "agriculture_half")), stringsAsFactors = FALSE)
  #   results_pest_control_df <- rbind(results_pest_control_df, setNames(data.frame(t(unlist(stability_pest_cotrol_result)), stringsAsFactors = FALSE), c("metrics", "forest_avg", "forest_half", "edge_avg", "edge_half", "agriculture_avg", "agriculture_half")), stringsAsFactors = FALSE)
  # }
# 
#   pest_control_name <- ifelse(is_all, {is_all <<- FALSE; "stability_potential_pest_control"}, "stability_confirmed_pest_control")
#   write.csv(results_pest_control_df, file = paste0(pest_control_name, ".csv"), row.names = FALSE)
# 
#   forest_canopy_row <- results_pest_control_df %>%
#     filter(metrics == "forest_canopy")
# 
#   # 将相应的数据更新到 basic_network_analysis 数据框中
#   basic_network_analysis <- basic_network_analysis %>%
#     mutate(
#       auc = case_when(
#         network_name == "Forest Interior" ~ forest_canopy_row$forest_avg,
#         network_name == "Forest Edge" ~ forest_canopy_row$edge_avg,
#         network_name == "Agriculture" ~ forest_canopy_row$agriculture_avg
#       ),
#       half_life = case_when(
#         network_name == "Forest Interior" ~ forest_canopy_row$forest_half,
#         network_name == "Forest Edge" ~ forest_canopy_row$edge_half,
#         network_name == "Agriculture" ~ forest_canopy_row$agriculture_half
#       )
#     )
#   
#   
#   results_pest_control_long <- results_pest_control_df %>%
#     pivot_longer(cols = -metrics, names_to = "variable", values_to = "value") %>%
#     mutate(value = as.numeric(value))
# 
#   for (metric in unique(results_pest_control_df$metrics)) {
#     p <- plot_auc_half_life(results_pest_control_long, metric)
#     ggsave(paste0("plot_", metric, ".png"), plot = p, width = 8, height = 6)
#   }
#   
#   plot_network_metric(basic_network_analysis)
#   plot_species_metric(species_result)
#   
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

results <- statistic_abundance(data_individual_in_site_type, bird_meta)

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


# input targeted list of analysis and web
venn_calculate <- function(sequence, bird_insect_web) {
  pest_species <- insect_role_all$Species[insect_role_all$Role == "Pest"]
  pests_in_df <- intersect(rownames(bird_insect_web), pest_species)
  num_all_insect <- length(pests_in_df)
  # Filter the bird_insect_web based on the sequence provided
  filtered_data <- bird_insect_web[, sequence, drop = FALSE]
  filtered_data <- filtered_data[apply(filtered_data, 1, function(row) any(row != 0)), ]
  
  # Number of birds (columns) and insects (rows)
  num_birds <- ncol(filtered_data)
  num_insects <- nrow(filtered_data)
  
  # Identify pests in the filtered data
  pest_species <- insect_role_all$Species[insect_role_all$Role == "Pest"]
  pests_in_df <- intersect(rownames(filtered_data), pest_species)
  num_pests <- length(pests_in_df)
  
  # Count non-forest birds
  non_forest_bird_columns <- colnames(filtered_data) %in% bird_forest_coverage$BirdLife.Name.CORRECT[bird_forest_coverage$Bird.Type == "Non-Forest bird"]
  num_non_forest_birds <- sum(non_forest_bird_columns)
  
  # Proportion of pests
  pest_proportion <- num_pests / num_insects
  
  # Print the results
  cat("Pest:", num_pests, "Non-pest:", num_insects - num_pests, "Proportion:",  num_pests / num_all_insect,"Proportion:",  num_insects / num_all_insect, "\n")
  cat("Insect", num_insects, "Bird:", num_birds ,"\n")
  cat("Forest Bird:", num_birds - num_non_forest_birds, "Non-forest Bird:", num_non_forest_birds ,"\n")
  cat("\n")
}


# Initialize an empty list to store the processed networks
data_venn <- list()
for(name in names(data_in_site_type)){
  print(name)
  data_merged <- data_in_site_type[[name]]
  
  # Filter valid columns based on bird names in bird_forest_coverage
  valid_columns <- bird_forest_coverage$BirdLife.Name.CORRECT %in% colnames(data_merged)
  data_merged_reordered <- data_merged[, bird_forest_coverage$BirdLife.Name.CORRECT[valid_columns]]
  
  # Add insect species as a column and convert to row names
  data_merged_reordered <- cbind(Insect_species=data_merged[["Insect_species"]], data_merged_reordered)
  bird_insect_web <- data_merged_reordered %>% column_to_rownames(var = 'Insect_species')
  
  # Convert to binary (presence/absence) matrix
  bird_insect_web[bird_insect_web != 0] <- 1
  
  # Remove rows that are entirely zero
  bird_insect_web <- bird_insect_web[rowSums(bird_insect_web != 0) > 0, ]
  
  # Remove columns that are entirely zero
  bird_insect_web <- bird_insect_web[, colSums(bird_insect_web != 0) > 0]
  
  # Save the processed network into the data_venn list
  data_venn[[name]] <- bird_insect_web
  
  insect_species_in_web <- rownames(bird_insect_web)
  
  # Extract pest species from insect_role dataframe
  pest_species <- insect_role$Species[insect_role$Role == "Pest"]
  
  # Determine which of the insect species in the bird_insect_web are pests
  pests_in_web <- intersect(insect_species_in_web, pest_species)
  
  # Calculate the number of pests and non-pests
  num_pests <- length(pests_in_web)
  num_non_pests <- length(insect_species_in_web) - num_pests
  
  # Output the results
  cat("Number of pests:", num_pests, "\n")
  cat("Number of non-pests:", num_non_pests, "\n")
}


# Assuming data_venn contains the three data frames named "A", "B", "C"
A <- colnames(data_venn[["Forest Interior"]])
B <- colnames(data_venn[["Forest Edge"]])
C <- colnames(data_venn[["Agriculture"]])

# Calculate the set operations
A_only <- setdiff(A, union(B, C))
B_only <- setdiff(B, union(A, C))
C_only <- setdiff(C, union(A, B))

A_and_B <- intersect(A, B)
A_and_C <- intersect(A, C)
B_and_C <- intersect(B, C)

A_B_C <- Reduce(intersect, list(A, B, C))

for(name in names(data_venn)) {
  cat("Processing web:", name, "\n")
  web <- data_venn[[name]]
  venn_calculate(A_B_C, web)
}

venn_calculate(A_only, data_venn[["Forest Interior"]])
union(B_only, setdiff(A_and_B, A_B_C))
venn_calculate(, data_venn[["Forest Interior"]])

# start simulate the forest degrate
venn_calculate(A, data_venn[["Forest Interior"]])
venn_calculate(setdiff(A_and_B, A_B_C), data_venn[["Forest Interior"]])
venn_calculate(setdiff(A_and_C, A_B_C), data_venn[["Forest Interior"]])
venn_calculate(A_B_C, data_venn[["Forest Interior"]])
venn_calculate(A_only, data_venn[["Forest Interior"]])


venn_calculate(B, data_venn[["Forest Edge"]])
venn_calculate(setdiff(A_and_B, A_B_C), data_venn[["Forest Edge"]])
# venn_calculate(setdiff(A_and_C, A_B_C), data_venn[["Forest Edge"]])
venn_calculate(A_B_C, data_venn[["Forest Edge"]])

venn_calculate(setdiff(B_and_C, A_B_C), data_venn[["Forest Edge"]])
venn_calculate(setdiff(B_only, A_B_C), data_venn[["Forest Edge"]])

venn_calculate(C, data_venn[["Agriculture"]])
venn_calculate(A_B_C, data_venn[["Agriculture"]])
venn_calculate(setdiff(B_and_C, A_B_C), data_venn[["Agriculture"]])

venn_calculate(setdiff(A_and_C, A_B_C), data_venn[["Agriculture"]])
venn_calculate(C_only, data_venn[["Agriculture"]])

# start simulate the forest degrate
venn_calculate(A, data_venn[["Forest Interior"]]) # a forest
venn_calculate(A_and_B, data_venn[["Forest Edge"]])# a b edge
venn_calculate(setdiff(B, A_and_B), data_venn[["Forest Edge"]])# b only b and c - abc edge
venn_calculate(B, data_venn[["Forest Edge"]])# b  edge
venn_calculate(B_and_C, data_venn[["Agriculture"]])# b c agriculture
venn_calculate(setdiff(C, B_and_C), data_venn[["Agriculture"]])# c only a and c - abc agriculture
venn_calculate(C, data_venn[["Agriculture"]])# c agriculture

# start simulate the forest degrate 2#
venn_calculate(A, data_venn[["Forest Interior"]]) # a forest
venn_calculate(A_and_B, data_venn[["Forest Edge"]])# a b edge
venn_calculate(A_B_C, data_venn[["Agriculture"]]) # a forest
venn_calculate(setdiff(C, A_B_C), data_venn[["Agriculture"]]) # a forest
venn_calculate(C, data_venn[["Agriculture"]])# c agriculture

# start simulate the forest degrate 2#
venn_calculate(A, data_venn[["Forest Interior"]]) # a forest
venn_calculate(A_and_C, data_venn[["Agriculture"]])# a b edge
venn_calculate(setdiff(C, A_and_B), data_venn[["Agriculture"]]) # a forest
venn_calculate(C, data_venn[["Agriculture"]])# c agriculture

venn_calculate(A_and_C, data_venn[["Agriculture"]])# a b edge
venn_calculate(A_and_C, data_venn[["Forest Interior"]])# a b edge

combined_union <- union(agri_edge_minus_common, agri_fore_minus_common)
combined_union <- union(combined_union, common_birds_3_intersect)
agri_only <- setdiff(bird_3_union$`Agriculture`, agri_fore_minus_common)



venn_calculate(A_B_C, data_venn[["Forest Interior"]])
venn_calculate(A_B_C, data_venn[["Forest Edge"]])
venn_calculate(A_B_C, data_venn[["Agriculture"]])

venn_calculate(A_only, data_venn[["Forest Interior"]])
venn_calculate(B_only, data_venn[["Forest Edge"]])
venn_calculate(C_only, data_venn[["Agriculture"]])

venn_calculate(setdiff(A_and_B, A_B_C), data_venn[["Forest Interior"]])
venn_calculate(setdiff(A_and_C, A_B_C), data_venn[["Forest Interior"]])

venn_calculate(setdiff(A_and_B, A_B_C), data_venn[["Forest Edge"]])
venn_calculate(setdiff(B_and_C, A_B_C), data_venn[["Forest Edge"]])

venn_calculate(setdiff(B_and_C, A_B_C), data_venn[["Agriculture"]])
venn_calculate(setdiff(A_and_C, A_B_C), data_venn[["Agriculture"]])


venn_calculate(A_and_B, data_venn[["Forest Interior"]])
venn_calculate(A_and_C, data_venn[["Forest Interior"]])

venn_calculate(A_and_B, data_venn[["Forest Edge"]])
venn_calculate(B_and_C, data_venn[["Forest Edge"]])

venn_calculate(B_and_C, data_venn[["Agriculture"]])
venn_calculate(A_and_C, data_venn[["Agriculture"]])



venn_calculate(union(setdiff(B_and_C, A_B_C), C_only), data_venn[["Agriculture"]])
venn_calculate(union(B_and_C, C_only), data_venn[["Agriculture"]])


