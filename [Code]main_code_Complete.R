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
library(igraph)


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


calculate_modularity <- function(bird_insect_web){
  # 将 bird_insect_web 转换为 igraph 对象
  g <- graph_from_incidence_matrix(bird_insect_web, weighted = TRUE)
  # 使用适合二部图的社区检测算法
  # `cluster_fast_greedy` 不适用于加权二部图，所以这里使用 `cluster_louvain`
  community <- cluster_louvain(g)
  # 计算模块度
  modularity_value <- modularity(community)
  return(modularity_value)
  # # 如果想获取每个节点的社区
  # membership(community)
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
extinction_curve <- function(matrix_network, insect_roles, repeat_count = 10, extinction_sequence = NULL) {
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
  return(list(mean = mean_remaining_pests, error = sd_remaining_pests))
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

# # function to visualize the extinction curve of pest control
# save_extinction_curve <- function(extinction_data, plot_name) {
#   time_stamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
#   
#   # 合并所有灭绝曲线成一个数据框
#   combined_df <- bind_rows(lapply(extinction_data, function(data) {
#     data.frame(
#       Step = 0:(length(data$curve) - 1),
#       Remaining_Pests = data$curve,
#       Sequence = data$name
#     )
#   }))
#   
#   # 找到每个序列的第一个点和最后一个点
#   reference_lines <- combined_df %>%
#     group_by(Sequence) %>%
#     summarize(
#       First_Step = min(Step),
#       Last_Step = max(Step),
#       First_Remaining_Pests = first(Remaining_Pests),
#       Last_Remaining_Pests = last(Remaining_Pests)
#     )
#   
#   # 绘制结果
#   p <- ggplot(combined_df, aes(x = Step, y = Remaining_Pests, color = Sequence)) +
#     geom_line() +
#     geom_point() +
#     geom_segment(data = reference_lines, 
#                  aes(x = First_Step, xend = Last_Step, 
#                      y = First_Remaining_Pests, yend = Last_Remaining_Pests), 
#                  linetype = "dashed", color = "black") +
#     labs(title = paste0("Extinction Curves (Pest under control) in ", plot_name, " Site"),
#          x = "Extinction Step",
#          y = "Pests under control") +
#     theme_minimal() +
#     scale_color_manual(values = rainbow(length(extinction_data))) +
#     theme(legend.title = element_blank(),
#           legend.position = "right",
#           legend.box.margin = margin(0, 10, 0, 0))
#   
#   # 设置图形的宽度和高度，扩展图形的画布
#   ggsave(paste0("extinction_curves_", plot_name, "_", time_stamp, ".jpg"), p,
#          width = 7.5, height = 5, units = "in")
# }


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
    geom_ribbon(aes(ymin = Mean_Remaining_Pests - Error_Remaining_Pests, ymax = Mean_Remaining_Pests + Error_Remaining_Pests), alpha = 0.2, color = NA) +
    annotate("segment", x = 0, xend = 1, y = 1, yend = 0, color = "grey", linetype = "dashed") +
    labs(title = "Resilience of natural pest control under predator extinction scenarios",
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
# data_in_site_type <- list(
#   # "Total" = combine_data(bird_meta, data_RAA, data), 
#   "Forest Interior" = combine_data(bird_meta, data_RAA_forest, data), 
#   "Forest Edge" = combine_data(bird_meta, data_RAA_edge, data), 
#   "Agriculture" = combine_data(bird_meta, data_RAA_agriculture, data)
#   )
data_in_site_type <- list(
  "Forest Interior" = combine_data(bird_meta, data_RAA_forest_total, data), 
  "Agriculture" = combine_data(bird_meta, data_RAA_agriculture, data)
)
results_list <- list()


# Main loop for data in different site type
is_all <- TRUE
for(insect_role in list(insect_role_all, insect_role_confirmed)){
  results_pest_control_df <- data.frame(matrix(ncol = 9, nrow = 0), stringsAsFactors = FALSE) # make empty data frame to store stability pest control result
  names(results_pest_control_df) <- c("metrics", "total_avg", "total_half", "forest_avg", "forest_half", "edge_avg", "edge_half", "agriculture_avg", "agriculture_half")
  pest_distribution <- list()
  
  for(name in names(data_in_site_type)){
    print(name)
    data_merged <- data_in_site_type[[name]]
    valid_columns <- bird_forest_coverage$BirdLife.Name.CORRECT %in% colnames(data_merged)
    data_merged_reordered <- data_merged[, bird_forest_coverage$BirdLife.Name.CORRECT[valid_columns]]
    data_merged_reordered <- cbind(Insect_species=data_merged[["Insect_species"]], data_merged_reordered)
    bird_insect_web <- data_merged_reordered %>% column_to_rownames(var = 'Insect_species')
    
    # result_modularity <- calculate_modularity(bird_insect_web)
    # print(paste("modularity", result_modularity))
    # result_nested <- nested(bird_insect_web, method = "NODF")
    # print(paste("nestedness", result_nested))
    # species_specialization <- specieslevel(bird_insect_web, level = "higher", index = "d")$`d'`
    # result_max_specialization <- max(species_specialization, na.rm = TRUE)
    # print(paste("max_specialization", result_max_specialization))
    # result_mean_specialization <- mean(species_specialization, na.rm = TRUE)
    # print(paste("mean_specialization", result_mean_specialization))
    # result_median_specialization <- median(species_specialization, na.rm = TRUE)
    # print(paste("median_specialization", result_median_specialization))
    # result_min_specialization <- min(species_specialization, na.rm = TRUE)
    # print(paste("min_specialization", result_min_specialization))
    
    # # statistics pest distribution
    # pest_species <- insect_role$Species[insect_role$Role == "Pest"]
    # pest_distribution[[length(pest_distribution) + 1]] <- sum(pest_species %in% rownames(bird_insect_web))

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
    # static_data <- cbind(
    #   c("Pest", "Non-pest", "Natural enemy", "Vector of animal disease"),
    #   rbind(count_nonzero(bird_insect_web_by_role$Pest, forest_birds, non_forest_birds),
    #         count_nonzero(bird_insect_web_by_role$`Non-pest`, forest_birds, non_forest_birds),
    #         count_nonzero(bird_insect_web_by_role$`Natural enemy`, forest_birds, non_forest_birds),
    #         count_nonzero(bird_insect_web_by_role$`Vector of animal disease`, forest_birds, non_forest_birds)))
    # colnames(static_data) <- c("Insect Role", "Forest birds only", "Shared", "Non-forest birds only")
    # static_df <- as.data.frame(static_data)
    # static_df[, -1] <- lapply(static_df[, -1], as.numeric)
    # df_percentages <- static_df %>%
    #   mutate(across(-`Insect Role`, ~ . / rowSums(static_df[, -1]))) %>%
    #   mutate(`Insect Role` = factor(`Insect Role`, levels = rev(`Insect Role`)))
    # df_plot <- melt(df_percentages, id.vars = "Insect Role", variable.name = "Bird Type",
    #                 value.name = "proportion")
    # save_bar_chart_insect_role(df_plot, name)

    # # Save the node and edge data for Gephi
    # nodes_o <- paste0("nodes_", gsub(" ", "_", name), ".csv")
    # edges_o <- paste0("edges_", gsub(" ", "_", name), ".csv")
    # convert_data_for_gephi(data_merged, bird_forest_coverage, nodes_o, edges_o)
    
    # # Visualize the extinction curve in different metrics
    # result_1 <- extinction_curve(bird_insect_web, insect_role)
    # forest_dependence_frequency <- bird_forest_coverage %>% arrange(desc(ForestDependency)) %>% pull(BirdLife.Name.CORRECT)
    # result_2 <- extinction_curve(bird_insect_web, insect_role, forest_dependence_frequency)
    # forest_dependence_distance <- bird_forest_coverage %>% arrange(average.Distance.to.Forest.km) %>% pull(BirdLife.Name.CORRECT)
    # result_3 <- extinction_curve(bird_insect_web, insect_role, )
    # forest_dependence_canopy <- bird_forest_coverage %>% arrange(desc(Bird.Forest.Cover.200m)) %>% pull(BirdLife.Name.CORRECT)
    # result_4 <- extinction_curve(bird_insect_web, insect_role, )
    # 
    # extinction_data <- list(
    #   list(curve = result_1, name = "Average Random(1000 repeats)"),
    #   list(curve = result_2, name = "Forest Dependence(Frequency)"),
    #   list(curve = result_3, name = "Forest Dependence(Distance)"),
    #   list(curve = result_4, name = "Forest Dependence(Canopy)"),
    #   list(curve = result_4, name = "Forest Dependence(Canopy)")
    # )
    # save_extinction_curve(extinction_data, name)

    
    # Visualize the extinction curve in different site type
    repeat_time <- 100
    results_list[["random"]][[name]] <- extinction_curve(bird_insect_web, insect_role, repeat_time, )
    abundance <- bird_forest_coverage %>% arrange(Bird.Rarity.Score) %>% pull(BirdLife.Name.CORRECT)
    results_list[["abundance"]][[name]] <- extinction_curve(bird_insect_web, insect_role, repeat_time, abundance)
    forest_dependence_distance <- bird_forest_coverage %>% arrange(average.Distance.to.Forest.km) %>% pull(BirdLife.Name.CORRECT)
    results_list[["distance_to_forest"]][[name]] <- extinction_curve(bird_insect_web, insect_role, repeat_time, forest_dependence_distance)
    forest_dependence_canopy <- bird_forest_coverage %>% arrange(desc(Bird.Forest.Cover.200m)) %>% pull(BirdLife.Name.CORRECT)
    results_list[["forest_canopy"]][[name]] <- extinction_curve(bird_insect_web, insect_role, repeat_time, forest_dependence_canopy)
    
    # plot_extinction_curve(results_list[["random"]][[name]])
    # plot_extinction_curve(results_list[["abundance"]][[name]])
    # plot_extinction_curve(results_list[["distance_to_forest"]][[name]])
    # plot_extinction_curve(results_list[["forest_canopy"]][[name]])
    
    # # Find the Correlation between the order and the bird body metrics
    # abundance <- bird_forest_coverage %>% arrange(Bird.Rarity.Score) %>% pull(BirdLife.Name.CORRECT)
    # for(bird_body_metric in c("Beak.Length.Culmen", "Beak.Length.Nares", "Beak.Width", "Beak.Depth", "Tarsus.Length", "Wing.Length", "Kipps.Distance", "Secondary1", "Hand.wing.Index", "Tail.Length", "Gape.Size", "Weight...g.")){
    #   metric_sequence <- extract_beak_length(abundance, bird_meta, "BirdLife.Name.CORRECT", bird_body_metric)
    #   print(bird_body_metric)
    #   check_correlation(metric_sequence)
    # }
    # forest_dependence_distance <- bird_forest_coverage %>% arrange(average.Distance.to.Forest.km) %>% pull(BirdLife.Name.CORRECT)
    # for(bird_body_metric in c("Beak.Length.Culmen", "Beak.Length.Nares", "Beak.Width", "Beak.Depth", "Tarsus.Length", "Wing.Length", "Kipps.Distance", "Secondary1", "Hand.wing.Index", "Tail.Length", "Gape.Size", "Weight...g.")){
    #   metric_sequence <- extract_beak_length(forest_dependence_distance, bird_meta, "BirdLife.Name.CORRECT", bird_body_metric)
    #   print(bird_body_metric)
    #   check_correlation(metric_sequence)
    # }
    # forest_dependence_canopy <- bird_forest_coverage %>% arrange(desc(Bird.Forest.Cover.200m)) %>% pull(BirdLife.Name.CORRECT)
    # for(bird_body_metric in c("Beak.Length.Culmen", "Beak.Length.Nares", "Beak.Width", "Beak.Depth", "Tarsus.Length", "Wing.Length", "Kipps.Distance", "Secondary1", "Hand.wing.Index", "Tail.Length", "Gape.Size", "Weight...g.")){
    #   metric_sequence <- extract_beak_length(forest_dependence_canopy, bird_meta, "BirdLife.Name.CORRECT", bird_body_metric)
    #   print(bird_body_metric)
    #   check_correlation(metric_sequence)
    # }
  }

  for(metric in names(results_list)){
    stability_pest_cotrol_result <- list(c(metric))
    save_extinction_curve_by_metric(results_list[[metric]], gsub("_", " ", metric))
    # for(sites in names(results_list[[metric]])){
    #   normalize <- function(x) {
    #     return ((x - min(x)) / (max(x) - min(x)))
    #   }
    #   normalized_sequence <- normalize(results_list[[metric]][[sites]])
    #   stability_pest_cotrol_sum <- sum(normalized_sequence)/length(results_list[[metric]][[sites]])
    #   stability_pest_cotrol_result[[length(stability_pest_cotrol_result) + 1]] <- stability_pest_cotrol_sum
    #   half_life_proportion <- which(normalized_sequence < 0.5)[1] / length(normalized_sequence) # half-life stability
    #   stability_pest_cotrol_result[[length(stability_pest_cotrol_result) + 1]] <- half_life_proportion*100
    # }
    # results_pest_control_df <- rbind(results_pest_control_df, setNames(data.frame(t(unlist(stability_pest_cotrol_result)), stringsAsFactors = FALSE), c("metrics", "total_avg", "total_half", "forest_avg", "forest_half", "edge_avg", "edge_half", "agriculture_avg", "agriculture_half")), stringsAsFactors = FALSE)
  }
  # pest_control_name <- ifelse(is_all, {is_all <<- FALSE; "stability_potential_pest_control"}, "stability_confirmed_pest_control")
  # write.csv(results_pest_control_df, file = paste0(pest_control_name, ".csv"), row.names = FALSE)
  
  # print(pest_distribution) # actually same in different sites confirmed 40 potential 158
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


























