# 得到三交
bird_3_intersect <- list()
bird_2_intersect_agri_edge <- list()
bird_2_intersect_agri_fore <- list()
bird_2_intersect_edge_fore <- list()

for(name in names(data_in_site_type)){
  print(name)
  data_merged <- data_in_site_type[[name]]
  valid_columns <- bird_forest_coverage$BirdLife.Name.CORRECT %in% colnames(data_merged)
  data_merged_reordered <- data_merged[, bird_forest_coverage$BirdLife.Name.CORRECT[valid_columns]]
  data_merged_reordered <- cbind(Insect_species=data_merged[["Insect_species"]], data_merged_reordered)
  bird_insect_web_weight <- data_merged_reordered %>% column_to_rownames(var = 'Insect_species')
  bird_insect_web <- bird_insect_web_weight
  bird_insect_web[bird_insect_web != 0] <- 1
  
  # Append the column names as a list element named by `name`
  bird_3_intersect <- append(bird_3_intersect, setNames(list(colnames(bird_insect_web)), name))
  
  # Handle pairwise intersections
  if(name == "Forest Interior"){
    bird_2_intersect_agri_fore <- append(bird_2_intersect_agri_fore, setNames(list(colnames(bird_insect_web)), "Forest Interior"))
    bird_2_intersect_edge_fore <- append(bird_2_intersect_edge_fore, setNames(list(colnames(bird_insect_web)), "Forest Interior"))
  } else if(name == "Forest Edge") {
    bird_2_intersect_agri_edge <- append(bird_2_intersect_agri_edge, setNames(list(colnames(bird_insect_web)), "Forest Edge"))
    bird_2_intersect_edge_fore <- append(bird_2_intersect_edge_fore, setNames(list(colnames(bird_insect_web)), "Forest Edge"))
  } else { # Assuming the only other option is "Agriculture"
    bird_2_intersect_agri_edge <- append(bird_2_intersect_agri_edge, setNames(list(colnames(bird_insect_web)), "Agriculture"))
    bird_2_intersect_agri_fore <- append(bird_2_intersect_agri_fore, setNames(list(colnames(bird_insect_web)), "Agriculture"))
  }
}

# Calculate the intersection for all three networks
common_birds_3_intersect <- Reduce(intersect, bird_3_intersect)

# Calculate pairwise intersections
common_birds_2_intersect_agri_edge <- Reduce(intersect, bird_2_intersect_agri_edge)
common_birds_2_intersect_agri_fore <- Reduce(intersect, bird_2_intersect_agri_fore)
common_birds_2_intersect_edge_fore <- Reduce(intersect, bird_2_intersect_edge_fore)

# Subtract the intersection of all three networks from the pairwise intersections
agri_edge_minus_common <- setdiff(common_birds_2_intersect_agri_edge, common_birds_3_intersect)
agri_fore_minus_common <- setdiff(common_birds_2_intersect_agri_fore, common_birds_3_intersect)
edge_fore_minus_common <- setdiff(common_birds_2_intersect_edge_fore, common_birds_3_intersect)



# 筛选三交
insect_control_by_bird_3_intersect <- list()
for(name in names(data_in_site_type)){
  print(name)
  data_merged <- data_in_site_type[[name]]
  valid_columns <- bird_forest_coverage$BirdLife.Name.CORRECT %in% colnames(data_merged)
  data_merged_reordered <- data_merged[, bird_forest_coverage$BirdLife.Name.CORRECT[valid_columns]]
  data_merged_reordered <- cbind(Insect_species=data_merged[["Insect_species"]], data_merged_reordered)
  bird_insect_web_weight <- data_merged_reordered %>% column_to_rownames(var = 'Insect_species')
  bird_insect_web <- bird_insect_web_weight
  bird_insect_web[bird_insect_web != 0] <- 1
  
  filtered_data <- bird_insect_web[, common_birds, drop = FALSE]
  filtered_data <- filtered_data[apply(filtered_data, 1, function(row) any(row != 0)), ]
  insect_control_by_bird_3_intersect <- append(insect_control_by_bird_3_intersect, setNames(list(rownames(bird_insect_web)), name))
}
all_identical <- identical(insect_control_by_bird_3_intersect[[1]], insect_control_by_bird_3_intersect[[2]]) &&
  identical(insect_control_by_bird_3_intersect[[1]], insect_control_by_bird_3_intersect[[3]])

cat("Are all elements in the list identical (except for their names)? ", all_identical, "\n")


insect_control_by_bird_2_intersect_agri_edge <- list()
for(name in names(data_in_site_type)){
  print(name)
  data_merged <- data_in_site_type[[name]]
  valid_columns <- bird_forest_coverage$BirdLife.Name.CORRECT %in% colnames(data_merged)
  data_merged_reordered <- data_merged[, bird_forest_coverage$BirdLife.Name.CORRECT[valid_columns]]
  data_merged_reordered <- cbind(Insect_species=data_merged[["Insect_species"]], data_merged_reordered)
  bird_insect_web_weight <- data_merged_reordered %>% column_to_rownames(var = 'Insect_species')
  bird_insect_web <- bird_insect_web_weight
  bird_insect_web[bird_insect_web != 0] <- 1
  
  if (name == "Agriculture" | name == "Forest Edge"){
    filtered_data <- bird_insect_web[, agri_edge_minus_common, drop = FALSE]
    filtered_data <- filtered_data[apply(filtered_data, 1, function(row) any(row != 0)), ]
    insect_control_by_bird_2_intersect_agri_edge <- append(insect_control_by_bird_2_intersect_agri_edge, setNames(list(rownames(bird_insect_web)), name))
  }
}

all_identical <- identical(insect_control_by_bird_2_intersect_agri_edge[[1]], insect_control_by_bird_2_intersect_agri_edge[[2]])

cat("Are all elements in the list identical (except for their names)? ", all_identical, "\n")


insect_control_by_bird_2_intersect_agri_fore <- list()
for(name in names(data_in_site_type)){
  print(name)
  data_merged <- data_in_site_type[[name]]
  valid_columns <- bird_forest_coverage$BirdLife.Name.CORRECT %in% colnames(data_merged)
  data_merged_reordered <- data_merged[, bird_forest_coverage$BirdLife.Name.CORRECT[valid_columns]]
  data_merged_reordered <- cbind(Insect_species=data_merged[["Insect_species"]], data_merged_reordered)
  bird_insect_web_weight <- data_merged_reordered %>% column_to_rownames(var = 'Insect_species')
  bird_insect_web <- bird_insect_web_weight
  bird_insect_web[bird_insect_web != 0] <- 1
  
  if (name == "Agriculture" | name == "Forest Interior"){
    filtered_data <- bird_insect_web[, agri_fore_minus_common, drop = FALSE]
    filtered_data <- filtered_data[apply(filtered_data, 1, function(row) any(row != 0)), ]
    insect_control_by_bird_2_intersect_agri_fore <- append(insect_control_by_bird_2_intersect_agri_fore, setNames(list(rownames(bird_insect_web)), name))
  }
}

all_identical <- identical(insect_control_by_bird_2_intersect_agri_fore[[1]], insect_control_by_bird_2_intersect_agri_fore[[2]])

cat("Are all elements in the list identical (except for their names)? ", all_identical, "\n")


insect_control_by_bird_2_intersect_edge_fore <- list()
for(name in names(data_in_site_type)){
  print(name)
  data_merged <- data_in_site_type[[name]]
  valid_columns <- bird_forest_coverage$BirdLife.Name.CORRECT %in% colnames(data_merged)
  data_merged_reordered <- data_merged[, bird_forest_coverage$BirdLife.Name.CORRECT[valid_columns]]
  data_merged_reordered <- cbind(Insect_species=data_merged[["Insect_species"]], data_merged_reordered)
  bird_insect_web_weight <- data_merged_reordered %>% column_to_rownames(var = 'Insect_species')
  bird_insect_web <- bird_insect_web_weight
  bird_insect_web[bird_insect_web != 0] <- 1
  
  if (name == "Forest Interior" | name == "Forest Interior"){
    filtered_data <- bird_insect_web[, edge_fore_minus_common, drop = FALSE]
    filtered_data <- filtered_data[apply(filtered_data, 1, function(row) any(row != 0)), ]
    insect_control_by_bird_2_intersect_edge_fore <- append(insect_control_by_bird_2_intersect_edge_fore, setNames(list(rownames(bird_insect_web)), name))
  }
}

all_identical <- identical(insect_control_by_bird_2_intersect_edge_fore[[1]], insect_control_by_bird_2_intersect_edge_fore[[2]])

cat("Are all elements in the list identical (except for their names)? ", all_identical, "\n")
