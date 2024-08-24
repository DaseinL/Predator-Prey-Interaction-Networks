library(ggplot2)
library(dplyr)
library(tidyr)

# 定义数据
data <- data.frame(
  Site = c("Total", "Forest Interior", "Forest Edge", "Agriculture"),
  all_pest_modularity = c(0.640051477933599, 0.667435574771394, 0.726440864928892, 0.567560517252051),
  all_pest_nestedness = c(10.1818004374983, 3.6071500929944, 3.58264704253433, 3.31449449887201),
  confirmed_pest_modularity = c(0.637124127250994, 0.666495917913814, 0.72328298542654, 0.567572384044992),
  confirmed_pest_nestedness = c(10.1818004374983, 3.6071500929944, 3.58264704253433, 3.31449449887201)
)

data <- data.frame(
  Site = c("Forest Interior", "Forest Edge", "Agriculture"),
  all_pest_modularity = c(0.667435574771394, 0.726440864928892, 0.567560517252051),
  all_pest_nestedness = c(3.6071500929944, 3.58264704253433, 3.31449449887201),
  confirmed_pest_modularity = c(0.666495917913814, 0.72328298542654, 0.567572384044992),
  confirmed_pest_nestedness = c( 3.6071500929944, 3.58264704253433, 3.31449449887201)
)

# 转换数据为长格式
data_long <- data %>%
  pivot_longer(cols = -Site, names_to = "Metric", values_to = "Value") %>%
  separate(Metric, into = c("PestType", "Measure"), sep = "_", extra = "merge") %>%
  mutate(Metric = paste(PestType, Measure, sep = "_"))

# 绘制折线图
ggplot(data_long, aes(x = Site, y = Value, group = Metric, color = Metric)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  labs(title = "Modularity and Nestedness Across Sites",
       x = "Site",
       y = "Value",
       color = "Metric") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# function to draw an extinction curve
extinction_curve <- function(matrix_network, insect_roles, extinction_sequence = NULL) {
  matrix_network <- matrix_network[rowSums(matrix_network) != 0, colSums(matrix_network) != 0, drop = FALSE]
  
  if (is.null(extinction_sequence)) {
    extinction_sequence <- sample(colnames(matrix_network))
    repeat_count <- 10
  }else{
    repeat_count <- 1
  }
  extinction_sequence <- extinction_sequence[extinction_sequence %in% colnames(matrix_network)]
  remaining_pests_list <- rep(0, length(extinction_sequence) + 1)
  
  for (t in 1:repeat_count){
    print(paste0(t/repeat_count*100, "%"))
    i <- 1
    matrix_network_copy <- matrix_network
    matrix_network_copy <- matrix_network_copy[rowSums(matrix_network_copy) != 0, drop = FALSE]
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

