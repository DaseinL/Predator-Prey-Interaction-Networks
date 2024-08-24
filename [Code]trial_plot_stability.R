library(ggplot2)
library(tidyr)
library(dplyr)

# 数据框结构（confirmed pest）
# data <- data.frame(
#   metrics = c("random", "abundance", "distance_to_forest", "forest_canopy"),
#   total_avg = c(0.617361111, 0.604513889, 0.663888889, 0.653472222),
#   total_half = c(79.16666667, 75, 79.16666667, 75),
#   forest_avg = c(0.598648649, 0.737837838, 0.501351351, 0.460810811),
#   forest_half = c(70.27027027, 81.08108108, 48.64864865, 40.54054054),
#   edge_avg = c(0.695175439, 0.509868421, 0.589912281, 0.60745614),
#   edge_half = c(84.21052632, 60.52631579, 63.15789474, 60.52631579),
#   agriculture_avg = c(0.575256107, 0.501182033, 0.657210402, 0.654846336),
#   agriculture_half = c(65.95744681, 57.44680851, 70.21276596, 70.21276596)
# )

# all pest
data <- data.frame(
  metrics = c("random", "abundance", "distance_to_forest", "forest_canopy"),
  total_avg = c(0.617352321, 0.68204993, 0.619462025, 0.607419128),
  total_half = c(76.38888889, 79.16666667, 72.22222222, 73.61111111),
  forest_avg = c(0.732303732, 0.715894466, 0.582046332, 0.559202059),
  forest_half = c(89.18918919, 81.08108108, 64.86486486, 64.86486486),
  edge_avg = c(0.685016112, 0.638829216, 0.610902256, 0.602577873),
  edge_half = c(84.21052632, 71.05263158, 65.78947368, 65.78947368),
  agriculture_avg = c(0.601364914, 0.622641509, 0.640305098, 0.638297872),
  agriculture_half = c(65.95744681, 76.59574468, 70.21276596, 70.21276596),
  stringsAsFactors = FALSE
)

# 转换数据为长格式
data_long <- data %>%
  pivot_longer(cols = -metrics, names_to = "variable", values_to = "value")

# 绘图函数
plot_metric <- function(data, metric_name) {
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
    scale_y_continuous(name = "Value (x100 for average)", sec.axis = sec_axis(~ ., name = "Half-life (%)")) +
    scale_fill_manual(values = custom_colors) +
    scale_color_manual(values = c("Half-life (%)" = "blue")) +
    labs(title = paste("Metrics:", metric_name), x = "Variable", y = "Value") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom",
          panel.background = element_rect(fill = "white", color = NA),  # 设置面板背景为白色
          plot.background = element_rect(fill = "white", color = NA))   # 设置绘图背景为白色
  
  return(p)
}

# 绘制并保存图表
for (metric in unique(data$metrics)) {
  p <- plot_metric(data_long, metric)
  ggsave(paste0("plot_", metric, ".png"), plot = p, width = 8, height = 6)
}