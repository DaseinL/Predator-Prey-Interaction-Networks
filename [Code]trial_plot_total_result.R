# 使用指定的配色方案
colors <- c(
  "modularity" = "#E4572E",    # 深橙色
  "nestedness" = "#F3A712",    # 柔和的橙黄色
  "specialization" = "#A8C686", # 亮黄色
  "diameter" = "#76B041",       # 柔和的蓝色
  "auc" = "#4B8BBE",           # 深蓝色
  "half_life" = "#30638E"      # 柔和的绿色
)


# # 使用指定的配色方案
# colors <- c(
#   "modularity" = "#C02C27",    # 深橙色
#   "nestedness" = "#FF6F6B",    # 柔和的橙黄色
#   "specialization" = "#F4A183", # 亮黄色
#   "diameter" = "#8073B2",       # 柔和的蓝色
#   "auc" = "#8BBAD1",           # 深蓝色
#   "half_life" = "#2166AC"      # 柔和的绿色
# )
# 
# # 
# # 使用指定的配色方案
# colors <- c(
#   "modularity" = "#D55E00",    # 深橙色
#   "nestedness" = "#E69F00",    # 柔和的橙黄色
#   "specialization" = "#F0E442", # 亮黄色
#   "diameter" = "#009E73",       # 柔和的蓝色
#   "auc" = "#56B4E9",           # 深蓝色
#   "half_life" = "#0072B2"      # 柔和的绿色
# )



# 加载所需包
library(ggplot2)
library(dplyr)
library(tidyr)

# # 假设 basic_network_analysis 是你的数据框
# basic_network_analysis <- data.frame(
#   network_name = c("Forest Interior", "Forest Edge", "Agriculture"),
#   modularity = c(0.667422405487886, 0.727493182718954, 0.567355525293877),
#   nestedness = c(12.6195036509549, 14.7428720013782, 10.9018143168708),
#   specialization = c(0.72373348158306, 0.722501846401411, 0.592048196098917),
#   diameter = c(0.0882352941176471, 0.142857142857143, 0.0681818181818182),
#   auc = c(0.673850908119658, 0.688870222007722, 0.663559632378932),
#   half_life = c(78.3783783783784, 81.578947368421, 76.5957446808511)
# )

# 将数据从宽格式转换为长格式
df <- basic_network_analysis %>%
  gather(key = "Metric", value = "Value", -network_name) %>%
  mutate(LandUseType = case_when(
    network_name == "Forest Interior" ~ "Forest Interior",
    network_name == "Forest Edge" ~ "Forest Edge",
    network_name == "Agriculture" ~ "Agriculture"
  ))

# 将 LandUseType 设置为因子，并指定顺序为 森林、边缘、农田
df$LandUseType <- factor(df$LandUseType, levels = c("Forest Interior", "Forest Edge", "Agriculture"))

# 将 Metric 设置为因子，并指定小项的显示顺序
df$Metric <- factor(df$Metric, levels = c("modularity", "nestedness", "specialization", "diameter", "auc", "half_life"))

df$Value <- as.numeric(df$Value)
# 数据标准化：按每个指标在不同用地类型内的值进行标准化
df <- df %>%
  group_by(Metric) %>%
  mutate(StandardizedValue = (Value - mean(Value)) / sd(Value)) %>%
  ungroup()

# 绘制条形图，并在 y=0 处添加辅助线
p <- ggplot(df, aes(x = LandUseType, y = StandardizedValue, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +  # 条形图没有边框
  geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 0.7) +  # 添加 y=0 的辅助线
  scale_fill_manual(values = colors) +  # 使用指定的颜色
  theme_minimal() +
  labs(title = "Changes in Metrics Across Different Land Use Types",
       x = "Land Use Type",
       y = "Standardized Value",
       fill = "Metrics") +
  theme(
    axis.line = element_line(color = "black"),  # 添加x和y轴线
    panel.grid = element_blank(),  # 移除网格线
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.ticks = element_line(color = "black")  # 添加轴刻度
  )

# 显示图表
print(p)
