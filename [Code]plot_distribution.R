library(dplyr)
library(ggplot2)
library(tidyr)

# 过滤数据，只保留 phylum 为 "Arthropoda" 的数据
arthropoda_data <- data %>%
  filter(phylum == "Arthropoda")

# 处理 species 列，将 "sp." 视为空
arthropoda_data <- arthropoda_data %>%
  mutate(species = ifelse(str_detect(species, "\\bsp\\.$"), NA, species))

# 统计每一列的非空数据
non_empty_counts <- arthropoda_data %>%
  summarise(
    class_count = sum(!is.na(class) & class != ""),
    order_count = sum(!is.na(order) & order != ""),
    family_count = sum(!is.na(family) & family != ""),
    genus_count = sum(!is.na(genus) & genus != ""),
    species_count = sum(!is.na(species) & species != "")
  )

# 将统计结果转换为适合绘图的格式
counts_long <- non_empty_counts %>%
  pivot_longer(cols = everything(), names_to = "Category", values_to = "Count")

# 调整 Category 顺序，以确保图中的顺序从上至下为 纲、目、科、属、种
counts_long$Category <- factor(counts_long$Category, levels = rev(c("class_count", "order_count", "family_count", "genus_count", "species_count")))

# 绘制条形图
ggplot(counts_long, aes(x = Category, y = Count)) +
  geom_bar(stat = "identity", fill = "#3c77b0", width = 0.5) +  # 设置条形宽度
  coord_flip() +
  geom_text(aes(label = Count), hjust = -0.1) +  # 添加数值标签，hjust 调整位置
  scale_x_discrete(labels = c("class_count" = "Class", "order_count" = "Order", "family_count" = "Family", "genus_count" = "Genus", "species_count" = "Species")) +
  labs(title = "Non-empty Data Count for Arthropoda Categories",
       x = "",
       y = "Non-empty Count") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1),
    panel.grid.major = element_blank(),  # 取消主网格线
    panel.grid.minor = element_blank(),  # 取消次网格线
    axis.line = element_line(color = "black")  # 增加坐标轴
  )