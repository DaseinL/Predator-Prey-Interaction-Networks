if (!requireNamespace("BiMat", quietly = TRUE)) install.packages("BiMat")

library(BiMat)

# 读取CSV文件。假设文件名为'bipartite_matrix.csv'，
# 其中行表示虫子，列表示鸟类，每个格子代表在这类鸟的粪便中平均发现的虫的次数。
bipartite_data <- read.csv("processed_data.csv", row.names = 1)

# 由于数据已经是加权的，我们直接将其转换为二部网络
network <- as.bipartite_adjacency_matrix(bipartite_data)

# 创建二部网络对象，注意告知是加权网络
bp_network <- new_bipartite_network(network, is_weighted = TRUE)

# 运行模块化分析
modularity_analysis <- compute_modules(bp_network)

# 获取模块结果
modules <- modularity_analysis$modules

# 打印模块信息
print(modules)

# 可视化网络和模块
plot(modularity_analysis)

