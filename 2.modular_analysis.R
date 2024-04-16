# 加载igraph包
library(igraph)

# 读取CSV文件
data <- read.csv("./merged_data.csv", header = TRUE, row.names = 1)

# 构建网络
net <- graph.data.frame(data, directed = FALSE)

# 模块化分析
mem <- cluster_leading_eigen(net)
membership <- membership(mem)

# 可视化网络和模块
V(net)$color <- rainbow(max(membership))[membership]
V(net)$frame.color <- "white"
V(net)$size <- 5
plot(net, vertex.label.cex = 0.7, layout = layout_with_fr)
legend("bottomleft", legend = unique(membership), pch = 21,
       col = rainbow(max(membership)), pt.bg = rainbow(max(membership)),
       pt.cex = 2, bty = "n", ncol = 3)
