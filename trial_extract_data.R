row_num <- 500
csv_name <- paste("first", row_num, "rows.csv", sep="_")
bipartie_png_name <- paste("bipartie_image_of_first", row_num, "rows.png", sep='_')

# get data
data <- read.csv("/Users/daseinl/Downloads/009MRes Denian Li Project/Data/Ghana_fwh2_consensus_RAW.csv")

# filter the bat and empty column
columns_to_keep <- c(names(data)[1], grep("^G", names(data), value = TRUE)) # keep only interaction
data <- data[, columns_to_keep] 
columns_to_remove <- grep("^GB", names(data), value = TRUE) # get rid of bat
data <- subset(data, select = -match(columns_to_remove, names(data)))

first_rows <- head(data, row_num)
data_clean <- first_rows[, colSums(is.na(data)) < nrow(data)]
write.csv(data_clean, csv_name, row.names = FALSE)

cat(paste0("The first ", row_num, " rows and ", ncol(data_clean), " columns have been saved to first_", row_num, "_rows.csv\n"))



# start visualizing a bipartie graph
if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
if (!requireNamespace("ggraph", quietly = TRUE)) install.packages("ggraph")

library(igraph)
library(ggraph)

interaction_matrix <- read.csv(csv_name)

edge_list <- which(interaction_matrix > 0, arr.ind = TRUE)
edges <- data.frame(from = rownames(interaction_matrix)[edge_list[, 1]],
                    to = colnames(interaction_matrix)[edge_list[, 2]])

g <- graph_from_data_frame(edges, directed = FALSE)

V(g)$type <- V(g)$name %in% rownames(interaction_matrix)
V(g)$color <- ifelse(V(g)$type, "skyblue", "lightgreen")

png(bipartie_png_name, width = 1024, height = 768)
plot(g, layout = layout.bipartite, vertex.label.color = "black", vertex.size = 2,
     vertex.label.cex = 0.001, edge.arrow.size = 0.5, edge.color = rgb(0, 0, 1, alpha = 0.01),
     vertex.frame.color = "gray", vertex.label.dist = 1.5)
dev.off()

# modular analysis

