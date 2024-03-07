data <- read.csv("/Users/daseinl/Downloads/009MRes Denian Li Project/Data/Ghana_fwh2_consensus_RAW.csv")
bird_meta <- read.csv("/Users/daseinl/Downloads/009MRes Denian Li Project/Data/MistNetting_wrangled_fullcoor.csv")

if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
library(dplyr)

# Get column named "OTU" and start only with G
columns_to_keep <- grep("^G[0-9]+|^OTU$", names(data), value = TRUE)
data_filtered <- data[, columns_to_keep]

# Combine column using the meta data to spacies
species_mapping <- setNames(bird_meta$BirdLife.Name.CORRECT, bird_meta$Specimen_number)
data_merged <- data.frame(OTU = data_filtered$OTU)
unique_species <- unique(species_mapping)
for(species in unique_species) {
  specimen_numbers <- names(species_mapping[species_mapping == species]) # bird code
  columns_to_merge <- names(data_filtered) %in% specimen_numbers
  if(sum(columns_to_merge) > 0) { # get averange value of each spacies
    average_data <- rowMeans(data_filtered[, columns_to_merge, drop = FALSE], na.rm = TRUE)
    data_merged[[species]] <- average_data
  }
}

write.csv(data_merged, "processed_data.csv", row.names = FALSE)

if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
if (!requireNamespace("ggraph", quietly = TRUE)) install.packages("ggraph")

library(igraph)
library(ggraph)

edge_list <- which(data_merged > 0, arr.ind = TRUE)
edges <- data.frame(from = rownames(data_merged)[edge_list[, 1]],
                    to = colnames(data_merged)[edge_list[, 2]])

g <- graph_from_data_frame(edges, directed = FALSE)

V(g)$type <- V(g)$name %in% rownames(data_merged)
V(g)$color <- ifelse(V(g)$type, "skyblue", "lightgreen")

png("bipartie_graph.png", width = 1024, height = 768)
plot(g, layout = layout.bipartite, vertex.label.color = "black", vertex.size = 2,
     vertex.label.cex = 0.001, edge.arrow.size = 0.5, edge.color = rgb(0, 0, 1, alpha = 0.003),
     vertex.frame.color = "gray", vertex.label.dist = 1.5)
dev.off()
