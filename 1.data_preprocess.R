rm(list=ls())
setwd(dir="/Users/daseinl/Coding/Predator-Prey-Interaction-Networks")
getwd()

if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")

library(dplyr)
library(tidyr)
library(tidyverse) # version 2.0.0

data <- read.csv("/Users/daseinl/Downloads/009MRes Denian Li Project/Data/Ghana_fwh2_consensus_RAW.csv")
bird_meta <- read.csv("/Users/daseinl/Downloads/009MRes Denian Li Project/Data/MistNetting_wrangled_fullcoor.csv")

# The code to filter arthropods and ambiguous species by Ariane
arthropoda_and_species_only <- data %>%
  dplyr::filter(phylum == "Arthropoda") %>% # only keep arthropods
  dplyr::filter(order != "Mesostigmata"&  # exclude arachnids which are mites and therefore not dietary items
                  order != "Sarcoptiformes"&
                  order != "Trombidiformes") %>%
  filter(!is.na(species)) %>%  # only keep OTUs identified to genus or species level
  filter(!is.na(genus)) %>%
  distinct()

# Get column named "OTU" and start only with G
columns_to_keep <- grep("^G[0-9]+|^OTU$", names(arthropoda_and_species_only), value = TRUE)
data_filtered <- arthropoda_and_species_only[, columns_to_keep]

# The RRA code to eliminate bias from small data by Ariane
num_col <- ncol(data_filtered)
data_RAA <- data_filtered %>%
  dplyr::select(1:num_col) %>%
  column_to_rownames(var = "OTU") %>%
  mutate(across(1:num_col - 1, ~ ifelse(. < 0.003 * sum(.), 0, .))) %>% # apply RRA filtering threshold (here 0.3% of total reads)
  filter(rowSums(select(., 1:num_col - 1)) != 0) %>%  # filter out OTUs that do not occur following the RRA filtering
  rownames_to_column(var = "OTU")

write.csv(data_RAA, "preprocessed_data.csv", row.names = FALSE)


# Combine column using the meta data to species
species_mapping <- setNames(bird_meta$BirdLife.Name.CORRECT, bird_meta$Specimen_number)
data_merged <- data.frame(OTU = data_RAA$OTU)
unique_species <- unique(species_mapping)
for(species in unique_species) {
  specimen_numbers <- names(species_mapping[species_mapping == species]) # bird code
  columns_to_merge <- names(data_RAA) %in% specimen_numbers
  if(sum(columns_to_merge) > 0) { # get average value of each species
    average_data <- rowSums(data_RAA[, columns_to_merge, drop = FALSE], na.rm = TRUE)
    data_merged[[species]] <- average_data
  }
}

# produce a table containing bird forest coverage rate
bird_forest_coverage <- bird_meta %>% 
  group_by(BirdLife.Name.CORRECT) %>% 
  summarise(average.Distance.to.Forest.km = mean(Distance.to.Forest.km))

# sort bird coverage and net information
bird_forest_coverage_sorted <- bird_forest_coverage[order(bird_forest_coverage$average.Distance.to.Forest.km), ]
valid_columns <- bird_forest_coverage_sorted$BirdLife.Name.CORRECT %in% colnames(data_merged)
data_merged_reordered <- data_merged[, bird_forest_coverage_sorted$BirdLife.Name.CORRECT[valid_columns]]
data_merged_reordered <- cbind(OTU=data_merged[["OTU"]], data_merged_reordered)

write.csv(data_merged_reordered, "./merged_data.csv", row.names = FALSE)
write.csv(bird_forest_coverage_sorted, "./bird_forest_coverage.csv", row.names = FALSE)


if (!requireNamespace("permute", quietly = TRUE)) install.packages("permute")
if (!requireNamespace("lattice", quietly = TRUE)) install.packages("lattice")
if (!requireNamespace("vegan", quietly = TRUE)) install.packages("vegan")
if (!requireNamespace("statnet.common", quietly = TRUE)) install.packages("statnet.common")
if (!requireNamespace("sna", quietly = TRUE)) install.packages("sna")
if (!requireNamespace("bipartite", quietly = TRUE)) install.packages("bipartite")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")


# start load the required packages of bipartite package
library(permute) 
library(lattice)
library(vegan)
library(statnet.common)
library(sna)
library(bipartite) # version 2.18
library(ggplot2) # version 3.4.2l


### Bipartite network visualization

# load in the data ### 
sample_completeness_data <- read_csv('./bird_forest_coverage.csv', show_col_types = FALSE)

bird_insect_web <- read_csv("./merged_data.csv", show_col_types = FALSE)
bird_insect_web <- bird_insect_web %>% column_to_rownames(var = 'OTU')

# function for coloring the forest bird yellow ###
colour_sel_1 <- function(plant_spp, first_col){
  col_names <- colnames(plant_spp)
  colours <- if_else(col_names %in% first_col,"#238E23", "darkgoldenrod2", 'red')
  return(colours)
}

first_col_moth_t1 <- sample_completeness_data %>% filter(average.Distance.to.Forest.km < 0)
first_col<- unique(first_col_moth_t1$BirdLife.Name.CORRECT)

jpeg(paste0("bipartite_plot.jpg"), width = 2000, height = 1000, units = 'px', res = 300) 
plotweb(bird_insect_web, col.high = colour_sel_1(plant_spp = bird_insect_web, first_col= first_col),
        text.rot=90, 
        col.low= 'grey',
        bor.col.interaction = NA,
        #low.lab.dis = 0.01,
        #high.lab.dis = 0.3,
        #low.spacing = 0.02,
        #high.spacing = 0.8,
        labsize = 0.5,
        low.lablength = 0,
        # high.lablength = 0,
        col.interaction = colour_sel_1(plant_spp = bird_insect_web, first_col= first_col),
        y.lim=c(-0.6,2.5),
        method = 'normal',
        #add = TRUE
)
dev.off()

#### making the barcharts that are part of figure 1 ####
# Stacked + percent
# fig_1_bar_data <- read.table("./figure_1_stacked_barchart.csv", header = T, sep = ',')
# 
# ggplot(fig_1_bar_data, aes(fill=factor(Insect, levels=c("Bee", "Shared", "Moth")), y=Number,
#                            x=Time)) +
#   geom_bar(position="fill", stat = 'identity') +
#   ylab('Preportion of plant species (%)') +
#   labs(x = "Time", fill = "Insect") +
#   scale_fill_manual(values = c("dodgerblue4", 'grey40', "darkgoldenrod2")) +
#   theme_bw()


## Module analysis

module_result <- computeModules(t(bird_insect_web), method="Beckett", deleteOriginalFiles = TRUE,
               steps = 1000000, tolerance = 1)
# metaComputeModules(module_result, method="Beckett", deleteOriginalFiles = TRUE,
#                    steps = 1000000, tolerance = 1)
# print(module_result$modularity)

jpeg("bipartite_matrix.jpg", width = 5000, height = 700, units = 'px', res = 300) 
visweb(t(bird_insect_web))
dev.off()

jpeg("bipartite_matrix_module.jpg", width = 5000, height = 700, units = 'px', res = 300) 
visweb(module_result@moduleWeb)
dev.off()

jpeg(paste0("bipartite_module_plot.jpg"), width = 2000, height = 1000, units = 'px', res = 300) 
plotModuleWeb(module_result, plotModules = TRUE,
              rank = FALSE, weighted = TRUE, displayAlabels = TRUE,
              displayBlabels = TRUE, labsize = 0.5, xlabel = "", ylabel = "",
              square.border = "white", fromDepth = 0, upToDepth = -1)
dev.off()

library(ggplot2)

# 假设bird_insect_web为虫和鸟之间的交互矩阵，sample_completeness_data为鸟类数据
# 判断鸟类是否为森林鸟类
is_forest_bird <- sample_completeness_data$average.Distance.to.Forest.km == 0

# 标记森林鸟类和非森林鸟类
forest_birds <- names(is_forest_bird)[is_forest_bird]
non_forest_birds <- names(is_forest_bird)[!is_forest_bird]

# 计算每种虫子被森林鸟类和非森林鸟类捕捉的情况
forest_only <- rowSums(bird_insect_web[, forest_birds]) > 0 & rowSums(bird_insect_web[, non_forest_birds]) == 0
non_forest_only <- rowSums(bird_insect_web[, non_forest_birds]) > 0 & rowSums(bird_insect_web[, forest_birds]) == 0
both_types <- rowSums(bird_insect_web[, forest_birds]) > 0 & rowSums(bird_insect_web[, non_forest_birds]) > 0

# 准备数据用于绘图
data_for_plot <- data.frame(
  Category = c("Forest Only", "Non-Forest Only", "Both"),
  Count = c(sum(forest_only), sum(non_forest_only), sum(both_types))
)

# 绘制比例图
p <- ggplot(data_for_plot, aes(x = "", y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("Forest Only" = "#238E23", "Non-Forest Only" = "darkgoldenrod2", "Both" = "grey")) +
  theme_void() +
  labs(title = "Insect Capture by Bird Type", fill = "Category")

ggsave("proportion_plot.png", plot = p, width = 10, height = 6, units = "in")
