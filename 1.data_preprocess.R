data <- read.csv("/Users/daseinl/Downloads/009MRes Denian Li Project/Data/Ghana_fwh2_consensus_RAW.csv")
bird_meta <- read.csv("/Users/daseinl/Downloads/009MRes Denian Li Project/Data/MistNetting_wrangled_fullcoor.csv")

if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
library(dplyr)

# Get column named "OTU" and start only with G
columns_to_keep <- grep("^G[0-9]+|^OTU$", names(data), value = TRUE)
data_filtered <- data[, columns_to_keep]

write.csv(data_filtered, "preprocessed_data.csv", row.names = FALSE)

# Combine column using the meta data to spacies
species_mapping <- setNames(bird_meta$BirdLife.Name.CORRECT, bird_meta$Specimen_number)
data_merged <- data.frame(OTU = data_filtered$OTU)
unique_species <- unique(species_mapping)
for(species in unique_species) {
  specimen_numbers <- names(species_mapping[species_mapping == species]) # bird code
  columns_to_merge <- names(data_filtered) %in% specimen_numbers
  if(sum(columns_to_merge) > 0) { # get averange value of each spacies
    average_data <- rowSums(data_filtered[, columns_to_merge, drop = FALSE], na.rm = TRUE)
    data_merged[[species]] <- average_data
  }
}

write.csv(data_merged, "merged_data.csv", row.names = FALSE)

# produce a table containing bird forest coverage rate
bird_forest_coverage <- bird_meta %>% 
  group_by(BirdLife.Name.CORRECT) %>% 
  summarise(Avg_Forest = mean(Forest.Cover.200m, na.rm = TRUE))

write.csv(bird_forest_coverage, "bird_forest_coverage.csv", row.names = FALSE)





# if (!requireNamespace("bipartite", quietly = TRUE)) install.packages("bipartite")

# library(bipartite)

library(tidyverse) # version 2.0.0
# start load the required packages of bipartite package
library(permute) 
library(lattice)
library(vegan)
library(statnet.common)
library(sna)
library(bipartite) # version 2.18
library(ggplot2) # version 3.4.2l


### Bipartite network visualisation


# load in the data ### 
sample_completeness_data <- read_csv('./bird_forest_coverage.csv')
# (bipartite webs only works with integer matrices)

t1 <- read_csv("./merged_data.csv")
t1 <- t1 %>% column_to_rownames(var = 'OTU')

# these are binary webs because reads cannot be correlated to visits.../ number of pollen grains
# function for coloring the moth nodes yellow ###
colour_sel_1 <- function(plant_spp, first_col){
  col_names <- colnames(plant_spp)
  colours <- if_else(col_names %in% first_col,"darkgoldenrod2", "dodgerblue4", 'red')
  return(colours)
}

# Figure 1a Early summer
first_col_moth_t1 <- sample_completeness_data %>% filter(Avg_Forest >= 50)
first_col<- unique(first_col_moth_t1$BirdLife.Name.CORRECT)

png(filename = "bipartite_plot.png", width = 1000, height = 800, res = 300)

plotweb(t1, col.high = colour_sel_1(plant_spp = t1, first_col= first_col),
        text.rot=90, 
        col.low= 'grey',
        bor.col.interaction = NA,
        #low.lab.dis = 0.01,
        #high.lab.dis = 0.3,
        #low.spacing = 0.02,
        #high.spacing = 0.8,
        labsize = 1.5,
        low.lablength = 0,
        high.lablength = 0,
        col.interaction = colour_sel_1(plant_spp = t1, first_col= first_col),
        y.lim=c(-0.6,2.5),
        method = 'cca',
        #add = TRUE
)

dev.off()

#### making the barcharts that are part of figure 1 ####
# Stacked + percent
# fig_1_bar_data <- read.table("./figure_1_stacked_barchart.csv", header = T, sep = ',')

# ggplot(fig_1_bar_data, aes(fill=factor(Insect, levels=c("Bee", "Shared", "Moth")), y=Number,
#                            x=Time)) + 
#   geom_bar(position="fill", stat = 'identity') +
#   ylab('Preportion of plant species (%)') +
#   labs(x = "Time", fill = "Insect") +
#   scale_fill_manual(values = c("dodgerblue4", 'grey40', "darkgoldenrod2")) +
#   theme_bw()
