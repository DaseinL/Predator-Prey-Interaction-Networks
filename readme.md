# Predator-Prey Interaction Networks Analysis

This project analyzes predator-prey interaction networks between birds and arthropods across different habitat types (forest interior, forest edge, and agricultural areas). It focuses on understanding how forest-dependent birds contribute to pest control services in tropical ecosystems.

## Features
- Processes metabarcoding diet data
- Filters species by taxonomic group and abundance thresholds
- Analyzes bipartite network structure
- Quantifies pest control resilience through extinction simulations
- Generates interactive network visualizations
- Compares ecological metrics across land use types

## Requirements
- R 4.0+
- Required packages:
    ```r
    install.packages(c("dplyr", "tidyr", "tidyverse", "vegan", 
                   "bipartite", "ggplot2", "igraph", "iNEXT"))
    ```

## Workflow
### 1. Data Preparation
```r
# Load core datasets
data <- read.csv("./path/to/fwh_consensus.csv")            # Diet composition data
bird_meta <- read.csv("./path/to/mistnetting_data.csv")    # Bird metadata
insect_meta <- read.csv("./path/to/pest_annotations.csv")  # Arthropod trait data
```

### 2. Data Filtering
``` r
# Filter arthropod data
data_filtered <- filter_by_species(data)        # Remove non-arthropods and mites
data_RAA <- filter_by_RAA(data_filtered)        # Apply 0.3% abundance threshold

# Subset by habitat type
forest_data <- filter_by_site(data_RAA, bird_meta, "Forest Interior")
edge_data <- filter_by_site(data_RAA, bird_meta, "Forest Edge")
agri_data <- filter_by_site(data_RAA, bird_meta, "Agriculture")
```

### 3. Network Construction
```r
# Create interaction matrices
combined_data <- combine_data(bird_meta, filtered_data, full_dataset)

# Convert to bipartite format
bird_insect_web <- combined_data %>% 
  column_to_rownames("Insect_species") %>%
  as.matrix()
```

### 4. Network Analysis
```r
# Calculate key metrics
web_features <- web_feature(bird_insect_web, bird_forest_coverage, insect_roles)

# Visualize network structure
save_forest_category_bipartite(bird_forest_coverage, insect_roles, bird_insect_web, "Forest_Network")
```

### 5. Resilience Analysis
```r
# Simulate bird extinctions
extinction_results <- extinction_curve(bird_insect_web, insect_roles, 
                                      repetitions=1000, 
                                      sequence=extinction_order)

# Plot resilience curves
plot_extinction_curve(extinction_results)
```

## Key Functions
### Data Processing
- filter_by_species(): Filters arthropod species and removes ambiguous taxa
- filter_by_RAA(): Applies relative abundance threshold (0.3% of total reads)
- combine_data(): Merges sequencing data with ecological metadata

### Network Analysis
- web_feature(): Calculates pest-forest bird interaction strengths
- calculate_modularity(): Quantifies network modular structure
- statistic_richness(): Compares forest/non-forest bird diversity

### Visualization
save_forest_category_bipartite(): Generates colored bipartite networks

plot_species_metric(): Creates boxplots of species-level metrics

plot_network_metric(): Visualizes landscape-scale network properties

## Data Inputs
- Diet Data: CSV file containing OTU read counts per bird specimen
- Bird Metadata:
    - Species traits
    - Forest dependency scores
    - Spatial coordinates
- Arthropod Metadata:
    - Pest/non-pest classifications
    - Taxonomic information
    - Ecological roles

## Outputs
- Interaction networks (Gephi format)
- Resilience curves (JPEG)
- Metric comparisons (CSV)
- Network visualizations (JPEG)
- Statistical summaries (CSV)

## Reuse Potential
This framework can be adapted for:
- Food web studies in other ecosystems
- Landscape-scale conservation planning
- Ecological network robustness assessments
- Biodiversity-ecosystem function research