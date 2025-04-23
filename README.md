# Salmon and Canda Goose Distribution Analysis in British Columbia

## Overview
This R project analyzes the spatial distribution of salmon and canada goose observations in British Columbia and their relationship with **data centers** and environmental covariates. The analysis includes point pattern analysis, spatial transformation, quadrat testing, and Poisson point process modeling to understand factors influencing salmon distribution.


## Folder structure

```
.
├── data/
│   ├── salmon.csv                # Salmon Occurrence records (GBIF download)
│   ├── geese_canada_2010_2024.csv# Canada Goose Occurrence records (GBIF download)
│   ├── bc_data_centers.geojson   # Locations of data centres in BC
│   └── BC_Covariat.Rda           # R data file containing environmental covariates for BC
├── 1_preprocessing.R             # The preprocess script
├── 2_modelling_salmon.R          # The analysis script 
├── 2_modelling_goose.R           # The analysis script 
└── README.md                     
```


## Requirements
The following R packages are required:
- `spatstat` - For point pattern analysis
- `sf` - For spatial feature manipulation
- `sp` - For spatial data handling
- `terra` - For raster manipulation
- `raster` - For raster data processing


## Results
The analysis examines how salmon and canada goose distribution patterns relate to various environmental factors with special attention to the potential influence of data centers. Key findings include:

1. **Spatial Distribution Patterns**:
   - Significant spatial clustering of species detected (quadrat test p-value < 2.2e-16)
   - Comparison of species distribution between 2010 and 2024 showing population changes
   - Density maps highlighting population hotspots across British Columbia

2. **Environmental Relationships**:
   - Quantified relationships between species density and environmental covariates using rhohat functions
   - Analysis of distance to data centers as a potential influencing factor
   - Complex relationships with Human Footprint Index, elevation, water proximity, and forest cover

3. **Statistical Models**:
   - Final model includes quadratic terms for environmetal covariates
   - Significant improvement of full model over null model (based on AIC comparisons)
   - Statistical significance of data center distances in predicting species distribution

## Usage
To run the complete analysis:
1. Ensure all required packages are installed
2. Run the scripts in order:
   ```
   Rscript 1_preprocessing.R
   Rscript 2_modelling_salmon.R
   Rscript 2_modelling_goose.R
   ```
3. Alternatively, open and run scripts section by section in RStudio
