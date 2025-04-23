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
├── scripts/
│   ├── 1_preprocessing.R         # The preprocess script
│   ├── 2_modelling_salmon.R      # The analysis script 
│   └── 2_modelling_goose.R       # The analysis script 
└── README.md                     # You are here
```


## Requirements
The following R packages are required:
- `spatstat` - For point pattern analysis
- `sf` - For spatial feature manipulation
- `sp` - For spatial data handling
- `terra` - For raster manipulation
- `raster` - For raster data processing


## Results
The analysis examines how salmon and canada goose distribution patterns relate to various environmental factors with special attention to the potential influence of data centers.

