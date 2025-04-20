################################################################################
# Salmon and Goose Distribution Analysis in British Columbia
# 
# This script analyzes the spatial distribution of goose observations in British
# Columbia and their relationship with data centers
#
# Required libraries: spatstat, bcmaps, sf
################################################################################

# Load necessary libraries
library(spatstat)  # For point pattern analysis
library(bcmaps)    # For BC boundary data
library(sf)        # For spatial feature manipulation

################################################################################
# 1. Data Import and Preprocessing
################################################################################

# Function to load and preprocess data
load_data <- function(file_path, sep = ",", country_code = NULL) {
  # Import data
  data <- read.csv(file_path)
  
  # Filter out observations with missing coordinates
  data <- data[!is.na(data$decimalLongitude) & !is.na(data$decimalLatitude), ]
  
  # Filter by country if specified
  if (!is.null(country_code)) {
    data <- data[data$countryCode == country_code, ]
  }
  
  # Reset row names for clean data frame
  rownames(data) <- NULL
  
  return(data)
}

# Function to filter data by year
filter_by_year <- function(data, year) {
  year_data <- data[data$year == year, ]
  year_data <- year_data[!is.na(year_data$decimalLongitude) & !is.na(year_data$decimalLatitude), ]
  rownames(year_data) <- NULL
  return(year_data)
}

# Function to create point pattern (ppp) objects
create_ppp <- function(x, y, window = NULL) {
  if (is.null(window)) {
    window <- owin(xrange = range(x), yrange = range(y))
  }
  return(ppp(x = x, y = y, window = window))
}

# Import data centers
data_centers <- st_read("data/bc_data_centers.geojson")

# Import BC Covariates
load("data/BC_Covariat.Rda")

################################################################################
# 2. Spatial Transformation and BC Window Creation
################################################################################

# Get BC boundary
bc <- bc_bound()

# Create BC window in Albers projection
bc_window_albers <- as.owin(bc)

# Function to transform points to BC Albers projection and filter to BC only
transform_points_to_bc <- function(data, bc_window) {
  # Create simple feature from points (WGS84)
  points_sf <- st_as_sf(data.frame(x = data$decimalLongitude, y = data$decimalLatitude), 
                        coords = c("x", "y"), crs = 4326)
  
  # Transform to BC Albers projection
  points_albers <- st_transform(points_sf, 3005)
  
  # Extract coordinates
  coords_albers <- st_coordinates(points_albers)
  x_albers <- coords_albers[, "X"]
  y_albers <- coords_albers[, "Y"]
  
  # Create point pattern in BC Albers projection
  all_points_ppp <- ppp(x = x_albers, y = y_albers, window = bc_window)
  
  # Extract points within BC
  inside_bc <- which(inside.owin(x_albers, y_albers, w = bc_window))
  bc_points_x <- x_albers[inside_bc]
  bc_points_y <- y_albers[inside_bc]
  bc_points_ppp <- ppp(x = bc_points_x, y = bc_points_y, window = bc_window)
  
  return(list(
    all_points_ppp = all_points_ppp,
    bc_points_ppp = bc_points_ppp,
    x_albers = x_albers,
    y_albers = y_albers,
    bc_points_x = bc_points_x,
    bc_points_y = bc_points_y
  ))
}

# Transform data centers to BC Albers projection
data_centers_albers <- st_transform(data_centers, 3005)
data_center_coords <- st_coordinates(data_centers_albers)


################################################################################
# 3. Analysis Functions
################################################################################

# Function to create density plot
create_density_plot <- function(data_ppp, title) {
  lambda_u_hat <- density(data_ppp)
  plot(lambda_u_hat, main=title)
  plot(data_ppp, pch = 16, cex = 0.7, add = TRUE)
  plot(data_ppp, pch = 16, cex = 0.3, cols = "white", add = TRUE)
  return(lambda_u_hat)
}


# Function to perform quadrat test
perform_quadrat_test <- function(data_ppp, nx = 4, ny = 4, plot = TRUE) {
  Q <- quadratcount(data_ppp, nx = nx, ny = ny)
  
  if (plot) {
    plot(intensity(Q, image = TRUE))
    plot(data_ppp, pch = 16, use.marks = FALSE, cols = "#046C9A", add = TRUE)
    plot(Q, col = "white", add = TRUE)
  }
  
  test_result <- quadrat.test(Q)
  return(list(Q = Q, test_result = test_result))
}

# Function to calculate distance to nearest data center
calculate_center_distances <- function(points_sf, data_centers_albers) {
  distances <- st_distance(points_sf, data_centers_albers)
  min_distances <- apply(distances, 1, min)
  return(min_distances)
}


# Function to create distance visualization
create_distance_visualization <- function(point_ppp, min_distances, data_center_coords, title) {
  ppp_marked <- point_ppp %mark% as.numeric(min_distances)
  
  # Create smooth image of distances
  dist_im <- Smooth(ppp_marked, sigma = 20000, dimyx = c(300, 300), at = "pixels")
  
  # Plot
  plot(dist_im, main = title)
  plot(point_ppp, add = TRUE, pch = 16, cex = 0.6)
  points(data_center_coords[, "X"], data_center_coords[, "Y"], 
         col = "red", pch = 4, cex = 1.2)
  
  return(dist_im)
}


# Function to create feature visualization
create_feature_plot <- function(data, title, data_ppp = NA){
  plot(data, main = title)
  plot(data_ppp, pch = 16, cex = 0.7, add = TRUE)
  plot(data_ppp, pch = 16, cex = 0.3, cols = "white", add = TRUE)
  points(data_center_ppp$x, data_center_ppp$y, col = "green", pch = 4, cex = 1.2)
}



######################################
# I. EDA
######################################

######################################
# goose
######################################

# Import data
data <- load_data("data/geese_canada_2010_2024.csv")

# Check temporal range
cat("Temporal range of data:", min(data$year, na.rm = TRUE), 
    "to", max(data$year, na.rm = TRUE), "\n")

year_data <- filter_by_year(data, 2024)

# Transform points to BC Albers and filter to BC only
transformed_points <- transform_points_to_bc(year_data, bc_window_albers)
inside_bc <- transformed_points$bc_points_ppp
summary(inside_bc)
cat("The goose population in the year 2024: ", npoints(inside_bc)) #1743

data_2010 <- filter_by_year(data, 2010)
trans <- transform_points_to_bc(data_2010, bc_window_albers)
inside_bc_2010 <- trans$bc_points_ppp
cat("The goose population in the year 2010: ", npoints(inside_bc_2010)) #6232

# Visualize the density
create_density_plot(inside_bc_2010, "Goose Distribution in 2010")
create_density_plot(inside_bc, "Goose Distribution in 2024")


# Intensity
intensity(inside_bc_2010) #6.616102e-09
intensity(inside_bc) #1.850428e-09


# Function to perform quadrat test
quadrat_result <- perform_quadrat_test(inside_bc, nx=4, ny=4)
quadrat_result$test_result # p-value < 2.2e-16

######################################
# Distance from data centers
######################################

# Convert to sf for distance calculation
points_sf <- st_as_sf(
  data.frame(x = inside_bc$x, y = inside_bc$y),
  coords = c("x", "y"), crs = 3005
)
min_distances <- calculate_center_distances(points_sf, data_centers_albers) 

inside_bc_marked <- inside_bc %mark% as.numeric(min_distances)

# Create distance map
data_center_ppp <- ppp(
  x = data_center_coords[, "X"],
  y = data_center_coords[, "Y"],
  window = bc_window_albers
)
dist_im <- distmap(data_center_ppp)

# Plot distance map
plot(dist_im, main = "Distance to nearest data center (m)")
plot(inside_bc_marked, add = TRUE, pch = 16, cex = 0.6)
points(data_center_ppp$x, data_center_ppp$y, col = "red", pch = 4, cex = 1.2)
plot(Window(inside_bc_marked), add = TRUE)

######################################
# HFI (Human Footprint Index)
######################################
hfi <- DATA$HFI
create_feature_plot(hfi, "Human Footprint Index", inside_bc)

######################################
# Elevation
######################################
elev <- DATA$Elevation
create_feature_plot(elev, "Elevation (m)", inside_bc)

######################################
# Distance from water
######################################
water <- DATA$Dist_Water
create_feature_plot(water, "Distance water", inside_bc)

######################################
# Forest
######################################
forest <- DATA$Forest
create_feature_plot(forest, "Forest Cover (%)", inside_bc)


######################################
# II. Modelling
######################################
# Estimate intensity as function of distance
rho_distance <- rhohat(inside_bc, dist_im, confidence = 0.95)
plot(rho_distance, main = "Goose Density vs Distance to Data Centers")

rho_forest <- rhohat(inside_bc_marked, forest, confidence = 0.95)
plot(rho_forest, main = "Goose Density vs Forest Cover")

rho_water <- rhohat(inside_bc_marked, water, confidence = 0.95)
plot(rho_water, main = "Goose Density vs Water Distance", legend="topright")

rho_hfi <- rhohat(inside_bc_marked, hfi, confidence = 0.95)
plot(rho_hfi, main = "Goose Density vs HFI")

rho_elev <- rhohat(inside_bc_marked, elev, confidence = 0.95)
plot(rho_elev, main = "Goose Density vs Elevation")

# Colinearity -> not too strong
cor.im(dist_im, hfi, elev, forest, water, use="pairwise.complete.obs")



fit_final<- ppm(
  inside_bc ~ dist_im + hfi + elev + forest + water +  I(forest^2),
  data = list(
    dist_im = dist_im,
    hfi = hfi,  
    elev = elev,  
    water = water,  
    forest = forest  
  ),
  control = list(maxit = 1000)
)
summary(fit_final)

plot(fit_final, log = TRUE, se=FALSE, superimpose=FALSE, n=200)
plot(inside_bc,
     pch = 16,
     cex = 0.6,
     add = TRUE)
plot(inside_bc,
     pch = 16,
     cex = 0.3,
     cols = "white",
     add = TRUE)

## Model selection
fit_null <- ppm(inside_bc ~ 1)
# summary(fit_null)
AIC(fit_final)
AIC(fit_null)
anova(fit_final, fit_null)
quadrat.test(fit_final, nx = 4, ny = 4) # p-value < 2.2e-16

#Visualise residuals --HFI might have some NA values, so there a error message
res <- residuals(fit_final)
plot(res,
     cols = "transparent")

#Calculate the partial residuals as a function of elevation
par_res_dist <- parres(fit_final, "dist_im")
plot(par_res_dist,
     legend = FALSE,
     lwd = 2,
     main = "",
     xlab = "Distance")

par_res_hfi <- parres(fit_final, "hfi")
plot(par_res_hfi,
     legend = FALSE,
     lwd = 2,
     main = "",
     xlab = "HFI")

par_res_elev <- parres(fit_final, "elev")
plot(par_res_elev,
     legend = FALSE,
     lwd = 2,
     main = "",
     xlab = "Elevation")

par_res_water <- parres(fit_final, "water")
plot(par_res_water,
     legend = FALSE,
     lwd = 2,
     main = "",
     xlab = "Water Distance")

par_res_forest <- parres(fit_final, "forest")
plot(par_res_forest,
     legend = FALSE,
     lwd = 2,
     main = "",
     xlab = "Forest Cover")
