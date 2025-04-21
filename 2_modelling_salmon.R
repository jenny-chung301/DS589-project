################################################################################
# Salmon and Goose Distribution Analysis in British Columbia
# 
# This script analyzes the spatial distribution of salmon observations in British
# Columbia and their relationship with data centers
#
# Required libraries: spatstat, bcmaps, sf
################################################################################

# Load necessary libraries
library(spatstat)  # For point pattern analysis
library(sf)        # For spatial feature manipulation

################################################################################
# 1. Data Import and Preprocessing
################################################################################

# Function to load and preprocess data
load_data <- function(file_path, sep = ",", country_code = NULL) {
  # Import data
  data <- read.csv(file_path, sep = sep)
  
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

# Get a Coordinate Reference System (CRS)
target_crs <- proj4string(DATA$Window)
win <- as.owin(st_as_sf(DATA$Window))

# Function to transform points to BC Albers projection and filter to BC only
transform_points_to_bc <- function(data, crs=target_crs, window=win) {
  # Create simple feature from points (WGS84)
  points_sf <- st_as_sf(data.frame(x = data$decimalLongitude, y = data$decimalLatitude), 
                        coords = c("x", "y"), crs = 4326)
  
  
  # Transform to BC Albers projection
  points_albers <- st_transform(points_sf, crs = crs)
  
  # Extract coordinates
  coords_albers <- st_coordinates(points_albers)
  x_albers <- coords_albers[, "X"]
  y_albers <- coords_albers[, "Y"]
  
  # Create point pattern in BC Albers projection
  all_points_ppp <- ppp(x = x_albers, y = y_albers, window = window)
  
  # Extract points within BC
  inside_bc <- which(inside.owin(x_albers, y_albers, w = window))
  bc_points_x <- x_albers[inside_bc]
  bc_points_y <- y_albers[inside_bc]
  bc_points_ppp <- ppp(x = bc_points_x, y = bc_points_y, window = window)
  
  return(list(
    all_points_ppp = all_points_ppp,
    bc_points_ppp = bc_points_ppp,
    x_albers = x_albers,
    y_albers = y_albers,
    bc_points_x = bc_points_x,
    bc_points_y = bc_points_y
  ))
}


# Transform data centers to ppp
data_centers_albers <- st_transform(data_centers, target_crs)
data_center_coords <- st_coordinates(data_centers_albers)
data_center_ppp <- ppp(
  x = data_center_coords[, "X"],
  y = data_center_coords[, "Y"],
  window = win
)
dc_dist <- distmap(data_center_ppp)

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
  # points(data_center_ppp$x, data_center_ppp$y, col = "green", pch = 4, cex = 1.2)
}



######################################
# I. EDA
######################################

######################################
# Salmon
######################################

# Import data
data <- load_data("data/salmon.csv")

# Check temporal range
cat("Temporal range of data:", min(data$year, na.rm = TRUE), 
    "to", max(data$year, na.rm = TRUE), "\n")


######################################
# Comparison with 2010
######################################

year_data_2010 <- filter_by_year(data, 2010)

# Transform points to BC Albers and filter to BC only
transformed_points_2010 <- transform_points_to_bc(year_data_2010)
inside_bc_2010 <- transformed_points_2010$bc_points_ppp
summary(inside_bc_2010)
cat("The salmon population: ", npoints(inside_bc_2010))

# Visualize the density
create_density_plot(inside_bc_2010, "Salman Distribution")

######################################
# 2024
######################################
year_data <- filter_by_year(data, 2024)

# Transform points to BC Albers and filter to BC only
transformed_points <- transform_points_to_bc(year_data)
inside_bc <- transformed_points$bc_points_ppp
summary(inside_bc)
cat("The salmon population: ", npoints(inside_bc))

# Visualize the density
create_density_plot(inside_bc, "Salman Distribution")

# Intensity
intensity(inside_bc)

# Function to perform quadrat test
quadrat_result <- perform_quadrat_test(inside_bc, nx=4, ny=4)
quadrat_result$test_result # p-value < 2.2e-16


######################################
# Distance from data centers
######################################

# Plot distance map
plot(dist_im, main = "Distance to nearest data center (m)")
points(data_center_ppp$x, data_center_ppp$y, col = "red", pch = 4, cex = 1.2)
plot(Window(inside_bc), add = TRUE)


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
k_inside_bc <- Kest(inside_bc)

plot(k_inside_bc,
     main = "",
     lwd = 2)
######################################
# Evaluating Clustering under an Inhomogeneous Intensity Assumption
######################################
lambda_salmon_pos <- density(inside_bc,
                             sigma=bw.ppl,
                             positive=TRUE)
pcf_salmon_inhom <- envelope(inside_bc,
                             pcfinhom,
                             simulate = expression(rpoispp(lambda_salmon_pos)),
                             rank = 1,
                             nsim = 19)
#par(mfrow = c(1,2))
plot(pcf_salmon_inhom,
     xlim = c(0,20000),
     main = "",
     lwd = 2)


######################################
# Relationship with covariates
######################################

# Estimate intensity as function of distance
rho_dc <- rhohat(inside_bc, dc_dist, confidence = 0.95)
plot(rho_dc, main = "Salmon Density vs Distance to Data Centers")

rho_forest <- rhohat(inside_bc, forest, confidence = 0.95)
plot(rho_forest, main = "Salmon Density vs Forest Cover")

rho_water <- rhohat(inside_bc_marked, water, confidence = 0.95)
plot(rho_water, main = "Salmon Density vs Water Distance")

rho_hfi <- rhohat(inside_bc_marked, hfi, confidence = 0.95)
plot(rho_hfi, main = "Salmon Density vs HFI")

rho_elev <- rhohat(inside_bc_marked, elev, confidence = 0.95)
plot(rho_elev, main = "Salmon Density vs Elevation")

# Colinearity -> not too strong
cor.im(dist_im, hfi, elev, forest, water, use="pairwise.complete.obs")

# 
# fit <- ppm(
#   inside_bc ~ dist_im + hfi + elev+ water + forest + I(forest^2) + I(elev^2),
#   covariates = list(
#     dist_im = dist_im,
#     hfi = hfi,  
#     elev = elev,  
#     water = water,  
#     forest = forest  
#   ),
#   control = list(maxit = 5000)
# )
# summary(fit)

fit_final<- ppm(
  inside_bc ~ dc_dist + hfi + elev + forest + water + I(hfi^2)  + I(elev^2),
  data = list(
    dc_dist = dc_dist,
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
AIC(fit_final)
AIC(fit_null)
anova(fit_final, fit_null)
quadrat.test(fit_final, nx = 4, ny = 4) # p-value < 2.2e-16

#Visualise residuals
res <- residuals(fit_final)
plot(res,
     cols = "transparent")
## without hfi, interpolate on hfi,  

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
