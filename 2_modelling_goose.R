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
library(sp)
library(terra)
library(raster)

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

#Calculate the partial residuals as a function of elevation
create_partial_residual_plt <- function(model, covariate, xlab){
  par_res <- parres(model, covariate)
  plot(par_res,
       legend = FALSE,
       lwd = 2,
       main = "",
       xlab = xlab)
}

######################################
# I. EDA
######################################

######################################
# Goose
######################################

# Import data
data <- load_data("data/geese_canada_2010_2024.csv")

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
cat("The Goose population: ", npoints(inside_bc_2010))

# Visualize the density
create_density_plot(inside_bc_2010, "Goose Distribution")

######################################
# 2024
######################################
year_data <- filter_by_year(data, 2024)

# Transform points to BC Albers and filter to BC only
transformed_points <- transform_points_to_bc(year_data)
inside_bc <- transformed_points$bc_points_ppp
summary(inside_bc)
cat("The goose population: ", npoints(inside_bc))

# Visualize the density
create_density_plot(inside_bc, "Goose Distribution")

# Intensity
intensity(inside_bc)

# Function to perform quadrat test
quadrat_result <- perform_quadrat_test(inside_bc, nx=4, ny=4)
quadrat_result$test_result # p-value < 2.2e-16

lambda_goose_pos <- density(inside_bc,
                             sigma=bw.ppl,
                             positive=TRUE)
pcf_goose_inhom <- envelope(inside_bc,
                             pcfinhom,
                             simulate = expression(rpoispp(lambda_goose_pos)),
                             rank = 1,
                             nsim = 19)
par(mfrow = c(1,2))
plot(pcf_goose_inhom,
     xlim = c(0,10000),
     main = "",
     lwd = 2,
     legendargs = list(cex = 0.5))


######################################
# Distance from data centers
######################################

# Plot distance map
plot(dc_dist, main = "Distance to nearest data center (m)")
points(data_center_ppp$x, data_center_ppp$y, col = "red", pch = 4, cex = 1.2)
plot(Window(inside_bc), add = TRUE)


######################################
# HFI (Human Footprint Index)
######################################
hfi <- DATA$HFI
create_feature_plot(hfi, "Human Footprint Index", inside_bc)

hfi_raster <- rast(hfi)
hfi_median <- global(hfi_raster, fun = median, na.rm = TRUE)[[1]]

# Replace all NA values with the median
hfi_filled <- subst(hfi_raster, NA, hfi_median)

#Convert to raster::RasterLayer 
hfi_raster_filled <- raster(hfi_filled)

#crop to match analysis window
win_bbox <- as.rectangle(inside_bc$window)
hfi_crop <- crop(hfi_raster_filled,
                 extent(win_bbox$xrange[1], win_bbox$xrange[2],
                        win_bbox$yrange[1], win_bbox$yrange[2]))

#Convert to spatstat's im object (manual)
hfi_matrix <- as.matrix(hfi_crop)
x_coords <- seq(from = xmin(hfi_crop), to = xmax(hfi_crop), length.out = ncol(hfi_matrix))
y_coords <- seq(from = ymin(hfi_crop), to = ymax(hfi_crop), length.out = nrow(hfi_matrix))
hfi_matrix_flipped <- hfi_matrix[nrow(hfi_matrix):1, ]
hfi_im <- im(hfi_matrix_flipped, xcol = x_coords, yrow = y_coords)
plot(hfi_im)

######################################
# Elevation
######################################
elev <- DATA$Elevation
create_feature_plot(elev, "Elevation (m)", inside_bc)
#convert spatstat im to terra raster (if needed)
elev_raster <- rast(elev)
elev_median <- global(elev_raster, fun = median, na.rm = TRUE)[[1]]

# replace all NA values with the median
elev_filled <- subst(elev_raster, NA, elev_median)

# Convert to raster::RasterLayer (for as.im compatibility)
elev_raster_filled <- raster(elev_filled)

# crop to match analysis window
win_bbox <- as.rectangle(inside_bc$window)
elev_crop <- crop(elev_raster_filled,
                 extent(win_bbox$xrange[1], win_bbox$xrange[2],
                        win_bbox$yrange[1], win_bbox$yrange[2]))

# Convert to spatstat's im object (manual method – most reliable)
elev_matrix <- as.matrix(elev_crop)
x_coords <- seq(from = xmin(elev_crop), to = xmax(elev_crop), length.out = ncol(elev_matrix))
y_coords <- seq(from = ymin(elev_crop), to = ymax(elev_crop), length.out = nrow(elev_matrix))
elev_matrix_flipped <- elev_matrix[nrow(elev_matrix):1, ]
elev_im <- im(elev_matrix_flipped, xcol = x_coords, yrow = y_coords)

plot(elev_im)
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
rho_dc <- rhohat(inside_bc, dc_dist, confidence = 0.95)
plot(rho_dc, main = "Goose Density vs Distance to Data Centers")

rho_forest <- rhohat(inside_bc, forest, confidence = 0.95)
plot(rho_forest, main = "Goose Density vs Forest Cover")

rho_water <- rhohat(inside_bc, water, confidence = 0.95)
plot(rho_water, main = "Goose Density vs Water Distance")

rho_hfi <- rhohat(inside_bc, hfi_im, confidence = 0.95)
plot(rho_hfi, main = "Goose Density vs HFI")

rho_elev <- rhohat(inside_bc, elev, confidence = 0.95)
plot(rho_elev, main = "Goose Density vs Elevation")

# Colinearity -> not too strong
cor.im(dc_dist, hfi_im, elev, forest, water, use="pairwise.complete.obs")

fit_final<- ppm(
  inside_bc ~ dc_dist + hfi_im + elev_im + forest + water+I(elev_im^2),
  data = list(
    dc_dist = dc_dist,
    hfi_im = hfi_im,  
    elev_im = elev_im,  
    water = water,  
    forest = forest  
  ),
  control = list(maxit = 10000)
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

# Partial residual plot
create_partial_residual_plt(fit_final, "dc_dist", "Distance to Data Center")
create_partial_residual_plt(fit_final, "hfi_im", "HFI")
create_partial_residual_plt(fit_final, "elev_im", "Elevation")
create_partial_residual_plt(fit_final, "forest", "Forest Cover")
create_partial_residual_plt(fit_final, "water", "Distance to Water")
