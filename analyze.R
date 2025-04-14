################################################################################
# Salmon Distribution Analysis in British Columbia
# 
# This script analyzes the spatial distribution of salmon observations in British
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

# Function to load and preprocess salmon observation data
load_data <- function(file_path, country_code = NULL) {
  # Import data
  data <- read.csv(file_path, sep = "\t")
  
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

# Import data
salmon_data_ca <- load_data("data/salmon_0009832-250402121839773.csv", "CA")

# Check temporal range
cat("Temporal range of data:", min(salmon_data_ca$year, na.rm = TRUE), 
    "to", max(salmon_data_ca$year, na.rm = TRUE), "\n")

# Import data centers
data_centers <- st_read("data/bc_data_centers.geojson")


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

# Function to calculate distance to nearest data center
calculate_center_distances <- function(salmon_points_sf, data_centers_albers) {
  distances <- st_distance(salmon_points_sf, data_centers_albers)
  min_distances <- apply(distances, 1, min)
  return(min_distances)
}

# Function to create distance visualization
create_distance_visualization <- function(salmon_ppp, min_distances, data_center_coords, title) {
  salmon_marked <- salmon_ppp %mark% as.numeric(min_distances)
  
  # Create smooth image of distances
  dist_im <- Smooth(salmon_marked, sigma = 20000, dimyx = c(300, 300), at = "pixels")
  
  # Plot
  plot(dist_im, main = title)
  plot(salmon_ppp, add = TRUE, pch = 16, cex = 0.6)
  points(data_center_coords[, "X"], data_center_coords[, "Y"], 
         col = "red", pch = 4, cex = 1.2)
  
  return(dist_im)
}

# Function to perform quadrat test
perform_quadrat_test <- function(salmon_ppp, nx = 4, ny = 4, plot = TRUE) {
  Q <- quadratcount(salmon_ppp, nx = nx, ny = ny)
  
  if (plot) {
    plot(salmon_ppp, pch = 16, use.marks = FALSE, cols = "#046C9A")
    plot(Q, col = "red", add = TRUE)
    plot(intensity(Q, image = TRUE))
  }
  
  test_result <- quadrat.test(Q)
  return(list(Q = Q, test_result = test_result))
}

# Function to create density plot
create_density_plot <- function(salmon_ppp) {
  lambda_u_hat <- density(salmon_ppp)
  plot(lambda_u_hat)
  plot(salmon_ppp, pch = 16, cex = 0.7, add = TRUE)
  plot(salmon_ppp, pch = 16, cex = 0.3, cols = "white", add = TRUE)
  return(lambda_u_hat)
}

# Function to detect hotspots
detect_hotspots <- function(salmon_ppp, alpha = 0.01) {
  R <- bw.ppl(salmon_ppp)
  LR <- scanLRTS(salmon_ppp, r = R)
  
  plot(LR)
  pvals <- eval.im(pchisq(LR, df = 1, lower.tail = FALSE))
  
  # Filter p-values
  pvals_filtered <- eval.im(ifelse(pvals < alpha, pvals, NA))
  
  # Plot the filtered p-value image
  plot(pvals_filtered, main = paste0("Filtered p-values (p < ", alpha, ")"))
  
  return(list(LR = LR, pvals = pvals, pvals_filtered = pvals_filtered))
}

# Function to analyze relationship with distance to data centers
analyze_distance_relationship <- function(salmon_ppp, data_center_ppp) {
  # Create distance map
  dist_im <- distmap(data_center_ppp)
  
  # Plot distance map
  plot(dist_im, main = "Distance to nearest data center (m)")
  plot(salmon_ppp, add = TRUE, pch = 16, cex = 0.6)
  points(data_center_ppp$x, data_center_ppp$y, col = "red", pch = 4, cex = 1.2)
  plot(Window(salmon_ppp), add = TRUE)
  
  # Estimate intensity as function of distance
  rho_model <- rhohat(salmon_ppp, dist_im, confidence = 0.95)
  plot(rho_model, main = "Relationship between salmon density and distance to data centers")
  
  # Optional: create binned model
  breaks <- seq(0, max(dist_im$v), by = 1000)
  rho_model_binned <- rhohat(salmon_ppp, dist_im, breaks = breaks)
  plot(rho_model_binned, main = "Binned relationship (1km intervals)")
  
  return(list(dist_im = dist_im, rho_model = rho_model, rho_model_binned = rho_model_binned))
}


################################################################################
# 4. Analysis for 2010 Data
################################################################################

analyze_year <- function(year) {
  cat("\n===============================================\n")
  cat("Analyzing salmon data for year", year, "\n")
  cat("===============================================\n")
  
  # Filter data by year
  year_data <- filter_by_year(salmon_data_ca, year)
  cat("Number of observations in", year, ":", nrow(year_data), "\n")
  
  # Transform points to BC Albers and filter to BC only
  transformed_points <- transform_points_to_bc(year_data, bc_window_albers)
  salmon_inside_bc <- transformed_points$bc_points_ppp
  
  cat("Number of observations inside BC:", npoints(salmon_inside_bc), "\n")
  
  # Plot points inside BC
  plot(salmon_inside_bc, main = paste("Salmon points inside BC only (", year, ")"), 
       pch = 16, cex = 0.7)
  
  # Convert to sf for distance calculation
  salmon_points_sf <- st_as_sf(
    data.frame(x = salmon_inside_bc$x, y = salmon_inside_bc$y),
    coords = c("x", "y"), crs = 3005
  )
  
  # Calculate distances to data centers
  min_distances <- calculate_center_distances(salmon_points_sf, data_centers_albers)
  salmon_points_sf$min_center_dist <- min_distances
  
  # Create distance visualization
  dist_im <- create_distance_visualization(
    salmon_inside_bc, min_distances, data_center_coords,
    paste("Distance to nearest data center (", year, ")")
  )
  
  # Point pattern analysis
  cat("\nIntensity (points per unit area):", intensity(salmon_inside_bc), "\n")
  
  # Quadrat test
  quadrat_results <- perform_quadrat_test(salmon_inside_bc)
  cat("Quadrat test p-value:", quadrat_results$test_result$p.value, "\n")
  
  # Density estimation
  lambda_u_hat <- create_density_plot(salmon_inside_bc)
  
  # Hotspot detection
  hotspots <- detect_hotspots(salmon_inside_bc)
  
  # Create data center point pattern
  data_center_ppp <- ppp(
    x = data_center_coords[, "X"],
    y = data_center_coords[, "Y"],
    window = bc_window_albers
  )
  
  # Analyze relationship with distance to data centers
  distance_analysis <- analyze_distance_relationship(salmon_inside_bc, data_center_ppp)
  
  return(list(
    salmon_inside_bc = salmon_inside_bc,
    dist_im = dist_im,
    lambda_u_hat = lambda_u_hat,
    hotspots = hotspots,
    distance_analysis = distance_analysis
  ))
}



################################################################################
# 5. Run Analysis for Multiple Years
################################################################################

# Run analysis for 2010
results_2010 <- analyze_year(2010)

# Run analysis for 2024
results_2024 <- analyze_year(2024)

# Compare intensities
cat("\n===============================================\n")
cat("Comparison of years\n")
cat("===============================================\n")
cat("Intensity 2010:", intensity(results_2010$salmon_inside_bc), "\n")
cat("Intensity 2024:", intensity(results_2024$salmon_inside_bc), "\n")

####
E_parks <- envelope(results_2010$salmon_inside_bc,
                    Kest,
                    correction="border",
                    rank = 1,
                    nsim = 19, # 0.05
)
plot(E_parks,
     main = "",
     lwd = 2)
E_parks_2024 <- envelope(results_2024$salmon_inside_bc,
                         Kest,
                         correction="border",
                         rank = 1,
                         nsim = 19, # 0.05
)
plot(E_parks_2024,
     main = "",
     lwd = 2)
