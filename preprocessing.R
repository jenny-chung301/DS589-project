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


# load the dataset
data <- read.csv("data/salmon_0009832-250402121839773.csv", sep="\t")
data <- data[!is.na(data$decimalLongitude) & !is.na(data$decimalLatitude), ]

# get CA data
data_ca <- data[data$countryCode == "CA",]

max(data_ca$year, na.rm=TRUE) # 2025
min(data_ca$year, na.rm=TRUE) # 1800


# Get 2010 data
data_2010 <- data_ca[data_ca$year == 2010,]
data_2010 <- data_2010[!is.na(data_2010$decimalLongitude) & !is.na(data_2010$decimalLatitude), ]
rownames(data_2010) <- NULL

# Create a ppp object
x_2010 <- data_2010$decimalLongitude
y_2010 <- data_2010$decimalLatitude
win_2010 <- owin(xrange = range(x_2010), yrange = range(y_2010))
salmon_ppp_2010 <- ppp(x = x_2010, y = y_2010, window = win_2010)
print(salmon_ppp_2010)
summary(salmon_ppp_2010)
plot(salmon_ppp_2010)


# Get a BC window
# Our dataset: EPSG:4326 coordinates
# BC: BC Albers(EPSG:3005)  coordinates
bc <- bc_bound()
print(st_crs(bc))
points_sf <- st_as_sf(data.frame(x = x_2010, y = y_2010), 
                      coords = c("x", "y"), 
                      crs = 4326)

# Converted salmon locations to BC Albers
points_albers <- st_transform(points_sf, 3005)

# Get coordinates
coords_albers <- st_coordinates(points_albers)
x_albers <- coords_albers[, "X"]
y_albers <- coords_albers[, "Y"]

bc_window_albers <- as.owin(bc)
salmon_ppp_2010_albers <- ppp(x = x_albers, y = y_albers, window = bc_window_albers)

# Plot the salmon inside and outside BC 
plot(salmon_ppp_2010_albers, main="Salmon points in BC Albers projection (2010)")

# Extract points within BC
inside_bc <- which(inside.owin(x_albers, y_albers, w = bc_window_albers))
bc_points_x <- x_albers[inside_bc]
bc_points_y <- y_albers[inside_bc]
salmon_inside_bc <- ppp(x = bc_points_x, y = bc_points_y, window = bc_window_albers)

plot(salmon_inside_bc, main="Salmon points inside BC only (2010)",pch=16, cex=0.7)


##############################
# Load the data center dataset
data_centers <- st_read("data/bc_data_centers.geojson")
print(st_crs(data_centers))

# Change the coords to BC 
data_centers_albers <- st_transform(data_centers, 3005) 

salmon_points_sf <- st_as_sf(
  data.frame(
    x = salmon_inside_bc$x,
    y = salmon_inside_bc$y
  ),
  coords = c("x", "y"),
  crs = 3005
)

# Get distances between salmon and a data center
distances <- st_distance(salmon_points_sf, data_centers_albers)
min_distances <- apply(distances, 1, min)

salmon_points_sf$min_center_dist <- min_distances
salmon_marked <- salmon_ppp_2010_albers %mark% as.numeric(min_distances)

# Visualize the distance
dist_im <- Smooth(salmon_marked, sigma=20000, dimyx=c(300,300), at="pixels")
# dist_im <- distmap(salmon_marked)
model <- ppm(salmon_inside_bc ~ dist_im, covariates=list(dist_im=dist_im))

plot(dist_im, main="Distance to nearest data center")
plot(salmon_inside_bc, add=TRUE, pch=16, cex=0.6)
data_center_coords <- st_coordinates(data_centers_albers)
points(data_center_coords[, "X"], data_center_coords[, "Y"], col = "red", pch = 4, cex = 1.2)




###########
# Analysis#
###########

# Point pattern anlaysis - homogeneous vs. inhomogeneous
npoints(salmon_inside_bc) / area(Window(salmon_inside_bc))
intensity(salmon_inside_bc)

Q4by4 <- quadratcount(salmon_inside_bc, nx = 4, ny = 4)

plot(salmon_inside_bc, pch = 16, use.marks = F, cols ="#046C9A")
plot(Q4by4, col= "red", add= T)

plot(intensity(Q4by4, image= T))

quadrat.test(Q4by4) # 2.2e-16

lambda_u_hat <- density(salmon_inside_bc)
plot(lambda_u_hat)
plot(salmon_inside_bc, pch = 16, cex=0.7, add=T)
plot(salmon_inside_bc, pch = 16, cex=0.3, cols="white", add=T)

#hot spot
R <- bw.ppl(salmon_inside_bc)
LR <- scanLRTS(salmon_inside_bc, r = R)
plot(LR)
pvals <- eval.im(pchisq(LR, df= 1, lower.tail=FALSE))
# Create a new image where only p-values < 0.01 are kept; others become NA
pvals_filtered <- eval.im(ifelse(pvals < 0.01, pvals, NA))

# Plot the filtered p-value image
plot(pvals_filtered, main = "Filtered p-values (p < 0.01)")


# rho
# Crate a im object
data_center_ppp <- ppp(
  x = data_center_coords[, "X"],
  y = data_center_coords[, "Y"],
  window = bc_window_albers
)
dist_im <- distmap(data_center_ppp)

## check the rejection
rejected_points <- attr(data_center_ppp, "rejects")
plot(rejected_points)


plot(dist_im, main="Distance to nearest data center (m)")
plot(salmon_inside_bc, add=TRUE, pch=16, cex=0.6)
points(data_center_coords[, "X"], data_center_coords[, "Y"], col = "red", pch = 4, cex = 1.2)
plot(Window(salmon_inside_bc), add=TRUE)

# rhohat
rho_model <- rhohat(salmon_inside_bc, dist_im, confidence=0.95)
plot(rho_model, main="Relationship between salmon density and distance to data centers")

# the closer you are to a data center, the stronger (higher) the salmon intensity,

# Optional
breaks <- seq(0, max(dist_im$v), by=1000)
rho_model_binned <- rhohat(salmon_inside_bc, dist_im, breaks=breaks)
plot(rho_model_binned, main="Binned relationship (1km intervals)")

## Poisson

Q24by24 <- quadratcount(salmon_inside_bc, nx = 24, ny = 24)
mean(Q24by24)
hist(Q24by24)

R_Poisson <- rpois(n = 24^2, lambda = mean(Q24by24))
hist(R_Poisson)



##########
# For 2024
# Get 2010 data
data_2024 <- data_ca[data_ca$year == 2024,]
data_2024 <- data_2024[!is.na(data_2024$decimalLongitude) & !is.na(data_2024$decimalLatitude), ]
rownames(data_2024) <- NULL

# Create a ppp object
x_2024 <- data_2024$decimalLongitude
y_2024 <- data_2024$decimalLatitude
win_2024 <- owin(xrange = range(x_2024), yrange = range(y_2024))
salmon_ppp_2024 <- ppp(x = x_2024, y = y_2024, window = win_2024)
print(salmon_ppp_2024)
summary(salmon_ppp_2024)
plot(salmon_ppp_2024)


points_sf_2024 <- st_as_sf(data.frame(x = x_2024, y = y_2024), 
                      coords = c("x", "y"), 
                      crs = 4326)

# Converted salmon locations to BC Albers
points_albers_2024 <- st_transform(points_sf_2024, 3005)

# Get coordinates
coords_albers_2024 <- st_coordinates(points_albers_2024)
x_albers_2024 <- coords_albers_2024[, "X"]
y_albers_2024 <- coords_albers_2024[, "Y"]

bc_window_albers <- as.owin(bc)
salmon_ppp_2010_albers <- ppp(x = x_albers_2024, y = y_albers_2024, window = bc_window_albers)

# Extract points within BC
inside_bc <- which(inside.owin(x_albers_2024, y_albers_2024, w = bc_window_albers))
bc_points_x <- x_albers_2024[inside_bc]
bc_points_y <- y_albers_2024[inside_bc]
salmon_inside_bc <- ppp(x = bc_points_x, y = bc_points_y, window = bc_window_albers)

plot(salmon_inside_bc, main="Salmon points inside BC only (2024)",pch=16, cex=0.7)



salmon_points_sf <- st_as_sf(
  data.frame(
    x = salmon_inside_bc$x,
    y = salmon_inside_bc$y
  ),
  coords = c("x", "y"),
  crs = 3005
)

# Get distances between salmon and a data center
distances <- st_distance(salmon_points_sf, data_centers_albers)
min_distances <- apply(distances, 1, min)

salmon_points_sf$min_center_dist <- min_distances
salmon_marked <- salmon_ppp_2010_albers %mark% as.numeric(min_distances)

# Visualize the distance
dist_im <- Smooth(salmon_marked, sigma=20000, dimyx=c(300,300), at="pixels")
# dist_im <- distmap(salmon_marked)
model <- ppm(salmon_inside_bc ~ dist_im, covariates=list(dist_im=dist_im))

plot(dist_im, main="Distance to nearest data center")
plot(salmon_inside_bc, add=TRUE, pch=16, cex=0.6)
data_center_coords <- st_coordinates(data_centers_albers)
points(data_center_coords[, "X"], data_center_coords[, "Y"], col = "red", pch = 4, cex = 1.2)


lambda_u_hat <- density(salmon_inside_bc)
plot(lambda_u_hat)
plot(salmon_inside_bc, pch = 16, cex=0.7, add=T)
plot(salmon_inside_bc, pch = 16, cex=0.3, cols="white", add=T)
