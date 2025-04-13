library(spatstat)
library(bcmaps)
library(sf)

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
model <- ppm(salmon_inside_bc ~ dist_im, covariates=list(dist_im=dist_im))

plot(dist_im, main="Distance to nearest data center")
plot(salmon_inside_bc, add=TRUE, pch=16, cex=0.6)
data_center_coords <- st_coordinates(data_centers_albers)
points(data_center_coords[, "X"], data_center_coords[, "Y"], col = "red", pch = 4, cex = 1.2)

