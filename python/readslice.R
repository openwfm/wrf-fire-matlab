# Load required packages
library(ggplot2)
library(raster)

# Read in the NetCDF file
ds <- brick("smoke.nc")

# Extract the latitude, longitude, and PM2.5 data
lat <- ds$lat
lon <- ds$lon
pm25 <- ds$pm25

# Convert the PM2.5 data to a raster object
pm25_raster <- raster(t(pm25))
extent(pm25_raster) <- extent(lon[1], lon[length(lon)], lat[length(lat)], lat[1])

# Plot the raster object using ggplot2
ggplot() + 
  geom_raster(data = as.data.frame(pm25_raster), aes(x = lon, y = lat, fill = value)) + 
  scale_fill_gradient(low = "white", high = "red") + 
  coord_fixed() + 
  theme_bw()

