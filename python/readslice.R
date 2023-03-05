# Load required packages
library(ggplot2)
library(raster)
library(ncdf4)

ncfile <- nc_open("smoke.nc")

# Read the variables from the NetCDF file
pm25 <- ncvar_get(ncfile, "pm25")
lon <- ncvar_get(ncfile, "XLONG")
lat <- ncvar_get(ncfile, "XLAT")

# Close the NetCDF file
nc_close(ncfile)

print("Create a data frame with the latitude, longitude, and pm25 values")
df <- data.frame(lat = as.vector(lat), lon = as.vector(lon), pm25 = as.vector(pm25))

print("Create a scatter plot of the pm25 values at the lat/lon points")

dev.new()
plot(df$lon, df$lat, col = df$pm25)


