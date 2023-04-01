# Load required packages
library(ncdf4)  # For reading netCDF files
library(ggplot2)
library(dplyr)

# Open input netCDF file
nc_file <- nc_open("ts_smoke.nc")

# read the variables
lon_t <- ncvar_get(nc_file, "XLONG")
lat_t <- ncvar_get(nc_file, "XLAT")
times <- ncvar_get(nc_file, "Times")
pm25_t <- ncvar_get(nc_file, "pm25")
cat("Data dimensions are",dim(pm25_t),"\n")


# Loop over each time step and visualize pm25
for (i in seq_along(times)) {
  # Select pm25 data for the current time step
  pm25 <- pm25_t[,,i]
  lon <- lon_t[,,i]
  lat <- lat_t[,,i]
  
  # Create a data frame from the lon, lat, and pm25 arrays
  df <- data.frame(
    lon = as.vector(t(matrix(lon, nrow=nrow(lon), ncol=ncol(lon)))),
    lat = as.vector(matrix(lat, nrow=nrow(lat), ncol=ncol(lat))),
    pm25 = as.vector(matrix(pm25, nrow=nrow(pm25), ncol=ncol(pm25)))
  )

  # Create a ggplot2 raster plot with longitude and latitude as x and y, and pm25 as the fill color
  p <- ggplot(df, aes(x=lon, y=lat, fill=pm25)) +
    geom_raster() +
    scale_fill_gradientn(colours = c("blue", "yellow", "red")) +
    labs(x="Longitude", y="Latitude", fill="PM2.5")

  # Add a title that includes the times[i] value
  p + ggtitle(paste("PM2.5 at time ", times[i]))
  
  # force display
  print(p)

  # Wait for user to press enter before plotting the next time step
  readline("Press enter to continue...")
}

# Close netCDF file
nc_close(nc_file)


