library(ncdf4)

# Open the netCDF file
nc <- nc_open("wrfout.nc")

# Get the variable you need
pm25 <- ncvar_get(nc, "tr17_1", start = c(1, 1, 1, nctime), count = c(-1, -1, -1, 1))
lon <- ncvar_get(nc, "XLONG")
lat <- ncvar_get(nc, "XLAT")

# Close the netCDF file
nc_close(nc)

# Create a new netCDF file with the variable you extracted
nc_new <- nc_create("smoke.nc")
ncvar_def(nc_new, "pm25", "double", dim = list("south_north", "west_east", "bottom_top", "time"))
ncvar_put(nc_new, "pm25", pm25)
ncvar_def(nc_new, "lon", "double", dim = list("south_north", "west_east"))
ncvar_put(nc_new, "lon", lon)
ncvar_def(nc_new, "lat", "double", dim = list("south_north", "west_east"))
ncvar_put(nc_new, "lat", lat)
nc_close(nc_new)

