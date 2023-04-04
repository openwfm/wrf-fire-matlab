# Load required packages
library(ncdf4)  # For reading netCDF files

print_variable_info <- function(x) {
    cat("Variable:", deparse(substitute(x)),"class", class(x),"type",typeof(x),
    "dimension",dim(x))
    if(is.numeric(x)) {
      cat(" min",min(x),"max",max(x),"\n")
    } else {
      # The array contains non-numeric values
      cat("\n")
    }
}

# Open input netCDF file
nc_file <- nc_open("ts_smoke.nc")

# read the variables
lon <- ncvar_get(nc_file, "XLONG")
lat <- ncvar_get(nc_file, "XLAT")
times <- ncvar_get(nc_file, "Times")
pm25 <- ncvar_get(nc_file, "pm25")

print_variable_info(lon)
print_variable_info(lat)
print_variable_info(pm25)
print_variable_info(times)

nc_close(nc_file)


