import xarray as xr
import numpy as np
import sys, os

pm25_vars = []
times_vars = []
lon_vars = []
lat_vars = []
new = sys.argv[-1]

for file_path in sys.argv[1:-1]:

    print('Reading NetCDF file',file_path)
    ds = xr.open_dataset(file_path)

    # read lowest layer
    pm25_vars.append(ds["tr17_1"].isel(bottom_top=0))

    # Get the other variables, Times, longitude, latitude
    times_vars.append(ds["Times"])
    lon_vars.append(ds["XLONG"])
    lat_vars.append(ds["XLAT"])

print('Concatenating the data variables along the time dimension')
print('pm25')
pm25 = xr.concat(pm25_vars, dim="Time")
print('times')
times = xr.concat(times_vars, dim="Time")
print('lon')
lon = xr.concat(lon_vars, dim="Time")
print('lat')
lat = xr.concat(lat_vars, dim="Time")

del pm25_vars, times_vars, lon_vars, lat_vars 

ds_new = xr.Dataset({'pm25':pm25,'Times':times,'XLONG':lon,'XLAT':lat})

print('Creating new NetCDF file',new)
ds_new.to_netcdf(new)
    
    
    
