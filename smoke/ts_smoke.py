import netCDF4 as nc4
import numpy as np
from numpy import concatenate as cat
import sys, pickle

times_vars = []
pm25_vars = []
lon_vars = []
lat_vars = []
argv = sys.argv
argv=['wrfout.nc','ts_smoke.nc']
new_path = argv[-1]

for file_path in argv[0:-1]:
    print('Reading NetCDF file',file_path)
    d = nc4.Dataset(file_path,'r')
    # extract ESMF string times
    frames = [''.join(x) for x in d.variables['Times'][:].astype(str)]
    for i in frames:
       print(i)
    times_vars += [d.variables["Times"][:,:]]
    pm25_vars  += [d.variables["tr17_1"][:,0,:,:]]
    lon_vars   += [d.variables["XLONG"][:,:,:]]
    lat_vars   += [d.variables["XLAT"][:,:,:]]

# take the last metadata; should check if same except dim_size_time[0]
dtype_times=d.variables["Times"].dtype
dtype_vars=d.variables["XLONG"].dtype
dim_size_time = d.variables["Times"][:,:].shape
dim_name_vars = d.variables["XLONG"].dimensions
dim_size_vars = d.variables["XLONG"][:,:,:].shape

# concatenate variables 
times = cat(times_vars,0); del times_vars
pm25 = cat(pm25_vars,0); del pm25_vars
lon  = cat(lon_vars,0); del lon_vars
lat  = cat(lat_vars,0); del lat_vars

print('Writing NetCDF file',new_path)
new = nc4.Dataset(new_path, mode="w")

# create dimensions
new.createDimension('Time',None)
for dim_name, dim in d.dimensions.items():
    if dim_name in ['south_north', 'west_east','DateStrLen']:
        new.createDimension(dim_name, len(dim))
time_dim = new.dimensions['Time']
tstr_dim = new.dimensions['DateStrLen']
sn_dim = new.dimensions['south_north']
we_dim = new.dimensions['west_east']
var_dimensions = (time_dim, sn_dim, we_dim)

# create variables
times_var = new.createVariable('Times', dtype_times, (time_dim, tstr_dim))
lon_var = new.createVariable('XLONG',   dtype_vars, var_dimensions)
lat_var = new.createVariable('XLAT',    dtype_vars, var_dimensions)
pm25_var = new.createVariable('pm25',   dtype_vars, var_dimensions)

# copy values 
times_var[:]=times
lon_var[:]=lon
lat_var[:]=lat
pm25[:]=pm25

new.close()

