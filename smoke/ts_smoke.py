import netCDF4 as nc4
import numpy as np
from numpy import concatenate as cat
import sys, pickle

times = []
pm25_vars = []
lon_vars = []
lat_vars = []
new = sys.argv[-1]

# file_path = "wrfout.nc"
for file_path in sys.argv[1:-1]:

    print('Reading NetCDF file',file_path)
    d = nc4.Dataset(file_path,'r')
    # extract ESMF string times
    frames = [''.join(x) for x in d.variables['Times'][:].astype(str)]
    for i in frames:
       print(i)
    times     += frames 
    pm25_vars += [d.variables["tr17_1"][:,0,:,:]]
    lon_vars  += [d.variables["XLONG"][:,:,:]]
    lat_vars  += [d.variables["XLAT"][:,:,:]]

pm25 = cat(pm25_vars,0); del pm25_vars
lon  = cat(lon_vars,0); del lon_vars
lat  = cat(lat_vars,0); del lat_vars

