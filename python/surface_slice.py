import xarray as xr

file='wrfout.nc'
new ='smoke.nc'

print('Open the original NetCDF file',file)
ds = xr.open_dataset(file)

print('Read 2d slices from the last fram')
pm25 = ds['tr17_1'][-1,0,:,:]
lon  = ds['XLONG'][-1,:,:]
lat  = ds['XLAT'][-1,:,:]

print('Create a new dataset with 2d slices')
new_ds = xr.Dataset({'pm25':pm25, 'lat':lon, 'lat':lat})

# Remove the time dimension from the new dataset
new_ds = new_ds.squeeze(drop=True)
if 'XTIME' in new_ds:
    new_ds = new_ds.drop_vars('XTIME')

print('Create new NetCDF file',new,'and write the dataseet to it')
new_ds.to_netcdf(new)



