program netcdf_test

use module_utils
use module_netcdf

!*** local
character(len=128)::filename
integer::ncid


!*** executable
filename = "wrf.nc"

call ncopen(filename,nf90_nowrite,ncid)

call ncclose(ncid)

end program netcdf_test
