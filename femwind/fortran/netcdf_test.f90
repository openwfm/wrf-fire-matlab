program netcdf_test

use module_utils
use module_netcdf
use module_wrfout

!*** local
character(len=128)::filename
integer::ncid,frame=1,frame_fmw=1
real, pointer, dimension(:,:,:)::u0_fmw,v0_fmw,w0_fmw

!*** executable
filename = "wrf.nc"

call read_initial_wind(filename,u0_fmw,v0_fmw,w0_fmw,1,frame=1)

end program netcdf_test
