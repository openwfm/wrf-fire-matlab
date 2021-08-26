program netcdf_test

use module_utils
use module_netcdf

!*** local
character(len=128)::filename
integer::ncid,chsum0,istep
real, pointer, dimension(:,:,:)::u0,v0,w0
integer::sr(2)

!*** executable
filename = "wrf.nc"

call ncopen(filename,nf90_nowrite,ncid)

chsum0 = netcdf_read_int(ncid,"CHSUM0_FMW")
print *,"CHSUM0=",chsum0

sr=(/10,10/)  ! to strip at i j ends

call netcdf_read_array_wrf(ncid,u0,"U0_FMW",1,sr)
print *,"got array ikj shape ",shape(u0)

call ncclose(ncid)

end program netcdf_test
