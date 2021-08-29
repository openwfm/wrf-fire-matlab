program netcdf_test

use module_utils
use module_netcdf

!*** local
character(len=128)::filename
integer::ncid,chsum0,istep,chsum1
real, pointer, dimension(:,:,:)::u0,v0,w0
integer::sr(2)

!*** executable
filename = "wrf.nc"
istep = 1

call ncopen(filename,nf90_nowrite,ncid)

chsum0 = netcdf_read_int_wrf(ncid,"CHSUM0_FMW",istep)
print *,"CHSUM0=",chsum0

sr=(/10,10/)  ! to strip at i j ends

call netcdf_read_array_wrf(ncid,u0,"U0_FMW",istep,sr)
call netcdf_read_array_wrf(ncid,v0,"V0_FMW",istep,sr)
call netcdf_read_array_wrf(ncid,w0,"W0_FMW",istep,sr)

call ncclose(ncid)

chsum1 = get_chsum(u0)
chsum1 = ieor(chsum1,get_chsum(v0))
chsum1 = ieor(chsum1,get_chsum(w0))

print *,"CHSUM0 computed=",chsum1," difference ",chsum1-chsum0


end program netcdf_test
