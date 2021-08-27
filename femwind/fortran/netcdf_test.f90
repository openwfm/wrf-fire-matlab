program netcdf_test

use module_utils
use module_netcdf

!*** local
character(len=128)::filename
integer::ncid,chsum0,istep,chsum0_fmw
real, pointer, dimension(:,:,:)::u0,v0,w0
integer::sr(2)

!*** executable
filename = "wrf.nc"

call ncopen(filename,nf90_nowrite,ncid)

chsum0 = netcdf_read_int(ncid,"CHSUM0_FMW")
print *,"CHSUM0=",chsum0

sr=(/10,10/)  ! to strip at i j ends

call netcdf_read_array_wrf(ncid,u0,"U0_FMW",1,sr)
call netcdf_read_array_wrf(ncid,v0,"V0_FMW",1,sr)
call netcdf_read_array_wrf(ncid,w0,"W0_FMW",1,sr)

call ncclose(ncid)

chsum0_fmw = get_chsum(u0)
print *,"chsum0_fmw ", chsum0_fmw
chsum0_fmw = ieor(chsum0_fmw,get_chsum(v0))
print *,"chsum0_fmw ", chsum0_fmw
chsum0_fmw = ieor(chsum0_fmw,get_chsum(w0))
print *,"chsum0_fmw ", chsum0_fmw

print *,"CHSUM0 computed=",chsum0_fmw," difference ",chsum0_fmw-chsum0


end program netcdf_test
