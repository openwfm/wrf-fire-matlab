program netcdf_test

use module_utils
use module_netcdf

!*** local
character(len=128)::filename
integer::ncid,chsum0,istep,chsum0_fmw
real, pointer, dimension(:,:,:)::u0_fmw,v0_fmw,w0_fmw
integer::sr(2)

!*** executable
filename = "wrf.nc"
istep = 2

call ncopen(filename,nf90_nowrite,ncid)

chsum0 = netcdf_read_int_wrf(ncid,"CHSUM0_FMW",istep)
print *,"CHSUM0=",chsum0

sr=(/10,10/)  ! to strip at i j ends

call netcdf_read_array_wrf(ncid,u0_fmw,"U0_FMW",istep,sr)
call netcdf_read_array_wrf(ncid,v0_fmw,"V0_FMW",istep,sr)
call netcdf_read_array_wrf(ncid,w0_fmw,"W0_FMW",istep,sr)

call ncclose(ncid)

chsum0_fmw = get_chsum(u0_fmw)
print *,"chsum0_fmw ", chsum0_fmw
chsum0_fmw = ieor(chsum0_fmw,get_chsum(v0_fmw))
print *,"chsum0_fmw ", chsum0_fmw
chsum0_fmw = ieor(chsum0_fmw,get_chsum(w0_fmw))
print *,"chsum0_fmw ", chsum0_fmw
print *,"bound u0_fmw",lbound(u0_fmw),ubound(u0_fmw)
print *,"bound v0_fmw",lbound(v0_fmw),ubound(v0_fmw)
print *,"bound w0_fmw",lbound(w0_fmw),ubound(w0_fmw)
print *,"corner values",u0_fmw(1,1,1),v0_fmw(1,1,1),w0_fmw(1,1,1)

print *,"CHSUM0 computed=",chsum0_fmw," difference ",chsum0_fmw-chsum0


end program netcdf_test
