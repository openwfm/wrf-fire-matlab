program netcdf_test

use module_utils
use module_netcdf
use module_wrfout

!*** local
character(len=128)::filename
integer::ncid,chsum0,frame,chsum0_fmw
real, pointer, dimension(:,:,:)::u0_fmw,v0_fmw,w0_fmw
integer::sr(2)

!*** executable
filename = "wrf.nc"

call ncopen(filename,nf90_nowrite,ncid)

do frame = 1,9999

  if(read_initial_wind(ncid,u0_fmw,v0_fmw,w0_fmw,frame=frame).ne.0)then
    print *,'check sum does not agree'
    stop 1
  endif

enddo

call ncclose(ncid)

end program netcdf_test
