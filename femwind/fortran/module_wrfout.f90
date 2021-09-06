module module_wrfout

use module_netcdf

contains

integer function read_initial_wind(ncid,u0_fmw,v0_fmw,w0_fmw,frame)
implicit none
!*** arguments
integer, intent(in)::ncid  ! open file
real, pointer, intent(out), dimension(:,:,:)::u0_fmw,v0_fmw,w0_fmw
integer, intent(in), optional::frame
! return: 0=OK, 1=check sum does not match
!*** local
integer::istep=1,chsum0,sr(2),chsum0_fmw,ierr=0

!*** executable
if(present(frame))istep=frame 
if(netcdf_msglevel>=0)print *,"step ",istep

chsum0 = netcdf_read_int_wrf(ncid,"CHSUM0_FMW",istep)
if(netcdf_msglevel>=0)print *,"CHSUM0=",chsum0

!sr=(/10,10/)  ! to strip at i j ends
call get_sr(ncid,sr)

call netcdf_read_array_wrf(ncid,"U0_FMW",istep,sr,a3d=u0_fmw)
call netcdf_read_array_wrf(ncid,"V0_FMW",istep,sr,a3d=v0_fmw)
call netcdf_read_array_wrf(ncid,"W0_FMW",istep,sr,a3d=w0_fmw)

chsum0_fmw = get_chsum(u0_fmw)
if(netcdf_msglevel>=0)print *,"chsum0_fmw ", chsum0_fmw
chsum0_fmw = ieor(chsum0_fmw,get_chsum(v0_fmw))
if(netcdf_msglevel>=0)print *,"chsum0_fmw ", chsum0_fmw
chsum0_fmw = ieor(chsum0_fmw,get_chsum(w0_fmw))
if(netcdf_msglevel>=0)print *,"chsum0_fmw ", chsum0_fmw
!print *,"bound u0_fmw",lbound(u0_fmw),ubound(u0_fmw)
!print *,"bound v0_fmw",lbound(v0_fmw),ubound(v0_fmw)
!print *,"bound w0_fmw",lbound(w0_fmw),ubound(w0_fmw)
!print *,"corner values",u0_fmw(1,1,1),v0_fmw(1,1,1),w0_fmw(1,1,1)

if(netcdf_msglevel>=0)print *,"CHSUM0 computed=",chsum0_fmw," difference ",chsum0_fmw-chsum0
if (chsum0_fmw.ne.chsum0)then
   if(netcdf_msglevel>=0)print *,'chsum does not agree'
   ierr = 1 
endif

read_initial_wind = ierr

end function read_initial_wind

end module module_wrfout
