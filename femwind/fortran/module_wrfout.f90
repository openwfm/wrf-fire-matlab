module module_wrfout

use module_netcdf
use IFPORT    ! intel fortran only

contains

subroutine read_initial_wind(filename,u0_fmw,v0_fmw,w0_fmw,frame0_fmw,frame)
implicit none
!*** purpose
! read wrf data, cycle until frame_fmw matches and chsum is correct
!*** arguments
character(len=*), intent(in)::filename  ! open file
real, pointer, intent(out), dimension(:,:,:)::u0_fmw,v0_fmw,w0_fmw
integer, intent(in)::frame0_fmw   ! frame number to expect
integer, intent(in), optional::frame ! the default frame in the file to read, default=1
! return: 0=OK, >0 timed out
!*** local
integer::istep=1,chsum0,sr(2),chsum0_fmw,ierr=0,maxtry=100,frame0_in,itry,ncid
character(len=256)::msg

!*** executable
if(present(frame))istep=frame 
if(netcdf_msglevel>=0)print *,"reading from file frame ",istep

if(netcdf_msglevel>=0)print *,"expecting time step frame0_fmw=",frame0_fmw
call ncopen(filename,nf90_nowrite,ncid)
do itry=1,maxtry
  frame0_in = netcdf_read_int_wrf(ncid,"FRAME0_FMW",istep)
  if(netcdf_msglevel>=0)print *,"try ",itry," got ",frame0_in
  if(frame0_fmw .eq. frame0_in)goto 1
  call ncclose(ncid)  
  call sleep(1)
  call ncopen(filename,nf90_nowrite,ncid)
enddo
write(msg,*)'timed out after ',maxtry,' tries waiting for frame ',frame0_fmw,' got ',frame0_in
call crash(trim(msg))
1 continue

call get_sr(ncid,sr) ! submesh refinement factors
do itry=1,maxtry
  chsum0_fmw = netcdf_read_int_wrf(ncid,"CHSUM0_FMW",istep)
  if(netcdf_msglevel>=0)print *,"read CHSUM0_FMW=",chsum0_fmw
  
  call netcdf_read_array_wrf(ncid,"U0_FMW",istep,sr,a3d=u0_fmw)
  call netcdf_read_array_wrf(ncid,"V0_FMW",istep,sr,a3d=v0_fmw)
  call netcdf_read_array_wrf(ncid,"W0_FMW",istep,sr,a3d=w0_fmw)
  
  chsum0 = get_chsum(u0_fmw)
  chsum0 = ieor(chsum0,get_chsum(v0_fmw))
  chsum0 = ieor(chsum0,get_chsum(w0_fmw))
  if(netcdf_msglevel>=0)print *," computed chsum0 ", chsum0
  call ncclose(ncid)  
  if (chsum0_fmw.eq.chsum0)goto 2
  call sleep(1)
  call ncopen(filename,nf90_nowrite,ncid)
enddo
call ncclose(ncid)  
write(msg,*)'timed out after ',maxtry,' tries waiting for correct check sum'
call crash(trim(msg))
2 continue

end subroutine read_initial_wind

end module module_wrfout
