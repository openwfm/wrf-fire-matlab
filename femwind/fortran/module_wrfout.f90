module module_wrfout
! *** wrf specific netcdf i/o

use module_netcdf
use IFPORT    ! intel fortran only

contains


subroutine netcdf_read_array_wrf(ncid,name,frame,sr,a1d,a2d,a3d)
    implicit none

!*** arguments
    integer, intent(in)::                ncid ! open netcdf file
    real, pointer, intent(out),optional:: a1d(:)  ! the array pointer; remember to deallocate when done with it
    real, pointer, intent(out),optional:: a2d(:,:)  ! the array pointer; remember to deallocate when done with it
    real, pointer, intent(out),optional:: a3d(:,:,:)  ! the array pointer; remember to deallocate when done with it
    character(LEN=*),intent(in):: name
    integer, intent(in), optional::                frame ! time step in the file
    integer, intent(in), dimension(2), optional::  sr ! strip to remove in i and j dimension

!*** local
    integer d1,d2,d3
    integer,dimension(:),allocatable::star,cnts,dims,ends
    integer::i,j,k,varid,ndims
    real,dimension(:,:),allocatable::a1
    real,dimension(:,:,:),allocatable::a2
    real,dimension(:,:,:,:),allocatable::a3
    integer:: istep=1, srf(2)
    character(len=256) msg

!*** executable
    if(present(frame))istep=frame

    d1 = l2i(present(a1d))
    d2 = l2i(present(a2d))
    d3 = l2i(present(a3d))

    if(d1+d2+d3.ne.1)call crash('netcdf_read_array_wrf: must have exactly one of a1d a2d a3d arguments')
    ndims = d1+2*d2+3*d3

    if(ndims.eq.1)then
        srf = 0
    elseif(present(sr))then
        srf = sr
    else
        call get_wrf_dims(ncid,srf)
    endif

    allocate(star(ndims+1))
    allocate(cnts(ndims+1))
    allocate(dims(ndims+1))
    allocate(ends(ndims+1))

    write(msg,*)"netcdf_read_array_wrf reading variable ",trim(name), &
          " dimension ",ndims," time step ",istep
    call message(msg)

    ! get idx
    call netcdf_var_info(ncid,name,dims,varid,netcdf_msglevel)
    write(msg,*)"got dimensions ",trim(name),dims
    call message(msg)
    if(dims(ndims+1).eq.0)then
       write(msg,*)'netcdf_read_array_wrf: wrong dimension ndims=',ndims, ' 1dims=',dims
       call crash(msg)
    endif
     
    star   = 1
    star(ndims+1) = istep
    ends(1:ndims) = dims(1:ndims)
    ends(ndims+1) = istep
    cnts = ends - star + 1

    write(msg,*)"reading ",trim(name),star," to ",ends
    call message(msg)
 
    ! read from file
    select case (ndims)
    case(1) 
        allocate(a1(star(1):ends(1),star(2):ends(2)))
        allocate(a1d(star(1):ends(1)))
        call check(nf90_get_var(ncid, varid, a1, start = star, count = cnts),"nf90_get_var:"//trim(name))
        a1d = a1(:,istep)
        write(msg,*)"returning 1D array length ",shape(a1d)
    case(2) 
	allocate(a2(star(1):ends(1),star(2):ends(2),star(3):ends(3)))
        call check(nf90_get_var(ncid, varid, a2, start = star, count = cnts),"nf90_get_var:"//trim(name))
        ! postprocessing - strip
        write(msg,*)" stripping ",srf," at ij ends"
        call message(msg)
        allocate(a2d(star(1):ends(1)-srf(1),star(2):ends(2)-srf(2)))
        do j=star(2),ends(2)-srf(2)
            do i=star(1),ends(1)-srf(1)
                a2d(i,j) = a2(i,j,istep)
            enddo
        enddo
        write(msg,*)"returning array ij shape ",shape(a2d)
    case(3) 
        allocate(a3(star(1):ends(1),star(2):ends(2),star(3):ends(3),star(4):ends(4)))
        call check(nf90_get_var(ncid, varid, a3, start = star, count = cnts),"nf90_get_var:"//trim(name))
        ! transpose at -> a and strip
        write(msg,*)" stripping ",srf," at ij ends and transposing to ikj indexing order"
        call message(msg)
        allocate(a3d(star(1):ends(1)-srf(1),star(3):ends(3),star(2):ends(2)-srf(2)))
        do k=star(3),ends(3)
            do j=star(2),ends(2)-srf(2)
                do i=star(1),ends(1)-srf(1)
                    a3d(i,k,j) = a3(i,j,k,istep)
                enddo
            enddo
        enddo
        write(msg,*)"returning array ikj shape ",shape(a3d)
    case default
        call crash('wrong case ndims')
    end select
    call message(msg)
end subroutine netcdf_read_array_wrf

subroutine get_wrf_dims(ncid,sr,dims3d)
    implicit none
!*** arguments
    integer, intent(in)::ncid    ! open wrfout dataset
    integer, intent(out), dimension(2):: sr  ! fire "subgrid" refinement factors in x and y directions
    integer, intent(out), dimension(3), optional:: dims3d ! grid dimensions
!*** local
    integer, dimension(3) :: dims_atm, dims_fire   ! 2D + allow time dimension
    integer, dimension(4) :: dims4d   ! 3D + allow time dimension
    integer::varid, prints
    character(len=256)::msg
!*** executable
    prints = netcdf_msglevel
    dims_atm = 0
    dims_fire = 0
    call netcdf_var_info(ncid,"XLAT",dims_atm,varid,prints=prints)
    call netcdf_var_info(ncid,"FXLAT",dims_fire,varid,prints=prints)
    sr = dims_fire(1:2)/(dims_atm(1:2) + 1)
    write(msg,*)"get_wrf_dims: fire subgrid refinement factors are ",sr
    call message(msg)
    if(present(dims3d))then
        call netcdf_var_info(ncid,"PHB",dims4d,varid,prints=prints)
        write(msg,*)"get_wrf_dims: netcdf dimensions ",dims4d
        call message(msg)
        dims3d = dims4d(1:3)
    endif
end subroutine get_wrf_dims

subroutine write_fire_wind(filename,frame0_fmw,uf,vf,u_fmw,v_fmw,w_fmw,frame)
implicit none
!*** purpose
! write solution 
!*** arguments
character(len=*), intent(in)::filename  ! open file
real, intent(in), dimension(:,:)::uf,vf
real, intent(in), dimension(:,:,:),optional::u_fmw,v_fmw,w_fmw
integer, intent(in)::frame0_fmw   ! frame number 
integer, intent(in), optional::frame ! the default frame in the file
!*** local
integer::istep,chsum,ncid
character(len=256)::msg

!*** executable
istep=1
if(present(frame))istep=frame 
write(msg,*)"writing timestep ",frame0_fmw," to ", trim(filename), " frame ",istep
call message(msg)
call ncopen(filename,nf90_write,ncid)
call netcdf_write_2d(ncid,uf,"UF",istep)
call netcdf_write_2d(ncid,vf,"VF",istep)
!if(present(u_fmw))call netcdf_write_wrf_3d(ncid,u_fmw,"U_FMW",istep) ! diagnostics
!if(present(v_fmw))call netcdf_write_wrf_3d(ncid,v_fmw,"V_FMW",istep)
!if(present(w_fmw))call netcdf_write_wrf_3d(ncid,w_fmw,"W_FMW",istep)
chsum=ieor(get_chsum_2d(uf),get_chsum_2d(vf))
print *,'chsum=',chsum
call netcdf_write_int(ncid,chsum,"CHSUM_FMW")
call netcdf_write_int(ncid,frame0_fmw,"FRAME_FMW")
call ncclose(ncid)  

end subroutine write_fire_wind


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
integer::istep=1,chsum0,sr(2),chsum0_fmw,ierr=0,maxtry=200,frame0_in,itry,ncid
character(len=256)::msg

!*** executable
if(present(frame))istep=frame 
write(msg,*)"reading from file frame ",istep
call message(msg)
write(msg,*)"expecting time step frame0_fmw=",frame0_fmw
call message(msg)
call ncopen(filename,nf90_nowrite,ncid)
do itry=1,maxtry
  frame0_in = netcdf_read_int_wrf(ncid,"FRAME0_FMW",istep)
  write(msg,*)"try ",itry," got ",frame0_in," expecting",frame0_fmw
  call message(msg)
  if(frame0_fmw .eq. frame0_in)goto 1
  call ncclose(ncid)  
  if(frame0_fmw .eq. -99)then
       call message('received stop request frame=-99')
       stop
  endif
  call sleep(1)
  call ncopen(filename,nf90_nowrite,ncid)
enddo
write(msg,*)'timed out after ',maxtry,' tries waiting for frame ',frame0_fmw,' got ',frame0_in
call crash(trim(msg))
1 continue

call get_wrf_dims(ncid,sr) ! submesh refinement factors
do itry=1,maxtry
  chsum0_fmw = netcdf_read_int_wrf(ncid,"CHSUM0_FMW",istep)
  write(msg,*)"read CHSUM0_FMW=",chsum0_fmw
  call message(msg)
  
  call netcdf_read_array_wrf(ncid,"U0_FMW",istep,sr,a3d=u0_fmw)
  call netcdf_read_array_wrf(ncid,"V0_FMW",istep,sr,a3d=v0_fmw)
  call netcdf_read_array_wrf(ncid,"W0_FMW",istep,sr,a3d=w0_fmw)
  
  chsum0 = get_chsum(u0_fmw)
  chsum0 = ieor(chsum0,get_chsum(v0_fmw))
  chsum0 = ieor(chsum0,get_chsum(w0_fmw))
  write(msg,*)" computed chsum0 ", chsum0
  call message(msg)
  call ncclose(ncid)  
  if (chsum0_fmw.eq.chsum0)goto 2
  call sleep(1)
  call ncopen(filename,nf90_nowrite,ncid)
enddo
call ncclose(ncid)  
write(msg,*)'timed out after ',maxtry,' tries waiting for correct check sum'
call crash(trim(msg))
2 continue
write(msg,*)'success check sum match for time step ',frame0_fmw
call message(msg)

end subroutine read_initial_wind

end module module_wrfout
