program femwind_wrfout

use module_femwind
use module_utils
use module_common
use module_netcdf
use module_wrfout

implicit none

integer::write_debug = 0,call_femwind = 0

integer ::                          &
    ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
    ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
    ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts, jfte                          ! fire tile bounds

! declarations of A, msize included from module femwind

! variables read from wrfout
real, pointer:: u0_fmw(:,:,:), v0_fmw(:,:,:), w0_fmw(:,:,:), zsf(:,:), ht_fmw(:)

! variables allocated and computed here
real, pointer:: u0(:,:,:), v0(:,:,:), w0(:,:,:)
real, pointer:: u(:,:,:), v(:,:,:), w(:,:,:), uf(:,:), vf(:,:)

integer:: i, j, k, n(3), nx, ny, nz
integer::ncid,frame,sr(2),frame0_fmw,mframe=100,dims(3)
real:: rate,A(3,3),dx,dy,zx,zy
character(len=256)::filename, msg

!*** executable

filename = "wrf.nc" ! file to read from and write to
frame = 1           ! frame in the file 

!************************************
! read static data 
!************************************

call ncopen(filename,nf90_nowrite,ncid)

! horizontal height at cell centers
call netcdf_read_array_wrf(ncid,"ZSF",frame=frame,a2d=zsf) 

! height of the vertical layers above the terrain
call netcdf_read_array_wrf(ncid,"HT_FMW",frame=frame,a1d=ht_fmw)

!  fire subgrid refinement ratio, 3d grid dimensions cell centered
call get_wrf_dims(ncid,sr,dims)  
! horizontal mesh spacing
dx = netcdf_read_att(ncid,"DX")/sr(1)
dy = netcdf_read_att(ncid,"DY")/sr(2)

call ncclose(ncid)

!************************************
! process static data 
!************************************

!construct mg data 

params%A = reshape((/1.0,0.0,0.0, &
                     0.0,1.0,0.0, &
                     0.0,0.0,1.0 /), &
                   (/3,3/))

! finest level 1
mg(1)%nx = dims(1) + 1  ! dimensions are vertex centered
mg(1)%ny = dims(2) + 1  ! dimensions are vertex centered
mg(1)%nz = dims(3) + 1  ! dimensions are vertex centered

call get_mg_dims(mg(1), &                 !  set bounds compatible with WRF
    ifds, ifde, kfds,kfde, jfds, jfde,            & ! fire grid dimensions
    ifms, ifme, kfms,kfme, jfms, jfme,            &
    ifps, ifpe, kfps,kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts,kfte, jfts, jfte)

! allocations persistent over the whole run

allocate(mg(1)%X(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(mg(1)%Y(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(mg(1)%Z(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(u0(ifms:ifme,kfms:kfme,jfms:jfme)) ! vector components called u v W
allocate(v0(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(w0(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(u(ifms:ifme,kfms:kfme,jfms:jfme)) ! vector components called u v W
allocate(v(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(w(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(uf(ifts:ifte,jfts:jfte))
allocate(vf(ifts:ifte,jfts:jfte))
  write(msg,*)"shape(uf)=",shape(uf)
  call message(msg)
  write(msg,*)"shape(vf)=",shape(vf)
  call message(msg)


! X Y Z are corner based, upper bound larger by one
! first interpolate/extrapolate zsf to corners
do j=jfts,jfte+1
  do i=ifts,ifte+1
    if(i.eq.ifts)then
      zx = zsf(i,j) - 0.5*(zsf(i+1,j)-zsf(i,j))  ! extrapolate down
    elseif(i.eq.ifte+1)then
      zx = zsf(i,j) - 0.5*(zsf(i-1,j)-zsf(i,j))  ! extrapolate up 
    else
      zx = 0.5*(zsf(i,j) + zsf(i+1,j))           ! interpolate
    endif
    if(j.eq.jfts)then
      zy = zsf(i,j) - 0.5*(zsf(i,j+1)-zsf(i,j))  ! extrapolate down
    elseif(j.eq.jfte+1)then
      zy = zsf(i,j) - 0.5*(zsf(i,j-1)-zsf(i,j))  ! extrapolate up 
    else
      zy = 0.5*(zsf(i,j) + zsf(i,j+1))           ! interpolate
    endif
    mg(1)%Z(i,kfts,j) = 0.5*(zx + zy)
  enddo
enddo

! (X,Y) is the same uniform grid on all levels
do j=jfts,jfte+1
  do k=kfts,kfte+1
    do i=ifts,ifte+1
        mg(1)%X(i,k,j) = (i-1)*dx
        mg(1)%Y(i,k,j) = (j-1)*dy
    enddo
  enddo
enddo

! add top of levels above terrain to Z(:,kfts,:)
print *,'debug ht_fmw=',ht_fmw
do j=jfts,jfte+1
  do k=kfts+1,kfte+1
    do i=ifts,ifte+1
        mg(1)%Z(i,k,j) = mg(1)%Z(i,kfts,j) + ht_fmw(k - kfts)
    enddo
  enddo
enddo

! mesh spacing
mg(1)%dx = dx 
mg(1)%dy = dy
allocate(mg(1)%dz(mg(1)%nz-1))
do k=kfds,kfte
    mg(1)%dz(k)=mg(1)%Z(1,k+1,1)-mg(1)%Z(1,k,1)
enddo

deallocate(ht_fmw,zsf)  ! no longer needed

if(call_femwind.gt.0)then
call message('calling femwind_setup')
call femwind_setup(mg)    
call message('femwind_setup done')
endif

if(write_debug.gt.0)then
  ! write coordinates for debug/display, even if the caller knows
  call write_average_to_center(mg(1)%X,'X_c')
  call write_average_to_center(mg(1)%Y,'Y_c')
  call write_average_to_center(mg(1)%Z,'Z_c')
endif

!************************************
! loop over time steps 
!************************************

do frame0_fmw=1,mframe

  ! read initial velocity field, loop until delivered 
  call read_initial_wind(filename,u0_fmw,v0_fmw,w0_fmw,frame0_fmw,frame=1)
  
  ! copy the input data to tile sized bounds
  ! initial wind is at cell centers, indexing was already switched to ikj in reading
  u0(ifts:ifte,kfts:kfte,jfts:jfte) = u0_fmw
  v0(ifts:ifte,kfts:kfte,jfts:jfte) = v0_fmw
  w0(ifts:ifte,kfts:kfte,jfts:jfte) = w0_fmw
  
  ! save memory while solver is running
  deallocate(u0_fmw,v0_fmw,w0_fmw)

  ! compute/update mass consistent flow
  
if(call_femwind.gt.0)then
  write(*,'(a)')'calling femwind_solve'
  call femwind_solve(  mg,&
    ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
    ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
    ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte,                        & ! fire tile bounds
    u0, v0, w0,                                  & ! input arrays
    u, v, w,                                                  & ! output arrays
    rate)
  write(*,'(a)')'femwind_solve returned OK'
else
  ! interpolate to fire height
  call message('TESTING ONLY: copying lowest level of input to uf vf')
  write(msg,*)"shape(uf)=",shape(uf)
  call message(msg)
  write(msg,*)"shape(vf)=",shape(vf)
  call message(msg)
  uf = u0(ifts:ifte,1,jfts:jfte)
  vf = v0(ifts:ifte,1,jfts:jfte)
  call write_fire_wind(filename,uf,vf,frame0_fmw,frame=1)
endif

  if(write_debug.gt.0)then
    ! write output as is in 3D but with tight dimensions
    call write_array(u(ifts:ifte,kfts:kfte,jfts:jfte),'u')  
    call write_array(v(ifts:ifte,kfts:kfte,jfts:jfte),'v')  
    call write_array(w(ifts:ifte,kfts:kfte,jfts:jfte),'w')  
    call write_scalar(rate,'rate')
  endif

enddo

contains

subroutine write_average_to_center(F,name)
! average and return with tight bounds
implicit none
!*** arguments
real, intent(in)::F(ifms:ifme, kfms:kfme, jfms:jfme)
character(len=*)::name
!*** local
real, allocatable::Fm(:,:,:)
integer::ii,kk,jj
real::s
!*** executable
allocate(Fm(ifts:ifte,kfts:kfte,jfts:jfte))
do j=jfts,jfte
  do k=kfts,kfte
    do i=ifts,ifte
      s = 0.
      do ii=0,1
        do kk=0,1
          do jj=0,1
            s = s + F(i+ii,k+kk,j+jj) 
          enddo
        enddo
      enddo
      Fm(i,k,j) = s/8
    enddo
  enddo
enddo
call write_array(Fm,name)
end subroutine write_average_to_center 

end program femwind_wrfout
