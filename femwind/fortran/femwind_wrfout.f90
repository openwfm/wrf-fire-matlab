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
real, pointer:: u0_fmw(:,:,:), v0_fmw(:,:,:), w0_fmw(:,:,:), zsfe(:,:), ht_fmw(:), zsf_fmw(:,:), fwh_fmw(:,:)

! variables allocated and computed here
real, pointer:: u0(:,:,:), v0(:,:,:), w0(:,:,:), hth(:)
real, pointer:: u(:,:,:), v(:,:,:), w(:,:,:), uf(:,:), vf(:,:), wh(:,:)
integer, pointer:: kh(:,:)

integer:: i, j, k, n(3), nx, ny, nz, nfx, nfy, nfz
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

!  fire subgrid refinement ratio, 3d grid dimensions cell centered
call get_wrf_dims(ncid,sr,dims)  
! horizontal mesh spacing
dx = netcdf_read_att(ncid,"DX")/sr(1)
dy = netcdf_read_att(ncid,"DY")/sr(2)

! height of the vertical layers above the terrain
call netcdf_read_array_wrf(ncid,"HT_FMW",frame=frame,a1d=ht_fmw)

! fire mesh size in cells
nfx = dims(1)*sr(1)  
nfy = dims(2)*sr(2)  
nfz = size(ht_fmw,1)

! terrain height at cell centers
call netcdf_read_array_wrf(ncid,"ZSF",frame=frame,a2d=zsf_fmw) 

! fire wind height to interpolate to at cell centers
call netcdf_read_array_wrf(ncid,"FWH",frame=frame,a2d=fwh_fmw) 

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
mg(1)%nx = nfx + 1  ! dimensions are vertex centered
mg(1)%ny = nfy + 1  ! dimensions are vertex centered
mg(1)%nz = nfz + 1  ! dimensions are vertex centered

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
allocate(kh(ifts:ifte,jfts:jfte))
allocate(wh(ifts:ifte,jfts:jfte))
allocate(zsfe(ifts-1:ifte+1,jfts-1:jfte+1))

  write(msg,*)"shape(uf)=",shape(uf)
  call message(msg)
  write(msg,*)"shape(vf)=",shape(vf)
  call message(msg)

! X Y Z are corner based, upper bound larger by one

! copy interior
zsfe(ifts:ifte,jfts:jfte)=zsf_fmw
! extrapolate on sides
zsfe(ifts-1,jfts:jfte)=1.5*zsfe(ifts,jfts:jfte)-0.5*zsfe(ifts+1,jfts:jfte)
zsfe(ifte+1,jfts:jfte)=1.5*zsfe(ifte,jfts:jfte)-0.5*zsfe(ifte-1,jfts:jfte)
zsfe(ifts:ifte,jfts-1)=1.5*zsfe(ifts:ifte,jfts)-0.5*zsfe(ifts:ifte,jfts+1)
zsfe(ifts:ifte,jfte+1)=1.5*zsfe(ifts:ifte,jfte)-0.5*zsfe(ifts:ifte,jfte-1)
! extrapolate at domain corners
zsfe(ifts-1,jfts-1)=1.5*zsfe(ifts,jfts)-0.5*zsfe(ifts+1,jfts+1)
zsfe(ifte+1,jfts-1)=1.5*zsfe(ifte,jfts)-0.5*zsfe(ifte-1,jfts+1)
zsfe(ifts-1,jfte+1)=1.5*zsfe(ifts,jfte)-0.5*zsfe(ifts+1,jfte-1)
zsfe(ifte+1,jfte+1)=1.5*zsfe(ifte,jfte)-0.5*zsfe(ifte-1,jfte-1)
! interpolate to cell corners
mg(1)%Z(ifts:ifte+1,kfts,jfts:jfte+1) = 0.25 * (    &
    zsfe(ifts:ifte+1,jfts:jfte+1) +         &
    zsfe(ifts-1:ifte,jfts:jfte+1) +         &
    zsfe(ifts:ifte+1,jfts-1:jfte) +         &
    zsfe(ifts-1:ifte,jfts-1:jfte) )

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

! heights of wind nodes in cell centers
allocate(hth(0:kfte))
hth(0)=0
hth(1)=ht_fmw(1)*0.5
do k=2,kfte
  hth(k)=0.5*(ht_fmw(k-1)+ht_fmw(k))
enddo

! indices to interpolate to fwh_fmw
do j=jfts,jfte
  do i=ifts,ifte
    do k=1,kfte
      if (hth(k) .ge. fwh_fmw(i,j))then
        kh(i,j)=k
        wh(i,j) = (hth(k) - fwh_fmw(i,j))/(hth(k) - hth(k-1)) ! weight 0 if at the top
        ! write(msg,*)'i,j=',i,j,'k=',k,'hth=',hth(k),hth(k-1),'fwh=',fwh_fmw(i,j),'wh=',wh(i,j)
        ! call message(msg)
        goto 1
      endif
    enddo
    call crash('did not find interpolation layer') 
1   continue 
  enddo
enddo

deallocate(hth,ht_fmw,zsfe)  ! no longer needed

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
    u=u0
    v=v0
    w=w0
    write(*,'(a)')'femwind_solve skipped'
  endif

  ! interpolate u,v to fire height
  do j=jfts,jfte
    do i=ifts,ifte
       k=kh(i,j)
       if(k>1)then
         uf(i,j) = wh(i,j) * u(i,k,j) + (1.0 - wh(i,j) * u(i,k-1,j))
         vf(i,j) = wh(i,j) * v(i,k,j) + (1.0 - wh(i,j) * v(i,k-1,j))
       else
         uf(i,j) = u(i,k,j)
         vf(i,j) = v(i,k,j)
       endif
    enddo
  enddo
  call message('TESTING ONLY: copying lowest level of input to uf vf')
  write(msg,*)"shape(uf)=",shape(uf)
  call message(msg)
  write(msg,*)"shape(vf)=",shape(vf)
  call message(msg)
  uf = u0(ifts:ifte,1,jfts:jfte)
  vf = v0(ifts:ifte,1,jfts:jfte)
  call write_fire_wind(filename,uf,vf,frame0_fmw,frame=1)

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
