program femwind_wrfout

use module_femwind
use module_utils
use module_common
use module_netcdf
use module_wrfout

implicit none

integer ::                          &
    ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
    ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
    ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts, jfte                          ! fire tile bounds

! A, msize included from module femwind
real, pointer:: X(:,:,:),Y(:,:,:),Z(:,:,:), &
                u0(:,:,:), v0(:,:,:), w0(:,:,:),   &
                u(:,:,:), v(:,:,:), w(:,:,:),Kmat(:,:,:,:),A_m(:,:,:)

integer:: i, j, k, n(3)
integer:: ncid,frame=1
real:: rate
character(len=128)::filename 

!*** executable
filename = "wrf.nc"
call ncopen(filename,nf90_nowrite,ncid)

if(read_initial_wind(ncid,u0,v0,w0,frame=frame).ne.0)then
    print *,'check sum does not agree'
    stop 1
endif

return


!call read_array(A_m,'A_input')  ! matrices read from Matlab are _m
!call read_array(X_m,'X_input')
!call read_array(Y_m,'Y_input')
!call read_array(Z_m,'Z_input')

!params%A = reshape(A_m,(/3,3/))

n = shape(X)
mg(1)%nx = n(1)
mg(1)%ny = n(3)
mg(1)%nz = n(2)

call get_mg_dims(mg(1), &
    ifds, ifde, kfds,kfde, jfds, jfde,            & ! fire grid dimensions
    ifms, ifme, kfms,kfme, jfms, jfme,            &
    ifps, ifpe, kfps,kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts,kfte, jfts,jfte)

allocate(mg(1)%X(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(mg(1)%Y(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(mg(1)%Z(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(u0(ifms:ifme,kfms:kfme,jfms:jfme)) ! vector components called u v W
allocate(v0(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(w0(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(u(ifms:ifme,kfms:kfme,jfms:jfme)) ! vector components called u v W
allocate(v(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(w(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(Kmat(ifms:ifme,kfms:kfme,jfms:jfme,msize)) ! stifness matrix

! copy the input data to tile sized bounds
! X Y Z are corner based, upper bound larger by one
!do j=jfts,jfte+1
!  do k=kfts,kfte+1
!    do i=ifts,ifte+1
!        mg(1)%X(i,k,j) = X_m(i,k,j)
!        mg(1)%Y(i,k,j) = Y_m(i,k,j)
!        mg(1)%Z(i,k,j) = Z_m(i,k,j)
!    enddo
!  enddo
!enddo

mg(1)%dx = mg(1)%X(2,1,1)-mg(1)%X(1,1,1)
mg(1)%dy = mg(1)%Y(1,1,2)-mg(1)%Y(1,1,1)
allocate(mg(1)%dz(mg(1)%nz-1))
do k=kfds,kfte
    mg(1)%dz(k)=mg(1)%Z(1,k+1,1)-mg(1)%Z(1,k,1)
enddo

write(*,'(a)')'calling femwind_setup'
call femwind_setup(mg)    
write(*,'(a)')'femwind_setup returned OK'

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

! write output as is in 3D but with tight dimensions
call write_array(u(ifts:ifte,kfts:kfte,jfts:jfte),'u')  
call write_array(v(ifts:ifte,kfts:kfte,jfts:jfte),'v')  
call write_array(w(ifts:ifte,kfts:kfte,jfts:jfte),'w')  
call write_scalar(rate,'rate')

end program femwind_wrfout
