program femwind_test

use module_femwind
use module_utils

implicit none

type(mg_type):: mg(max_levels)  ! the main multigrid structure

real, pointer:: X_m(:,:,:),Y_m(:,:,:),Z_m(:,:,:), &
                u0_m(:,:,:), v0_m(:,:,:), w0_m(:,:,:),   &
                u_m(:,:,:), v_m(:,:,:), w_m(:,:,:)

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
real:: rate

call read_array(A_m,'A')  ! matrices read from Matlab are _m
call read_array(X_m,'X')
call read_array(Y_m,'Y')
call read_array(Z_m,'Z')
call read_array(u0_m, 'u0')
call read_array(v0_m, 'v0')
call read_array(w0_m, 'w0')

A = reshape(A_m,(/3,3/))

n = shape(X_m)
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
do j=jfts,jfte+1
  do k=kfts,kfte+1
    do i=ifts,ifte+1
        mg(1)%X(i,k,j) = X_m(i,k,j)
        mg(1)%Y(i,k,j) = Y_m(i,k,j)
        mg(1)%Z(i,k,j) = Z_m(i,k,j)
    enddo
  enddo
enddo

mg(1)%dx = mg(1)%X(2,1,1)-mg(1)%X(1,1,1)
mg(1)%dy = mg(1)%Y(1,1,2)-mg(1)%Y(1,1,1)
allocate(mg(1)%dz(mg(1)%nz-1))
do k=kfds,kfte
    mg(1)%dz(k)=mg(1)%Z(1,k+1,1)-mg(1)%Z(1,1,1)
enddo

write(*,'(a)')'calling femwind_setup'
call femwind_setup(mg)    
write(*,'(a)')'femwind_setup returned OK'

! u is midpoint based
do j=jfts,jfte
  do k=kfts,kfte
    do i=ifts,ifte
        u0(i,k,j) = u0_m(i,k,j)
        v0(i,k,j) = v0_m(i,k,j)
        w0(i,k,j) = w0_m(i,k,j)
    enddo
  enddo
enddo


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

end program femwind_test
