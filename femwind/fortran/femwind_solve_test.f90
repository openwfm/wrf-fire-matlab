program femwind_solve_test

use module_femwind_solve
use module_io_matlab

implicit none

real, pointer:: X_m(:,:,:),Y_m(:,:,:),Z_m(:,:,:), &
                u0_m(:,:,:), v0_m(:,:,:), w0_m(:,:,:),   &
                u_m(:,:,:), v_m(:,:,:), w_m(:,:,:)

integer ::                          &
    ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
    ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
    ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts, jfte                          ! fire tile bounds

real, pointer:: A(:,:,:), X(:,:,:),Y(:,:,:),Z(:,:,:), &
                u0(:,:,:), v0(:,:,:), w0(:,:,:),   &
                u(:,:,:), v(:,:,:), w(:,:,:)

integer:: i, j, k, n(3)

call read_array(A,'A')  ! matrices read from Matlab are _m
call read_array(X_m,'X')
call read_array(Y_m,'Y')
call read_array(Z_m,'Z')
call read_array(u0_m, 'u0')
call read_array(v0_m, 'v0')
call read_array(w0_m, 'w0')

n = shape(X)

ifts = 1        ! tile is defined in cells not vertices
ifte = n(1)-1
jfts = 1
jfte = n(3)-1
kfts = 1
kfte = n(2)-1
ifms = ifts-1  ! at least one larger
ifme = ifte+2
jfms = jfts-1
jfme = jfte+2
kfms = kfts-1
kfme = kfte+2
ifds = ifts
ifde = ifte
jfds = jfts
jfde = jfte
kfds = kfts
kfde = kfte


allocate(X(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(Y(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(Z(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(u0(ifms:ifme,kfms:kfme,jfms:jfme)) ! vector components called u v W
allocate(v0(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(w0(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(u(ifms:ifme,kfms:kfme,jfms:jfme)) ! vector components called u v W
allocate(v(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(w(ifms:ifme,kfms:kfme,jfms:jfme))

! copy the input data to tile sized bounds
! X Y Z are corner based, upper bound larger by one
do j=jfts,jfte+1
  do k=kfts,kfte+1
    do i=ifts,ifte+1
        X(i,k,j) = X_m(i,k,j)
        Y(i,k,j) = Y_m(i,k,j)
        Z(i,k,j) = Z_m(i,k,j)
    enddo
  enddo
enddo

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
call femwind_solve(  &
  ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
  ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
  ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire patch bounds
  ifts, ifte, kfts, kfte, jfts,jfte,                        & ! fire tile bounds
  A, X , Y, Z, u0, v0, w0,                                  & ! input arrays
  u, v, w )                                                   ! output arrays

! write output as is in 3D but with tight dimensions
call write_array(u(ifts:ifte,kfts:kfte,jfts:jfte),'u')  
call write_array(v(ifts:ifte,kfts:kfte,jfts:jfte),'v')  
call write_array(w(ifts:ifte,kfts:kfte,jfts:jfte),'w')  

end program femwind_solve_test
