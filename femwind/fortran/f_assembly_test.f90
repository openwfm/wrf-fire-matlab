program f_assembly_test

use module_f_assembly
use module_utils

implicit none

real, pointer:: F(:,:,:),Xmat(:,:,:),Ymat(:,:,:),Zmat(:,:,:), &
                Xu0mat(:,:,:), Yu0mat(:,:,:), Zu0mat(:,:,:)  ! fortran is not case sensitive

integer :: F_dim, x_dim(3),u_dim(3),                          &
    ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
    ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
    ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts, jfte                          ! fire tile bounds

real, pointer:: A(:,:,:), X(:,:,:),Y(:,:,:),Z(:,:,:), &
                Xu0(:,:,:), Yu0(:,:,:), Zu0(:,:,:),   &
                Fm(:,:,:)    ! fortran is not case sensitive
integer:: i, j, k
integer:: fsize(3)

call read_array(A,'A')
call read_array(X,'X')
call read_array(Y,'Y')
call read_array(Z,'Z')
call read_array(Xu0, 'Xu0')
call read_array(Yu0, 'Yu0')
call read_array(Zu0, 'Zu0')

x_dim = shape(X)
u_dim = shape(Xu0)

ifts = 1
ifte = x_dim(1)-1
jfts = 1
jfte = x_dim(3)-1
kfts = 1
kfte = x_dim(2)-1
ifms = ifts-1
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


allocate(Xmat(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(Ymat(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(Zmat(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(Xu0mat(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(Yu0mat(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(Zu0mat(ifms:ifme,kfms:kfme,jfms:jfme))

! copy the input data to tile sized bounds
! X Y Z are corner based, upper bound larger by one
do j=jfts,jfte+1
  do k=kfts,kfte+1
    do i=ifts,ifte+1
        Xmat(i,k,j) = X(i,k,j)
        Ymat(i,k,j) = Y(i,k,j)
        Zmat(i,k,j) = Z(i,k,j)
    enddo
  enddo
enddo

! u is midpoint based
do j=jfts,jfte
  do k=kfts,kfte
    do i=ifts,ifte
        Xu0mat(i,k,j) = Xu0(i,k,j)
        Yu0mat(i,k,j) = Yu0(i,k,j)
        Zu0mat(i,k,j) = Zu0(i,k,j)
    enddo
  enddo
enddo




F_dim = x_dim(1)*x_dim(2)*x_dim(3)
allocate(F(ifms:ifme, kfms:kfme, jfms:jfme))
F = 0.

write(*,'(a)')'calling ndt_f_assembly'
call f_assembly(  &
  ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
  ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
  ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire patch bounds
  ifts, ifte, kfts, kfte, jfts,jfte,                        & ! fire tile bounds
  A, Xmat , Ymat, Zmat, Xu0mat, Yu0mat, Zu0mat,             &
  F)

! write output as is in 3D but with tight dimensions
call write_array(F(ifts:ifte+1,kfts:kfte+1,jfts:jfte+1),'F')  

end program f_assembly_test
