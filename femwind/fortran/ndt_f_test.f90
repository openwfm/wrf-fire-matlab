program ndt_f_test

use module_ndt_f_assembly
use module_hexa
use module_io_matlab

implicit none

real, pointer:: F(:),Xmat(:,:,:),Ymat(:,:,:),Zmat(:,:,:), Xu0mat(:,:,:), Yu0mat(:,:,:), Zu0mat(:,:,:)  ! fortran is not case sensitive

integer :: F_dim, x_dim(3),u_dim(3),                          &
    ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
    ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
    ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts, jfte,                       &  ! fire tile bounds
    iuds, iude, kuds, kude, juds, jude,                       &
    iums, iume, kums, kume, jums, jume

real, pointer:: A(:,:,:), X(:,:,:),Y(:,:,:),Z(:,:,:), Xu0(:,:,:), Yu0(:,:,:), Zu0(:,:,:), aflags(:,:,:)    ! fortran is not case sensitive
integer:: iflags, i, j, k

call read_array(A,'A')
call read_array(aflags,'iflags')
call read_array(X,'X')
call read_array(Y,'Y')
call read_array(Z,'Z')
call read_array(Xu0, 'Xu0')
call read_array(Yu0, 'Yu0')
call read_array(Zu0, 'Zu0')

iflags = aflags(1,1,1) 
x_dim = shape(X)

ifts = 1
ifte = x_dim(1)
jfts = 1
jfte = x_dim(3)
kfts = 1
kfte = x_dim(2)
ifms = ifts
ifme = ifte
jfms = jfts
jfme = jfte
kfms = kfts
kfme = kfte

iuds = 1
iude = u_dim(1)
juds = 1
jude = u_dim(3)
kuds = 1
kude = u_dim(2)
iums = iuds
iume = iude
jums = juds
jume = jude
kums = kuds
kume = kude


allocate(Xmat(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(Ymat(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(Zmat(ifms:ifme,kfms:kfme,jfms:jfme))

! copy the input data to tile sized bounds
do j=jfts,jfte
  do k=kfts,kfte
    do i=ifts,ifte
        Xmat(i,k,j) = X(i,k,j)
        Ymat(i,k,j) = Y(i,k,j)
        Zmat(i,k,j) = Z(i,k,j)
    enddo
  enddo
enddo





F_dim = x_dim(1)*x_dim(2)*x_dim(3)
allocate(F(1:F_dim))
F = 0.

write(*,'(a)')'calling ndt_f_assembly'
call ndt_f_assembly(  &
  ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
  ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
  ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire patch bounds
  ifts, ifte, kfts, kfte, jfts,jfte,                        & ! fire tile bounds
  A, Xmat , Ymat, Zmat, Xu0, Yu0, Zu0, iflags, F, F_dim)

!F_m = reshape(F,(/F_dim,1,1/))

call write_array_nd(F,(/F_dim,1,1/),'Fvec')

end program ndt_f_test
