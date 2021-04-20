program ndt_f_test

use module_ndt_f_assembly
use module_hexa
use module_io_matlab

implicit none

real, pointer:: F(:,:,:),Xmat(:,:,:),Ymat(:,:,:),Zmat(:,:,:), &
                Xu0mat(:,:,:), Yu0mat(:,:,:), Zu0mat(:,:,:)  ! fortran is not case sensitive

integer :: F_dim, x_dim(3),u_dim(3),                          &
    ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
    ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
    ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts, jfte,                       &  ! fire tile bounds
    iuds, iude, kuds, kude, juds, jude,                       &
    iums, iume, kums, kume, jums, jume

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
allocate(Xu0mat(iums:iume,kums:kume,jums:jume))
allocate(Yu0mat(iums:iume,kums:kume,jums:jume))
allocate(Zu0mat(iums:iume,kums:kume,jums:jume))
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

do j=juds,jude
  do k=kuds,kude
    do i=iuds,iude
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
call ndt_f_assembly(  &
  ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
  ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
  ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire patch bounds
  ifts, ifte, kfts, kfte, jfts,jfte,                        & ! fire tile bounds
  iums, iume, kums, kume, jums, jume,                       &
  A, Xmat , Ymat, Zmat, Xu0mat, Yu0mat, Zu0mat,             &
  F)

!F_m = reshape(F,(/F_dim,1,1/))
allocate(Fm(ifts:ifte,kfts:kfte,jfts:jfte))
do j=jfts,jfte
  do k=kfts,kfte
    do i=ifts,ifte
        Fm(i,k,j)=F(i,k,j)
    enddo
  enddo
enddo

fsize = (/ifte-ifts+1,kfte-kfts+1,jfte-jfts+1/)

call write_array_nd(reshape(Fm,(/product(fsize)/)),fsize,'Fvec')

end program ndt_f_test
