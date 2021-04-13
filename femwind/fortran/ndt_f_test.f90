program ndt_f_test

use module_ndt_f_assembly
use module_hexa
use module_io_matlab

implicit none

real, pointer:: F(:), F_m(:,:,:)  ! fortran is not case sensitive

integer :: F_dim, x_dim(3),                                       &
    ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
    ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
    ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte                           ! fire tile bounds
    

real, pointer:: A(:,:,:), X(:,:,:),Y(:,:,:),Z(:,:,:), u0(:,:,:), aflags(:,:,:)    ! fortran is not case sensitive
integer:: iflags

call read_array(A,'A')
call read_array(aflags,'iflags')
call read_array(X,'X')
call read_array(Y,'Y')
call read_array(Z,'Z')

iflags = aflags(1,1,1) 
x_dim = shape(X)

ifts = 1
ifte = x_dim(1)
jfts = 1
jfte = x_dim(2)
kfts = 1
kfte = x_dim(3)
ifms = ifts-1
ifme = ifte+1
jfms = jfts-1
jfme = jfte+1
kfms = kfts-1
kfme = kfte+1

F_dim = x_dim(1)*x_dim(2)*x_dim(3)
allocate(F(1:F_dim))
F = 0.

write(*,'(a)')'calling ndt_f_assembly'
call ndt_f_assembly(  &
  ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
  ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
  ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire patch bounds
  ifts, ifte, kfts, kfte, jfts,jfte,                        & ! fire tile bounds
  A, X ,Y ,Z, iflags, F, F_dim)

F_m = reshape(F,(/F_dim,1,1/))

call write_array_nd(F,(/F_dim,1,1/),'Fvec')

end program ndt_f_test
