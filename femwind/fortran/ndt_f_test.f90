program ndt_f_test

use module_ndt_f_assembly
use module_io_matlab

implicit none

real, pointer:: F(:)  ! fortran is not case sensitive

integer :: F_dim ,                                            &
    ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
    ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
    ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte                           ! fire tile bounds
    

real, pointer:: A(:,:,:), X(:,:,:), u0(:,:,:), aflags(:,:,:)    ! fortran is not case sensitive
integer:: iflags(3,1,1)

call read_array(A,'A')
call read_array(X,'X')
call read_array(u0,'u0')
call read_array(aflags,'iflags')

iflags = reshape(aflags,(/3,1,1/)) 

ifts = 1
ifte = size(X,1)
jfts = 1
jfte = size(X,2)
kfts = 1
kfte = size(X,3)
ifms = ifts-1
ifme = ifte+1
jfms = jfts-1
jfme = jfte+1
kfms = kfts-1
kfme = kfte+1

F_dim = size(X,1)*size(X,2)*size(X,3)
allocate(F(F_dim))
F = 0.

write(*,'(a)')'calling ndt_f_assembly'
call ndt_f_assembly(  &
  ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
  ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
  ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire patch bounds
  ifts, ifte, kfts, kfte, jfts,jfte,                        & ! fire tile bounds
  A, X ,X ,X, u0, iflags, F, F_dim)

call write_array(F,'Fvec')

end program ndt_f_test
