program vec_boundary_conditions_test

use module_boundary_conditions   ! testing only
use module_utils ! to read and write matrices as text files from matlab

implicit none

real, pointer:: F(:,:,:), &  ! fortran is not case sensitive
                F_m(:,:,:)
integer :: n(3)

integer ::  &
    ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
    ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
    ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire memory bounds
    ifts, ifte, kfts, kfte, jfts,jfte                            ! fire tile bounds
integer :: i,j,k,jx

! read input arrays in ikj index ordering and tight bounds
call read_array(F_m,'F')
n = shape(F_m)

ifts = 1
ifte = n(1)-1
kfts = 1
kfte = n(2)-1
jfts = 1
jfte = n(3)-1
ifds=ifts
ifde=ifte
kfds=kfts
kfde=kfte
jfds=jfts
jfde=jfte
ifms = ifts-1
ifme = ifte+2
jfms = jfts-1
jfme = jfte+2
kfms = kfts-1
kfme = kfte+2

! allocate a little bigger with zeros in extra areas
allocate(F(ifms:ifme,kfms:kfme,jfms:jfme))

! copy the input data 
F(1:n(1),1:n(2),1:n(3))=F_m

write(*,'(a)')'calling vec_boundary_conditions'
call vec_boundary_conditions(  &
  ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
  ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
  ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire memory bounds
  ifts, ifte, kfts, kfte, jfts,jfte,                        & ! fire tile bounds
  F)

! copy the output data 
F_m=F(1:n(1),1:n(2),1:n(3))

call write_array(F_m,'Fb')

end program vec_boundary_conditions_test
