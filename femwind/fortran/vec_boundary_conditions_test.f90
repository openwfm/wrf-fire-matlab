program vec_boundary_conditions_test

use module_boundary_conditions   ! testing only
use module_io_matlab ! to read and write matrices as text files from matlab

implicit none

real, pointer:: lambda(:,:,:), &  ! fortran is not case sensitive
                lambda_m(:,:,:)
integer :: n(3)

integer ::  &
    ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
    ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
    ifts, ifte, kfts, kfte, jfts,jfte                            ! fire tile bounds
integer :: i,j,k,jx

! read input arrays in ikj index ordering and tight bounds
call read_array(lambda_m,'lambda')
n = shape(lambda_m)

ifts = 1
ifte = n(1)
kfts = 1
kfte = n(2)
jfts = 1
jfte = n(3)
ifds=ifts
ifde=ifte
kfds=kfts
kfde=kfte
jfds=jfts
jfde=jfte
ifms = ifts-1
ifme = ifte+1
jfms = jfts-1
jfme = jfte+1
kfms = kfts-1
kfme = kfte+1

! allocate a little bigger with zeros in extra areas
allocate(lambda(ifms:ifme,kfms:kfme,jfms:jfme))

! copy the input data 
lambda(ifts:ifte,kfts:kfte,jfts:jfte)=lambda_m

write(*,'(a)')'calling vec_boundary_conditions'
call vec_boundary_conditions(  &
  ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
  ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
  ifts, ifte, kfts, kfte, jfts,jfte,                        & ! fire tile bounds
  lambda)

! copy the output data 
lambda_m=lambda(ifts:ifte,kfts:kfte,jfts:jfte)

call write_array(lambda_m,'lambdab')

end program ndt_boundary_conditions_test
