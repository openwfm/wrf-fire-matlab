program ndt_boundary_conditions_test

use module_boundary_conditions   ! testing only
use module_io_matlab ! to read and write matrices as text files from matlab

implicit none

real, pointer:: kmat(:,:,:,:), &  ! fortran is not case sensitive
                kmat_m(:,:,:,:),a(:)
integer :: s(4),n(3)

integer :: msize, &
    ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
    ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
    ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte                            ! fire tile bounds
integer :: i,j,k,jx

! read input arrays in ikj index ordering and tight bounds
call read_array_nd(a,s,'kmat')
allocate(kmat_m(s(1),s(2),s(3),s(4)))
kmat_m = reshape(a,s)
n = s(1:3)

ifts = 1
ifte = n(1)
jfts = 1
jfte = n(2)
kfts = 1
kfte = n(3)
msize = s(4)
if(msize.ne.14)call crash('msize must be 14')
ifms = ifts-1
ifme = ifte+1
jfms = jfts-1
jfme = jfte+1
kfms = kfts-1
kfme = kfte+1

! allocate a little bigger with zeros in extra areas
allocate(kmat(ifms:ifme,kfms:kfme,jfms:jfme,1:msize))
kmat = 0.

! copy the input data 
do j=jfts,jfte
  do k=kfts,kfte
    do i=ifts,ifte
      do jx = 1,msize
        kmat(i,k,j,jx) = kmat_m(i,j,k,jx)
      enddo
    enddo
  enddo
enddo
           
write(*,'(a)')'calling ntd_boundary_conditions'
call ndt_boundary_conditions(  &
  ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
  ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
  ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire patch bounds
  ifts, ifte, kfts, kfte, jfts,jfte,                        & ! fire tile bounds
  kmat)

! copy the output data 
do j=jfts,jfte
  do k=kfts,kfte
    do i=ifts,ifte
      do jx = 1,msize
        kmat_m(i,k,j,jx) = kmat(i,j,k,jx)
      enddo
    enddo
  enddo
enddo

call write_array_nd(reshape(kmat_m,(/product(s)/)),s,'kmat')

end program ndt_boundary_conditions_test
