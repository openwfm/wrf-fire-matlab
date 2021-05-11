program sweeps_test

use module_sweeps
use module_utils

implicit none

real, pointer:: Kmat(:,:,:,:), x(:,:,:),F(:,:,:), &  ! fortran is not case sensitive
                Kmat_m(:,:,:,:), x_m(:,:,:), F_m(:,:,:)
real, pointer::a(:)
integer :: s(4),n(3)
real::reldif

integer :: msize, &
    ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
    ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
    ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte                            ! fire tile bounds
integer :: i,j,k,jx

! read input arrays in ijk index ordering and tight bounds
call read_array_nd(a,s,'Kmat')
allocate(Kmat_m(s(1),s(2),s(3),s(4)))
Kmat_m = reshape(a,s)
call read_array_nd(a,n,'Fmat')
allocate(F_m(n(1),n(2),n(3)))
F_m = reshape(a,n)
call read_array_nd(a,n,'x_sweeps')
allocate(x_m(n(1),n(2),n(3)))
x_m = reshape(a,n)

if (s(1).ne.n(1).or.s(2).ne.n(2).or.s(3).ne.n(3))call crash('sweeps_test: inconsistent size Kmat and x')

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
allocate(Kmat(ifms:ifme,kfms:kfme,jfms:jfme,1:msize))
allocate(   F(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(   x(ifms:ifme,kfms:kfme,jfms:jfme))
Kmat = 0.
F = 0.
x = 0.


! copy the input data
do j=jfts,jfte
  do k=kfts,kfte
    do i=ifts,ifte
      do jx = 1,msize
        Kmat(i,k,j,jx) = Kmat_m(i,j,k,jx)
      enddo
      x(i,k,j)=x_m(i,j,k)
      F(i,k,j)=F_m(i,j,k)
    enddo
  enddo
enddo

write(*,'(a)')'calling sweeps'
call sweeps(  &
  ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
  ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
  ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire patch bounds
  ifts, ifte, kfts, kfte, jfts,jfte,                        & ! fire tile bounds
  Kmat, F, x, reldif)

write(*,'(a,3i8)')'copying the output data to array size ',n
allocate(x_m(1:n(1),1:n(2),1:n(3)))
do j=jfts,jfte
  do k=kfts,kfte
    do i=ifts,ifte
      x_m(i,j,k)=x(i,k,j)
    enddo
  enddo
enddo

call write_array_nd(reshape(x_m,(/product(n)/)),n,'x_sweeps')

end program sweeps_test
