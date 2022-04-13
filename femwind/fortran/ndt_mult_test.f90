program ndt_mult_test

use module_ndt_mult   ! testing only
use module_utils ! to read and write matrices as text files from matlab

implicit none

real, pointer:: kmat(:,:,:,:), u(:,:,:), y(:,:,:), &  ! fortran is not case sensitive
                kmat_m(:,:,:,:), u_m(:,:,:), y_m(:,:,:),  r(:,:,:)
real, pointer::a(:)
integer :: s(4),n(3)

integer :: msize, &
    ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
    ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
    ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte                            ! fire tile bounds
integer :: i,j,k,jx
real:: siz,relres

! read input arrays in ikj index ordering and tight bounds
call read_array_nd(a,s,'kmat')
allocate(kmat_m(s(1),s(2),s(3),s(4)))
kmat_m = reshape(a,s)
call read_array_nd(a,n,'u')
allocate(u_m(n(1),n(2),n(3)))
u_m = reshape(a,n)

if (s(1).ne.n(1).or.s(2).ne.n(2).or.s(3).ne.n(3))call crash('ndt_mult_test: inconsistent size kmat and u')

ifts = 1
ifte = n(1) - 1
kfts = 1
kfte = n(2) - 1
jfts = 1
jfte = n(3) - 1

ifds = ifts
ifde = ifte
jfds = jfts
jfde = jfte
kfds = kfts
kfde = kfte

msize = s(4)
if(msize.ne.14)call crash('msize must be 14')
ifms = ifts-1
ifme = ifte+2
kfms = kfts-1
kfme = kfte+2
jfms = jfts-1
jfme = jfte+2

! allocate a little bigger with zeros in extra areas
allocate(Kmat(ifms:ifme,kfms:kfme,jfms:jfme,1:msize))
allocate(   u(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(   y(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(   r(ifms:ifme,kfms:kfme,jfms:jfme))
kmat = 0.
u = 0.
y = 0.

print *,'copying the input data'
Kmat(1:s(1),1:s(2),1:s(3),1:s(4))=Kmat_m(1:s(1),1:s(2),1:s(3),1:s(4))
u(1:s(1),1:s(2),1:s(3))=u_m(1:s(1),1:s(2),1:s(3))
           
write(*,'(a)')'ntd_mult now computing r = y - Kmat*u, calling with y=0'
call ndt_mult(  &
  ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
  ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
  ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire patch bounds
  ifts, ifte, kfts, kfte, jfts,jfte,                        & ! fire tile bounds
  kmat, u, y, r, siz, relres)

write(*,'(a,3i8)')'copying the resulting -y=K*u to array size ',n
allocate(y_m(1:n(1),1:n(2),1:n(3)))
y_m(1:s(1),1:s(2),1:s(3))=-r(1:s(1),1:s(2),1:s(3))

call write_array_nd(reshape(y_m,(/product(n)/)),n,'y')

end program ndt_mult_test
