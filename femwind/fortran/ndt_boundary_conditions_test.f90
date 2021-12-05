program ndt_boundary_conditions_test

use module_boundary_conditions   ! testing only
use module_utils ! to read and write matrices as text files from matlab

implicit none

real, pointer:: Kmat(:,:,:,:), &  ! fortran is not case sensitive
                Kmat_m(:,:,:,:),a(:)
integer :: s(4)

integer :: msize, &
    ifds, ifde, kfds, kfde, jfds, jfde,         & ! fire grid dimensions
    ifms, ifme, kfms, kfme, jfms, jfme,         &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,         & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts, jfte
integer :: i,j,k,jx

! read input arrays in ikj index ordering and tight bounds
call read_array_nd(a,s,'K')
allocate(Kmat_m(s(1),s(2),s(3),s(4)))
Kmat_m = reshape(a,s)

ifts = 1
ifte = s(1) - 1
kfts = 1
kfte = s(2) - 1
jfts = 1
jfte = s(3) - 1
ifds = 1
ifde = s(1) - 1
kfds = 1
kfde = s(2) - 1
jfds = 1
jfde = s(3) - 1
msize = s(4)
if(msize.ne.14)call crash('msize must be 14')
ifms = ifts-1
ifme = ifte+2
jfms = jfts-1
jfme = jfte+2
kfms = kfts-1
kfme = kfte+2

! allocate a little bigger with zeros in extra areas
print *,'allocating Kmat size',ifms,ifme,kfms,kfme,jfms,jfme,1,msize
allocate(Kmat(ifms:ifme,kfms:kfme,jfms:jfme,1:msize))
Kmat = 0.

print *,'copying the input matrix'
Kmat(1:s(1),1:s(2),1:s(3),1:s(4))=Kmat_m(1:s(1),1:s(2),1:s(3),1:s(4))
!do j=1,s(1)
!  do k=1,s(2)
!    do i=1,s(3)
!      do jx = 1,s(4)
!        Kmat(i,k,j,jx) = Kmat_m(i,k,j,jx)
!      enddo
!    enddo
!  enddo
!enddo
           
write(*,'(a)')'calling ntd_boundary_conditions'
call ndt_boundary_conditions(  &
    ifds, ifde, kfds, kfde, jfds, jfde,         & ! fire grid dimensions
    ifms, ifme, kfms, kfme, jfms, jfme,         &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,         & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts, jfte,         &    
    Kmat)

print *,'copying the output matrix'
Kmat_m(1:s(1),1:s(2),1:s(3),1:s(4))=Kmat(1:s(1),1:s(2),1:s(3),1:s(4))
!do j=1,s(1)
!  do k=1,s(2)
!    do i=1,s(3)
!      do jx = 1,s(4)
!        Kmat_m(i,k,j,jx) = Kmat(i,k,j,jx)
!      enddo
!    enddo
!  enddo
!enddo

call write_array_nd(reshape(Kmat_m,(/product(s)/)),s,'Kb')

end program ndt_boundary_conditions_test
