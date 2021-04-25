program ndt_assembly_test

use module_ndt_assembly
use module_hexa
use module_utils

implicit none

real:: Amat(3,3)
real, pointer:: Xmat(:,:,:),Ymat(:,:,:), Zmat(:,:,:), Kmat(:,:,:,:)
real, pointer::a1(:), a2(:), a3(:), a4(:), X(:,:,:), Y(:,:,:), Z(:,:,:), Km(:,:,:,:)
integer ::n1(2),n2(3),n3(3), n4(3), m(3)

integer :: msize, &
    ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
    ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
    ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte,                        & ! fire tile bounds
    iats, iate, jats, jate, iams,iame, 			      &	        ! Amat bounds
    jams, jame

integer :: i,j,k,jx
integer :: iflags = 1                !Flags to construct K in hexa module
integer:: ksize(4)                               ! Global Stifness Matrix Dimensions in Matlab

call read_array_nd(a1,n1,'A') !Recovering X-Matrix and dimension of X matrix
if (n1(1).ne.3.or.n1(2).ne.3)then
    call crash('A must be 3 by 3')
    stop
endif

Amat = reshape(a1,n1)
    
! read input arrays in ikj index ordering and tight bounds
call read_array(X,'X')        !Recovering X-Matrix and dimension of X matrix
call read_array(Y,'Y')        !Recovering X-Matrix and dimension of X matrix
call read_array(Z,'Z')        !Recovering X-Matrix and dimension of X matrix

n2 = shape(X)

ifts = 1
ifte = n2(1)
jfts = 1 
jfte = n2(3)
kfts = 1
kfte = n2(2)
msize = 14

ifms = ifts
ifme = ifte
jfms = jfts
jfme = jfte
kfms = kfts
kfme = kfte


allocate(Kmat(ifms:ifme,kfms:kfme,jfms:jfme, 1:msize))

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

!write(*,'(a)')'calling ndt_assembly'
call ndt_assembly(  &
  ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
  ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
  ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire patch bounds
  ifts, ifte, kfts, kfte, jfts,jfte,                        & ! fire tile bounds
  Amat,Xmat,Ymat,Zmat, iflags, Kmat)

!write(*,'(a,3i8)')'copying the output data to array size ',n2,msize
allocate(Km(ifts:ifte,kfts:kfte,jfts:jfte, 1:msize))
do j=jfts,jfte
  do k=kfts,kfte
    do i=ifts,ifte
        do jx = 1,msize
      	    Km(i,k,j,jx)=Kmat(i,k,j,jx)
        enddo
    enddo
  enddo
enddo


ksize = (/ifte-ifts+1,kfte-kfts+1,jfte-jfts+1,msize/)
call write_array_nd(reshape(Km,(/product(ksize)/)),ksize,'K')
!***print *,(/product(ksize)/)


end program ndt_assembly_test
