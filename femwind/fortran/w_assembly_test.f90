program w_assembly_test

use module_w_assembly
use module_hexa
use module_io_matlab

!Purpose: Create Arrays of Wind Vector Component Values at Center Points of Spatial Grid
!In:
!A Coefficient Matrix size 3X3, symmetric positive definite
!u0, v0, w0  Initial wind speed values in x,y,z direction at grid cell centers
!iflags1 iflags1 = 1 returns Kloc and Jg from hexa, iflags2 = 2 returns Floc and Jg from hexa
!iflags2 iflags2 =1  indicates add initial wind to calculated wind
!out:
!U,V,W Computed wind values in x,y,z direction 

implicit none

real:: Amat(3,3)
real, pointer:: u0mat(:,:,:),v0mat(:,:,:), w0mat(:,:,:), Umat(:,:,:), &
                Vmat(:,:,:), Wmat(:,:,:), u0(:,:,:), v0(:,:,:), w0(:,:,:), &
                U(:,:,:), V(:,:,:),W(:,:,:)
integer, pointer :: lambda(:) 
real, pointer::a1(:)
integer ::n1(2),n2(3),  m(3)

integer :: msize, &
    ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
    ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
    ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte,                        & ! fire tile bounds
    iats, iate, jats, jate, iams,iame, 			      &	        ! Amat bounds
    jams, jame

integer :: i,j,k,jx
integer :: iflags1 = 3                 !Set iflags=1 to construct K in hexa module, iflags = 3 to construct Jg
integer :: iflags2 = 1
integer :: usize(3)
call read_array_nd(a1,n1,'A') !Recovering X-Matrix and dimension of X matrix
if (n1(1).ne.3.or.n1(2).ne.3)then
    call crash('A must be 3 by 3')
    stop
endif

Amat = reshape(a1,n1)

! read input arrays in ikj index ordering and tight bounds
call read_array(u0,'u0')        !Recovering Inital Wind Arrays
call read_array(v0,'v0')        
call read_array(w0,'w0')        

n2 = shape(u0)

allocate(lambda(product(n2)))
lambda = 0.

ifts = 1
ifte = n2(1)
jfts = 1 
jfte = n2(3)
kfts = 1
kfte = n2(2)

ifms = ifts - 1
ifme = ifte + 1
jfms = jfts - 1
jfme = jfte + 1
kfms = kfts - 1
kfme = kfte + 1



allocate(u0mat(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(v0mat(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(w0mat(ifms:ifme,kfms:kfme,jfms:jfme))

allocate(Umat(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(Vmat(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(Wmat(ifms:ifme,kfms:kfme,jfms:jfme))

! copy the input data to tile sized bounds
do j=jfts,jfte
  do k=kfts,kfte
    do i=ifts,ifte
        u0mat(i,k,j) = u0(i,k,j)
	v0mat(i,k,j) = v0(i,k,j)
	w0mat(i,k,j) = w0(i,k,j)	
    enddo
  enddo
enddo

!write(*,'(a)')'calling w_assembly'
call w_assembly(  &
  ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
  ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
  ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire patch bounds
  ifts, ifte, kfts, kfte, jfts,jfte,                        & ! fire tile bounds
  Amat, u0,v0,w0, lambda, iflags1, iflags2,n2,          & !Input from femwind, u0, v0, w0 
    U,V,W)                                  

!write(*,'(a,3i8)')'copying the output data to array size ',n2,msize

allocate(U(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(V(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(W(ifms:ifme,kfms:kfme,jfms:jfme))
!keep track of indexing
do j=jfts,jfte
  do k=kfts,kfte
    do i=ifts,ifte
      	    U(i,k,j)=Umat(i,k,j)
            V(i,k,j)=Vmat(i,k,j)
            W(i,k,j)=Wmat(i,k,j)
    enddo
  enddo
enddo


usize = (/ifte-ifts+1,kfte-kfts+1,jfte-jfts+1/)
call write_array_nd(reshape(U,(/product(usize)/)),usize,'U')
call write_array_nd(reshape(V,(/product(usize)/)),usize,'V')
call write_array_nd(reshape(W,(/product(usize)/)),usize,'W')
!***print *,product(ksize)


end program w_assembly_test
