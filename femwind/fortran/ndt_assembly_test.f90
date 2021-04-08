program ndt_assembly_test

use module_ndt_assembly
use module_hexa
use module_io_matlab

implicit none

real, pointer:: Amat(:,:), Xmat(:,:,:),Ymat(:,:,:), Zmat(:,:,:), Kmat(:,:,:,:),&
		Kmat_m(:,:,:,:), Amat_m(:,:),Xmat_m(:,:,:),Ymat_m(:,:,:), Zmat_m(:,:,:)
real, pointer::a1(:), a2(:), a3(:), a4(:) 
integer ::n1(2),n2(3),n3(3), n4(3), m(3)




integer :: msize, &
    ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
    ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
    ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte,                        & ! fire tile bounds
    iats, iate, jats, jate, iams,iame, 			      &	        ! Amat bounds
    jams, jame

integer :: i,j,k,jx
integer :: iflags(3) = reshape((/1,0,1/),(/3/))                 !Flags to construct K in hexa module
integer:: ksize(4)                                              ! Global Stifness Matrix Dimensions in Matlab

call read_array_nd(a1,n1,'A') !Recovering X-Matrix and dimension of X matrix
allocate(Amat_m(n1(1),n1(2))) !Do this again for Y, Z
Amat_m = reshape(a1,n1)
iats = 1
iate = n1(1)
jats = 1
jate = n1(2)

iams = iats
iame = iate
jams = jats 
jame = jate 

allocate(Amat(iams:iame,jams:jame))

do j=jats,jate
    do i=iats,iate
        Amat(i,j) = Amat_m(i,j)	
  enddo
enddo
	

! read input arrays in ijk index ordering and tight bounds
call read_array_nd(a2,n2,'X')        !Recovering X-Matrix and dimension of X matrix
allocate(Xmat_m(n2(1),n2(2),n2(3)))  !Do this again for Y, Z
Xmat_m = reshape(a2,n2)

call read_array_nd(a3,n3,'Y') !Recovering X-Matrix and dimension of X matrix
allocate(Ymat_m(n3(1),n3(2),n3(3))) 
Ymat_m = reshape(a3,n3)

call read_array_nd(a4,n4,'Z') !Recovering X-Matrix and dimension of X matrix
allocate(Zmat_m(n4(1),n4(2),n4(3))) 
Zmat_m = reshape(a4,n4)


ifts = 1
ifte = n2(1)
jfts = 1 
jfte = n2(2)
kfts = 1
kfte = n2(3)
msize = 14

ifms = ifts
ifme = ifte
jfms = jfts
jfme = jfte
kfms = kfts
kfme = kfte

ksize = (/ifte,jfte,kfte,msize/)


! allocate a little bigger with zeros in extra areas
allocate(Kmat(ifms:ifme,kfms:kfme,jfms:jfme, 1:msize))

allocate(Xmat(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(Ymat(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(Zmat(ifms:ifme,kfms:kfme,jfms:jfme))


Kmat = 0.
Xmat = 9999.
Ymat = 9999.
Zmat = 9999.


! copy the input data 
do j=jfts,jfte
  do k=kfts,kfte
    do i=ifts,ifte
        Xmat(i,k,j) = Xmat_m(i,j,k)
	Ymat(i,k,j) = Ymat_m(i,j,k)
	Zmat(i,k,j) = Zmat_m(i,j,k)	
    enddo
  enddo
enddo

write(*,'(a)')'calling ndt_assembly'
call ndt_assembly(  &
   ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
  ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
  ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire patch bounds
  ifts, ifte, kfts, kfte, jfts,jfte,                        & ! fire tile bounds
  Amat,Xmat,Ymat,Zmat, iflags, Kmat)

write(*,'(a,3i8)')'copying the output data to array size ',n2,msize
allocate(Kmat_m(ifts:ifte,jfts:jfte,kfts:kfte, 1:msize))
do j=jfts,jfte
  do k=kfts,kfte
    do i=ifts,ifte
        do jx = 1,msize
      	    Kmat_m(i,j,k,jx)=Kmat(i,k,j,jx)
        enddo
    enddo
  enddo
enddo


call write_array_nd(reshape(Kmat_m,(/product(ksize)/)),ksize,'K')
!***print *,(/product(ksize)/)


end program ndt_assembly_test
