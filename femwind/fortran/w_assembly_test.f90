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
real, pointer:: u0mat(:,:,:),v0mat(:,:,:), w0mat(:,:,:), Umat(:,:,:),              &
                Vmat(:,:,:), Wmat(:,:,:), u0(:,:,:), v0(:,:,:), w0(:,:,:),         &
                U(:,:,:), V(:,:,:),W(:,:,:)                                         ! Calculated final windFinal 
                Xmat(:,:,:),Ymat(:,:,:),Zmat(:,:,:), X(:,:,:),Y(:,:,:),Z(:,:,:)
real:: lambda(120) 
real, pointer::a1(:), a2(:)
integer ::dim_A(2),dim_lam(1),u_dim(3), x_dim(3)

integer :: msize, &
    ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
    ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
    ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte,                        & ! fire tile bounds
    iats, iate, jats, jate, iams,iame, 			              & ! Amat bounds
    jams, jame,                                               &
    iuds, iude, kuds, kude, juds, jude,                       & ! Wind tile and and memory bounds
    iums, iume, kums, kume, jums, jume

integer :: i,j,k,jx
integer :: aflags(2) = (/3,1/)                 !Set iflags=1 to construct K in hexa module, iflags = 3 to construct Jg
!integer :: iflags2 = 1
integer :: usize(3)
call read_array_nd(a1,dim_A,'A') !Recovering X-Matrix and dimension of X matrix
if (n1(1).ne.3.or.n1(2).ne.3)then
    call crash('A must be 3 by 3')
    stop
endif

Amat = reshape(a1,dim_A)

call read_array_nd(a2,dim_lam, 'lambda')


lambda = reshape(a2,dim_lam)

! read input arrays in ikj index ordering and tight bounds
call read_array(X,'X')
call read_array(Y,'Y')
call read_array(Z,'Z')

call read_array(u0,'u0')        !Recovering Inital Wind Arrays
call read_array(v0,'v0')        
call read_array(w0,'w0')        

x_dim = shape(X)
u_dim = shape(u0)


ifts = 1
ifte = x_dim(1)
jfts = 1
jfte = x_dim(3)
kfts = 1
kfte = x_dim(2)
ifms = ifts - 1
ifme = ifte + 1
jfms = jfts - 1 
jfme = jfte + 1
kfms = kfts - 1
kfme = kfte + 1 

iuds = 1
iude = u_dim(1)
juds = 1
jude = u_dim(3)
kuds = 1
kude = u_dim(2)
iums = iuds-1
iume = iude+1
jums = juds-1
jume = jude+1
kums = kuds-1
kume = kude+1


allocate(Xmat(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(Ymat(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(Zmat(ifms:ifme,kfms:kfme,jfms:jfme))

allocate(u0mat(iums:iume,kums:kume,jums:jume))
allocate(v0mat(iums:iume,kums:kume,jums:jume))
allocate(w0mat(iums:iume,kums:kume,jums:jume))




do j=jfts,jfte
  do k=kfts,kfte
    do i=ifts,ifte
        Xmat(i,k,j) = X(i,k,j)
        Ymat(i,k,j) = Y(i,k,j)
        Zmat(i,k,j) = Z(i,k,j)
    enddo
  enddo
enddo

! copy the input data to tile sized bounds
do j=jums,jume
  do k=kums,kume
    do i=iums,iume
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
  iuds, iude, kuds, kude, juds, jude,                       & ! Wind tile and and memory bounds
  iums, iume, kums, kume, jums, jume 
  A, u0,v0,w0, lambda, aflags, n,     & !Input from femwind, u0, v0, w0
  X, Y, Z,                                      & !Spatial Grid Data     
  U,V,W)                                    

!write(*,'(a,3i8)')'copying the output data to array size ',n2,msize

!allocate(U(ifms:ifme,kfms:kfme,jfms:jfme))
!allocate(V(ifms:ifme,kfms:kfme,jfms:jfme))
!allocate(W(ifms:ifme,kfms:kfme,jfms:jfme))
!keep track of indexing
!do j=jfts,jfte
  !do k=kfts,kfte
    !do i=ifts,ifte
      	    !Umat(i,k,j)=U(i,k,j)
            !Vmat(i,k,j)=V(i,k,j)
            !Wmat(i,k,j)=W(i,k,j)
    !enddo
  !enddo
1enddo


usize = (/ifte-ifts+1,kfte-kfts+1,jfte-jfts+1/)
call write_array_nd(reshape(U,(/product(usize)/)),usize,'U')
call write_array_nd(reshape(V,(/product(usize)/)),usize,'V')
call write_array_nd(reshape(W,(/product(usize)/)),usize,'W')
!***print *,product(ksize)


end program w_assembly_test
