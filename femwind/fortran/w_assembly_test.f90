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

real, pointer:: u0mat(:,:,:),v0mat(:,:,:), w0mat(:,:,:), Umat(:,:,:),              &
                Vmat(:,:,:), Wmat(:,:,:), u0(:,:,:), v0(:,:,:), w0(:,:,:),         &
                U(:,:,:), V(:,:,:),W(:,:,:),lambda(:,:,:), lambdamat(:,:,:),       &  ! Calculated final windFinal 
                Xmat(:,:,:),Ymat(:,:,:),Zmat(:,:,:), X(:,:,:),Y(:,:,:),Z(:,:,:),   &
               

real, pointer :: a1(:), a2(:)
integer :: n1(2),lambda_dim(3),u_dim(3), x_dim(3)

integer :: , &
    ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
    ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
    ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte,                        & ! fire tile bounds
    iats, iate, jats, jate, iams,iame, 			      & ! Amat bounds
    jams, jame,                                               &
    iuds, iude, kuds, kude, juds, jude,                       & ! Wind tile and and memory bounds
    iums, iume, kums, kume, jums, jume

integer :: i,j,k,jx
integer :: aflags(2) = (/3,1/)                 !Set iflags=1 to construct K in hexa module, iflags = 3 to construct Jg
!integer :: iflags2 = 1
integer :: usize(3)


!call read_array_nd(a,n,'u')
!allocate(u_m(n(1),n(2),n(3)))
!u_m = reshape(a,n)

call read_array(lambda, 'lambda')
lambda_dim = shape(lambda)


! read input arrays in ikj index ordering and tight bounds
call read_array(X,'X')
call read_array(Y,'Y')
call read_array(Z,'Z')

call read_array(u0,'u0')        !Recovering Inital Wind Arrays
call read_array(v0,'v0')        
call read_array(w0,'w0')        

x_dim = shape(X)
u_dim = shape(u0)

!Checking that dimensions of lambda and are X arrays are consistent
if (product(x_dim).ne.dim_lam) then
    call crash('Length of lambda must equal the product of dimensions of X')
    stop
endif

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

allocate(lambdamat(ifms:ifme, kfms:kfme, jfms:jfme))


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
        lambdamat(i,k,j) = lambda(i,k,j)
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
  u0mat,v0mat,w0mat, lambda, aflags,     & !Input from femwind, u0, v0, w0
  dim_lam,x_dim, Xmat, Ymat, Zmat,                                      & !Spatial Grid Data     
  U,V,W)                                    

!write(*,'(a,3i8)')'copying the output data to array size ',n2,msize

allocate(Umat(1:u_dim(1),1:u_dim(3),1:u_dim(2)))
allocate(Vmat(1:u_dim(1),1:u_dim(3),1:u_dim(2)))
allocate(Wmat(1:u_dim(1),1:u_dim(3),1:u_dim(2)))
keep track of indexing
do j=jfts,jfte
  do k=kfts,kfte
    do i=ifts,ifte
      	    Umat(i,k,j)=U(i,k,j)
            Vmat(i,k,j)=V(i,k,j)
            Wmat(i,k,j)=W(i,k,j)
    enddo
  enddo
enddo


usize = (/ifte-ifts+1,kfte-kfts+1,jfte-jfts+1/)
call write_array_nd(reshape(U,(/product(usize)/)),usize,'U')
call write_array_nd(reshape(V,(/product(usize)/)),usize,'V')
call write_array_nd(reshape(W,(/product(usize)/)),usize,'W')


end program w_assembly_test
