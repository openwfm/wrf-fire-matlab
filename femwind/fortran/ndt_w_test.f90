program w_assembly_test

use module_w_assembly
use module_hexa
use module_io_matlab

!Purpose: Create Arrays of Wind Vector Component Values at Center Points of Spatial Grid
!In:
!A Coefficient Matrix size 3X3, symmetric positive definite
!u0, v0, w0  Initial wind speed values in x,y,z direction at grid cell centers
!X,Y,Z       3-D Physical Location in the Mesh Grid
!iflags1 iflags1 = 1 returns Kloc and Jg from hexa, iflags2 = 2 returns Floc and Jg from hexa
!iflags2 iflags2 =1  indicates add initial wind to calculated wind
!out:
!U,V,W Computed wind values in x,y,z direction 

implicit none

real, pointer:: u0mat(:,:,:),v0mat(:,:,:), w0mat(:,:,:), Umat(:,:,:),              &
                Vmat(:,:,:), Wmat(:,:,:), u0(:,:,:), v0(:,:,:), w0(:,:,:),         &
                U(:,:,:), V(:,:,:),W(:,:,:),lambda(:,:,:), lambdamat(:,:,:),       &  ! Calculated final windFinal 
                Xmat(:,:,:),Ymat(:,:,:),Zmat(:,:,:), X(:,:,:),Y(:,:,:),Z(:,:,:)
               

real, pointer :: a1(:), a2(:)
integer :: n1(2),lambda_dim(3),u_dim(3), x_dim(3)
real :: Amat(3,3)
integer :: &
    ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
    ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
    ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte,                        & ! fire tile bounds
    iuts, iute, kuts, kute, juts, jute,                       &
    iums, iume, kums, kume, jums, jume
                              
  

integer :: i,j,k,jx
integer :: aflags(2) = (/3,1/)                 !Set iflags=1 to construct K in hexa module, iflags = 3 to construct Jg
!integer :: iflags2 = 1
integer :: usize(3)

call read_array_nd(a1,n1,'A') !Recovering X-Matrix and dimension of X matrix
if (n1(1).ne.3.or.n1(2).ne.3)then
    call crash('A must be 3 by 3')
    stop
endif

Amat = reshape(a1,n1)

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
if (x_dim(1).ne.lambda_dim(1).or.x_dim(2).ne.lambda_dim(2).or.x_dim(3).ne.lambda_dim(3)) then
    call crash('Lambda dimensions must equal the dimensions of X')
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

iuts = 1
iute = u_dim(1)
juts = 1
jute = u_dim(3)
kuts = 1
kute = u_dim(2)
iums = iuts - 1
iume = iute + 1
jums = juts - 1
jume = jute + 1
kums = kuts - 1
kume = kute + 1



allocate(lambdamat(ifts:ifte, kfts:kfte, jfts:jfte))


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
do j=juts,jute
  do k=kuts,kute
    do i=iuts,iute
    u0mat(i,k,j) = u0(i,k,j)
	v0mat(i,k,j) = v0(i,k,j)
	w0mat(i,k,j) = w0(i,k,j)	
    enddo
  enddo
enddo

allocate(U(iums:iume,kums:kume,jums:jume))
allocate(V(iums:iume,kums:kume,jums:jume))
allocate(W(iums:iume,kums:kume,jums:jume))

U = 0.
V = 0.
W = 0.

!write(*,'(a)')'calling w_assembly'
call w_assembly(  &
  ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
  ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
  ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire patch bounds
  ifts, ifte, kfts, kfte, jfts,jfte,                        & ! fire tile bounds
  iuts, iute, kuts, kute, juts, jute,                       &
  iums, iume, kums, kume, jums, jume,                       &
  lambdamat,u0, v0, w0,                                        & !Input from femwind, u0, v0, w0
  Amat, X, Y, Z,                                            & !Spatial Grid Data     
  U,V,W)                                    

!write(*,'(a,3i8)')'copying the output data to array size ',n2,msize

allocate(Umat(iuts:iute,kuts:kute,juts:jute))
allocate(Vmat(iuts:iute,kuts:kute,juts:jute))
allocate(Wmat(iuts:iute,kuts:kute,juts:jute))
!keep track of indexing
do j=juts,jute
  do k=kuts,kute
    do i=iuts,iute
      	    Umat(i,k,j)=U(i,k,j)
            Vmat(i,k,j)=V(i,k,j)
            Wmat(i,k,j)=W(i,k,j)
    enddo
  enddo
enddo


usize = (/iute-iuts+1,kute-kuts+1,jute-juts+1/)
call write_array_nd(reshape(Umat,(/product(usize)/)),usize,'U')
call write_array_nd(reshape(Vmat,(/product(usize)/)),usize,'V')
call write_array_nd(reshape(Wmat,(/product(usize)/)),usize,'W')


end program w_assembly_test
