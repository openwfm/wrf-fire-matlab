program w_assembly_test

use module_w_assembly
use module_hexa
use module_utils

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
                Xmat(:,:,:),Ymat(:,:,:),Zmat(:,:,:), X(:,:,:),Y(:,:,:),Z(:,:,:),A_m(:,:,:)
               

real, pointer :: a1(:), a2(:)
integer :: n1(2),lambda_dim(3),u_dim(3), x_dim(3)
real :: Amat(3,3)
integer :: &
    ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
    ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
    ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte                           ! fire tile bounds
                              
  

integer :: i,j,k,jx
integer :: usize(3)

call read_array(A_m,'A_test')
Amat = reshape(A_m,(/3,3/)) 

call read_array(lambda, 'lambda_test')
lambda_dim = shape(lambda)

print *, 'lambda array has shape', shape(lambda)

! read input arrays in ikj index ordering and tight bounds
call read_array(X,'X_test')
call read_array(Y,'Y_test')
call read_array(Z,'Z_test')

call read_array(u0,'u0_test')        !Recovering Inital Wind Arrays
call read_array(v0,'v0_test')        
call read_array(w0,'w0_test')        

x_dim = shape(X)
u_dim = shape(u0)

!Checking that dimensions of lambda and are X arrays are consistent
if (x_dim(1).ne.lambda_dim(1).or.x_dim(2).ne.lambda_dim(2).or.x_dim(3).ne.lambda_dim(3)) then
    call crash('Lambda dimensions must equal the dimensions of X')
    stop
endif

ifts = 1
ifte = x_dim(1)-1
jfts = 1
jfte = x_dim(3)-1
kfts = 1
kfte = x_dim(2)-1
ifms = ifts-1
ifme = ifte+2
jfms = jfts-1
jfme = jfte+2
kfms = kfts-1
kfme = kfte+2
ifds = ifts
ifde = ifte
jfds = jfts
jfde = jfte
kfds = kfts
kfde = kfte

allocate(lambdamat(ifts:ifte+1, kfts:kfte+1, jfts:jfte+1))



allocate(Xmat(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(Ymat(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(Zmat(ifms:ifme,kfms:kfme,jfms:jfme))

allocate(u0mat(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(v0mat(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(w0mat(ifms:ifme,kfms:kfme,jfms:jfme))

! copy the input data to tile sized bounds
! X Y Z are corner based, upper bound larger by one
do j=jfts,jfte+1
  do k=kfts,kfte+1
    do i=ifts,ifte+1
        Xmat(i,k,j) = X(i,k,j)
        Ymat(i,k,j) = Y(i,k,j)
        Zmat(i,k,j) = Z(i,k,j)
        lambdamat(i,k,j) = lambda(i,k,j)
    enddo
  enddo
enddo

do j=jfts,jfte
  do k=kfts,kfte
    do i=ifts,ifte
        u0mat(i,k,j) = u0(i,k,j)
        v0mat(i,k,j) = v0(i,k,j)
        w0mat(i,k,j) = w0(i,k,j)
    enddo
  enddo
enddo

print *, 'u0 is', u0


allocate(U(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(V(ifms:ifme,kfms:kfme,jfms:jfme))
allocate(W(ifms:ifme,kfms:kfme,jfms:jfme))

U = 0.
V = 0.
W = 0.

!write(*,'(a)')'calling w_assembly'
call w_assembly(  &
  ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
  ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
  ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire patch bounds
  ifts, ifte, kfts, kfte, jfts,jfte,                        & ! fire tile bounds
  lambdamat,u0mat, v0mat, w0mat, Amat,                      & !Input from femwind, u0, v0, w0,Spatial Grid Data
  Xmat, Ymat, Zmat,                                         &
  U,V,W)                                    

!write(*,'(a,3i8)')'copying the output data to array size ',n2,msize

allocate(Umat(ifts:ifte,kfts:kfte,jfts:jfte))
allocate(Vmat(ifts:ifte,kfts:kfte,jfts:jfte))
allocate(Wmat(ifts:ifte,kfts:kfte,jfts:jfte))
!keep track of indexing
do j=jfts,jfte
  do k=kfts,kfte
    do i=ifts,ifte
      	    Umat(i,k,j)=U(i,k,j)
            Vmat(i,k,j)=V(i,k,j)
            Wmat(i,k,j)=W(i,k,j)
    enddo
  enddo
enddo
!print *, 'Shape of Umat is', shape(Umat)

call write_array(Umat,'U_test')
call write_array(Vmat,'V_test')
call write_array(Wmat,'W_test')


end program w_assembly_test
