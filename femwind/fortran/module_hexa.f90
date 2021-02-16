module module_hexa   ! testing only

contains

subroutine hexa(A,X,u0,Kloc,Floc,Jg,iflags)
! purpose: create local stiffness matrix etc for hexa 3d element
! in:
!   A   coefficient matrix size 3x3, symmetric positive definite
!   X   nodes coordinates size 3x8, one each column is one node
!   u0   column vector of input wind at the center of the element
!   iflags  iflags(1)>0 compute Kloc, iflags(2)>0 compute Floc, iflags(3)>0 compute Jg  
! out:
!   Kloc   local stiffness matrix
!   Floc   local divergence load vector
!   Jg     gradient at center of function with values V is V'*Jg          

implicit none

!*** arguments

real, intent(in):: A(3,3), X(3,8), u0(3)    ! fortran is not case sensitive
integer, intent(in)::iflags(3)
real, intent(out):: Kloc(8,8), Floc(8), Jg(8,3)
!, Jg(3)

!*** local variables
!real :: Nb = 8
!real :: Ng
real :: g = 0.5773502691896257
real :: ib(8,3), s(8,3), gradf(8,3), Jg_tmp(8,3)
real :: Jx(3,3), qdash(3,3), Q(3,3), R(3,3), Q_tran(3,3), R_inv(3,3)
real :: Jg_tran(3,8), A_tmp(3,8)
real :: K_at_s(8,8)
real :: tmp_mat(8)
integer :: i,j,k,m
real :: detJx = 0
real :: tmp = 0
real :: vol = 0
!*** executable

! temporary, assign someting 
! to prevent compiler warning about unassigned variables
!Kloc = 0 
!Floc = 0
!Jg = 0
!s = 0
!i=0
!j=0
!k=0
!m=0

ib = reshape((/-1,-1,-1,-1,1,1,1,1,-1,-1,1,1,-1,-1,1,1,-1,1,-1,1,-1,1,-1,1/),shape(ib))
!Still in testing process, so will most likely break

s= g*ib

if (iflags(3) > 0) then

!gradf piece
do i=1,8
do j = 1,8
do k = 1,3
gradf(j,k) = ((1+ib(j,1)*s(i,1))*(1+ib(j,2)*s(i,2))*(1+ib(j,3)*s(i,3)))/8
end do
end do
end do

Jx = matmul(X,gradf)

!!!Start QR Decomp!!!        
!Get first column of Q
tmp = 0
do j = 1,3
qdash(j,1) = Jx(j,1)
tmp = tmp + qdash(j,1)*qdash(j,1)
end do

R(1,1) = sqrt(tmp)

do j = 1,3
Q(j,1) = qdash(j,1)/R(1,1)
R(1,2) = R(1,2) + Q(j,1)*Jx(j,2)
end do

!Get second column of Q
tmp = 0
do j =1,3
qdash(j,2) = Jx(j,2)-R(1,2)*Q(j,1)
tmp = tmp + qdash(j,2)*qdash(j,2)
end do

R(2,2) = sqrt(tmp)

do j =1,3
Q(j,2) = qdash(j,2)/R(2,2)
R(1,3) = R(1,3) + Q(j,1)*Jx(j,3)
R(2,3) = R(2,3) + Q(j,2)*Jx(j,3)
end do

!Get third column of Q
tmp = 0
do j =1,3
qdash(j,3) = Jx(j,3) - R(1,3)*Q(j,1) - R(2,3)*Q(j,2)
tmp = tmp + qdash(j,3)*qdash(j,3)
end do

R(3,3) = sqrt(tmp)

do j = 1,3
Q(j,3) = qdash(j,3)/R(3,3)
end do
!!!End QR Decomp!!!

        
detJx = abs(R(1,1)*R(2,2)*R(3,3))

do j = 1,3
do k = 1,3
Q_tran(k,j) = Q(j,k)
R_inv(k,j) = R(j,k)/(detJx)
end do
end do


Jg_tmp=matmul(gradf,R_inv)
Jg = matmul(Jg_tmp,Q_tran)
end if

!check to calc Kloc
if (iflags(1) > 0) then
do j = 1,8
do k = 1,3
Jg_tran(k,j) = Jg(j,k)*detJx
end do
enddo

A_tmp = matmul(A, Jg_tran)
K_at_s = matmul(Jg,A_tmp)
Kloc = Kloc-K_at_s
end if

!check to calc Floc
if (iflags(2) > 0) then
vol = detJx*8
do j = 1,8
tmp = 0
do m = 1,3
tmp = tmp+Jg(j,k)*u0(k)
end do
tmp_mat(j) = tmp*vol        
end do
Floc = Floc-tmp_mat
end if    
end subroutine hexa

end module module_hexa
