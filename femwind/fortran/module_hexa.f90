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
integer, intent(in)::iflags
real, intent(out):: Kloc(8,8), Floc(8), Jg(8,3)
!*** local variables
!real, parameter :: g = 0.5773502691896257
real,dimension(9,3),save :: ib = reshape((/-1,-1,-1,-1,1,1,1,1,0,-1,-1,1,1,-1,-1,1,1,0,-1,1,-1,1,-1,1,-1,1,0/),(/9,3/))
real,dimension(9,3),save :: s = reshape((/-0.5773502691896257,-0.5773502691896257,-0.5773502691896257,-0.5773502691896257&
,0.5773502691896257,0.5773502691896257,0.5773502691896257,0.5773502691896257,0.0,-0.5773502691896257,-0.5773502691896257&
,0.5773502691896257,0.5773502691896257,-0.5773502691896257,-0.5773502691896257,0.5773502691896257,0.5773502691896257&
,0.0,-0.5773502691896257,0.5773502691896257,-0.5773502691896257,0.5773502691896257,-0.5773502691896257,0.5773502691896257&
,-0.5773502691896257,0.5773502691896257,0.0/),(/9,3/))
real :: gradf(8,3), Jg_tmp(8,3)
real :: Jx(3,3), Jx_inv(3,3)
real :: Jg_tran(3,8), A_tmp(3,8)
real :: K_at_s(8,8)
real :: tmp_mat(8)
real :: u0_tmp(3)
integer :: i,j,k,m
real :: detJx = 0
real :: tmp = 0
real :: vol = 0
Kloc = 0.
Floc = 0.
Jg = 0.
!*** executable

!Calculate Jg loop
if (iflags .eq.  3) then

do i=1,9

!Calculate gradf
!first column of gradf
do j = 1,8
gradf(j,1) = (ib(j,1)*(1+ib(j,2)*s(i,2))*(1+ib(j,3)*s(i,3)))/8
end do
!second coumn of gradf
do j = 1,8
gradf(j,2) = ((1+ib(j,1)*s(i,1))*ib(j,2)*(1+ib(j,3)*s(i,3)))/8
end do
!third column of gradf
do j = 1,8
gradf(j,3) = ((1+ib(j,1)*s(i,1))*(1+ib(j,2)*s(i,2))*ib(j,3))/8
end do

Jx = matmul(X,gradf)

!!!Compute Jx_inv!!!
detJx = (Jx(1,1)*Jx(2,2)*Jx(3,3) - Jx(1,1)*Jx(2,3)*Jx(3,2)-&
         Jx(1,2)*Jx(2,1)*Jx(3,3) + Jx(1,2)*Jx(2,3)*Jx(3,1)+&
         Jx(1,3)*Jx(2,1)*Jx(3,2) - Jx(1,3)*Jx(2,2)*Jx(3,1))

Jx_inv(1,1) = +(1/detJx) * (Jx(2,2)*Jx(3,3) - Jx(2,3)*Jx(3,2))
Jx_inv(2,1) = -(1/detJx) * (Jx(2,1)*Jx(3,3) - Jx(2,3)*Jx(3,1))
Jx_inv(3,1) = +(1/detJx) * (Jx(2,1)*Jx(3,2) - Jx(2,2)*Jx(3,1))
Jx_inv(1,2) = -(1/detJx) * (Jx(1,2)*Jx(3,3) - Jx(1,3)*Jx(3,2))
Jx_inv(2,2) = +(1/detJx) * (Jx(1,1)*Jx(3,3) - Jx(1,3)*Jx(3,1))
Jx_inv(3,2) = -(1/detJx) * (Jx(1,1)*Jx(3,2) - Jx(1,2)*Jx(3,1))
Jx_inv(1,3) = +(1/detJx) * (Jx(1,2)*Jx(2,3) - Jx(1,3)*Jx(2,2))
Jx_inv(2,3) = -(1/detJx) * (Jx(1,1)*Jx(2,3) - Jx(1,3)*Jx(2,1))
Jx_inv(3,3) = +(1/detJx) * (Jx(1,1)*Jx(2,2) - Jx(1,2)*Jx(2,1))

Jg = matmul(gradf,Jx_inv)

end do
endif !end for computing Jg

!Calculate Kloc loop
if (iflags .eq.  1) then

do i=1,9

!Calculate gradf
!first column of gradf
do j = 1,8
gradf(j,1) = (ib(j,1)*(1+ib(j,2)*s(i,2))*(1+ib(j,3)*s(i,3)))/8
end do
!second coumn of gradf
do j = 1,8
gradf(j,2) = ((1+ib(j,1)*s(i,1))*ib(j,2)*(1+ib(j,3)*s(i,3)))/8
end do
!third column of gradf
do j = 1,8
gradf(j,3) = ((1+ib(j,1)*s(i,1))*(1+ib(j,2)*s(i,2))*ib(j,3))/8
end do

Jx = matmul(X,gradf)

!!!Compute Jx_inv!!!
detJx = (Jx(1,1)*Jx(2,2)*Jx(3,3) - Jx(1,1)*Jx(2,3)*Jx(3,2)-&
         Jx(1,2)*Jx(2,1)*Jx(3,3) + Jx(1,2)*Jx(2,3)*Jx(3,1)+&
         Jx(1,3)*Jx(2,1)*Jx(3,2) - Jx(1,3)*Jx(2,2)*Jx(3,1))

Jx_inv(1,1) = +(1/detJx) * (Jx(2,2)*Jx(3,3) - Jx(2,3)*Jx(3,2))
Jx_inv(2,1) = -(1/detJx) * (Jx(2,1)*Jx(3,3) - Jx(2,3)*Jx(3,1))
Jx_inv(3,1) = +(1/detJx) * (Jx(2,1)*Jx(3,2) - Jx(2,2)*Jx(3,1))
Jx_inv(1,2) = -(1/detJx) * (Jx(1,2)*Jx(3,3) - Jx(1,3)*Jx(3,2))
Jx_inv(2,2) = +(1/detJx) * (Jx(1,1)*Jx(3,3) - Jx(1,3)*Jx(3,1))
Jx_inv(3,2) = -(1/detJx) * (Jx(1,1)*Jx(3,2) - Jx(1,2)*Jx(3,1))
Jx_inv(1,3) = +(1/detJx) * (Jx(1,2)*Jx(2,3) - Jx(1,3)*Jx(2,2))
Jx_inv(2,3) = -(1/detJx) * (Jx(1,1)*Jx(2,3) - Jx(1,3)*Jx(2,1))
Jx_inv(3,3) = +(1/detJx) * (Jx(1,1)*Jx(2,2) - Jx(1,2)*Jx(2,1))
Jg = matmul(gradf,Jx_inv)

!check to calc Kloc
if (i < 9) then
do j = 1,8
do k = 1,3
Jg_tran(k,j) = Jg(j,k)
end do
enddo

A_tmp = matmul(A, Jg_tran)
K_at_s = matmul(Jg,A_tmp)
Kloc = Kloc+(K_at_s*abs(detJx))
end if !end for computing Kloc
end do !end of outter most do loop
end if !end for computing Kloc


!Calculate Floc loop
if (iflags .eq. 2) then

do i=1,9

!Calculate gradf
!first column of gradf
do j = 1,8
gradf(j,1) = (ib(j,1)*(1+ib(j,2)*s(i,2))*(1+ib(j,3)*s(i,3)))/8
end do
!second coumn of gradf
do j = 1,8
gradf(j,2) = ((1+ib(j,1)*s(i,1))*ib(j,2)*(1+ib(j,3)*s(i,3)))/8
end do
!third column of gradf
do j = 1,8
gradf(j,3) = ((1+ib(j,1)*s(i,1))*(1+ib(j,2)*s(i,2))*ib(j,3))/8
end do

Jx = matmul(X,gradf)

!!!Compute Jx_inv!!!
detJx = (Jx(1,1)*Jx(2,2)*Jx(3,3) - Jx(1,1)*Jx(2,3)*Jx(3,2)-&
         Jx(1,2)*Jx(2,1)*Jx(3,3) + Jx(1,2)*Jx(2,3)*Jx(3,1)+&
         Jx(1,3)*Jx(2,1)*Jx(3,2) - Jx(1,3)*Jx(2,2)*Jx(3,1))

Jx_inv(1,1) = +(1/detJx) * (Jx(2,2)*Jx(3,3) - Jx(2,3)*Jx(3,2))
Jx_inv(2,1) = -(1/detJx) * (Jx(2,1)*Jx(3,3) - Jx(2,3)*Jx(3,1))
Jx_inv(3,1) = +(1/detJx) * (Jx(2,1)*Jx(3,2) - Jx(2,2)*Jx(3,1))
Jx_inv(1,2) = -(1/detJx) * (Jx(1,2)*Jx(3,3) - Jx(1,3)*Jx(3,2))
Jx_inv(2,2) = +(1/detJx) * (Jx(1,1)*Jx(3,3) - Jx(1,3)*Jx(3,1))
Jx_inv(3,2) = -(1/detJx) * (Jx(1,1)*Jx(3,2) - Jx(1,2)*Jx(3,1))
Jx_inv(1,3) = +(1/detJx) * (Jx(1,2)*Jx(2,3) - Jx(1,3)*Jx(2,2))
Jx_inv(2,3) = -(1/detJx) * (Jx(1,1)*Jx(2,3) - Jx(1,3)*Jx(2,1))
Jx_inv(3,3) = +(1/detJx) * (Jx(1,1)*Jx(2,2) - Jx(1,2)*Jx(2,1))

Jg = matmul(gradf,Jx_inv)

vol = abs(detJx)*8
u0_tmp = u0*vol
do j = 1,8
tmp = 0
do k = 1,3
tmp = tmp+Jg(j,k)*u0_tmp(k)
end do
tmp_mat(j) = tmp        
end do
Floc = Floc-tmp_mat
end do
end if !end for computing Floc

end subroutine hexa

end module module_hexa
