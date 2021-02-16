program hexa_test

use module_hexa   ! testing only
use module_io_matlab ! to read and write matrices as text files from matlab

implicit none

real:: A(3,3), X(3,8), u0(3)    ! fortran is not case sensitive
integer::iflags(3)
real:: Kloc(8,8), Floc(8), Jg(8,3)

call hexa(A,X,u0,Kloc,Floc,Jg,iflags)
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

end program hexa_test
