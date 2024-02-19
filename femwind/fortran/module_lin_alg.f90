module module_lin_alg

use module_utils

contains

subroutine inv3(M, Inv_M)
!Purpose: Calculate the inverse of a 3X3 matrix
!In:
!M  3X3 matrix
!Out:
!M_inv Inverse of M 
implicit none
real,intent(in), dimension(3,3):: M
real,intent(out), dimension(3,3):: Inv_M

!Local Variables
real, dimension(3,3)::M_T
real :: det_M, det_M_inv

!!Compute Inverse of M
det_M =   M(1,1)*(M(2,2)*M(3,3)-M(3,2)*M(2,3)) -  &
          M(1,2)*(M(1,2)*M(3,3) - M(3,1)*M(2,3))+ &
          M(1,3)*(M(2,1)*M(2,2) - M(3,1)*M(3,2))

if(abs(det_M).lt.10*tiny(det_M))then
        print *,'det_M=',det_M
        call crash('Inv3: The matrix is numerically singular')
    endif

det_M_inv = 1./det_M

M_T   =  transpose(M)

Inv_M(1,1) = det_M_inv*(M_T(2,2)*M_T(3,3) - M_T(2,3)*M_T(3,2))
Inv_M(1,2) = det_M_inv*(M_T(2,1)*M_T(3,3) - M_T(3,1)*M_T(2,3))
Inv_M(1,3) = det_M_inv*(M_T(2,1)*M_T(3,2) - M_T(2,2)*M_T(3,1))
Inv_M(2,1) = det_M_inv*(M_T(1,2)*M_T(3,3) - M_T(3,2)*M_T(1,3))
Inv_M(2,2) = det_M_inv*(M_T(1,1)*M_T(3,3) - M_T(3,1)*M_T(1,3))
Inv_M(2,3) = det_M_inv*(M_T(1,1)*M_T(3,2) - M_T(3,1)*M_T(1,2))
Inv_M(3,1) = det_M_inv*(M_T(1,2)*M_T(2,3) - M_T(1,3)*M_T(2,2))
Inv_M(3,2) = det_M_inv*(M_T(1,1)*M_T(2,3) - M_T(2,1)*M_T(1,3))
Inv_M(3,3) = det_M_inv*(M_T(1,1)*M_T(2,2) - M_T(1,2)*M_T(2,1))

end subroutine inv3
end module module_lin_alg 
