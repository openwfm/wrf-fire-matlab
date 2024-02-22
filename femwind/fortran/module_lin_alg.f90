module module_lin_alg

use module_utils

contains

subroutine inv3(M, M_inv)
!Purpose: Calculate the inverse of a 3X3 matrix
!In:
!M  3X3 matrix
!Out:
!M_inv Inverse of M 
implicit none
real,intent(in), dimension(3,3):: M
real,intent(out), dimension(3,3):: M_inv

!Local Variables

real :: detM


    ! Compute M_inv
    detM  = (M(1,1)*M(2,2)*M(3,3) - M(1,1)*M(2,3)*M(3,2)-&
             M(1,2)*M(2,1)*M(3,3) + M(1,2)*M(2,3)*M(3,1)+&
             M(1,3)*M(2,1)*M(3,2) - M(1,3)*M(2,2)*M(3,1))

    if(abs(detM).lt.10.0*tiny(detM))then
        print *,'detM=',detM
        call crash('The matrix is numerically singular')
    endif

    detM = 1.0/detM

    M_inv(1,1) = +detM * (M(2,2)*M(3,3) - M(2,3)*M(3,2))
    M_inv(2,1) = -detM * (M(2,1)*M(3,3) - M(2,3)*M(3,1))
    M_inv(3,1) = +detM * (M(2,1)*M(3,2) - M(2,2)*M(3,1))
    M_inv(1,2) = -detM * (M(1,2)*M(3,3) - M(1,3)*M(3,2))
    M_inv(2,2) = +detM * (M(1,1)*M(3,3) - M(1,3)*M(3,1))
    M_inv(3,2) = -detM * (M(1,1)*M(3,2) - M(1,2)*M(3,1))
    M_inv(1,3) = +detM * (M(1,2)*M(2,3) - M(1,3)*M(2,2))
    M_inv(2,3) = -detM * (M(1,1)*M(2,3) - M(1,3)*M(2,1))
    M_inv(3,3) = +detM * (M(1,1)*M(2,2) - M(1,2)*M(2,1))

end subroutine inv3
end module module_lin_alg 
