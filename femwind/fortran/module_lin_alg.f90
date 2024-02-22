module module_lin_alg

use module_utils

contains

subroutine inv3(M, Minv)
!Purpose: Calculate the inverse of a 3X3 matrix
!In:
!M  3X3 matrix
!Out:
!M_inv Inverse of M 
implicit none
real,intent(in), dimension(:,:):: M
real,intent(out), dimension(:,:):: Minv

!Local Variables
real :: det_M, det_M_inv

!!Compute Inverse of M
det_M =   M(1,1)*(M(2,2)*M(3,3) - M(3,2)*M(2,3)) &
        - M(1,2)*(M(2,1)*M(3,3) - M(3,1)*M(2,3)) &
        + M(1,3)*(M(2,1)*M(3,2) - M(3,1)*M(2,2))

! print *,'det_M=',det_M

if(abs(det_M).lt.10*tiny(det_M))then
        call print_matrix('M',M)
        print *,'det_M=',det_M
        call crash('Inv3: The matrix is numerically singular')
endif

det_M_inv = 1./det_M

! Compute the inverse of M (Minv)

Minv(1,1) =  (M(2,2)*M(3,3) - M(2,3)*M(3,2)) * det_M_inv
Minv(1,2) = -(M(1,2)*M(3,3) - M(1,3)*M(3,2)) * det_M_inv
Minv(1,3) =  (M(1,2)*M(2,3) - M(1,3)*M(2,2)) * det_M_inv
Minv(2,1) = -(M(2,1)*M(3,3) - M(2,3)*M(3,1)) * det_M_inv
Minv(2,2) =  (M(1,1)*M(3,3) - M(1,3)*M(3,1)) * det_M_inv
Minv(2,3) = -(M(1,1)*M(2,3) - M(1,3)*M(2,1)) * det_M_inv
Minv(3,1) =  (M(2,1)*M(3,2) - M(2,2)*M(3,1)) * det_M_inv
Minv(3,2) = -(M(1,1)*M(3,2) - M(1,2)*M(3,1)) * det_M_inv
Minv(3,3) =  (M(1,1)*M(2,2) - M(1,2)*M(2,1)) * det_M_inv

end subroutine inv3

subroutine print_matrix(name,A)
character(len=1)::name
real :: A(:,:)
integer i,j
print *,'Matrix ',trim(name),lbound(A,1),':',ubound(A,1),' by ', &
                          lbound(A,2),':',ubound(A,2)
do i=lbound(A,1),ubound(A,1)
    print *,(A(i,j),j=lbound(A,2),ubound(A,2))
enddo
end subroutine print_matrix



end module module_lin_alg 
