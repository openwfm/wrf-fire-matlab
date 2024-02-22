program lin_alg_test
use module_lin_alg
implicit none
integer, parameter :: n=3
real, dimension(n,n) :: A, Ainv, Id, res
integer i,j

A  = reshape([1., 2., 3., 4., 5., 6., 7., 8., 10.],[3, 3])
Id = reshape([1., 0., 0., 0., 1., 0., 0., 0., 1.],[3, 3])

!do i=1,n
!    do j=1,n
!        call random_number(A(i,j))
!        if(i .eq. j)then
!            id(i,j)=1.
!        else
!            id(i,j)=0.
!        endif
!    enddo
!enddo

call print_matrix('A',A)
call inv3( A, Ainv)
call print_matrix('Ainv',Ainv)
res = matmul(A, Ainv) 
call print_matrix('A*Ainv',res)
res = res - Id
call print_matrix('A*Ainv - I',res)

print *,'lin_alg_test: error=',maxval(abs(res))
end program lin_alg_test