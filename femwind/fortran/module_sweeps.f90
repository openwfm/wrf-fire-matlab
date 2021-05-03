module module_sweeps
use module_utils

contains

subroutine sweeps(ifds, ifde, kfds,kfde, jfds, jfde,               &     ! fire grid dimensions
                  ifms, ifme, kfms,kfme, jfms, jfme,               &
                  ifps, ifpe, kfps, kfpe, jfps, jfpe,              &     ! fire patch bounds
                  ifts, ifte, kfts, kfte, jfts,jfte,               &
                  Kmat, F, x)                                            !input and output matrices/vectors

implicit none

!*** arguments

integer, intent(in)::                                              &
    ifds, ifde, kfds, kfde, jfds, jfde,                            &     ! fire domain bounds
    ifms, ifme, kfms, kfme, jfms, jfme,                            &     ! fire memory bounds
    ifps, ifpe, kfps, kfpe, jfps, jfpe,                            &     ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte                                    ! fire tile bounds


integer, parameter:: msize = 14
real, intent(in), dimension(ifms:ifme,kfms:kfme,jfms:jfme,msize):: Kmat  ! global stiffness matrix
real, intent(in), dimension(ifms:ifme,kfms:kfme,jfms:jfme):: F           ! global load vector
real,intent(out), dimension(ifms:ifme,kfms:kfme,jfms:jfme):: x           ! output vector

!*** local

integer:: i, j, k, r1, r2 

!*** executable
 
call write_array(F(ifts: ifte, kfts: kfte, jfts:jfte),'F_sweeps_in')
call write_array(x(ifts: ifte, kfts: kfte, jfts:jfte),'x_sweeps_in')

do r1 = 1,2
    do r2 = 1,2
        do i = r1,ifte,2
           do j =r2,jfte,2 
              do k = 1,kfte
                  x(i,k,j)=                                         &
                            (1/Kmat(i  ,k  ,j  , 1))*               &
                            ( F(i,k,j) -                            &
                            (                                       &
                             Kmat(i-1,k-1,j-1,14)*x(i-1,k-1,j-1) +  &
                             Kmat(i  ,k-1,j-1,13)*x(i  ,k-1,j-1) +  &
                             Kmat(i+1,k-1,j-1,12)*x(i+1,k-1,j-1) +  &
                             Kmat(i-1,k-1,j  ,11)*x(i-1,k-1,j  ) +  &
                             Kmat(i  ,k-1,j  ,10)*x(i  ,k-1,j  ) +  &
                             Kmat(i+1,k-1,j  , 9)*x(i+1,k-1,j  ) +  &
                             Kmat(i-1,k-1,j+1, 8)*x(i-1,k-1,j+1) +  &
                             Kmat(i  ,k-1,j+1, 7)*x(i  ,k-1,j+1) +  &
                             Kmat(i+1,k-1,j+1, 6)*x(i+1,k-1,j+1) +  &
                             Kmat(i-1,k  ,j-1, 5)*x(i-1,k  ,j-1) +  &
                             Kmat(i  ,k  ,j-1, 4)*x(i  ,k  ,j-1) +  &
                             Kmat(i+1,k  ,j-1, 3)*x(i+1,k  ,j-1) +  &
                             Kmat(i-1,k  ,j  , 2)*x(i-1,k  ,j  ) +  &
                             Kmat(i  ,k  ,j  , 2)*x(i+1,k  ,j  ) +  &
                             Kmat(i  ,k  ,j  , 3)*x(i-1,k  ,j+1) +  &
                             Kmat(i  ,k  ,j  , 4)*x(i  ,k  ,j+1) +  &
                             Kmat(i  ,k  ,j  , 5)*x(i+1,k  ,j+1) +  &
                             Kmat(i  ,k  ,j  , 6)*x(i-1,k+1,j-1) +  &
                             Kmat(i  ,k  ,j  , 7)*x(i  ,k+1,j-1) +  &
                             Kmat(i  ,k  ,j  , 8)*x(i+1,k+1,j-1) +  &
                             Kmat(i  ,k  ,j  , 9)*x(i-1,k+1,j  ) +  &
                             Kmat(i  ,k  ,j  ,10)*x(i  ,k+1,j  ) +  &
                             Kmat(i  ,k  ,j  ,11)*x(i+1,k+1,j  ) +  &
                             Kmat(i  ,k  ,j  ,12)*x(i-1,k+1,j+1) +  &
                             Kmat(i  ,k  ,j  ,13)*x(i  ,k+1,j+1) +  &
                             Kmat(i  ,k  ,j  ,14)*x(i+1,k+1,j+1)    &
                            )                                       &
                            )
              end do
           end do  
        end do   
    end do
end do

call write_array(x(ifts: ifte, kfts: kfte, jfts:jfte),'x_sweeps_out')
stop
end subroutine sweeps

end module module_sweeps
