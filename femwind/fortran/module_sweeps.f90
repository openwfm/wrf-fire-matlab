module module_sweeps

contains

subroutine sweeps(ifds, ifde, kfds,kfde, jfds, jfde,                       & ! fire grid dimensions
    ifms, ifme, kfms,kfme, jfms, jfme,            &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte,             &
    K, x_out, x_in, F)

implicit none

!*** arguments

integer, intent(in)::                             &
    ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
    ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
    ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte                            ! fire tile bounds


implicit none



integer, parameter:: msize = 14
real, intent(in), dimension(ifms:ifme,kfms:kfme,jfms:jfme,msize):: K  ! global stiffness matrix
real, intent(in), dimension(ifms:ifme,kfms:kfme,jfms:jfms):: F
real,intent(in),  dimension(ifms:ifme,kfms:kfme,jfms:jfme):: x_in          ! input vector
real,intent(out), dimension(ifms:ifme,kfms:kfme,jfms:jfme):: x_out          ! output vector

integer:: i, j, k, r1, r2 

do r1 = 1,2
    do r2 = 1,2
        do i = rb1,ifte,2
           do j =rb2,jfte,2 
              do k = 1,kfte
                  x_out(i,k,j)= &
      ! x(i,k,j)-&
       (1/K(i  ,k  ,j  , 1)*x_in(i  ,k  ,j  ))* &

       (F(i,k,j) - &
       ( &
        K(i-1,k-1,j-1,14)*x_in(i-1,k-1,j-1) +  &
        K(i  ,k-1,j-1,13)*x_in(i  ,k-1,j-1) +  &
        K(i+1,k-1,j-1,12)*x_in(i+1,k-1,j-1) +  &
        K(i-1,k-1,j  ,11)*x_in(i-1,k-1,j  ) +  &
        K(i  ,k-1,j  ,10)*x_in(i  ,k-1,j  ) +  &
        K(i+1,k-1,j  , 9)*x_in(i+1,k-1,j  ) +  &
        K(i-1,k-1,j+1, 8)*x_in(i-1,k-1,j+1) +  &
        K(i  ,k-1,j+1, 7)*x_in(i  ,k-1,j+1) +  &
        K(i+1,k-1,j+1, 6)*x_in(i+1,k-1,j+1) +  &
        K(i-1,k  ,j-1, 5)*x_in(i-1,k  ,j-1) +  &
        K(i  ,k  ,j-1, 4)*x_in(i  ,k  ,j-1) +  &
        K(i+1,k  ,j-1, 3)*x_in(i+1,k  ,j-1) +  &
        K(i-1,k  ,j  , 2)*x_in(i-1,k  ,j  ) +  &
        K(i  ,k  ,j  , 2)*x_in(i+1,k  ,j  ) +  &
        K(i  ,k  ,j  , 3)*x_in(i-1,k  ,j+1) +  &
        K(i  ,k  ,j  , 4)*x_in(i  ,k  ,j+1) +  &
        K(i  ,k  ,j  , 5)*x_in(i+1,k  ,j+1) +  &
        K(i  ,k  ,j  , 6)*x_in(i-1,k+1,j-1) +  &
        K(i  ,k  ,j  , 7)*x_in(i  ,k+1,j-1) +  &
        K(i  ,k  ,j  , 8)*x_in(i+1,k+1,j-1) +  &
        K(i  ,k  ,j  , 9)*x_in(i-1,k+1,j  ) +  &
        K(i  ,k  ,j  ,10)*x_in(i  ,k+1,j  ) +  &
        K(i  ,k  ,j  ,11)*x_in(i+1,k+1,j  ) +  &
        K(i  ,k  ,j  ,12)*x_in(i-1,k+1,j+1) +  &
        K(i  ,k  ,j  ,13)*x_in(i  ,k+1,j+1) +  &
        K(i  ,k  ,j  ,14)*x_in(i+1,k+1,j+1)
        ) &
       ! - F(i,k,j)
        )

              end do
           end do  
        end do   
    end do
end do

end subroutine sweeps

end module module_sweeps
