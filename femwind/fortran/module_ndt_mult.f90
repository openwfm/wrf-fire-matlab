module module_ndt_mult

contains

subroutine ndt_mult(                              &
    ifds, ifde, jfds, jfde,                       & ! fire grid dimensions
    ifms, ifme, kfms,kfme, jfms, jfme,            &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts,ifte, kfts, kfte, jfts,jfte,             &
    kmat, u, y)

implicit none

!*** arguments

integer, intent(in)::                             &
    ifds, ifde, jfds, jfde,                       & ! fire domain bounds
    ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
    ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire patch bounds
    ifts,ifte, kfts, kfte, jfts,jfte                            ! fire tile bounds

integer, parameter:: msize = 14
real, intent(in), dimension(ifms:ifme,kfms:kfme,jfms:jfme,msize):: kmat  ! global stiffness matrix
real,intent(in),  dimension(ifms:ifme,kfms:kfme,jfms:jfme):: u          ! input vector 
real,intent(out), dimension(ifms:ifme,kfms:kfme,jfms:jfme):: y          ! output vector 

!*** local

integer:: i,j,k  

!** executable

do j=jfts,jfte
  do k=kfts,kfte
    do i=ifts,ifte
      y(i,k,j)= &
        kmat(i-1,k-1,j-1,14)*u(i-1,k-1,j-1) +  &
        kmat(i  ,k-1,j-1,13)*u(i  ,k-1,j-1) +  &
        kmat(i+1,k-1,j-1,12)*u(i+1,k-1,j-1) +  &
        kmat(i-1,k-1,j  ,11)*u(i-1,k-1,j  ) +  &
        kmat(i  ,k-1,j  ,10)*u(i  ,k-1,j  ) +  &
        kmat(i+1,k-1,j  , 9)*u(i+1,k-1,j  ) +  &
        kmat(i-1,k-1,j+1, 8)*u(i-1,k-1,j+1) +  &
        kmat(i  ,k-1,j+1, 7)*u(i  ,k-1,j+1) +  &
        kmat(i+1,k-1,j+1, 6)*u(i+1,k-1,j+1) +  &
        kmat(i-1,k  ,j-1, 5)*u(i-1,k  ,j-1) +  &
        kmat(i  ,k  ,j-1, 4)*u(i  ,k  ,j-1) +  &
        kmat(i+1,k  ,j-1, 3)*u(i+1,k  ,j-1) +  &
        kmat(i-1,k  ,j  , 2)*u(i-1,k  ,j  ) +  &
        kmat(i  ,k  ,j  , 1)*u(i  ,k  ,j  ) +  &
        kmat(i  ,k  ,j  , 2)*u(i+1,k  ,j  ) +  &
        kmat(i  ,k  ,j  , 3)*u(i-1,k  ,j+1) +  &
        kmat(i  ,k  ,j  , 4)*u(i  ,k  ,j+1) +  &
        kmat(i  ,k  ,j  , 5)*u(i+1,k  ,j+1) +  &
        kmat(i  ,k  ,j  , 6)*u(i-1,k+1,j-1) +  &
        kmat(i  ,k  ,j  , 7)*u(i  ,k+1,j-1) +  &
        kmat(i  ,k  ,j  , 8)*u(i+1,k+1,j-1) +  &
        kmat(i  ,k  ,j  , 9)*u(i-1,k+1,j  ) +  &
        kmat(i  ,k  ,j  ,10)*u(i  ,k+1,j  ) +  &
        kmat(i  ,k  ,j  ,11)*u(i+1,k+1,j  ) +  &
        kmat(i  ,k  ,j  ,12)*u(i-1,k+1,j+1) +  &
        kmat(i  ,k  ,j  ,13)*u(i  ,k+1,j+1) +  &
        kmat(i  ,k  ,j  ,14)*u(i+1,k+1,j+1) 
    enddo
  enddo
enddo
          
end subroutine ndt_mult

end module module_ndt_mult


