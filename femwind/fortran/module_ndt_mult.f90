module module_ndt_mult

use module_utils

contains

subroutine ndt_mult(                              &
    ifds, ifde, kfds,kfde, jfds, jfde,                       & ! fire grid dimensions
    ifms, ifme, kfms,kfme, jfms, jfme,            &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte,             &
    kmat, lambda, y, r, siz, relres)

implicit none

!*** purpose: compute r = y - Kglo*lambda and 2-norm of r

!*** arguments

integer, intent(in)::                             &
    ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
    ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
    ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte                            ! fire tile bounds

integer, parameter:: msize = 14
real, intent(in), dimension(ifms:ifme,kfms:kfme,jfms:jfme,msize):: kmat  ! global stiffness matrix
real,intent(in),  dimension(ifms:ifme,kfms:kfme,jfms:jfme):: lambda, y          ! input vectors 
real,intent(out), dimension(ifms:ifme,kfms:kfme,jfms:jfme):: r          ! output vector 
real, intent(out):: siz, relres
!*** local

integer:: i,j,k,ie,je,ke  
real:: t, res

!** executable
ie = snode(ifte,ifde,+1)
je = snode(jfte,jfde,+1)
ke = snode(kfte,kfde,+1)

siz = 0.
res = 0.

do j=jfts,je
  do k=kfts,ke
    do i=ifts,ie
      t = y(i,k,j) - ( &
        kmat(i-1,k-1,j-1,14)*lambda(i-1,k-1,j-1) +  &
        kmat(i  ,k-1,j-1,13)*lambda(i  ,k-1,j-1) +  &
        kmat(i+1,k-1,j-1,12)*lambda(i+1,k-1,j-1) +  &
        kmat(i-1,k-1,j  ,11)*lambda(i-1,k-1,j  ) +  &
        kmat(i  ,k-1,j  ,10)*lambda(i  ,k-1,j  ) +  &
        kmat(i+1,k-1,j  , 9)*lambda(i+1,k-1,j  ) +  &
        kmat(i-1,k-1,j+1, 8)*lambda(i-1,k-1,j+1) +  &
        kmat(i  ,k-1,j+1, 7)*lambda(i  ,k-1,j+1) +  &
        kmat(i+1,k-1,j+1, 6)*lambda(i+1,k-1,j+1) +  &
        kmat(i-1,k  ,j-1, 5)*lambda(i-1,k  ,j-1) +  &
        kmat(i  ,k  ,j-1, 4)*lambda(i  ,k  ,j-1) +  &
        kmat(i+1,k  ,j-1, 3)*lambda(i+1,k  ,j-1) +  &
        kmat(i-1,k  ,j  , 2)*lambda(i-1,k  ,j  ) +  &  ! K(I,J)*x(J)   I=(i-1,j,k)  storing upper triangle of K only, K(I,J) stored in another row
        kmat(i  ,k  ,j  , 1)*lambda(i  ,k  ,j  ) +  &  ! K(I,I)*x(I)   I=(i,j,k)
        kmat(i  ,k  ,j  , 2)*lambda(i+1,k  ,j  ) +  &  ! K(I,J)*x(J)   J=(i+1,j,k)
        kmat(i  ,k  ,j  , 3)*lambda(i-1,k  ,j+1) +  &  ! etc rest of the row I in upper triangle
        kmat(i  ,k  ,j  , 4)*lambda(i  ,k  ,j+1) +  &
        kmat(i  ,k  ,j  , 5)*lambda(i+1,k  ,j+1) +  &
        kmat(i  ,k  ,j  , 6)*lambda(i-1,k+1,j-1) +  &
        kmat(i  ,k  ,j  , 7)*lambda(i  ,k+1,j-1) +  &
        kmat(i  ,k  ,j  , 8)*lambda(i+1,k+1,j-1) +  &
        kmat(i  ,k  ,j  , 9)*lambda(i-1,k+1,j  ) +  &
        kmat(i  ,k  ,j  ,10)*lambda(i  ,k+1,j  ) +  &
        kmat(i  ,k  ,j  ,11)*lambda(i+1,k+1,j  ) +  &
        kmat(i  ,k  ,j  ,12)*lambda(i-1,k+1,j+1) +  &
        kmat(i  ,k  ,j  ,13)*lambda(i  ,k+1,j+1) +  &
        kmat(i  ,k  ,j  ,14)*lambda(i+1,k+1,j+1) )
      r(i,k,j) = t
      siz = max(siz,abs(y(i,k,j)))
      res = max(res,abs(t))
    enddo
  enddo
enddo

relres = res/max(tiny(siz),siz)

         
end subroutine ndt_mult

end module module_ndt_mult


