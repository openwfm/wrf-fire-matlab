module module_boundary_conditions

use module_utils

contains

subroutine ndt_boundary_conditions(                              &
    ifds, ifde, kfds, kfde, jfds, jfde,         & ! fire grid dimensions
    ifms, ifme, kfms, kfme, jfms, jfme,         &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,         & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts, jfte,         &    
    kmat)

implicit none

!*** arguments

integer, intent(in)::                             &
    ifds, ifde, kfds, kfde, jfds, jfde,         & ! fire grid dimensions
    ifms, ifme, kfms, kfme, jfms, jfme,         &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,         & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts, jfte             
integer, parameter:: msize = 14
real, intent(inout), dimension(ifms:ifme,kfms:kfme,jfms:jfme,msize):: kmat  ! global stiffness matrix

!*** local

real:: s
integer:: i,j,k,ie,je,ke

!** executable

! loop upper bounds
ie = snode(ifte,ifde,+1)
je = snode(jfte,jfde,+1)
ke = snode(kfte,kfde,+1)

! scale
s=0.

do j=jfts,je
  do k=kfts,ke
    do i=ifts,ie
       s=max(s,abs(kmat(i  ,k  ,j  , 1)))
    enddo
  enddo
enddo

do j=jfts,je
  do k=kfts,ke
    do i=ifts,ie
      ! not efficient but will be executed only once
      if(i.eq.ifds.or.i.eq.ifde+1.or.j.eq.jfds.or.j.eq.jfde+1.or.k.eq.kfde+1)then
        ! replace the row/col (i,k,j) by scaled identity
        kmat(i-1,k-1,j-1,14)=0.
        kmat(i  ,k-1,j-1,13)=0.
        kmat(i+1,k-1,j-1,12)=0.
        kmat(i-1,k-1,j  ,11)=0.
        kmat(i  ,k-1,j  ,10)=0.
        kmat(i+1,k-1,j  , 9)=0.
        kmat(i-1,k-1,j+1, 8)=0.
        kmat(i  ,k-1,j+1, 7)=0.
        kmat(i+1,k-1,j+1, 6)=0.
        kmat(i-1,k  ,j-1, 5)=0.
        kmat(i  ,k  ,j-1, 4)=0.
        kmat(i+1,k  ,j-1, 3)=0.
        kmat(i-1,k  ,j  , 2)=0.
        kmat(i  ,k  ,j  , 1)=s
        kmat(i  ,k  ,j  , 2)=0.
        kmat(i  ,k  ,j  , 3)=0.
        kmat(i  ,k  ,j  , 4)=0.
        kmat(i  ,k  ,j  , 5)=0.
        kmat(i  ,k  ,j  , 6)=0.
        kmat(i  ,k  ,j  , 7)=0.
        kmat(i  ,k  ,j  , 8)=0.
        kmat(i  ,k  ,j  , 9)=0.
        kmat(i  ,k  ,j  ,10)=0.
        kmat(i  ,k  ,j  ,11)=0.
        kmat(i  ,k  ,j  ,12)=0.
        kmat(i  ,k  ,j  ,13)=0.
        kmat(i  ,k  ,j  ,14)=0.
      endif
    enddo
  enddo
enddo
          
end subroutine ndt_boundary_conditions

subroutine vec_boundary_conditions(               &
    ifds, ifde, kfds, kfde, jfds, jfde,         & ! fire grid dimensions
    ifms, ifme, kfms, kfme, jfms, jfme,         &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,         & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts, jfte,         &    
    F)

implicit none

!*** arguments

integer, intent(in)::                             &
    ifds, ifde, kfds, kfde, jfds, jfde,         & ! fire grid dimensions
    ifms, ifme, kfms, kfme, jfms, jfme,         &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,         & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts, jfte    
integer, parameter:: msize = 14
real, intent(inout), dimension(ifms:ifme,kfms:kfme,jfms:jfme):: F  ! corner-based scalar field

!*** local

integer:: i,j,k,ie,je,ke

!** executable

! loop upper bounds
ie = snode(ifte,ifde,+1)
je = snode(jfte,jfde,+1)
ke = snode(kfte,kfde,+1)

do j=jfts,je
  do k=kfts,ke
    do i=ifts,ie
      ! not efficient, change later 
      if(i.eq.ifds.or.i.eq.ifde+1.or.j.eq.jfds.or.j.eq.jfde+1.or.k.eq.kfde+1)then
        F(i,k,j)=0.
      endif
    enddo
  enddo
enddo
          
end subroutine vec_boundary_conditions
end module module_boundary_conditions


