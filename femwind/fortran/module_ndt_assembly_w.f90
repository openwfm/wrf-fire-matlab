module module_ndt_assembly
       use module_hexa

contains

subroutine ndt_assembly(                              &
    ifds, ifde, kfds,kfde, jfds, jfde,            & ! fire grid dimensions
    ifms, ifme, kfms,kfme, jfms, jfme,            &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte,            &
    A, u0,v0,w0, lambda, iflags1,iflags2,          & !Input from femwind, u0, v0, w0 
    U,V,W)                                          !U,V,W  
!Purpose: Create Arrays of Wind Vector Component Values at Center Points of Spatial Grid
!In:
!A Coefficient Matrix size 3X3, symmetric positive definite
!u0, v0, w0  Initial wind speed values in x,y,z direction at grid cell centers
!iflags1 iflags =1  indicates add initial wind to calculated wind
!iflags2 iflags2 = 1 returns Kloc and Jg from hexa, iflags2 = 2 returns Floc and Jg from hexa
!out:
! U,V,W Computed wind values in x,y,z direction 

implicit none


!*** arguments

integer, intent(in)::                     &
    ifds, ifde, kfds,kfde, jfds, jfde,            & ! fire grid dimensions
    ifms, ifme, kfms,kfme, jfms, jfme,            &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte            




real, intent(in), dimension(3,3):: A
real, intent(in), dimension(ifms:ifme, kfms:kfme, jfms: jfme):: u0,v0,w0!spatial grid
integer, intent(in) :: iflags1
!Input for hexa
integer, intent(in)::iflags2 = 1.

real, intent(out), dimension(ifms:ifme, kfms:kfme, jfms:jfme,1:msize)::U,V,W


!*** local

integer:: ie1, ie2, ie3, ic1, ic2, ic3, iloc, i &
          kloc, k1, k2, k3
integer, dimension(8):: kglo                     !global index
real ::  Kloc(8,8), Floc(8), Jg(8,3)
real ::  Xloc(3,8), u0loc(3)
real :: grad(3)
         
!*** u0loc is an input for module_hexa, but is not used to construct K. Do I need to define this?
!*** integer, dimension(3,1,1), save ::iflags = reshape((/1,0,1/),(/3,1,1/)) !define iflags to construct JG and Kloc in hexa



!** executable

do ie2=jfts,jfte
    do ie3=kfts, kfte
        do ie1=ifts, ifte
            do ic2=0,1
                do ic3=0,1
                    do ic1=0,1
                          kloc = 1 + ic1 + 2*(ic2 + 2*ic3)
                          k1 = ie1 + ic1
                          k2 = ie2 + ic2 
                          k3 = ie3 + ic3
                          kglo(kloc) = k1 + n(1)*((k2-1) + n(2)*(k3-1))  
                    enddo
                enddo
            enddo
            call hexa(A,Xloc,u0loc,Kloc,Floc,Jg,iflags2)
            grad =transpose(lambda(kglo))*Jg
            call inv(grad, grad_inv)
            grad = matmult(grad_inv,A)
            U(ie1, ie2, ie3)=grad(1)
            V(ie1, ie2, ie3)=grad(2)
            W(ie1, ie2, ie3)=grad(3)
            
        enddo
    enddo
enddo

if.iflags1.eq.1
             U = U + u0
             V = V + v0
             W = W + w0
end if
end subroutine ndt_assembly
end module  module_ndt_assembly
        

                            
                        
                        
            














