module module_w_assembly
       use module_hexa
       use module_lin_alg 

contains

subroutine w_assembly(                              &
    ifds, ifde, kfds,kfde, jfds, jfde,            & ! fire grid dimensions
    ifms, ifme, kfms,kfme, jfms, jfme,            &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte,            &
    A, u0,v0,w0, lambda, iflags1, iflags2,n2,     & !Input from femwind, u0, v0, w0
    X, Y, Z,                                      & !Spatial Grid Data     
    U,V,W)                                          !U,V,W  
!Purpose: Create Arrays of Wind Vector Component Values at Center Points of Spatial Grid
!In:
!A Coefficient Matrix size 3X3, symmetric positive definite
!u0, v0, w0  Initial wind speed values in x,y,z direction at grid cell centers
!iflags1 iflags1 = 1 returns Kloc and Jg from hexa, iflags2 = 2 returns Floc and Jg from hexa
!iflags2 iflags2 =1  indicates add initial wind to calculated wind
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
real, intent(in), dimension(ifms:ifme, kfms:kfme, jfms: jfme):: u0,v0,w0
integer, intent(in) :: iflags1, iflags2
integer, intent(in) :: n2(3)
integer, intent(in) :: lambda(product(n2))
!Input for hexa


real, intent(out), dimension(ifms:ifme, kfms:kfme, jfms:jfme)::U,V,W


!*** local

integer:: ie1, ie2, ie3, ic1, ic2, ic3, iloc, i, &
          kloci, k1, k2, k3
real ::  Kloc(8,8), Floc(8), Jg(8,3)
real ::  Xloc(3,8), u0loc(3)
real :: grad(3)
real :: kgloi(8)
real :: A_inv(3,3)
Xloc = 99999.
Jg = 0.
Kloc = 0.
Floc = 0.
grad = 0.

        
!*** u0loc is an input for module_hexa, but is not used to construct K. Do I need to define this?
!** executable

do ie2=jfts,jfte-1
    do ie3=kfts, kfte-1
        do ie1=ifts, ifte-1
            do ic2=0,1
                do ic3=0,1
                    do ic1=0,1
                        iloc=1+ic1+2*(ic2+2*ic3);  !local index of the node in the element
                        Xloc(1,iloc)=X(ie1 + ic1, ie3 + ic3, ie2 + ic2)
                        Xloc(2,iloc)=Y(ie1 + ic1, ie3 + ic3, ie2 + ic2)
                        Xloc(3,iloc)=Z(ie1 + ic1, ie3 + ic3, ie2 + ic2)
                        kloci = 1 + ic1 + 2*(ic2 + 2*ic3)
                        k1 = ie1 + ic1
                        k2 = ie2 + ic2 
                        k3 = ie3 + ic3
                        kgloi(kloci) = k1 + n2(1)*((k2-1) + n2(3)*(k3-1))  
                    enddo
                enddo
            enddo
            call hexa(A,Xloc,u0loc,Kloc,Floc,Jg,iflags1)
            grad = matmul(transpose(Jg),lambda(kgloi))
            call Inv3(A, A_inv)
            grad = matmul(transpose(A_inv),grad)
            U(ie1, ie2, ie3)=grad(1)
            V(ie1, ie2, ie3)=grad(2)
            W(ie1, ie2, ie3)=grad(3) 
        enddo
    enddo
enddo

if (iflags2.eq.1) then
             U = U + u0
             V = V + v0
             W = W + w0
end if
end subroutine w_assembly
end module  module_w_assembly
        

                            
                        
                        
            














