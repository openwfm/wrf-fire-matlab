module module_w_assembly
       use module_hexa
       use module_lin_alg 

contains

subroutine w_assembly(                            &
    ifds, ifde, kfds,kfde, jfds, jfde,            & ! fire grid dimensions
    ifms, ifme, kfms,kfme, jfms, jfme,            &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte,            &
    lambda,u0, v0, w0, A,  X, Y, Z,               & !Input from femwind, u0, v0, w0, Spatial Grid Data
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

integer, intent(in)::                             &
    ifds, ifde, kfds,kfde, jfds, jfde,            & ! fire grid dimensions
    ifms, ifme, kfms,kfme, jfms, jfme,            &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte            
     
real, intent(in), dimension(ifms:ifme,kfms:kfme,jfms:jfme):: lambda

!Input for hexa
real, intent(in) :: A(3,3)
real, intent(in), dimension(ifms:ifme, kfms:kfme, jfms:jfme)::X,Y,Z,u0, v0, w0 

!Output
real, intent(out), dimension(ifms:ifme, kfms:kfme, jfms:jfme)::U,V,W

!*** local

integer:: ie1, ie2, ie3, ic1, ic2, ic3, iloc, i, k1, k2, k3
real ::  Kloc(8,8), Floc(8), Jg(8,3)
real ::  Xloc(3,8), u0loc(3)
real :: grad(3)
real :: lambda_loc(8)
real :: A_inv(3,3)
real :: vol

!*** executable

lambda_loc = 0.
Xloc = 99999.
Jg = 0.
Kloc = 0.
Floc = 0.
grad = 0.
u0loc =0.
       
!*** u0loc is an input for module_hexa, but is not used to construct K. Do I need to define this?
!** executable

!print *, 'u0 vector is', u0
do ie2=jfts,jfte
    do ie3=kfts, kfte
        do ie1=ifts, ifte
            ! constant part
            do ic2=0,1
                do ic3=0,1
                    do ic1=0,1
                        Xloc(1,iloc)=X(ie1 + ic1, ie3 + ic3, ie2 + ic2)
                        Xloc(2,iloc)=Y(ie1 + ic1, ie3 + ic3, ie2 + ic2)
                        Xloc(3,iloc)=Z(ie1 + ic1, ie3 + ic3, ie2 + ic2)
                    enddo
                enddo
            enddo
            u0loc(1) = u0(ie1,ie3,ie2)
            u0loc(2) = v0(ie1, ie3,ie2)
            u0loc(3) = w0(ie1, ie3,ie2)
            !fine print *, 'local lambda is', lambda_loc 
            !print* , 'Xloc is', Xloc
            !print* , 'u0loc is', u0loc(1)  
            call hexa(A,Xloc,u0loc,Kloc,Floc,Jg,vol,3)
            !print*, 'Jg is', Jg
            !print*, shape(jg)
            !*** end of constant part
            do ic2=0,1
                do ic3=0,1
                    do ic1=0,1
                        iloc=1+ic1+2*(ic2+2*ic3);  !local index of the node in the element
                        k1 = ie1 + ic1
                        k2 = ie2 + ic2 
                        k3 = ie3 + ic3
                        lambda_loc(iloc) = lambda(k1,k3,k2)  
                    enddo
                enddo
            enddo
            grad = matmul(transpose(Jg),lambda_loc)
            !not ok print *,'Grad before multiplication by A_inv is', grad
  
            call Inv3(A, A_inv)
            grad = matmul(transpose(A_inv),grad)
            ! Not ok print *,'Grad after multiplication by A_inv is', grad
            
            U(ie1, ie3, ie2)=grad(1)+ u0(ie1, ie3, ie2)
            V(ie1, ie3, ie2)=grad(2)+ v0(ie1, ie3, ie2)
            W(ie1, ie3, ie2)=grad(3)+ w0(ie1, ie3, ie2)
            
        enddo
       !print *, 'lambda array', lambda
    enddo
enddo
!print *,'Shape of U', shape(U)

end subroutine w_assembly
end module  module_w_assembly
        

                            
                        
                        
            














