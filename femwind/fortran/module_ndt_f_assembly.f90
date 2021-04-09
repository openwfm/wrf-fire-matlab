module module_ndt_f_assembly
       use module_hexa
contains
subroutine ndt_f_assembly(                        &
    ifds, ifde, kfds,kfde, jfds, jfde,            & ! fire grid dimensions
    ifms, ifme, kfms,kfme, jfms, jfme,            &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte,            &
    A, X,Y,Z, iflags,                         & !Input from femwind, U0, V0, W0 not used in hexa to construct Kloc or JG
    F, F_dim)                                             !Global load vector output  

implicit none


!*** arguments

integer, intent(in)::                             &
    ifds, ifde, kfds,kfde, jfds, jfde,            & ! fire grid dimensions
    ifms, ifme, kfms,kfme, jfms, jfme,            &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte            




real, intent(in), dimension(3,3):: A
real, intent(in), dimension(ifms:ifme, kfms:kfme, jfms: jfme):: X,Y,Z!spatial grid
!Input for hexa
integer, intent(in)::iflags

integer,intent(in)::F_dim

real,intent(out), dimension(1:F_dim)::F

!*** local

integer:: ie1, ie2, ie3, ic1, ic2, ic3, iloc, k1, k2, k3, i
real :: Kloc(8,8), Floc(8), Jg(8,3)
real ::  Xloc(3,8), u0(3) 
real :: kglo(8)
!*** u0loc is an input for module_hexa, but is not used to construct K. Do I need to define this?
!*** integer, dimension(3,1,1), save ::iflags = reshape((/1,0,1/),(/3,1,1/)) !define iflags to construct JG and Kloc in hexa


Xloc = 9999.
F = 0.
kglo = 0.

!** executable
do ie3=jfts,jfte -1
    do ie2=kfts, kfte -1
        do ie1=ifts, ifte -1
            do ic3=0,1
                do ic2=0,1
                    do ic1=0,1
                        iloc=1+ic1+2*(ic2+2*ic3);  !local index of the node in the element
                        k1 = ie1 + ic1; 
                        k2 = ie2+ic2;
                        k3 = ie3+ic3;
                        kglo(iloc) = k1+ifte*((k2-1)+jfte*(k3-1));
                            Xloc(1,iloc)=X(ie1 + ic1, ie2 + ic2, ie3 + ic3)
                            Xloc(2,iloc)=Y(ie1 + ic1, ie2 + ic2, ie3 + ic3)
                            Xloc(3,iloc)=Z(ie1 + ic1, ie2 + ic2, ie3 + ic3)

                    enddo
                enddo
            enddo
            call hexa(A,Xloc,u0,Kloc,Floc,Jg,iflags)
            do i = 1,8
                F(kglo(i)) = F(kglo(i)) + Floc(i)
            enddo
        enddo
    enddo
enddo

end subroutine ndt_f_assembly
end module  module_ndt_f_assembly
