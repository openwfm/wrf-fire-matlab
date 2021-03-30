module module_ndt_assembly
contains
subroutine ndt_assembly(                              &
    ifds, ifde, kfds,kfde, jfds, jfde,            & ! fire grid dimensions
    ifms, ifme, kfms,kfme, jfms, jfme,            &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte,            &
    A, X, msize, u0, iflags,                                     & !Input from femwind
    Kmat)                                            & !Global stiffness matrix output  

implicit none


!*** arguments

integer, intent(in)::
    ifds, ifde, kfds,kfde, jfds, jfde,                       & ! fire grid dimensions
    ifms, ifme, kfms,kfme, jfms, jfme,            &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte,            &




integer, parameter:: msize = 14
integer, intent(in)::nn                             !product of fire domain bounds

real, intent(in), dimension(3,ifms:ifme, kfms:kfme, jfms: jfme)::X !spatial grid
!Input for hexa
real, intent(in):: A(3,3), Xloc(3,8), u0(3)    
integer, intent(in)::iflags(3)
real, intent(out):: Kloc(8,8), Floc(8), Jg(8,3)

real, intent(out), dimension(ifds:ifde, kfds:kfde, jfds:jfde,msize)::Kmat


!*** local

integer:: ie1, ie2, ie3, ic1, ic2, ic3, iloc, i

real, intent(in), dimension(3,8)::Xloc

!** executable
z=0.0
do ie3=jfts,jfter -1
    do ie2=kfts, kfte -1
        do ie1=ifts, ifte -1
            do ic3=0,1
                do ic2=0,1
                    do ic1=0,1
                        iloc=1+ic1+2*(ic2+2*ic3);  !local index of the node in the element
                        do i=1,3
                            Xloc(i,iloc)=X(i,ie1 + ic1, ie2 + ic2, ie3 + ic3)
                        enddo
                    enddo
                enddo
            enddo
            call hexa(A,Xloc,u0,iflags,Kloc,Floc,Jg)
            Kmat(i,k,j, 1) = Kmat(i,k,j, 1) + Kloc( 1, 1)+ Kloc( 2, 2)+ Kloc( 3, 3) + Kloc( 4, 4)&
                             + Kloc( 5, 5)+ Kloc( 6, 6) + Kloc( 7, 7) + Kloc( 8, 8)             
            Kmat(i,k,j, 2) = Kmat(i,k,j, 2)+Kloc( 8, 7) + Kloc( 7, 8) + Kloc( 6, 5) + Kloc( 5, 6)&
                             + Kloc( 4, 3) + Kloc( 3, 4) + Kloc( 1, 2) + Kloc( 2, 1) 
            Kmat(i,k,j, 3) = Kmat(i,k,j, 3)+Kloc( 7, 6) + Kloc( 6, 7) + Kloc( 3, 2) + Kloc( 2, 3)
            Kmat(i,k,j, 4) = Kmat(i,k,j, 4)+Kloc( 8, 6) + Kloc( 7, 5) + Kloc( 6, 8) + Kloc( 5, 7)&
                             + Kloc( 4, 2) + Kloc( 3, 1) + Kloc( 1, 3) + Kloc( 2, 4)
            Kmat(i,k,j, 5) = Kmat(i,k,j, 5) + Kloc( 8, 5) + Kloc( 5, 8) + Kloc( 4, 1) + Kloc( 1, 4) 
            Kmat(i,k,j, 6) = Kmat(i,k,j, 6 + Kloc( 5, 4) + Kloc( 4, 5) 
            Kmat(i,k,j, 7) = Kmat(i,k,j, 7) + Kloc( 6, 4) + Kloc( 5, 3) + Kloc( 4, 6) + Kloc( 3, 5)
            Kmat(i,k,j, 8) = Kmat(i,k,j, 8) + Kloc( 6, 3) + Kloc( 3, 6)
            Kmat(i,k,j, 9) = Kmat(i,k,j, 9) + Kloc( 7, 4) + Kloc( 5, 2) + Kloc( 4, 7) + Kloc( 2, 5)
            Kmat(i,k,j,10) = Kmat(i,k,j,10) + Kloc( 8, 4) + Kloc( 7, 3) + Kloc( 6, 2) + Kloc( 5, 1)&
                             Kloc( 4, 8) + Kloc( 3, 7) + Kloc( 2, 6) + Kloc( 1, 5)
            Kmat(i,k,j,11) = Kmat(i,k,j,11) + Kloc( 8, 3) + Kloc( 6, 1) + Kloc( 3, 8) + Kloc( 1, 6) 
            Kmat(i,k,j,12) = Kmat(i,k,j,12) + Kloc( 7, 2) + Kloc( 2, 7) 
            Kmat(i,k,j,13) = Kmat(i,k,j,13)+ Kloc( 8, 2) +Kloc( 2, 8)+ Kloc( 7, 1) +Kloc( 1, 7)
            Kmat(i,k,j,14) = Kmat(i,k,j,14)+ Kloc( 8, 1) + Kloc( 1, 8)
        enddo
    enddo
enddo

end subroutine ndt_assembly
end module ndt_assembly
        

                            
                        
                        
            














