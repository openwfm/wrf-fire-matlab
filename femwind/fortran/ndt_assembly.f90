module module_ndt_assembly
contains
subroutine ndt_assembly(                              &
    ifds, ifde, kfds,kfde, jfds, jfde,            & ! fire grid dimensions
    ifms, ifme, kfms,kfme, jfms, jfme,            &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte,            &
    A, X, u0, iflags,                                     & !Input from femwind
    K)                                             !Global stiffness matrix output  

implicit none


!*** arguments

integer, intent(in)::
    ifds, ifde, kfds,kfde, jfds, jfde,                       & ! fire grid dimensions
    ifms, ifme, kfms,kfme, jfms, jfme,            &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte,            




integer, integer(in):: msize = 14
integer, intent(in):: nn                             !product of fire domain bounds

real, intent(in), dimension(3,ifms:ifme, kfms:kfme, jfms: jfme):: X !spatial grid
!Input for hexa
real, intent(in):: A(3,3), Xloc(3,8), u0(3)    
integer, intent(in)::iflags(3)
real, intent(out):: Kloc(8,8), Floc(8), Jg(8,3)

real, intent(out), dimension(ifds:ifde, kfds:kfde, jfds:jfde,1:msize)::K


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
            K(ie1  ,ie2  ,ie3  , 1) =   K(ie1  ,ie2  ,ie3  , 1) + Kloc( 1,  1)           
            K(ie1  ,ie2  ,ie3  , 4) =   K(ie1  ,ie2  ,ie3  , 4) + Kloc( 1,  3)
            K(ie1  ,ie2  ,ie3  , 5) =   K(ie1  ,ie2  ,ie3  , 5) + Kloc( 1,  4)
            K(ie1  ,ie2  ,ie3  ,10) =   K(ie1  ,ie2  ,ie3  ,10) + Kloc( 1,  5)
            K(ie1  ,ie2  ,ie3  ,11) =   K(ie1  ,ie2  ,ie3  ,11) + Kloc( 1,  6)
            K(ie1  ,ie2  ,ie3  ,13) =   K(ie1  ,ie2  ,ie3  ,13) + Kloc( 1,  7)
            K(ie1  ,ie2  ,ie3  ,14) =   K(ie1  ,ie2  ,ie3  ,14) + Kloc( 1,  8)
            K(ie1+1,ie2  ,ie3  , 1) =   K(ie1+1,ie2  ,ie3  , 1) + Kloc( 2,  2)
            K(ie1+1,ie2  ,ie3  , 3) =   K(ie1+1,ie2  ,ie3  , 3) + Kloc( 2,  3)
            K(ie1+1,ie2  ,ie3  , 4) =   K(ie1+1,ie2  ,ie3  , 4) + Kloc( 2,  4)
            K(ie1+1,ie2  ,ie3  , 9) =   K(ie1+1,ie2  ,ie3  , 9) + Kloc( 2,  5)
            K(ie1+1,ie2  ,ie3  ,10) =   K(ie1+1,ie2  ,ie3  ,10) + Kloc( 2,  6)
            K(ie1+1,ie2  ,ie3  ,12) =   K(ie1+1,ie2  ,ie3  ,12) + Kloc( 2,  7)
            K(ie1+1,ie2  ,ie3  ,13) =   K(ie1+1,ie2  ,ie3  ,13) + Kloc( 2,  8)
            K(ie1  ,ie2+1,ie3  , 1) =   K(ie1  ,ie2+1,ie3  , 1) + Kloc( 3,  3)
            K(ie1  ,ie2+1,ie3  , 2) =   K(ie1  ,ie2+1,ie3  , 2) + Kloc( 3,  4)
            K(ie1  ,ie2+1,ie3  , 7) =   K(ie1  ,ie2+1,ie3  , 7) + Kloc( 3,  5)
            K(ie1  ,ie2+1,ie3  , 8) =   K(ie1  ,ie2+1,ie3  , 8) + Kloc( 3,  6)
            K(ie1  ,ie2+1,ie3  ,10) =   K(ie1  ,ie2+1,ie3  ,10) + Kloc( 3,  7)
            K(ie1  ,ie2+1,ie3  ,11) =   K(ie1  ,ie2+1,ie3  ,11) + Kloc( 3,  8)
            K(ie1+1,ie2+1,ie3  , 1) =   K(ie1+1,ie2+1,ie3  , 1) + Kloc( 4,  4)
            K(ie1+1,ie2+1,ie3  , 6) =   K(ie1+1,ie2+1,ie3  , 6) + Kloc( 4,  5)
            K(ie1+1,ie2+1,ie3  , 7) =   K(ie1+1,ie2+1,ie3  , 7) + Kloc( 4,  6)
            K(ie1+1,ie2+1,ie3  , 9) =   K(ie1+1,ie2+1,ie3  , 9) + Kloc( 4,  7)
            K(ie1+1,ie2+1,ie3  ,10) =   K(ie1+1,ie2+1,ie3  ,10) + Kloc( 4,  8)
            K(ie1  ,ie2  ,ie3+1, 1) =   K(ie1  ,ie2  ,ie3+1, 1) + Kloc( 5,  5)
            K(ie1  ,ie2  ,ie3+1, 2) =   K(ie1  ,ie2  ,ie3+1, 2) + Kloc( 5,  6)
            K(ie1  ,ie2  ,ie3+1, 4) =   K(ie1  ,ie2  ,ie3+1, 4) + Kloc( 5,  7)
            K(ie1  ,ie2  ,ie3+1, 5) =   K(ie1  ,ie2  ,ie3+1, 5) + Kloc( 5,  8)
            K(ie1+1,ie2  ,ie3+1, 1) =   K(ie1+1,ie2  ,ie3+1, 1) + Kloc( 6,  6)
            K(ie1+1,ie2  ,ie3+1, 3) =   K(ie1+1,ie2  ,ie3+1, 3) + Kloc( 6,  7)
            K(ie1+1,ie2  ,ie3+1, 4) =   K(ie1+1,ie2  ,ie3+1, 4) + Kloc( 6,  8)
            K(ie1  ,ie2+1,ie3+1, 1) =   K(ie1  ,ie2+1,ie3+1, 1) + Kloc( 7,  7)
            K(ie1  ,ie2+1,ie3+1, 2) =   K(ie1  ,ie2+1,ie3+1, 2) + Kloc( 7,  8)
            K(ie1+1,ie2+1,ie3+1, 1) =   K(ie1+1,ie2+1,ie3+1, 1) + Kloc( 8,  8)
	    
        enddo
    enddo
enddo

end subroutine ndt_assembly
end module_ndt_assembly
        

                            
                        
                        
            














