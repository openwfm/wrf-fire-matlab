module module_ndt_assembly
       use module_hexa
contains
subroutine ndt_assembly(                              &
    ifds, ifde, kfds,kfde, jfds, jfde,            & ! fire grid dimensions
    ifms, ifme, kfms,kfme, jfms, jfme,            &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte,            &
    A, X,Y,Z, iflags,                   	  & !Input from femwind, U0, V0, W0 not used in hexa to construct Kloc or JG
    K)                                             !Global stiffness matrix output  

implicit none


!*** arguments

integer, intent(in)::    			  &
    ifds, ifde, kfds,kfde, jfds, jfde,            & ! fire grid dimensions
    ifms, ifme, kfms,kfme, jfms, jfme,            &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte            




integer, parameter:: msize=14
real, intent(in), dimension(3,3):: A
real, intent(in), dimension(ifms:ifme, kfms:kfme, jfms: jfme):: X,Y,Z!spatial grid
!Input for hexa
integer, intent(in), dimension(3,1,1)::iflags

real, intent(out), dimension(ifms:ifme, kfms:kfme, jfms:jfme,1:msize)::K


!*** local

integer:: ie1, ie2, ie3, ic1, ic2, ic3, iloc
real ::  Kloc(8,8), Floc(8), Jg(8,3)
real ::  Xloc(3,8), u0loc(3)
    	
!*** u0loc is an input for module_hexa, but is not used to construct K. Do I need to define this?
!*** integer, dimension(3,1,1), save ::iflags = reshape((/1,0,1/),(/3,1,1/)) !define iflags to construct JG and Kloc in hexa



!** executable
Xloc = 99999.
K = 0.
!print *, 'A=', A(:,:)
!print *, 'X=(:,1,:)', X(:,1,:)
!print *, 'X=(:,2,:)', X(:,2,:)
!print *, 'X=(:,3,:)', X(:,3,:)

do ie2=jfts,jfte -1
    do ie3=kfts, kfte -1
        do ie1=ifts, ifte -1
            do ic2=0,1
                do ic3=0,1
                    do ic1=0,1
                        iloc=1+ic1+2*(ic2+2*ic3);  !local index of the node in the element
                            Xloc(1,iloc)=X(ie1 + ic1, ie3 + ic3, ie2 + ic2)
                            Xloc(2,iloc)=Y(ie1 + ic1, ie3 + ic3, ie2 + ic2)
                            Xloc(3,iloc)=Z(ie1 + ic1, ie3 + ic3, ie2 + ic2)
				
                    enddo
                enddo
            enddo
	    print *,'Xloc(1,8)=', Xloc(:,:)
            call hexa(A,Xloc,u0loc,Kloc,Floc,Jg,iflags)
	    !print *, 'Kloc(1,1)=',Kloc(1,1)
	    !print *, ie1, ie2, ie3, K(1,1,2,1)
            K(ie1  ,ie3  ,ie2  , 1) =   K(ie1  ,ie3  ,ie2  , 1) + Kloc( 1,  1) 
            K(ie1  ,ie3  ,ie2  , 2) =   K(ie1  ,ie3  ,ie2  , 2) + Kloc( 1,  2) 
            K(ie1  ,ie3  ,ie2  , 4) =   K(ie1  ,ie3  ,ie2  , 4) + Kloc( 1,  3) 
            K(ie1  ,ie3  ,ie2  , 5) =   K(ie1  ,ie3  ,ie2  , 5) + Kloc( 1,  4) 
            K(ie1  ,ie3  ,ie2  ,10) =   K(ie1  ,ie3  ,ie2  ,10) + Kloc( 1,  5) 
            K(ie1  ,ie3  ,ie2  ,11) =   K(ie1  ,ie3  ,ie2  ,11) + Kloc( 1,  6) 
            K(ie1  ,ie3  ,ie2  ,13) =   K(ie1  ,ie3  ,ie2  ,13) + Kloc( 1,  7) 
            K(ie1  ,ie3  ,ie2  ,14) =   K(ie1  ,ie3  ,ie2  ,14) + Kloc( 1,  8) 
            K(ie1+1,ie3  ,ie2  , 1) =   K(ie1+1,ie3  ,ie2  , 1) + Kloc( 2,  2) 
            K(ie1+1,ie3  ,ie2  , 3) =   K(ie1+1,ie3  ,ie2  , 3) + Kloc( 2,  3) 
            K(ie1+1,ie3  ,ie2  , 4) =   K(ie1+1,ie3  ,ie2  , 4) + Kloc( 2,  4) 
            K(ie1+1,ie3  ,ie2  , 9) =   K(ie1+1,ie3  ,ie2  , 9) + Kloc( 2,  5) 
            K(ie1+1,ie3  ,ie2  ,10) =   K(ie1+1,ie3  ,ie2  ,10) + Kloc( 2,  6) 
            K(ie1+1,ie3  ,ie2  ,12) =   K(ie1+1,ie3  ,ie2  ,12) + Kloc( 2,  7) 
            K(ie1+1,ie3  ,ie2  ,13) =   K(ie1+1,ie3  ,ie2  ,13) + Kloc( 2,  8) 
            K(ie1  ,ie3  ,ie2+1, 1) =   K(ie1  ,ie3  ,ie2+1, 1) + Kloc( 3,  3) 
            K(ie1  ,ie3  ,ie2+1, 2) =   K(ie1  ,ie3  ,ie2+1, 2) + Kloc( 3,  4) 
            K(ie1  ,ie3  ,ie2+1, 7) =   K(ie1  ,ie3  ,ie2+1, 7) + Kloc( 3,  5) 
            K(ie1  ,ie3  ,ie2+1, 8) =   K(ie1  ,ie3  ,ie2+1, 8) + Kloc( 3,  6) 
	 K(ie1  ,ie3  ,ie2+1,10) =   K(ie1  ,ie3  ,ie2+1,10) + Kloc( 3,  7) 
	 K(ie1  ,ie3  ,ie2+1,11) =   K(ie1  ,ie3  ,ie2+1,11) + Kloc( 3,  8) 
	 K(ie1+1,ie3  ,ie2+1, 1) =   K(ie1+1,ie3  ,ie2+1, 1) + Kloc( 4,  4) 
	 K(ie1+1,ie3  ,ie2+1, 6) =   K(ie1+1,ie3  ,ie2+1, 6) + Kloc( 4,  5) 
	 K(ie1+1,ie3  ,ie2+1, 7) =   K(ie1+1,ie3  ,ie2+1, 7) + Kloc( 4,  6) 
	 K(ie1+1,ie3  ,ie2+1, 9) =   K(ie1+1,ie3  ,ie2+1, 9) + Kloc( 4,  7) 
	 K(ie1+1,ie3  ,ie2+1,10) =   K(ie1+1,ie3  ,ie2+1,10) + Kloc( 4,  8) 
	 K(ie1  ,ie3+1,ie2  , 1) =   K(ie1  ,ie3+1,ie2  , 1) + Kloc( 5,  5) 
	 K(ie1  ,ie3+1,ie2  , 2) =   K(ie1  ,ie3+1,ie2  , 2) + Kloc( 5,  6) 
	 K(ie1  ,ie3+1,ie2  , 4) =   K(ie1  ,ie3+1,ie2  , 4) + Kloc( 5,  7) 
	 K(ie1  ,ie3+1,ie2  , 5) =   K(ie1  ,ie3+1,ie2  , 5) + Kloc( 5,  8) 
	 K(ie1+1,ie3+1,ie2  , 1) =   K(ie1+1,ie3+1,ie2  , 1) + Kloc( 6,  6) 
	 K(ie1+1,ie3+1,ie2  , 3) =   K(ie1+1,ie3+1,ie2  , 3) + Kloc( 6,  7) 
	 K(ie1+1,ie3+1,ie2  , 4) =   K(ie1+1,ie3+1,ie2  , 4) + Kloc( 6,  8) 
	 K(ie1  ,ie3+1,ie2+1, 1) =   K(ie1  ,ie3+1,ie2+1, 1) + Kloc( 7,  7) 
	 K(ie1  ,ie3+1,ie2+1, 2) =   K(ie1  ,ie3+1,ie2+1, 2) + Kloc( 7,  8) 
	 K(ie1+1,ie3+1,ie2+1, 1) =   K(ie1+1,ie3+1,ie2+1, 1) + Kloc( 8,  8) 
	    !print *, ie1, ie2, ie3, K(1,1,2,1)	    

        enddo
    enddo
enddo

end subroutine ndt_assembly
end module  module_ndt_assembly
        

                            
                        
                        
            














