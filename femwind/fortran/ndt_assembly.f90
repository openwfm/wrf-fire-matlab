module module_ndt_assembly_test
contains
subroutine ndt_assembly(                              &
    ifds, ifde, kfds,kfde, jfds, jfde,            & ! fire grid dimensions
    ifms, ifme, kfms,kfme, jfms, jfme,            &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte,            &
    a, x, msize,                                      & !Input from femwind
    k)                                            & !Global stiffness matrix output  

implicit none


!*** arguments

integer, intent(in)::
    ifds, ifde, kfds,kfde, jfds, jfde,                       & ! fire grid dimensions
    ifms, ifme, kfms,kfme, jfms, jfme,            &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte,            &



real, intent(in), dimension(3,3):: a
real, intent(in), dimension(3,ifms:ifme, kfms:kfme, jfms: jfme)::x !spatial coordinate  matrix 
integer, parameter:: msize = 14
integer, intent(in)::nn                             !product of fire domain bounds

real, intent(in), dimension(3,8)::xloc
real, intent(in), dimension(nn,1)::fmat
real, dimension
real, intent(in), dimension(3,1)::z


real, intent(out), dimension(ifds:ifde, kfds:kfde, jfds:jfde,msize)::kmat


!*** local

integer:: ie1, ie2, ie3, ic1, ic2, ic3, iloc, i

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
                            xloc(i,iloc)=x(i,ie1 + ic1, ie2 + ic2, ie3 + ic3)
                        enddo
                    enddo
                enddo
            enddo
            call hexa(A,X,u0,Kloc,Floc,Jg,iflags)
            for ic3=0:1
                for ic2=0:1
                    for ic1=0:!
                        for kc3=0:1
                            for kc2=0:1
                                for kc1=0:1
                                    iloc=1+ic1+2*(ic2+2*ic3); % index in the local element matrix
                                    kloc=1+kc1+2*(kc2+2*kc3); % index in the local element matrix
                                    % global index triple of node i
                                    i1=ie1+ic1;
                                    i2=ie2+ic2;
                                    i3=ie3+ic3;
                                    % relative position of k vs. i
                                    j1 = kc1-ic1;
                                    j2 = kc2-ic2;
                                    j3 = kc3-ic3;
                                    % index triple of row m of K where the entry (i,k) is stored
                                    % in fortran we won't have the 2+ because
                                    % the array t will be indexed -1:1
                                    % storing in this row only 
                                    m3 = i3+t(3,2+j1,2+j2,2+j3);
                                    m2 = i2+t(2,2+j1,2+j2,2+j3);
                                    m1 = i1+t(1,2+j1,2+j2,2+j3);
                                    % index of the matrix entry (i,k) in K(m,:) 
                                    jx=     t(4,2+j1,2+j2,2+j3);
                                    % add entry of the local matrix 
                                    % this row only, no duplicates if triangle
                                    if m1==i1 && m2 == i2 && m3 == i3
                                           K(m1,m2,m3,jx) = K(m1,m2,m3,jx) + Kloc(iloc,kloc);
                                 enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo            
        enddo
    enddo
enddo
!Creating outlines of ndt_storage_table subroutine.
subroutine ndt_storage_table(m,t)
real, intent(in), dimension()::t
integer, intent(in)::j1, j2, j3

    select case(m)
        case(27)

        case(14)


end
real function g(j1, j2, j3)
g = 1 + (j1+1) + 3*(j2+1+3*(j3+1))
return
end


        

                            
                        
                        
            














