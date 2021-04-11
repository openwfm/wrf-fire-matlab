module module_coarsening

use module_io_matlab, only: crash

contains

subroutine prolongation(   &
    ifds, ifde, kfds,kfde, jfds, jfde,            & ! fire grid dimensions
    ifms, ifme, kfms,kfme, jfms, jfme,            &
    ifts, ifte, kfts, kfte, jfts,jfte,            &
    ifcds, ifcde, kfcds,kfcde, jfcds,jfcde,       & ! coarse grid domain
    ifcms, ifcme, kfcms,kfcme, jfcms,jfcme,       & ! coarse grid dimensions
    ifcts, ifcte, kfcts,kfcte, jfcts,jfcte,       & ! coarse grid tile 
    u,uc,cr_x,cr_y,cl_z,X,Y,Z)

! Multiply by the prolongation matrix
! In:
!   uc      coarse grid vector
!   cr_x, cr_y  coarsening factor in horizontal directions x and y
!   cl_z    1D array, indices of coarse grid in the z directions
!   X,Y,Z   grid coordinates 
! Out:
!   u      fine grid vector 
  
implicit none
!*** arguments

integer, intent(in)::                             & 
    ifds, ifde, kfds,kfde, jfds,jfde,             & ! fire grid dimensions
    ifms, ifme, kfms,kfme, jfms,jfme,             &
    ifts, ifte, kfts, kfte, jfts,jfte

integer, intent(in)::                             &  
    ifcds, ifcde, kfcds,kfcde, jfcds,jfcde,       & ! coarse grid domain
    ifcms, ifcme, kfcms,kfcme, jfcms,jfcme,       & ! coarse grid dimensions
    ifcts, ifcte, kfcts,kfcte, jfcts,jfcte          ! coarse grid tile 

real, intent(in), dimension(ifms:ifme,kfms:kfme,jfms:jfme):: X,Y,Z !spatial grid
real, intent(in), dimension(ifcms:ifcme,kfcms:kfcme,jfcms:jfcme):: uc ! coarse vector

integer, intent(in):: cr_x, cr_y, &       ! coarsening factors in the horizonal directions
    cl_z(kfcts:kfcte)                      ! indices of coarse grid in the vertical direction

real, intent(out), dimension(ifms:ifme,kfms:kfme,jfms:jfme):: u ! fine grid interpolant

!*** local
integer :: i,j,k,ic,jc,kc,ifc,jfc,kfc,ifs,ife,jfs,jfe,kfs,kfe
real:: qi,qj,qk

!*** executable

! zero the output to ready for contributions
do j=jfts,jfte
    do k=kfts,kfte
        do i=ifts,ifte
            u(i,k,j)=0.
        enddo
    enddo
enddo

if ((cl_z(kfcts) .ne. kfts) .or. (cl_z(kfcte) .ne. kfte)) then
    call crash('vertical corsening must include all domain')
endif

do kc=kfcts,kfcte           ! loop over coarse layers    
    kfc=cl_z(kc);          ! the fine grid number of the coarse layer
    if (kfc>kfts) then
        kfs=cl_z(kc-1)+1   ! from above previous coarse     
    else
        kfs=kfc            ! itself, there is no previous fine layer
    endif
    if (kfc<kfde) then
        kfe=cl_z(kc+1)-1   ! up to under next layer
    else
        kfe=kfc            ! itself, there is no next layer
    endif
    !print *,'vertical coarse layer ',kc,' at ',kfc,' contributes to layers ',kfs,':',kfe 
    do k=kfs,kfe
        ! really needs to be made fine point oriented and draw from coarse points
        do jc=jfcts,jfcte          
            jfc=cr_y*(jc-jfcds)+jfds ! fine grid index of the coarse point
            jfs=max(jfc-cr_y+1,jfds) ! start support
            jfe=min(jfc+cr_y-1,jfde) ! end support
            if (jfc > jfde) then  ! after end of domain not divisible by coarse ratio
               jfc = jfde
	    elseif (jfc > jfde - cr_y .and. jfc < jfde)then ! just before end of domain not divisible by coarse ratio
               jfe = jfde -1
            endif
            !print *,'coarse y ',jc,' at j ',jfc,' contributes to ',jfs,':',jfe
                do ic=ifcts,ifcte
                ifc=cr_y*(ic-ifcds)+ifds ! fine grid index of the coarse point
                ifs=max(ifc-cr_y+1,ifds) ! start support
                ife=min(ifc+cr_y-1,ifde) ! end support
                if (ifc > ifde) then  ! after end of domain not divisible by coarse ratio
                   ifc = ifde
    	        elseif (ifc > ifde - cr_x .and. ifc < ifde)then ! just before end of domain not divisible by coarse ratio
                   ife = ifde -1
                endif
                !print *,'coarse x ',ic,' at i ',ifc,' contributes to ',ifs,':',ife
                do j=jfs,jfe
                    do i=ifs,ife
                        if (i>ifc) then 
                            qi=(X(i,k,j)-X(ife+1,k,j))/(X(ifc,k,j)-X(ife+1,k,j))
                        elseif (i<ifc) then 
                            qi=(X(i,k,j)-X(ifs-1,k,j))/(X(ifc,k,j)-X(ifs-1,k,j))
                        else
                            qi=1.
                        endif
                        if (j>jfc) then 
                            qj=(Y(i,k,j)-Y(i,k,jfe+1))/(Y(i,k,jfc)-Y(i,k,jfe+1))
                        elseif (j<jfc) then 
                            qj=(Y(i,k,j)-Y(i,k,jfs-1))/(Y(i,k,jfc)-Y(i,k,jfs-1))
                        else
                            qj=1.
                        endif
                        if (k>kfc) then 
                            qk=(Z(i,k,j)-Z(i,kfe+1,j))/(Z(i,kfc,j)-Z(i,kfe+1,j))
                        elseif (k<kfc) then 
                            qk=(Z(i,k,j)-Z(i,kfs-1,j))/(Z(i,kfc,j)-Z(i,kfs-1,j))
                        else
                            qk=1.
                        endif
                        u(i,k,j) = u(i,k,j) + qi*qk*qj*uc(ic,kc,jc);
                    enddo
                enddo
            enddo
        enddo
    enddo
enddo

end subroutine prolongation

end module module_coarsening

