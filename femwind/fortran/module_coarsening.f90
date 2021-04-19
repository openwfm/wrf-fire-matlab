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
    u,uc,cr_x,cr_y,icl_z,X,Y,Z)

! Multiply by the prolongation matrix
! In:
!   uc      coarse grid vector
!   cr_x, cr_y  coarsening factor in horizontal directions x and y
!   icl_z    1D array, indices of coarse grid in the z directions
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
    icl_z(kfcts:kfcte)                      ! indices of coarse grid in the vertical direction

real, intent(out), dimension(ifms:ifme,kfms:kfme,jfms:jfme):: u ! fine grid interpolant

!*** local
integer :: i,j,k,ic,jc,kc,ifc,jfc,kfc,ifs,ife,jfs,jfe,kfs,kfe
real:: qi,qj,qk

!*** executable

! zero the output to ready for contributions
u(ifts:ifte,kfts:kfte,jfts:jfte) = 0.

if ((icl_z(kfcts) .ne. kfts) .or. (icl_z(kfcte) .ne. kfte)) then
    call crash('vertical corsening must include all domain')
endif

do kc=kfcts,kfcte           ! loop over coarse layers    
    kfc=icl_z(kc);          ! the fine grid number of the coarse layer
    if (kfc>kfts) then
        kfs=icl_z(kc-1)+1   ! from above previous coarse     
    else
        kfs=kfc            ! itself, there is no previous fine layer
    endif
    if (kfc<kfde) then
        kfe=icl_z(kc+1)-1   ! up to under next layer
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

subroutine restriction(   &
    ifds, ifde, kfds,kfde, jfds, jfde,            & ! fire grid dimensions
    ifms, ifme, kfms,kfme, jfms, jfme,            &
    ifts, ifte, kfts, kfte, jfts,jfte,            &
    ifcds, ifcde, kfcds,kfcde, jfcds,jfcde,       & ! coarse grid domain
    ifcms, ifcme, kfcms,kfcme, jfcms,jfcme,       & ! coarse grid dimensions
    ifcts, ifcte, kfcts,kfcte, jfcts,jfcte,       & ! coarse grid tile 
    uc,u,cr_x,cr_y,icl_z,X,Y,Z)

! Multiply by the prolongation matrix transpose
! In:
!   u      fine grid vector 
!   cr_x, cr_y  coarsening factor in horizontal directions x and y
!   icl_z    1D array, indices of coarse grid in the z directions
!   X,Y,Z   grid coordinates 
! Out:
!   uc      coarse grid vector
  
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
real, intent(out), dimension(ifcms:ifcme,kfcms:kfcme,jfcms:jfcme):: uc ! coarse vector

integer, intent(in):: cr_x, cr_y, &       ! coarsening factors in the horizonal directions
    icl_z(kfcts:kfcte)                      ! indices of coarse grid in the vertical direction

real, intent(in), dimension(ifms:ifme,kfms:kfme,jfms:jfme):: u ! fine grid interpolant

!*** local
integer :: i,j,k,ic,jc,kc,ifc,jfc,kfc,ifs,ife,jfs,jfe,kfs,kfe
real:: qi,qj,qk

!*** executable

! zero the output to ready for contributions
uc(ifcts:ifcte,kfcts:kfcte,jfcts:jfcte) = 0.

if ((icl_z(kfcts) .ne. kfts) .or. (icl_z(kfcte) .ne. kfte)) then
    call crash('vertical corsening must include all domain')
endif

do kc=kfcts,kfcte           ! loop over coarse layers    
    kfc=icl_z(kc);          ! the fine grid number of the coarse layer
    if (kfc>kfts) then
        kfs=icl_z(kc-1)+1   ! from above previous coarse     
    else
        kfs=kfc            ! itself, there is no previous fine layer
    endif
    if (kfc<kfde) then
        kfe=icl_z(kc+1)-1   ! up to under next layer
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
                        uc(ic,kc,jc) = uc(ic,kc,jc) + u(i,k,j)*qi*qk*qj
                    enddo
                enddo
            enddo
        enddo
    enddo
enddo

end subroutine restriction

subroutine coarsening_icl(cr_x,cr_y,icl_z,dx,dy,dz)
! decide on coarsening
! in:
!   dx,dy       mesh spacings, scalar
!   dz          verticl element size, vector
!   params      structure
! out:
!   hzc         horizontal coarsening factors in directions 1 and 2
!   icl_z        coarse indices in direction 3
implicit none

!*** arguments
real,intent(in)::dx,dy,dz(:)
real,intent(out)::cr_x,cr_y
real,pointer,intent(out):icl_z(:)

!*** local


!*** executable
    if ~isvector(dz),
          error('dz must be a vector')
    end
    dz = dz(:)';  % make sure dz is a row
    % add 0 to htt bottom
    dxy=min(dx,dy);  % horizontal step
    n3 = length(dz)+1;
    % decide on horizontal coarsening factor
    crit=(dz(1)/dxy)/params.a(3);
    if crit > params.minaspect
        hzc=[2,2]; % future proofing if they are different 
    else
        hzc=[1,1];
    end
    hzcavg=sqrt(hzc(1)*hzc(2)); 
    fprintf('horizontal coarsening factor %g %g because weighted height=%g\n',...
        hzc, crit)
    icl3=zeros(1,n3); % allocate max
    lcl=1; % last coarse level
    icl3(1)=lcl;
    nc3=0;
    for i=1:n3
        newlcl=lcl+1; % next coarse level by 1
        if lcl+2 <= n3 
            crit = ((dz(lcl)+dz(lcl+1))/(2*dxy*hzcavg/2))/params.a(3);
            if crit < params.maxaspect  
                newlcl=lcl+2; % next coarse level by 2
            end
        end
        lcl = newlcl;
        if lcl <= n3
            icl3(i+1)=lcl;
        else % at the top already
            nc3=i;
            icl3 = icl3(1:i);
            break
        end
    end     
    if nc3==0
        error('number of coarse layers is 0')
    end
    disp(['vertical coarse layers ',num2str(icl3)])
    hg=[0,cumsum(dz)];
    hgc=hg(icl3);
    disp(['heights above terrain ',num2str(hg)])
    disp(['coarse heights above terrain ',num2str(hgc)])
end


! 
end subroutine coarsening_icl
end module module_coarsening

