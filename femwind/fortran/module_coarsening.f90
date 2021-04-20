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

subroutine coarsening_icl(cr_x,cr_y,icl_z,dx,dy,dz,A,minaspect,maxaspect)
! decide on coarsening
! in:
!   dx,dy       mesh spacings, scalar
!   dz          verticl element size, vector
!   A           matrix size (3,3), only diagonal used
! out:
!   cr_x,cr_y   horizontal coarsening factors in directions 1 and 2
!   icl_z       coarse indices in direction 3
implicit none

!*** arguments
real,intent(in):: dx,dy,dz(:),A(:,:)
real,intent(in):: minaspect,maxaspect
integer,intent(out):: cr_x,cr_y
integer,pointer,intent(out):: icl_z(:)

!*** local
real:: dxy,crit,hzcavg,arat
integer, allocatable:: icl3(:)
integer:: nc3,newlcl,lcl,i,n3


!*** executable
    dxy=min(dx,dy)  ! horizontal step
    n3 = size(dz)+1 !
    arat = A(3,3)/min(A(1,1),A(2,2))  ! scaled vertical penalty
    ! decide on horizontal coarsening factors
    crit=(dz(1)/dxy)/arat
    if (crit > minaspect) then
        cr_x = 2 
        cr_y = 2 
    else
        cr_x = 1 
        cr_y = 1 
    endif
    hzcavg=sqrt(real(cr_x*cr_y)); 
    print *,'horizontal coarsening factors ',cr_x,cr_y, ' because weighted height is ', crit
    allocate(icl3(n3+1))
    icl3=0
    lcl=1 ! last coarse level
    icl3(1)=lcl
    nc3=0
    do i=1,n3
        newlcl=lcl+1  ! next coarse level by 1
        if (lcl+2 <= n3) then
            crit = ((dz(lcl)+dz(lcl+1))/(dxy*hzcavg))/arat
            if (crit < maxaspect ) then
                newlcl=lcl+2 ! next coarse level by 2
            endif
        endif
        lcl = newlcl;
        if (lcl <= n3) then
            icl3(i+1)=lcl
        else ! at the top already
            nc3=i
            allocate(icl_z(i))
            icl_z = icl3(1:i)
            exit
        endif
    enddo     
    deallocate(icl3)
    if (nc3==0) then
        call crash('coarsening_icl: number of coarse layers is 0')
    endif
    print *,'vertical coarse layers ',icl_z

end subroutine coarsening_icl

subroutine hzc2icl(nc_x, nc_y, icl_x, icl_y, cr_x, cr_y, n_x, n_y)
! compute horizontal coarse index vectors 
!
! In:
!   cr_x, cr_y coarsening factors in x and y
!   n_x, n_y   fine mesh size in x and y
! Out:
!   nc_x, nc_y      coarse mesh size
!   icl_x, icl_y    fine indices of coarse grid lines

implicit none

!*** arguments
integer, intent(in)::cr_x, cr_y, n_x, n_y
integer, intent(out):: nc_x, nc_y
integer, intent(out), pointer, dimension(:) :: icl_x, icl_y

!*** local
integer::i

nc_x = (n_x + 2)/2
nc_y = (n_y + 2)/2

allocate(icl_x(nc_x))
allocate(icl_y(nc_y))

do i=1,n_x,2
    icl_x(i)=1+2*(i-1)  ! every second node
enddo
icl_x(nc_x)=n_x ! last coarse is last fine

do i=1,n_y,2
    icl_y(i)=1+2*(i-1)  ! every second node
enddo
icl_y(nc_y)=n_y ! last coarse is last fine

end subroutine hzc2icl

end module module_coarsening

