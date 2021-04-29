! **** TESTING ***
! called from femwind_test.f90 compiled as femwind_test.exe
! femwind_test.exe is executed from femwind_fortran.m if params.run_fortran
! femwind_fortran.m is called from femwind_main.m   
! run femwind_main from e.g. femwind_rate_test or femwind_wrfout_test



module module_femwind
       use module_utils
       use module_f_assembly
       use module_ndt_assembly
       use module_w_assembly
       use module_coarsening
       use module_boundary_conditions

! parameters
integer, parameter::max_levels=20
integer, parameter::msize=14
real:: minaspect=1./3.,maxaspect=3.
real::A(3,3)=reshape((/1., 0., 0.,  0., 1., 0.,  0., 0., 1./),(/3, 3/))

integer::nlevels                                ! number of levels

type mg_type
    real, dimension(3,3):: A                        ! penalty weight matrix
    real, pointer, dimension(:,:,:):: X, Y, Z       ! grid vertices
    real, pointer, dimension(:,:,:,:):: Kglo        ! global stiffness matrix
    
    real:: dx,dy                                    ! horizontal spacing, scalar 
    real, pointer, dimension(:)::dz                 ! vertical spacing of the layers
    integer::nx,ny,nz,nn                            ! mesh size in vertices
  
    integer::level

    integer:: cr_x, cr_y                            ! coarsening factors
    integer, pointer, dimension(:):: icl_x, icl_y, icl_z  ! coarsening indices

    integer:: &
    ifds, ifde, kfds, kfde, jfds, jfde,           & ! fire grid dimensions
    ifms, ifme, kfms, kfme, jfms, jfme,           & ! memory dimensions
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts, jfte              ! tile dimensions
 
end type
       
contains

subroutine get_mg_dims(mg_level, &
    ifds, ifde, kfds, kfde, jfds, jfde,           & ! fire grid dimensions
    ifms, ifme, kfms, kfme, jfms, jfme,           & ! memory dimensions
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts, jfte)              ! tile dimensions
implicit none
type(mg_type),intent(in)::mg_level
integer, intent(out)::  &
    ifds, ifde, kfds, kfde, jfds, jfde,           & ! fire grid dimensions
    ifms, ifme, kfms, kfme, jfms, jfme,           & ! memory dimensions
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts, jfte              ! tile dimensions
 
    ifds = 1
    ifde = mg_level%nx-1  ! dimension in elements
    jfds = 1
    jfde = mg_level%ny-1 
    kfds = 1
    kfde = mg_level%nz-1

    ifts = ifds
    ifte = ifde
    jfts = jfds
    jfte = jfde
    kfts = kfds
    kfte = kfde

    ifps = ifds
    ifpe = ifde
    jfps = jfds
    jfpe = jfde
    kfps = kfds
    kfpe = kfde

    ifms = ifps - 1
    ifme = ifpe + 2  ! need +1 for vertex grid, and +1 for zero strip 
    jfms = jfps - 1
    jfme = jfpe + 2
    kfms = kfps - 1
    kfme = kfpe + 2

end subroutine get_mg_dims

subroutine write_3d(mg,l,var,name,v)
! write array in text file for debugging 
implicit none

!*** arguments
type(mg_type),intent(in)::mg(:)    !multigrid structure
integer,intent(in)::l              !level
real, intent(in)::var(:,:,:)       !the array
character(len=*),intent(in)::name  
character(len=1),intent(in)::v     !'v' or 'e'

!*** local
integer:: &
    ifds, ifde, kfds, kfde, jfds, jfde,           & ! fire grid dimensions
    ifms, ifme, kfms, kfme, jfms, jfme,           & ! memory dimensions
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts, jfte              ! tile dimensions
integer::ie,ke,je
character(len=2)::lc

!*** executable

call get_mg_dims(mg(l), &
    ifds, ifde, kfds,kfde, jfds, jfde,            & ! fire grid dimensions
    ifms, ifme, kfms,kfme, jfms, jfme,            &
    ifps, ifpe, kfps,kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts,kfte, jfts,jfte)            

write(lc,'(i2.2)')l

if(v.eq.'v')then
    ie = snode(ifte,ifde,+1)  ! vertex based
    je = snode(jfte,jfde,+1)
    ke = snode(kfte,kfde,+1)
elseif(v.eq.'e')then
    ie=ifte
    je=jfte
    ke=kfte
else
    call crash('write_3d: bad v')
endif
 
call write_tight(var) ! write dereferenced

contains

subroutine write_tight(varr)
    real, intent(in)::varr(ifms:ifme,kfms:kfme,jfms:jfme)
    ! for element/cell/midpoint based arrays, like fmw winds
    ! open(iu,file='varr.dat',form='unformatted',status='unknown')
    ! write(iu)varr(ifts:ifte,kfts:kfte,jfts:jfte)
    ! close(iu)
    call write_array(varr(ifts:ie,kfts:ke,jfts:je),name//lc)
end subroutine write_tight

end subroutine write_3d

subroutine print_mg_dims(mg,l)
type(mg_type),intent(in)::mg(:)  ! multigrid level
    print *,'level ',mg(l)%level
    if(associated(mg(l)%X))print *,'X memory shape ',shape(mg(l)%X)
    if(associated(mg(l)%Y))print *,'Y memory shape ',shape(mg(l)%Y)
    if(associated(mg(l)%Z))print *,'Z memory shape ',shape(mg(l)%Z)
    print *,'dx=',mg(l)%dx,' dy=',mg(l)%dy
    if(associated(mg(l)%dz))print *,'dz=',mg(l)%dz,' shape ',shape(mg(l)%dz)
    print *,'grid size nx=',mg(l)%nx,' ny=',mg(l)%ny,' nz=',mg(l)%nz, &
       ' total ndof=',mg(l)%nn
end subroutine print_mg_dims


subroutine femwind_setup(mg)
implicit none
! set up the mg_level structure
! input: mg(1)%X,Y,Z,dx,dy,dz already set
!*** arguments
type(mg_type),intent(inout)::mg(:)  ! multigrid level
!*** local

integer::k,l,m,nzc
real::s
integer::   &
    ifds, ifde, kfds, kfde, jfds, jfde,           & ! fire grid dimensions
    ifms, ifme, kfms, kfme, jfms, jfme,           & ! memory dimensions
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts, jfte,           & ! tile dimensions
    ifcds, ifcde, kfcds,kfcde, jfcds,jfcde,       & ! coarse grid domain
    ifcms, ifcme, kfcms,kfcme, jfcms,jfcme,       & ! coarse grid dimensions
    ifcps, ifcpe, kfcps,kfcpe, jfcps,jfcpe,       & ! coarse grid dimensions
    ifcts, ifcte, kfcts,kfcte, jfcts,jfcte          ! coarse grid tile


!*** executable

! decide on the grid coarsening and compute scalars and 1D index arrays

    mg(1)%nn = mg(1)%nx *  mg(1)%ny * mg(1)%nz
    mg(1)%level = 1

    print *,'femwind_setup received'
    call print_mg_dims(mg,1)
    
    nlevels = max_levels-1


    do l=1,nlevels

        ! get horizontal coarsening ratios and vertical coarse list

        print *,'multigrid level ',l,' grid size ', mg(l)%nx, mg(l)%ny, mg(l)%nz
        ! decide about the next level
        call coarsening_icl(mg(l)%cr_x,mg(l)%cr_y,mg(l)%icl_z, &
                            mg(l)%dx,mg(l)%dy,mg(l)%dz,A,minaspect,maxaspect)

        ! get horizontal coarsening lists

        call coarsening_hzc2icl(mg(l)%icl_x, mg(l)%icl_y, &
                                mg(l)%cr_x,mg(l)%cr_y,&
                                mg(l)%nx, mg(l)%ny)

        ! update coarse mesh scalars

        mg(l+1)%nx = size(mg(l)%icl_x)
        mg(l+1)%ny = size(mg(l)%icl_y)
        mg(l+1)%nz = size(mg(l)%icl_z)
        mg(l+1)%nn = mg(l+1)%nx *  mg(l+1)%ny * mg(l+1)%nz
        mg(l+1)%dx = mg(l)%dx/mg(l)%cr_x
        mg(l+1)%dy = mg(l)%dy/mg(l)%cr_y
        mg(l+1)%level = l+1

        call print_mg_dims(mg,l+1)

        if (mg(l+1)%nn >= mg(l)%nn .or. l == max_levels) then
            print *,'stopping at ',l,' levels because coarse ndof ',mg(l+1)%nn, &
                '>= fine ', mg(l)%nn,' or l= ',l,' = ',max_levels 
            nlevels = l 
            exit ! the loop
        endif

        ! coarsen dz

        allocate(mg(l+1)%dz(mg(l+1)%nz-1))
	do k=1,mg(l+1)%nz - 1 
            ! print *,'computing coarse dz of layer ',k
            s = 0.
            do m=mg(l)%icl_z(k),mg(l)%icl_z(k+1)-1
                ! print *,'adding fine dz ',m,' equal ',mg(l)%dz(m)
                s = s + mg(l)%dz(m)
            enddo
            mg(l+1)%dz(k)=s
            ! print *,'height dz of coarse layer',k,' is ',s
        enddo
			

    enddo
    print *,nlevels,' levels'

    ! get coarse grid coordinates 
    do l=1,nlevels-1 
        call get_mg_dims(mg(l), &
            ifds, ifde, kfds,kfde, jfds, jfde,            & ! fire grid dimensions
            ifms, ifme, kfms,kfme, jfms, jfme,            &
            ifps, ifpe, kfps,kfpe, jfps, jfpe,           & ! fire patch bounds
            ifts, ifte, kfts,kfte, jfts,jfte)            
        call get_mg_dims(mg(l+1), &
            ifcds, ifcde, kfcds,kfcde, jfcds, jfcde,            & ! fire grid dimensions
            ifcms, ifcme, kfcms,kfcme, jfcms, jfcme,            &
            ifcps, ifcpe, kfcps,kfcpe, jfcps, jfcpe,           & ! fire patch bounds
            ifcts, ifcte, kfcts,kfcte, jfcts,jfcte)            
        allocate(mg(l+1)%X(ifms: ifme, kfms: kfme, jfms: jfme))
        allocate(mg(l+1)%Y(ifms: ifme, kfms: kfme, jfms: jfme))
        allocate(mg(l+1)%Z(ifms: ifme, kfms: kfme, jfms: jfme))
        call coarsening_grid(l, &
            ifds, ifde, kfds, kfde, jfds, jfde,           & ! fire grid dimensions
            ifms, ifme, kfms, kfme, jfms, jfme,           & ! memory dimensions
            ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
            ifts, ifte, kfts, kfte, jfts, jfte,           & ! tile dimensions
            ifcds, ifcde, kfcds,kfcde, jfcds,jfcde,       & ! coarse grid domain
            ifcms, ifcme, kfcms,kfcme, jfcms,jfcme,       & ! coarse grid dimensions
            ifcps, ifcpe, kfcps,kfcpe, jfcps,jfcpe,       & ! coarse grid dimensions
            ifcts, ifcte, kfcts,kfcte, jfcts,jfcte,       & ! coarse grid tile
            mg(l)%icl_x, mg(l)%icl_y, mg(l)%icl_z,        & ! indices of coarse grid
            mg(l)%X, mg(l)%Y, mg(l)%Z,                    & ! fine grid coordinates
            mg(l+1)%X, mg(l+1)%Y, mg(l+1)%Z)                ! coarse grid coordinates
    enddo
 
    do l=1,nlevels
        call write_3d(mg,l,mg(l)%X,'X','v')
        call write_3d(mg,l,mg(l)%Y,'Y','v')
        call write_3d(mg,l,mg(l)%Z,'Z','v')
    enddo


    ! assemble the stiffness matrices 
    do l=1,nlevels 
        call get_mg_dims(mg(l), &
            ifds, ifde, kfds,kfde, jfds, jfde,            & ! fire grid dimensions
            ifms, ifme, kfms,kfme, jfms, jfme,            &
            ifps, ifpe, kfps,kfpe, jfps, jfpe,           & ! fire patch bounds
            ifts, ifte, kfts,kfte, jfts,jfte)            
        allocate(mg(l)%Kglo(ifms:ifme, kfms:kfme, jfms:jfme, msize))
        call ndt_assembly(                              &
            ifds, ifde, kfds,kfde, jfds, jfde,            & ! fire grid dimensions
            ifms, ifme, kfms,kfme, jfms, jfme,            &
            ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
            ifts, ifte, kfts, kfte, jfts,jfte,            &
            A, mg(l)%X,mg(l)%Y,mg(l)%Z, 1,                & 
            mg(l)%Kglo)        
    enddo

end subroutine femwind_setup

subroutine femwind_solve(mg,      &
     ifds, ifde, kfds,kfde, jfds, jfde,            & ! fire grid dimensions
     ifms, ifme, kfms,kfme, jfms, jfme,            &
     ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
     ifts, ifte, kfts, kfte, jfts,jfte,            &
     u0, v0, w0, u, v, w, rate)                  ! outputs 

implicit none

!*** arguments

type(mg_type),intent(in)::mg(:)  ! multigrid levels

integer, intent(in)::                             &
    ifds, ifde, kfds,kfde, jfds, jfde,            & ! fire grid dimensions
    ifms, ifme, kfms,kfme, jfms, jfme,            &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte            

real, intent(in), dimension(ifms:ifme, kfms:kfme, jfms:jfme):: u0, v0, w0   !initial wind vector at midpoints
real,intent(out), dimension(ifms:ifme, kfms:kfme, jfms:jfme)::u, v, w       !mass consistent wind at midpoints
real,intent(out):: rate

!*** local
real, dimension(ifms:ifme, kfms:kfme, jfms:jfme):: f  
real, dimension(ifts:ifte+1, kfts:kfte+1, jfts:jfte+1):: f_f ! testing only

!*** executable

! F = div(u0)
! F = f_assembly_fortran(A,X,U0,lambda,params);
call f_assembly(                        &
    ifds, ifde, kfds, kfde, jfds, jfde,           & ! fire grid dimensions
    ifms, ifme, kfms, kfme, jfms, jfme,           &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte,            &
    A, mg(1)%X, mg(1)%Y, mg(1)%Z, u0, v0, w0,                       & !Input from femwind, U0, V0, W0 not used in hexa to construct Kloc or JG
    F)                                             !Global load vector output  

call vec_boundary_conditions(                              &
    ifds, ifde, kfds,kfde, jfds, jfde,                       & ! fire grid dimensions
    ifms, ifme, kfms,kfme, jfms, jfme,            &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte,             &
    F)


f_f=f(ifts:ifte+1, kfts:kfte+1, jfts:jfte+1) ! testing only
call write_array(f_f,'F_f') ! testing only

! initialize for now
u=0.
v=0.
w=0.
rate =0.

end subroutine femwind_solve
end module  module_femwind
