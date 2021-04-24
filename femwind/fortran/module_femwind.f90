module module_femwind
       use module_f_assembly
       use module_ndt_assembly
       use module_w_assembly
       use module_coarsening

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
    
    real:: dx,dy  
    real, pointer, dimension(:)::dz                  ! spacing of the layers
    integer::nx,ny,nz,nn                            ! mesh size in vertices

    integer:: cr_x, cr_y                            ! coarsening factors
    integer, pointer:: icl_z(:)                     ! coarsening indices in z

    integer:: &
    ifds, ifde, kfds, kfde, jfds, jfde,           & ! fire grid dimensions
    ifms, ifme, kfms, kfme, jfms, jfme,           & ! memory dimensions
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts, jfte              ! tile dimensions
 
end type

type(mg_type):: mg(max_levels)
       
contains

subroutine get_mg_dims(mg, &
    ifds, ifde, kfds, kfde, jfds, jfde,           & ! fire grid dimensions
    ifms, ifme, kfms, kfme, jfms, jfme,           & ! memory dimensions
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts, jfte)              ! tile dimensions
implicit none
type(mg_type),intent(in)::mg
integer, intent(out)::  &
    ifds, ifde, kfds, kfde, jfds, jfde,           & ! fire grid dimensions
    ifms, ifme, kfms, kfme, jfms, jfme,           & ! memory dimensions
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts, jfte              ! tile dimensions
 
    ifds = mg%ifds
    ifde = mg%ifde
    jfds = mg%jfds
    jfde = mg%jfde
    call crash('get_mg_dims not finished')

end subroutine get_mg_dims


subroutine femwind_setup(mg)
! set up the mg_level structure
! input: mg(1)%X,Y,Z,dx,dy,dz already set

type(mg_type),intent(inout)::mg(max_levels)  ! multigrid level
integer ::l,                                      & ! level
    ifds, ifde, kfds,kfde, jfds, jfde,            & ! fire grid dimensions
    ifms, ifme, kfms,kfme, jfms, jfme,            &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte            

!*** decide on the grid coarsening
    mg(1)%nx=size(mg(1)%X,1)
    mg(1)%ny=size(mg(1)%X,3)
    mg(1)%nz=size(mg(1)%X,2)
    mg(1)%nn =  mg(1)%nx *  mg(1)%ny *  mg(1)%nz
    nlevels = max_levels-1
    do l=1,nlevels
        ! get horizontal coarsening ratios and vertical coarse list
        print *,'multigrid level ',l,' grid size ', mg(l)%nx, mg(l)%ny, mg(l)%nz
        allocate(mg(l)%Kglo(ifms:ifme, kfms:kfme, jfms:jfme, msize))
        ! decide about the next level
        call coarsening_icl(mg(l)%cr_x,mg(l)%cr_y,mg(l)%icl_z, &
                            mg(l)%dx,mg(l)%dy,mg(l)%dz,A,minaspect,maxaspect)
        if (mg(l+1)%nn >= mg(l)%nn) then
            levels = l 
            exit ! the loop
        endif
        call crash('allocate and compute coarse X Y Z dz here')
        call coarsening_grid   ! set mg(l+1)%X,Y,Z,dx,dy,dz,nx,ny,nz,nn
    enddo
    print *,nlevels,' levels'

    ! assemble the matrices 
    do l=1,nlevels 
        call get_mg_dims(mg(l), &
            ifds, ifde, kfds,kfde, jfds, jfde,            & ! fire grid dimensions
            ifms, ifme, kfms,kfme, jfms, jfme,            &
            ifps, ifpe, kfps,kfpe, jfps, jfpe,           & ! fire patch bounds
            ifts, ifte, kfts,kfte, jfts,jfte)            
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

!*** executable

! F = div(u0)
! F = f_assembly_fortran(A,X,U0,lambda,params);
call f_assembly(                        &
    ifds, ifde, kfds, kfde, jfds, jfde,           & ! fire grid dimensions
    ifms, ifme, kfms, kfme, jfms, jfme,           &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte,            &
    A, mg(1)%X, mg(1)%Y, mg(1)%Z, u0, v0, w0,                       & !Input from femwind, U0, V0, W0 not used in hexa to construct Kloc or JG
    f)                                             !Global load vector output  

! initialize for now
u=0.
v=0.
w=0.
rate =0.

end subroutine femwind_solve
end module  module_femwind