! **** HOW TO TEST ****
! run femwind_fortran_rate_test
!
! details:
! called from femwind_test.f90 compiled as femwind_test.exe
! femwind_test.exe is executed from femwind_fortran.m if params.run_fortran
! femwind_fortran.m is called from femwind_main.m   
! run femwind_main from e.g. femwind_rate_test or femwind_wrfout_test
! femwind_fortran_rate_test sets the params and calls femwind_rate_test


module module_femwind
   use module_utils
   use module_f_assembly
   use module_ndt_assembly
   use module_w_assembly
   use module_coarsening
   use module_boundary_conditions
   use module_common

! parameters
integer, parameter::max_levels=20
real:: minaspect=1./3.,maxaspect=3.
real::A(3,3)=reshape((/1., 0., 0.,  0., 1., 0.,  0., 0., 1./),(/3, 3/))

integer::nlevels                                ! number of levels
       
contains

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

        ! decide on coarsening

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

    ! allocate grid variables
    do l=1,nlevels
        call allocate_mg_level(mg(l))
    enddo

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
            ifds, ifde, kfds, kfde, jfds, jfde,         & ! fire grid dimensions
            ifms, ifme, kfms, kfme, jfms, jfme,         &
            ifps, ifpe, kfps, kfpe, jfps, jfpe,         & ! fire patch bounds
            ifts, ifte, kfts, kfte, jfts, jfte)            
        allocate(mg(l)%Kglo(ifms:ifme, kfms:kfme, jfms:jfme, msize))
        call ndt_assembly(                              &
            ifds, ifde, kfds, kfde, jfds, jfde,         & ! fire grid dimensions
            ifms, ifme, kfms, kfme, jfms, jfme,         &
            ifps, ifpe, kfps, kfpe, jfps, jfpe,         & ! fire patch bounds
            ifts, ifte, kfts, kfte, jfts, jfte,         &      
            A, mg(l)%X,mg(l)%Y,mg(l)%Z, 1,              & 
            mg(l)%Kglo)        
        call ndt_boundary_conditions(                   &
            ifds, ifde, kfds, kfde, jfds, jfde,         & ! fire grid dimensions
            ifms, ifme, kfms, kfme, jfms, jfme,         &
            ifps, ifpe, kfps, kfpe, jfps, jfpe,         & ! fire patch bounds
            ifts, ifte, kfts, kfte, jfts, jfte,         &      
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

integer, intent(inout)::                             &
    ifds, ifde, kfds,kfde, jfds, jfde,            & ! fire grid dimensions
    ifms, ifme, kfms,kfme, jfms, jfme,            &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte            

real, intent(in), dimension(ifms:ifme, kfms:kfme, jfms:jfme):: u0, v0, w0   !initial wind vector at midpoints
real,intent(out), dimension(ifms:ifme, kfms:kfme, jfms:jfme)::u, v, w       !mass consistent wind at midpoints
real,intent(out):: rate

!*** local

!*** executable

! F = div(u0)
! F = f_assembly_fortran(A,X,U0,lambda,params);

call f_assembly(                        &
    ifds, ifde, kfds, kfde, jfds, jfde,           & ! fire grid dimensions
    ifms, ifme, kfms, kfme, jfms, jfme,           &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts, jfte,           &
    A, mg(1)%X, mg(1)%Y, mg(1)%Z, u0, v0, w0,     &                    	
    mg(1)%F)                                        !Global load vector output  

call vec_boundary_conditions(                              &
    ifds, ifde, kfds, kfde, jfds, jfde,           & ! fire grid dimensions
    ifms, ifme, kfms, kfme, jfms, jfme,           &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts, jfte,           &
    mg(1)%F)

! initialize for now
u=0.
v=0.
w=0.
rate =0.

end subroutine femwind_solve
end module  module_femwind
