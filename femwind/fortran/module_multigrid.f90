module module_multigrid
   use module_common
   use module_utils
   use module_f_assembly
   use module_ndt_assembly
   use module_w_assembly
   use module_ndt_mult
   use module_coarsening
   use module_boundary_conditions
   use module_sweeps

integer::nlevels                                ! number of levels
       
contains

subroutine multigrid_setup(mg)
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

    print *,'multigrid_setup received'
    call print_mg_dims(mg,1)
    
    do l=1,max_levels

        ! decide on coarsening

        print *,'multigrid level ',l,' grid size ', mg(l)%nx, mg(l)%ny, mg(l)%nz
        ! decide about the next level
        call coarsening_icl(mg(l)%cr_x,mg(l)%cr_y,mg(l)%icl_z,  &
                            mg(l)%dx,mg(l)%dy,mg(l)%dz,         &
                            params%A,params%minaspect,params%maxaspect)

        ! get horizontal coarsening lists

        call coarsening_hzc2icl(mg(l)%icl_x, mg(l)%icl_y, &
                                mg(l)%cr_x,mg(l)%cr_y,&
                                mg(l)%nx, mg(l)%ny)

        ! update coarse mesh scalars

        mg(l+1)%nx = size(mg(l)%icl_x)
        mg(l+1)%ny = size(mg(l)%icl_y)
        mg(l+1)%nz = size(mg(l)%icl_z)
        mg(l+1)%nn = mg(l+1)%nx *  mg(l+1)%ny * mg(l+1)%nz
        mg(l+1)%dx = mg(l)%dx * mg(l)%cr_x
        mg(l+1)%dy = mg(l)%dy * mg(l)%cr_y
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
        if(params%debug_level >=l)call write_3d(mg,l,mg(l)%X,'X','v')
        if(params%debug_level >=l)call write_3d(mg,l,mg(l)%Y,'Y','v')
        if(params%debug_level >=l)call write_3d(mg,l,mg(l)%Z,'Z','v')
    enddo


    ! assemble the stiffness matrices 
    do l=1,nlevels 
        call get_mg_dims(mg(l),                         &
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
            params%A, mg(l)%X,mg(l)%Y,mg(l)%Z, 1,              & 
            mg(l)%Kglo)        
        call ndt_boundary_conditions(                   &
            ifds, ifde, kfds, kfde, jfds, jfde,         & ! fire grid dimensions
            ifms, ifme, kfms, kfme, jfms, jfme,         &
            ifps, ifpe, kfps, kfpe, jfps, jfpe,         & ! fire patch bounds
            ifts, ifte, kfts, kfte, jfts, jfte,         &      
            mg(l)%Kglo)        
    enddo

end subroutine multigrid_setup


recursive subroutine multigrid_cycle(mg,l,rate)  
implicit none
 
!*** arguments

type(mg_type),intent(inout)::mg(:)  ! multigrid levels
integer, intent(in)::l           ! level
real, intent(out)::rate

!*** local

integer::   &
    ifds, ifde, kfds, kfde, jfds, jfde,           & ! fire grid dimensions
    ifms, ifme, kfms, kfme, jfms, jfme,           & ! memory dimensions
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts, jfte,           & ! tile dimensions
    ifcds, ifcde, kfcds,kfcde, jfcds,jfcde,       & ! coarse grid domain
    ifcms, ifcme, kfcms,kfcme, jfcms,jfcme,       & ! coarse grid dimensions
    ifcps, ifcpe, kfcps,kfcpe, jfcps,jfcpe,       & ! coarse grid dimensions
    ifcts, ifcte, kfcts,kfcte, jfcts,jfcte          ! coarse grid tile

integer::it, maxit, nit, nsmooth, cycles
logical::coarse,coarsest
real:: tol, norm_f, norm_res
character(len=10)::it_kind

!*** executable

call get_mg_dims(mg(l),                         &
    ifds, ifde, kfds, kfde, jfds, jfde,         & ! fire grid dimensions
    ifms, ifme, kfms, kfme, jfms, jfme,         &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,         & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts, jfte)
if(l<nlevels) then 
    call get_mg_dims(mg(l+1),                         &
    ifcds, ifcde, kfcds,kfcde, jfcds,jfcde,       & ! coarse grid domain
    ifcms, ifcme, kfcms,kfcme, jfcms,jfcme,       & ! coarse grid dimensions
    ifcps, ifcpe, kfcps,kfcpe, jfcps,jfcpe,       & ! coarse grid dimensions
    ifcts, ifcte, kfcts,kfcte, jfcts,jfcte)         ! coarse grid tile
!*** temporary fix, meshes are in vertices not elements here.
!*** should be fixed using snode in all subroutines called from here 
    ifcte = ifcte + 1
    jfcte = jfcte + 1
    kfcte = kfcte + 1
endif

!*** temporary fix, meshes are in vertices not elements here.
!*** should be fixed using snode in all subroutines called from here 
ifte = ifte + 1
jfte = jfte + 1
kfte = kfte + 1

! initialize
cycles = 0

! decide on number of iterations
coarsest = l.eq.nlevels
if(coarsest)then
    ! coarsest level
    maxit = params%coarsest_iter
    nsmooth = 0
elseif(l == 1)then
    ! top level
    maxit = params%maxit
    nsmooth = params%nsmooth
else  
    ! some coarse level in between
    maxit = params%maxit_coarse
    nsmooth = params%nsmooth_coarse
endif

mg(l)%lambda=0.                             ! initial solution 0

if(params%check_relres)then
    norm_f = norm2( mg(l)%f(ifts: ifte, kfts: kfte, jfts: jfte))  !
    tol = norm_f * params%restol
    print *,'starting mg cycle level ',l,' norm_f ',norm_f,' residual tolerance ',tol
endif


do it=1,maxit
    coarse = mod(it,nsmooth+1)==0 .and. l < nlevels
    print *,'level=',l,' it=',it,' coarse=',coarse
    if(coarse)then                                      ! coarse correction
        if(params%debug_level >=l)call  write_array(mg(l)%lambda(ifts: ifte, kfts: kfte, jfts:jfte),'cc_lambda_in')
        it_kind='coarse correction'
        ! compute residual residual = f - Kglo*lambda
        call ndt_mult(  &
            ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
            ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
            ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire patch bounds
            ifts, ifte, kfts, kfte, jfts,jfte,                        & ! fire tile bounds
            mg(l)%Kglo, mg(l)%lambda, mg(l)%f, mg(l)%res)
        if(params%debug_level >=l)call  write_array(mg(l)%res(ifts: ifte, kfts: kfte, jfts:jfte),'cc_res')
        ! restriction: f_coarse = R*residual      
        call restriction(   &
            ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
            ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
            ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire patch bounds
            ifts, ifte, kfts, kfte, jfts,jfte,                        & ! fire tile boundss                
            ifcds, ifcde, kfcds,kfcde, jfcds,jfcde,       & ! coarse grid domain
            ifcms, ifcme, kfcms,kfcme, jfcms,jfcme,       & ! coarse grid dimensions
            ifcps, ifcpe, kfcps,kfcpe, jfcps,jfcpe,       & ! coarse grid dimensions
            ifcts, ifcte, kfcts,kfcte, jfcts,jfcte,       & ! coarse grid tile                
            mg(l+1)%f,mg(l)%res,                          &
            mg(l)%cr_x,mg(l)%cr_y,mg(l)%icl_z,            &
            mg(l)%X,mg(l)%Y,mg(l)%Z)
        if(params%debug_level >=l)call  write_array(mg(l+1)%f(ifcts: ifcte, kfcts: kfcte, jfcts:jfcte),'cc_f_c')
        if(params%debug_level >=l)call pause
        !if (l.eq.2) call crash('stopping here')
        ! call self on level l+1
        call multigrid_cycle(mg,l+1,rate)

        ! prolongation lambda = lambda + P*lambda_coarse

        if(params%debug_level >=l)call  write_array(mg(l+1)%lambda(ifcts: ifcte+1, kfcts: kfcte+1, jfcts:jfcte+1),'cc_lambda_c')
        call prolongation(   &
            ifds, ifde, kfds, kfde, jfds, jfde,           & ! fire grid dimensions
            ifms, ifme, kfms, kfme, jfms, jfme,           & ! memory dimensions
            ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
            ifts, ifte, kfts, kfte, jfts, jfte,           & ! tile dimensions
            ifcds, ifcde, kfcds,kfcde, jfcds,jfcde,       & ! coarse grid domain
            ifcms, ifcme, kfcms,kfcme, jfcms,jfcme,       & ! coarse grid dimensions
            ifcps, ifcpe, kfcps,kfcpe, jfcps,jfcpe,       & ! coarse grid dimensions
            ifcts, ifcte, kfcts,kfcte, jfcts,jfcte,       & ! coarse grid tile
            mg(l)%lambda,mg(l+1)%lambda,                  &
            mg(l)%cr_x,mg(l)%cr_y,mg(l)%icl_z,            &
            mg(l)%X,mg(l)%Y,mg(l)%Z)
    else
        it_kind = 'smoothing'
        if(coarsest)it_kind='coarsest'
        if(params%debug_level >=l)call  write_array(mg(l)%f(ifts: ifte, kfts: kfte, jfts:jfte),'sweep_f_in')
        if(params%debug_level >=l)call  write_array(mg(l)%lambda(ifts: ifte, kfts: kfte, jfts:jfte),'sweep_lambda_in')
        call sweeps( &
            ifds, ifde, kfds, kfde, jfds, jfde,           & ! fire grid dimensions
            ifms, ifme, kfms, kfme, jfms, jfme,           & ! memory dimensions
            ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
            ifts, ifte, kfts, kfte, jfts, jfte,           & ! tile dimensions                  
            mg(l)%Kglo, mg(l)%f, mg(l)%lambda) 
        if(params%debug_level >=l)call  write_array(mg(l)%lambda(ifts: ifte, kfts: kfte, jfts:jfte),'sweep_lambda_out')
        if((.not.coarsest.or.it==maxit).and.(params%debug_level >=l))call pause
    endif

    if(params%check_relres)then
        call ndt_mult(  &
            ifds, ifde, kfds, kfde, jfds, jfde,                       & ! fire domain bounds
            ifms, ifme, kfms, kfme, jfms, jfme,                       & ! fire memory bounds
            ifps, ifpe, kfps, kfpe, jfps, jfpe,                       & ! fire patch bounds
            ifts, ifte, kfts, kfte, jfts,jfte,                        & ! fire tile bounds
            mg(l)%Kglo, mg(l)%lambda, mg(l)%f, mg(l)%res)
        if(params%debug_level >=l)call  write_array(mg(l+1)%res(ifcts: ifcte, kfcts: kfcte, jfcts:jfcte),'it_res')
        norm_res = norm2( mg(l)%res(ifts: ifte, kfts: kfte, jfts: jfte))   ! residual norm
        if (mod(it,params%nsmooth+1)==params%nsmooth)then
            cycles=cycles+1;
        endif
        print *,'cycle ',cycles,' level ',l, 'avg rate ',rate
        rate = (norm_res/norm_f)**(1d0/cycles);
        print *,'level ',l,' iteration ',it,' ',it_kind,' residual ',norm_res, &
            ' norm_f ',norm_f,' rate ',rate
        if(norm_res < tol)exit
    endif
    if(params%debug_level >=l)call  write_array(mg(l)%lambda(ifts: ifte, kfts: kfte, jfts:jfte),'lambda_out')

        
enddo

end subroutine multigrid_cycle

end module  module_multigrid
