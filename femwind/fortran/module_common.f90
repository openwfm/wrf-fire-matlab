module module_common 

use module_utils

integer, parameter::msize=14
integer, parameter::max_levels=5

! method parameters
type params_type
    real:: minaspect=1./3.,maxaspect=3.
    real:: A(3,3)=reshape((/1., 0., 0.,  0., 1., 0.,  0., 0., 1./),(/3, 3/))
    integer:: coarsest_iter=100 ! to solve the coarse problem 
    integer:: maxit=50     ! total iterations
    integer:: nsmooth=3    ! smoothing iterations before correction
    integer:: nsmooth_coarse=2 ! smoothing iterations on coarse levels
    integer:: maxit_coarse=8 ! on levels>1: 2 smoothing, coarse, 2 smoothing, coarse, 2 smoothing
    integer:: debug_level=-1 ! multigrid level to debug up to
    integer:: msglevel=1 ! type of messages to print  
    integer:: print_level=2 ! multigrid level to print messages up to
    !logical:: check_relres=.true.  ! check residual after every base iteration if to continue
    logical:: check_relres=.false.
    real:: restol=1e-5 ! residual tolerance to stop iterating
    ! real:: restol=1e-6 ! residual tolerance to stop iterating
    !logical:: check_reldif=.false.  ! check difference after every base iteration if to continue
    logical:: check_reldif=.true.  ! check difference after every base iteration if to continue
    ! real:: diftol_finest=1e-6 ! keep iterating while iterates change by this
    real:: diftol_finest=1e-5 ! keep iterating while iterates change by this
    real:: diftol_coarse=1e-1 ! keep iterating while iterates change by this on corse levels
    real:: diftol_coarsest=1e-2 ! keep iterating while iterates change by this on corsest level
end type

type(params_type)::params

! multigrid structure

type mg_type
    real, dimension(3,3):: A                        ! penalty weight matrix
    real, pointer, dimension(:,:,:):: X, Y, Z       ! grid vertices
    real, pointer, dimension(:,:,:)::        &      !
         f, lambda, res                             ! multigrid variables                                       
    
    real, pointer, dimension(:,:,:,:):: Kglo        ! global stiffness matrix
    
    real:: dx,dy                                    ! horizontal spacing, scalar 
    real, pointer, dimension(:)::dz                 ! vertical spacing of the layers
    integer::nx,ny,nz,nn                            ! mesh size in terms of vertices
  
    integer::level

    integer:: cr_x, cr_y                            ! coarsening factors
    integer, pointer, dimension(:):: icl_x, icl_y, icl_z  ! coarsening indices
 
end type

type(mg_type):: mg(max_levels+1)  ! the main multigrid structure

contains

subroutine read_params
    real:: minaspect=1./3.,maxaspect=3.
    real:: A(3,3)=reshape((/1., 0., 0.,  0., 1., 0.,  0., 0., 1./),(/3, 3/))
    integer:: coarsest_iter=100 ! to solve the coarse problem 
    integer:: maxit=50     ! total iterations
    integer:: nsmooth=3    ! smoothing iterations before correction
    integer:: nsmooth_coarse=2 ! smoothing iterations on coarse levels
    integer:: maxit_coarse=8 ! on levels>1: 2 smoothing, coarse, 2 smoothing, coarse, 2 smoothing
    integer:: debug_level=-1 ! multigrid level to debug up to
    integer:: msglevel=1 ! type of messages to print  
    integer:: print_level=2 ! multigrid level to print messages up to
    !logical:: check_relres=.true.  ! check residual after every base iteration if to continue
    logical:: check_relres=.false.
    real:: restol=1e-6 ! residual tolerance to stop iterating
    !logical:: check_reldif=.false.  ! check difference after every base iteration if to continue
    logical:: check_reldif=.true.  ! check difference after every base iteration if to continue
    real:: diftol_finest=1e-6 ! keep iterating while iterates change by this
    real:: diftol_coarse=1e-1 ! keep iterating while iterates change by this on corse levels
    real:: diftol_coarsest=1e-2 ! keep iterating while iterates change by this on corsest level

    namelist/params_nml/ & 
      minaspect, &
      A, &
      coarsest_iter, &
      maxit, &
      nsmooth, &
      nsmooth_coarse, &
      maxit_coarse, &
      debug_level, &
      msglevel, &
      print_level, &
      check_relres, &
      restol, &
      check_reldif, &
      diftol_finest, &
      diftol_coarse, &
      diftol_coarsest

    write(*,a)'Default parameters:'
    write(*,nml=params_nml)
    open(8,file='params.nml',status='old',err=9)
    read(8,nml=params_nml)
    close(8)
     params%minaspect=minaspect
     params%A=A 
     params%coarsest_iter=coarsest_iter
     params%maxit=maxit
     params%nsmooth=nsmooth
     params%nsmooth_coarse=nsmooth_coarse
     params%maxit_coarse=maxit_coarse
     params%debug_level=debug_level
     params%msglevel=msglevel
     params%print_level=print_level
     params%check_relres=check_relres
     params%restol=restol
     params%check_reldif=check_reldif
     params%diftol_finest=diftol_finest
     params%diftol_coarse=diftol_coarse
     params%diftol_coarsest=diftol_coarsest
    write(*,a)'Using modified parameters:'
    write(*,nml=params_nml)
    return
9   write(*,a)'Cannot open file params.nml, using defaults'
    end subroutine read_params
    

subroutine allocate_mg_level(mg)    !allocate one multigrid level
type(mg_type),intent(inout)::mg     !multigrid structure
!*** local
integer:: &
    ifds, ifde, kfds, kfde, jfds, jfde,           & ! fire grid dimensions
    ifms, ifme, kfms, kfme, jfms, jfme,           & ! memory dimensions
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts, jfte              ! tile dimensions

!*** executable
call get_mg_dims(mg, &
    ifds, ifde, kfds,kfde, jfds, jfde,            & ! fire grid dimensions
    ifms, ifme, kfms,kfme, jfms, jfme,            &
    ifps, ifpe, kfps,kfpe, jfps, jfpe,            & ! fire patch bounds
    ifts, ifte, kfts,kfte, jfts,jfte)      

if (.not.associated(mg%X))allocate(mg%X(ifms: ifme, kfms: kfme, jfms: jfme))
if (.not.associated(mg%Y))allocate(mg%Y(ifms: ifme, kfms: kfme, jfms: jfme))
if (.not.associated(mg%Z))allocate(mg%Z(ifms: ifme, kfms: kfme, jfms: jfme))
allocate(mg%F(ifms: ifme, kfms: kfme, jfms: jfme))
allocate(mg%lambda(ifms: ifme, kfms: kfme, jfms: jfme))
allocate(mg%res(ifms: ifme, kfms: kfme, jfms: jfme))

end subroutine allocate_mg_level

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
character(len=256)::msg
 
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
write(msg,*)"get_mg_dims: domain",ifds, ifde, kfds, kfde, jfds, jfde
call message(msg)
write(msg,*)"get_mg_dims: memory",ifms, ifme, kfms, kfme, jfms, jfme
call message(msg)
write(msg,*)"get_mg_dims: patch ",ifps, ifde, kfps, kfpe, jfps, jfpe
call message(msg)
write(msg,*)"get_mg_dims: tile  ",ifts, ifte, kfts, kfte, jfts, jfte
call message(msg)

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
    print *,'grid size nx=',mg(l)%nx,' ny=',mg(l)%ny,' nz=',mg(l)%nz, &
       ' total ndof=',mg(l)%nn
    print *,'dx=',mg(l)%dx,' dy=',mg(l)%dy
    if(associated(mg(l)%dz))then
        print *,'dz shape',shape(mg(l)%dz)
        if(product(shape(mg(l)%dz))>0) then
            print *,'size ',size(mg(l)%dz)
            print *,'dz=',mg(l)%dz
        endif   
    endif
end subroutine print_mg_dims

end module module_common 
