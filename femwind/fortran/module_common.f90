module module_common 

use module_utils

! declarations of multigrid structure
integer, parameter::msize=14



type mg_type
    real, dimension(3,3):: A                        ! penalty weight matrix
    real, pointer, dimension(:,:,:):: X, Y, Z       ! grid vertices
    real, pointer, dimension(:,:,:)::        &      !
         F, lambda, res                             ! multigrid variables                                       
    
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
