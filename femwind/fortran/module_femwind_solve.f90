module module_femwind_solve
       use module_f_assembly
       
contains

subroutine femwind_setup(                         &
    ifds, ifde, kfds, kfde, jfds, jfde,           & ! fire grid dimensions
    ifms, ifme, kfms, kfme, jfms, jfme,           &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte,            &
    A, X, Y, Z,                                   & !  inputs
    K)                                              ! outputs 

integer, intent(in)::                             &
    ifds, ifde, kfds,kfde, jfds, jfde,            & ! fire grid dimensions
    ifms, ifme, kfms,kfme, jfms, jfme,            &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte            

integer, parameter::msize=14
real, intent(in), dimension(3,3):: A
real, intent(in), dimension(ifms:ifme, kfms:kfme, jfms:jfme):: X,Y,Z        !spatial grid at corners
real,intent(out), dimension(ifms:ifme, kfms:kfme, jfms:jfme, msize):: K     !stiffness matrix


end subroutine femwind_setup

subroutine femwind_solve(                         &
    ifds, ifde, kfds, kfde, jfds, jfde,           & ! fire grid dimensions
    ifms, ifme, kfms, kfme, jfms, jfme,           &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte,            &
    A, X, Y, Z, u0, v0, w0,                       & !  inputs
    u, v, w, rate)                                  ! outputs 

implicit none

!*** arguments

integer, intent(in)::                             &
    ifds, ifde, kfds,kfde, jfds, jfde,            & ! fire grid dimensions
    ifms, ifme, kfms,kfme, jfms, jfme,            &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte            

real, intent(in), dimension(3,3):: A
real, intent(in), dimension(ifms:ifme, kfms:kfme, jfms:jfme):: X,Y,Z, &     !spatial grid at corners
                                                               u0, v0, w0   !initial wind vector at midpoints
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
    A, X, Y, Z, u0, v0, w0,                       & !Input from femwind, U0, V0, W0 not used in hexa to construct Kloc or JG
    f)                                             !Global load vector output  

! initialize for now
u=0.
v=0.
w=0.
rate =0.
end subroutine femwind_solve
end module  module_femwind_solve
