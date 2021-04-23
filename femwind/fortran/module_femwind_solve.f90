module module_femwind_solve
       use module_f_assembly
       
contains

subroutine femwind_solve(                         &
    ifds, ifde, kfds, kfde, jfds, jfde,           & ! fire grid dimensions
    ifms, ifme, kfms, kfme, jfms, jfme,           &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte,            &
    A, X, Y, Z, u0, v0, w0,                       & !  inputs
    u, v, w)                                        ! outputs 

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

!*** local
real, dimension(ifms:ifme, kfms:kfme, jfms:jfme):: f  

!*** executable

! F = div(u0)
call f_assembly(                        &
    ifds, ifde, kfds, kfde, jfds, jfde,           & ! fire grid dimensions
    ifms, ifme, kfms, kfme, jfms, jfme,           &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte,            &
    A, X, Y, Z, u0, v0, w0,                       & !Input from femwind, U0, V0, W0 not used in hexa to construct Kloc or JG
    f)                                             !Global load vector output  

end subroutine femwind_solve
end module  module_femwind_solve
