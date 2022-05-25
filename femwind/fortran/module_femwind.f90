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
   use module_multigrid
   use module_common
   use module_utils
   use module_f_assembly
   use module_ndt_assembly
   use module_w_assembly
   ! use module_ndt_mult
   ! use module_coarsening
   use module_boundary_conditions
   ! use module_sweeps

contains

subroutine femwind_setup(mg)
implicit none
! set up the mg_level structure
! input: mg(1)%X,Y,Z,dx,dy,dz already set
!*** arguments
type(mg_type),intent(inout)::mg(:)  ! multigrid level
!*** local

!*** executable

call multigrid_setup(mg)

end subroutine femwind_setup

subroutine femwind_solve(mg,      &
     ifds, ifde, kfds,kfde, jfds, jfde,            & ! fire grid dimensions
     ifms, ifme, kfms,kfme, jfms, jfme,            &
     ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
     ifts, ifte, kfts, kfte, jfts,jfte,            &
     u0, v0, w0,                                   & ! input 
     u, v, w, rate)                                  ! output 

implicit none

!*** arguments

type(mg_type),intent(inout)::mg(:)  ! multigrid levels

integer, intent(in)::                             &
    ifds, ifde, kfds,kfde, jfds, jfde,            & ! fire grid dimensions
    ifms, ifme, kfms,kfme, jfms, jfme,            &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts,jfte            

real, intent(in), dimension(ifms:ifme, kfms:kfme, jfms:jfme):: u0, v0, w0   !initial wind vector at midpoints
real,intent(out), dimension(ifms:ifme, kfms:kfme, jfms:jfme):: u, v, w      !mass consistent wind at midpoints
real,intent(out):: rate

!*** local

!*** executable

! f = div(u0)
! f = f_assembly_fortran(A,X,U0,lambda,params);

print *,'femwind solve start'
print *,'calling f_assembly'
call f_assembly(                                  &
    ifds, ifde, kfds, kfde, jfds, jfde,           & ! fire grid dimensions
    ifms, ifme, kfms, kfme, jfms, jfme,           &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts, jfte,           &
    params%A, mg(1)%X, mg(1)%Y, mg(1)%Z,          &
    u0, v0, w0,                                   &                    	
    mg(1)%f)                                        !Global load vector output  

print *,'calling vec_boundary_conditions'
call vec_boundary_conditions(                              &
    ifds, ifde, kfds, kfde, jfds, jfde,           & ! fire grid dimensions
    ifms, ifme, kfms, kfme, jfms, jfme,           &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts, jfte,           &
    mg(1)%f)

print *,'calling multigrid_cycle'
call multigrid_cycle(mg,1,rate)                        ! start multigrid from level 1

if(params%debug_level >=0)call  write_array(mg(1)%lambda(ifts: ifte+1, kfts: kfte+1, jfts:jfte+1),'lambda_sol')

if(params%debug_level >=0)call  write_array(u0(ifts: ifte, kfts: kfte, jfts:jfte),'u0_sol')
if(params%debug_level >=0)call  write_array(v0(ifts: ifte, kfts: kfte, jfts:jfte),'v0_sol')
if(params%debug_level >=0)call  write_array(w0(ifts: ifte, kfts: kfte, jfts:jfte),'w0_sol')
if(params%debug_level >=0)call  write_array(mg(1)%X(ifts: ifte+1, kfts: kfte+1, jfts:jfte+1),'X_sol')
if(params%debug_level >=0)call  write_array(mg(1)%Y(ifts: ifte+1, kfts: kfte+1, jfts:jfte+1),'Y_sol')
if(params%debug_level >=0)call  write_array(mg(1)%Z(ifts: ifte+1, kfts: kfte+1, jfts:jfte+1),'Z_sol')

print *,'calling w_assembly'
call w_assembly(                                  &
    ifds, ifde, kfds, kfde, jfds, jfde,           & ! fire grid dimensions
    ifms, ifme, kfms, kfme, jfms, jfme,           &
    ifps, ifpe, kfps, kfpe, jfps, jfpe,           & ! fire patch bounds
    ifts, ifte, kfts, kfte, jfts, jfte,           & 
    mg(1)%lambda, u0, v0, w0,                     &
    params%A, mg(1)%X, mg(1)%Y, mg(1)%Z,          & !Input from femwind, u0, v0, w0, Spatial Grid Data
    u, v, w)                                        ! final output  

print *,'end femwind_solve'

if(params%debug_level >=0)call  write_array(u(ifts: ifte, kfts: kfte, jfts:jfte),'u_sol')
if(params%debug_level >=0)call  write_array(v(ifts: ifte, kfts: kfte, jfts:jfte),'v_sol')
if(params%debug_level >=0)call  write_array(w(ifts: ifte, kfts: kfte, jfts:jfte),'w_sol')

end subroutine femwind_solve


end module  module_femwind
