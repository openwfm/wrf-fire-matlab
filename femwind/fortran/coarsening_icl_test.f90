program coarsening_icl_test

use module_coarsening   
use module_utils ! to read and write matrices as text files from matlab

implicit none

real, pointer, dimension(:,:,:):: icl_z_m,dx_m,dy_m,dz_m,A_m,minaspect_m,maxaspect_m ! to read from files
real, dimension(1,1,1)::cr_x_m,cr_y_m
real:: A(3,3), dx,dy,minaspect,maxaspect ! to pass on
real, pointer::dz(:)
integer, pointer:: icl_z(:)
integer:: cr_x,cr_y

integer::n3

! read and copy
call read_array(dx_m,'dx')
dx=dx_m(1,1,1)
call read_array(dy_m,'dy')
dy=dy_m(1,1,1)
call read_array(dz_m,'dz')
n3 = size(dz_m)
allocate(dz(n3))
dz=reshape(dz_m,(/n3/))
call read_array(A_m,'A')
A=reshape(A_m,(/3,3/))
call read_array(minaspect_m,'minaspect')
minaspect=minaspect_m(1,1,1)
call read_array(maxaspect_m,'maxaspect')
maxaspect=maxaspect_m(1,1,1)

call coarsening_icl(cr_x,cr_y,icl_z,dx,dy,dz,A,minaspect,maxaspect)

cr_x_m(1,1,1)=cr_x
call write_array(cr_x_m,'cr_x')
cr_y_m(1,1,1)=cr_y
call write_array(cr_y_m,'cr_y')
call write_array(reshape(real(icl_z),(/size(icl_z),1,1/)),'icl_z')

end program coarsening_icl_test
