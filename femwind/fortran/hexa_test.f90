program hexa_test

use module_hexa   ! testing only
use module_utils ! to read and write matrices as text files from matlab
use module_lin_alg

implicit none

real, pointer:: A(:,:,:), X(:,:,:), u0(:,:,:), aflag(:,:,:)    ! fortran is not case sensitive
integer:: iflag
!real:: A(3,3), X(3,8), u0(3)    ! fortran is not case sensitive
real, pointer:: Kloc(:,:,:), Floc(:,:,:), Jg(:,:,:) ! for convenience, make all arrays 3D
real:: vol

print *,'hexa_test.exe start:'
allocate(Kloc(8,8,1))
allocate(Floc(8,1,1))
allocate(Jg(8,3,1)) ! for convenience, make all arrays 3D

call read_array(A,'A')
call read_array(X,'X')
call read_array(u0,'u0')
call read_array(aflag,'iflags')


iflag = aflag(1,1,1)

call print_matrix('A',A(:,:,1))
call print_matrix('X',X(:,:,1))
call print_matrix('u0',u0(:,:,1))
write(*,*)'iflag=',iflag

write(*,*)'calling hexa'
call hexa(A(:,:,1),X(:,:,1),u0(:,:,1),Kloc(:,:,1),Floc(:,:,1),Jg(:,:,1),vol,iflag)
write(*,*)'hexa output'

call print_matrix('Kloc',Kloc(:,:,1))
call print_matrix('Floc',Floc(:,:,1))
call print_matrix('Jg',Jg(:,:,1))
write(*,*)'vol=',vol

call write_array(Kloc,'Kloc')
call write_array(Floc,'Floc')
call write_array(Jg,'Jg')

print *,'hexa_test.exe end:'

end program hexa_test
