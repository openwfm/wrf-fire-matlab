program hexa_test

use module_hexa   ! testing only
use module_utils ! to read and write matrices as text files from matlab

implicit none

real, pointer:: A(:,:,:), X(:,:,:), u0(:,:,:), aflag(:,:,:)    ! fortran is not case sensitive
integer:: iflag
!real:: A(3,3), X(3,8), u0(3)    ! fortran is not case sensitive
real, pointer:: Kloc(:,:,:), Floc(:,:,:), Jg(:,:,:) ! for convenience, make all arrays 3D
real:: vol
allocate(Kloc(8,8,1))
allocate(Floc(8,1,1))
allocate(Jg(8,3,1)) ! for convenience, make all arrays 3D

call read_array(A,'A')
call read_array(X,'X')
call read_array(u0,'u0')
call read_array(aflag,'iflags')


iflag = aflag(1,1,1)

write(*,*)'A=',A
write(*,*)'X=',X
write(*,*)'u0=',u0
write(*,*)'iflag=',iflag

write(*,*)'calling hexa'
call hexa(A,X,u0,Kloc,Floc,Jg,vol,iflag)
write(*,*)'hexa output'

write(*,*)'Kloc=',Kloc
write(*,*)'Floc=',Floc
write(*,*)'Jg=',Jg
write(*,*)'vol=',vol

call write_array(Kloc,'Kloc')
call write_array(Floc,'Floc')
call write_array(Jg,'Jg')

end program hexa_test
