module module_io_matlab

contains

subroutine crash(msg)
! terminate with error
character(len=*),intent(in)::msg
write(*,*)'FATAL:',trim(msg)
stop
end subroutine crash

subroutine read_array(a,name)

!*** purpose read array from text file
implicit none

!*** arguments
real, pointer, intent(out):: a(:,:,:)  ! the array pointer; remember to deallocate when done with it
character(len=*),intent(in)::name! file name root, .txt will be added

!*** internal
integer:: n(3),iu=8,j

!*** executable
open(iu,file=trim(name)//'.txt',form='formatted',status='old')
read(iu,*)j
if(j.ne.456)call crash('read_array: wrong magic number')
read(iu,*)j
if(j.ne.3)call crash('read_array: must have 3 dimensions')
read(iu,*)n
1 format('reading matrix size ',3i5,' from file ',a)
write(*,1)n,trim(name)//'.txt'
allocate(a(n(1),n(2),n(3)))
read(iu,*)a
close(iu)
end subroutine read_array

subroutine write_array(a,name)
!*** purpose write array to text file
!*** arguments
real, pointer, intent(in):: a(:,:,:)  
character(len=*),intent(in)::name! file name root, .txt will be added

!*** internal
integer:: n(3),iu=8,i1,i2,i3

!*** executable
iu=8
open(iu,file=trim(name)//'.txt',form='formatted',status='unknown')
write(iu,*)456
write(iu,*)3
n=shape(a)
write(iu,*)n(1)
write(iu,*)n(2)
write(iu,*)n(3)
1 format('writing matrix size ',3i5,' to file ',a)
write(*,1)n,trim(name)//'.txt.'
do i3=1,n(3)
do i2=1,n(2)
do i1=1,n(1)
write(iu,*)a(i1,i2,i3)
enddo
enddo
enddo
close(iu)
end subroutine write_array

subroutine read_array_nd(a,s,name)
!*** purpose read nd array from text file
!   usage
!   integer :: s(k)  ! k is constant
!   real, pointer :: a(:)
!   call read_array_nd(a,s,'file')
!   target_array = reshape(a,s)
!   
implicit none

!*** arguments
integer, intent(out):: s(:)  ! the array shape pointer
real, pointer, intent(out):: a(:)  ! the array data pointer
character(len=*),intent(in)::name! file name root, .txt will be added

!*** internal
integer:: iu=8,n,j,sn(1)

!*** executable
open(iu,file=trim(name)//'.txt',form='formatted',status='old')
read(iu,*)j
if(j.ne.456)call crash('read_array: wrong magic number')
read(iu,*)j
sn = shape(s)
if(j.ne.sn(1))call crash('read_array: wrong number of dimensions')
if(j.lt.1.or.j.gt.7)call crash('read_array: must have 1 to 7 dimensions')

read(iu,*)s
1 format('file ',a,' reading matrix size ',7i8)
write(*,1)trim(name)//'.txt.',s
allocate(a(product(s)))
read(iu,*)a
close(iu)
end subroutine read_array_nd

subroutine write_array_nd(a,s,name)
!*** purpose write array to text file
!*** usage
!   s = shape(a)
!   call write_array_nd(reshape(a,(/product(s)\)),'file')

implicit none
!*** arguments
integer, intent(in):: s(:)  
real, intent(in):: a(:)  
character(len=*),intent(in)::name! file name root, .txt will be added

!*** internal
integer::n(1),m(1),i,iu
!*** executable
iu=8
1 format('writing ',a,' matrix size ',7i5)
write(*,1)trim(name)//'.txt',s
n = shape(s)
m = shape(a)
if (product(s).ne.m(1))call crash('write_array_nd: wrong size of a')
open(iu,file=trim(name)//'.txt',form='formatted',status='unknown')
write(iu,*)456
write(iu,*)n
do i=1,n(1)
    write(iu,*)s(i)
enddo
do i=1,m(1)
    write(iu,*)a(i)
enddo
close(iu)
end subroutine write_array_nd

end module module_io_matlab
