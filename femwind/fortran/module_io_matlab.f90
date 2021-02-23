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
write(*,1)n,trim(name)//'.txt.'
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

end module module_io_matlab
