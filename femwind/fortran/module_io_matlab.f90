module module_io_matlab

contains

subroutine crash(msg)
! terminate with error
character(len=*),intent(in)::msg
write(*,*)'FATAL:',trim(msg)
stop
end subroutine crash

subroutine read_array(n1,n2,n3,a,name)

!*** purpose read array from text file
! the dimensions must be known and agree with the file
implicit none

!*** arguments
integer, intent(in):: n1,n2,n3   ! dimensions
real, intent(out):: a(n1,n2,n3)  ! the array
character(len=*),intent(in)::name! file name root, .txt will be added

!*** internal
integer:: i1,i2,i3,iu=10,j
real*8:: x

!*** executable
open(unit=iunit,file=trim(name)//'.txt.',form='formatted',status='old',action='read')
read(iu,*)j
if(j.ne.456)call crash('read_array: wrong magic number')
read(iu,*)j
if(j.ne.3)call crash('read_array: must have 3 dimensions')

1 format('reading matrix size ',3i5,' from file ',a)
write(*,1)n1,n2,n3,trim(name)//'.txt.'

read(iu,*)i1
if(i1.ne.n1)call crash('read_array: n1 does not match')
read(iu,*)i2
if(i2.ne.n2)call crash('read_array: n2 does not match')
read(iu,*)i1
if(i3.ne.n3)call crash('read_array: n3 does not match')
do i3=1,n3
   do i2=1,n2
       do i1=1,n1   ! first index changing fastest
           read(iu,*) a(i1,i2,i3)
       enddo
   enddo
enddo
close(iu)
end subroutine read_array

subroutine write_array(n1,n2,n3,a,name)
!*** purpose write array to text file
integer, intent(in):: n1,n2,n3
real, intent(in):: a(n1,n2,n3)
character(len=*),intent(in)::name




end subroutine write_array

end module module_io_matlab
