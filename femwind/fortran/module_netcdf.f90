module module_netcdf

use netcdf
use module_utils

integer::netcdf_msglevel = 0

contains

! from https://github.com/openwfm/wrf-fire/blob/master/standalone/wrf_netcdf.F
subroutine ncopen(filename,mode,ncid)
!*** purpose: open netcdf file wrapper with an informative error message
implicit none
!*** arguments
character(len=*), intent(in):: filename
integer, intent(in)::mode
integer, intent(out):: ncid
!*** executable
call check(nf90_open(trim(filename),mode,ncid),"Cannot open file "//trim(filename))
print *,"Opened netcdf file ",trim(filename)," as ",ncid," mode ",mode
end subroutine ncopen

subroutine ncclose(ncid)
!*** purpose: open netcdf file wrapper with informative error message
implicit none
!*** arguments
integer, intent(in):: ncid
print *,"Closing netcdf file ",ncid
call check(nf90_close(ncid),"Cannot close netcdf file ")
end subroutine ncclose

real function netcdf_read_att(ncid,name)
! read real global scalar attribute
!*** arguments
integer, intent(in)::ncid
character(len=*), intent(in)::name
!*** local
real::value
integer::xtype,len,attnum
character(len=256)::msg
!*** executable
call check(nf90_inquire_attribute(ncid, nf90_global, trim(name), xtype, len, attnum),"nf90_inquire_attribute")
if(xtype.ne.nf90_float.or.len.ne.1)then
   write(msg,*)"netcdf_read_att can read only float scalar but ",trim(name)," has xtype=",xtype," len=",len
   call crash(trim(msg))
endif
call check(nf90_get_att(ncid,nf90_global, trim(name), value),"nf90_get_att")
write(msg,*)"netcdf_read_att returning ",trim(name),"=",value
call message(msg)
netcdf_read_att = value
end function netcdf_read_att

integer function netcdf_read_int_wrf(ncid,name,istep)
implicit none
!*** Read one integer 
!*** arguments
    integer, intent(in)::ncid                      ! open netcdf file
    character(LEN=*),intent(in)::name              ! variable name
    integer, intent(in)::istep                     ! index in unlimited dimension (timestep number)
!*** local
    integer::ia(1)                     ! variable to store
    integer::ierr,varid
!*** executable
        print *,"netcdf_read_int_wrf reading variable ",trim(name)," time step ",istep
        call check(nf90_inq_varid(ncid, trim(name), varid), &
            "netcdf_read_int_wrf/nf90_inq_varid:"//trim(name))
        call check(nf90_get_var(ncid, varid, ia, start = (/istep/), count = (/1/)), &
             "netcdf_read_int_wrf/nf90_get_var:"//trim(name))
        print *,"netcdf_read_int_wrf:", trim(name), " = ",ia
        netcdf_read_int_wrf = ia(1)
end function netcdf_read_int_wrf

subroutine netcdf_write_int(ncid,ia,varname)
implicit none
!*** arguments
    integer, intent(in)::                         &
    ncid,                                         & ! open netcdf file
    ia                                              ! variable to write
    character(LEN=*),intent(in):: varname
!*** local
    integer::varid,ival
    character(len=256)::msg
!*** executable
    write(msg,*)'netcdf_write_int: varname=',varname,' value=',ia
    call message(msg)
        call check(nf90_inq_varid(ncid, trim(varname), varid), &
            "netcdf_write_int/nf90_inq_varid:"//trim(varname))
        ival = ia
        call check(nf90_put_var(ncid, varid, ival), &
            "netcdf_write_int/nf90_put_var:"//trim(varname))
end subroutine netcdf_write_int

subroutine netcdf_write_array(ncid,a,name)
    implicit none

!*** arguments
    integer, intent(in)::ncid                    ! open netcdf file
    real,intent(in),dimension(:,:,:)::a  
    character(LEN=*),intent(in):: name

!*** local
    integer,dimension(4)::star,cnts
    integer::i,j,k,varid,ends(4),dims(4),n(3)
    real,dimension(:,:,:,:),allocatable::at
    character(len=256) msg

    ! get idx
    n=shape(a)
    call netcdf_var_info(ncid,name,dims,varid,netcdf_msglevel)
    star   = (/1,1,1,1/)
    ends   = (/dims(1),dims(2),dims(3),1/)
    ends   = min(ends,dims)
    cnts = ends - star + 1
    
    ! transpose a -> at
    allocate(at(star(1):ends(1),star(2):ends(2),star(3):ends(3),1))
    do k=star(3),ends(3)
        do j=star(2),ends(2)
            do i=star(1),ends(1)
                at(i,j,k,1)=a(i,k,j)
            enddo
        enddo
    enddo

    if(netcdf_msglevel>=0) &
      write(msg,*)"writing ",trim(name),n(1),star(1),ends(1),n(2),star(2),ends(2),n(3),star(3),ends(3)
    call message(msg)
 
    ! write to file
    call check(nf90_put_var(ncid, varid, at, start = star, count = cnts),"nf90_put_var:"//trim(name))

    deallocate(at)
    
end subroutine netcdf_write_array


subroutine netcdf_write_2d(ncid,a,name,iframe)
    implicit none
!*** purpose
!   write a 2d array to netcdf file

!*** arguments
    integer, intent(in)::ncid                    ! open netcdf file
    real,intent(in),dimension(:,:)::a  
    character(LEN=*),intent(in):: name
    integer, intent(in)::iframe                  ! time frame to write in 

!*** local
    integer,dimension(3)::star,cnts
    integer::i,j,k,varid,ends(3),dims(3),n(2)
    character(len=256) msg

    ! get idx
    n=shape(a)
    call netcdf_var_info(ncid,name,dims,varid,netcdf_msglevel)
    write(msg,*)"array ",trim(name)," shape ",n," NetCDF dimensions ",dims
    call message(msg)
    
    if(dims(1).lt.n(1).or.dims(2).lt.n(2))call crash("array shape too large")
    star   = (/1,1,iframe/)
    ends   = (/n(1),n(2),iframe/)
    if(iframe.gt.dims(3))call crash('netcdf_write_2d: frame not in file')
    cnts = ends - star + 1
    
    write(msg,*)"writing ",trim(name)," from ",star," to ",ends
    call message(msg)
 
    ! write to file
    call check(nf90_put_var(ncid, varid, a, start = star, count = cnts),"nf90_put_var:"//trim(name))

end subroutine netcdf_write_2d

integer function l2i(l)
    implicit none
    logical, intent(in)::l
    if(l)then
        l2i = 1
    else
        l2i = 0
    endif
end function l2i


subroutine netcdf_var_info(ncid,varname,dims,varid,prints)
    implicit none
!*** arguments
    integer, intent(in)::ncid
    character(len=*)::varname
    integer,intent(out)::dims(:),varid
    integer,intent(in),optional::prints 
!*** local
    integer, parameter::mdims = 256
    integer:: xtype, ndims, natts, dimids(mdims),i,j,attnum
    integer :: len
    character(len=nf90_max_name):: name
    integer:: values_int(mdims)
    real:: values_real(mdims)
    character(len=mdims):: values_char
    character(LEN=256):: filename, msg
    logical::verbose=.true.

    if(present(prints)) verbose = prints>0

    call check(nf90_inq_varid(ncid,trim(varname),varid),"nf90_inq_varid"//trim(varname))
    call check(nf90_inquire_variable(ncid, varid, name, xtype, ndims, dimids, natts),"nf90_inquire_variable")
    if(ndims>mdims)call crash("netcdf_var_info: increase mdims")
    if(ndims>size(dims))call crash("netcdf_var_info: dims too short")
    if(verbose)then
        write(msg,*)"variable ",trim(name), " xtype",xtype, "ndims",ndims, "natts",natts
        call message(msg)
    endif
    do i=1,ndims
        call check(nf90_inquire_dimension(ncid, dimids(i), name, len),"nf90_inquire_dimension")
        dims(i)=len
        if(verbose)then
            write(msg,*)"dimension ",i,trim(name)," length",len
            call message(msg)
        endif
    enddo
    if(.not.verbose)return
    do i=1,natts
        attnum = i
        call check(nf90_inq_attname(ncid, varid, attnum, name),"nf90_inq_attname")
        call check(nf90_inquire_attribute(ncid, varid, trim(name), xtype, len, attnum),"nf90_inquire_attribute")
        if(len>mdims)call crash("netcdf_var_info: increase mdims")
        !write(msg,*)"attribute ",i,trim(name),' type',xtype
        !call message(msg)
        select case (xtype) 
            case (nf90_char)
                call check(nf90_get_att(ncid, varid, trim(name), values_char),"nf90_get_att")
                write(msg,*)"attribute ",i,trim(name)," type ",xtype," values",len," : ",trim(values_char)
            case (nf90_int,nf90_short,nf90_ushort,nf90_uint,nf90_int64,nf90_uint64)
                call check(nf90_get_att(ncid, varid, trim(name), values_int),"nf90_get_att")
                write(msg,*)"attribute ",i,trim(name)," type ",xtype," values",len," : ",(values_int(j),j=1,len)
            case (nf90_float,nf90_double)
                call check(nf90_get_att(ncid, varid, trim(name), values_real),"nf90_get_att")
                write(msg,*)"attribute ",i,trim(name)," type ",xtype," values",len," : ",(values_real(j),j=1,len)
            case default
                write(msg,*)'attribute type ',xtype,' not supported'
        end select
        call message(msg)
    enddo
end subroutine netcdf_var_info

subroutine check(ierr,errmsg)
    implicit none
    integer, intent(in)::ierr
    character(len=*), intent(in)::errmsg
    character(len=256)msg
    if(ierr.ne.0)then
        write(msg,"(a,a,i6,1x,a)")trim(errmsg)," error",ierr,trim(nf90_strerror(ierr))
        call crash(trim(msg))
    endif
end subroutine check

end module module_netcdf
