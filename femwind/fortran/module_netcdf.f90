module module_netcdf

use netcdf
use module_utils

contains

subroutine netcdf_read_int(ncid,ia,varname)
implicit none
!*** arguments
    integer, intent(in)::ncid                      ! open pnetcdf file
    integer, intent(inout)::ia                     ! variable to store
    character(LEN=*),intent(in)::varname           ! variable name
!*** local
    integer::ierr,varid
    character(len=256)::msg
!*** executable
        call check(nf90_inq_varid(ncid, trim(varname), varid), &
            "netcdf_read_int/nf90_inq_varid:"//trim(varname))
        call check(nf90_get_var(ncid, varid, ia), &
            "netcdf_read_int/nf90_get_var:"//trim(varname))
        write(msg,*)'netcdf_read_int: varname=',varname,' value=',ia
end subroutine netcdf_read_int

subroutine netcdf_write_int(ncid,ia,varname)
implicit none
!*** arguments
    integer, intent(in)::                         &
    ncid,                                         & ! open pnetcdf file
    ia                                              ! variable to write
    character(LEN=*),intent(in):: varname
!*** local
    integer::varid,ival
    character(len=256)::msg
!*** executable
    write(msg,*)'netcdf_write_int: varname=',varname,' value=',ia
    call message(msg,level=0)
        call check(nf90_inq_varid(ncid, trim(varname), varid), &
            "netcdf_write_int/nf90_inq_varid:"//trim(varname))
        ival = ia
        call check(nf90_put_var(ncid, varid, ival), &
            "netcdf_write_int/nf90mpi_put_var:"//trim(varname))
end subroutine netcdf_write_int

subroutine pnetcdf_write_arr(ncid,                &
    ids,ide, kds,kde, jds,jde,                    & ! atm grid dimensions
    ims,ime, kms,kme, jms,jme,                    &
    ips,ipe, kps,kpe, jps,jpe,                    &
    a,name)
    implicit none

!*** arguments
    integer, intent(in)::                         &
    ncid,                                         & ! open pnetcdf file
    ids,ide, kds,kde, jds,jde,                    & ! atm grid dimensions
    ims,ime, kms,kme, jms,jme,                    &
    ips,ipe, kps,kpe, jps,jpe                     
    real,intent(in),dimension(ims:ime,kms:kme,jms:jme)::a  
    character(LEN=*),intent(in):: name

!*** local
    integer(kind=MPI_OFFSET_KIND),dimension(4)::star,cnts
    integer::i,j,k,varid,ends(4),dims(4)
    real,dimension(:,:,:,:),allocatable::at
    character(len=256) msg

    ! get idx
    call pnetcdf_var_info(ncid,name,dims,varid,1)
    star   = (/ips,jps,kps,1/)
    ends   = (/ipe,jpe,kpe,1/)
    ends   = min(ends,dims)
    ! at end of domain, extend patch by one
    if (ends(1).eq.dims(1)-1)ends(1)=dims(1) 
    if (ends(2).eq.dims(2)-1)ends(2)=dims(2)
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

    write(msg,*)"writing ",trim(name),star(1),ends(1),star(2),ends(2),star(3),ends(3)
    call message(msg)
 
    ! write to file
    call check(nf90mpi_put_var_all(ncid, varid, at, start = star, count = cnts),"nf90mpi_put_var:"//trim(name))

    deallocate(at)
    
end subroutine pnetcdf_write_arr

subroutine pnetcdf_read_arr(ncid,                &
    ids,ide, kds,kde, jds,jde,                    & ! atm grid dimensions
    ims,ime, kms,kme, jms,jme,                    &
    ips,ipe, kps,kpe, jps,jpe,                    &
    a,name)
    implicit none

!*** arguments
    integer, intent(in)::                         &
    ncid,                                         & ! open pnetcdf file
    ids,ide, kds,kde, jds,jde,                    & ! atm grid dimensions
    ims,ime, kms,kme, jms,jme,                    &
    ips,ipe, kps,kpe, jps,jpe                     
    real,intent(out),dimension(ims:ime,kms:kme,jms:jme)::a  
    character(LEN=*),intent(in):: name

!*** local
    integer(kind=MPI_OFFSET_KIND),dimension(4)::star,cnts
    integer::i,j,k,varid,ends(4),dims(4)
    real,dimension(:,:,:,:),allocatable::at
    character(len=256) msg

    ! get idx
    call pnetcdf_var_info(ncid,name,dims,varid,1)
    star   = (/ips,jps,kps,1/)
    ends   = (/ipe,jpe,kpe,1/)
    ends   = min(ends,dims)
    ! at end of domain, extend patch by one
    if (ends(1).eq.dims(1)-1)ends(1)=dims(1) 
    if (ends(2).eq.dims(2)-1)ends(2)=dims(2)
    cnts = ends - star + 1

    write(msg,*)"reading ",trim(name),star(1),ends(1),star(2),ends(2),star(3),ends(3)
    call message(msg)
 
    ! read from file
    allocate(at(star(1):ends(1),star(2):ends(2),star(3):ends(3),1))
    call check(nf90mpi_get_var_all(ncid, varid, at, start = star, count = cnts),"nf90mpi_get_var:"//trim(name))
    
    ! transpose at -> a
    do k=star(3),ends(3)
        do j=star(2),ends(2)
            do i=star(1),ends(1)
                a(i,k,j) = at(i,j,k,1)
            enddo
        enddo
    enddo

    deallocate(at)
    
end subroutine pnetcdf_read_arr

subroutine pnetcdf_var_info(ncid,varname,dims,varid,prints)
    implicit none
!*** arguments
    integer, intent(in)::ncid
    character(len=*)::varname
    integer,intent(out)::dims(:),varid
    integer,intent(in),optional::prints 
!*** local
    integer, parameter::mdims = 256
    integer:: xtype, ndims, natts, dimids(mdims),i,j,attnum
    integer(kind=MPI_OFFSET_KIND) :: len
    character(len=nf90_max_name):: name
    integer:: values_int(mdims)
    real:: values_real(mdims)
    character(len=mdims):: values_char
    character(LEN=256):: filename, msg
    logical::verbose=.true.

    if(present(prints)) verbose = prints>0
     

    call check(nf90mpi_inq_varid(ncid,trim(varname),varid),"nf90mpi_inq_varid"//trim(varname))
    call check(nf90mpi_inquire_variable(ncid, varid, name, xtype, ndims, dimids, natts),"nf90mpi_inquire_variable")
    if(ndims>mdims)call crash("pnetcdf_var_info: increase mdims")
    if(ndims>size(dims))call crash("pnetcdf_var_info: dims too short")
    if(verbose)then
        write(msg,*)"variable ",trim(name), " xtype",xtype, "ndims",ndims, "natts",natts
        call message(msg)
    endif
    do i=1,ndims
        call check(nf90mpi_inquire_dimension(ncid, dimids(i), name, len),"nf90mpi_inquire_dimension")
        dims(i)=len
        if(verbose)then
            write(msg,*)"dimension ",i,trim(name)," length",len
            call message(msg)
        endif
    enddo
    if(.not.verbose)return
    do i=1,natts
        attnum = i
        call check(nf90mpi_inq_attname(ncid, varid, attnum, name),"nf90mpi_inq_attname")
        call check(nf90mpi_inquire_attribute(ncid, varid, trim(name), xtype, len, attnum),"nf90mpi_inquire_attribute")
        if(len>mdims)call crash("pnetcdf_var_info: increase mdims")
        !write(msg,*)"attribute ",i,trim(name),' type',xtype
        !call message(msg)
        select case (xtype) 
            case (nf90_char)
                call check(nf90mpi_get_att(ncid, varid, trim(name), values_char),"nf90mpi_get_att")
                write(msg,*)"attribute ",i,trim(name)," type ",xtype," values",len," : ",trim(values_char)
            case (nf90_int,nf90_short,nf90_ushort,nf90_uint,nf90_int64,nf90_uint64)
                call check(nf90mpi_get_att(ncid, varid, trim(name), values_int),"nf90mpi_get_att")
                write(msg,*)"attribute ",i,trim(name)," type ",xtype," values",len," : ",(values_int(j),j=1,len)
            case (nf90_float,nf90_double)
                call check(nf90mpi_get_att(ncid, varid, trim(name), values_real),"nf90mpi_get_att")
                write(msg,*)"attribute ",i,trim(name)," type ",xtype," values",len," : ",(values_real(j),j=1,len)
            case default
                write(msg,*)'attribute type ',xtype,' not supported'
        end select
        call message(msg)
    enddo
end subroutine pnetcdf_var_info

subroutine check(ierr,errmsg)
    implicit none
    integer, intent(in)::ierr
    character(len=*), intent(in)::errmsg
    character(len=256)msg
    if(ierr.ne.0)then
        write(msg,"(a,a,i6,1x,a)")trim(errmsg)," error",ierr,trim(nf90mpi_strerror(ierr))
        call crash(trim(msg))
    endif
end subroutine check
#endif

module module_netcdf
