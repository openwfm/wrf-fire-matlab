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

integer function l2i(l)
    implicit none
    logical, intent(in)::l
    if(l)then
        l2i = 1
    else
        l2i = 0
    endif
end function l2i

subroutine netcdf_read_array_wrf(ncid,name,frame,sr,a1d,a2d,a3d)
    implicit none

!*** arguments
    integer, intent(in)::                ncid ! open netcdf file
    real, pointer, intent(out),optional:: a1d(:)  ! the array pointer; remember to deallocate when done with it
    real, pointer, intent(out),optional:: a2d(:,:)  ! the array pointer; remember to deallocate when done with it
    real, pointer, intent(out),optional:: a3d(:,:,:)  ! the array pointer; remember to deallocate when done with it
    character(LEN=*),intent(in):: name
    integer, intent(in), optional::                frame ! time step in the file
    integer, intent(in), dimension(2), optional::  sr ! strip to remove in i and j dimension

!*** local
    integer d1,d2,d3
    integer,dimension(:),allocatable::star,cnts,dims,ends
    integer::i,j,k,varid,ndims
    real,dimension(:,:),allocatable::a1
    real,dimension(:,:,:),allocatable::a2
    real,dimension(:,:,:,:),allocatable::a3
    integer:: istep=1, srf(2)
    character(len=256) msg

!*** executable
    if(present(frame))istep=frame

    d1 = l2i(present(a1d))
    d2 = l2i(present(a2d))
    d3 = l2i(present(a3d))

    if(d1+d2+d3.ne.1)call crash('netcdf_read_array_wrf: must have exactly one of a1d a2d a3d arguments')
    ndims = d1+2*d2+3*d3

    if(ndims.eq.1)then
        srf = 0
    elseif(present(sr))then
        srf = sr
    else
        call get_sr(ncid,srf)
    endif

    allocate(star(ndims+1))
    allocate(cnts(ndims+1))
    allocate(dims(ndims+1))
    allocate(ends(ndims+1))

    write(msg,*)"netcdf_read_array_wrf reading variable ",trim(name), &
          " dimension ",ndims," time step ",istep
    call message(msg)

    ! get idx
    call netcdf_var_info(ncid,name,dims,varid,netcdf_msglevel)
    write(msg,*)"got dimensions ",trim(name),dims
    call message(msg)
    if(dims(ndims+1).eq.0)then
       write(msg,*)'netcdf_read_array_wrf: wrong dimension ndims=',ndims, ' 1dims=',dims
       call crash(msg)
    endif
     
    star   = 1
    star(ndims+1) = istep
    ends(1:ndims) = dims(1:ndims)
    ends(ndims+1) = istep
    cnts = ends - star + 1

    write(msg,*)"reading ",trim(name),star," to ",ends
    call message(msg)
 
    ! read from file
    select case (ndims)
    case(1) 
        allocate(a1(star(1):ends(1),star(2):ends(2)))
        allocate(a1d(star(1):ends(1)))
        call check(nf90_get_var(ncid, varid, a1d, start = star, count = cnts),"nf90_get_var:"//trim(name))
        a1d = a1(:,istep)
        write(msg,*)"returning 1D array length ",shape(a1d)
    case(2) 
	allocate(a2(star(1):ends(1),star(2):ends(2),star(3):ends(3)))
        call check(nf90_get_var(ncid, varid, a2, start = star, count = cnts),"nf90_get_var:"//trim(name))
        ! postprocessing - strip
        write(msg,*)" stripping ",srf," at ij ends"
        call message(msg)
        allocate(a2d(star(1):ends(1)-srf(1),star(2):ends(2)-srf(2)))
        do j=star(2),ends(2)-srf(2)
            do i=star(1),ends(1)-srf(1)
                a2d(i,j) = a2(i,j,istep)
            enddo
        enddo
        write(msg,*)"returning array ij shape ",shape(a2d)
    case(3) 
        allocate(a3(star(1):ends(1),star(2):ends(2),star(3):ends(3),star(4):ends(4)))
        call check(nf90_get_var(ncid, varid, a3, start = star, count = cnts),"nf90_get_var:"//trim(name))
        ! transpose at -> a and strip
        write(msg,*)" stripping ",srf," at ij ends and transposing to ikj indexing order"
        call message(msg)
        allocate(a3d(star(1):ends(1)-srf(1),star(3):ends(3),star(2):ends(2)-srf(2)))
        do k=star(3),ends(3)
            do j=star(2),ends(2)-srf(2)
                do i=star(1),ends(1)-srf(1)
                    a3d(i,k,j) = a3(i,j,k,istep)
                enddo
            enddo
        enddo
        write(msg,*)"returning array ikj shape ",shape(a3d)
    case default
        call crash('wrong case ndims')
    end select

    call message(msg)

end subroutine netcdf_read_array_wrf

subroutine get_sr(ncid,sr)
    implicit none
!*** arguments
    integer, intent(in)::ncid    ! open wrfout dataset
    integer, intent(out), dimension(2):: sr  ! fire "subgrid" refinement factors in x and y directions
!*** local
    integer, dimension(3) :: dims_atm, dims_fire   ! 2D + allow time dimension
    integer::varid, prints
!*** executable
    prints = netcdf_msglevel
    dims_atm = 0
    dims_fire = 0
    call netcdf_var_info(ncid,"XLAT",dims_atm,varid,prints=prints)
    call netcdf_var_info(ncid,"FXLAT",dims_fire,varid,prints=prints)
    sr = dims_fire(1:2)/(dims_atm(1:2) + 1)
    if(prints >= 0) print *,"get_sr: fire subgrid refinement factors are ",sr
end subroutine get_sr

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
