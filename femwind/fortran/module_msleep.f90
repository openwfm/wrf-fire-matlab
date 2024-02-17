module module_msleep

contains

subroutine call_msleep(m)
    ! sleep m miliseconds
    use, intrinsic :: iso_c_binding, only: c_long

    implicit none

    interface
        ! Correctly declare the C function as a subroutine in Fortran
        subroutine msleep(msec) bind(C, name="msleep")
            import :: c_long
            integer(c_long), value :: msec
        end subroutine msleep
    end interface

    integer :: m                    ! argument
    integer(c_long) :: msec_c_long  ! local

    ! Call the C subroutine
    msec_c_long = int(m, kind=c_long)
    call msleep(msec_c_long) ! Sleep for 1000 milliseconds or 1 second
end subroutine call_msleep

end module module_msleep
