module timing
    use mod_types, only: wp => dp
    implicit none

    ! start and end of system time
    integer(wp) :: start, end

    ! rate at which time is polled
    real(wp) :: rate

end module timing