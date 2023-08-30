module timing
    use mod_types, only: wp => dp
    implicit none

    ! needs to be of type INT8 otherwise the timing accuracy only goes to around E-2
    integer(wp) :: start, end

    real(wp) :: rate

end module timing