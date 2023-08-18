module timing
    use mod_types, only: wp => dp
    implicit none

    ! Declares global timing variables used throughout multiple subroutines.

    ! start time, end time
    real(wp) :: start, end
end module timing