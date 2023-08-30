module grid_vars
    use mod_types, only: wp => dp
    implicit none

    ! grid arrays
    real(wp), allocatable :: xi(:, :), et(:, :), xn(:, :), yn(:, :)

    ! node counters
    integer(wp) :: in, jn
    ! domain parameters
    integer(wp) :: yn_min, yn_max, fn_beg, fn_end

end module grid_vars