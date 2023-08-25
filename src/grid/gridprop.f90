module gridprop
    use mod_types, only: wp => dp
    implicit none

    ! arrays
    real(wp), allocatable :: xn(:, :), yn(:, :), xi(:, :), et(:, :)

    ! nodes
    integer(wp) :: in, jn
    ! domain
    integer(wp) :: yn_min, yn_max, fn_beg, fn_end
    ! timing
    integer(wp) :: start, end

    ! timing
    real(wp) :: rate
    

end module gridprop