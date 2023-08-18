module gridprop
    use mod_types, only: wp => dp
    implicit none

    integer :: in, jn, n_max, in_max, jn_max, xy_bounds(4), fn_bounds(2), yn_min, yn_max, fn_beg, fn_end
    real(wp), allocatable :: xn(:, :), yn(:, :), xi(:, :), et(:, :)
    real(wp) :: start, end, tol

end module gridprop