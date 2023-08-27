module input
    use mod_types, only: wp => dp
    implicit none

    ! maximum number of iterations for the elliptic Gauss-Seidel method
    integer(wp), parameter :: max_iter = 1e5
    ! maximum number of nodes in the x-direction
    integer(wp), parameter :: in_max   = 51
    ! maximum number of nodes in the y-direction
    integer(wp), parameter :: jn_max   = 21

    ! tolerance for the elliptic Gauss-Seidel method
    real(wp), parameter :: tol          = 10e-10
    ! x-limits for the real domain
    real(wp), parameter :: x_bounds(2)  = (/0, 5/)
    ! y-limits for the real domain
    real(wp), parameter :: y_bounds(2)  = (/0, 1/)
    ! x-limits between which the functions are applied
    real(wp), parameter :: fn_bounds(2) = (/2, 3/)

    ! output algebraic grid file
    logical, parameter :: save_alg = .true.
    ! output elliptic grid file
    logical, parameter :: save_elp = .true.

    ! save location for algebraic grid file
    character(len=19), parameter :: alg_path = '../../data/xy_alg.x'
    ! save location for elliptic grid file
    character(len=19), parameter :: elp_path = '../../data/xy_elp.x'

end module input