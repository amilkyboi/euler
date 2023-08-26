module input
    use mod_types,  only: wp => dp
    use free_vars, only: p_inf
    implicit none

    ! maximum number of RK4 iterations
    integer, parameter :: n_iter = 5000

    ! freestream mach number
    real(wp), parameter :: mach_inf = 0.3_wp
    ! second-order switch magnitude
    real(wp), parameter :: nu2      = 0.0_wp
    ! fourth-order switch magnitude
    real(wp), parameter :: nu4      = 0.002_wp
    ! Courant–Friedrichs–Lewy number
    real(wp), parameter :: cfl      = 1.0_wp
    ! exit pressure
    real(wp), parameter :: p_exit   = 0.1_wp * p_inf

    ! output file name for saving plots
    character(len=19), parameter :: plt_str = '../../data/soln.plt'
    ! output file name for saving residual
    character(len=19), parameter :: res_str = '../../data/resd.csv'

end module input
