module input
    use mod_types,  only: wp => dp
    use free_vars, only: p_inf
    implicit none

    ! switch & exit pressure values that have worked well in the past for reference
    ! mach 0.3: nu2 = 0.0_wp, nu4 = 0.002_wp, p_exit = 1.0_wp * p_inf
    ! mach 0.5: nu2 = 0.0_wp, nu4 = 0.003_wp, p_exit = 1.0_wp * p_inf
    ! mach 0.7: nu2 = 0.2_wp, nu4 = 0.003_wp, p_exit = 0.8_wp * p_inf

    ! maximum number of RK4 iterations
    integer, parameter :: max_iter = 5000

    ! freestream mach number
    real(wp), parameter :: mach_inf  = 0.7_wp
    ! second-order switch magnitude
    real(wp), parameter :: nu2       = 0.2_wp
    ! fourth-order switch magnitude
    real(wp), parameter :: nu4       = 0.003_wp
    ! exit pressure
    real(wp), parameter :: p_exit    = 0.8_wp * p_inf
    ! inlet angle of attack
    real(wp), parameter :: inlet_aoa = 0.0_wp
    ! maximum residual tolerance
    real(wp), parameter :: res_tol   = 10e-8_wp
    ! Courant–Friedrichs–Lewy number
    real(wp), parameter :: cfl       = 1.0_wp

    ! output file name for saving plots
    character(len=19), parameter :: plt_str = '../../data/soln.plt'
    ! output file name for saving residual
    character(len=19), parameter :: res_str = '../../data/res.csv'

end module input
