module input
    use mod_types, only: wp => dp
    implicit none

    ! All relevant input parameters specified by the user.

    ! freestream mach number
    real(wp), parameter :: mach_inf = 0.7_wp
    ! maximum number of RK4 iterations
    integer, parameter :: n_iter = 5000
    ! second-order switch magnitude
    real(wp), parameter :: nu2 = 0.4_wp
    ! fourth-order switch magnitude
    real(wp), parameter :: nu4 = 0.003_wp
    ! Courant–Friedrichs–Lewy number
    real(wp), parameter :: cfl = 1.0_wp
    ! output file name for saving plots
    character(len=30), parameter :: plt_str = '../../data/soln.plt'
    ! output file name for saving residual
    character(len=30), parameter :: res_str = '../../data/res.csv'
end module input
