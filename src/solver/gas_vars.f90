module gas_vars
    use mod_types, only: wp => dp
    implicit none

    ! declares all variables related to the working fluid

    ! ratio of specific heats
    real(wp), parameter :: gamma            = 1.4_wp
    ! gamma - 1
    real(wp), parameter :: gamm1            = 0.4_wp
    ! 1 / gamma
    real(wp), parameter :: over_gamma       = 1.0_wp / 1.4_wp
    ! 1 / (gamma - 1)
    real(wp), parameter :: over_gamm1       = 1.0_wp / 0.4_wp
    ! 1 / (gamma * (gamma - 1))
    real(wp), parameter :: over_gtgm1       = 1.0_wp / 0.56_wp
    ! gamma / (gamma - 1)
    real(wp), parameter :: gamma_over_gamm1 = 1.4_wp / 0.4_wp

end module gas_vars