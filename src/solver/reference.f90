module reference
    use mod_types, only: wp => dp
    use gas_vars,  only: gamma, gammam1
    implicit none

    ! Declares all variables that serve as standard reference values for non-dimensionalization and
    ! re-dimensionalization.

    ! reference length, 1 m
    real(wp), parameter :: l_ref = 1.0_wp
    ! reference specific gas constant for air, 287 J/(kg*K)    
    real(wp), parameter :: R_ref = 287.0_wp
    ! reference temperature, 298 K
    real(wp), parameter :: T_ref = 298.0_wp
    ! reference speed of sound, m/s
    real(wp), parameter :: a_ref = sqrt(gamma * 287.0_wp * 298.0_wp)
    ! reference specific heat capacity at constant volume, J/(kg*K)
    real(wp), parameter :: cv_ref = 287.0_wp / gammam1

end module reference