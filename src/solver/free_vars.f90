module free_vars
    use mod_types, only: wp => dp
    use gas_vars,   only: gamma
    implicit none

    ! Declares all variables that involve the freestream flow conditions.

    ! freestream stagnation pressure
    real(wp), parameter :: p_inf = 1.0_wp / gamma

end module free_vars