module free_vars
    use mod_types, only: wp => dp
    use gas_vars,   only: gamma
    implicit none

    ! Declares all variables that involve the freestream flow conditions.

    ! freestream stagnation pressure
    real(wp), parameter :: p_inf = 1.0_wp / gamma
    ! freestream x-velocity, freestream y-velocity, freestream Riemann1 invariant
    real(wp) :: u_inf, v_inf, riem1_inf
end module free_vars