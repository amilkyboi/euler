subroutine bcwall()
    use mod_types, only: wp => dp
    use flowprop,  only: r, u, v, vel, p, E, T, c, s, mach
    use reference, only: a_ref, cv_ref, T_ref, R_ref
    use gridprop,  only: ic, ic_max, jc_max
    use gasprop,   only: gamma, gammam1
    use fluxes,    only: q
    use functions
    implicit none

    ! Applies wall boundary conditions to all velocity terms in the state vector and the two flux
    ! vectors. Ensures that the velocity at the wall is tangent to the surface, thereby enforcing
    ! the no-penetration boundary condition.

    real(wp) :: nn(2), ns(2), dotn, dots

    do ic = -1, ic_max + 2
        ! upward facing normal of ghost cells on bottom wall
        nn = normal(ic, 0, 1)

        dotn = q(ic, 1, 2) * nn(1) + q(ic, 1, 3) * nn(2)

        ! first layer of ghost cells on the bottom
        ! update state vector
        q(ic, 0, 1) = q(ic, 1, 1)
        q(ic, 0, 2) = q(ic, 1, 2) - 2 * dotn * nn(1)
        q(ic, 0, 3) = q(ic, 1, 3) - 2 * dotn * nn(2)
        q(ic, 0, 4) = q(ic, 1, 4)

        ! update tracked quantities
        r(ic, 0) = q(ic, 0, 1)
        u(ic, 0) = q(ic, 0, 2) / q(ic, 0, 1)
        v(ic, 0) = q(ic, 0, 3) / q(ic, 0, 1)
        E(ic, 0) = q(ic, 0, 4) / q(ic, 0, 1)

        vel(ic, 0) = sqrt(u(ic, 0)**2 + v(ic, 0)**2)
        T(ic, 0) = ((E(ic, 0) * a_ref**2 - 0.5 * (vel(ic, 0) * a_ref)**2) / cv_ref) / T_ref
        c(ic, 0) = sqrt(gamma * R_ref * T(ic, 0) * T_ref) / a_ref
        mach(ic, 0) = vel(ic, 0) / c(ic, 0)
        p(ic, 0) = gammam1 * r(ic, 0) * (E(ic, 0) - 0.5 * vel(ic, 0)**2)
        s(ic, 0) = p(ic, 0) / r(ic, 0)**gamma

        ! second layer of ghost cells on the bottom
        ! update state vector
        q(ic, -1, 1) = q(ic, 2, 1)
        q(ic, -1, 2) = q(ic, 2, 2) - 2 * dotn * nn(1)
        q(ic, -1, 3) = q(ic, 2, 3) - 2 * dotn * nn(2)
        q(ic, -1, 4) = q(ic, 2, 4)

        ! update tracked quantities
        r(ic, -1) = q(ic, -1, 1)
        u(ic, -1) = q(ic, -1, 2) / q(ic, -1, 1)
        v(ic, -1) = q(ic, -1, 3) / q(ic, -1, 1)
        E(ic, -1) = q(ic, -1, 4) / q(ic, -1, 1)

        vel(ic, -1) = sqrt(u(ic, -1)**2 + v(ic, -1)**2)
        T(ic, -1) = ((E(ic, -1) * a_ref**2 - 0.5 * (vel(ic, -1) * a_ref)**2) / cv_ref) / T_ref
        c(ic, -1) = sqrt(gamma * R_ref * T(ic, -1) * T_ref) / a_ref
        mach(ic, -1) = vel(ic, -1) / c(ic, -1)
        p(ic, -1) = gammam1 * r(ic, -1) * (E(ic, -1) - 0.5 * vel(ic, -1)**2)
        s(ic, -1) = p(ic, -1) / r(ic, -1)**gamma

        ! bottom facing normal of ghost cells on top wall
        ns = normal(ic, jc_max+1, 2)

        dots = q(ic, jc_max, 2) * ns(1) + q(ic, jc_max, 3) * ns(2)

        ! first layer of ghost cells on the top
        ! update state vector
        q(ic, jc_max+1, 1) = q(ic, jc_max, 1)
        q(ic, jc_max+1, 2) = q(ic, jc_max, 2) - 2 * dots * ns(1)
        q(ic, jc_max+1, 3) = q(ic, jc_max, 3) - 2 * dots * ns(2)
        q(ic, jc_max+1, 4) = q(ic, jc_max, 4)

        ! update tracked quantities
        r(ic, jc_max+1) = q(ic, jc_max+1, 1)
        u(ic, jc_max+1) = q(ic, jc_max+1, 2) / q(ic, jc_max+1, 1)
        v(ic, jc_max+1) = q(ic, jc_max+1, 3) / q(ic, jc_max+1, 1)
        E(ic, jc_max+1) = q(ic, jc_max+1, 4) / q(ic, jc_max+1, 1)

        vel(ic, jc_max+1) = sqrt(u(ic, jc_max+1)**2 + v(ic, jc_max+1)**2)
        T(ic, jc_max+1) = ((E(ic, jc_max+1) * a_ref**2 - 0.5 * (vel(ic, jc_max+1) * a_ref)**2) / cv_ref) / T_ref
        c(ic, jc_max+1) = sqrt(gamma * R_ref * T(ic, jc_max+1) * T_ref) / a_ref
        mach(ic, jc_max+1) = vel(ic, jc_max+1) / c(ic, jc_max+1)
        p(ic, jc_max+1) = gammam1 * r(ic, jc_max+1) * (E(ic, jc_max+1) - 0.5 * vel(ic, jc_max+1)**2)
        s(ic, jc_max+1) = p(ic, jc_max+1) / r(ic, jc_max+1)**gamma

        ! second layer of ghost cells on the top
        ! update state vector
        q(ic, jc_max+2, 1) = q(ic, jc_max-1, 1)
        q(ic, jc_max+2, 2) = q(ic, jc_max-1, 2) - 2 * dots * ns(1)
        q(ic, jc_max+2, 3) = q(ic, jc_max-1, 3) - 2 * dots * ns(2)
        q(ic, jc_max+2, 4) = q(ic, jc_max-1, 4)

        ! update tracked quantities
        r(ic, jc_max+2) = q(ic, jc_max+2, 1)
        u(ic, jc_max+2) = q(ic, jc_max+2, 2) / q(ic, jc_max+2, 1)
        v(ic, jc_max+2) = q(ic, jc_max+2, 3) / q(ic, jc_max+2, 1)
        E(ic, jc_max+2) = q(ic, jc_max+2, 4) / q(ic, jc_max+2, 1)

        vel(ic, jc_max+2) = sqrt(u(ic, jc_max+2)**2 + v(ic, jc_max+2)**2)
        T(ic, jc_max+2) = ((E(ic, jc_max+2) * a_ref**2 - 0.5 * (vel(ic, jc_max+2) * a_ref)**2) / cv_ref) / T_ref
        c(ic, jc_max+2) = sqrt(gamma * R_ref * T(ic, jc_max+2) * T_ref) / a_ref
        mach(ic, jc_max+2) = vel(ic, jc_max+2) / c(ic, jc_max+2)
        p(ic, jc_max+2) = gammam1 * r(ic, jc_max+2) * (E(ic, jc_max+2) - 0.5 * vel(ic, jc_max+2)**2)
        s(ic, jc_max+2) = p(ic, jc_max+2) / r(ic, jc_max+2)**gamma
    end do
end subroutine bcwall