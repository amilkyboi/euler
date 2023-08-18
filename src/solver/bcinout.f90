subroutine bcinout()
    use mod_types,  only: wp => dp
    use flowprop,   only: r, p, u, v, vel, c, s, E, T, mach
    use freestream, only: p_inf, riem1_inf
    use reference,  only: a_ref, cv_ref, T_ref
    use gridprop,   only: jc, ic_max, jc_max
    use gasprop,    only: gamma, gammam1
    use input,      only: mach_inf
    use fluxes,     only: q
    use functions
    implicit none

    ! Applies inlet and outlet boundary conditions to the state vector at ic = 1, and i = ic_max
    ! respectively. 

    real(wp) :: riem1, riem2_1, riem2_2, alpha_inlet, p0_inf, T0_inf, T_inf

    ! inlet boundary conditions
    alpha_inlet = 0.0_wp

    ! freestream stagnation temperature
    T_inf = 1.0_wp
    T0_inf = T_inf * (1.0_wp + 0.5 * gammam1 * mach_inf**2)

    ! freestream stagnation pressure
    p0_inf = p_inf * (1.0_wp + 0.5 * gammam1 * mach_inf**2)**(gamma/gammam1)
    

    ! 1. find riem1 at freestream and calculate riem1 at i = 1
    riem1_inf = mach_inf + 2 / gammam1
    riem1 = riem1_inf

    do jc = -1, jc_max + 2
        ! 2. find riem2 at i = 2 and calculate riem2 at i = 1
        riem2_2 = vel(2, jc) - (2 * c(2, jc)) / gammam1
        riem2_1 = riem2_2

        ! 3. find velocity at i = 1
        vel(1, jc) = 0.5 * (riem1 + riem2_1)

        ! 4. find u and v at i = 1
        u(1, jc) = vel(1, jc) * cos(alpha_inlet)
        v(1, jc) = vel(1, jc) * sin(alpha_inlet)

        ! 5. find speed of sound at i = 1
        c(1, jc) = 0.25 * gammam1 * (riem1 - riem2_1)

        ! 6. find mach number at i = 1
        mach(1, jc) = vel(1, jc) / c(1, jc)

        ! 7. find static pressure at i = 1
        p(1, jc) = p0_inf / (1.0_wp + 0.5 * gammam1 * mach(1, jc)**2)**(gamma/gammam1)

        ! 8. find density at i = 1
        r(1, jc) = (gamma * p(1, jc)) / c(1, jc)**2

        ! 9. find energy at i = 1
        E(1, jc) = p(1, jc) / (r(1, jc) * gammam1) + 0.5 * vel(1, jc)**2

        ! 10. find T and s
        T(1, jc) = T0_inf / (1.0_wp + 0.5 * gammam1 * mach(1, jc)**2)
        s(1, jc) = p(1, jc) / r(1, jc)**gamma

        ! 10. update state vector components
        q(1, jc, 1) = r(1, jc)
        q(1, jc, 2) = r(1, jc) * u(1, jc)
        q(1, jc, 3) = r(1, jc) * v(1, jc)
        q(1, jc, 4) = r(1, jc) * E(1, jc)
    end do

    do jc = -1, jc_max + 2
        ! outlet boundary conditions
        ! 1. static pressure at i = i_max
        p(ic_max, jc) = p(ic_max-1, jc)
        ! if (mach(ic_max, jc) >= 1.0_wp) then
        !     ! supersonic case
        !     p(ic_max, jc) = p(ic_max-1, jc)
        ! else
        !     ! subsonic case
        !     p(ic_max, jc) = p_inf
        ! end if

        ! 2. entropy at i = i_max
        s(ic_max, jc) = s(ic_max-1, jc)

        ! 3. riem1 at i = i_max
        riem1 = vel(ic_max-1, jc) + (2 * c(ic_max-1, jc)) / gammam1

        ! 4. y-velocity at i = i_max
        v(ic_max, jc) = v(ic_max-1, jc)

        ! 5. other variables at i = i_max
        r(ic_max, jc) = (p(ic_max, jc) / s(ic_max, jc))**(1 / gamma)
        c(ic_max, jc) = sqrt((gamma * p(ic_max, jc)) / r(ic_max, jc))
        vel(ic_max, jc) = riem1 - (2 * c(ic_max, jc)) / gammam1
        u(ic_max, jc) = sqrt(vel(ic_max, jc)**2 - v(ic_max, jc)**2)
        E(ic_max, jc) = p(ic_max, jc) / (r(ic_max, jc) * gammam1) + 0.5 * vel(ic_max, jc)**2
        T(ic_max, jc) = ((E(ic_max, jc) * a_ref**2 - 0.5 * (vel(ic_max, jc) * a_ref)**2) / cv_ref) / T_ref

        ! 6. update state vector components
        q(ic_max, jc, 1) = r(ic_max, jc)
        q(ic_max, jc, 2) = r(ic_max, jc) * u(ic_max, jc)
        q(ic_max, jc, 3) = r(ic_max, jc) * v(ic_max, jc)
        q(ic_max, jc, 4) = r(ic_max, jc) * E(ic_max, jc)
    end do

end subroutine bcinout