subroutine bcinout()
    use mod_types,  only: wp => dp
    use flowvars,   only: dens, pres, xvel, yvel, vmag, vsnd, entr, enrg, temp, mach
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

    integer :: k

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
        riem2_2 = vmag(2, jc) - (2 * vsnd(2, jc)) / gammam1
        riem2_1 = riem2_2

        ! 3. find velocity at i = 1
        vmag(1, jc) = 0.5 * (riem1 + riem2_1)

        ! 4. find u and v at i = 1
        xvel(1, jc) = vmag(1, jc) * cos(alpha_inlet)
        yvel(1, jc) = vmag(1, jc) * sin(alpha_inlet)

        ! 5. find speed of sound at i = 1
        vsnd(1, jc) = 0.25 * gammam1 * (riem1 - riem2_1)

        ! 6. find mach number at i = 1
        mach(1, jc) = vmag(1, jc) / vsnd(1, jc)

        ! 7. find static pressure at i = 1
        pres(1, jc) = p0_inf / (1.0_wp + 0.5 * gammam1 * mach(1, jc)**2)**(gamma/gammam1)

        ! 8. find density at i = 1
        dens(1, jc) = (gamma * pres(1, jc)) / vsnd(1, jc)**2

        ! 9. find energy at i = 1
        enrg(1, jc) = pres(1, jc) / (dens(1, jc) * gammam1) + 0.5 * vmag(1, jc)**2

        ! 10. find T and s
        temp(1, jc) = T0_inf / (1.0_wp + 0.5 * gammam1 * mach(1, jc)**2)
        entr(1, jc) = pres(1, jc) / dens(1, jc)**gamma

        ! 10. update state vector components
        q(1, jc, 1) = dens(1, jc)
        q(1, jc, 2) = dens(1, jc) * xvel(1, jc)
        q(1, jc, 3) = dens(1, jc) * yvel(1, jc)
        q(1, jc, 4) = dens(1, jc) * enrg(1, jc)
    end do

    do jc = -1, jc_max + 2
        ! outlet boundary conditions
        ! 1. static pressure at i = i_max
        if (mach(ic_max, jc) >= 1.0_wp) then
            ! supersonic case
            pres(ic_max, jc) = pres(ic_max - 1, jc)

            q(ic_max + 1, jc, 4) = pres(ic_max, jc) / gammam1 + 0.5 * q(ic_max + 1, jc, 1) * vmag(ic_max + 1, jc)**2
        else
            ! subsonic case
            pres(ic_max, jc) = p_inf

            q(ic_max + 1, jc, 4) = p_inf / gammam1 + 0.5 * q(ic_max + 1, jc, 1) * vmag(ic_max + 1, jc)**2
        end if

        ! 2. entropy at i = i_max
        entr(ic_max, jc) = entr(ic_max-1, jc)

        ! 3. riem1 at i = i_max
        riem1 = vmag(ic_max-1, jc) + (2 * vsnd(ic_max-1, jc)) / gammam1

        ! 4. y-velocity at i = i_max
        yvel(ic_max, jc) = yvel(ic_max-1, jc)

        ! 5. other variables at i = i_max
        dens(ic_max, jc) = (pres(ic_max, jc) / entr(ic_max, jc))**(1 / gamma)
        vsnd(ic_max, jc) = sqrt((gamma * pres(ic_max, jc)) / dens(ic_max, jc))
        vmag(ic_max, jc) = riem1 - (2 * vsnd(ic_max, jc)) / gammam1
        xvel(ic_max, jc) = sqrt(vmag(ic_max, jc)**2 - yvel(ic_max, jc)**2)
        enrg(ic_max, jc) = pres(ic_max, jc) / (dens(ic_max, jc) * gammam1) + 0.5 * vmag(ic_max, jc)**2
        temp(ic_max, jc) = ((enrg(ic_max, jc) * a_ref**2 - 0.5 * (vmag(ic_max, jc) * a_ref)**2) / cv_ref) / T_ref

        ! 6. update state vector components
        q(ic_max, jc, 1) = dens(ic_max, jc)
        q(ic_max, jc, 2) = dens(ic_max, jc) * xvel(ic_max, jc)
        q(ic_max, jc, 3) = dens(ic_max, jc) * yvel(ic_max, jc)
        q(ic_max, jc, 4) = dens(ic_max, jc) * enrg(ic_max, jc)

        do k = 1, 3
            q(ic_max + 1, jc, k) = 2 * q(ic_max, jc, k) - q(ic_max - 1, jc, k)
        end do
    end do

end subroutine bcinout