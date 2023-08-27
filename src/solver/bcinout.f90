subroutine bcinout()
    use mod_types, only: wp => dp
    use flow_vars, only: dens, pres, xvel, yvel, vmag, vsnd, entr, enrg, temp, mach
    use gas_vars,  only: gamma, gamm1, over_gamma, over_gamm1, gamma_over_gamm1
    use input,     only: mach_inf, inlet_aoa, p_exit
    use reference, only: a_ref, cv_ref, T_ref
    use grid_vars, only: jc, ic_max, jc_max
    use free_vars, only: p_inf
    use flux_vars, only: q
    use functions
    implicit none

    ! applies inlet and outlet boundary conditions to the state vector at ic = 1, and i = ic_max

    ! temporary counter
    integer :: k

    ! Riemann invariants
    real(wp) :: riem1_inf, riem1_1, riem2_1, riem2_2, stag_pres_inf, stag_temp_inf, temp_inf

    ! freestream stagnation temperature
    temp_inf      = 1.0_wp
    stag_temp_inf = temp_inf * (1.0_wp + 0.5 * gamm1 * mach_inf**2)

    ! freestream stagnation pressure
    stag_pres_inf = p_inf * (1.0_wp + 0.5 * gamm1 * mach_inf**2)**gamma_over_gamm1

    ! inlet bc's
    ! 1. find riem1 at freestream and calculate riem1 at i = 1
    riem1_inf = mach_inf + 2 * over_gamm1
    riem1_1   = riem1_inf

    do jc = -1, jc_max + 2
        ! 2. find riem2 at i = 2 and calculate riem2 at i = 1
        riem2_2 = vmag(2, jc) - (2 * vsnd(2, jc)) * over_gamm1
        riem2_1 = riem2_2

        ! 3. find velocity at i = 1
        vmag(1, jc) = 0.5 * (riem1_1 + riem2_1)

        ! 4. find u and v at i = 1
        xvel(1, jc) = vmag(1, jc) * cos(inlet_aoa)
        yvel(1, jc) = vmag(1, jc) * sin(inlet_aoa)

        ! 5. find speed of sound at i = 1
        vsnd(1, jc) = 0.25 * gamm1 * (riem1_1 - riem2_1)

        ! 6. find mach number at i = 1
        mach(1, jc) = vmag(1, jc) / vsnd(1, jc)

        ! 7. find static pressure at i = 1
        pres(1, jc) = stag_pres_inf / (1.0_wp + 0.5 * gamm1 * mach(1, jc)**2)**gamma_over_gamm1

        ! 8. find density at i = 1
        dens(1, jc) = (gamma * pres(1, jc)) / vsnd(1, jc)**2

        ! 9. find energy at i = 1
        enrg(1, jc) = pres(1, jc) / (dens(1, jc) * gamm1) + 0.5 * vmag(1, jc)**2

        ! 10. find T and s
        temp(1, jc) = stag_temp_inf / (1.0_wp + 0.5 * gamm1 * mach(1, jc)**2)
        entr(1, jc) = pres(1, jc) / dens(1, jc)**gamma

        ! 10. update state vector components
        q(1, jc, 1) = dens(1, jc)
        q(1, jc, 2) = dens(1, jc) * xvel(1, jc)
        q(1, jc, 3) = dens(1, jc) * yvel(1, jc)
        q(1, jc, 4) = dens(1, jc) * enrg(1, jc)
    end do

    ! outlet bc's
    do jc = -1, jc_max + 2
        ! 1. static pressure at i = i_max

        ! check if the flow inside the last portion of the domain is supersonic
        if (mach(ic_max - 1, jc) >= 1.0_wp) then
            ! supersonic case
            pres(ic_max, jc) = pres(ic_max - 1, jc)
        else
            ! subsonic case
            pres(ic_max, jc) = p_exit
        end if

        ! 2. entropy at i = i_max
        entr(ic_max, jc) = entr(ic_max-1, jc)

        ! 3. riem1 at i = i_max
        riem1_1 = vmag(ic_max-1, jc) + (2 * vsnd(ic_max-1, jc)) * over_gamm1

        ! 4. y-velocity at i = i_max
        yvel(ic_max, jc) = yvel(ic_max-1, jc)

        ! 5. other variables at i = i_max
        dens(ic_max, jc) = (pres(ic_max, jc) / entr(ic_max, jc))**over_gamma
        vsnd(ic_max, jc) = sqrt((gamma * pres(ic_max, jc)) / dens(ic_max, jc))
        vmag(ic_max, jc) = riem1_1 - (2 * vsnd(ic_max, jc)) * over_gamm1
        xvel(ic_max, jc) = sqrt(vmag(ic_max, jc)**2 - yvel(ic_max, jc)**2)
        enrg(ic_max, jc) = pres(ic_max, jc) / (dens(ic_max, jc) * gamm1) + 0.5 * vmag(ic_max, jc)**2
        temp(ic_max, jc) = ((enrg(ic_max, jc) * a_ref**2 - 0.5 * &
                            (vmag(ic_max, jc) * a_ref)**2) / cv_ref) / T_ref
        mach(ic_max, jc) = vmag(ic_max, jc) / vsnd(ic_max, jc)

        ! 6. update state vector components
        q(ic_max, jc, 1) = dens(ic_max, jc)
        q(ic_max, jc, 2) = dens(ic_max, jc) * xvel(ic_max, jc)
        q(ic_max, jc, 3) = dens(ic_max, jc) * yvel(ic_max, jc)
        q(ic_max, jc, 4) = dens(ic_max, jc) * enrg(ic_max, jc)

        ! assert the zero-gradient bc (asserts that there is essentially no change in the state
        ! vector between the exit and the first ghost cell)

        ! TODO: this isn't really true, need to explore ways to increase the accuracy here
        do k = 1, 4
            q(ic_max + 1, jc, k) = q(ic_max, jc, k)
        end do
    end do

end subroutine bcinout