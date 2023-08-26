subroutine bcwall()
    use mod_types, only: wp => dp
    use flow_vars, only: dens, xvel, yvel, vmag, pres, enrg, temp, vsnd, entr, mach
    use reference, only: a_ref, cv_ref, T_ref, R_ref
    use grid_vars, only: ic, ic_max, jc_max
    use gas_vars,  only: gamma, gammam1
    use flux_vars, only: q
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
        dens(ic, 0) = q(ic, 0, 1)
        xvel(ic, 0) = q(ic, 0, 2) / q(ic, 0, 1)
        yvel(ic, 0) = q(ic, 0, 3) / q(ic, 0, 1)
        enrg(ic, 0) = q(ic, 0, 4) / q(ic, 0, 1)

        vmag(ic, 0) = sqrt(xvel(ic, 0)**2 + yvel(ic, 0)**2)
        temp(ic, 0) = ((enrg(ic, 0) * a_ref**2 - 0.5 * (vmag(ic, 0) * a_ref)**2) / cv_ref) / T_ref
        vsnd(ic, 0) = sqrt(gamma * R_ref * temp(ic, 0) * T_ref) / a_ref
        mach(ic, 0) = vmag(ic, 0) / vsnd(ic, 0)
        pres(ic, 0) = gammam1 * dens(ic, 0) * (enrg(ic, 0) - 0.5 * vmag(ic, 0)**2)
        entr(ic, 0) = pres(ic, 0) / dens(ic, 0)**gamma

        ! second layer of ghost cells on the bottom
        ! update state vector
        q(ic, -1, 1) = q(ic, 2, 1)
        q(ic, -1, 2) = q(ic, 2, 2) - 2 * dotn * nn(1)
        q(ic, -1, 3) = q(ic, 2, 3) - 2 * dotn * nn(2)
        q(ic, -1, 4) = q(ic, 2, 4)

        ! update tracked quantities
        dens(ic, -1) = q(ic, -1, 1)
        xvel(ic, -1) = q(ic, -1, 2) / q(ic, -1, 1)
        yvel(ic, -1) = q(ic, -1, 3) / q(ic, -1, 1)
        enrg(ic, -1) = q(ic, -1, 4) / q(ic, -1, 1)

        vmag(ic, -1) = sqrt(xvel(ic, -1)**2 + yvel(ic, -1)**2)
        temp(ic, -1) = ((enrg(ic, -1) * a_ref**2 - 0.5 * (vmag(ic, -1) * a_ref)**2) / cv_ref) / T_ref
        vsnd(ic, -1) = sqrt(gamma * R_ref * temp(ic, -1) * T_ref) / a_ref
        mach(ic, -1) = vmag(ic, -1) / vsnd(ic, -1)
        pres(ic, -1) = gammam1 * dens(ic, -1) * (enrg(ic, -1) - 0.5 * vmag(ic, -1)**2)
        entr(ic, -1) = pres(ic, -1) / dens(ic, -1)**gamma

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
        dens(ic, jc_max+1) = q(ic, jc_max+1, 1)
        xvel(ic, jc_max+1) = q(ic, jc_max+1, 2) / q(ic, jc_max+1, 1)
        yvel(ic, jc_max+1) = q(ic, jc_max+1, 3) / q(ic, jc_max+1, 1)
        enrg(ic, jc_max+1) = q(ic, jc_max+1, 4) / q(ic, jc_max+1, 1)

        vmag(ic, jc_max+1) = sqrt(xvel(ic, jc_max+1)**2 + yvel(ic, jc_max+1)**2)
        temp(ic, jc_max+1) = ((enrg(ic, jc_max+1) * a_ref**2 - 0.5 * (vmag(ic, jc_max+1) * a_ref)**2) / cv_ref) / T_ref
        vsnd(ic, jc_max+1) = sqrt(gamma * R_ref * temp(ic, jc_max+1) * T_ref) / a_ref
        mach(ic, jc_max+1) = vmag(ic, jc_max+1) / vsnd(ic, jc_max+1)
        pres(ic, jc_max+1) = gammam1 * dens(ic, jc_max+1) * (enrg(ic, jc_max+1) - 0.5 * vmag(ic, jc_max+1)**2)
        entr(ic, jc_max+1) = pres(ic, jc_max+1) / dens(ic, jc_max+1)**gamma

        ! second layer of ghost cells on the top
        ! update state vector
        q(ic, jc_max+2, 1) = q(ic, jc_max-1, 1)
        q(ic, jc_max+2, 2) = q(ic, jc_max-1, 2) - 2 * dots * ns(1)
        q(ic, jc_max+2, 3) = q(ic, jc_max-1, 3) - 2 * dots * ns(2)
        q(ic, jc_max+2, 4) = q(ic, jc_max-1, 4)

        ! update tracked quantities
        dens(ic, jc_max+2) = q(ic, jc_max+2, 1)
        xvel(ic, jc_max+2) = q(ic, jc_max+2, 2) / q(ic, jc_max+2, 1)
        yvel(ic, jc_max+2) = q(ic, jc_max+2, 3) / q(ic, jc_max+2, 1)
        enrg(ic, jc_max+2) = q(ic, jc_max+2, 4) / q(ic, jc_max+2, 1)

        vmag(ic, jc_max+2) = sqrt(xvel(ic, jc_max+2)**2 + yvel(ic, jc_max+2)**2)
        temp(ic, jc_max+2) = ((enrg(ic, jc_max+2) * a_ref**2 - 0.5 * (vmag(ic, jc_max+2) * a_ref)**2) / cv_ref) / T_ref
        vsnd(ic, jc_max+2) = sqrt(gamma * R_ref * temp(ic, jc_max+2) * T_ref) / a_ref
        mach(ic, jc_max+2) = vmag(ic, jc_max+2) / vsnd(ic, jc_max+2)
        pres(ic, jc_max+2) = gammam1 * dens(ic, jc_max+2) * (enrg(ic, jc_max+2) - 0.5 * vmag(ic, jc_max+2)**2)
        entr(ic, jc_max+2) = pres(ic, jc_max+2) / dens(ic, jc_max+2)**gamma
    end do
end subroutine bcwall