subroutine flux()
    use flowvars, only: dens, xvel, yvel, enrg, vmag, pres, temp, vsnd, entr, mach
    use reference, only: R_ref, T_ref, a_ref, cv_ref
    use gridprop,  only: ic, jc, ic_max, jc_max
    use gasprop,   only: gamma, gammam1
    use fluxes,    only: q, f, g
    implicit none

    ! Calculates the f and g fluxes based on the state vector values. This is called after each step
    ! of the Runge-Kutta method to update the solution.

    do ic = -1, ic_max + 2
        do jc = -1, jc_max + 2
            ! update flow properties
            dens(ic, jc) = q(ic, jc, 1)
            xvel(ic, jc) = q(ic, jc, 2) / q(ic, jc, 1)
            yvel(ic, jc) = q(ic, jc, 3) / q(ic, jc, 1)
            enrg(ic, jc) = q(ic, jc, 4) / q(ic, jc, 1)
            vmag(ic, jc) = sqrt(xvel(ic, jc)**2 + yvel(ic, jc)**2)
            pres(ic, jc) = gammam1 * dens(ic, jc) * (enrg(ic, jc) - 0.5 * vmag(ic, jc)**2)
            temp(ic, jc) = ((enrg(ic, jc) * a_ref**2 - 0.5 * (vmag(ic, jc) * a_ref)**2) / cv_ref) / T_ref
            vsnd(ic, jc) = sqrt(gamma * R_ref * temp(ic, jc) * T_ref) / a_ref
            entr(ic, jc) = pres(ic, jc) / dens(ic, jc)**gamma
            mach(ic, jc) = vmag(ic, jc) / vsnd(ic, jc)

            ! update fluxes
            f(ic, jc, 1) =  q(ic, jc, 2)
            f(ic, jc, 2) = (q(ic, jc, 2)**2 / q(ic, jc, 1)) + pres(ic, jc)
            f(ic, jc, 3) = (q(ic, jc, 2) * q(ic, jc, 3)) / q(ic, jc, 1)
            f(ic, jc, 4) = (q(ic, jc, 4) + pres(ic, jc)) * (q(ic, jc, 2) / q(ic, jc, 1))

            g(ic, jc, 1) =  q(ic, jc, 3)
            g(ic, jc, 2) = (q(ic, jc, 3) * q(ic, jc, 2)) / q(ic, jc, 1)
            g(ic, jc, 3) = (q(ic, jc, 3)**2 / q(ic, jc, 1)) + pres(ic, jc)
            g(ic, jc, 4) = (q(ic, jc, 4) + pres(ic, jc)) * (q(ic, jc, 3) / q(ic, jc, 1))
        end do
    end do

end subroutine