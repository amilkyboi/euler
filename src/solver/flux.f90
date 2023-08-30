subroutine flux()
    use flow_vars, only: dens, xvel, yvel, enrg, vmag, pres, temp, vsnd, entr, mach
    use reference, only: R_ref, T_ref, a_ref, cv_ref
    use grid_vars, only: ic, jc, ic_max, jc_max
    use gas_vars,  only: gamma, gamm1
    use flux_vars, only: q, f, g
    implicit none

    ! calculates the f and g fluxes based on the state vector values

    ! TODO: see about accessing more variables directly from the state vector instead of going
    !       through the flow property arrays

    do ic = -1, ic_max + 2
        do jc = -1, jc_max + 2
            ! first update flow properties directly with the current state vector
            dens(ic, jc) = q(ic, jc, 1)
            xvel(ic, jc) = q(ic, jc, 2) / q(ic, jc, 1)
            yvel(ic, jc) = q(ic, jc, 3) / q(ic, jc, 1)
            enrg(ic, jc) = q(ic, jc, 4) / q(ic, jc, 1)

            ! other properties follow
            vmag(ic, jc) = sqrt(xvel(ic, jc)**2 + yvel(ic, jc)**2)
            pres(ic, jc) = gamm1 * dens(ic, jc) * (enrg(ic, jc) - 0.5 * vmag(ic, jc)**2)
            temp(ic, jc) = ((enrg(ic, jc) * a_ref**2 - 0.5 * (vmag(ic, jc) * a_ref)**2) / cv_ref) / T_ref
            vsnd(ic, jc) = sqrt(gamma * R_ref * temp(ic, jc) * T_ref) / a_ref
            entr(ic, jc) = pres(ic, jc) / dens(ic, jc)**gamma
            mach(ic, jc) = vmag(ic, jc) / vsnd(ic, jc)

            ! update f flux
            f(ic, jc, 1) =  q(ic, jc, 2)
            f(ic, jc, 2) = (q(ic, jc, 2)**2 / q(ic, jc, 1)) + pres(ic, jc)
            f(ic, jc, 3) = (q(ic, jc, 2) * q(ic, jc, 3)) / q(ic, jc, 1)
            f(ic, jc, 4) = (q(ic, jc, 4) + pres(ic, jc)) * (q(ic, jc, 2) / q(ic, jc, 1))

            ! update g flux
            g(ic, jc, 1) =  q(ic, jc, 3)
            g(ic, jc, 2) = (q(ic, jc, 3) * q(ic, jc, 2)) / q(ic, jc, 1)
            g(ic, jc, 3) = (q(ic, jc, 3)**2 / q(ic, jc, 1)) + pres(ic, jc)
            g(ic, jc, 4) = (q(ic, jc, 4) + pres(ic, jc)) * (q(ic, jc, 3) / q(ic, jc, 1))
        end do
    end do

end subroutine