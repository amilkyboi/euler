subroutine flux()
    use mod_types, only: wp => dp
    use flowprop, only: r, u, v, E, vel, p, T, c, s, mach
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
            r(ic, jc) = q(ic, jc, 1)
            u(ic, jc) = q(ic, jc, 2) / q(ic, jc, 1)
            v(ic, jc) = q(ic, jc, 3) / q(ic, jc, 1)
            E(ic, jc) = q(ic, jc, 4) / q(ic, jc, 1)
            vel(ic, jc) = sqrt(u(ic, jc)**2 + v(ic, jc)**2)
            p(ic, jc) = gammam1 * r(ic, jc) * (E(ic, jc) - 0.5 * vel(ic, jc)**2)
            T(ic, jc) = ((E(ic, jc) * a_ref**2 - 0.5 * (vel(ic, jc) * a_ref)**2) / cv_ref) / T_ref
            c(ic, jc) = sqrt(gamma * R_ref * T(ic, jc) * T_ref) / a_ref
            s(ic, jc) = p(ic, jc) / r(ic, jc)**gamma
            mach(ic, jc) = vel(ic, jc) / c(ic, jc)

            ! update fluxes
            f(ic, jc, 1) =  q(ic, jc, 2)
            f(ic, jc, 2) = (q(ic, jc, 2)**2 / q(ic, jc, 1)) + p(ic, jc)
            f(ic, jc, 3) = (q(ic, jc, 2) * q(ic, jc, 3)) / q(ic, jc, 1)
            f(ic, jc, 4) = (q(ic, jc, 4) + p(ic, jc)) * (q(ic, jc, 2) / q(ic, jc, 1))

            g(ic, jc, 1) =  q(ic, jc, 3)
            g(ic, jc, 2) = (q(ic, jc, 3) * q(ic, jc, 2)) / q(ic, jc, 1)
            g(ic, jc, 3) = (q(ic, jc, 3)**2 / q(ic, jc, 1)) + p(ic, jc)
            g(ic, jc, 4) = (q(ic, jc, 4) + p(ic, jc)) * (q(ic, jc, 3) / q(ic, jc, 1))
        end do
    end do

end subroutine