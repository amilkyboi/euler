subroutine initializer()
    use mod_types,  only: wp => dp
    use gridprop,   only: ic_max, jc_max, ic, jc, xn, yn, alpha
    use flowprop,   only: r, p, T, E, c, s, u, v, vel, mach
    use reference,  only: a_ref, R_ref, T_ref, cv_ref
    use gasprop,    only: gamma, gammam1
    use input,      only: mach_inf
    use fluxes,     only: q, f, g
    use functions
    use timing
    implicit none

    ! Initializes the state vector and the two flux vectors in each cell. The values are calculated
    ! based on conservation properties. The initial mach number for each cell is found based on the
    ! cell's height ratio compared to the inlet cell. This method provides a more accurate first
    ! guess than initializing the state vector and flux vectors to purely freestream values.

    integer :: an(2), bn(2), cn(2), dn(2)
    real(wp) :: cnst1, cnst2, mach_temp, alfa

    call cpu_time(start)

    ! commonly reused constants
    cnst1 = 1.0_wp / (gamma * gammam1)
    cnst2 = 1.0_wp / gamma

    do ic = -1, ic_max + 2
        do jc = -1, jc_max + 2
            ! mach number for each cell based on height ratio
            mach_temp = mach_inf! * (yn(1, jc_max) - yn(1, 1)) / &
                                !   (0.5 * ((yn(ic, jc_max) - yn(ic, 1)) + (yn(ic+1, jc_max) - yn(ic+1, 1))))

            an = ijnode(ic, jc, 1)
            bn = ijnode(ic, jc, 2)
            cn = ijnode(ic, jc, 3)
            dn = ijnode(ic, jc, 4)

            ! angle of attack of each cell, calculated as an average between top and bottom slopes
            alpha(ic, jc) = atan(0.5 * ((yn(dn(1), dn(2)) - yn(cn(1), cn(2)))/(xn(dn(1), dn(2)) - xn(cn(1), cn(2))) + & 
                                        (yn(bn(1), bn(2)) - yn(an(1), an(2)))/(xn(bn(1), bn(2)) - xn(an(1), an(2)))))
            alfa = alpha(ic, jc)

            ! initialize state vector
            q(ic, jc, 1) = 1.0_wp
            q(ic, jc, 2) = mach_temp * cos(alfa)
            q(ic, jc, 3) = mach_temp * sin(alfa)
            q(ic, jc, 4) = cnst1 + 0.5 * mach_temp**2

            ! store all relevant flow properties seprately in non-dimensional form
            ! rho, u, v, vel, E, T, c, mach, p, s
            r(ic, jc) = q(ic, jc, 1)
            u(ic, jc) = q(ic, jc, 2) / q(ic, jc, 1)
            v(ic, jc) = q(ic, jc, 3) / q(ic, jc, 1)
            E(ic, jc) = q(ic, jc, 4) / q(ic, jc, 1)

            vel(ic, jc) = sqrt(u(ic, jc)**2 + v(ic, jc)**2)
            T(ic, jc) = ((E(ic, jc) * a_ref**2 - 0.5 * (vel(ic, jc) * a_ref)**2) / cv_ref) / T_ref
            c(ic, jc) = sqrt(gamma * R_ref * T(ic, jc) * T_ref) / a_ref
            mach(ic, jc) = vel(ic, jc) / c(ic, jc)
            p(ic, jc) = gammam1 * r(ic, jc) * (E(ic, jc) - 0.5 * vel(ic, jc)**2)
            s(ic, jc) = p(ic, jc) / r(ic, jc)**gamma

            ! initialize fluxes
            f(ic, jc, 1) = mach_temp * cos(alfa)
            f(ic, jc, 2) = (mach_temp * cos(alfa))**2 + cnst2
            f(ic, jc, 3) = mach_temp**2 * cos(alfa) * sin(alfa)
            f(ic, jc, 4) = (cnst1 + 0.5 * mach_temp**2 + cnst2)*mach_temp*cos(alfa)

            g(ic, jc, 1) = mach_temp * sin(alfa)
            g(ic, jc, 2) = mach_temp**2 * sin(alfa) * cos(alfa)
            g(ic, jc, 3) = (mach_temp * sin(alfa))**2 + cnst2
            g(ic, jc, 4) = (cnst1 + 0.5 * mach_temp**2 + cnst2)*mach_temp*sin(alfa)
        end do
    end do

    call cpu_time(end)
    print *, 'subroutine initializer took ', end - start, ' seconds'

end subroutine initializer