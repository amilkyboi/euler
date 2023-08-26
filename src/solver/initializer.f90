subroutine initializer()
    use mod_types, only: wp => dp
    use flow_vars, only: dens, pres, temp, enrg, vsnd, entr, xvel, yvel, vmag, mach
    use gas_vars,  only: gamma, gammam1, over_gamma, over_gtgm1
    use grid_vars, only: ic_max, jc_max, ic, jc, xn, yn
    use reference, only: a_ref, R_ref, T_ref, cv_ref
    use input,     only: mach_inf
    use flux_vars, only: q, f, g
    use functions
    use timing
    implicit none

    ! initializes the state vector and the two flux vectors in each cell

    ! used to store node locations (in, jn)
    integer :: node_a(2), node_b(2), node_c(2), node_d(2)
    ! various temporary variables for calculating continuity throughout the domain
    real(wp) :: mach_continuity, cell_aoa, inlet_height, current_avg_domain_height

    call system_clock(start, rate)

    ! calculates the height of the inlet
    inlet_height = yn(1, jc_max) - yn(1, 1)

    do ic = -1, ic_max + 2
        ! iterates along the domain in the x-direction and gets the current average height of the
        ! domain
        current_avg_domain_height = 0.5 * ((yn(ic, jc_max) - yn(ic, 1)) + &
                                           (yn(ic+1, jc_max) - yn(ic+1, 1)))

        ! applies continuity by accelerating the flow when the domain converges

        ! TODO: this assumption only works for a subsonic inlet condition, should make a switch that
        !       turns it off when mach_inf >= 1.0_wp
        mach_continuity = mach_inf * inlet_height / current_avg_domain_height

        do jc = -1, jc_max + 2
            ! convert from cell notation (ic, jc) to node notation (in, jn)
            node_a = ijcell_to_ijnode(ic, jc, 1)
            node_b = ijcell_to_ijnode(ic, jc, 2)
            node_c = ijcell_to_ijnode(ic, jc, 3)
            node_d = ijcell_to_ijnode(ic, jc, 4)

            ! angle of attack of each cell, calculated as an average between the slope of the top
            ! wall and the slope of the bottom wall
            cell_aoa = atan(0.5 * ((yn(node_d(1), node_d(2)) - yn(node_c(1), node_c(2))) / &
                                   (xn(node_d(1), node_d(2)) - xn(node_c(1), node_c(2))) + & 
                                   (yn(node_b(1), node_b(2)) - yn(node_a(1), node_a(2))) / &
                                   (xn(node_b(1), node_b(2)) - xn(node_a(1), node_a(2)))))

            ! initialize the state vector with non-dimensional freestream values
            q(ic, jc, 1) = 1.0_wp
            q(ic, jc, 2) = mach_continuity * cos(cell_aoa)
            q(ic, jc, 3) = mach_continuity * sin(cell_aoa)
            q(ic, jc, 4) = over_gtgm1 + 0.5 * mach_continuity**2

            ! store all relevant flow properties seprately in non-dimensional form
            ! components directly obtainable from the state vector
            dens(ic, jc) = q(ic, jc, 1)
            xvel(ic, jc) = q(ic, jc, 2) / q(ic, jc, 1)
            yvel(ic, jc) = q(ic, jc, 3) / q(ic, jc, 1)
            enrg(ic, jc) = q(ic, jc, 4) / q(ic, jc, 1)

            ! other useful quantities
            vmag(ic, jc) = sqrt(xvel(ic, jc)**2 + yvel(ic, jc)**2)
            temp(ic, jc) = ((enrg(ic, jc) * a_ref**2 - 0.5 * (vmag(ic, jc) * a_ref)**2) / cv_ref) / T_ref
            vsnd(ic, jc) = sqrt(gamma * R_ref * temp(ic, jc) * T_ref) / a_ref
            mach(ic, jc) = vmag(ic, jc) / vsnd(ic, jc)
            pres(ic, jc) = gammam1 * dens(ic, jc) * (enrg(ic, jc) - 0.5 * vmag(ic, jc)**2)
            entr(ic, jc) = pres(ic, jc) / dens(ic, jc)**gamma

            ! f flux
            f(ic, jc, 1) =  mach_continuity * cos(cell_aoa)
            f(ic, jc, 2) = (mach_continuity * cos(cell_aoa))**2 + over_gamma
            f(ic, jc, 3) =  mach_continuity**2 * cos(cell_aoa) * sin(cell_aoa)
            f(ic, jc, 4) = (over_gtgm1 + 0.5 * mach_continuity**2 + over_gamma) * &
                            mach_continuity * cos(cell_aoa)

            ! g flux
            g(ic, jc, 1) =  mach_continuity * sin(cell_aoa)
            g(ic, jc, 2) =  mach_continuity**2 * sin(cell_aoa) * cos(cell_aoa)
            g(ic, jc, 3) = (mach_continuity * sin(cell_aoa))**2 + over_gamma
            g(ic, jc, 4) = (over_gtgm1 + 0.5 * mach_continuity**2 + over_gamma) * &
                            mach_continuity * sin(cell_aoa)
        end do
    end do

    call system_clock(end)
    print *, 'subroutine initializer took ', (end - start) / rate, ' seconds'

end subroutine initializer