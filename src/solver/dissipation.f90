subroutine dissipation()
    use mod_types, only: wp => dp
    use gridprop,  only: ic, jc, ic_max, jc_max
    use flowprop,  only: p, u, v, c
    use input,     only: nu2, nu4
    use fluxes,    only: q, dis
    use functions
    implicit none

    ! Calculates the dissipation term, the most critical part of the Jameson-Schmidt-Turkel scheme.

    real(wp), allocatable :: s2c_xi(:, :), s2c_et(:, :), s2f(:, :, :), s4f(:, :, :), eigf(:, :, :)
    real(wp) :: nn(2), ns(2), ne(2), nw(2), delta2p_xi, delta2p_et

    allocate(s2c_xi(0:ic_max+1, 0:jc_max+1), s2c_et(0:ic_max+1, 0:jc_max+1), s2f(ic_max, jc_max, 4), &
             s4f(ic_max, jc_max, 4), eigf(ic_max, jc_max, 4))

    ! second-order switches on cells ij
    do ic = 0, ic_max + 1
        do jc = 0, jc_max + 1
            delta2p_xi = p(ic+1, jc) - 2*p(ic, jc) + p(ic-1, jc)
            delta2p_et = p(ic, jc+1) - 2*p(ic, jc) + p(ic, jc-1)

            s2c_xi(ic, jc) = (nu2*abs(delta2p_xi))/(p(ic+1, jc) + 2*p(ic, jc) + p(ic-1, jc))
            s2c_et(ic, jc) = (nu2*abs(delta2p_et))/(p(ic+1, jc) + 2*p(ic, jc) + p(ic-1, jc))
        end do
    end do

    ! second-order and fourth-order switches on faces
    do ic = 1, ic_max
        do jc = 1, jc_max
            s2f(ic, jc, 1) = 0.5 * (s2c_et(ic, jc+1) + s2c_et(ic, jc)) ! north
            s2f(ic, jc, 2) = 0.5 * (s2c_et(ic, jc-1) + s2c_et(ic, jc)) ! south
            s2f(ic, jc, 3) = 0.5 * (s2c_xi(ic+1, jc) + s2c_xi(ic, jc)) ! east
            s2f(ic, jc, 4) = 0.5 * (s2c_xi(ic-1, jc) + s2c_xi(ic, jc)) ! west

            s4f(ic, jc, 1) = max(0.0_wp, nu4 - s2f(ic, jc, 1)) ! north
            s4f(ic, jc, 2) = max(0.0_wp, nu4 - s2f(ic, jc, 2)) ! south
            s4f(ic, jc, 3) = max(0.0_wp, nu4 - s2f(ic, jc, 3)) ! east
            s4f(ic, jc, 4) = max(0.0_wp, nu4 - s2f(ic, jc, 4)) ! west
        end do
    end do

    ! eigenvalues on faces
    do ic = 1, ic_max
        do jc = 1, jc_max
            ! normal vectors for each face
            nn = normal(ic, jc, 1) ! north
            ns = normal(ic, jc, 2) ! south
            ne = normal(ic, jc, 3) ! east
            nw = normal(ic, jc, 4) ! west

            eigf(ic, jc, 1) = abs(u(ic, jc)*nn(1) + v(ic, jc)*nn(2)) + c(ic, jc) ! north
            eigf(ic, jc, 2) = abs(u(ic, jc)*ns(1) + v(ic, jc)*ns(2)) + c(ic, jc) ! south
            eigf(ic, jc, 3) = abs(u(ic, jc)*ne(1) + v(ic, jc)*ne(2)) + c(ic, jc) ! east
            eigf(ic, jc, 4) = abs(u(ic, jc)*nw(1) + v(ic, jc)*nw(2)) + c(ic, jc) ! west
        end do
    end do

    ! dissipation
    do ic = 1, ic_max
        do jc = 1, jc_max
            dis(ic, jc, :) = ((s2f(ic, jc, 3)*length(ic, jc, 3)*eigf(ic, jc, 3)*(q(ic+1, jc, :) - q(ic, jc, :)) - &
                               s2f(ic, jc, 4)*length(ic, jc, 4)*eigf(ic, jc, 4)*(q(ic, jc, :) - q(ic-1, jc, :))) + &
                              (s2f(ic, jc, 1)*length(ic, jc, 1)*eigf(ic, jc, 1)*(q(ic, jc+1, :) - q(ic, jc, :)) - &
                               s2f(ic, jc, 2)*length(ic, jc, 2)*eigf(ic, jc, 2)*(q(ic, jc, :) - q(ic, jc-1, :)))) - &
                             ((s4f(ic, jc, 3)*length(ic, jc, 3)*eigf(ic, jc, 3)*((q(ic+2, jc, :) - &
                               2*q(ic+1, jc, :) + q(ic, jc, :)) - (q(ic+1, jc, :) - 2*q(ic, jc, :) + q(ic-1, jc, :))) - &
                               s4f(ic, jc, 4)*length(ic, jc, 4)*eigf(ic, jc, 4)*((q(ic+1, jc, :) - &
                               2*q(ic, jc, :) + q(ic-1, jc, :)) - (q(ic, jc, :) - 2*q(ic-1, jc, :) + q(ic-2, jc, :)))) + &
                              (s4f(ic, jc, 1)*length(ic, jc, 1)*eigf(ic, jc, 1)*((q(ic, jc+2, :) - &
                               2*q(ic, jc+1, :) + q(ic, jc, :)) - (q(ic, jc+1, :) - 2*q(ic, jc, :) + q(ic, jc-1, :))) - &
                               s4f(ic, jc, 2)*length(ic, jc, 2)*eigf(ic, jc, 2)*((q(ic, jc+1, :) - &
                               2*q(ic, jc, :) + q(ic, jc-1, :)) - (q(ic, jc, :) - 2*q(ic, jc-1, :) + q(ic, jc-2, :)))))
        end do
    end do

end subroutine dissipation