subroutine timestep()
    use mod_types, only: wp => dp
    use gridprop,  only: ic, jc, ic_max, jc_max, ar, dt_min, dt
    use reference, only: a_ref, l_ref
    use flowprop,  only: c, u, v
    use input,     only: cfl    
    use functions
    implicit none

    ! Calculate the time step needed for the RK4 iterations. Only considers the real cells.

    real(wp) :: nn(2), ns(2), ne(2), nw(2), eign, eigs, eige, eigw

    do ic = 1, ic_max
        do jc = 1, jc_max
            nn = normal(ic, jc, 1)
            ns = normal(ic, jc, 2)
            ne = normal(ic, jc, 3)
            nw = normal(ic, jc, 4)

            eign = abs(u(ic, jc)*nn(1) + v(ic, jc)*nn(2)) + c(ic, jc) ! north
            eigs = abs(u(ic, jc)*ns(1) + v(ic, jc)*ns(2)) + c(ic, jc) ! south
            eige = abs(u(ic, jc)*ne(1) + v(ic, jc)*ne(2)) + c(ic, jc) ! east
            eigw = abs(u(ic, jc)*nw(1) + v(ic, jc)*nw(2)) + c(ic, jc) ! west

            dt(ic, jc) = (2 * ar(ic, jc)) / (eign*length(ic, jc, 1) + &
                          eigs*length(ic, jc, 2) + eige*length(ic, jc, 3) + eigw*length(ic, jc, 4))
        end do
    end do

    dt_min = cfl * minval(dt)

end subroutine timestep