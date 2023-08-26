subroutine timestep()
    use mod_types, only: wp => dp
    use gridprop,  only: ic, jc, ic_max, jc_max, ar, dt_min, dt
    use flowvars,  only: vsnd, xvel, yvel
    use input,     only: cfl    
    use functions
    use timing
    implicit none

    ! Calculate the time step needed for the RK4 iterations. Only considers the real cells.

    real(wp) :: nn(2), ns(2), ne(2), nw(2), eign, eigs, eige, eigw

    call system_clock(start, rate)

    do ic = 1, ic_max
        do jc = 1, jc_max
            nn = normal(ic, jc, 1)
            ns = normal(ic, jc, 2)
            ne = normal(ic, jc, 3)
            nw = normal(ic, jc, 4)

            eign = abs(xvel(ic, jc)*nn(1) + yvel(ic, jc)*nn(2)) + vsnd(ic, jc) ! north
            eigs = abs(xvel(ic, jc)*ns(1) + yvel(ic, jc)*ns(2)) + vsnd(ic, jc) ! south
            eige = abs(xvel(ic, jc)*ne(1) + yvel(ic, jc)*ne(2)) + vsnd(ic, jc) ! east
            eigw = abs(xvel(ic, jc)*nw(1) + yvel(ic, jc)*nw(2)) + vsnd(ic, jc) ! west

            dt(ic, jc) = (2 * ar(ic, jc)) / (eign*length(ic, jc, 1) + &
                          eigs*length(ic, jc, 2) + eige*length(ic, jc, 3) + eigw*length(ic, jc, 4))
        end do
    end do

    dt_min = cfl * minval(dt)

    call system_clock(end)
    print *, 'subroutine timestep took ', (end - start) / rate, ' seconds'

end subroutine timestep