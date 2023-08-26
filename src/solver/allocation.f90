subroutine allocation()
    use gridprop, only: in_max, jn_max, ic_max, jc_max, ar, alpha, dt
    use flowprop, only: r, u, v, vel, E, T, c, mach, p, s
    use fluxes,   only: q, f, g, res, dis
    use timing
    implicit none

    call system_clock(start, rate)

    ! Allocates all arrays except for xn and yn.

    ! define the number of cells in the real domain (not including ghost cells)
    ic_max = in_max - 1
    jc_max = jn_max - 1

    ! allocate grid properties
    allocate(ar(ic_max, jc_max), alpha(-1:ic_max+2, -1:jc_max+2), dt(ic_max, jc_max))
    ! allocate state vector and fluxes
    allocate(q(-1:ic_max+2, -1:jc_max+2, 4), f(-1:ic_max+2, -1:jc_max+2, 4), &
             g(-1:ic_max+2, -1:jc_max+2, 4), dis(ic_max, jc_max, 4), &
             res(ic_max, jc_max, 4))
    ! allocate flow properties
    allocate(r(-1:ic_max+2, -1:jc_max+2), u(-1:ic_max+2, -1:jc_max+2), v(-1:ic_max+2, -1:jc_max+2), &
             vel(-1:ic_max+2, -1:jc_max+2), E(-1:ic_max+2, -1:jc_max+2), T(-1:ic_max+2, -1:jc_max+2), &
             c(-1:ic_max+2, -1:jc_max+2), mach(-1:ic_max+2, -1:jc_max+2), p(-1:ic_max+2, -1:jc_max+2), &
             s(-1:ic_max+2, -1:jc_max+2))

    call system_clock(end)
    print *, 'subroutine allocation took ', (end - start) / rate, ' seconds'

end subroutine allocation