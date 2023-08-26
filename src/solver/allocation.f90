subroutine allocation()
    use grid_vars, only: in_max, jn_max, ic_max, jc_max, area, dt
    use flow_vars, only: dens, xvel, yvel, vmag, enrg, temp, vsnd, mach, pres, entr
    use flux_vars,   only: q, f, g, res, dis
    use timing
    implicit none

    call system_clock(start, rate)

    ! Allocates all arrays except for xn and yn.

    ! define the number of cells in the real domain (not including ghost cells)
    ic_max = in_max - 1
    jc_max = jn_max - 1

    ! allocate grid properties
    allocate(area(-1:ic_max + 2, -1:jc_max + 2), dt(ic_max, jc_max))
    ! allocate state vector and fluxes
    allocate(q(-1:ic_max+2, -1:jc_max+2, 4), f(-1:ic_max+2, -1:jc_max+2, 4), &
             g(-1:ic_max+2, -1:jc_max+2, 4), dis(ic_max, jc_max, 4), &
             res(ic_max, jc_max, 4))
    ! allocate flow properties
    allocate(dens(-1:ic_max+2, -1:jc_max+2), xvel(-1:ic_max+2, -1:jc_max+2), yvel(-1:ic_max+2, -1:jc_max+2), &
             vmag(-1:ic_max+2, -1:jc_max+2), enrg(-1:ic_max+2, -1:jc_max+2), temp(-1:ic_max+2, -1:jc_max+2), &
             vsnd(-1:ic_max+2, -1:jc_max+2), mach(-1:ic_max+2, -1:jc_max+2), pres(-1:ic_max+2, -1:jc_max+2), &
             entr(-1:ic_max+2, -1:jc_max+2))

    call system_clock(end)
    print *, 'subroutine allocation took ', (end - start) / rate, ' seconds'

end subroutine allocation