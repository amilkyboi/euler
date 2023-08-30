subroutine allocation()
    use flow_vars, only: dens, xvel, yvel, vmag, enrg, temp, vsnd, mach, pres, entr
    use grid_vars, only: in_max, jn_max, ic_max, jc_max, area, dt
    use flux_vars, only: q, f, g, res, dis
    use timing
    implicit none

    call system_clock(start, rate)

    ! allocates all arrays except for xn and yn

    ! define the number of cells in the real domain (not including ghost cells)
    ic_max = in_max - 1
    jc_max = jn_max - 1

    ! cell areas
    allocate(area(-1:ic_max + 2, -1:jc_max + 2))
    ! time step, dependent on cell geometry
    allocate(dt(ic_max, jc_max))
    ! state vector (rho, rho*u, rho*v, rho*E)
    allocate(q(-1:ic_max + 2, -1:jc_max + 2, 4))
    ! f flux vector (rho*u, rho*u**2 + p, rho*u*v, rho*H*u)
    allocate(f(-1:ic_max + 2, -1:jc_max + 2, 4))
    ! g flux vector (rho*v, rho*u*v, rho*v**2 + p, rho*H*v)
    allocate(g(-1:ic_max + 2, -1:jc_max + 2, 4))
    ! dissipation
    allocate(dis(ic_max, jc_max, 4))
    ! residual
    allocate(res(ic_max, jc_max, 4))

    ! density
    allocate(dens(-1:ic_max + 2, -1:jc_max + 2))
    ! x-velocity
    allocate(xvel(-1:ic_max + 2, -1:jc_max + 2))
    ! y-velocity
    allocate(yvel(-1:ic_max + 2, -1:jc_max + 2))
    ! velocity magnitude (sqrt(u**2 + v**2))
    allocate(vmag(-1:ic_max + 2, -1:jc_max + 2))
    ! energy
    allocate(enrg(-1:ic_max + 2, -1:jc_max + 2))
    ! static temperature
    allocate(temp(-1:ic_max + 2, -1:jc_max + 2))
    ! speed of sound
    allocate(vsnd(-1:ic_max + 2, -1:jc_max + 2))
    ! mach number
    allocate(mach(-1:ic_max + 2, -1:jc_max + 2))
    ! static pressure
    allocate(pres(-1:ic_max + 2, -1:jc_max + 2))
    ! specific entropy
    allocate(entr(-1:ic_max + 2, -1:jc_max + 2))

    call system_clock(end)
    print *, 'subroutine allocation took ', (end - start) / rate, ' seconds'

end subroutine allocation