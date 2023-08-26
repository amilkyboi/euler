subroutine startup()
    use grid_vars, only: in_max, jn_max, in, jn, xn, yn
    use timing
    implicit none

    ! Reads the data created by the elliptic grid generator, allocates the xn and yn arrays, and
    ! stores grid data inside each vector.

    call system_clock(start, rate)

    open(8, file='../../data/xy_elp.x')
    read(8, *) in_max, jn_max

    allocate(xn(-1:in_max+2, -1:jn_max+2), yn(-1:in_max+2, -1:jn_max+2))

    read(8, *) ((xn(in, jn), in=1, in_max), jn=1, jn_max), ((yn(in, jn), in=1, in_max), jn=1, jn_max)
    close(8)
    
    print *, in_max, 'x', jn_max, ' grid'

    call system_clock(end)
    print *, 'subroutine startup took ', (end - start) / rate, ' seconds'

end subroutine startup