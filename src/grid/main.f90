program main
    use mod_types, only: wp => dp
    use gridprop
    implicit none

    n_max  = 100000
    in_max = 51
    jn_max = 11
    tol    = 10.0**(-10.0)

    allocate(xn(in_max, jn_max), yn(in_max, jn_max), xi(in_max, jn_max), et(in_max, jn_max))

    xy_bounds = (/0, 5, 0, 1/)
    fn_bounds = (/2, 3/)

    call algebraic
    
    call elliptic

    print *, 'done'

end program main