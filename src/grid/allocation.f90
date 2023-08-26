subroutine allocation()
    use gridprop, only: xi, et, xn, yn
    use input,    only: in_max, jn_max
    implicit none

    ! allocate the (xi, et) grid
    allocate(xi(in_max, jn_max))
    allocate(et(in_max, jn_max))

    ! allocate the (x, y) grid
    allocate(xn(in_max, jn_max))
    allocate(yn(in_max, jn_max))

end subroutine allocation