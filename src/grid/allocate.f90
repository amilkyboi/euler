subroutine allocate()
    use gridprop, only: xn, yn, xi, et
    use input,    only: in_max, jn_max
    implicit none

    allocate(xn(in_max, jn_max), yn(in_max, jn_max), xi(in_max, jn_max), et(in_max, jn_max))

end subroutine allocate