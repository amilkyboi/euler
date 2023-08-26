subroutine area()
    use gridprop,  only: ic, jc, ic_max, jc_max, xn, yn, ar
    use reference, only: l_ref
    use functions
    use timing
    implicit none

    ! Finds the area of each cell, including ghost cells. After the complete array is computed, the
    ! four ghost cells on each corner have their area removed since they are not needed for
    ! computation.

    integer :: an(2), bn(2), cn(2), dn(2)

    call system_clock(start, rate)

    do ic = 1, ic_max
        do jc = 1, jc_max
            an = ijnode(ic, jc, 1)
            bn = ijnode(ic, jc, 2)
            cn = ijnode(ic, jc, 3)
            dn = ijnode(ic, jc, 4)

            ar(ic, jc) = 0.5*((xn(cn(1), cn(2)) - xn(an(1), an(2)))*(yn(dn(1), dn(2)) - yn(bn(1), bn(2))) - &
                              (xn(dn(1), dn(2)) - xn(bn(1), bn(2)))*(yn(cn(1), cn(2)) - yn(an(1), an(2))))
        end do
    end do

    ! non-dimensionalize the entire area vector
    ar = ar / l_ref**2

    call system_clock(end)
    print *, 'subroutine area took ', (end - start) / rate, ' seconds'

end subroutine area