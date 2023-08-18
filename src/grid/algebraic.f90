subroutine algebraic()
    use gridprop, only: fn_beg, fn_end, yn_min, yn_max, fn_bounds, xy_bounds, in, jn, &
                        in_max, jn_max, xn, yn, xi, et, start, end
    use functions
    implicit none

    integer :: xn_min, xn_max

    call cpu_time(start)

    fn_beg = fn_bounds(1)
    fn_end = fn_bounds(2)
    xn_min = xy_bounds(1)
    xn_max = xy_bounds(2)
    yn_min = xy_bounds(3)
    yn_max = xy_bounds(4)

    do in = 1, in_max
        do jn = 1, jn_max
            ! Map from the (in, jn) grid to the (xi, et) grid
            xi(in, jn) = (in - 1)/dble(in_max - 1)
            et(in, jn) = (jn - 1)/dble(jn_max - 1)
            
            ! Map from the (xi, et) grid to the (xn, yn) grid
            xn(in, jn) = xn_min + xi(in, jn)*(xn_max - xn_min)

            ! Determine where the functions defining the top and bottom boundaries apply
            if ((xn(in, jn) <= fn_beg) .or. (xn(in, jn) >= fn_end)) then
                yn(in, jn) = yn_min + et(in, jn)*(yn_max - yn_min)
            else
                yn(in, jn) = y_btm(xn(in, jn)) + et(in, jn)*(y_top(xn(in, jn)) - y_btm(xn(in, jn)))
            end if
        end do
    end do

    call cpu_time(end)
    print *, 'subroutine algebraic took ', end - start, ' seconds'

end subroutine algebraic