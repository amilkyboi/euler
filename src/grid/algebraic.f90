subroutine algebraic()
    use mod_types, only: wp => dp
    use gridprop,  only: xi, et, xn, yn, in, jn, yn_min, yn_max, fn_beg, fn_end
    use input,     only: in_max, jn_max, x_bounds, y_bounds, fn_bounds, save_alg, alg_path
    use functions
    use timing
    implicit none

    ! minimum and maximum x-axis bounds of the domain
    integer(wp) :: xn_min, xn_max

    call system_clock(start, rate)

    ! variables for clarity
    fn_beg = fn_bounds(1)
    fn_end = fn_bounds(2)

    xn_min = x_bounds(1)
    xn_max = x_bounds(2)

    yn_min = y_bounds(1)
    yn_max = y_bounds(2)

    do in = 1, in_max
        do jn = 1, jn_max
            ! map from the (in, jn) grid to the (xi, et) grid

            ! ensure the division is carried out with double-precision real
            xi(in, jn) = (in - 1) / dble(in_max - 1)
            et(in, jn) = (jn - 1) / dble(jn_max - 1)
            
            ! map from the (xi, et) grid to the (xn, yn) grid
            xn(in, jn) = xn_min + xi(in, jn)*(xn_max - xn_min)

            ! determine where the functions defining the top and bottom boundaries apply
            if ((xn(in, jn) < fn_beg) .or. (xn(in, jn) > fn_end)) then
                yn(in, jn) = yn_min + et(in, jn) * (yn_max - yn_min)
            else
                yn(in, jn) = y_btm(xn(in, jn)) + &
                             et(in, jn) * (y_top(xn(in, jn)) - y_btm(xn(in, jn)))
            end if
        end do
    end do

    ! optionally save the algebraic mesh
    if (save_alg) then
        open (8, file=alg_path)
        write(8,*) in_max, jn_max
        write(8,*) ((xn(in,jn), in=1,in_max), jn=1,jn_max), ((yn(in,jn), in=1,in_max), jn=1,jn_max)
        close(8)
    end if

    call system_clock(end)
    print *, 'subroutine algebraic took ', (end - start) / rate, ' seconds'

end subroutine algebraic