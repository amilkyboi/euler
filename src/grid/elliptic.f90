subroutine elliptic()
    use mod_types, only: wp => dp
    use gridprop, only: fn_bounds, xy_bounds, yn_min, yn_max, fn_beg, fn_end, xn, yn, xi, et, &
                        in, jn, in_max, jn_max, n_max, tol, start, end
    use functions
    implicit none

    integer :: min_l, min_r
    integer :: n

    real(wp), allocatable :: dist_l(:), dist_r(:), xn_new(:, :), yn_new(:, :)
    real(wp) :: alp, bet, gam, d_xi, d_et, e_x, e_y

    call cpu_time(start)

    allocate(dist_l(in_max), dist_r(in_max), xn_new(in_max, jn_max), yn_new(in_max, jn_max))

    fn_beg = fn_bounds(1)
    fn_end = fn_bounds(2)
    yn_min = xy_bounds(3)
    yn_max = xy_bounds(4)

    ! Determine the step size in the xi and et directions
    d_xi = xi(2, 1) - xi(1, 1)
    d_et = et(1, 2) - et(1, 1)

    xn_new = xn
    yn_new = yn
    do n = 1, n_max
        do in = 2, in_max-1
            do jn = 2, jn_max-1
                alp = (1/(4*d_et**2))*((xn(in, jn+1) - xn(in, jn-1))**2 + (yn(in, jn+1) - yn(in, jn-1))**2)

                bet = (1/(4*d_xi*d_et))*((xn(in+1, jn) - xn(in-1, jn))*(xn(in, jn+1) - xn(in, jn-1)) + &
                      (yn(in+1, jn) - yn(in-1, jn))*(yn(in, jn+1) - yn(in, jn-1)))

                gam = (1/(4*d_xi**2))*((xn(in+1, jn) - xn(in-1,jn))**2 + (yn(in+1, jn) - yn(in-1, jn))**2)

                ! Calculate new (xn, yn) points using the Gauss-Seidel method
                xn_new(in, jn) = ((d_xi**2*d_et**2)/(2*(d_xi**2*gam + d_et**2*alp)))* &
                                ((alp/d_xi**2)*(xn(in+1, jn) + xn(in-1, jn)) + &
                                (bet/(2*d_xi*d_et))*(xn(in+1, jn-1) + xn(in-1, jn+1) - xn(in+1,jn+1) - xn(in-1,jn-1)) + &
                                (gam/d_et**2)*(xn(in, jn+1) + xn(in, jn-1)))

                yn_new(in, jn) = ((d_xi**2*d_et**2)/(2*(d_xi**2*gam + d_et**2*alp)))* &
                                ((alp/d_xi**2)*(yn(in+1, jn) + yn(in-1, jn)) + &
                                (bet/(2*d_xi*d_et))*(yn(in+1, jn-1) + yn(in-1, jn+1) - yn(in+1,jn+1) - yn(in-1,jn-1)) + &
                                (gam/d_et**2)*(yn(in, jn+1) + yn(in, jn-1)))
            end do
            ! Apply Neumann boundary conditions
            if ((xn(in, 1) <= fn_beg) .or. (xn(in, 1) >= fn_end)) then
                ! Ensures that the slope of the functions is considered outside of where the
                ! functions are actually applied. If this is not done, the grid squares near the
                ! beginning and end of the bumps are distorted and not squarelike. The yn direction
                ! is not considered since the Dirichlet boundary conditions determine the yn
                ! position when outside the domain of the functions.
                xn_new(in, 1)     = xn_new(in, 2)
                xn_new(in, jn_max) = xn_new(in, jn_max-1)
            else
                ! Ensures that the slope of the functions are considered where the functions are
                ! applied. Reorganizes the xn points to be normal to the surface of the bump. The yn
                ! points must also be considered so that the points do not neglect the BCs of the
                ! bump itself.
                xn_new(in, 1)     = xn_new(in, 2) + d_et*yp_btm(xn_new(in, 1))
                xn_new(in, jn_max) = xn_new(in, jn_max-1) - d_et*yp_top(xn_new(in, jn_max))
                yn_new(in, 1)     = y_btm(xn_new(in, 1))
                yn_new(in, jn_max) = y_top(xn_new(in, jn_max))
            end if
        end do
        
        ! Find the maximum error in both the xn and yn directions
        e_x = maxval(abs(xn_new - xn))
        e_y = maxval(abs(yn_new - yn))

        xn = xn_new
        yn = yn_new

        ! Break the loop if both xn and yn errors are below a set value
        if (e_x < tol .and. e_y < tol) then
            exit
        end if
    end do

    ! Find distance between the nodes on the first row (in.e. where jn = 1) to the beginning and end
    ! of the function. Since the domain is symmetric, the distances found on the bottom row are the
    ! same as those on the top row.
    do in = 1, in_max
        dist_l(in) = 2 - sqrt(xn(in, 1)**2 + yn(in, 1)**2)
        dist_r(in) = 3 - sqrt(xn(in, 1)**2 + yn(in, 1)**2)
    end do

    ! Find the location of the minimum value in the distance arrays.
    min_l = minloc(abs(dist_l), 1)
    min_r = minloc(abs(dist_r), 1)
    
    ! Use the index of the minimum position to replace 4 of the nodes in the mesh. Ensures that the
    ! functions boundary conditions are respected. If this is not done, the functions will not
    ! always begins and end at exactly xn = 2 and xn = 3 respectively.
    xn(min_l, jn_max) = fn_beg
    yn(min_l, jn_max) = yn_max
    xn(min_r, jn_max) = fn_end
    yn(min_r, jn_max) = yn_max
    xn(min_l, 1) = fn_beg
    yn(min_l, 1) = yn_min
    xn(min_r, 1) = fn_end
    yn(min_r, 1) = yn_min

    open (8, file='../../data/xyelliptic.x')
    write(8,*) in_max, jn_max
    write(8,*) ((xn(in,jn), in=1,in_max), jn=1,jn_max), ((yn(in,jn), in=1,in_max), jn=1,jn_max)
    close(8)

    call cpu_time(end)
    print *, 'subroutine elliptic took ', end - start, ' seconds'

end subroutine elliptic