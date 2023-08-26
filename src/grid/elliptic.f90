subroutine elliptic()
    use mod_types, only: wp => dp
    use gridprop,  only: xi, et, xn, yn, in, jn, yn_min, yn_max, fn_beg, fn_end
    use input,     only: max_iter, in_max, jn_max, tol, y_bounds, fn_bounds, save_elp, elp_path
    use functions
    use timing
    implicit none

    ! distance arrays for ensuring boundary conditions, temporary grid storage arrays
    real(wp), allocatable :: dist_left(:), dist_rght(:), xn_new(:, :), yn_new(:, :)

    ! minimum distance indicies, iterations
    integer(wp) :: min_dist_left_idx, min_dist_rght_idx, iter

    ! elliptic grid parameters, grid steps, and errors
    real(wp) :: alp, bet, gam, d_xi, d_et, err_x, err_y

    call system_clock(start, rate)

    ! distance arrays for ensuring boundary conditions
    allocate(dist_left(in_max))
    allocate(dist_rght(in_max))
    ! temporary grid storage arrays
    allocate(xn_new(in_max, jn_max))
    allocate(yn_new(in_max, jn_max))

    ! variables for clarity
    fn_beg = fn_bounds(1)
    fn_end = fn_bounds(2)

    yn_min = y_bounds(1)
    yn_max = y_bounds(2)

    ! determine the step size in the xi and et directions
    d_xi = xi(2, 1) - xi(1, 1)
    d_et = et(1, 2) - et(1, 1)

    xn_new = xn
    yn_new = yn

    do iter = 1, max_iter
        do in = 2, in_max-1
            do jn = 2, jn_max-1
                ! find the elliptic grid parameters on each iteration for each cell; yes, it's slow
                ! but the grid produced is more well-behaved and the extra time is worth it
                alp = (1/(4*d_et**2))*((xn(in, jn+1) - xn(in, jn-1))**2 + &
                           (yn(in, jn+1) - yn(in, jn-1))**2)

                bet = (1/(4*d_xi*d_et))*((xn(in+1, jn) - xn(in-1, jn))*(xn(in, jn+1) - &
                xn(in, jn-1)) + &
                (yn(in+1, jn) - yn(in-1, jn))*(yn(in, jn+1) - yn(in, jn-1)))

                gam = (1/(4*d_xi**2))*((xn(in+1, jn) - xn(in-1, jn))**2 + &
                                    (yn(in+1, jn) - yn(in-1, jn))**2)

                ! calculate new (xn, yn) points using the Gauss-Seidel method
                xn_new(in, jn) = ((d_xi**2*d_et**2)/(2*(d_xi**2*gam + d_et**2*alp)))* &
                                ((alp/d_xi**2)*(xn(in+1, jn) + xn(in-1, jn)) + &
                                (bet/(2*d_xi*d_et))*(xn(in+1, jn-1) + xn(in-1, jn+1) - xn(in+1,jn+1) - xn(in-1,jn-1)) + &
                                (gam/d_et**2)*(xn(in, jn+1) + xn(in, jn-1)))

                yn_new(in, jn) = ((d_xi**2*d_et**2)/(2*(d_xi**2*gam + d_et**2*alp)))* &
                                ((alp/d_xi**2)*(yn(in+1, jn) + yn(in-1, jn)) + &
                                (bet/(2*d_xi*d_et))*(yn(in+1, jn-1) + yn(in-1, jn+1) - yn(in+1,jn+1) - yn(in-1,jn-1)) + &
                                (gam/d_et**2)*(yn(in, jn+1) + yn(in, jn-1)))
            end do

            ! apply Neumann boundary conditions
            if ((xn(in, 1) < fn_beg) .or. (xn(in, 1) > fn_end)) then
                ! ensures that the slope of the functions is considered outside of where the
                ! functions are actually applied
                
                ! if this is not done, the grid squares near the beginning and end of the bumps are
                ! distorted and not squarelike

                ! the yn direction is not considered since the Dirichlet boundary conditions
                ! determine the yn position when outside the domain of the functions
                xn_new(in, 1)      = xn_new(in, 2)
                xn_new(in, jn_max) = xn_new(in, jn_max-1)
            else
                ! ensures that the slope of the functions are considered where the functions are
                ! applied

                ! reorganizes the xn points to be normal to the surface of the bump
                
                ! the yn points must also be considered so that the points do not neglect the BCs of
                ! the bump itself
                xn_new(in, 1)      = xn_new(in, 2)          + d_et * yp_btm(xn_new(in, 1))
                xn_new(in, jn_max) = xn_new(in, jn_max - 1) - d_et * yp_top(xn_new(in, jn_max))
                yn_new(in, 1)      = y_btm(xn_new(in, 1))
                yn_new(in, jn_max) = y_top(xn_new(in, jn_max))
            end if
        end do
        
        ! find the maximum error in both the xn and yn directions
        err_x = maxval(abs(xn_new - xn))
        err_y = maxval(abs(yn_new - yn))

        ! reset the grid arrays
        xn = xn_new
        yn = yn_new

        ! break the loop if tolerance is reached
        if (err_x < tol .and. err_y < tol) then
            exit
        end if
    end do

    ! find distance between the nodes on the first row (i.e. where jn = 1) to the beginning and end
    ! of the function

    ! since the domain is symmetric, the distances found on the bottom row are the same as those on
    ! the top row
    do in = 1, in_max
        dist_left(in) = fn_beg - sqrt(xn(in, 1)**2 + yn(in, 1)**2)
        dist_rght(in) = fn_end - sqrt(xn(in, 1)**2 + yn(in, 1)**2)
    end do

    ! find the location of the minimum value in the distance arrays
    min_dist_left_idx = minloc(abs(dist_left), 1)
    min_dist_rght_idx = minloc(abs(dist_rght), 1)
    
    ! use the index of the minimum position to replace 4 of the nodes in the mesh
    
    ! if this is not done, the functions will not always begin and end where specified
    xn(min_dist_left_idx, jn_max) = fn_beg
    yn(min_dist_left_idx, jn_max) = yn_max
    xn(min_dist_rght_idx, jn_max) = fn_end
    yn(min_dist_rght_idx, jn_max) = yn_max
    xn(min_dist_left_idx, 1) = fn_beg
    yn(min_dist_left_idx, 1) = yn_min
    xn(min_dist_rght_idx, 1) = fn_end
    yn(min_dist_rght_idx, 1) = yn_min

    ! optionally save the elliptic mesh
    if (save_elp) then
        open (8, file=elp_path)
        write(8,*) in_max, jn_max
        write(8,*) ((xn(in,jn), in=1,in_max), jn=1,jn_max), ((yn(in,jn), in=1,in_max), jn=1,jn_max)
        close(8)
    end if

    call system_clock(end)
    print *, 'subroutine elliptic took ', (end - start) / rate, ' seconds'

end subroutine elliptic