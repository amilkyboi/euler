subroutine rk4()
    use mod_types, only: wp => dp
    use grid_vars, only: ic, jc, ic_max, jc_max, area, dt_min
    use input,     only: max_iter, res_tol, res_str, frc_str
    use flux_vars, only: q, res, dis
    use reference, only: l_ref
    use flow_vars, only: pres
    use functions
    use timing
    implicit none

    ! performs numerical integration on the state vector, residual, and dissipation terms using a
    ! modified version of the Runge-Kutta 4 method

    ! local variables for iteration number, number of RK4 stages, and current stage
    integer :: iter, n_stages, stage

    ! multiplicative constants for the modified Runge-Kutta algorithm, error is the difference
    ! between the current and previous density residual, north-facing normal
    real(wp) :: a_vals(4), nn(2)
    ! temporary state vector storage for the RK4 algorithm and error analysis
    real(wp), allocatable :: q_old(:, :, :), res_dis_old(:, :, :), res_dis_err(:, :, :), fp_x(:), &
                             fp_y(:)

    call system_clock(start, rate)

    ! previous state and residual vectors
    allocate(q_old(ic_max, jc_max, 4))
    allocate(res_dis_old(ic_max, jc_max, 4))
    allocate(res_dis_err(ic_max, jc_max, 4))
    allocate(fp_x(ic_max))
    allocate(fp_y(ic_max))

    ! RK constants
    a_vals(1) = 0.25_wp
    a_vals(2) = 0.33333333333333333_wp
    a_vals(3) = 0.5_wp
    a_vals(4) = 1.0_wp

    ! number of RK stages
    n_stages = 4

    ! apply boundary conditions before anything else to ensure domain conformance
    call bcinout
    call bcwall
    call flux
    call residual
    call dissipation

    ! set previous state and residual vectors
    q_old       = q
    res_dis_old = res - dis

    ! file that contains iteration number and residual quantities for plotting
    open(1, file=res_str)
    write(1,'(*(g0.15,:,","))') 'iter', 'rho', 'rho_u', 'rho_v', 'rho_e', 'err_rho', 'err_rho_u', &
                                'err_rho_v', 'err_rho_e'

    do iter = 1, max_iter

        ! dissipation must be called before the state vector is updated in the RK algorithm
        call dissipation

        ! RK4
        do stage = 1, n_stages
            do ic = 1, ic_max
                do jc = 1, jc_max
                    q(ic, jc, :) = q_old(ic, jc, :) - a_vals(stage) * dt_min * (res(ic, jc, :) - dis(ic, jc, :)) / area(ic, jc)
                end do
            end do

            ! apply boundary conditions before calculating flux and residual
            call bcinout
            call bcwall
            call flux
            call residual

        end do

        ! track residual minus dissipation as the error quantity of interest
        res_dis_err = res_dis_old - (res - dis)

        ! reset tracking vectors
        q_old       = q
        res_dis_old = res - dis

        ! write maximum residuals across the entire domain for each quantity along with the maximum
        ! error between the current and previous iteration
        write(1,'(*(g0.15,:,","))') iter, maxval(abs(res(:, :, 1) - dis(:, :, 1))), &
                                          maxval(abs(res(:, :, 2) - dis(:, :, 2))), &
                                          maxval(abs(res(:, :, 3) - dis(:, :, 3))), &
                                          maxval(abs(res(:, :, 4) - dis(:, :, 4))), &
                                          maxval(abs(res_dis_err(:, :, 1))), &
                                          maxval(abs(res_dis_err(:, :, 2))), &
                                          maxval(abs(res_dis_err(:, :, 3))), &
                                          maxval(abs(res_dis_err(:, :, 4)))

        ! if the tolerance threshold is met before the maximum number of iterations is reached, the
        ! algorithm is complete
        if (maxval(abs(res_dis_err)) <= res_tol) then
            exit
        end if
    end do

    close(1)

    ! apply boundary conditions one final time
    call bcinout
    call bcwall

    open(2, file=frc_str)
    write(2,'(*(g0.15,:,","))') 'dist', 'fp_x', 'fp_y'

    ! force on the bottom wall
    ! F_x = -p * dA * n_x
    ! F_y = -p * dA * n_y
    ! dA  = len_cell * len_reference
    do ic = 1, ic_max
        ! upward pointing normal from the top face of the first bottom ghost cell
        nn = normal(ic, 0, 1)

        ! note that the magnitude of the force will decrease as the grid is refined since the length
        ! of each cell will get smaller
        fp_x(ic) = - pres(ic, 1) * length(ic, 0, 1) * l_ref * nn(1)
        fp_y(ic) = - pres(ic, 1) * length(ic, 0, 1) * l_ref * nn(2)

        write(2,'(*(g0.15,:,","))') xn(ic, 1), fp_x(ic), fp_y(ic)
    end do

    close(2)

    print *, iter
    print *, maxval(abs(res_dis_err))

    call system_clock(end)
    print *, 'subroutine rk4 took ', (end - start) / rate, ' seconds'

end subroutine rk4