subroutine rk4()
    use mod_types, only: wp => dp
    use grid_vars,  only: ic, jc, ic_max, jc_max, area, dt_min
    use input,     only: n_iter, res_str
    use flux_vars,    only: q, res, dis
    use functions
    use timing
    implicit none

    ! Performs numerical integration on the state vector, residual, and dissipation terms using a
    ! modified version of the Runge-Kutta 4 method. First, a maximum time step is found for each
    ! cell, after which the smallest value is chosen as the global delta-t. The RK method is then
    ! applied to the state vector to update it in time.

    integer :: iter, nstages, i
    real(wp) :: a_vals(4), err
    real(wp), allocatable :: qold(:, :, :)

    call system_clock(start, rate)

    allocate(qold(ic_max, jc_max, 4))

    a_vals(1) = 0.25_wp
    a_vals(2) = 0.3333333333333333_wp
    a_vals(3) = 0.5_wp
    a_vals(4) = 1.0_wp

    nstages = 4

    open(1, file=res_str)

    write(1,'(*(g0.15,:,","))') 'iter', 'rho', 'rho_u', 'rho_v', 'rho_e'

    call bcinout
    call bcwall
    call flux
    call residual
    call dissipation

    qold = q

    do iter = 1, n_iter
        do i = 1, nstages
            do ic = 1, ic_max
                do jc = 1, jc_max
                    q(ic, jc, :) = qold(ic, jc, :) - a_vals(i)*dt_min*(res(ic, jc, :) - dis(ic, jc, :))/area(ic, jc)
                end do
            end do

            call bcinout
            call bcwall
            call flux
            call residual

        end do

        call dissipation

        err = maxval(abs(qold - q))
        qold = q
        ! print *, iter
        write(1,'(*(g0.15,:,","))') iter, maxval(abs(res(:, :, 1))), maxval(abs(res(:, :, 2))), maxval(abs(res(:, :, 3))), &
                                          maxval(abs(res(:, :, 4)))

        if (err <= 10.0_wp**(-8)) then
            exit
        end if
    end do

    call bcinout
    call bcwall

    close(1)

    call system_clock(end)
    print *, 'subroutine rk4 took ', (end - start) / rate, ' seconds'

end subroutine rk4