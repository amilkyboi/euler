subroutine rk4()
    use mod_types, only: wp => dp
    use gridprop,  only: ic, jc, ic_max, jc_max, ar, dt_min
    use flowprop,  only: p
    use input,     only: n_iter, res_str, frc_str
    use fluxes,    only: q, res, dis
    use reference, only: l_ref
    use functions
    implicit none

    ! Performs numerical integration on the state vector, residual, and dissipation terms using a
    ! modified version of the Runge-Kutta 4 method. First, a maximum time step is found for each
    ! cell, after which the smallest value is chosen as the global delta-t. The RK method is then
    ! applied to the state vector to update it in time.

    integer :: iter, nstages, i
    real(wp) :: a_vals(4), err, nn(2)
    real(wp), allocatable :: qold(:, :, :), fp_x(:), fp_y(:)

    allocate(qold(ic_max, jc_max, 4), fp_x(ic_max), fp_y(ic_max))

    a_vals(1) = 0.25_wp
    a_vals(2) = 0.3333333333333333_wp
    a_vals(3) = 0.5_wp
    a_vals(4) = 1.0_wp

    nstages = 4

    open(1, file=res_str)

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
                    q(ic, jc, :) = qold(ic, jc, :) - a_vals(i)*dt_min*(res(ic, jc, :) - dis(ic, jc, :))/ar(ic, jc)
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

    open(2, file=frc_str)

    do ic = 1, ic_max
        if (xn(ic, 1) >= 2.0_wp .and. xn(ic, 1) <= 3.0_wp) then
            nn = normal(ic, 0, 1) ! upward facing normal from the bottom wall

            ! force due to pressure in x-direction
            fp_x(ic) = p(ic, 1) * l_ref * length(ic, 1, 2) * nn(1)

            ! force due to pressure in y-direction
            fp_y(ic) = p(ic, 1) * l_ref * length(ic, 1, 2) * nn(2)

            write(2,'(*(g0.15,:,","))') xn(ic, 1), fp_x(ic), fp_y(ic)
        end if
    end do

    close(2)

end subroutine rk4