subroutine ghost()
    use grid_vars, only: in_max, jn_max, in, jn, xn, yn
    use timing
    implicit none

    ! adds two ghost cells to all four sides of the computational domain

    call system_clock(start, end)

    ! TODO: creates the cells based on the 5x1 geometry only, need to change this by adding an input
    !       variable or something

    ! add front two columns of ghost cells
    do in = -1, 0
        do jn = 1, jn_max
            xn(in, jn) = (5.d0 / (in_max - 1)) * (in - 1)
            yn(in, jn) = yn(1, jn)
        end do
    end do

    ! add back two columns of ghost cells
    do in = in_max + 1, in_max + 2
        do jn = 1, jn_max
            xn(in, jn) = xn(in_max, jn) + (5.d0 / (in_max - 1)) * (in - in_max)
            yn(in, jn) = yn(1, jn)
        end do
    end do

    ! add top two rows of ghost cells
    do in = -1, in_max + 2
        do jn = jn_max + 1, jn_max + 2
            xn(in, jn) = xn(in, 1)
            yn(in, jn) = yn(in, jn_max) + (yn(in, jn_max) - yn(in, 1)) * &
                                          (1.d0 / (jn_max - 1)) * (jn - jn_max)
        end do
    end do

    ! add bottom two rows of ghost cells
    do in = -1, in_max + 2
        do jn = -1, 0
            xn(in, jn) = xn(in, 1)
            yn(in, jn) = yn(in, 1) + (yn(in, jn_max) - yn(in, 1)) * (1.d0 / (jn_max - 1)) * (jn - 1)
        end do
    end do

    call system_clock(end)
    print *, 'subroutine ghost took ', (end - start) / rate, ' seconds'

end subroutine ghost