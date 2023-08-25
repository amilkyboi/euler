subroutine ghost()
    use gridprop, only: in_max, jn_max, in, jn, xn, yn
    use timing
    implicit none

    ! Adds two ghost cells to all four sides of the computational domain.

    call cpu_time(start)

    ! front
    do in = -1, 0
        do jn = 1, jn_max
            xn(in, jn) = (5.d0/(in_max - 1))*(in - 1)
            yn(in, jn) = yn(1, jn)
        end do
    end do

    ! back
    do in = in_max+1, in_max+2
        do jn = 1, jn_max
            xn(in, jn) = xn(in_max, jn) + (5.d0/(in_max - 1))*(in - in_max)
            yn(in, jn) = yn(1, jn)
        end do
    end do

    ! top
    do in = -1, in_max+2
        do jn = jn_max+1, jn_max+2
            xn(in, jn) = xn(in, 1)
            yn(in, jn) = yn(in, jn_max) + (yn(in, jn_max) - yn(in, 1))*(1.d0/(jn_max - 1))*(jn - jn_max)
        end do
    end do

    ! bottom
    do in = -1, in_max+2
        do jn = -1, 0
            xn(in, jn) = xn(in, 1)
            yn(in, jn) = yn(in, 1) + (yn(in, jn_max) - yn(in, 1))*(1.d0/(jn_max - 1))*(jn - 1)
        end do
    end do

    ! write the extended grid to a separate output file for verification
    ! open (8, file='../../data/xyghost.x')
    ! write(8,*) in_max+4, jn_max+4
    ! write(8,*) ((xn(in,jn), in=-1,in_max+2), jn=-1,jn_max+2), ((yn(in,jn), in=-1,in_max+2), jn=-1,jn_max+2)
    ! close(8)

    call cpu_time(end)
    print *, 'subroutine ghost took ', end - start, ' seconds'

end subroutine ghost