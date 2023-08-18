program main
    use timing
    implicit none

    ! Calls all relevant subroutines and executes the Euler solver.

    call cpu_time(start)

    ! reads the elliptic grid file, finds in_max and jn_max, allocates and fills xn and yn
    call startup

    ! allocates all arrays except for xn and yn that have already been allocated by startup
    call allocation

    ! adds two ghost cells to each side of the domain
    call ghost

    ! calculates the area matrix
    call area

    ! initial state vector and flux
    call initializer

    call timestep

    call rk4

    ! plot the solution
    call plt

    call cpu_time(end)
    print *, 'total program took ', end - start, ' seconds'
    print *, 'done'

end program main