program main
    implicit none

    ! calls all relevant subroutines and executes the Euler solver

    call startup

    call allocation

    call ghost

    call area

    call initializer

    call timestep

    call rk4

    call plt

end program main