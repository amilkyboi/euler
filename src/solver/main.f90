program main
    implicit none

    call startup

    call allocation

    call ghost

    call get_area

    call initializer

    call timestep

    call rk4

    call plt

end program main