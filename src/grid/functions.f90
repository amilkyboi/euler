module functions
    use mod_types, only: wp => dp
    implicit none

    real(wp), parameter :: pi = 4.0*atan(1.0)

    contains

    ! function defining the top of the domain, y_top(x)
    real(wp) function y_top(x)
        real(wp), intent(in) :: x

        y_top = 1 - 0.1 * sin(pi * (x - 2))

    end function y_top

    ! y'_top(x)
    real(wp) function yp_top(x)
        real(wp), intent(in) :: x

        yp_top = -0.1 * pi * cos(pi * (x - 2))

    end function yp_top

    ! function defining the bottom of the domain, y_btm(x)
    real(wp) function y_btm(x)
        real(wp), intent(in) :: x

        y_btm = 0.1 * sin(pi * (x - 2))

    end function y_btm

    ! y'_btm(x)
    real(wp) function yp_btm(x)
        real(wp), intent(in) :: x

        yp_btm = 0.1 * pi * cos(pi * (x - 2))

    end function yp_btm

end module functions