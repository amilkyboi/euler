module flowprop
    use mod_types, only: wp => dp
    implicit none

    ! Declares all variables directly related to the flow itself.

    ! density, x-velocity, y-velocity, velocity, specific energy, temperature, speed of sound,
    ! mach number, static pressure, and entropy
    real(wp), allocatable :: r(:, :), u(:, :), v(:, :), vel(:, :), E(:, :), T(:, :), c(:, :), &
                             mach(:, :), p(:, :), s(:, :)
end module flowprop