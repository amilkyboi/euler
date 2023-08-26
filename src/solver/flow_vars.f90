module flow_vars
    use mod_types, only: wp => dp
    implicit none

    ! declares all variables directly related to the flow

    ! density, x-velocity, y-velocity, velocity magnitude, specific energy,
    ! static temperature, speed of sound, mach number, static pressure, and specific entropy
    real(wp), allocatable :: dens(:, :), xvel(:, :), yvel(:, :), vmag(:, :), enrg(:, :), &
                             temp(:, :), vsnd(:, :), mach(:, :), pres(:, :), entr(:, :)
end module flow_vars