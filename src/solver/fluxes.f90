module fluxes
    use mod_types, only: wp => dp
    implicit none

    ! Declares all variables that involve the state vector or flux.

    ! cell state vector, cell f flux, cell g flux, cell residual, cell dissipation
    real(wp), allocatable :: q(:, :, :), f(:, :, :), g(:, :, :), res(:, :, :), dis(:, :, :)

end module fluxes