module grid_vars
    use mod_types, only: wp => dp
    implicit none

    ! Declares all variables directly related to the geometry of the grid.

    ! i of node, j of node, i of cell, j of cell, max i node, max j node, max i cell, max j cell
    integer :: in, jn, ic, jc, in_max, jn_max, ic_max, jc_max

    ! node x-locations, node y-locations, cell area, cell angle of attack
    real(wp), allocatable :: xn(:, :), yn(:, :), area(:, :), dt(:, :)
    ! minimum time step
    real(wp) :: dt_min

end module grid_vars
