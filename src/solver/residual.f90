subroutine residual()
    use mod_types, only: wp => dp
    use grid_vars,  only: ic, jc, ic_max, jc_max, xn, yn
    use flux_vars,    only: f, g, res
    use functions
    implicit none

    ! Calculates the residual term of the main equation using the f and g fluxes on the four faces
    ! of each cell.

    integer(wp) :: an(2), bn(2), cn(2), dn(2)

    do ic = 1, ic_max
        do jc = 1, jc_max
            ! grab the four node locations for each cell
            an = ijcell_to_ijnode(ic, jc, 1)
            bn = ijcell_to_ijnode(ic, jc, 2)
            cn = ijcell_to_ijnode(ic, jc, 3)
            dn = ijcell_to_ijnode(ic, jc, 4)

            res(ic, jc, :) = 0.5 * ((f(ic, jc, :) + f(ic+1, jc, :)) * (yn(cn(1), cn(2)) - yn(bn(1), bn(2))) &
                           -        (g(ic, jc, :) + g(ic+1, jc, :)) * (xn(cn(1), cn(2)) - xn(bn(1), bn(2))) &
                           +        (f(ic, jc, :) + f(ic, jc+1, :)) * (yn(dn(1), dn(2)) - yn(cn(1), cn(2))) &
                           -        (g(ic, jc, :) + g(ic, jc+1, :)) * (xn(dn(1), dn(2)) - xn(cn(1), cn(2))) &
                           +        (f(ic, jc, :) + f(ic-1, jc, :)) * (yn(an(1), an(2)) - yn(dn(1), dn(2))) &
                           -        (g(ic, jc, :) + g(ic-1, jc, :)) * (xn(an(1), an(2)) - xn(dn(1), dn(2))) &
                           +        (f(ic, jc, :) + f(ic, jc-1, :)) * (yn(bn(1), bn(2)) - yn(an(1), an(2))) &
                           -        (g(ic, jc, :) + g(ic, jc-1, :)) * (xn(bn(1), bn(2)) - xn(an(1), an(2))))
        end do
    end do
end subroutine residual