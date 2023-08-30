subroutine residual()
    use grid_vars, only: ic, jc, ic_max, jc_max, xn, yn
    use flux_vars, only: f, g, res
    use functions
    implicit none

    ! calculates the residual term using the f and g fluxes on the four faces of each cell

    ! used to store node locations (in, jn)
    integer :: node_a(2), node_b(2), node_c(2), node_d(2)

    do ic = 1, ic_max
        do jc = 1, jc_max
            ! convert from cell notation (ic, jc) to node notation (in, jn)
            node_a = ijcell_to_ijnode(ic, jc, 1)
            node_b = ijcell_to_ijnode(ic, jc, 2)
            node_c = ijcell_to_ijnode(ic, jc, 3)
            node_d = ijcell_to_ijnode(ic, jc, 4)

            ! uses f and g fluxes on each cell face for a total of 8 terms
            res(ic, jc, :) = 0.5 * ((f(ic, jc, :) + f(ic + 1, jc, :)) * &
                                    (yn(node_c(1), node_c(2)) - yn(node_b(1), node_b(2))) - &
                                    (g(ic, jc, :) + g(ic + 1, jc, :)) * &
                                    (xn(node_c(1), node_c(2)) - xn(node_b(1), node_b(2))) + &
                                    (f(ic, jc, :) + f(ic, jc + 1, :)) * &
                                    (yn(node_d(1), node_d(2)) - yn(node_c(1), node_c(2))) - &
                                    (g(ic, jc, :) + g(ic, jc + 1, :)) * &
                                    (xn(node_d(1), node_d(2)) - xn(node_c(1), node_c(2))) + &
                                    (f(ic, jc, :) + f(ic - 1, jc, :)) * &
                                    (yn(node_a(1), node_a(2)) - yn(node_d(1), node_d(2))) - &
                                    (g(ic, jc, :) + g(ic - 1, jc, :)) * &
                                    (xn(node_a(1), node_a(2)) - xn(node_d(1), node_d(2))) + &
                                    (f(ic, jc, :) + f(ic, jc - 1, :)) * &
                                    (yn(node_b(1), node_b(2)) - yn(node_a(1), node_a(2))) - &
                                    (g(ic, jc, :) + g(ic, jc - 1, :)) * &
                                    (xn(node_b(1), node_b(2)) - xn(node_a(1), node_a(2))))
        end do
    end do

end subroutine residual