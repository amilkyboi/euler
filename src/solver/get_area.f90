subroutine get_area()
    use grid_vars, only: ic, jc, ic_max, jc_max, xn, yn, area
    use reference, only: l_ref
    use functions
    use timing
    implicit none

    ! finds the area of each cell, including ghost cells

    ! used to store node locations (in, jn)
    integer :: node_a(2), node_b(2), node_c(2), node_d(2)

    call system_clock(start, rate)

    ! loop over all cells, including ghost cells
    do ic = -1, ic_max + 2
        do jc = -1, jc_max + 2
            ! convert from cell notation (ic, jc) to node notation (in, jn)
            node_a = ijcell_to_ijnode(ic, jc, 1)
            node_b = ijcell_to_ijnode(ic, jc, 2)
            node_c = ijcell_to_ijnode(ic, jc, 3)
            node_d = ijcell_to_ijnode(ic, jc, 4)

            ! simple area formula using the cross product of the two distance vectors
            area(ic, jc) = 0.5 * ((xn(node_c(1), node_c(2)) - xn(node_a(1), node_a(2))) * &
                                (yn(node_d(1), node_d(2)) - yn(node_b(1), node_b(2))) - &
                                (xn(node_d(1), node_d(2)) - xn(node_b(1), node_b(2))) * &
                                (yn(node_c(1), node_c(2)) - yn(node_a(1), node_a(2))))
        end do
    end do

    ! non-dimensionalize the entire area vector with respect to the reference length
    area = area / l_ref**2

    call system_clock(end)
    print *, 'subroutine area took ', (end - start) / rate, ' seconds'

end subroutine get_area