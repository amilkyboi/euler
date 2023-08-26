module functions
    use mod_types, only: wp => dp
    use gas_vars,  only: gammam1
    use grid_vars, only: xn, yn
    implicit none

    contains

    function ijcell_to_ijnode(ic, jc, node)
        ! converts from cell notation (ic, jc) to node notation (in, jn)

        integer, intent(in) :: ic, jc, node
        integer :: ijcell_to_ijnode(2)

        if (node == 1) then
            ! southwest corner, location a
            ijcell_to_ijnode = [ic    , jc    ]
        else if (node == 2) then
            ! southeast corner, location b
            ijcell_to_ijnode = [ic + 1, jc    ]
        else if (node == 3) then
            ! northeast corner, location c
            ijcell_to_ijnode = [ic + 1, jc + 1]
        else if (node == 4) then
            ! northwest corner, location d
            ijcell_to_ijnode = [ic    , jc + 1]
        else
            error stop
        end if

    end function ijcell_to_ijnode

    function length(ic, jc, dir)
        ! calculates the length of a specified cell face given cell indices

        integer, intent(in) :: ic, jc, dir
        integer :: node_a(2), node_b(2), node_c(2), node_d(2)

        ! float length
        real(wp) :: length

        ! convert from cell indices (ic, jc) to node indices (in, jn)
        node_a = ijcell_to_ijnode(ic, jc, 1)
        node_b = ijcell_to_ijnode(ic, jc, 2)
        node_c = ijcell_to_ijnode(ic, jc, 3)
        node_d = ijcell_to_ijnode(ic, jc, 4)

        if (dir == 1) then
            ! north face
            length = sqrt((xn(node_c(1), node_c(2)) - xn(node_d(1), node_d(2)))**2 + &
                          (yn(node_c(1), node_c(2)) - yn(node_d(1), node_d(2)))**2)
        else if (dir == 2) then
            ! south face
            length = sqrt((xn(node_b(1), node_b(2)) - xn(node_a(1), node_a(2)))**2 + &
                          (yn(node_b(1), node_b(2)) - yn(node_a(1), node_a(2)))**2)
        else if (dir == 3) then
            ! east face
            length = sqrt((xn(node_c(1), node_c(2)) - xn(node_b(1), node_b(2)))**2 + &
                          (yn(node_c(1), node_c(2)) - yn(node_b(1), node_b(2)))**2)
        else if (dir == 4) then
            ! west face
            length = sqrt((xn(node_d(1), node_d(2)) - xn(node_a(1), node_a(2)))**2 + &
                          (yn(node_d(1), node_d(2)) - yn(node_a(1), node_a(2)))**2)
        else
            error stop
        end if

    end function length

    function normal(ic, jc, dir)
        ! finds the direction cosines in the x-direction and y-direction for the outward pointing
        ! normal vector of a specified cell face

        integer, intent(in) :: ic, jc, dir
        integer :: node_a(2), node_b(2), node_c(2), node_d(2)

        ! normal vector
        real(wp) :: normal(2)

        ! convert from cell indices (ic, jc) to node indices (in, jn)
        node_a = ijcell_to_ijnode(ic, jc, 1)
        node_b = ijcell_to_ijnode(ic, jc, 2)
        node_c = ijcell_to_ijnode(ic, jc, 3)
        node_d = ijcell_to_ijnode(ic, jc, 4)

        if (dir == 1) then
            ! North
            normal = [(yn(node_d(1), node_d(2)) - yn(node_c(1), node_c(2))) / length(ic, jc, dir), &
                     -(xn(node_d(1), node_d(2)) - xn(node_c(1), node_c(2))) / length(ic, jc, dir)]
        else if (dir == 2) then
            ! South
            normal = [(yn(node_b(1), node_b(2)) - yn(node_a(1), node_a(2))) / length(ic, jc, dir), &
                     -(xn(node_b(1), node_b(2)) - xn(node_a(1), node_a(2))) / length(ic, jc, dir)]
        else if (dir == 3) then
            ! East
            normal = [(yn(node_c(1), node_c(2)) - yn(node_b(1), node_b(2))) / length(ic, jc, dir), &
                     -(xn(node_c(1), node_c(2)) - xn(node_b(1), node_b(2))) / length(ic, jc, dir)]
        else if (dir == 4) then
            ! West
            normal = [(yn(node_a(1), node_a(2)) - yn(node_d(1), node_d(2))) / length(ic, jc, dir), &
                     -(xn(node_a(1), node_a(2)) - xn(node_d(1), node_d(2))) / length(ic, jc, dir)]
        else
            error stop
        end if

    end function normal

end module functions