module functions
    use mod_types, only: wp => dp
    use gasprop,   only: gammam1
    use gridprop,  only: xn, yn
    implicit none

    ! Various useful functions.

    contains
    function ijnode(ic, jc, node)
        ! Converts from cell notation to node notation. Outputs an array of inode and jnode.

        integer, intent(in) :: ic, jc, node
        integer :: ijnode(2)

        if (node == 1) then
            ! SW / A
            ijnode = [ic  , jc  ]
        else if (node == 2) then
            ! SE / B
            ijnode = [ic+1, jc  ]
        else if (node == 3) then
            ! NE / C
            ijnode = [ic+1, jc+1]
        else if (node == 4) then
            ! NW / D
            ijnode = [ic  , jc+1]
        end if
    end function ijnode

    function length(ic, jc, dir)
        ! Calculates the length of a specified cell face.

        integer, intent(in) :: ic, jc, dir
        integer :: an(2), bn(2), cn(2), dn(2)
        real(wp) :: length
        
        an = ijnode(ic, jc, 1)
        bn = ijnode(ic, jc, 2)
        cn = ijnode(ic, jc, 3)
        dn = ijnode(ic, jc, 4)

        if (dir == 1) then
            ! North
            length = sqrt((xn(cn(1), cn(2)) - xn(dn(1), dn(2)))**2 + (yn(cn(1), cn(2)) - yn(dn(1), dn(2)))**2)
        else if (dir == 2) then
            ! South
            length = sqrt((xn(bn(1), bn(2)) - xn(an(1), an(2)))**2 + (yn(bn(1), bn(2)) - yn(an(1), an(2)))**2)
        else if (dir == 3) then
            ! East
            length = sqrt((xn(cn(1), cn(2)) - xn(bn(1), bn(2)))**2 + (yn(cn(1), cn(2)) - yn(bn(1), bn(2)))**2)
        else if (dir == 4) then
            ! West
            length = sqrt((xn(dn(1), dn(2)) - xn(an(1), an(2)))**2 + (yn(dn(1), dn(2)) - yn(an(1), an(2)))**2)
        else
            error stop
        end if
    end function length

    function normal(ic, jc, dir)
        ! Outputs the direction cosines in the x-direction and y-direction for the outward pointing
        ! normal vector of a specified cell face.

        integer, intent(in) :: ic, jc, dir
        integer :: an(2), bn(2), cn(2), dn(2)
        real(wp) :: normal(2)
        
        an = ijnode(ic, jc, 1)
        bn = ijnode(ic, jc, 2)
        cn = ijnode(ic, jc, 3)
        dn = ijnode(ic, jc, 4)

        if (dir == 1) then
            ! North
            normal = [(yn(dn(1), dn(2)) - yn(cn(1), cn(2)))/length(ic, jc, dir), &
                     -(xn(dn(1), dn(2)) - xn(cn(1), cn(2)))/length(ic, jc, dir)]
        else if (dir == 2) then
            ! South
            normal = [(yn(bn(1), bn(2)) - yn(an(1), an(2)))/length(ic, jc, dir), &
                     -(xn(bn(1), bn(2)) - xn(an(1), an(2)))/length(ic, jc, dir)]
        else if (dir == 3) then
            ! East
            normal = [(yn(cn(1), cn(2)) - yn(bn(1), bn(2)))/length(ic, jc, dir), &
                     -(xn(cn(1), cn(2)) - xn(bn(1), bn(2)))/length(ic, jc, dir)]
        else if (dir == 4) then
            ! West
            normal = [(yn(an(1), an(2)) - yn(dn(1), dn(2)))/length(ic, jc, dir), &
                     -(xn(an(1), an(2)) - xn(dn(1), dn(2)))/length(ic, jc, dir)]
        end if
    end function normal

end module functions