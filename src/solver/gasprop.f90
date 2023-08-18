module gasprop
    use mod_types, only: wp => dp
    implicit none

    ! Declares all variables related to the gas properties of the flow.

    ! ratio of specific heats
    real(wp), parameter :: gamma = 1.4_wp
    ! ratio of specific heats minus 1
    real(wp), parameter :: gammam1 = 0.4_wp
end module gasprop