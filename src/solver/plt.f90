subroutine plt()
   use mod_types, only: wp => dp
   use grid_vars, only: in, jn, ic, jc, in_max, jn_max, ic_max, jc_max, xn, yn
   use input,     only: plt_str
   use flow_vars, only: mach
   use flux_vars, only: q
   use timing
   implicit none

   integer :: i1, i2, i3, i4
   real(wp), allocatable :: qn(:, :, :), mn(:, :)

   call system_clock(start, rate)

   allocate(qn(in_max, jn_max, 4), mn(in_max, jn_max))

   do in = 1, in_max
      do jn = 1, jn_max
         qn(in, jn, :) = 0.25 * (q(in, jn, :) + q(in-1, jn, :) + q(in, jn-1, :) + q(in-1, jn-1, :))
         mn(in, jn) = 0.25 * (mach(in, jn) + mach(in-1, jn) + mach(in, jn-1) + mach(in-1, jn-1))
      end do
   end do

   open(33, file=plt_str)
   write(33, *) 'variables = "x", "y", "rho", "rho_u", "rho_v", "rho_e", "mach"'
   write(33, *) "zone f=fepoint, et=quadrilateral, n=" , in_max*jn_max, ", e=", ic_max*jc_max

   do jn = 1, jn_max
      do in = 1, in_max
         write(33, *) xn(in,jn), yn(in,jn), qn(in,jn,1), qn(in,jn,2), qn(in,jn,3), qn(in,jn,4), mn(in,jn)
      end do
   end do

   ! connectivity
   do jc = 1, jc_max
      do ic = 1, ic_max
         i1 = ic + (jc-1) * (ic_max+1)
         i2 = i1 + 1
         i3 = i2 + ic_max + 1
         i4 = i3 - 1
         write(33, *) i1, i2, i3, i4
      end do
   end do
   close(33)

   call system_clock(end)
   print *, 'subroutine plt took ', (end - start) / rate, ' seconds'

end subroutine plt