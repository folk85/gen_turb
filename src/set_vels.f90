!@brief Calculate velocities
subroutine set_vels()
  USE prec_mod
  USE tmp_mod
  implicit none
  integer :: i, j

  do i= 1, ncell
    do j = 1, 3
      u(j,i) = 2.0d0 * SUM(ac_m(j,1:nmodes) * COS( b_m(j,1:nmodes)*xp(j,i) + &
        c_m(j,1:nmodes)))  
    enddo
  enddo
  RETURN
end subroutine set_vels