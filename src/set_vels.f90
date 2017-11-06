!@brief Calculate velocities
subroutine set_vels()
  USE prec_mod
  USE tmp_mod
  implicit none
  integer :: i, j, k
  real(prec) :: summ

  do j = 1, 3
    do i= 1, ncell
      summ = 0.0d0
      do k= 1, nmodes
        summ = summ + ac_m(j,k) * COS( DOT_PRODUCT(b_m(1:3,k),xp(1:3,i)) + c_m(k))
      enddo   
      u(j,i) = 2.0d0 * summ  
    enddo
  enddo
  RETURN
end subroutine set_vels