!@brief Calculate velocities in 3D space ()
subroutine set_vels()
  USE prec_mod
  USE tmp_mod
  implicit none
  integer :: i, j, k
  real(prec) :: summ, dtmp

  do j = 1, 3
    do i= 1, ncell
      summ = 0.0d0
      do k= 1, nmodes
        dtmp = DOT_PRODUCT(b_m(1:3,k),xp(1:3,i)) + c_m(k)
        summ = summ + ac_m(j,k) * COS( dtmp) + as_m(j,k) * SIN( dtmp) 
      enddo   
      u(j,i) = 2.0d0 * summ  
    enddo
  enddo
  RETURN
end subroutine set_vels


!@brief Calculate velocities in 4Dim. space (XYZ+Time)
subroutine set_vels_time_space()
  USE prec_mod
  USE tmp_mod
  implicit none
  integer :: i, j, k, icell, itime, ilog
  real(prec) :: summ, dtmp

  ! WRITE(*,*) "Switch Time coefficiet in set_vels_time_space !!!"

  ilog = 0
  do j = 1, 3
    do i = 1, ncell * ntimes
      summ = 0.0d0
      icell = MOD(i-1 , ncell) + 1
      itime = INT((i-1) / ncell) + 1
      if (itime /= ilog) THEN
        write(*,*) " In Cycle time",itime, icell
        ilog = itime
      END IF
      do k = 1, nmodes * ntimes
        dtmp = DOT_PRODUCT(b_m(1:3,k),xp(1:3,i)) + c_m(k) * dtime(itime)
        ! dtmp = DOT_PRODUCT(b_m(1:3,k),xp(1:3,i)) + c_m(k)
        summ = summ + ac_m(j,k) * COS( dtmp) + as_m(j,k) * SIN( dtmp) 
      enddo   
      u(j,i) = 2.0d0 * summ  
    enddo
  enddo
  RETURN
end subroutine set_vels_time_space