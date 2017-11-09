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
  integer :: i, j, k, icell, it, ilog, nkall
  real(prec) :: summ, dtmp

  ! WRITE(*,*) "Switch Time coefficiet in set_vels_time_space !!!"

  nkall = nmodes * ntimes
  ilog = 0
  do i= 1, ncell * ntimes
    icell = MOD(i-1 , ncell) + 1
    it = INT((i-1) / ncell) + 1
    if (it /= ilog) THEN
      write(*,*) " In Cycle time",it, icell
      ilog = it
    END IF
    cs_m(1:nkall) = b_m(1,1:nkall)*xp(1,icell) +&
                     b_m(2,1:nkall)*xp(2,icell) + &
                     b_m(3,1:nkall)*xp(3,icell) + &
                     c_m(1:nkall) * dtim(it)
    do j= 1, 3
      u(j,i) = 2.0d0 * SUM(ac_m(j,1:nkall)*COS(cs_m(1:nkall)) &
        + as_m(j,1:nkall)*SIN(cs_m(1:nkall)))
    enddo
  enddo
  RETURN
end subroutine set_vels_time_space