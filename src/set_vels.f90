!>@brief Calculate velocities in 3D space ()
subroutine set_vels()
  USE prec_mod
  USE comm1, ONLY: xp, u
  USE tmp_mod
  implicit none
  integer :: i, j, k
  real(prec) :: summ, dtmp

  do j = 1, 3
    do i= 1, tcell
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


!>@brief Calculate velocities in 4Dim. space (XYZ+Time)
subroutine set_vels_time_space()
  USE prec_mod
  USE comm1, ONLY: xp, u
  USE tmp_mod
  implicit none
  integer :: i, j, k, icell, it, ilog, nkall
  real(prec) :: summ, dtmp

  ! WRITE(*,*) "Switch Time coefficiet in set_vels_time_space !!!"

  nkall = nmodes * ntimes
  ilog = 0
  do i= 1, tcell * ntimes
    icell = MOD(i-1 , tcell) + 1
    it = INT((i-1) / tcell) + 1
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


!>@brief Calculate velocities in 4Dim. space (XYZ+Time)
subroutine set_vels_at_time(dtim_in)
  
  USE comm0, ONLY: ncell
  USE comm1, ONLY: xp, u
  USE prec_mod
  USE fire_mpi_var
  USE tmp_mod
  
  implicit none
  
  real(prec), intent(in) :: dtim_in !< current time step
  
  integer :: i, j, k, icell, it, ilog, nkall, il
  real(prec) :: summ, dtmp

  ! WRITE(*,*) "Switch Time coefficiet in set_vels_time_space !!!"

  nkall = nmodes * ntimes
  ilog = 0

  ilog = INT(ncell**(0.66666667))
  il = INT(ncell**(1.d0/3.0d0))
  ! cycle over cells
  do i= 1, ncell !!!!* ntimes
    icell = MOD(i-1 , tcell) + 1
    ! it = INT((i-1) / tcell) + 1
    ! if (it /= ilog) THEN
    !   write(*,*) " In Cycle time",it, icell
    !   ilog = it
    ! END IF
    k = MOD(icell-1,ilog)
    if (mpi_master .and. k.eq.0) then
      write(*,'(2(a,i))') "Gen vels : ",INT(icell/ilog),"/",il
    endif
    cs_m(1:nkall) =  b_m(1,1:nkall)*xp(1,icell) + &
                     b_m(2,1:nkall)*xp(2,icell) + &
                     b_m(3,1:nkall)*xp(3,icell) + &
                     c_m(1:nkall) * dtim_in
    do j= 1, 3
      u(j,i) = 2.0d0 * SUM( ac_m(j,1:nkall) * COS(cs_m(1:nkall)) &
        + as_m(j,1:nkall) * SIN(cs_m(1:nkall)) )
    enddo
  enddo
  RETURN
end subroutine set_vels_at_time


!>@brief Calculate velocities in 4Dim. space (XYZ+Time)
subroutine set_vels_at_space_time(dtim_in,dxx,vels)
  
  USE comm0, ONLY: ncell
  USE comm1, ONLY: xp, u
  USE prec_mod
  USE fire_mpi_var
  USE tmp_mod
  
  implicit none
  
  real(prec), intent(in) :: dtim_in !< current time step
  real(prec),dimension(1:3), intent(inout) :: dxx !< current time step
  real(prec),dimension(1:3), intent(inout) :: vels !< current time step
  
  integer :: i, j, k, icell, it, ilog, nkall, il
  real(prec) :: summ, dtmp

  ! WRITE(*,*) "Switch Time coefficiet in set_vels_time_space !!!"

  nkall = nmodes * ntimes

  cs_m(1:nkall) =  b_m(1,1:nkall)*dxx(1) + &
                   b_m(2,1:nkall)*dxx(2) + &
                   b_m(3,1:nkall)*dxx(3) + &
                   c_m(1:nkall) * dtim_in
  do j= 1, 3
    vels(j) = 2.0d0 * SUM( ac_m(j,1:nkall) * COS(cs_m(1:nkall)) &
      + as_m(j,1:nkall) * SIN(cs_m(1:nkall)) )
  enddo
  RETURN
end subroutine set_vels_at_space_time
