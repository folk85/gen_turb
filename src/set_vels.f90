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
subroutine set_vels_time_space(icase)
  USE prec_mod
  USE comm1, ONLY: xp, u
  USE tmp_mod
  implicit none
  integer, OPTIONAL :: icase
  integer :: i, j, k, icell, it, ilog, nkall, ic
  real(prec) :: summ, dtmp
  real(prec),dimension(1:3,1:tcell) :: u_tmp

  ! WRITE(*,*) "Switch Time coefficiet in set_vels_time_space !!!"

  ! if (PRESENT(icase)) then
  !   ic = icase
  ! else
  !   ic = 4
  ! endif
  do i = 1, ntimes
    icell = tcell
    ! write(*,'(a,4i7)') " In Cycle time new_vec", i , ntimes,tcell, icell
    CALL set_vels_at_space_time2(dtim(i),icell,xp(1:3,1:icell),u_tmp(1:3,1:tcell))
    u(1:3,tcell*(i-1)+1:tcell*i) = u_tmp(1:3,1:tcell)
  end do
  return

  ! IF (ic .eq. 4) THEN
    nkall = nkmod
    ilog = 0
    do i= 1, tcell * ntimes
      icell = MOD(i-1 , tcell) + 1
      it = INT((i-1) / tcell) + 1
      if (it /= ilog) THEN
        write(*,*) " In Cycle time new",it, icell
        ilog = it
      END IF
      
      ! Call the same routine as in USEINI
      CALL set_vels_at_space_time(dtim(it),xp(1:3,icell),u(1:3,i))
      ! CALL set_vels_at_space_time2(dtim(it),1,xp(1:3,icell),u(1:3,i))

    enddo
  ! ELSE IF (ic .eq. 5) THEN
  !   ! Call the same routine as in USEINI
  !   CALL set_vels_at_space_time(dtim(it),xp(1:3,icell),u(1:3,i))
  ! END IF !(ic .eq. 5) THEN

  RETURN
end subroutine set_vels_time_space


!>@brief Calculate velocities in 4Dim. space (XYZ+Time)
subroutine set_vels_at_time(dtim_in)
  
  USE comm0, ONLY: ncell,nsp,nep
  USE comm1, ONLY: xp, u
  USE prec_mod
  USE fire_mpi_var
  USE tmp_mod
  
  implicit none
  
  real(prec), intent(in) :: dtim_in !< current time step
  
  integer :: i, j, k, icell, it, ilog, nkall, il
  real(prec) :: summ, dtmp

  ! WRITE(*,*) "Switch Time coefficiet in set_vels_time_space !!!"

  nkall = nkmod
  ilog = 0

  ilog = INT(ncell**(0.66666667))
  il = INT(ncell**(1.d0/3.0d0))
  ! cycle over cells
  do i= nsp(1), nep(1) !!!!* ntimes
    icell = MOD(i-1 , tcell) + 1
    ! it = INT((i-1) / tcell) + 1
    ! if (it /= ilog) THEN
    !   write(*,*) " In Cycle time",it, icell
    !   ilog = it
    ! END IF
    k = MOD(icell-1,ilog)
    if (mpi_master .and. k.eq.0) then
      write(*,'(2(a20,i7))') "Gen vels : ",INT(icell/ilog),"/",il
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

  ! CALL exchng(u,3,1)
  ! CALL exchng(u,3,2)
  ! CALL exchng(u,3,3)
  RETURN
end subroutine set_vels_at_time


!>@brief Calculate velocities in 4Dim. space (XYZ+Time)
subroutine set_vels_at_space_time(dtim_in,dxx,vels)
  
  USE comm0, ONLY: ncell
  ! USE comm1, ONLY: xp, u
  USE prec_mod
  USE fire_mpi_var
  USE tmp_mod
  
  implicit none
  
  real(prec), intent(in) :: dtim_in !< current time step
  real(prec),dimension(1:3), intent(in) :: dxx !< current time step
  real(prec),dimension(1:3), intent(out) :: vels !< current time step
  
  integer :: i, j, k, icell, it, ilog, nkall, il
  real(prec),dimension(1:3) :: dv_cos!< current time step
  real(prec) :: summ, dtmp

  ! 1 - set number of modes in all times
  nkall = nkmod

  ! 2 - COS predicted vars
  cs_m(1:nkall) = c_m(1:nkall)

  ! 3 - Calc $ K \cdot x + w t$
  CALL dgemm('N','N',1,nkall,3,1.0d0,dxx,1,b_m,3,dtim_in,cs_m,1)
  ! write(*,'(a,8es13.5)')"A1: ", cs_m(1:8)
  ! 4 - calc sin and cos
  CALL vdsincos(nkall,cs_m(1),sin_m(1),cos_m(1))
  ! 5 - Calc COS  like $AC \cdot \cos{ CS }$
  CALL dgemm('N','N',3,1,nkall,1.0d0,ac_m,3,cos_m,nkall,0.0d0,dv_cos,3)
  ! 6 - Calc SIN  like $AS \cdot \sin{ CS }$ and add it to prev val.
  CALL dgemm('N','N',3,1,nkall,1.0d0,as_m,3,sin_m,nkall,1.0d0,dv_cos,3)
  vels(1:3) = dv_cos(1:3)
  RETURN
end subroutine set_vels_at_space_time


!>@brief Calculate velocities in 4Dim. space (XYZ+Time)
subroutine set_vels_at_space_time2(dtim_in,inel,dxx,vels)
  
  USE comm0, ONLY: ncell
  ! USE comm1, ONLY: xp, u
  USE prec_mod
  USE fire_mpi_var
  USE tmp_mod
  
  implicit none
  
  real(prec), intent(in) :: dtim_in !< current time step
  integer, intent(in) :: inel !< current time step
  real(prec),dimension(1:3,1:inel), intent(in) :: dxx !< current time step
  real(prec),dimension(1:3,1:inel), intent(out) :: vels !< current time step
  !------
  integer :: i, j, k, icell, it, ilog, nkall, nel

  real(prec),dimension(1:3,1:inel) :: dv_cos!< current time step
  ! real(prec),dimension(1:nkmod,1:inel) :: dcs!< current time step
  real(prec),dimension(1:nkmod,1:inel) :: dcos_v!< current time step
  real(prec),dimension(1:nkmod,1:inel) :: dsin_v!< current time step
  real(prec),dimension(1:nkmod,1:inel) :: dcs1 !< current time step
  ! real(prec),dimension(1:nmodes*ntimes*inel) :: dcos_v!< current time step
  ! real(prec),dimension(1:nmodes*ntimes*inel) :: dsin_v!< current time step
  real(prec) :: summ, dtmp

  ! 1 - set number of modes in all times
  nkall = nkmod
  nel = inel
  ! dcs(1:nkmod,1:inel) = 0.0d0

  do i = 1, nel
    dcs1(1:nkall,i) = c_m(1:nkall)
  enddo
  ! 2 - COS predicted vars
  ! cs_m(1:nkall) = c_m(1:nkall)

  ! 3 - Calc $ K \cdot x + w t$
  ! CALL dgemm('T','N',nkall,nel,3,1.0d0,b_m,3,dxx,nel,dtim_in,cs_m,nel)
  ! CALL dgemm('N','N',nel,nkall,3,1.0d0,dxx,nel,b_m,3,dtim_in,cs_m,nel)
  CALL dgemm('T','N',nkall,nel,3,1.0d0,b_m,3,dxx,3,dtim_in,dcs1,nkall)
  ! call gemm(b_m, dxx, cs_m ,transa='T',transb='N')
  ! write(*,'(a,8es13.5)')"A2: ", dcs(1:8,1)
  ! write(*,'(a,8es13.5)')"A2: ", cs_m(1:8)
  ! 4 - calc sin and cos
  CALL vdsincos(nkall*nel,dcs1(1,1),dsin_v(1,1),dcos_v(1,1))
  ! do i = 1, inel
  !   j = nkall * (i-1)
  !   CALL vdsincos(nkall,dcs(j),dsin_v(j),dcos_v(j))
  ! enddo
  ! 5 - Calc COS  like $AC \cdot \cos{ CS }$
  ! CALL dgemm('N','N',3,1,nkall,1.0d0,ac_m,3,cos_m,nkall,0.0d0,dv_cos,3)
  CALL dgemm('N','N',3,nel,nkall,1.0d0,ac_m,3,dcos_v,nkall,0.0d0,dv_cos,3)
  ! 6 - Calc SIN  like $AS \cdot \sin{ CS }$ and add it to prev val.
  CALL dgemm('N','N',3,nel,nkall,1.0d0,as_m,3,dsin_v,nkall,1.0d0,dv_cos,3)
  vels(1:3,1:inel) = dv_cos(1:3,1:inel)


  RETURN
end subroutine set_vels_at_space_time2
