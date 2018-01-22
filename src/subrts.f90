
!----------------------------------------------------------------------
!>@brief Allocate modes arrays
!!
subroutine tmp_alloc(icase)
  USE tmp_mod
  implicit none
  integer, optional :: icase
  integer :: i, ic

  ! ic = 4
  ! ! IF (PRESENT(icase)) ic = icase
  ! WRITE(*,*) "ICASE " , icase
  ! ic = icase

  ! IF (ic .eq. 5) THEN
  !   i = nmodes
  ! ELSE
  !   i = nmodes * ntimes  
  ! END IF

  i = nmodes
  ic = 0
  IF (ALLOCATED(c_m)) ic = SIZE(c_m)
  !
  ! 
  IF (ic.NE.i) THEN
    IF (ALLOCATED(ac_m)  ) DEALLOCATE(ac_m)  
    IF (ALLOCATED(as_m)  ) DEALLOCATE(as_m)  
    IF (ALLOCATED(b_m)   ) DEALLOCATE(b_m)   
    IF (ALLOCATED(c_m)   ) DEALLOCATE(c_m)   
    IF (ALLOCATED(dphi_m)) DEALLOCATE(dphi_m)
    IF (ALLOCATED(cs_m)  ) DEALLOCATE(cs_m)  
    IF (ALLOCATED(cos_m) ) DEALLOCATE(cos_m) 
    IF (ALLOCATED(sin_m) ) DEALLOCATE(sin_m) 
    IF (ALLOCATED(dtim)  ) DEALLOCATE(dtim)  
    IF (ALLOCATED(dRtang)) DEALLOCATE(dRtang)
    IF (ALLOCATED(dRlong)) DEALLOCATE(dRlong)
    IF (ALLOCATED(dRtime)) DEALLOCATE(dRtime)
  END IF !(ic.NE.i) THEN
  
  
  ALLOCATE(ac_m(1:3,1:i))
  ALLOCATE(as_m(1:3,1:i))
  ALLOCATE(b_m(1:3,1:i))
  ALLOCATE(c_m(1:i))
  ALLOCATE(dphi_m(1:3,1:i))
  ALLOCATE(cs_m(1:i))
  ALLOCATE(cos_m(1:i))
  ALLOCATE(sin_m(1:i))

  ac_m(:,:) = 0.0d0
  as_m(:,:) = 0.0d0
  b_m(:,:) = 0.0d0
  c_m(:) = 0.0d0
  dphi_m(:,:) = 0.0d0
  cs_m(:) = 0.0d0
  cos_m(:) = 0.0d0
  sin_m(:) = 0.0d0


  !-----
  ! Allocate timesteps
  ALLOCATE(dtim(1:ntimes))
  dtim(:) = 0.0d0


  !-----
  ! Allocate Vectors for correlations
  ALLOCATE(dRtime(0:ntimes))
  ALLOCATE(dRlong(0:int(tcell**0.333333)))
  ALLOCATE(dRtang(0:int(tcell**0.333333)))
  dRtime(:) = 0.0d0
  dRlong(:) = 0.0d0
  dRtang(:) = 0.0d0
  RETURN
end subroutine tmp_alloc

!----------------------------------------------------------------------
!>@brief Initialize the COMM1 module elements
!!
subroutine comm1_alloc()
  USE tmp_mod
  use comm1
    implicit none
    integer :: ikey
    integer :: i
!-----
!  Use temporary variables for coordinates and velocities in cells
  i = tcell * ntimes
  ALLOCATE(xp(1:3,1:tcell))
  ALLOCATE(u(1:3,1:i))

  xp(:,:) = 0.0d0
  u(:,:) = 0.0d0
    
end subroutine comm1_alloc


!>@brief set own seed for random generator
subroutine random_seed_user()
  IMPLICIT NONE
  ! ----- variables for portable seed setting -----
  INTEGER :: i_seed
  INTEGER, DIMENSION(:), ALLOCATABLE :: a_seed
  INTEGER, DIMENSION(1:8) :: dt_seed
  ! ----- end of variables for seed setting -----

  ! write(*,*) "Use the same RANDOM_SEED(1) to reproduce results"
  ! ! CALL RANDOM_SEED(1)
  ! RETURN

  ! ----- Set up random seed portably -----
  CALL RANDOM_SEED(size=i_seed)
  ALLOCATE(a_seed(1:i_seed))
  CALL RANDOM_SEED(get=a_seed)
  CALL DATE_AND_TIME(values=dt_seed)
  a_seed(i_seed)=dt_seed(8); a_seed(1)=dt_seed(8)*dt_seed(7)*dt_seed(6)
  CALL RANDOM_SEED(put=a_seed)
  DEALLOCATE(a_seed)
  ! ----- Done setting up random seed -----
end subroutine random_seed_user


!
!>@brief Random Sample from normal (Gaussian) distribution
!
FUNCTION rand_normal(mean,stdev) RESULT(c)
  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: PI=3.141592653589793238462d0
  DOUBLE PRECISION :: mean,stdev,c,temp(2), r, theta
  IF(stdev <= 0.0d0) THEN

    WRITE(*,*) "Standard Deviation must be +ve"
  ELSE
    CALL RANDOM_NUMBER(temp)
    r = (-2.0d0 * LOG(temp(1)))**5.0d-1
    theta = 2.0d0*PI*temp(2)
    c = mean + stdev * r * SIN(theta)
  END IF
END FUNCTION

!
!>@brief Random Sample from normal (Gaussian) distribution
!
SUBROUTINE rand_normal_sub(nel,mean,stdev,c)
  USE prec_mod
  USE constants, ONLY: f_zero, f_pi
  IMPLICIT NONE

  integer, intent(IN) :: nel !< number of elements
  real(prec), intent(IN) :: mean !< mean Value
  real(prec), intent(IN) :: stdev !< deviation
  real(prec), dimension(1:nel),intent(OUT) :: c !< return an array
  real(prec), dimension(1:nel*2) :: temp !< temporary
  real(prec) :: r, theta
  integer :: i
  IF(stdev <= 0.0d0) THEN

    WRITE(*,*) "Standard Deviation must be +ve"
  ELSE
    
    CALL random_seed_user()

    CALL RANDOM_NUMBER(temp(1:nel*2))
    
    do i= 1, nel
      r = (-2.0d0 * LOG(temp(i+i-1)))**5.0d-1
      theta = 2.0d0 * f_pi * temp(i+i)
      c(i) = mean + stdev * r * SIN(theta)
    enddo
  END IF
  RETURN
END SUBROUTINE rand_normal_sub

!----------------------------------------------------------------------
!>@brief Get correlations Vectors by time and space
!----------------------------------------------------------------------
subroutine get_cor()
  USE tmp_mod
  USE prec_mod
  USE comm1
  implicit none
  integer :: i, j, k, isc,iec,jsc,jec, itot, l, i1,j1
  real(prec) :: dtmp, rxx, ryy, rzz, rxy, rxz, ryx, ryz, rzx, rzy

  do k = 0, ntimes - 1
    dRtime(k) = 0.0d0
    ! write(*,'(6a5)') "i" , 'ntimes','isc','iec','jsc','jec'
    itot = ntimes - k
    do i= 1, itot
      isc = tcell * (i-1) + 1
      iec = tcell * (i)
      jsc = tcell * (i+k-1) + 1
      jec = tcell * (i+k)
      ! write(*,'(6i7)') i , ntimes,isc,iec,jsc,jec
      do j= 1, 3 !ntimes - i + 1
        dtmp = SUM(u(j,isc:iec)*u(j,jsc:jec)) 
        dtmp = dtmp / (3.0d0 * DBLE(itot) * DBLE(tcell))
        dRtime(k) = dRtime(k) + dtmp 
      enddo
    enddo
    ! if (k > 0) dRtime(k) = dRtime(k) / dRtime(0)
  end do !k = 1, ntimes
    
  do l = 0, min(nx,ny) - 1
    dRlong(l) = 0.0d0
    dRtang(l) = 0.0d0
    rxx = 0.0d0
    rxy = 0.0d0
    rxz = 0.0d0
    ryy = 0.0d0
    ryx = 0.0d0
    ryz = 0.0d0
    rzz = 0.0d0
    rzx = 0.0d0
    rzy = 0.0d0
    do i = 1, nx
      do j = 1, ny
        do k = 1, nz
          i1 = (i-1) * ny * nz + nz * (j-1) + k
          if (l + i .LT. nx) THEN
            j1 = (i-1+l) * ny * nz + nz * (j-1) + k
            rxx = rxx + SUM(u(1,i1::tcell)*u(1,j1::tcell))/DBLE(ntimes)
            rxy = rxy + SUM(u(2,i1::tcell)*u(2,j1::tcell))/DBLE(ntimes)
            rxz = rxz + SUM(u(3,i1::tcell)*u(3,j1::tcell))/DBLE(ntimes)
          end if !(k1 + i .LE. nx) THEN
          if (l + j .LT. ny) THEN
            j1 = (i-1) * ny * nz + nz * (j-1+l) + k
            ryy = ryy + SUM(u(2,i1::tcell)*u(2,j1::tcell))/DBLE(ntimes)
            ryx = ryx + SUM(u(1,i1::tcell)*u(1,j1::tcell))/DBLE(ntimes)
            ryz = ryz + SUM(u(3,i1::tcell)*u(3,j1::tcell))/DBLE(ntimes)
          end if !(k1 + i .LE. nx) THEN
          if (l + k .LT. nz) THEN
            j1 = (i-1) * ny * nz + nz * (j-1) + k+l
            rzz = rzz + SUM(u(3,i1::tcell)*u(3,j1::tcell))/DBLE(ntimes)
            rzx = rzx + SUM(u(1,i1::tcell)*u(1,j1::tcell))/DBLE(ntimes)
            rzy = rzy + SUM(u(2,i1::tcell)*u(2,j1::tcell))/DBLE(ntimes)
          end if !(k1 + i .LE. nx) THEN
        enddo
      enddo
    enddo
    dRlong(l) = (rxx / DBLE(nx-l)+ryy / DBLE(ny-l)+rzz / DBLE(nz-l)) / 3.0d0
    dRtang(l) = dRtang(l) + (rxy + rxz) / DBLE(nx-l) / 6.0d0
    dRtang(l) = dRtang(l) + (ryx + ryz) / DBLE(ny-l) / 6.0d0
    dRtang(l) = dRtang(l) + (rzx + rzy) / DBLE(nz-l) / 6.0d0
  enddo !k1 = 0, min(nx,ny)

    

end subroutine get_cor


subroutine read_coef()
  USE tmp_mod
  implicit none
  integer :: ifile
  logical :: lop
  CHARACTER(len=128) :: sname
  CHARACTER(len=256) :: txt1,txt2
  integer :: i,j

  ifile  = 212
  INQUIRE(ifile,OPENED= lop)
  IF(.NOT. lop) THEN
    OPEN(ifile,FILE=TRIM(sname))
    ! READ(ifile,*) sname
    READ(ifile,*) nmodes
  END IF ! itst.eq.1

  ! Reallocate variables
  CALL tmp_alloc()
  !
  ! Read all coefficients from file
  DO i = 1, nmodes
    read(ifile,*) j, ac_m(1:3,i), as_m(1:3,i), b_m(1:3,i), c_m(i)
  END DO !i = 1, nmodes
  CLOSE(ifile)
    
end subroutine read_coef


!----------------------------------------------------------------------
!>@brief Write config files 
subroutine write_coef()
  USE tmp_mod
  implicit none
  integer :: ifile
  logical :: lop
  CHARACTER(len=128) :: sname
  CHARACTER(len=256) :: txt1,txt2
  integer :: i,j

  ifile  = 212
  INQUIRE(ifile,OPENED= lop)
  IF(.NOT. lop) THEN
    WRITE(sname,'(2a)') 'user_coef','.dat'
    OPEN(ifile,FILE=TRIM(sname))
    ! WRITE(ifile,'(2a)') 'mode  P_cos_1 P_cos_2 P_cos_3 Q_cos_1', &
    ! ' Q_cos_2 Q_cos_3 K_1 K_2 K_3 time_freq'
  END IF ! itst.eq.1
  write(ifile,'(i7)') nmodes
  DO i = 1, nmodes
    write(ifile,'(i6,10es17.9)')i, ac_m(1:3,i), as_m(1:3,i), b_m(1:3,i), c_m(i)
  END DO !i = 1, nmodes
  CLOSE(ifile)
    
end subroutine write_coef