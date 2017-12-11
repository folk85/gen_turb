
!----------------------------------------------------------------------
!>@brief Allocate modes arrays
!!
subroutine tmp_alloc()
  USE tmp_mod
  implicit none
  integer :: i

  i = nmodes * ntimes  
  ALLOCATE(ac_m(1:3,1:i))
  ALLOCATE(as_m(1:3,1:i))
  ALLOCATE(b_m(1:3,1:i))
  ALLOCATE(c_m(1:i))
  ALLOCATE(dphi_m(1:3,1:i))
  ALLOCATE(cs_m(1:i))

  ac_m(:,:) = 0.0d0
  as_m(:,:) = 0.0d0
  b_m(:,:) = 0.0d0
  c_m(:) = 0.0d0
  dphi_m(:,:) = 0.0d0
  cs_m(:) = 0.0d0

!-----
!  Use temporary variables for coordinates and velocities in cells
  ! i = tcell * ntimes
  ! ALLOCATE(xp(1:3,1:tcell))
  ! ALLOCATE(u(1:3,1:i))

  ! xp(:,:) = 0.0d0
  ! u(:,:) = 0.0d0

!-----
! Allocate timesteps
  ALLOCATE(dtim(1:ntimes))
  dtim(:) = 0.0d0

  RETURN
end subroutine tmp_alloc

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

