!@brief  main porgram algorithm 
!!
!@todo 
!!- Make time dependence
!!- Extend for 2D space 
!!- Integrate to FIRE via `useini`
program gen_flow_saad
  USE prec_mod
  USE tmp_mod
  implicit none
  real(prec) :: dlt !< Time duration
  real(prec) :: dlx !< Length by X-corrdinate
  real(prec) :: dly !< Length by Y-corrdinate
  real(prec) :: dlz !< Length by Z-corrdinate
  integer :: nt !< number of elements by time
  integer :: nx !< number of elements by X-corrdinate
  integer :: ny !< number of elements by Y-corrdinate
  integer :: nz !< number of elements by Z-corrdinate

  ! integer :: nmodes !< number of Modes

  real(prec) :: dsigma    !< define Deviation of velocities
  real(prec) :: dlength   !< define Integral Length Scale 
  real(prec) :: dtau      !< define Integral Time Scale 
  real(prec) :: std      !< define Integral Time Scale 
  real(prec) :: dmean      !< define Integral Time Scale 
!-----
  integer :: i, i1
  integer :: j
  integer :: k
  integer :: ion                            !< number of unitfor I/O
  integer, dimension(3)  :: nels            !< 
  real(prec), dimension(3)  :: dels         !<
  real(prec), dimension(3)  :: tmp3         !< temporary
  real(prec), dimension(:,:),allocatable  :: dx_i         !< temporary
  real(prec), dimension(3)  :: std_i         !<
  real(prec), dimension(3)  :: dmean_i         !<
  real(prec), dimension(:), allocatable  :: dtmp   !< temporary array
  logical :: ltest = .FALSE.   !< Set test run  .TRUE.
  real(prec) :: rand_normal !< RNG function
!-----
  integer :: icase = 4                       !< Set the case, which we want to run
  ! real(prec), allocatable :: u(:,:,:,:)

  ! write(*,*) "Check cross_product _ initial"

  ! CALL check_cross()

  ! STOP
  if (ltest) then
    write(*,*) "Generate Random Numbers by Gauss distribution"
    k = 1000
    ALLOCATE(dtmp(1:k))
    do i = 100, 1000, 100
      CALL rand_normal_sub(i,0.0d0,1.0d0,dtmp(1:i))
      dmean = SUM(dtmp(1:i)) / DBLE(i)
      std = SQRT(SUM((dtmp(1:i)-dmean)**2)/ DBLE(i))
      write(*,"(i13,2es13.5)") i, dmean, std
    END DO
    OPEN(UNIT=123,FILE='tests/gauss.dat')
    do i= 1, k
      write(123,'(3es13.5)') dtmp(i), rand_normal(0.0d0,1.0d0)
    enddo
    CLOSE(123)

    STOP

  endif

  write(*,*) "Welcome into th program"
  ! set  time duration
  dlt = 1.0d-2
  ! set space Length
  dlx = 1.0d-1
  dly = dlx
  dlz = dlx
  ! set number of nodes
  nx = 64
  ny = nx
  nz = nx

  !set number of timesteps
  nt = 64
  ntimes = nt

  !set number of Modes
  nmodes = 100
  ncell  = nx * ny * nz
  CALL tmp_alloc()

  !define coordinates
  ALLOCATE(dx_i(1:3,1:nx))
  do i = 1, nx
    dx_i(1,i) = dlx * DBLE(2*i-1) * 5.d-1 / DBLE(nx)
    dx_i(2,i) = dly * DBLE(2*i-1) * 5.d-1 / DBLE(nx)
    dx_i(3,i) = dlz * DBLE(2*i-1) * 5.d-1 / DBLE(nx)
  end do
  do i = 1, nx
    do j = 1, ny
      do k = 1, nz
        i1 = (i-1) * ny * nz + nz * (j-1) + k
        xp(1,i1) = dx_i(1,i)
        xp(2,i1) = dx_i(2,j)
        xp(3,i1) = dx_i(3,k)
      enddo
    enddo
  enddo

  !define timesteps
  do i= 1, ntimes
    dtim(i) = DBLE(i) * dlt / DBLE(ntimes)
  enddo
  ! set starting index for the coefficients in modes
  in_time = 1

  !set the Integral values
  dlength = 1.0d-1
  dsigma = 1.0d+1
  dtau = dlength / dsigma

  !generate arrays 
  dels(1) = dlx
  dels(2) = dly
  dels(3) = dlz

  nels(1) = nx
  nels(2) = ny
  nels(3) = nz

  
  write(*,*) "Generate  the modes"
  if (icase == 3) THEN
    ! generate fields
    CALL gen_flow_3d(dels, nels, dsigma, dlength, dtau)


    write(*,*) "Calculate the velocities"
    CALL set_vels()

    do i=1, 3
      dmean_i(i) = SUM(u(i,:)) / DBLE(ncell)
      std_i(i) = SQRT(SUM((u(i,:)-dmean_i(i))**2) / DBLE(ncell))
      write(*,'(i13,2es13.5)') i, dmean_i(i),std_i(i)
    enddo
  ELSE IF (icase == 4) THEN
    write(*,*) "work in 3D-space + Time"

    do i = 1, ntimes
      in_time = ntimes * (i-1) + 1
    ! generate fields
      CALL gen_flow_3d(dels, nels, dsigma, dlength, dtau)
    end do

    write(*,*) "Calculate the velocities"
    CALL set_vels_time_space()

    write(*,*) "Calcs mean and deviation"
    do i=1, 3
      dmean_i(i) = SUM(u(i,:)) / DBLE(ncell * ntimes)
      std_i(i) = SQRT(SUM((u(i,:)-dmean_i(i))**2) / DBLE(ncell*ntimes))
      write(*,'(i13,2es13.5)') i, dmean_i(i),std_i(i)
    enddo
  END IF

  write(*,*) "Store the array"
  !store the field
  ion = 121
  OPEN(ion,file='store.dat')
  ! write(ion,'(3i5)') nx, ny, nz
  do i=1, ncell * ntimes
    write(ion,'(3es13.5)') u(1:3,i)
  enddo
  CLOSE(ion)
  write(*,*) "Exit the program"

end program gen_flow_saad
!----------------------------------------------------------------------


!----------------------------------------------------------------------
!@brief Allocate modes arrays
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
  i = ncell * ntimes
  ALLOCATE(xp(1:3,1:ncell))
  ALLOCATE(u(1:3,1:i))

  xp(:,:) = 0.0d0
  u(:,:) = 0.0d0

!-----
! Allocate timesteps
  ALLOCATE(dtim(1:ntimes))
  dtim(:) = 0.0d0

  RETURN
end subroutine tmp_alloc

!@brief set own seed for random generator
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
!@brief Random Sample from normal (Gaussian) distribution
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
!@brief Random Sample from normal (Gaussian) distribution
!
SUBROUTINE rand_normal_sub(nel,mean,stdev,c)
  USE prec_mod
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


