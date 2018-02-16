!>@brief  main porgram algorithm 
!!
!>@todo 
!!- Make time dependence
!!- Extend for 2D space 
!!- Integrate to FIRE via `useini`
program gen_flow_saad
  USE prec_mod
  USE tmp_mod
  USE comm0, ONLY: ncell, nsp,nep
  USE comm1, ONLY: u, xp
  implicit none
  real(prec) :: dlt !< Time duration
  real(prec) :: dlx !< Length by X-corrdinate
  real(prec) :: dly !< Length by Y-corrdinate
  real(prec) :: dlz !< Length by Z-corrdinate
  integer :: nt !< number of elements by time
  ! integer :: nx !< number of elements by X-corrdinate
  ! integer :: ny !< number of elements by Y-corrdinate
  ! integer :: nz !< number of elements by Z-corrdinate

  ! integer :: nmodes !< number of Modes

  real(prec) :: dsigma    !< define Deviation of velocities
  real(prec) :: dlength   !< define Integral Length Scale 
  real(prec) :: dtau      !< define Integral Time Scale 
  real(prec) :: std       !< @var define Integral Time Scale 
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
  integer :: icase = 5                       !< Set the case, which we want to run

  interface
  subroutine tmp_alloc(icase)
  USE tmp_mod
  implicit none
  integer, optional :: icase
  end subroutine tmp_alloc

  subroutine set_vels_time_space(icase)
    USE prec_mod
    USE comm1, ONLY: xp, u
    USE tmp_mod
    implicit none
    integer, OPTIONAL :: icase
  end subroutine set_vels_time_space

  end interface
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
  dlt = 1.0d-6
  ! set space Length
  dlx = 2.0d-3
  dly = dlx
  dlz = dlx
  ! set number of nodes
  nx = 40
  ny = nx
  nz = nx

  !set number of timesteps
  nt = 10
  ntimes = nt

  !set number of Modes
  nmodes = 10000
  tcell  = nx * ny * nz
  CALL tmp_alloc()

  !set VALs from FIRE
  ncell = tcell
  nsp(1) = 1
  nsp(1) = ncell
  CALL comm1_alloc()

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
  dlength = 1.0d-2
  dsigma = 1.0d+1
  dtau = dlength / dsigma

  !generate arrays 
  dels(1) = dlx * 2.5d+1
  dels(2) = dly * 2.5d+1
  dels(3) = dlz * 2.5d+1

  nels(1) = nx * 25
  nels(2) = ny * 25
  nels(3) = nz * 25


  !-----------------------------------------------------------------

  ! write(*,*) "Write out time correlation"
  ! OPEN(121,FILE='store.dat',ACTION="READ")
  ! do i=1,tcell * ntimes
  !   read(121,'(3es13.5)') u(1:3,i)
  ! enddo
  ! CALL get_cor()
  ! write(*,'(4es13.5)') 0.0d0, dRtime(0), 0.0d0,dRlong(0),dRtang(0)
  ! do i = 1, ntimes
  ! if (i .lt. nx)write(*,'(5es13.5)') dtim(i), dRtime(i),dx_i(1,i),dRlong(i),dRtang(l)
  ! if (i .ge. nx)write(*,'(4es13.5)') dtim(i), dRtime(i),0.0d0,0.0d0
  ! enddo

  ! STOP
  !-----------------------------------------------------------------

  
  write(*,*) "Generate  the modes"
  if (icase == 3) THEN
    !-----
    ! Generate modes accoring to Kraichnan article
    !
    ! generate fields
    CALL gen_flow_3d(dels, nels, dsigma, dlength, dtau)


    write(*,*) "Calculate the velocities"
    CALL set_vels()

  ELSE IF (icase == 4) THEN
    !-----
    ! Generate modes accoring to time and space modes
    !
    write(*,*) "work in 3D-space + Time"

    do i = 1, ntimes
      in_time = nmodes * (i-1) + 1
    ! generate fields
      CALL gen_flow_3d(dels, nels, dsigma, dlength, dtau)
    end do

    write(*,*) "Calculate the velocities"
    CALL set_vels_time_space(icase)


  ELSE IF (icase == 5) THEN
    !-----
    ! Generate modes accoring to Kraichnan article for different times
    !
    write(*,*) "work in 3D-space + Time. Generate only modes"

    ! generate fields
    CALL gen_flow_3d(dels, nels, dsigma, dlength, dtau)

    write(*,*) "Calculate the velocities"
    CALL set_vels_time_space(icase)

  END IF

  write(*,*) "Calcs mean and deviation"
  do i=1, 3
    dmean_i(i) = SUM(u(i,:)) / DBLE(tcell * ntimes)
    std_i(i) = SQRT(SUM((u(i,:)-dmean_i(i))**2) / DBLE(tcell*ntimes))
    write(*,'(i13,2es13.5)') i, dmean_i(i),std_i(i)
  enddo

  write(*,*) "Write out time correlation"
  CALL get_cor()
  ion = 121
  OPEN(ion,file='store_corr.dat')
  write(ion,'(5es13.5)') 0.0d0, dRtime(0), 0.0d0,dRlong(0),dRtang(0)
  do i = 1, MAX(ntimes,nx)-1
    if (i .lt. MIN(nx,ntimes))write(ion,'(5es13.5)') dtim(i), dRtime(i),dx_i(1,i),dRlong(i),dRtang(i)
    if (i.ge.nx .and. i.lt.ntimes)write(ion,'(5es13.5)') dtim(i), dRtime(i),0.0d0,0.0d0,0.0d0
    if (i.lt.nx .and. i.ge.ntimes)write(ion,'(5es13.5)') 0.0d0, 0.0d0,dx_i(1,i),dRlong(i),dRtang(i)
  enddo
  CLOSE(ion)


  write(*,*) "Store the array"
  !store the field
  ion = 121
  OPEN(ion,file='store.dat')
  ! write(ion,'(3i5)') nx, ny, nz
  do i=1, tcell * ntimes
    write(ion,'(3es13.5)') u(1:3,i)
  enddo
  CLOSE(ion)
  write(*,*) "Exit the program"

end program gen_flow_saad
!----------------------------------------------------------------------


