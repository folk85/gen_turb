!----------------------------------------------------------------------
!@brief Generate coefficients for modes using current geometry, 
!! mesh resolution and number of required wave-modes
!!
subroutine gen_flow_3d(dls,nels,dsigma,dlength,dtau)
  USE prec_mod
  USE constants, ONLY: f_zero, f_pi
  USE tmp_mod
  implicit none
  integer, dimension(1:3), intent(IN) :: nels   !< number of modes
  real(prec), dimension(3), intent(IN) :: dls !< Physical boundaries
  real(prec), intent(IN) :: dsigma    !< define Deviation of velocities
  real(prec), intent(IN) :: dlength   !< define Integral Length Scale 
  real(prec), intent(IN) :: dtau      !< define Integral Time Scale 
  integer :: ix !< number of elements by X-corrdinate
  integer :: iy !< number of elements by Y-corrdinate
  integer :: iz !< number of elements by Z-corrdinate
  real(prec) :: dx      !< cell width by X-corrdinate
  real(prec) :: dy      !< cell width by Y-corrdinate
  real(prec) :: dz      !< cell width by Z-corrdinate
  real(prec),allocatable :: dkm(:)      !< wave Length of m-th mode
  real(prec) :: kmin    !< minimal wave Length
  real(prec) :: kmax    !< maximal wave Length
  real(prec) :: kwid    !< interval between wave modes
  real(prec),allocatable :: dqm(:)      !< amplitude of 

  ! real(prec) :: set_eturb     !< external function
  ! real(prec) ::set_unit_vector !< external function
  ! real(prec) ::cross_product !< external function
  ! external :: set_eturb
  ! external :: set_unit_vector
  ! external :: cross_product

  real(prec),dimension(1:3) :: tmp3     !< temporary vector 3 - elements
  real(prec),allocatable :: vtmp(:)     !< temporary vector
  ! real(prec),allocatable :: vtmp2(:,:)     !< temporary vector
  real(prec),allocatable :: dkun_i(:,:)     !< Unit vectors of wave number
  real(prec),allocatable :: dk_i(:,:)     !< wave number vector -direc
  real(prec),allocatable :: dsim_i(:,:)     !< unit direction vector in COS
  real(prec),allocatable :: dksi_i(:,:)     !< unit direction vector in SIN
  real(prec),allocatable :: dpsi(:)     !< random value in [0,2Pi)

  real(prec) :: dw_freq !< wave-time frequency multiplier
  real(prec) :: dtmp 

  integer :: i  !< temporary index
  integer :: j  !< temporary index
  integer :: ist  !< Starting index
  integer :: ien  !< Ending index
  real(prec) :: rand_normal
  ! real(prec) :: dtmp1
  INTERFACE
    function set_eturb(dk, dl_in, dsigma_in) result(de)
      USE prec_mod
      implicit none
      real(prec), intent(IN) :: dk !< wave number
      real(prec), intent(IN), optional :: dl_in !< Integral Length Scale (set by default 0.1 m)
      real(prec), intent(IN), optional :: dsigma_in !< RMS of velocities (set by default 10 m/s)
      real(prec) :: de !< return value
   end function set_eturb
  end INTERFACE
  !generate random seed
  ! CALL CPU_TIME(time)
  ! CALL RANDOM_SEED(PUT=time)
  CALL random_seed_user()
  
  ! allocatable wave number array
  if(.NOT.ALLOCATED(dkm)) ALLOCATE(dkm(1:nmodes))
  dkm(:) = 0.0d0

  ! allocatable amplitude array
  if(.NOT.ALLOCATED(dqm)) ALLOCATE(dqm(1:nmodes))
  dqm(:) = 0.0d0

  ! Allocate unit vectors of wave number
  if(.NOT.ALLOCATED(dkun_i)) ALLOCATE(dkun_i(1:3,1:nmodes))
  dkun_i(:,:) = 0.0d0

  ! allocatable wave number vector
  if(.NOT.ALLOCATED(dk_i)) ALLOCATE(dk_i(1:3,1:nmodes))
  dk_i(:,:) = 0.0d0

  ! allocatable wave number vector in COS
  if(.NOT.ALLOCATED(dsim_i)) ALLOCATE(dsim_i(1:3,1:nmodes))
  dsim_i(:,:) = 0.0d0

  ! allocatable wave number vector in SIN
  if(.NOT.ALLOCATED(dksi_i)) ALLOCATE(dksi_i(1:3,1:nmodes))
  dksi_i(:,:) = 0.0d0

  ! allocatable random value
  if(.NOT.ALLOCATED(dpsi)) ALLOCATE(dpsi(1:nmodes))
  dpsi(:) = 0.0d0

  ix = nels(1)
  iy = nels(2)
  iz = nels(3)

  ! set cell width
  dx = dls(1) / REAL(nels(1),prec)
  dy = dls(2) / REAL(nels(2),prec)
  dz = dls(3) / REAL(nels(3),prec)

  ! 2-3 -  define minimal wave Length
  kmin = 0.0d0
  kmax = 0.0d0
  do i=1, 3
    kmin = MAX(kmin,2.0d0 * f_pi / dls(i))
    kmax = MAX(kmax,2.0d0*f_pi * DBLE(nels(i)) / dls(i))
  enddo

  ! 4 - Generate a list of M equidistant modes:
  !   $$k_m = k_o + \frac{k_{max}-k_o}{M} (m - 1) $$
  kwid = (kmax - kmin) / DBLE(nmodes)
  ! set waves for modes
  do i = 1, nmodes
    dkm(i) = kmin + kwid * DBLE(i - 1)
    !5 - generate amplitudes by equation $q_m = \sqrt{E(k_m)\Delta k}$
    !       Get SQRT below
    dqm(i) = set_eturb(dkm(i),dl_in=dlength,dsigma_in=dsigma) * kwid
  enddo

  ! 6 - 7 Calculate wave number vector
  ! ALLOCATE(vtmp(1:3,1:nmodes))
  if(.NOT.ALLOCATED(vtmp)) ALLOCATE(vtmp(1:nmodes*4))
  
  CALL RANDOM_NUMBER(vtmp(:))

  do i = 1, nmodes
    ! generate unit vector in 3D space
    ! tmp3(1:3) = set_unit_vector(vtmp(1+3*(i-1):3*i))
    j = 3 * i - 2
    ! CALL set_unit_vector_sub(vtmp(j:j+1),tmp3(1:3))
    CALL set_unit_vector_sub(vtmp(j:j+2),tmp3(1:3))
    dkun_i(1:3,i) = tmp3(1:3) !set_unit_vector(vtmp(1:2,i))
    dk_i(1:3,i) = tmp3(1:3) !set_unit_vector(vtmp(1:2,i))
    ! dk_i(1,i) = 2.0d0 * SIN(5.0d-1 * dx * dkm(i) * tmp3(1)) / dx
    ! dk_i(2,i) = 2.0d0 * SIN(5.0d-1 * dy * dkm(i) * tmp3(2)) / dy
    ! dk_i(3,i) = 2.0d0 * SIN(5.0d-1 * dz * dkm(i) * tmp3(3)) / dz
    !----------
    ! Set another distribution
    dk_i(1,i) = 1.0d0 * SIN(1.0d-0 * dx * dkm(i) * tmp3(1)) / dx
    dk_i(2,i) = 1.0d0 * SIN(1.0d-0 * dy * dkm(i) * tmp3(2)) / dy
    dk_i(3,i) = 1.0d0 * SIN(1.0d-0 * dz * dkm(i) * tmp3(3)) / dz
  enddo

  ! 8 - Define random unity vectors
  CALL RANDOM_NUMBER(vtmp(1:nmodes*4))

 ! 9 - generate (unit) direction vector
  do i= 1, nmodes
    CALL set_unit_vector_sub(vtmp(4*i-3),tmp3(1:3))
    !-----------------------------
    ! Set addtitional conditions for Periodic BND
    !-----------------------------
    ! tmp3(1) = SIN(5.0d-1 * (dls(1) - dx) * dkm(i) * dkun_i(1,i))
    ! tmp3(2) = SIN(5.0d-1 * (dls(2) - dy) * dkm(i) * dkun_i(2,i))
    ! tmp3(3) = SIN(5.0d-1 * (dls(3) - dz) * dkm(i) * dkun_i(3,i))
    ! tmp3(2) = SIGN(1.0d0,vtmp(4*i-2)-5.0d-1) * SQRT(ABS(1.0d0-tmp3(1)**2-tmp3(3)**2))
    !!!!!!!
    !
    CALL cross_product_sub(tmp3(1:3), dk_i(1:3,i),dsim_i(1:3,i))

    ! 9b - set for second unit vector in SIN coefficients
    CALL set_unit_vector_sub(vtmp(4*i-1),tmp3(1:3))
    CALL cross_product_sub(tmp3(1:3), dk_i(1:3,i),dksi_i(1:3,i))
    ! CALL cross_product_sub(dsim_i(1:3,i), dk_i(1:3,i),dksi_i(1:3,i))
  enddo

  ! 10 - Generate random value
  ! CALL RANDOM_NUMBER(dpsi(1:nmodes))
  ! dw_freq = 1.0d0
  ! dpsi(:) = dpsi(:) * 2.0d0 * f_pi * dw_freq
  !! According to @cite{juve_1999.pdf} use Gauss distribution function
  do i=1, nmodes
    dtmp = dkm(i) * dsigma
    dpsi(i) = rand_normal(dtmp,dtmp)
  enddo

  ! 11 - Calculate velocities in every point

  ! write(*,*) "Generate Spectrum profile"
  ! OPEN(UNIT=123,FILE='tests/spectr.dat')
  ! do i= 1, nmodes
  !   ! dtmp = 10. + 10.*(i-1)
  !   ! dtmp1 = set_eturb(dtmp,dlength,dsigma)
  !   ! write(*,'(2es13.5)')dtmp, dtmp1
  !   ! if (dtmp /= dtmp) write(*,*) "Some errors"
  !   ! write(*,'(3es13.5)') dkm(i),set_eturb(dkm(i),dlength,dsigma), dqm(i)
  !   write(123,'(3es13.5)') dkm(i),set_eturb(dkm(i),dl_in=dlength,dsigma_in=dsigma), dqm(i)
  ! enddo
  ! CLOSE(123)

  !fill the coefficients in from tmp_mod
  ist = in_time
  ien = in_time + nmodes - 1

  do i = 1, 3
    CALL rand_normal_sub(nmodes,0.0d0,1.0d0,vtmp(1:nmodes))
    ! ac_m(i,ist:ien) = dsim_i(i,:) * SQRT(dqm(:))
    ac_m(i,ist:ien) = dsim_i(i,:) * SQRT(dqm(:) * vtmp(1:nmodes)**2 * 1.0d0  &
      /(1.0d0+vtmp(1:nmodes)**2)/DBLE(ntimes))* SIGN(1.0d0,vtmp(1:nmodes))
    as_m(i,ist:ien) = dksi_i(i,:) * SQRT(dqm(:) *1.0d0                       &
      /(1.0d0+vtmp(1:nmodes)**2)/DBLE(ntimes))* SIGN(1.0d0,vtmp(1:nmodes))
    b_m(i,ist:ien)  = dkun_i(i,:) * dkm(:)
  enddo

  ! WRITE(*,*) "Switch Time coefficiet in `gen_flow_3d` !!!"

  ! do i= 1, nmodes
  !   j = i + ist - 1
  !   c_m(j)  = dpsi(i) * dsigma * dkm(i)
  ! enddo
  c_m(ist:ien) = dpsi(:) !* dkm(:) * dsigma /3.0d0/4.0d0
  ! c_m(:)  = dpsi(:)

  write(*,*) "END: Generate Spectrum profile"

  ! CALL RANDOM_NUMBER(vels(1:3,1:ix,1:iy,1:iz))

  RETURN
end subroutine gen_flow_3d
!----------------------------------------------------------------------
