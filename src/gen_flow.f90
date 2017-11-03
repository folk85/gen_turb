!@brief  main porgram algorithm 
!!
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

  integer :: nmodes !< number of Modes

  real(prec) :: dsigma    !< define Deviation of velocities
  real(prec) :: dlength   !< define Integral Length Scale 
  real(prec) :: dtau      !< define Integral Time Scale 
!-----
  integer :: i
  integer :: j
  integer :: k
  integer :: ion                            !< number of unitfor I/O
  integer, dimension(3)  :: nels            !< 
  real(prec), dimension(3)  :: dels         !<
  real(prec), allocatable :: u(:,:,:,:)

  ! write(*,*) "Check cross_product _ initial"

  ! CALL check_cross()

  ! STOP

  write(*,*) "Welcome into th program"
  ! set  time duration
  dlt = 0.0d0
  ! set space Length
  dlx = 2.0d0 * f_pi * 1.0d-1
  dly = dlx
  dlz = dlx
  ! set number of nodes
  nx = 64
  ny = nx
  nz = nx

  !set number of Modes
  nmodes = 1000
  CALL tmp_alloc(nmodes)

  !set timesteps
  nt = 1

  !set the Integral values
  dlength = 7.0d-2
  dsigma = 1.0d+1
  dtau = dlength / dsigma

  !generate arrays 
  dels(1) = dlx
  dels(2) = dly
  dels(3) = dlz

  nels(1) = nx
  nels(2) = ny
  nels(3) = nz

  ! allocate the velocities
  ALLOCATE(u(1:3,1:nx,1:ny,1:nz))
  u(:,:,:,:) = 0.0d0

  write(*,*) "Generate  the flow"
  ! generate fields
  CALL gen_flow_3d(dels, nels,nmodes, dsigma,dlength,dtau)

  write(*,*) "Store the array"
  !store the field
  ion = 121
  OPEN(ion,file='store.dat')
  write(ion,'(3i5)') nx, ny, nz
  do i=1, nz
    do j=1, ny
      do k=1, nx
        write(ion,'(3es13.5)') u(1:3,k,j,i)
      enddo
    enddo
  enddo
  CLOSE(ion)
  write(*,*) "Exit the program"

end program gen_flow_saad
!----------------------------------------------------------------------



!----------------------------------------------------------------------
!@brief Generate coefficients for modes using current geometry, 
!! mesh resolution and number of required wave-modes
!!
subroutine gen_flow_3d(dls,nels,nms, dsigma,dlength,dtau)
  USE prec_mod
  USE tmp_mod
  implicit none
  integer, dimension(1:3), intent(IN) :: nels   !< number of modes
  real(prec), dimension(3), intent(IN) :: dls !< Physical boundaries
  integer, intent(IN) :: nms !< number of modes
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

  real(prec) :: set_eturb     !< external function
  ! real(prec) ::set_unit_vector !< external function
  ! real(prec) ::cross_product !< external function
  ! external :: set_eturb
  ! external :: set_unit_vector
  ! external :: cross_product

  real(prec),dimension(1:3) :: tmp3     !< temporary vector 3 - elements
  real(prec),allocatable :: vtmp(:)     !< temporary vector
  ! real(prec),allocatable :: vtmp2(:,:)     !< temporary vector
  real(prec),allocatable :: dk_i(:,:)     !< wave number vector -direc
  real(prec),allocatable :: dsim_i(:,:)     !< unit direction vector
  real(prec),allocatable :: dpsi(:)     !< random value in [0,2Pi)


  integer :: i  !< temporary index
  integer :: j  !< temporary index
  ! real(prec) :: dtmp
  ! real(prec) :: dtmp1

  ! allocatable wave number array
  if(ALLOCATED(dkm)) DEALLOCATE(dkm)
  ALLOCATE(dkm(1:nms))
  dkm(:) = 0.0d0

  ! allocatable amplitude array
  if(ALLOCATED(dqm)) DEALLOCATE(dqm)
  ALLOCATE(dqm(1:nms))
  dqm(:) = 0.0d0

  ! allocatable wave number vector
  if(ALLOCATED(dk_i)) DEALLOCATE(dk_i)
  ALLOCATE(dk_i(1:3,1:nms))
  dk_i(:,:) = 0.0d0

  ! allocatable wave number vector
  if(ALLOCATED(dsim_i)) DEALLOCATE(dsim_i)
  ALLOCATE(dsim_i(1:3,1:nms))
  dsim_i(:,:) = 0.0d0

  ! allocatable random value
  if(ALLOCATED(dpsi)) DEALLOCATE(dpsi)
  ALLOCATE(dpsi(1:nms))
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

  ! 4 - define inverals and 
  kwid = (kmax - kmin) / DBLE(nms)
  ! set waves for modes
  do i = 1, nms
    dkm(i) = kmin + kwid * DBLE(i - 1)
    !5 - generate amplitudes by equation $q_m = \sqrt{E(k_m)\Delta k}$
    dqm(i) = SQRT(set_eturb(dkm(i),dlength,dsigma) * kwid)
  enddo

  ! 6 - 7 Calculate wave number vector
  ! ALLOCATE(vtmp(1:3,1:nms))
  ALLOCATE(vtmp(1:nms*3))
  CALL RANDOM_NUMBER(vtmp(:))
  do i = 1, nms
    ! generate unit vector in 3D space
    ! tmp3(1:3) = set_unit_vector(vtmp(1+3*(i-1):3*i))
    j = 3 * i - 2
    ! CALL set_unit_vector_sub(vtmp(j:j+1),tmp3(1:3))
    CALL set_unit_vector_sub(vtmp(j:j+2),tmp3(1:3))
    dk_i(1:3,i) = tmp3(1:3) !set_unit_vector(vtmp(1:2,i))
    dk_i(1,i) = 2.0d0 * SIN(5.0d-1 * dx * dkm(i) * tmp3(1)) / dx
    dk_i(2,i) = 2.0d0 * SIN(5.0d-1 * dy * dkm(i) * tmp3(2)) / dy
    dk_i(3,i) = 2.0d0 * SIN(5.0d-1 * dz * dkm(i) * tmp3(3)) / dz
  enddo

  ! 8 - Define random unity vectors
  CALL RANDOM_NUMBER(vtmp(1:nms*3))

 ! 9 - generate (unit) direction vector
  do i= 1, nms
    ! tmp3(1:3) = set_unit_vector(vtmp(1+3*(i-1):3*i))
    CALL set_unit_vector_sub(vtmp(1+3*(i-1):3*i),tmp3(1:3))
    ! dsim_i(1:3,i) = cross_product(tmp3(1:3), dk_i(1:3,i))
    CALL cross_product_sub(tmp3(1:3), dk_i(1:3,i),dsim_i(1:3,i))
  enddo

  ! 10 - Generate random value
  CALL RANDOM_NUMBER(dpsi(1:nms))
  dpsi(:) = dpsi(:) * 2.0d0 * f_pi

  ! 11 - Calculate velocities in every point

  write(*,*) "Generate Spectrum profile"
  OPEN(UNIT=123,FILE='tests/spectr.dat')
  do i= 1, nms
    ! dtmp = 10. + 10.*(i-1)
    ! dtmp1 = set_eturb(dtmp,dlength,dsigma)
    ! write(*,'(2es13.5)')dtmp, dtmp1
    ! if (dtmp /= dtmp) write(*,*) "Some errors"
    ! write(*,'(3es13.5)') dkm(i),set_eturb(dkm(i),dlength,dsigma), dqm(i)
    write(123,'(3es13.5)') dkm(i),set_eturb(dkm(i),dlength,dsigma), dqm(i)
  enddo
  CLOSE(123)

  !fill the coefficients in from tmp_mod
  do i = 1, 3
    ac_m(i,:) = dsim_i(i,:) * dqm(:)
    b_m(i,:)  = dk_i(i,:) * dkm(:)
    c_m(i,:)  = dpsi(:)
  enddo

  write(*,*) "END: Generate Spectrum profile"

  ! CALL RANDOM_NUMBER(vels(1:3,1:ix,1:iy,1:iz))

  RETURN
end subroutine gen_flow_3d


!----------------------------------------------------------------------
!@brief Allocate modes arrays
!!
subroutine tmp_alloc(nels)
  USE tmp_mod
  implicit none
  integer, intent(IN) :: nels  !< number of modes
  
  ALLOCATE(ac_m(1:3,1:nels))
  ALLOCATE(as_m(1:3,1:nels))
  ALLOCATE(b_m(1:3,1:nels))
  ALLOCATE(c_m(1:3,1:nels))
  ALLOCATE(dphi_m(1:3,1:nels))

  ac_m(:,:) = 0.0d0
  as_m(:,:) = 0.0d0
  b_m(:,:) = 0.0d0
  c_m(:,:) = 0.0d0
  dphi_m(:,:) = 0.0d0
  RETURN
end subroutine tmp_alloc