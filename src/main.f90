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

  ! integer :: nmodes !< number of Modes

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
  CALL gen_flow_3d(dels, nels, dsigma, dlength, dtau)

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