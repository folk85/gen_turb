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
  integer :: i, i1
  integer :: j
  integer :: k
  integer :: ion                            !< number of unitfor I/O
  integer, dimension(3)  :: nels            !< 
  real(prec), dimension(3)  :: dels         !<
  real(prec), dimension(3)  :: tmp3         !< temporary
  ! real(prec), allocatable :: u(:,:,:,:)

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
  ncell  = nx * ny * nz
  CALL tmp_alloc()

  !define coordinates
  do i = 1, nx
    tmp3(1) = dlx * (2*i-1)*5.d-1/ nx
    do j = 1, ny
      tmp3(2) = dly * (2*j-1)*5.d-1/ ny
      do k = 1, nz
        tmp3(3) = dlz * (2*k-1)*5.d-1/ nz
        i1 = (i-1) * ny * nz + nz * (j-1) + k
        xp(1:3,i1) = tmp3(1:3)
      enddo
    enddo
  enddo
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

  write(*,*) "Generate  the modes"
  ! generate fields
  CALL gen_flow_3d(dels, nels, dsigma, dlength, dtau)


  write(*,*) "Calculate the velocities"
  CALL set_vels()

  write(*,*) "Store the array"
  !store the field
  ion = 121
  OPEN(ion,file='store.dat')
  write(ion,'(3i5)') nx, ny, nz
  do i=1, ncell
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
  
  ALLOCATE(ac_m(1:3,1:nmodes))
  ALLOCATE(as_m(1:3,1:nmodes))
  ALLOCATE(b_m(1:3,1:nmodes))
  ALLOCATE(c_m(1:3,1:nmodes))
  ALLOCATE(dphi_m(1:3,1:nmodes))

  ac_m(:,:) = 0.0d0
  as_m(:,:) = 0.0d0
  b_m(:,:) = 0.0d0
  c_m(:,:) = 0.0d0
  dphi_m(:,:) = 0.0d0

!-----
!  Use temporary variables for coordinates and velocities in cells
  ALLOCATE(xp(1:3,1:ncell))
  ALLOCATE(u(1:3,1:ncell))

  xp(:,:) = 0.0d0
  u(:,:) = 0.0d0
  RETURN
end subroutine tmp_alloc