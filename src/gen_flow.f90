!@breaf  main porgram algorithm 
!!
program gen_flow_saad
  USE prec_mod
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

  write(*,*) "Welcome into th program"
  ! set  time duration
  dlt = 0.0d0
  ! set space Length
  dlx = 2.0d0 * f_pi * 1.0d-2
  dly = dlx
  dlz = dlx
  ! set number of nodes
  nx = 64
  ny = nx
  nz = nx
  !set number of Modes
  nmodes = nx 

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
  CALL gen_flow_3d(u, dels, nels,nmodes)

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


subroutine gen_flow_3d(vels,dls,nels,nms)
  USE prec_mod
  implicit none
  integer, dimension(1:3), intent(IN) :: nels   !< number of modes
  real(prec), dimension(1:3,1:nels(1),1:nels(2),1:nels(3)), intent(INOUT) :: vels !< in/Out velocities
  real(prec), dimension(3), intent(IN) :: dls !< Physical boundaries
  integer, intent(IN) :: nms !< number of modes
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

  integer :: i,j,k

  if(ALLOCATED(dkm)) DEALLOCATE(dkm)
  ALLOCATE(dkm(1:nms))
  dkm(:) = 0.0d0

  ix = nels(1)
  iy = nels(2)
  iz = nels(3)

  ! set cell width
  dx = dls(1) / REAL(nels(1),prec)
  dy = dls(2) / REAL(nels(2),prec)
  dz = dls(3) / REAL(nels(3),prec)

  ! define minimal wave Length
  kmin = 0.0d0
  kmax = 0.0d0
  do i=1, 3
    kmin = MAX(kmin,2.0d0 * f_pi / dls(i))
    kmax = MAX(kmax,2.0d0*f_pi * nels(i) / dls(i))
  enddo

  ! define inverals and 
  kwid = (kmax - kmin) / nms
  ! set waves for modes
  do i = 1, nms
    dkm(i) = kmin + kwid * (i - 1)
  enddo


  CALL RANDOM_NUMBER(vels(1:3,1:ix,1:iy,1:iz))

  RETURN
end subroutine gen_flow_3d

! !@brief Module with definition of precision
! module prec_mod
! implicit none

!   integer, parameter    :: prec = kind(1.0d0)  !< define precision mode
!   real(prec), parameter :: f_zero 1.0_prec !< Pi in Fire used precision. 
!   real(prec), parameter :: f_pi 3.141592653589793238462643383279502884197_prec !< Pi in Fire used precision. 
! end module prec_mod