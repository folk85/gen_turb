!@brief generate random unit vector
!!  Use following formula for XY plane
!!  $$ e_x = r \cdot cos(\phi) $$
!!  $$ e_y = r \cdot sin(\phi) $$
!! For case in 3D space we change $r$ to $sin$ and $cos$
!!  $$ e_x = sin(\theta) \cdot cos(\phi) $$
!!  $$ e_y = sin(\theta)  \cdot sin(\phi) $$
!!  $$ e_y = cos(\theta) $$
!! We have to define $r$ non-negative, therefore define 
!!  $$ \theta = arccos(\tau) $$ 
!!  $$ \tau = \textit{N}(-1,1) $$ 
!! The angle $\phi \in \textit{N}(-2\pi;2\pi) $, where 
!! $\textit{N}(a,b)$ is uniform distribution in range $[a,b)$ 
!------------------------------------------------------------------
pure function set_unit_vector(dk) result(dv)
  USE prec_mod
  implicit none
  real*8,dimension(*), intent(IN)  :: dk !< initial angles in polar system of coordinatesgeneration of random numbers in [0,1)
  real*8,dimension(1:3) :: dv !< Return vector 
  real(prec) :: dtheta !< angle [0,Pi]
  real(prec) :: dtau !< temporal number  [-1, 1)
  real(prec) :: dphi !< angle [-Pi,Pi)

  
  ! CALL RANDOM_NUMBER(dk)

  dphi = dk(1) * 2.0d0 * f_pi
  dtau = dk(2) * 2.0d0 - 1.0d0
  dtheta = ACOS(dtau)
  dv(1) = SIN(dtheta) * COS(dphi)
  dv(2) = SIN(dtheta) * SIN(dphi)
  dv(3) = COS(dtheta)
  RETURN
end function set_unit_vector

!----------------------------------------------------------------------
!@brief Calculate normalized cross product
!!
function cross_product(a, b) !result(r)
  USE prec_mod
  implicit none
  real*8, dimension(3)              :: cross_product !< result vector
  real(prec), dimension(3) , intent(IN) :: a !< first vector
  real(prec), dimension(3) , intent(IN) :: b !< second vector
  real(prec), dimension(3)              :: r !< result vector
  real(prec) :: d !< result vector
  integer :: i

  write(*,*) "Calc cross product"
  do i=1, 3
    write(*,'(i10,2es10.2)') i, a(i), b(i)
  enddo

  r(1) = a(2) * b(3) - a(3) * b(2)
  r(2) = a(3) * b(1) - a(1) * b(3)
  r(3) = a(1) * b(2) - a(2) * b(1)

  write(*,*) "|r|"

  d = 1.0d0 !SQRT(SUM(r(1:3)**2))

  write(*,'(a,3es10.2)') "result",r(1:3)

  write(*,*) "Noramlize"

  do i=1,3
    cross_product(i)  = r(i) / (d + small)
  enddo
  ! cross_product(1:3) = r(1:3) !/ (d + small)
  RETURN     
end function cross_product

!----------------------------------------------------------------------
!@brief Calculate normalized cross product
!!
pure subroutine cross_product_sub(a, b, cr_p) !result(r)
  USE prec_mod
  implicit none
  real(prec), dimension(3) , intent(IN) :: a !< first vector
  real(prec), dimension(3) , intent(IN) :: b !< second vector
  real(prec), dimension(3) , intent(OUT) :: cr_p !< second vector
  real(prec), dimension(3)              :: r !< result vector
  real(prec) :: d !< result vector

  r(1) = a(2) * b(3) - a(3) * b(2)
  r(2) = a(3) * b(1) - a(1) * b(3)
  r(3) = a(1) * b(2) - a(2) * b(1)

  d = SQRT(SUM(r(1:3)**2))

  cr_p(1:3) = r(1:3) / (d + small)

  RETURN     
end subroutine cross_product_sub


subroutine check_cross()
  USE prec_mod
    implicit none
    real(prec),dimension(3) :: i,j,k
    real(prec),dimension(3) :: res
    real*8 :: cross_product
    integer :: m

    i(1:3) = (/ 1. , 0. , 0. /)
    j(1:3) = (/ 0. , 1. , 0. /)
    k(1:3) = (/ 0. , 0. , 1. /)
    do m=1, 3
      write(*,'(i10,3es10.2)') m, i(m),j(m),k(m)
    enddo

    write(*,*) "Call cross_product"
    ! res = cross_product(i(1:3),j(1:3))
    CALL cross_product_sub(i,j,res)
    write(*,'(3es10.2)') res
    CALL cross_product_sub(j,k,res)
    write(*,'(3es10.2)') res
    CALL cross_product_sub(k,i,res)
    write(*,'(3es10.2)') res

    return
end subroutine check_cross