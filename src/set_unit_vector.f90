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
function set_unit_vector() result(dv)
  USE prec_mod
  implicit none
  real*8,dimension(1:3) :: dv !< Return vector 
  real*8,dimension(2) :: dk !< generation of random numbers in [0,1)
  real(prec) :: dtheta !< angle [0,Pi]
  real(prec) :: dtau !< temporal number  [-1, 1)
  real(prec) :: dphi !< angle [-Pi,Pi)

  
  CALL RANDOM_NUMBER(dk)

  dphi = dk(1) * 2.0d0 * f_pi
  dtau = dk(2) * 2.0d0 - 1.0d0
  dtheta = acos(dtau)
  dv(1) = SIN(dtheta) * COS(dphi)
  dv(2) = SIN(dtheta) * SIN(dphi)
  dv(3) = COS(dtheta)
  RETURN
end function set_unit_vector