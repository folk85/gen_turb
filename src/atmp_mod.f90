!@brief Module considers arrays to caclculate vels
!!   $$\textbf{u} = 2 \sum_{m=1}^M Ac_m cos(\textbf(B_m}\cdot \textbf{x} + C_m) $$
!! 
!! here we define angle `dphi_m` = $\phi_m = \textbf(B_m}\cdot \textbf{x} + C_m$
!!
module tmp_mod
  USE prec_mod
  implicit none
  real(prec),dimension(:,:), allocatable :: ac_m   !< coefficient amplitude with COS
  real(prec),dimension(:,:), allocatable :: as_m   !< coefficient amplitude with SIN
  real(prec),dimension(:,:), allocatable :: b_m    !< coefficient vector K inside SIN and COS
  real(prec),dimension(:,:), allocatable :: c_m    !< random numbers inside SIN and COS
  real(prec),dimension(:,:), allocatable :: dphi_m    !< caclculated angle for SIN and COS
contains
    
end module tmp_mod