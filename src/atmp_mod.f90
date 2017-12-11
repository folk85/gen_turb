!@brief Module considers arrays to caclculate vels
!!   $$\textbf{u} = 2 \sum_{m=1}^M Ac_m cos(\textbf(B_m}\cdot \textbf{x} + C_m) $$
!! 
!! here we define angle `dphi_m` = $\phi_m = \textbf(B_m}\cdot \textbf{x} + C_m$
!!
module tmp_mod
  USE prec_mod
  implicit none
  integer :: nmodes                                   !< number of Modes
  integer :: ntimes                                   !< number of timesteps
  real(prec), dimension(:,:), allocatable :: ac_m     !< coefficient amplitude with COS
  real(prec), dimension(:,:), allocatable :: as_m     !< coefficient amplitude with SIN
  real(prec), dimension(:,:), allocatable :: b_m      !< coefficient vector K inside SIN and COS
  real(prec), dimension(:), allocatable :: c_m        !< random numbers inside SIN and COS
  real(prec), dimension(:,:), allocatable :: dphi_m   !< caclculated angle for SIN and COS
  real(prec), dimension(:), allocatable :: cs_m       !< Internal temporary coefficient

!-----
!  Use temporary variables for coordinates and velocities in cells
  integer :: tcell                                    !< number of cells
  ! real(prec), dimension(:,:), allocatable :: xp       !< Cell Cartesian coordinates
  ! real(prec), dimension(:,:), allocatable :: u        !< velocities in Cells
  integer :: in_time                                  !< Current ttimestep. Used in ge_flow
  real(prec), dimension(:), allocatable :: dtim       !< times in array of timesteps
end module tmp_mod