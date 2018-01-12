!>@brief (from FIRE ) involve module COMM1 from FIRE 
!!
module comm1
use prec_mod, only : prec
implicit none
  real(prec), dimension(:,:), allocatable :: xp       !< Cell Cartesian coordinates
  real(prec), dimension(:,:), allocatable :: u        !< velocities in Cells
    
end module comm1