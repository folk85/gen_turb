!>@brief (from FIRE ) involve module COMM0 from FIRE 
!!
module comm0

  USE prec_mod, ONLY: prec
  implicit none
  integer :: ncell                 !< number of cells (in FIRE notation)
  integer, dimension(1:1) :: nsp  ! (from FIRE) number of starting cell in domain
  integer, dimension(1:1) :: nep  ! (from FIRE) number of ending cell in domain
  integer :: I_USEINI
  integer :: itst
  real(prec) :: time
end module comm0