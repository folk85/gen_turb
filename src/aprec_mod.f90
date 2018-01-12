!@brief Module with definition of precision
module prec_mod
implicit none

  integer, parameter    :: prec = kind(1.0d0)  !< define precision mode
  ! real(prec), parameter :: f_zero = 1.0_prec !< Pi in Fire used precision. 
  ! real(prec), parameter :: f_pi = 3.141592653589793238462643383279502884197_prec !< Pi in Fire used precision. 
  ! real(prec), parameter :: small = 1.0d-30 !< Small value
end module prec_mod