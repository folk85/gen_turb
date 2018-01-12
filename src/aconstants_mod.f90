!@brief Module with definition of precision
module constants
  USE prec_mod, ONLY: prec
  implicit none
  real(prec), parameter :: f_zero = 1.0_prec !< Pi in Fire used precision. 
  real(prec), parameter :: f_pi = 3.141592653589793238462643383279502884197_prec !< Pi in Fire used precision. 
end module constants