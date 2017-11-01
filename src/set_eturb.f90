!@brief  Calculate $E(k)$  by Karman formula
!@detail Для генерации изотропной турбулентности используем метод, предложенный в работе Rogallo [1]. 
!!
!! Функция Кармана спектральной плотности энергии:
!!  $${ E(k) = \frac{\sigma^2  L }{2 \pi} \frac{(\alpha k L)^4}{(1+(\alpha k L)^2)^{17/6}} }$$,
!! где _k_ - волновая частота, $\sigma = \sqrt{u'^2+v'^2+w'^2}$ - среднеквадратичное отклонение; 
!!  $\alpha$ - коэффициент ??? ($\alpha\approx 1.339$). _L_ - масштаб турбулентности.
!!
function set_eturb(dk, dl_in, dsigma_in) result(de)
  USE prec_mod
  implicit none
  real(prec), intent(IN) :: dk !< wave number
  real(prec), intent(IN), optional :: dl_in !< Integral Length Scale (set by default 0.1 m)
  real(prec), intent(IN), optional :: dsigma_in !< RMS of velocities (set by default 10 m/s)
  real(prec) :: de !< return value

  real(prec), save :: dl_def = 1.0d-1
  real(prec), save :: dsig_def = 1.0d+1
  real(prec) :: dl
  real(prec) :: dsig
  real(prec) :: dalpha
  real(prec) :: dtmp

! check if presented
  dl = dl_def ! [m]
  if(present(dl_in)) dl = dl_in

  dsig = dsig_def ! [m/s]
  if(present(dsigma_in)) dsig = dsigma_in

  ! write(*,*) "set_eturb: dl and dsig", dl,dsig

  dalpha = 1.339d0
  dtmp = (dalpha * dl * dk ) ** 2

  de = dl * 5.5d+1 / f_pi / 2.7d+1
  de = de * (dtmp * dsig)**2
  de = de / (1.d0 + dtmp)**(1.7d+1/6.0d0)

  return
end function set_eturb
