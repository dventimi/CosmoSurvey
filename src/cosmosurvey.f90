program main
  use M_kracken
  use utils
  use quadrature
  use constants
  use cosmology
  use survey3
  implicit none
  real :: z1, z2, l1, l2, dOmega, flux_lim
  real :: w, Omega_Lambda_0, Omega_M_0, sigma_8, tau_0, lambda_0, alpha, beta, sig_tau
  real :: sig_lambda, phi, psi, rho, chi, gamma
  real, dimension(:), pointer :: dN
  integer :: Nz, Nl
  type(CosmoParams) :: theta
  character(len=255) :: cmd

  cmd = &
       &' -h .false. -z1 0.01 -z2 2.5 -nz 5 -l1 44.0 -l2 46.0 -nl 5 '//&
       &' -do 12.0 -fl 1.25e-13 -w -1.0 -om 0.30 -s8 0.8 -g 0.55 '
  call kracken ('cmd', cmd)
  if (lget('cmd_h')) then
     print *, ''
     print *, 'Usage:'
     print *, ''
     print *, '  cosmosurvey [[OPTION]]'
     print *, ''
     print *, 'Options (Default):'
     print *, '  -h                Print help information'
     print *, '  -z1  (0.01)       Redshift start'
     print *, '  -z2  (2.50)       Redshift end'
     print *, '  -nz  (5)          Number of redshift bins'
     print *, '  -l1  (44.0)       log10(L1) luminosity start, [L1] = ergs/s'
     print *, '  -l2  (46.0)       log10(L2) luminosity end, [L2] = ergs/s'
     print *, '  -nl  (5)          Number of luminosity bins'
     print *, '  -do  (12.0)       Survey solid angle dOmega, [dOmega] = steradians'
     print *, '  -fl  (1.25e-13)   Survey flux limit, [fl] = ergs/s/cm^2'
     print *, '  -w   (-1.0)       Dark Energy equation of state parameter w'
     print *, '  -om  (0.30)       Matter density parameter Omega_M'
     print *, '  -s8  (0.80)       Power spectrum normalization sigma_8'
     print *, '  -g   (0.55)       Linder growth index gamma'
     print *, ''
  else
     theta = f_lcdm()
     theta%w = rget('cmd_w')
     theta%Omega_M_0 = rget('cmd_om')
     theta%sigma_8 = rget('cmd_s8')
     theta%gamma = rget('cmd_g')
     dN => f_N(theta, &
          &rget('cmd_z1'), &
          &rget('cmd_z2'), &
          &rget('cmd_l1'), &
          &rget('cmd_l2'), &
          &iget('cmd_nz'), &
          &iget('cmd_nl'), &
          &rget('cmd_do'), &
          &rget('cmd_fl'))
     print *, dN
     deallocate(dN)
  end if
end program main

