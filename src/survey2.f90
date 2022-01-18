module survey2
  use cosmology
  use structure
  implicit none
  private
  public :: setenv, f_dN_dz, f_dN_dV
  real :: Mlim_G, z_G
contains
  subroutine setenv(Mlim)
    real, intent(in) :: Mlim
    Mlim_G = Mlim
  end subroutine setenv

  real function f_Mlim (z)
    real, intent(in) :: z
    f_Mlim = Mlim_G
  end function f_Mlim

  real function f_integrand (mu, argv)
    real, intent(in) :: mu
    real, intent(in), dimension(:), optional :: argv
    real :: z, M
    M = 10.0**mu
    z = z_G
    f_integrand = f_dN_dV_dmu(M, z)
  end function f_integrand

  real function f_dN_dV (z)
    real, intent(in) :: z
    real :: mu1
    mu1 = log10(f_Mlim(z))
    z_G = z
    f_dN_dV = f_qag(f_integrand, mu1, log10(1e16), 1, [0.])
  end function f_dN_dV

  real function f_dN_dz (z, dOmega)
    real, intent(in) :: z, dOmega
    f_dN_dz = f_dN_dV(z) * f_dV_dOmega_dz(z) * dOmega * (c/H100)**3
  end function f_dN_dz
end module survey2
