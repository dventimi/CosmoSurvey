module survey3
  use quadrature
  use constants
  use cosmology
  use structure
  implicit none
  private
  public :: f_dN, f_N
  real, parameter :: tol = 1e-3
  integer, parameter :: scale = 3
  real :: z_G, l1_G, l2_G
contains
  ! Probability density function P(lambda|mu)
  real function f_pdf(lambda, mu, z)
    real, intent(in) :: lambda, mu, z
    real :: expected_lambda, sig_lambda
    real, dimension(2, 2) :: C
    C = f_Cov(z)
    sig_lambda = sqrt(C(2, 2))
    expected_lambda = f_lambda_of_mu(mu, z)
    f_pdf = 1./(sqrt(2.*pi)*sig_lambda)*exp(-0.5*(lambda - expected_lambda)**2/sig_lambda**2)
  end function f_pdf

  ! Number density per unit log(M) per unit log(L)
  subroutine dN_dV_dl_dmu (ndim, u, ncomp, f)
    integer, intent(in) :: ndim, ncomp
    double precision :: u(ndim), f(ncomp)
    real :: P, dN_dV_dmu, mu1, mu2, dmu, mu, C(2, 2), sig_l, dl, l
    C = f_Cov(z_G)
    sig_l = sqrt(C(2, 2))
    dl = l2_G - l1_G
    l = u(1)*dl + l1_G
    mu1 = f_mu_of_lambda(l - scale*sig_l, z_G)
    mu2 = f_mu_of_lambda(l + scale*sig_l, z_G)
    dmu = mu2 - mu1
    mu = u(2)*dmu + mu1
    P = f_pdf(l, mu, z_G)*2*scale*sig_l
    dN_dV_dmu = f_dN_dV_dmu(10.**mu, z_G)
    f = dN_dV_dmu * P
  end subroutine dN_dV_dl_dmu

  ! Number per comoving volume (number density)
  real function f_dN_dV (z, l1, l2, flux)
    real, intent(in) :: z, l1, l2, flux
    real :: intval, abserr
    integer :: neval
    z_G = z
    l1_G = max(f_luminosity(flux, z), l1)
    l2_G = l2
    if (l2_G.gt.l1_G) then
       call cuba_cuhre(dN_dV_dl_dmu, 2, intval, 1e-3, abserr, neval)
    else
       intval = 0.
    end if
    f_dN_dV = intval
  end function f_dN_dV

  ! Number per solid angle per unit redshift
  real function f_dN_dOmega_dz (z, argv)
    real, intent(in) :: z
    real, intent(in), dimension(:), optional :: argv
    real :: l1, l2, flux
    l1 = argv(1)
    l2 = argv(2)
    flux = argv(3)
    f_dN_dOmega_dz = f_dV_dOmega_dz(z) * f_dN_dV(z, l1, l2, flux)
  end function f_dN_dOmega_dz

  ! Number per solid angle
  real function f_dN_dOmega (z1, z2, l1, l2, flux)
    real, intent(in) :: z1, z2, l1, l2, flux
    f_dN_dOmega = f_qag(f_dN_dOmega_dz, z1, z2, 10, [l1, l2, flux])
  end function f_dN_dOmega

  ! Number (in a given redshift-luminosity bin)
  real function f_dN (z1, z2, l1, l2, dOmega, flux)
    real, intent(in) :: z1, z2, l1, l2, dOmega, flux
    f_dN = f_dN_dOmega(z1, z2, l1, l2, flux)*dOmega*(c/H100)**3
  end function f_dN

  ! Survey (vector of cluster numbers in redshift-luminosity bins)
  function f_N(theta, z1, z2, l1, l2, nz, nl, dOmega, flux_lim)
    type(CosmoParams), intent(in) :: theta
    real, intent(in) :: z1, z2, l1, l2, dOmega, flux_lim
    integer, intent(in) :: nz, nl
    real, dimension(:), pointer :: f_N
    real :: dz, dl
    integer :: i, j, idx
    theta_G = theta
    dz = (z2 - z1)/real(Nz)
    dl = (l2 - l1)/real(Nl)
    allocate(f_N(Nz*Nl))
    do i = 1, Nz
       do j = 1, Nl
          idx = (i-1)*Nl + j
          f_N(idx) = f_dN(z1+(i-1)*dz, z1+i*dz, l1+(j-1)*dl, l1+j*dl, dOmega, flux_lim)
       end do
    end do
  end function f_N
end module survey3
