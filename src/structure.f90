! Model for growth of structure
module structure
  use quadrature
  use constants
  use cosmology
  implicit none
contains
  ! The integrand used in the calculation of the growth factor G(z).
  function f_Integrand_Linder(z, argv) result(val)
    real, intent(in) :: z
    real, intent(in), dimension(:), optional :: argv
    real :: val
    val = (f_Omega_M(z)**theta_G%gamma - 1.)/(1.+z)
  end function f_Integrand_Linder

  ! Calculation of the growth factor
  function f_G(z) result(G)
    real, intent(in) :: z
    real :: G, G0
    G = exp(f_qagi(f_Integrand_Linder, z, 1))/(1.+z)
  end function f_G

  ! Shape-governing function sigma(M), save for a normalization factor (fit function)
  function f_ln_inv_sigma(M, z) result(ln_inv_sigma)
    real, intent(in) :: M, z
    real :: ln_15, a, b, x, g, ln_g, ln_inv_sigma
    ln_15 = -1*log(theta_G%sigma_8-0.3)
    a = 0.281
    b = 0.0123
    x = log(M/1d15)
    g = f_g(z)/f_g(0.0)
    ln_g = -1*log(g)
    ln_inv_sigma = ln_15 + a*x + b*x*x + ln_g
  end function f_ln_inv_sigma

  function f_alpha_eff (M) result(alpha_eff)
    real, intent(in) :: M
    real :: a, b, x, alpha_eff
    a = 0.281
    b = 0.0123
    x = log(M/1d15)
    alpha_eff = a + 2*b*x
  end function f_alpha_eff

  function f_mass_fraction (ln_inv_sigma) result(mass_fraction)
    real, intent(in) :: ln_inv_sigma
    real :: a, b, eps, mass_fraction
    a = 0.22
    b = 0.73
    eps = 3.86
    mass_fraction = a*exp(-1*(abs(ln_inv_sigma + b))**eps)
  end function f_mass_fraction

  function f_dN_dV_dmu (M, z) result(dN_dV_dm)
    real, intent(in) :: M, z
    real :: dN_dV_dm
    dN_dV_dm = f_Omega_M(z)*rho_crit*(M_scale/M)*f_alpha_eff(M)*f_mass_fraction(f_ln_inv_sigma(M, z))*log(10.0)
  end function f_dN_dV_dmu

  ! T(M) mass-observable relation
  function f_tau_of_mu(mu, z) result(tau_of_mu)
    real, intent(in) :: mu, z
    real :: tau_of_mu, E
    E = f_E(z)
    tau_of_mu = theta_G%tau_0 + theta_G%alpha*mu + theta_G%alpha*log(E)
  end function f_tau_of_mu

  ! Slope of the T(M) relation
  function f_dtau_dmu() result(dmu_dtau)
    real :: dmu_dtau
    dmu_dtau = theta_G%alpha
  end function f_dtau_dmu

  ! M(T) mass-observable relation (it's an inverse function)
  function f_mu_of_tau(tau, z) result(mu_of_tau)
    real, intent(in) :: tau, z
    real :: mu_of_tau, E
    E = f_E(z)
    mu_of_tau = (tau - theta_G%tau_0 - theta_G%alpha*log(E))/theta_G%alpha
  end function f_mu_of_tau

  ! L(M) mass-observable relation
  function f_lambda_of_mu(mu, z) result(lambda_of_mu)
    real, intent(in) :: mu, z
    real :: lambda_of_mu, E
    E = f_E(z)
    lambda_of_mu = theta_G%lambda_0 + (mu-log10(M_scale))*theta_G%beta + 2.*log10(E) + log10(L_scale)
  end function f_lambda_of_mu

  function f_mu_of_lambda(lambda, z) result(mu_of_lambda)
    real, intent(in) :: lambda, z
    real :: mu_of_lambda, E
    E = f_E(z)
    mu_of_lambda = (lambda - log10(1d45) - theta_G%lambda_0 - 2.*log10(E))/theta_G%beta + log10(M_scale)
  end function f_mu_of_lambda

  ! Get a covariance matrix in T & L
  function f_Cov(z) result(C)
    real, intent(in) :: z
    real, dimension(2, 2) :: C
    real :: stau, slambda, rho
    stau = theta_G%sig_tau*(1. + z)**theta_G%phi
    slambda = theta_G%sig_lambda*(1. + z)**theta_G%psi
    rho = theta_G%rho
    C(1,1) = stau*stau
    C(1,2) = rho*stau*slambda
    C(2,1) = rho*slambda*stau
    C(2,2) = slambda*slambda
  end function f_Cov

  ! Get the determinant of a covariance matrix C
  function f_C_det(C_in) result(C_det)
    real, dimension(:, :), intent(in) :: C_in
    real, dimension(size(C_in, 1), size(C_in, 1)) :: C
    real :: C_det
    C = C_in
    C_det = C(1,1)*C(2,2) - C(1,2)*C(2,1)
  end function f_C_det

  ! Calculates the luminosity of a source at redshift z
  function f_luminosity(flux, z) result(lum)
    real, intent(in) :: flux, z
    real(kind=selected_real_kind(8,75)) :: lum, D_L, L, T, f
    D_L = f_D_L(z)*(c/H100)*cm_per_Mpc
    L = 4.*pi*D_L**2*flux
    lum = log10(L)
  end function f_luminosity
end module structure
