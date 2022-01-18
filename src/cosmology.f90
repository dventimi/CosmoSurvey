! Cosmological model, mainly with functions for calculating geometry
module cosmology
  use quadrature
  use constants
  implicit none
  type CosmoParams
     real :: w               !Dark energy equation of state parameter
     real :: Omega_M_0       !Matter parameter
     real :: Omega_B_0       !Baryon density parameter
     real :: sigma_8         !Matter power spectrum normalization
     real :: tau_0           !M-T relationship normalization
     real :: lambda_0        !L-M relationship normalization
     real :: alpha           !M-T relationship slope
     real :: beta            !L-M relationship slope
     real :: sig_tau         !Mass-observable scatter parameter
     real :: sig_lambda      !Mass-observable scatter parameter
     real :: phi             !Mass-observable scatter parameter
     real :: psi             !Mass-observable scatter parameter
     real :: rho             !Mass-observable scatter parameter
     real :: chi             !Mass-observable scatter parameter
     real :: gamma           !Linder growth function parameter
  end type CosmoParams
  ! Model parameters are stored in a global data structure,
  ! a global variable of type CosmoParams, which is a user-defined type.
  ! We use a global variable partly because it's the most convenient way
  ! to make parameter data available to the various integrands.
  type(CosmoParams) :: theta_G
  !$OMP THREADPRIVATE (theta_G)
contains
  ! Function to return fiducial set of parameters
  ! (Note:  in a Minuit setting, these typically will be overwritten)
  ! These are the standard LCDM parameters
  function f_lcdm() result(theta)
    type(CosmoParams) :: theta
    theta%w = -1.0              !1
    theta%Omega_M_0 = 0.30      !3 
    theta%Omega_B_0 = 0.02      !4
    theta%sigma_8 = 0.10        !5 
    theta%tau_0 = 0.            !6
    theta%lambda_0 = 0.         !7
    theta%alpha = 2./3.         !8
    theta%beta = 1.8            !9
    theta%sig_tau = 0.1         !10
    theta%sig_lambda = 0.1      !11
    theta%phi = 0.              !12
    theta%psi = 0.              !13
    theta%rho = 0.              !14
    theta%chi = 0.              !15
    theta%gamma = 0.55          !16
  end function f_lcdm

  ! Function to return the priors associated with the lcdm param vals
  function f_lcdm_priors() result(priors)
    real, dimension(16,2) :: priors
    priors(1,:) = [-1.02, -0.70]
    priors(2,:) = [0.65, 0.80]
    priors(3,:) = [0.24, 0.36]
    priors(4,:) = [0.00, 0.10]
    priors(5,:) = [0.80, 0.90]
    priors(6,:) = [-0.50, 0.50]
    priors(7,:) = [-0.50, 0.50]
    priors(8,:) = [2./3*0.9, 2./3/0.9]
    priors(9,:) = [1.8*0.9, 1.8/0.9]
    priors(10,:) = [0.00, 1.00]
    priors(11,:) = [0.00, 1.00]
    priors(12,:) = [0.00, 1.00]
    priors(13,:) = [0.00, 1.00]
    priors(14,:) = [0.00, 1.00]
    priors(15,:) = [0.00, 1.00]
    priors(16,:) = [0.50, 0.60]
  end function f_lcdm_priors

  function f_lcdm_names() result(names)
    character(len=17), dimension(16) :: names
    names = [&
          &'{/=20 w}        ', &   !1
          &'{/Symbol=20 W_L}', &   !2
          &'{/Symbol=20 W}_m', &   !3
          &'{/Symbol=20 W}_b', &   !4
          &'{/Symbol=20 s}_8', &   !5
          &'{/Symbol=20 t}_0', &   !6
          &'{/Symbol=20 l}_0', &   !7
          &'{/Symbol=20 a}  ', &   !8
          &'{/Symbol=20 b}  ', &   !9
          &'{/Symbol=20 s_t}', &   !10
          &'{/Symbol=20 s_l}', &   !11
          &'{/Symbol=20 f}  ', &   !12
          &'{/Symbol=20 y}  ', &   !13
          &'{/Symbol=20 r}  ', &   !14
          &'{/Symbol=20 c}  ', &   !15
          &'{/Symbol=20 g}  ']     !16
  end function f_lcdm_names

  ! E(z) = H(z)/H_0
  function f_E(z) result(E)
    real, intent(in) :: z
    real :: E, Omega_0, Omega_kappa_0, Omega_lambda_0
    Omega_lambda_0 = 1.0 - theta_G%Omega_M_0
    E = sqrt(theta_G%Omega_M_0*(1.+z)**3 + &
         &Omega_lambda_0*(1.+z)**(3*(1.+theta_G%w)))
  end function f_E

  ! dr/dz where r is comoving distance
  function f_dr_dz(z, argv) result(dr_dz)
    real, intent(in) :: z
    real, intent(in), dimension(:), optional :: argv
    real :: E, dr_dz
    E = f_E(z)
    dr_dz = 1./E
  end function f_dr_dz

  ! Comoving distance r between z1 and z2
  function f_r(z1, z2) result(r)
    real, intent(in) :: z1, z2
    real :: r, Omega_kappa_0, Omega_0, x, rho, R_kappa
    r = f_qag(f_dr_dz, z1, z2, 10)
  end function f_r

  ! ! Radius of curvature
  ! function f_R_kappa() result(R_kappa)
  !   real :: R_kappa, Omega_0
  !   Omega_0 = theta_G%Omega_M_0 + theta_G%Omega_lambda_0
  !   R_kappa = sqrt(1./abs((Omega_0 - 1.)))
  ! end function f_R_kappa

  ! dV/(dz*dOmega) where V is comoving volume element
  real function f_dV_dOmega_dz (z, argv)
    real, intent(in) :: z
    real, intent(in), dimension(:), optional :: argv
    real :: r, E
    r = f_r(0., z)
    E = f_E(z)
    f_dV_dOmega_dz = 1./E*r**2
  end function f_dV_dOmega_dz

  ! Comoving volume between z1 and z2, in solid angle dOmega
  function f_dV(z1, z2, dOmega) result(dV)
    real, intent(in) :: z1, z2, dOmega
    real :: dV
    dV = f_qag(f_dV_dOmega_dz, z1, z2, 10)*dOmega
  end function f_dV

  ! Matter parameter Omega_M(z) as a function of redshift
  function f_Omega_M(z) result(Omega_M)
    real, intent(in) :: z
    real :: Omega_M, E
    E = f_E(z)
    Omega_M = theta_G%Omega_M_0*(1.+z)**3/E**2
  end function f_Omega_M

  ! Dark energy parameter Omega_lambda(z) as a function of redshift
  function f_Omega_lambda(z) result(Omega_lambda)
    real, intent(in) :: z
    real :: Omega_lambda, E
    E = f_E(z)
    Omega_lambda = (1. - theta_G%Omega_M_0)/E**2
  end function f_Omega_lambda

  ! Distance measure out to redshift z
  function f_D(z) result(D)
    real, intent(in) :: z
    real :: D, Omega_K_0
    D = f_r(0., z)
  end function f_D

  ! Luminosity distance out to redshift z
  function f_D_L(z) result(D_L)
    real, intent(in) :: z
    real :: D_L
    D_L = (1.+z)*f_D(z)
  end function f_D_L
end module cosmology
