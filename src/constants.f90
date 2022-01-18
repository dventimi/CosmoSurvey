! Global constants
module constants
  implicit none
  integer, parameter :: bignum = selected_real_kind(10, 100)
  real, parameter :: rho_crit = 2.77466274901e-4       !Msun/Mpc^3
  real, parameter :: H100 = 3.24077928534e-18          !s^-1
  real, parameter :: c = 9.71561187787e-15             !Mpc / s
  real, parameter :: cm_per_Mpc = 3.0856775854e24      !cm / Mpc
  real, parameter :: pi = 4.*atan(1.0)                 !pi
  real, parameter :: M_scale = 1e15                    !Msun
  real(kind=bignum), parameter :: L_scale = 1d45      !ergs/s
  real, parameter :: M_8 = 0.595070983403 !In units of M15
end module constants

