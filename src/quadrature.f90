! *****************************************************************************
! * WRAPPER FUNCTIONS FOR THE INTEGRATION ROUTINES
! *****************************************************************************
module quadrature
  use quadpack
  implicit none
  double precision, parameter :: eps_G = 1d-3 !desired integration accuracy
  integer, parameter :: maxpts = 1000000 !maximum allowed invocations of integrand
contains
! *****************************************************************************
! * MULTI-DIMENSIONAL INTEGRATORS
! *****************************************************************************
  ! Calls CUBA's fastest non-deterministic (Monte Carlo) integrator
  subroutine cuba_divonne(Integrand, ndim, val, err, n)
    interface
       subroutine Integrand(ndim, x, ncomp, f)
         integer :: ndim, ncomp
         double precision x(ndim), f(ncomp)
       end subroutine Integrand
    end interface
    integer, intent(in) :: ndim
    real, intent(out) :: val, err
    integer, intent(out) :: n
    integer :: nregions, neval, fail
    double precision :: integral(ndim), error(ndim), prob(ndim)
    call divonne(ndim, 1, Integrand,&
         &    eps_G, 1d-12, 0, 0, maxpts,&
         &    47, 1, 1, 5,&
         &    0d0, 10d0, 0.25d0,&
         &    0, ndim, 0, 0, 0,&
         &    nregions, neval, fail, integral, error, prob)
    val = integral(1)
    err = error(1)
    n = neval
  end subroutine cuba_divonne

  ! Calls CUBA's deterministic integrator
  subroutine cuba_cuhre(Integrand, ndim, val, tol, err, n)
    interface
       subroutine Integrand(ndim, x, ncomp, f)
         integer, intent(in) :: ndim, ncomp
         double precision x(ndim), f(ncomp)
       end subroutine Integrand
    end interface
    integer, intent(in) :: ndim
    real, intent(out) :: val, err
    real, intent(in) :: tol
    integer, intent(out) :: n
    double precision :: eps
    integer :: nregions, neval, fail
    double precision :: integral(ndim), error(ndim), prob(ndim)
    eps = tol
    call cuhre(ndim, 1, Integrand,&
         &    eps, 1d-12, 0, 0, maxpts,&
         &    0, &
         &    nregions, neval, fail, integral, error, prob)
    val = integral(1)
    err = error(1)
    n = neval
  end subroutine cuba_cuhre

! *****************************************************************************
! * ONE-DIMENSIONAL INTEGRATORS
! *****************************************************************************
  ! Manual quadrature, which doesn't use Netlib.  It Just adds up boxes.
  function f_quad(f_dy_dx, x1, x2, N, argv) result(val)
    real, intent(in) :: x1, x2
    integer, intent(in) :: N
    real, intent(in), dimension(:), optional :: argv
    interface
       real function f_dy_dx(x, argv)
         real, intent(in) :: x
         real, intent(in), dimension(:), optional :: argv
       end function f_dy_dx
    end interface
    real abserr,epsabs,epsrel,f,val,work
    integer ier,neval
    real :: dx
    integer :: i
    real, dimension(N) :: dy_dx
    epsabs = 0.0
    epsrel = eps_G
    dx = (x2 - x1)/(N - 1)
    dy_dx(1) = 0.5*f_dy_dx(x1, argv)
    dy_dx(N) = 0.5*f_dy_dx(x2, argv)
    do i = 2, N - 1
       dy_dx(i) = f_dy_dx(x1 + (i - 1)*dx, argv)
    end do
    val = sum(dy_dx)*dx
  end function f_quad

 ! Uses Netlib's adaptive quadrature algorithm 'qags'
 function f_qags(f_dy_dx, x1, x2, N, argv) result(val)
   real, intent(in) :: x1, x2
   integer, intent(in) :: N
   real, intent(in), dimension(:), optional :: argv
   interface
      real function f_dy_dx(x, argv)
        real, intent(in) :: x
        real, intent(in), dimension(:), optional :: argv
      end function f_dy_dx
   end interface
   real abserr,epsabs,epsrel,f,val,work
   integer ier,iwork,last,lenw,limit,neval
   dimension iwork(100),work(400)
   epsabs = 0.0
   epsrel = eps_G
   limit = N
   lenw = limit*4
   call qags(f_dy_dx,x1,x2,epsabs,epsrel,val,abserr,neval,ier,argv)
 end function f_qags

 ! Uses Netlib's adaptive quadrature algorithm 'qag'
 function f_qag(f_dy_dx, x1, x2, N, argv) result(val)
   real, intent(in) :: x1, x2
   integer, intent(in) :: N
   real, intent(in), dimension(:), optional :: argv
   interface
      real function f_dy_dx(x, argv)
        real, intent(in) :: x
        real, intent(in), dimension(:), optional :: argv
      end function f_dy_dx
   end interface
   real abserr,epsabs,epsrel,f,val,work
   integer ier,iwork,last,lenw,limit,neval
   dimension iwork(100),work(400)
   epsabs = 0.0
   epsrel = eps_G
   limit = N
   lenw = limit*4
   call qag(f_dy_dx,x1,x2,epsabs,epsrel,1,val,abserr,neval,ier,argv)
 end function f_qag

 ! Uses Netlib's quadrature algorithm for open (infinite) intervals, 'qagi'
 function f_qagi(f_dy_dx, x1, inf, argv) result(val)
   real, intent(in) :: x1
   integer, intent(in) :: inf
   real, intent(in), dimension(:), optional :: argv
   interface
      real function f_dy_dx(x, argv)
        real, intent(in) :: x
        real, intent(in), dimension(:), optional :: argv
      end function f_dy_dx
   end interface
   real abserr,epsabs,epsrel,f,val,work
   integer ier,iwork,last,lenw,limit,neval
   dimension iwork(100),work(400)
   epsabs = 0.0e0
   epsrel = eps_G
   limit = 100
   lenw = limit*4
   call qagi(f_dy_dx,x1,inf,epsabs,epsrel,val,abserr,neval,ier,argv)
 end function f_qagi

 ! Uses's Netlib's non-adaptive integrator 'qng'
 function f_qng(f_dy_dx, x1, x2, N, argv) result(val)
   real, intent(in) :: x1, x2
   integer, intent(in) :: N
   real, intent(in), dimension(:), optional :: argv
   interface
      real function f_dy_dx(x, argv)
        real, intent(in) :: x
        real, intent(in), dimension(:), optional :: argv
      end function f_dy_dx
   end interface
   real abserr,epsabs,epsrel,f,val,work
   integer ier,neval
   epsabs = 0.0
   epsrel = eps_G
   call qng(f_dy_dx,x1,x2,epsabs,epsrel,val,abserr,neval,ier,argv)
 end function f_qng
end module quadrature
