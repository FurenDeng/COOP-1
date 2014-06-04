module coop_function_mod
  use coop_constants_mod
  use coop_basicutils_mod
  use coop_arguments_mod
  implicit none
#include "constants.h"

  private

  public:: coop_function

  type coop_function
     private
     COOP_INT method
     logical xlog, ylog
     COOP_INT n
     COOP_REAL xmin, xmax, dx, fleft, fright, slopeleft, sloperight
     COOP_REAL, dimension(:), allocatable::f, f2
   contains
     procedure::set_boundary => coop_function_set_boundary
     procedure::init => coop_function_initialize
     procedure::eval => coop_function_evaluate
     procedure::integrate => coop_function_integrate
     procedure::derivative => coop_function_derivative
     procedure::free => coop_function_free
  end type coop_function


  interface coop_function
     module procedure coop_function_constructor
  end interface coop_function

contains

  subroutine coop_function_set_boundary(this, fleft, fright, slopeleft, sloperight)
    class(coop_function)::this
    COOP_REAL,optional::fleft, fright, slopeleft, sloperight
    if(present(fleft))this%fleft = fleft
    if(present(fright))this%fright = fright
    if(present(slopeleft))this%slopeleft = slopeleft
    if(present(sloperight))this%sloperight = sloperight
  end subroutine coop_function_set_boundary

  function coop_function_constructor(f, xmin, xmax, xlog, ylog, args) result(cf)
    external f
    COOP_REAL f, xmin, xmax, dx, lnxmin, lnxmax
    logical,optional::xlog
    logical,optional::ylog
    type(coop_arguments),optional::args
    COOP_REAL_ARRAY::y
    COOP_INT i
    type(coop_function) :: cf
    if(present(xlog))then
       cf%xlog = xlog
    else
       cf%xlog = .false.
    endif
    if(present(ylog))then
       cf%ylog = ylog
    else
       cf%ylog = .false.
    endif
    if(cf%xlog)then
       lnxmin = log(xmin)
       lnxmax = log(xmax)
       dx = (lnxmax - lnxmin)/(coop_default_array_size - 1)
       if(present(args))then
          y(1) = f(xmin, args)
          y(coop_default_array_size) = f(xmax, args)
          !$omp parallel do
          do i=2, coop_default_array_size-1
             y(i) = f(exp(lnxmin + dx*(i-1)), args)
          enddo
          !$omp end parallel do
       else
          y(1) = f(xmin)
          y(coop_default_array_size) = f(xmax)
          !$omp parallel do
          do i=2, coop_default_array_size-1
             y(i) = f(exp(lnxmin + dx*(i-1)))
          enddo
          !$omp end parallel do
       endif
    else
       dx = (xmax - xmin)/(coop_default_array_size - 1)
       if(present(args))then
          y(1) = f(xmin, args)
          y(coop_default_array_size) = f(xmax, args)
          !$omp parallel do
          do i=2, coop_default_array_size-1
             y(i) = f(xmin + dx*(i-1), args)
          enddo
          !$omp end parallel do
       else
          y(1) = f(xmin)
          y(coop_default_array_size) = f(xmax)
          !$omp parallel do
          do i=2, coop_default_array_size-1
             y(i) = f(xmin + dx*(i-1))
          enddo
          !$omp end parallel do
       endif
    endif
    call cf%init(coop_default_array_size, xmin, xmax, y, method = COOP_INTERPOLATE_SPLINE, xlog = cf%xlog, ylog = cf%ylog)
  end function coop_function_constructor

  subroutine coop_function_free(this)
    class(coop_function):: this
    if(allocated(this%f))deallocate(this%f)
    if(allocated(this%f2))deallocate(this%f2)
  end subroutine coop_function_free

  subroutine coop_function_initialize(this, n, xmin, xmax, f, method, fleft, fright, slopeleft, sloperight, chebyshev_order, xlog, ylog)
    class(coop_function):: this
    logical, optional::xlog, ylog
    COOP_INT,intent(IN):: n
    COOP_REAL,intent(IN):: xmin, xmax, f(n)
    COOP_REAL, optional::fleft, fright, slopeleft, sloperight
    COOP_INT, optional::method
    COOP_INT, optional::chebyshev_order
    call this%free()
    if(present(xlog))then
       this%xlog = xlog
    else
       this%xlog = .false.
    endif
    if(present(ylog))then
       this%ylog = ylog
    else
       this%ylog = .false.
    endif
    if(n .le. 1 .or. xmin .eq. xmax)then
       write(*,*) "coop function cannot be initialized for xmin = xmax"
       stop 
    endif
    if(present(method))then
       this%method = method
    else
       this%method = COOP_INTERPOLATE_LINEAR
    endif
    if(this%method .eq. COOP_INTERPOLATE_CHEBYSHEV)then
       if(present(chebyshev_order))then
          this%n = chebyshev_order + 1
       else
          this%n = coop_default_chebyshev_fit_order + 1
       endif
    else
       this%n = n
    endif
    allocate(this%f(this%n), this%f2(this%n))
    if(this%method .ne. COOP_INTERPOLATE_CHEBYSHEV)then
       if(this%ylog)then
          this%f = log(f)
       else
          this%f = f
       endif
    endif
    if(this%xlog)then
       this%xmin = log(xmin)
       this%xmax = log(xmax)
    else
       this%xmin = xmin
       this%xmax = xmax
    endif
    this%dx = (this%xmax - this%xmin)/(this%n-1)
    if(present(fleft))then
       this%fleft = fleft
    else
       this%fleft = f(1)
    endif
    if(present(fright))then
       this%fright = fright
    else
       this%fright = f(n)
    endif
    if(this%ylog)then
       this%fleft = log(this%fleft)
       this%fright = log(this%fright)
    endif
    if(present(slopeleft))then
       this%slopeleft = slopeleft
    else
       this%slopeleft= 0.d0
    endif
    if(present(sloperight))then
       this%sloperight = sloperight
    else
       this%sloperight = 0.d0
    endif
    select case(this%method)
    case(COOP_INTERPOLATE_LINEAR)
       this%f2 = 0.
    case(COOP_INTERPOLATE_QUDRATIC)
       this%f2(2:n-1) = this%f(3:n) + this%f(1:n-2) - 2.*this%f(2:n-1)
       this%f2(1) = this%f2(2)
       this%f2(n) = this%f2(n-1)
       this%f2 = this%f2/6.
    case(COOP_INTERPOLATE_SPLINE)
       call coop_spline_uniform(this%n, this%f, this%f2)
    case(COOP_INTERPOLATE_CHEBYSHEV)
       call coop_chebfit_uniform(n, f, this%n, this%f)
       call coop_chebfit_derv(n, this%xmin, this%xmax, this%f, this%f2)
    end select
  end subroutine coop_function_initialize

  function coop_function_evaluate(this, x) result(f)
    class(coop_function):: this
    COOP_REAL x, f, a, b, xdiff
    COOP_INT l, r
    if(this%xlog)then
       xdiff =  log(x) - this%xmin
    else
       xdiff = x - this%xmin
    endif
    b = xdiff/this%dx + 1.d0
    l = floor(b)
    if(l .lt. 1)then
       f = this%fleft + this%slopeleft*xdiff
    elseif(l.ge. this%n)then
       if(this%xlog)then
          xdiff =  log(x) - this%xmax
       else
          xdiff = x - this%xmax
       endif
       f = this%fright + this%sloperight*xdiff
    else
       select case(this%method)
       case(COOP_INTERPOLATE_LINEAR)
          b = b - l
          f = this%f(l) * (1.d0-b) + this%f(l+1) * b 
       case(COOP_INTERPOLATE_QUDRATIC, COOP_INTERPOLATE_SPLINE)
          b = b - l
          r = l + 1
          a = 1. - b
          f = this%f(l) * a + this%f(r) * b + this%f2(l) * (a**2-1.)*a + this%f2(r)*(b**2-1.)*b
       case(COOP_INTERPOLATE_CHEBYSHEV)
          call coop_chebeval(this%n, this%xmin, this%xmax, this%f, xdiff+this%xmin, f)
       case default
          stop "UNKNOWN interpolation method"
       end select
    endif
    if(this%ylog)then
       f = exp(f)
    endif
  end function coop_function_evaluate

  
  !!simple integration 
  !!to be optimized
  function coop_function_integrate(this, a, b) result(integral)
    class(coop_function)::this
    integer,parameter::nsteps = 8192
    COOP_REAL a, b, integral, dx, y1, y2, ym, x1, x2, xm, lna, lnb, dxby2
    integer i
    if(this%xlog)then
       lna = log(a)
       lnb = log(b)
       dx = (lnb - lna)/nsteps
       dxby2 = dx/2.d0
       x1 = lna
       y1 = this%eval(exp(x1))*exp(x1)
       integral = y1
       do i=1, nsteps-1
          xm = x1 + dxby2
          x2 = x1 + dx
          ym = this%eval(exp(xm))*exp(xm)
          y2 = this%eval(exp(x2))*exp(x2)
          integral = integral + (2.d0*y2+ym*4.d0)
          x1 = x2
          y2 = y1
       enddo
       xm = x1 + dxby2
       x2 = x1 + dx
       ym = this%eval(exp(xm))*exp(xm)
       y2 = this%eval(exp(x2))*exp(x2)
       integral = (integral + (y2+ym*4.d0))*(dx/6.d0)
    else
       dx = (b - a)/nsteps
       dxby2 = dx/2.d0
       x1 = a
       y1 = this%eval(x1)
       integral = y1
       do i=1, nsteps-1
          xm = x1 + dxby2
          x2 = x1 + dx
          ym = this%eval(xm)
          y2 = this%eval(x2)
          integral = integral + (2.d0*y2+ym*4.d0)
          x1 = x2
          y2 = y1
       enddo
       xm = x1 + dxby2
       x2 = x1 + dx
       ym = this%eval(xm)
       y2 = this%eval(x2)
       integral = (integral + (y2+ym*4.d0))*(dx/6.d0)
    endif
  end function coop_function_integrate


!!to be optimized
  function coop_function_derivative(this, x) result(fp)
    class(coop_function)::this
    COOP_REAL x, fp, f, dx
    select case(this%method)
    case(COOP_INTERPOLATE_CHEBYSHEV)
       if(this%xlog)then
          call coop_chebeval(this%n, this%xmin, this%xmax, this%f2, log(x), fp)
          fp = fp/x
       else
          call coop_chebeval(this%n, this%xmin, this%xmax, this%f2, x, fp)
       endif
    case default
       dx = this%dx/4.d0
       if(this%xlog)then
          if(this%ylog)then
             fp =(log(this%eval(x*exp(dx)))-log(this%eval(x*exp(-dx))))/(2.d0*dx*x)
          else
             fp = (this%eval(x*exp(dx))-this%eval(x*exp(-dx)))/(2.d0*dx*x)
          endif
       else
          if(this%ylog)then
             fp =(log(this%eval(x+dx))-log(this%eval(x-dx)))/(2.d0*dx)
          else
             fp =(this%eval(x+dx)-this%eval(x-dx))/(2.d0*dx)
          endif
       endif
    end select
    if(this%ylog)then
       fp = fp * this%eval(x)
    endif
  end function coop_function_derivative


end module coop_function_mod