module tmp
#include "constants.h"      
  use coop_wrapper_utils
  implicit none
contains

  function F(x)
    COOP_REAL::x, F
    if(x.gt. 0.03d0)then
       F = sqrt(1.d0+ x**(-3)) - log(x**1.5d0 + sqrt(1.d0+x**3))/x**3
    else
       F = x**1.5*(2.d0/3.d0 - x**3/5.d0)
    endif
  end function F

  function f2byx(x) !!f(x)^2/x
    COOP_REAL:: x, f2byx
    if(x.gt. 0.03d0)then
       F2byx = (sqrt(1.d0+ x**(-3)) - log(x**1.5d0 + sqrt(1.d0+x**3))/x**3)**2/x
    else
       F2byx = (x*(2.d0/3.d0 - x**3/5.d0))**2
    endif
  end function f2byx

  function fit(x)
    COOP_REAL::x, fit
    COOP_REAL,parameter::x1 = 1.49d0, x2 = 1.51d0
    if(x.lt. x1)then
       fit = x**3*(4.d0/27.d0)/(1.d0+0.3d0*x**3)*(1.d0+x**4.95*0.02)
    elseif(x .gt. x2)then
       fit = log(x+0.6d0/x**0.9) - 0.37123
    else
       fit = ((x1**3*(4.d0/27.d0)/(1.d0+0.3d0*x1**3)*(1.d0+x1**4.95*0.02))*(1.51d0-x) + (log(1.51d0+0.6d0/1.51d0**0.9) - 0.37123)*(x-1.49d0))/0.02d0
    endif
  end function fit
  
end module tmp

program Test
#include "constants.h"    
  use coop_wrapper_utils
  use tmp
  implicit none
  COOP_INT,parameter::n = 151
  COOP_REAL::x(n), y(n)
  COOP_INT::i
  call coop_set_uniform(n, x, 0.2d0, 2.d0)
  do i=1, n
     y(i) = coop_integrate(f2byx, 0.d0, x(i), 1.d-12)
     print*, x(i), y(i)/fit(x(i))-1.d0
  enddo
  
  
end program Test



