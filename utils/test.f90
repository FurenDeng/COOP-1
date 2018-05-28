module tmp
  use coop_wrapper_utils
#include "constants.h"      
contains
  function f1(x)
    COOP_REAL::f1,x
    f1 = 1.d0/(1.d0+x**4)
  end function f1

  function f2(x)
    COOP_REAL::f2,x
    f2 = f1(x)*x
  end function f2
end module tmp

program Test
#include "constants.h"    
  use tmp
  use coop_wrapper_utils
  implicit none

  print*, coop_integrate(f2,0.d0,300.d0)/coop_integrate(f1, 0.d0, 300.d0)
end program Test



