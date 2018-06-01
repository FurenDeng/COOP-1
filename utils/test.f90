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

  subroutine doublearr(x)
    COOP_REAL,dimension(:)::x
    x = x*2
  end subroutine doublearr
end module tmp



program Test
#include "constants.h"    
  use tmp
  use coop_wrapper_utils
  implicit none
  COOP_REAL:: x(2,3)
  x(1,:) = 1
  x(2,:) = 2
  call doublearr(x(2,:))
  write(*,*) x(1,:)
  write(*,*) x(2,:)
end program Test



