program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_REAL::n,s
  read(*,*) n,s
  print*, coop_IncompleteGamma(n/2.d0, s/2.d0)/exp(log_gamma(n/2.d0))

end program Test  
