program Test
#include "constants.h"    
  use coop_wrapper_utils
  implicit none
  type(coop_nn)::nn
  COOP_REAL::err1, err2
  COOP_INT::i, j
  call nn%full_init( (/ 192, 32, 8, 4 /) )

  call random_number(nn%layers(1)%v)
  nn%true_out = (/ 1.d0, 2.d0, 3.d0, 4.d0 /)
  do i=0, 1000
     call nn%walk(step = 0.1d0)
     if(mod(i,50).eq.0) print*, nn%Err(), nn%best_Err
  enddo
end program Test



