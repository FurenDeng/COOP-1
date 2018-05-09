program shells
  use coop_hnn_mod
  implicit none
#include "constants.h"
  type(coop_hnn):: hnn
  type(coop_healpix_maps)::map
  COOP_INT::i
  call  map%read("zetaproj/gp_meanchi43_sigmachi1_300_TE.fits", nmaps_wanted = 1)

  call hnn%init( nmaps = (/ 10,  20,  20,  5,  5, 1 /), nside = (/ 256,  256, 64,  64,  16, 16 /), input = map, delta_ell = 20 )
  hnn%true_out = 1.
  do i=0, 100
     call hnn%walk(step = 0.1d0)
     print*, hnn%Err(), hnn%best_Err
  enddo
  


end program shells
