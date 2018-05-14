program shells
  use coop_hnn_mod
  implicit none
#include "constants.h"
  type(coop_hnn):: hnn
  type(coop_healpix_maps)::map
  COOP_INT::i
  call  map%read("zetaproj/gp_meanchi43_sigmachi1_300_TE.fits", nmaps_wanted = 1)

  call hnn%init( nmaps = (/ 5,  5,  3,  3,  1,  1/), nside = (/ 256,  64,  64, 16, 16, 1/), input = map, delta_ell = 20 )
  
  hnn%true_out = 12.
  do i=0, 500
     call hnn%walk(step = 0.1d0)
     print*, hnn%Err(), hnn%best_Err
  enddo
  
  do i=1, hnn%nlayers
     call hnn%layers(i)%v%write("hnn"//COOP_STR_OF(i)//".fits")
  enddo

end program shells
