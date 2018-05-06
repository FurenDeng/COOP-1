program shells
  use coop_healcnn_mod
  implicit none
#include "constants.h"
  type(coop_healcnn_layer)::layer
  call layer%init(nside_in = 256, nside_pooling = 128, nside_out = 128, nmaps_in = 2, nmaps_out = 3)
  print*, layer%h%nside, layer%h%nmaps, layer%h%lmax
  print*, layer%c(6)%nb
  print*, layer%c(6)%kernel(1, 0.99d0),  layer%c(6)%kernel(2, 0.99d0),  layer%c(6)%kernel(3, 0.99d0),  layer%c(6)%kernel(4, 0.99d0)
end program shells
