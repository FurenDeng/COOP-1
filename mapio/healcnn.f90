module coop_healpix_cnn1d_mod
  use coop_wrapper_utils
  use coop_sphere_mod
  use coop_fitsio_mod
  use coop_healpix_mod
  use coop_stacking_mod
  use coop_fitswrap_mod
  use coop_gstack_mod
  implicit none

#include "constants.h"

  type coop_healpix_cnn1d_layer
     COOP_INT::nin = 0
     COOP_INT::nout = 0
     type(coop_healpix_maps)::maps
     type(coop_healpix_maps)::dmaps
     logical::is_pooling = .false.
     COOP_INT::lmax = 200
     COOP_INT::lmin = 15
     COOP_INT::depth = 2
     COOP_REAL,dimension(:,:,:),allocatable::weight
     COOP_REAL,dimension(:,:,:),allocatable::dw     
     COOP_REAL,dimension(:),allocatable::bias
     COOP_REAL,dimension(:),allocatable::db
  end type coop_healpix_cnn1d_layer

  type coop_healpix_cnn1d_binary_selection
     COOP_INT::nlayers
     type(coop_healpix_cnn1d_layer),dimension(:),allocatable::layers
     COOP_REAL::threshold
  end type coop_healpix_cnn1d_binary_selection
  

contains

  
  
end module coop_healpix_cnn1d_mod


