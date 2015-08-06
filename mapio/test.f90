program test
  use coop_wrapper_utils
  use coop_healpix_mod
  use coop_fitswrap_mod
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::mask, dustmap
  COOP_SINGLE::mean, sigma
  COOP_REAL::xbar, rms, A
  COOP_INT, dimension(:), allocatable::listpix
  COOP_INT::nlist, i, j
  call mask%read("dust/lat30_mask_n1024.fits")
  call dustmap%read("dust/dust_i_n1024_15a.fits")
  nlist = count(mask%map(:,1).gt.0.5)
  mean = sum(dustmap%map(:,1)*mask%map(:,1))/nlist
  sigma = sqrt(sum((dustmap%map(:,1)-mean)**2*mask%map(:,1))/nlist)
  allocate(listpix(nlist))
  j = 0
  do i=0, mask%npix-1
     if(mask%map(i, 1) .gt. 0.5)then
        j = j + 1
        listpix(j) = i
     endif
  enddo
  call coop_fit_gaussian( dble(dustmap%map(listpix, 1)), 100, xbar, rms, A)
  print*, xbar, rms
  print*, mean, sigma
  call coop_asy_histogram( x = dble(dustmap%map(listpix, 1)), nbins = 100, xlabel = "$\ln I$", ylabel = "$dP/d\ln I$", filename = "lnI_hist.txt", fit_gaussian = .true.) 
end program test
