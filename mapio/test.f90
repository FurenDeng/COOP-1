program shells
  use coop_healcnn_mod
  implicit none
#include "constants.h"
  type(coop_healcnn)::cnn
  call cnn%init(nlayers = 3, map = "zetaproj/gp_meanchi43_sigmachi1_300_TE.fits", mask="zetaproj/mask.fits",  nmaps= (/ 1, 3, 1 /), nside = (/ 256, 32, 8 /) , nside_pooling = (/ 128, 32 /) )
  write(*,*) cnn%layers(1)%c(1)%lmin
  write(*,*) cnn%layers(1)%c(1)%lmax
  write(*,*) cnn%layers(1)%c(2)%lmin
  write(*,*) cnn%layers(1)%c(2)%lmax
  write(*,*) cnn%layers(1)%c(3)%lmin
  write(*,*) cnn%layers(1)%c(3)%lmax  
  call cnn%fp()
  call cnn%layers(1)%h%write("cnn1.fits")    
  call cnn%layers(2)%h%write("cnn2.fits")  
  call cnn%layers(3)%h%write("cnn3.fits")
end program shells
