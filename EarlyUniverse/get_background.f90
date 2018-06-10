program Test
  use coop_wrapper_utils
  use coop_lattice_fields_mod
  implicit none
#include "constants.h"
  type(coop_ode)::bg
  COOP_INT,parameter::nflds = 2
  COOP_REAL,parameter::stop_at_epsilon = 1.d0
  COOP_REAL::Hubble, phi_ini(nflds), dotphi_ini(nflds), y(nflds*2+1), dlna, eps, lasteps
  COOP_INT::i
  call bg%init(n=2*nflds+1, method=COOP_ODE_DVERK, tol=1.d-10)
  phi_ini = (/ 20.d0*coop_lattice_Mp, 0.d0 /)
  !!approximate, ignore kinetic energy  
  Hubble = sqrt(coop_lattice_fields_V(phi_ini))/(coop_sqrt3*coop_lattice_Mp)
  !!slow-roll approximation
  dotphi_ini  = coop_lattice_fields_dVdphi(phi_ini)/(-3.d0*Hubble)
  !!now iterate to get more accurate values
  Hubble = sqrt(coop_lattice_fields_V(phi_ini) + sum(dotphi_ini**2)/2.d0)/(coop_sqrt3*coop_lattice_Mp)
  dotphi_ini  = coop_lattice_fields_dVdphi(phi_ini)/(-3.d0*Hubble)
  Hubble = sqrt(coop_lattice_fields_V(phi_ini) + sum(dotphi_ini**2)/2.d0)/(coop_sqrt3*coop_lattice_Mp)
  y(1:nflds) = phi_ini
  y(nflds+1:2*nflds) = dotphi_ini
  y(2*nflds+1) = Hubble
  call bg%set_initial_conditions(xini = 0.d0, yini = y)
  eps = sum(bg%y(nflds+1:2*nflds)**2)/2.d0/bg%y(2*nflds+1)**2/coop_lattice_Mpsq
  do while( eps .lt. stop_at_epsilon*0.999d0 )
     dlna = 1.d0/(1.d0/abs(stop_at_epsilon - eps) + 1.d0 )     
     call bg%evolve(coop_lattice_background_eqs, bg%x+dlna)
     lasteps = eps
     eps =  sum(bg%y(nflds+1:2*nflds)**2)/2.d0/bg%y(2*nflds+1)**2/coop_lattice_Mpsq     
  enddo
  dlna = dlna * (stop_at_epsilon - eps)/(eps - lasteps)
  call bg%evolve(coop_lattice_background_eqs, bg%x+dlna)
  eps =  sum(bg%y(nflds+1:2*nflds)**2)/2.d0/bg%y(2*nflds+1)**2/coop_lattice_Mpsq     
  
  write(*,"(6G16.7)") eps, dlna, bg%y(1)/coop_lattice_Mp, bg%y(nflds+1)/coop_lattice_Mpsq, bg%y(2*nflds+1)/coop_lattice_Mp, (coop_lattice_fields_V(bg%y(1:nflds))+sum(bg%y(nflds+1:2*nflds)**2)/2.d0)/(3.d0*coop_lattice_Mpsq*bg%y(2*nflds+1)**2)-1.d0
end program Test
