program Test
  use coop_wrapper_utils
  use coop_lattice_fields_mod
  implicit none
#include "constants.h"
#include "lattice.h"
  type(coop_lattice_fields)::this
  COOP_INT::i
  COOP_REAL::dt, phi(2), pi(2), phi_sigma2(2), power_index(2)

  !!initial background; for serious calculations you should get it from inflation.
  phi  =  (/ coop_lattice_Mp, 0.d0 /)
  pi = (/ -coop_lattice_Mpsq*1.d-8, 0.d0 /)
  !!rms fluctuations; for serious calculations you should get it from inflation.
  phi_sigma2 = (/ coop_lattice_Mpsq*1.d-14,  coop_lattice_Mpsq*1.d-14 /)
  power_index = 0.d0
  !!initialize   
  call this%init( n = 32, LH = 10.d0, phi = phi, pi = pi, phi_sigma2 = phi_sigma2, power_index = power_index, use_conformal_time = .true.)
  !!set time step  
  dt = this%dx/40.d0
  !!set initial scale factor and Hubble parameter
  call this%set_pi_y()
  !!choose order for ODE solver (2, 4, or 6)
  this%ode_order = 6.d0
  write(*,"(4A16)") " a ", "E_K/E_tot", " E_G/E_tot ", "  8\pi G H^2/(3\rho)"
  write(*,"(4G16.7)") this%a, this%ke/(this%ke+this%ge+this%pe), this%ge/(this%ke+this%ge+this%pe), this%H**2*3.d0*coop_lattice_Mpsq / (this%ke + this%ge + this%pe) - 1.d0
  !!evolve the fields
  do i = 1, 100
     call this%evolve(dt, 50)
     call this%set_pi_y()
     write(*,"(4G16.7)") this%a, this%ke/(this%ke+this%ge+this%pe), this%ge/(this%ke+this%ge+this%pe), this%H**2*3.d0*coop_lattice_Mpsq / (this%ke + this%ge + this%pe) - 1.d0
  enddo
end program Test
