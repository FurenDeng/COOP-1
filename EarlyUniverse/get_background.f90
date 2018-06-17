program Test
  use coop_wrapper_utils
  use coop_lattice_fields_mod
  implicit none
#include "constants.h"

#define FIELDS bg%y(1:nflds)
#define DOTFIELDS bg%y(nflds+1:2*nflds)
#define LNA bg%y(2*nflds+1)  
#define HUBBLE bg%y(2*nflds+2)

  
  COOP_INT,parameter::max_nstep = 20000
  type(coop_ode)::bg

  !!number of fields in the model, see include/lattice_fields_model.h
  COOP_INT,parameter::nflds = 2

  !!defines the end of inflation: epsilon = 1, where epsilon = - \dot H / H^2
  COOP_REAL,parameter::epsilon_end = 1.d0

  !!other variables
  COOP_REAL::Hini
  COOP_REAL:: phi_ini(nflds), dotphi_ini(nflds), y(nflds*2+2), dt, fluk
  COOP_INT::i, j, nstep
  COOP_REAL,dimension(max_nstep)::lnH_array, eps_array, lna_array, lnH_array2, eps_array2
  COOP_REAL,dimension(:,:,:),allocatable::masses_array, masses_array2
  logical::eps_increase = .false.
  
  !!===============initialize ODE solver============================
  call bg%init(n=2*nflds+2, method=COOP_ODE_DVERK, tol=1.d-10)

  !!=================set inital conditions =========================================
  !!default initial field values; you can change it if your model is different
  phi_ini = 0.d0
  phi_ini(1) = 20.d0*coop_lattice_Mp  !!initial inflaton value
  
  !!approximate, ignore kinetic energy  
  Hini = sqrt(coop_lattice_fields_V(phi_ini))/(coop_sqrt3*coop_lattice_Mp)
  !!slow-roll approximation
  dotphi_ini  = coop_lattice_fields_dVdphi(phi_ini)/(-3.d0*Hini)
  !!now iterate to get more accurate values
  Hini = sqrt(coop_lattice_fields_V(phi_ini) + sum(dotphi_ini**2)/2.d0)/(coop_sqrt3*coop_lattice_Mp)
  dotphi_ini  = coop_lattice_fields_dVdphi(phi_ini)/(-3.d0*Hini)
  Hini = sqrt(coop_lattice_fields_V(phi_ini) + sum(dotphi_ini**2)/2.d0)/(coop_sqrt3*coop_lattice_Mp)
  
  y(1:nflds) = phi_ini
  y(nflds+1:2*nflds) = dotphi_ini
  y(2*nflds+1) = 0.d0  !!initial ln a = 0
  y(2*nflds+2) = Hini
  call bg%set_initial_conditions(xini = 0.d0, yini = y)


  !!======================evolution=================================================
  dt = 0.1/Hini
  nstep = 1
  lna_array(nstep) = LNA
  lnH_array(nstep) = log(HUBBLE)
  eps_array(nstep) =  coop_lattice_background_epsilon(nflds, bg%y)
  if(eps_array(nstep) .gt. 0.2d0) stop "Please set initial conditions during inflation."
  do while( eps_array(nstep) .lt. epsilon_end*0.9999d0 .and. nstep .lt. max_nstep-1)
     if(eps_increase .and. eps_array(nstep) .gt.  epsilon_end*0.5d0 )then
        dt = max(min(0.3 * dt * (epsilon_end-eps_array(nstep))/(eps_array(nstep) - eps_array(nstep-1)), 0.1d0/HUBBLE), 1.d-6/HUBBLE)
     else
        dt = 0.1/HUBBLE
     endif     
     call bg%evolve(coop_lattice_background_eqs, bg%x+dt)

     nstep = nstep + 1
     lna_array(nstep) = LNA
     lnH_array(nstep) = log(HUBBLE)
     eps_array(nstep) =  coop_lattice_background_epsilon(nflds, bg%y)
     eps_increase = (eps_array(nstep) .gt. eps_array(nstep-1) + 1.d-10)
     
     if(abs(coop_lattice_background_rhotot(nflds, bg%y) /(3.d0*coop_lattice_Mpsq*HUBBLE**2) - 1.d0) .gt. 1.d-4) stop "energy conservation failed: check if you have a typo in include/lattice_fields_mode.h"
     
  enddo
  dt = dt * (epsilon_end - eps_array(nstep))/(eps_array(nstep) - eps_array(nstep-1))
  call bg%evolve(coop_lattice_background_eqs, bg%x+dt)
  nstep = nstep + 1
  lna_array(nstep) = LNA
  lnH_array(nstep) = log(HUBBLE)
  eps_array(nstep) =  coop_lattice_background_epsilon(nflds, bg%y)

  lna_array = lna_array - LNA !!set lna = 0 at the end of inflation
  
  call coop_spline(nstep, lna_array(1:nstep), eps_array(1:nstep), eps_array2(1:nstep))
  call coop_spline(nstep, lna_array(1:nstep), lnH_array(1:nstep), lnH_array2(1:nstep))
  allocate(masses_array(nstep, nflds, nflds), masses_array2(nstep, nflds, nflds))
  do i=1, nflds
     do j=1, i
        call coop_spline(nstep, lna_array(1:nstep), masses_array(:, i, j), masses_array2(:, i, j))
     enddo
  enddo
  
  write(*,"(A)") "================================================================================="
  write(*,"(A, E13.3)") "Energy conservation check: relative error = ", coop_lattice_background_rhotot(nflds, bg%y) /(3.d0*coop_lattice_Mpsq*HUBBLE**2) - 1.d0
  write(*, "(A, E16.7)") "H/M_p at the end of inflation:", HUBBLE/coop_lattice_Mp  
  write(*,"(A)") "---------------------------------------------------------------------------------"    
  write(*, "(A)") "Field values (/M_p) at the end of inflation:"
  write(*, "("//COOP_STR_OF(nflds)//"E16.7)") FIELDS/coop_lattice_Mp
  write(*,"(A)") "---------------------------------------------------------------------------------"      
  write(*, "(A)") "Time derivatives of field values (/M_p^2) at the end of inflation:"
  write(*, "("//COOP_STR_OF(nflds)//"E16.7)") DOTFIELDS/coop_lattice_Mpsq
  write(*,"(A)") "---------------------------------------------------------------------------------"

contains

  subroutine interpolate(lna, H, eps, masses)
    COOP_REAL::lna, H, eps, masses(nflds, nflds)
    COOP_INT::i, j
    call coop_splint(nstep, lna_array(1:nstep), eps_array(1:nstep), eps_array2(1:nstep), lna, eps)
    call coop_splint(nstep, lna_array(1:nstep), lnH_array(1:nstep), lnH_array2(1:nstep), lna, H)
    do i=1, nflds
       do j=1, i
          call coop_splint(nstep, lna_array(1:nstep), masses_array(:, i, j), masses_array2(:, i, j), lna, masses(i, j))
          masses(j,i) = masses(i,j)
       enddo
    enddo
    H = exp(H)
  end subroutine interpolate

end program Test
