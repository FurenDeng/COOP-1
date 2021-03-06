program test
  use coop_wrapper_firstorder
  use coop_coupledDE_collapse_mod
  implicit none
#include "constants.h"
#if DO_COUPLED_DE
  type(coop_dictionary)::dict
  COOP_STRING::params_file
  type(coop_coupledDE_collapse_params)::params
  COOP_INT,parameter::na = 801
  COOP_REAL::zvir1
  COOP_REAL::F_pk, p_nu, e_nu
  COOP_INT::i
  COOP_REAL::a(na), x(12, na)
  type(coop_file)::fp
  COOP_STRING::output

  if(iargc().lt. 1)then
     write(*,*) "========================================================"
     write(*,*) "Syntax:"
     write(*,*) "./CDSolve params.ini"
     write(*,*) "========================================================"
     stop
  endif
  call coop_get_Input(1, params_file)
  call coop_load_dictionary(params_file, dict)
  call coop_dictionary_lookup(dict, "output", output)
  call params%init(dict, update_cosmology = .true.)
  call coop_set_uniform(na, a, 0.03d0, 1.d0)
  call params%get_solution(a, x)
  call fp%open(output)
  write(fp%unit, "(13A16)")  "# a             ", " x1 ", " x2 ", " x3 ", " dot x1 ", " dot x2 ", " dot x3 ", " phi ", " dot phi ", " aH ", " D/a ", " phi_bg ", " dot phi_bg"
  do i=1, na
     write(fp%unit, "(13E16.7)") a(i), x(:, i)
  enddo
  call fp%close()
  write(*,"(A)") "The solution is successfully written to "//trim(output)//"."
  zvir1 = params%zvir1()
  if(zvir1 .ge. 0.d0) then
     write(*,*) "z_collapse = ", zvir1
  else
     write(*,*) "Not collapsed."
  endif
#else
  stop "To use CDSolve you need to compile the COOP with DARK_ENERGY_MODEL  = COUPLED_DE in configure.in"
#endif
end program test
