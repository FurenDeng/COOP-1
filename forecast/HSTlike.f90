module coop_HSTlike_mod
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"

  !!from arxiv:1311.3461
  type coop_HST_object
     COOP_REAL::zeff = 0.04
     COOP_REAL::angconversion = 11425.8d0
     COOP_REAL:: H0 = 70.6
     COOP_REAL:: H0_err = 3.3
   contains
     procedure::LogLike => coop_HST_object_LogLike
  end type coop_HST_object

contains

  function coop_HST_object_LogLike(this, dlzeff) result(logLike)
    class(coop_HST_object)::this
    COOP_REAL Heff, loglike, dlzeff
    Heff = this%angconversion / dlzeff
    Loglike = ((Heff - this%H0)/this%H0_err)**2/2.d0
  end function coop_HST_object_LogLike
  
  
end module coop_HSTlike_mod
