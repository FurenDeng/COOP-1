program test
  use coop_wrapper_firstorder
  use coop_forecast_mod
  implicit none
#include "constants.h"
  COOP_UNKNOWN_STRING,parameter::action= "MCMC"
  logical::use_CMB = .true.
  logical::use_BAO = .true.
  logical::use_SN = .true.
  logical::use_HST = .true.
  type(coop_clik_object),target::pl(3)
  type(coop_HST_object),target::HSTlike
  type(coop_data_JLA),target:: jla
  type(coop_bao_object),target::bao(4)
  type(coop_cosmology_firstorder),target::cosmology

  type(coop_mcmc_params)::mcmc
  type(coop_data_pool)::pool
  COOP_STRING::inifile
  COOP_UNKNOWN_STRING, parameter::planckdata_path = "../data/cmb/"
  COOP_UNKNOWN_STRING, parameter::planckdata_path2 =   "/home/zqhuang/includes/planck13/data" 
  COOP_INT i
  COOP_INT,parameter::total_steps = 10000
  COOP_INT,parameter::update_freq = 500
  logical do_update_propose 
  COOP_REAL::loglike
  call coop_MPI_init()

  if(iargc().ge.1)then
     inifile = trim(coop_InputArgs(1))
  else
     write(*,*) "missing ini file"
     write(*,*) "Syntax: "
     write(*,*) "./DOCLIK myinis/xx.ini"
     stop 
  endif

  mcmc%cosmology => cosmology
  call mcmc%init(trim(inifile))
  mcmc%do_flush = .true.  !!do not buffer
  
  if(mcmc%do_fastslow .and. coop_MPI_Rank().eq.0)then
     write(*,*) "doing fast-slow MCMC"
     write(*,*) "number of varying fast parameters:"//COOP_STR_OF(mcmc%n_fast)
     write(*,*) "number of varying slow parameters:"//COOP_STR_OF(mcmc%n_slow)
     if(mcmc%n_slow .eq.0)then
        mcmc%fast_per_round = 100
     endif
  endif


  !!BAO
  if(use_BAO)then
     call bao(1)%init("../data/bao/sdss_6DF_bao.dataset")
     call bao(2)%init("../data/bao/sdss_MGS_bao.dataset")
     call bao(3)%init("../data/bao/sdss_DR11LOWZ_bao.dataset")
     call bao(4)%init("../data/bao/sdss_DR11CMASS_bao.dataset")     
     pool%BAO%baolike => bao
  endif

  !!HST
  if(use_HST) pool%HST%HSTlike => HSTlike

  !!supernova  
  if(use_SN)then
     call jla%read("../data/jla/jla.dataset")
     pool%SN_JLA%JLALike => jla
  endif

  if(use_CMB)then
     if(coop_file_exists(trim(planckdata_path)//"/CAMspec_v6.2TN_2013_02_26_dist.clik"))then
        call pl(1)%init(trim(planckdata_path)//"/CAMspec_v6.2TN_2013_02_26_dist.clik")
        call pl(2)%init(trim(planckdata_path)//"/commander_v4.1_lm49.clik")
        call pl(3)%init(trim(planckdata_path)//"/lowlike_v222.clik")
     else !!try another path
        call pl(1)%init(trim(planckdata_path2)//"/CAMspec_v6.2TN_2013_02_26_dist.clik")
        call pl(2)%init(trim(planckdata_path2)//"/commander_v4.1_lm49.clik")
        call pl(3)%init(trim(planckdata_path2)//"/lowlike_v222.clik")
        
     endif
     pool%CMB%cliklike => pl
  endif


  select case(action)
  case("TEST", "test")
     mcmc%params = mcmc%center
     mcmc%fullparams(mcmc%used) = mcmc%params
     call mcmc%get_lmax_from_data(pool)
     call mcmc%set_cosmology()
     loglike = pool%loglike(mcmc)
     write(*,*) "-ln(likelihood) = ", loglike
  case("MCMC", "mcmc")

     
     do_update_propose = (mcmc%settings%index("propose_matrix") .ne. 0)
     mcmc%do_memsave = do_update_propose

     !!do MCMC
     do i = 1, total_steps
        if(i.lt. 20 .or. mod(i, 10) .eq. 0) print*, "on Node ", coop_MPI_Rank(), ": step", i, " likelihood = ", mcmc%loglike
        if(do_update_propose)then
           if(mod(i, update_freq).eq.0)then
              call mcmc%update_propose()
           endif
           if(i .gt. total_steps/3)then
              do_update_propose = .false.
              call mcmc%chain%init()
              mcmc%do_memsave = .false.
           endif
        endif
        call mcmc%mcmc_step(pool)
     enddo
  case default
     print*, action
     stop "unknown action"
  end select
  call coop_MPI_finalize()
end program test
