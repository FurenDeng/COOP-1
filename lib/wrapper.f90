module coop_wrapper
  use coop_wrapper_firstorder
  use camb_mypp
  implicit none
#include "constants.h"

  type(coop_cosmology_firstorder):: coop_global_cosmology
  type(coop_species):: coop_global_baryon
  type(coop_species):: coop_global_cdm
  type(coop_species):: coop_global_radiation
  type(coop_species):: coop_global_massless_neutrinos
  type(coop_species):: coop_global_massive_neutrinos
  type(coop_species):: coop_global_de

  type(coop_arguments)::coop_global_cosmological_parameters, coop_global_param_index

  COOP_REAL,parameter::coop_pp_lnkmin = -9.22d0
  COOP_REAL::coop_pp_lnkmax = -0.3d0



  COOP_INT:: cosmomc_de_index = 8
  COOP_INT:: cosmomc_de2pp_num_params = 7

  COOP_INT:: cosmomc_de_model = COOP_DE_COSMOLOGICAL_CONSTANT
  COOP_INT:: cosmomc_pp_model = COOP_PP_STANDARD
  COOP_INT:: cosmomc_de_num_params = 2
  COOP_INT:: cosmomc_pp_num_params = 8
  COOP_INT:: cosmomc_pp_num_origin = 8
  COOP_INT::  coop_pp_nleft, coop_pp_nright
  COOP_REAL:: coop_pp_lnk_per_knot
  COOP_REAL:: coop_pp_scalar_lnkpivot = log(0.05d0)
  COOP_REAL:: coop_pp_tensor_lnkpivot =  log(0.05d0)
  logical ::cosmomc_pp_inflation_consistency = .true.
  COOP_INT, parameter::coop_pp_lmin = 2
  COOP_INT, parameter::coop_pp_lmax = 2850
  COOP_REAL::coop_pp_ells(coop_pp_lmin:coop_pp_lmax)  
  COOP_REAL::coop_pp_scalar_Cls(coop_num_cls, coop_pp_lmin:coop_pp_lmax)
  COOP_REAL::coop_pp_lensed_Cls(coop_num_cls, coop_pp_lmin:coop_pp_lmax)
  COOP_REAL::coop_pp_tensor_Cls(coop_num_cls, coop_pp_lmin:coop_pp_lmax)
  COOP_REAL::coop_pp_total_Cls(coop_num_cls, coop_pp_lmin:coop_pp_lmax)
  COOP_INT,parameter::coop_pp_n = 1024  
  COOP_REAL,dimension(coop_pp_n)::coop_pp_lnkMpc, coop_pp_lnps, coop_pp_lnpt, coop_pp_lnps2, coop_pp_lnpt2, coop_pp_lneps, coop_pp_lnV, coop_pp_phi, coop_pp_lnH
  COOP_INT::coop_pp_ipivot
 
  interface coop_setup_cosmology_from_cosmomc
     module procedure coop_setup_cosmology_from_cosmomc_s, coop_setup_cosmology_from_cosmomc_d
  end interface coop_setup_cosmology_from_cosmomc

contains

  subroutine coop_setup_cosmology_from_cosmomc_s(params, h, want_firstorder)
    real params(:)
    double precision, optional::h
    logical,optional::want_firstorder
    call COOP_COSMO_PARAMS%init(r = COOP_REAL_OF(params), i = (/ cosmomc_de_model, cosmomc_de_index, cosmomc_de_num_params, cosmomc_pp_model, cosmomc_de_index + cosmomc_de_num_params + cosmomc_de2pp_num_params, cosmomc_pp_num_params /), l = (/ cosmomc_pp_inflation_consistency /)  )

    
    if(COOP_INFLATION_CONSISTENCY)then
       COOP_NT = - COOP_AMP_RATIO / 8.d0
    endif
    COOP_COSMO%As = 1.d-10 * exp(COOP_LN10TO10AS)
    COOP_COSMO%ns = COOP_NS
    COOP_COSMO%nrun = COOP_NRUN
    COOP_COSMO%r =  COOP_AMP_RATIO
    COOP_COSMO%nt = COOP_NT
    COOP_COSMO%inflation_consistency = COOP_INFLATION_CONSISTENCY
    if(present(h))then
       if(present(want_firstorder))then
          call coop_setup_global_cosmology_with_h(COOP_REAL_OF(h), want_firstorder = want_firstorder)
       else
          call coop_setup_global_cosmology_with_h(COOP_REAL_OF(h))
       endif
    endif
  end subroutine coop_setup_cosmology_from_cosmomc_s

  subroutine coop_setup_cosmology_from_cosmomc_d(params, h, want_firstorder)
    doubleprecision params(:)
    double precision, optional::h
    logical,optional::want_firstorder    
    call COOP_COSMO_PARAMS%init(r = COOP_REAL_OF(params), i = (/ cosmomc_de_model, cosmomc_de_index, cosmomc_de_num_params, cosmomc_pp_model, cosmomc_de_index + cosmomc_de_num_params + cosmomc_de2pp_num_params, cosmomc_pp_num_params /), l = (/ cosmomc_pp_inflation_consistency /) )

    
    if(COOP_INFLATION_CONSISTENCY)then
       COOP_NT = - COOP_AMP_RATIO / 8.d0
    endif
    COOP_COSMO%As = 1.d-10 * exp(COOP_LN10TO10AS)
    COOP_COSMO%ns = COOP_NS
    COOP_COSMO%nrun = COOP_NRUN
    COOP_COSMO%r =  COOP_AMP_RATIO
    COOP_COSMO%nt = COOP_NT
    COOP_COSMO%inflation_consistency = COOP_INFLATION_CONSISTENCY
    
    if(present(h))then
       if(present(want_firstorder))then
          call coop_setup_global_cosmology_with_h(COOP_REAL_OF(h), want_firstorder = want_firstorder)          
       else
          call coop_setup_global_cosmology_with_h(COOP_REAL_OF(h))
       endif       
    endif
  end subroutine coop_setup_cosmology_from_cosmomc_d


  subroutine coop_print_info()
    write(*,*) "This COOP Version "//trim(coop_version)
    write(*,*) "Author: Zhiqi Huang"
  end subroutine coop_print_info

  subroutine coop_setup_global_cosmology_with_h(h, want_background, want_firstorder)
    COOP_REAL h
    type(coop_arguments) args
    COOP_INT l
    logical,optional::want_background, want_firstorder
    call COOP_COSMO%free()
    call COOP_COSMO%init(name = "COOP_GLOBAL_COSMOLOGY",  id = 0, h = h)
    if(h.le.0.d0)return  !!return for bad h
    coop_global_baryon = coop_baryon(COOP_OMEGABH2/h**2)
    call COOP_COSMO%add_species(coop_global_baryon)
    coop_global_radiation = coop_radiation(COOP_COSMO%Omega_radiation())
    call COOP_COSMO%add_species(coop_global_radiation)
    if(COOP_MNU .eq. 0.d0)then
       coop_global_massless_neutrinos = coop_neutrinos_massless(COOP_COSMO%Omega_massless_neutrinos())
       call COOP_COSMO%add_species(coop_global_massless_neutrinos)
    else
       coop_global_massless_neutrinos = coop_neutrinos_massless(COOP_COSMO%Omega_massless_neutrinos_per_species()*(COOP_COSMO%NNu()-1))
       call COOP_COSMO%add_species(coop_global_massless_neutrinos)  !!assuming one species
       coop_global_massive_neutrinos = coop_neutrinos_massive( &
            COOP_COSMO%Omega_nu_per_species_from_mnu_eV(COOP_MNU) ,&
            COOP_COSMO%Omega_massless_neutrinos_per_species())
       call COOP_COSMO%add_species(coop_global_massive_neutrinos)
    endif    
    select case(COOP_DE_MODEL)
    case(COOP_DE_COSMOLOGICAL_CONSTANT)
       coop_global_cdm = coop_cdm(COOP_OMEGACH2/h**2)
       call COOP_COSMO%add_species(coop_global_cdm)
       coop_global_de = coop_de_lambda( &
            COOP_COSMO%Omega_k()-COOP_OMEGAK &
            )
       call COOP_COSMO%add_species(coop_global_de)
       COOP_COSMO%de_genre = COOP_PERT_NONE                     
    case(COOP_DE_W0)
       coop_global_cdm = coop_cdm(COOP_OMEGACH2/h**2)
       call COOP_COSMO%add_species(coop_global_cdm)             
       coop_global_de = coop_de_w0( &
            COOP_COSMO%Omega_k()-COOP_OMEGAK, &
            COOP_COSMO_PARAMS%r(COOP_INDEX_DE) &
            )
       call COOP_COSMO%add_species(coop_global_de)
       COOP_COSMO%de_genre = COOP_PERT_NONE                     
    case(COOP_DE_W0WA)
       coop_global_cdm = coop_cdm(COOP_OMEGACH2/h**2)
       call COOP_COSMO%add_species(coop_global_cdm)                    
       coop_global_de = coop_de_w0wa( &
            COOP_COSMO%Omega_k()-COOP_OMEGAK, &
            COOP_COSMO_PARAMS%r(COOP_INDEX_DE), &
            COOP_COSMO_PARAMS%r(COOP_INDEX_DE+1) &
            )
       call COOP_COSMO%add_species(coop_global_de)
       COOP_COSMO%de_genre = COOP_PERT_NONE                     
    case(COOP_DE_QUINTESSENCE)
       coop_global_cdm = coop_cdm(COOP_OMEGACH2/h**2)
       call COOP_COSMO%add_species(coop_global_cdm)                    
       coop_global_de = coop_de_quintessence( &
            COOP_COSMO%Omega_k()-COOP_OMEGAK, &
            COOP_COSMO_PARAMS%r(COOP_INDEX_DE), &
            COOP_COSMO_PARAMS%r(COOP_INDEX_DE+1), &
            COOP_COSMO_PARAMS%r(COOP_INDEX_DE+2) &
            )
       call COOP_COSMO%add_species(coop_global_de)
       COOP_COSMO%de_genre = COOP_PERT_NONE
#if DO_COUPLED_DE       
    case(COOP_DE_COUPLED_QUINTESSENCE)
       stop "setup_global_cosmology_with_h: COUPLED_DE needs to be done"
#endif       
    case default
       stop "UNKNOWN DARK ENERGY MODEL"
    end select
    if(COOP_COSMO%h().le.0.d0)return  !!rejected model
    if(present(want_background))then
       if(want_background) &
            call coop_global_cosmology_setup_background()
    endif
    if(present(want_firstorder))then
       if(want_firstorder) &
            call coop_global_cosmology_setup_firstorder
    endif
  end subroutine coop_setup_global_cosmology_with_h


  subroutine coop_global_cosmology_setup_background()
    call COOP_COSMO%setup_background()        
  end subroutine coop_global_cosmology_setup_background

  subroutine coop_global_cosmology_setup_firstorder()
    type(coop_arguments) args
    COOP_INT l

    if(COOP_COSMO%need_setup_background)call coop_global_cosmology_setup_background()
    COOP_COSMO%optre = COOP_TAU
    call COOP_COSMO%set_xe()
    call COOP_COSMO%compute_source(0)
    if(COOP_COSMO%has_tensor)then
       call COOP_COSMO%compute_source(2)
    endif
    call coop_global_cosmology_compute_Cls()
  end subroutine coop_global_cosmology_setup_firstorder

    subroutine coop_global_cosmology_prepare_from_camb()
    type(coop_arguments) args
    COOP_INT l
    call coop_global_cosmology_setup_background()
    COOP_COSMO%optre = COOP_TAU
    call COOP_COSMO%set_xe()
  end subroutine coop_global_cosmology_prepare_from_camb



  subroutine coop_global_cosmology_compute_Cls(do_lensing)
    logical,save::init = .true.
    logical,optional::do_lensing
    COOP_INT::l
    if(init)then
       do l = coop_pp_lmin, coop_pp_lmax
          coop_pp_ells(l) = dble(l)
       enddo
       init = .false.
    endif
    call coop_setup_pp()
    call COOP_COSMO%source(0)%get_all_Cls(coop_pp_lmin, coop_pp_lmax, coop_pp_scalar_Cls)       
    if(COOP_COSMO%has_tensor)then
       call COOP_COSMO%source(2)%get_all_Cls(coop_pp_lmin, coop_pp_lmax, coop_pp_tensor_Cls)       
    else
       coop_pp_tensor_cls = 0.d0
    endif
    if(present(do_lensing))then
       if(do_lensing)then
          call coop_get_lensing_Cls(coop_pp_lmin, coop_pp_lmax, coop_pp_scalar_Cls, coop_pp_lensed_Cls)
       else
          coop_pp_lensed_Cls  = 0.d0
       endif
    else
       coop_pp_lensed_Cls  = 0.d0
    endif
    coop_pp_total_cls = coop_pp_scalar_cls + coop_pp_lensed_cls + coop_pp_tensor_cls
    
  end subroutine coop_global_cosmology_compute_Cls


  subroutine coop_setup_global_cosmology()
    COOP_REAL,parameter::hmin = 0.5d0, hmax = 0.9d0
    COOP_REAL hl, hr, hm, tl, tr, tm
    hl = hmin
    call coop_setup_global_cosmology_with_h(hl)
    tl = 100.d0*COOP_COSMO%cosmomc_theta()
    hr = hmax
    call coop_setup_global_cosmology_with_h(hr)
    tr = 100.d0*COOP_COSMO%cosmomc_theta()
    if(tl .lt. COOP_100THETA-1.d-8 .and. tr .gt. COOP_100THETA+1.d-8)then
       do while(hr - hl .gt. 0.0001)
          hm = (hl + hr)/2.
          call coop_setup_global_cosmology_with_h(hm)
          tm = 100.d0*COOP_COSMO%cosmomc_theta()
          if(tm .gt. COOP_100THETA+1.d-8)then
             hr  = hm
          elseif(tm .lt. COOP_100THETA - 1.d-8)then
             hl = hm
          else
             return
          endif
       enddo
    else
       write(*,*) hl, tl
       write(*,*) hr, tr
       write(*,*) "Bad cosmological parameters. Setting h = 0."
       call coop_setup_global_cosmology_with_h(COOP_REAL_OF(0.))
    endif
  end subroutine coop_setup_global_cosmology



  subroutine coop_setup_pp()
    COOP_REAL,dimension(:),allocatable::lnk, lnps, lnps2
    COOP_REAL  dlnk
    COOP_INT  coop_pp_nleft, coop_pp_nright, nknots, i
    type(coop_arguments)::args
    select case(COOP_PP_MODEL)
    case(COOP_PP_STANDARD)
       call coop_set_uniform(coop_pp_n, coop_pp_lnkMpc, coop_pp_lnkmin, coop_pp_lnkmax)
       coop_pp_lnps = COOP_LN10TO10AS - 10.d0*coop_ln10 + ( COOP_NS - 1.d0 ) * (coop_pp_lnkMpc - coop_pp_scalar_lnkpivot) + (COOP_NRUN/2.d0) *   (coop_pp_lnkMpc - coop_pp_scalar_lnkpivot) ** 2 + (COOP_NRUNRUN/6.d0) *  (coop_pp_lnkMpc - coop_pp_scalar_lnkpivot) ** 3
    case(COOP_PP_BUMP)
       call coop_set_uniform(coop_pp_n, coop_pp_lnkMpc, coop_pp_lnkmin, coop_pp_lnkmax)
       coop_pp_lnps = COOP_LN10TO10AS - 10.d0*coop_ln10 + ( COOP_NS - 1.d0 ) * (coop_pp_lnkMpc - coop_pp_scalar_lnkpivot) + (COOP_NRUN/2.d0) *   (coop_pp_lnkMpc - coop_pp_scalar_lnkpivot) ** 2 + (COOP_NRUNRUN/6.d0) *  (coop_pp_lnkMpc - coop_pp_scalar_lnkpivot) ** 3 + COOP_COSMO_PARAMS%r(COOP_INDEX_PP+cosmomc_pp_num_origin) * exp(-((coop_pp_lnkMpc - COOP_COSMO_PARAMS%r(COOP_INDEX_PP+cosmomc_pp_num_origin+1))/COOP_COSMO_PARAMS%r(COOP_INDEX_PP+cosmomc_pp_num_origin+2))**2/2)

    case(COOP_PP_SCAN_SPLINE)
       nknots =  COOP_NUM_PP - cosmomc_pp_num_origin + 1
       if(nknots .lt. 5) stop "You need at least 5 knots for scan_spline mode"
       coop_pp_nleft = nint((nknots - 1)* (coop_pp_scalar_lnkpivot-coop_pp_lnkmin) / (-coop_pp_lnkmin))
       coop_pp_nright = nknots - 1 - coop_pp_nleft
       dlnk = (coop_pp_scalar_lnkpivot-coop_pp_lnkmin)/coop_pp_nleft
       coop_pp_lnk_per_knot = dlnk
       coop_pp_lnkmax = coop_pp_lnkmin + (nknots-1)*dlnk
       call coop_set_uniform(coop_pp_n, coop_pp_lnkMpc, coop_pp_lnkmin, coop_pp_lnkmax)       
       allocate(lnk(nknots), lnps(nknots), lnps2(nknots))
       call coop_set_uniform(nknots, lnk, -dlnk*coop_pp_nleft, dlnk*coop_pp_nright)
       lnps(1:coop_pp_nleft) = COOP_COSMO_PARAMS%r(COOP_INDEX_PP+cosmomc_pp_num_origin:COOP_INDEX_PP+cosmomc_pp_num_origin+coop_pp_nleft -1)
       lnps(coop_pp_nleft+1) = 0.d0
       lnps(coop_pp_nleft+2:nknots) = COOP_COSMO_PARAMS%r(COOP_INDEX_PP+cosmomc_pp_num_origin+coop_pp_nleft:COOP_INDEX_PP+COOP_NUM_PP-1)
       call coop_spline(nknots, lnk, lnps, lnps2)
       do i=1, coop_pp_n
          call coop_splint(nknots, lnk, lnps, lnps2, coop_pp_lnkMpc(i) - coop_pp_scalar_lnkpivot, coop_pp_lnps(i))
       enddo
       !!modified to resolve the bump issue
       coop_pp_lnps = coop_pp_lnps + COOP_LN10TO10AS -10.d0*coop_ln10&
            + ( COOP_NS - 1.d0 ) * (coop_pp_lnkMpc - coop_pp_scalar_lnkpivot)
       deallocate(lnk, lnps,  lnps2)
    case(COOP_PP_SCAN_LINEAR)
       nknots =  COOP_NUM_PP - cosmomc_pp_num_origin + 1
       if(nknots .lt. 5) stop "You need at least 5 knots for scan_linear mode"
       coop_pp_nleft = nint((nknots - 1)* (coop_pp_scalar_lnkpivot-coop_pp_lnkmin) / (-coop_pp_lnkmin))
       coop_pp_nright = nknots - 1 - coop_pp_nleft
       dlnk = (coop_pp_scalar_lnkpivot-coop_pp_lnkmin)/coop_pp_nleft
       coop_pp_lnk_per_knot = dlnk       
       coop_pp_lnkmax = coop_pp_lnkmin + (nknots-1)*dlnk
       call coop_set_uniform(coop_pp_n, coop_pp_lnkMpc, coop_pp_lnkmin, coop_pp_lnkmax)       
       allocate(lnk(nknots), lnps(nknots), lnps2(nknots))
       call coop_set_uniform(nknots, lnk, -dlnk*coop_pp_nleft, dlnk*coop_pp_nright)
       lnps(1:coop_pp_nleft) = COOP_COSMO_PARAMS%r(COOP_INDEX_PP+cosmomc_pp_num_origin:COOP_INDEX_PP+cosmomc_pp_num_origin+coop_pp_nleft -1)
       lnps(coop_pp_nleft+1) = 0.d0
       lnps(coop_pp_nleft+2:nknots) = COOP_COSMO_PARAMS%r(COOP_INDEX_PP+cosmomc_pp_num_origin+coop_pp_nleft:COOP_INDEX_PP+COOP_NUM_PP-1)
       do i=1, coop_pp_n
          call coop_linear_interp(nknots, lnk, lnps, coop_pp_lnkMpc(i) - coop_pp_scalar_lnkpivot, coop_pp_lnps(i))
       enddo
       !!modified to resolve the bump issue
       coop_pp_lnps = coop_pp_lnps + COOP_LN10TO10AS - 10.d0*coop_ln10 &
            + ( COOP_NS - 1.d0 ) * (coop_pp_lnkMpc - coop_pp_scalar_lnkpivot)
       deallocate(lnk, lnps, lnps2)
    case(COOP_PP_GENERAL_SINGLE_FIELD)
       stop "not applied yet"
    end select
    call coop_spline(coop_pp_n, coop_pp_lnkMpc, coop_pp_lnps, coop_pp_lnps2)

    if(COOP_AMP_RATIO .gt. 0.d0)then
       select case(COOP_PP_MODEL)
       case(COOP_PP_STANDARD, COOP_PP_SCAN_SPLINE, COOP_PP_SCAN_LINEAR, COOP_PP_BUMP)
          coop_pp_lnpt = log(COOP_AMP_RATIO * coop_primordial_ps(exp(coop_pp_tensor_lnkpivot)))+(COOP_NT)*(coop_pp_lnkMpc - coop_pp_tensor_lnkpivot) + COOP_NTRUN*(coop_pp_lnkMpc - coop_pp_tensor_lnkpivot)**2
       case (COOP_PP_GENERAL_SINGLE_FIELD)
          write(*,*) "COOP_PP_GENERAL_SINGLE_FIELD not applied yet"
          stop          
       case default
          write(*,*) "COOP_PP_MODEL"//COOP_STR_OF(COOP_PP_MODEL)//"not applied yet"
          stop
       end select
       call coop_spline(coop_pp_n, coop_pp_lnkMpc, coop_pp_lnpt, coop_pp_lnpt2)
    else
       coop_pp_lnpt = -50.d0
       coop_pp_lnpt2 = 0.d0
    endif
    call COOP_COSMO%set_power(coop_pp_get_power, args)    
  end subroutine coop_setup_pp


  function coop_primordial_lnps(kMpc)  result(lnps)
    COOP_REAL kMpc, lnps
    call coop_splint(coop_pp_n, coop_pp_lnkMpc, coop_pp_lnps, coop_pp_lnps2, log(kMpc), lnps)
  end function coop_primordial_lnps


  function coop_primordial_ps(kMpc)  result(ps)
    COOP_REAL kMpc, ps
    ps = exp(coop_primordial_lnps(kMpc))
  end function coop_primordial_ps

  function coop_primordial_lnpt(kMpc)  result(lnpt)
    COOP_REAL kMpc, lnpt
    if(COOP_AMP_RATIO .eq. 0.d0)then
       lnpt = -50.d0
    else
       call coop_splint(coop_pp_n, coop_pp_lnkMpc, coop_pp_lnpt, coop_pp_lnpt2, log(kMpc), lnpt)
    endif
  end function coop_primordial_lnpt

  function coop_primordial_pt(kMpc)  result(pt)
    COOP_REAL kMpc, pt
    if(COOP_AMP_RATIO .eq. 0.d0)then
       pt = 0.d0
    else
       pt = exp(coop_primordial_lnpt(kMpc))
    endif
  end function coop_primordial_pt

  subroutine coop_pp_get_potential()
    integer i
    COOP_REAL, parameter::max_delta = 0.4
    COOP_REAL, parameter::max_lneps = log(0.4)
    COOP_REAL dlnk, fourdlnk, dphiby2(coop_pp_n), eps(coop_pp_n), delta(coop_pp_n)
    coop_pp_ipivot = coop_minloc(abs(coop_pp_lnkMpc - coop_pp_scalar_lnkpivot))
    dlnk = coop_pp_lnkMpc(2)-coop_pp_lnkMpc(1)
    fourdlnk = dlnk*4.d0
    do i=2, coop_pp_n - 1
       delta(i) = (coop_pp_lnps(i-1)-coop_pp_lnps(i+1))/fourdlnk
       if(abs(delta(i)) .gt. max_delta)then
          delta(i) = sign(max_delta, delta(i))
       endif
    enddo
    delta(1) = delta(2)
    delta(coop_pp_n) = delta(coop_pp_n  - 1)
    coop_pp_lneps = min(coop_pp_lnpt - coop_pp_lnps - log(16.d0), max_lneps) 
    eps = exp(coop_pp_lneps)
    coop_pp_lneps = min(max_lneps,  coop_pp_lneps - log( (1.d0 - (2.d0*(coop_ln2+coop_EulerC-1.d0))*eps ) /(1.d0-2.d0*eps+(2.d0*(2.d0-coop_EulerC-coop_ln2))*delta))) !!slow-roll correction
    eps = exp(coop_pp_lneps)
    coop_pp_lnH = log(coop_pi2/2.d0*exp(coop_pp_lnpt)/( 1.d0 - (2.d0*(coop_ln2+coop_EulerC-1.d0))*eps ))/2.d0
    coop_pp_lnV = 2.d0*coop_pp_lnH + log(3.d0*(1.d0-eps/3.d0))
    coop_pp_phi(1) = 0.d0
    dphiby2 = sqrt(2.d0*eps)/(1.d0-eps)*dlnk/2.d0
    do i=2, coop_pp_n 
       coop_pp_phi(i) = coop_pp_phi(i-1) + (dphiby2(i-1)+dphiby2(i))
    enddo
    coop_pp_phi = coop_pp_phi - coop_pp_phi(coop_pp_ipivot)
  end subroutine coop_pp_get_potential


  subroutine coop_pp_get_power(kbykpiv, ps, pt, cosmology, args)
    COOP_REAL kbykpiv, ps, pt, kMpc
    type(coop_cosmology_firstorder)::cosmology
    type(coop_arguments)::args
    kMpc = exp(coop_pp_scalar_lnkpivot)*kbykpiv
    ps = coop_primordial_ps(kMpc)
    pt = coop_primordial_pt(kMpc)
  end subroutine coop_pp_get_power


  !!in the coupled DE case, this gives the "effective dark energy" density ratio (@scale factor = a vs. @today).
  !!The "effective dark energy" includes the scalar field and non-conserved part of CDM (rho_CDM - rho_CDM_conserved),  where rho_CDM_conserved is extrapolated from early time when the coupling between CDM and DE can be ignored.
  function coop_global_cosmology_DE_Rhoa4_Ratio_eff(a) result(rhoa4)
    COOP_REAL a, rhoa4, omc_early
    COOP_REAL,parameter::a_early = 1.d-8    
    if(coop_global_de%genre .eq. COOP_SPECIES_COUPLED)then
       omc_early = coop_global_cdm%omega * coop_global_cdm%rhoa3_ratio(a_early) 
       rhoa4 = ( coop_global_de%rhoa3_ratio(a) * coop_global_de%omega + coop_global_cdm%rhoa3_ratio(a) * coop_global_cdm%omega - omc_early ) * a / (coop_global_de%omega +  coop_global_cdm%omega - omc_early) 
    else
       rhoa4 = coop_global_de%rhoa4_ratio(a)
    endif
  end function coop_global_cosmology_DE_Rhoa4_Ratio_eff
  

end module coop_wrapper
