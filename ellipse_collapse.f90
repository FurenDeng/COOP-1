module coop_ellipse_collapse_mod
!  use coop_wrapper_firstorder
  implicit none
!#include "constants.h"
!#define COOP_F_TYPE
!  integer,parameter::coop_short_int_length1 = 2
!  integer,parameter::coop_integer_length1 = kind(1)
!  integer,parameter::coop_real_length1 = kind(1.d0)  ! double precision
!  integer,parameter::coop_string_length1 = 1024
!  integer,parameter::coop_short_string_length1 = 32
#ifndef furen_version

#define furen_version

#endif

#ifndef normalize_zvir1

#define normalize_zvir1

#endif


#ifndef COOP_F_TYPE
#define COOP_F_TYPE

#define COOP_REAL real(8)
#define COOP_INT integer(4)
#define COOP_STRING character(len=1024)
#define COOP_SHORT_STRING character(len=32)
#define COOP_UNKNOWN_STRING character(len=*)
#define COOP_DEALLOC(x)  if(allocated(x))deallocate(x)
!should be delete
#define COOP_INTERPOLATE_LINEAR 1
#define COOP_INTERPOLATE_QUADRATIC 2
#define COOP_INTERPOLATE_SPLINE 3
#define COOP_INTERPOLATE_CHEBYSHEV 4
#define COOP_INTERPOLATE_POLYNOMIAL 5
#define COOP_INTERPOLATE_POWERLAW 6
#define COOP_INTERPOLATE_RATIONAL 7
#define COOP_INTERPOLATE_ZIGZAG 8
#define COOP_INTERPOLATE_NONUNIFORM 9

#endif


  !!by default everything is a function of scale factor a;  a = 1 today.
  !!physical time t is used in the definition of some variables, its unit = 1/H_0 

  !!this is a global accuracy parameter
  !!for high accuracy test you can use something like 1.e-4
  !!for normal runs you can use something ~ 1.e-3
#ifdef furen_version

  COOP_REAL, parameter::coop_ellipse_collapse_accuracy = 1.d-4

#else

  COOP_REAL, parameter::coop_ellipse_collapse_accuracy = 1.d-3

#endif

  !!set z_vir = this value if not collapsed
  COOP_REAL, parameter::coop_ellipse_collapse_bad_zvir = -1.d0
  COOP_REAL,parameter:: coop_pi = acos(-1.d0)
  COOP_REAL,parameter:: coop_pio2 = coop_pi/2.d0

  !import by furen
  type coop_2dfunction
    COOP_SHORT_STRING::name = "NoName2D"
    logical::initialized = .false.
    COOP_INT::nx = 0
    COOP_INT::ny = 0
    COOP_REAL::xmin=0.d0
    COOP_REAL::xmax=1.d0
    COOP_REAL::ymin=0.d0
    COOP_REAL::ymax=1.d0
    COOP_REAL::dx = 0.d0
    COOP_REAL::dy = 0.d0
    logical::xlog = .false.
    logical::ylog = .false.
    logical::zlog = .false.
    COOP_REAL,dimension(:,:),allocatable::f, fxx, fyy
  contains
    !    procedure::init => coop_2dfunction_init
    procedure::init_symmetric => coop_2dfunction_init_symmetric
    procedure::free => coop_2dfunction_free
    procedure::eval_bare => coop_2dfunction_evaluate_bare
    procedure::eval => coop_2dfunction_evaluate
  end type coop_2dfunction

  type coop_function
    COOP_SHORT_STRING::name="NoName"
    logical::initialized = .false.
    logical::is_zero = .false.
    COOP_INT::method = COOP_INTERPOLATE_LINEAR
    logical::xlog =.false.
    logical::ylog = .false.
    logical::check_boundary = .true.
    COOP_INT::n = 0
    COOP_INT::n_down = 0
    COOP_REAL::scale = 1.d0   !!f -> scale * f + shift
    COOP_REAL::shift = 0.d0
    COOP_REAL xmin, xmax, dx, fleft, fright, slopeleft, sloperight
    COOP_REAL, dimension(:), allocatable::f, f1, f2, f3
  contains
    !    procedure::set_boundary => coop_function_set_boundary
    procedure::init => coop_function_init
    !    procedure::mult_const => coop_function_mult_const
    !    procedure::add_const => coop_function_add_const
    procedure::init_polynomial => coop_function_init_polynomial
    !    procedure::init_rational => coop_function_init_rational
    !    procedure::init_powerlaw => coop_function_init_powerlaw
    !    procedure::init_NonUniform => coop_function_init_NonUniform
    !    procedure::init_zigzag => coop_function_init_zigzag
    procedure::eval => coop_function_evaluate
    procedure::eval_bare => coop_function_evaluate_bare !!without log scaling
    !    procedure::derivative_bare => coop_function_derivative_bare !!without log scaling
    !    procedure::derivative2_bare => coop_function_derivative2_bare !!without logc scaling
    !    procedure::integrate => coop_function_integrate
    !    procedure::derivative => coop_function_derivative
    !    procedure::derivative2 => coop_function_derivative2
    !    procedure::maxloc => coop_function_maxloc
    !    procedure::minloc => coop_function_minloc
    procedure::free => coop_function_free
    !    procedure::monotonic_solution => coop_function_monotonic_solution
  end type coop_function
  !end

  type coop_ellipse_collapse_params
     COOP_REAL,dimension(3)::lambda = (/ 0.d0, 0.d0, 0.d0 /) !!lambda's
     COOP_REAL::Omega_m = 0.3d0  !!fractional matter density
     !!dark energy EOS  w+ wa(1-a)
     COOP_REAL::w = -1.d0
     COOP_REAL::wa = 0.d0
     !!if set to nonzero, ignore w, wa and use HBK one-parameter parametrization (Huang, Bond, Kofman 2011, ApJ)
     COOP_REAL::epsilon_s = 0.d0
     COOP_REAL::Omega_r, Omega_de
     COOP_REAL::Omega_k = 0.d0
     COOP_REAL::h = 0.7d0
     COOP_REAL::T_CMB = 0.d0 !2.726  !!set to zero to turn off radiation effects, for comparison with peakpatch
     logical::is_spherical = .false.
     COOP_INT::num_ode_vars = 6
     COOP_REAL,dimension(3)::collapse_a_ratio = (/ 0.18d0, 0.18d0, 0.18d0 /)
     type(coop_2dfunction)::bprime
     type(coop_function)::Dbya, Dbyadot

     !defined by furen deng=============================================================================================
     type(coop_function)::chi_tab,omde_tab,t_tab,H_tab,HBK_w_tab
     COOP_REAL::amin = 1.d-7, amax = 1.0, astart = 1.d-3
     COOP_INT::nmax = 2000


#ifdef normalize_zvir1
     COOP_REAL::z_norm_zvir1 = 1100.d0, D_norm_zvir1
#endif
     !end===============================================================================================================

   contains
     procedure::free => coop_ellipse_collapse_params_free  !!subroutine; no argument; release the memory allocated for this object
     procedure::init => coop_ellipse_collapse_params_init  !!subroutine; this%init(Omega_m, w, Omega_k, h, F_pk, e_nu, p_nu);  initialize the object with these parameters
     procedure::get_bprime => coop_ellipse_collapse_params_get_bprime !!subroutine; this%get_bprime(x, bprime)  where x(1:3) is the input array x_1, x_2, x_3, bprime(1:3) is the output b'_1, b'_2, b'_3;
     procedure::Growth_D => coop_ellipse_collapse_params_Growth_D !! this%Growth_D(a) is a function that returns the growth factor D; input a is the scale factor
     procedure::Growth_H_D => coop_ellipse_collapse_params_Growth_H_D  !!this%Growth_H_D(a) is a function that returns d ln D/d t, where a is the scale factor
     procedure::dadt => coop_ellipse_collapse_params_aH  !!this%dadt(a) is a function that  returns da/dt
     procedure::ddotabya => coop_ellipse_collapse_params_ddotabya !!this%ddotabya is a function that returns ( \ddot a / a )
     procedure::set_initial_conditions => coop_ellipse_collapse_params_set_initial_conditions  !!this%set_initial_conditions(y) is a subroutine set the initial conditions for y = (x_1, x_2, x_3, d x_1/dt, d x_2/dt, d x_3/dt)
     procedure::evolve=> coop_ellipse_collapse_params_evolve  !!this%evolve(a, y, a_end) is a subroutine that evolves the current  y = (x_1, x_2, x_3, d x_1/dt, d x_2/dt, d x_3/dt) from current scale factor a to the scale factor a_end; the input is y at a; after the call a becomes a_end.
     procedure::get_solution => coop_ellipse_collapse_params_get_solution  !!this%get_solution(a_arr, x_arr) is a subroutine; input a_arr(1:n) is an array of scale factors in ascending order; return the result in x_arr(1:3, 1:n), where x_arr(1:3, i) is the solution of (x_1, x_2, x_3) at scale factor a_arr(i).
     procedure::zvir1 => coop_ellipse_collapse_params_zvir1  !!this%zvir1() is a function that returns zvir1
     procedure::set_growth_initial_conditions => coop_ellipse_collapse_params_set_growth_initial_conditions !!set initial conditions for the growth function solver
     procedure::make_tab => coop_make_table_for_pkp
     procedure::pkp_init => pkp_init_COOP
!     procedure::chi_t => ode_forw_chi_t
     procedure::init_tab => set_initial_condition_for_table
  end type coop_ellipse_collapse_params


contains

  subroutine coop_ellipse_collapse_params_set_initial_conditions(this, a_ini, y)
    !!set inital vector y = ( x_1, x_2,  x_3,  d x_1/dt, d x_2/dt, d x_3/dt ) at a = a_ini
    class(coop_ellipse_collapse_params)::this
    COOP_REAL,intent(IN)::a_ini
    COOP_REAL::y(:), D_ini,  dadt_ini, corr(3), suml, suml2
    D_ini = this%Growth_D(a_ini)
    dadt_ini = this%dadt(a_ini)
    suml = sum(this%lambda)
    suml2 = sum(this%lambda**2)
    corr = (1.2d0*this%lambda*suml + 0.6d0*suml2 + 11.d0/35.d0*(suml**2-suml2)-3.d0*this%lambda**2)/10.d0  !!this is for 2nd-order correction of initial conditions; in spherical case corr = 3/7 lambda^2
    y(1:3) = a_ini *(1.d0-this%lambda*D_ini - corr*D_ini**2)
    y(4:6) = y(1:3)*dadt_ini/a_ini - a_ini*(D_ini*this%Growth_H_D(a_ini)*(this%lambda + 2.d0*corr*D_ini))
  end subroutine coop_ellipse_collapse_params_set_initial_conditions

  subroutine coop_ellipse_collapse_odes(n, a, y, dyda, params)
  !!the ODE that evolves  y = ( x_1, x_2,  x_3,  d x_1/dt, d x_2/dt, d x_3/dt )
  !!other inputs: 
  !!n = 6 is the dimension of y
  !!a is the scale factor
  !!params is the object containing all the parameters and methods
  !! return dyda = d y/d a
    COOP_INT::n
    COOP_REAL::a, y(n), dyda(n), dadt, bprime(3), growthD, delta, dark_Energy_term, rhomby3, radiation_term, delta_plus_1
    type(coop_ellipse_collapse_params)::params
    COOP_REAL,parameter::eps = coop_ellipse_collapse_accuracy
    COOP_REAL::suppression_factor, arat(3), aeq
    dadt = params%dadt(a)
    radiation_term = - params%omega_r/a**4*2.d0
    !!dark energy contribution; ignore dark energy perturbations in wCDM
    if(params%epsilon_s .ne. 0)then
       aeq = (params%Omega_m/params%omega_de)**(1.d0/(3.d0-1.08*(1.d0-params%omega_m)*params%epsilon_s))  !!see Huang, Bond, Kofman, 2011 ApJ

       dark_Energy_term = - params%omega_de*fit_HBK_rho_DE_ratio(params%epsilon_s, aeq, a)*(1.d0+3.d0*HBK_w(params%epsilon_s, aeq, a))
    else
       dark_Energy_term =  - params%omega_de*a**(-3.d0*(1.d0+params%w+params%wa))*exp(-3.d0*params%wa*(1.d0-a))*(1.d0+3.d0*(params%w+params%wa*(1.d0-a)))
    endif
    if(params%is_spherical)then
       arat(1) = (y(1)/a/params%collapse_a_ratio(1) - 1.d0)/eps
       if(arat(1) .lt. -1.d0 .and. y(4) .lt. 0.d0)then  !!collapsed; freeze it
          dyda(1) = 0.d0
          dyda(4) = 0.d0
       else  !!still collapsing
          !!----------- equation for x_1 -------------------------
          dyda(1) = y(4) / dadt    !!d x_1/da = d x_1/dt /(da/dt)
          dyda(4) = y(1)/dadt/2.d0 * ( &   !! d( dx_1/dt)/da =(d^2 x_1/dt^2)/(da/dt)
               dark_Energy_term  &
               + radiation_term &
               -  params%Omega_m/(y(1))**3 &  !!matter contribution
               )
          !!-----------end of equation for x_1 -------------------------
          if(arat(1) .lt. 0.d0  .and. y(4) .lt. 0.d0)then  !!do suppression around x_i/a = fr_i so that the derivative is continuous; no need to change
             suppression_factor =  sin(coop_pio2*(1.d0+arat(1)))**4
             dyda(1) = dyda(1)*suppression_factor
             dyda(4) = dyda(4)*suppression_factor

          endif
       endif
       dyda(2:3) = dyda(1)
       dyda(5:6) = dyda(4)
    else
       rhomby3 = params%Omega_m/a**3 !!I am working in unit of H_0^2/(8\pi G)
       delta = a**3/(y(1)*y(2)*y(3))-1.d0
       call params%get_bprime(y(1:3), bprime)
       growthD = params%growth_D(a)
       arat = (y(1:3)/a/params%collapse_a_ratio - 1.d0)/eps
       if(arat(1) .lt. -1.d0  .and. y(4) .lt. 0.d0)then  !!collapsed; freeze it
          dyda(1) = 0.d0
          dyda(4) = 0.d0
       else  !!still collapsing
          !!----------- equation for x_1 -------------------------
          dyda(1) = y(4) / dadt    !!d x_1/da = d x_1/dt /(da/dt)
          dyda(4) = y(1)/dadt/2.d0 * ( &   !! d( dx_1/dt)/da =(d^2 x_1/dt^2)/(da/dt)
               dark_Energy_term  &
               + radiation_term &
               -  rhomby3*(1.d0 + delta *(1.d0+bprime(1)*1.5d0) + (3.d0*params%lambda(1)-sum(params%lambda))*growthD ) &  !!matter contribution
               )
          !!----------- end of equation for x_1 -------------------------
          if(arat(1) .lt. 0.d0)then !!do suppression around x_i/a = fr_i so that the derivative is continuous; this is for stability of the ode solver;
             suppression_factor =  sin(coop_pio2*(1.d0+arat(1)))**4
             dyda(1) = dyda(1)*suppression_factor
             dyda(4) = dyda(4)*suppression_factor
          endif
       endif
       if(arat(2).lt.-1.d0 .and. y(5) .lt. 0.d0 )then !!collapsed; freeze it
          dyda(2) = 0.d0
          dyda(5) = 0.d0
       else !!still collapsing
          !!----------- equation for x_2 -------------------------
          dyda(2) = y(5) / dadt    !!d x_1/da = d x_1/dt /(da/dt)
          dyda(5) = y(2)/dadt/2.d0 * ( &   !! d( dx_1/dt)/da =(d^2 x_1/dt^2)/(da/dt)
               dark_Energy_term  &
               + radiation_term &
               - rhomby3 *(1.d0 + delta *(1.d0+bprime(2)*1.5d0) + (3.d0*params%lambda(2)-sum(params%lambda))*growthD ) &  !!matter contribution
               )
          !!----------- end of equation for x_2 -------------------------
          if(arat(2) .lt. 0.d0)then !!do suppression around x_i/a = fr_i so that the derivative is continuous; this is for stability of the ode solver;
             suppression_factor =  sin(coop_pio2*(1.d0+arat(2)))**4
             dyda(2) = dyda(2)*suppression_factor
             dyda(5) = dyda(5)*suppression_factor
          endif
       endif
       if(arat(3).lt.-1.d0 .and. y(6) .lt. 0.d0)then !!collapsed; freeze it
          dyda(3) = 0.d0
          dyda(6) = 0.d0
       else !!still collapsing
          !!----------- equation for x_3 -------------------------
          dyda(3) = y(6) / dadt    !!d x_1/da = d x_1/dt /(da/dt)
          dyda(6) = y(3)/dadt/2.d0 * ( &   !! d( dx_1/dt)/da =(d^2 x_1/dt^2)/(da/dt)
               dark_Energy_term  &
               + radiation_term &
               -  rhomby3*(1.d0 + delta *(1.d0+bprime(3)*1.5d0) + (3.d0*params%lambda(3)-sum(params%lambda))*growthD ) &  !!matter contribution
               )
          !!----------- end of equation for x_3 -------------------------
          if(arat(3) .lt. 0.d0)then  !!do suppression around x_i/a = fr_i so that the derivative is continuous; this is for stability of the ode solver;
             suppression_factor =  sin(coop_pio2*(1.d0+arat(3)))**4
             dyda(3) = dyda(3)*suppression_factor
             dyda(6) = dyda(6)*suppression_factor
          endif
       endif
    endif
  end subroutine coop_ellipse_collapse_odes


!!y = ( D/a,  d (D/a)/dt )
!!set initial conditions at a_ini
  subroutine coop_ellipse_collapse_params_set_growth_initial_conditions(this, a_ini, y)
    class(coop_ellipse_collapse_params)::this
    COOP_REAL::a_ini, dadt, addbya, H
    COOP_REAL::y(2)
    dadt = this%dadt(a_ini)
    H = dadt/a_ini
    addbya = this%ddotabya(a_ini)
    y(1) = 1.d0
    y(2) =  - (2.d0*H**2+addbya-1.5d0*this%Omega_m/a_ini**3)*y(1)/4.d0/H
  end subroutine coop_ellipse_collapse_params_set_growth_initial_conditions


!!y = ( D/a,  d (D/a)/dt )
!!the equation in the standard wCDM case is
!!\ddot (aD) + 4 H (aD) + (2H^2 - \ddot a/ a - 3/2 Omega_m0 / a^3) (aD) = 0
!!n is the dimension of y (here = 2)
!!a is the scale factor
!!dyda returns dy/da
!!params is the object that contains all the parameters
  subroutine coop_ellipse_collapse_growth_ode(n, a, y, dyda, params)
    type(coop_ellipse_collapse_params)::params
    COOP_INT::n
    COOP_REAL::a, dadt, addbya, H
    COOP_REAL::y(2), dyda(2)
    dadt = params%dadt(a)
    H = dadt/a
    addbya = params%ddotabya(a)
    dyda(1) = y(2)/dadt  !! d(D/a)/ da  = d(D/a)/dt / (da/dt)
    dyda(2) = (-4.d0*H*y(2) - (2.d0*H**2+addbya-1.5d0*params%Omega_m/a**3)*y(1))/dadt
  end subroutine coop_ellipse_collapse_growth_ode



!!=====================You don't need to read anything below ==================

!!==================== methods of class coop_ellipse_collapse_params ============
  subroutine coop_ellipse_collapse_params_evolve(this, a, y, a_end)
  !!evolve y from a to a_end
  !!vector y = ( x_1, x_2,  x_3,  d x_1/dt, d x_2/dt, d x_3/dt )
  !!the object this contains all the parameters for the model
    class(coop_ellipse_collapse_params)::this
    COOP_REAL::a, y(this%num_ode_vars), a_end
    COOP_REAL,parameter::tol = max(1.d-10, coop_ellipse_collapse_accuracy*1.d-5)
    COOP_INT::ind
    COOP_REAL::c(24), w(this%num_ode_vars, 9)
    select type(this)
    type is(coop_ellipse_collapse_params)
       ind = 1
       call coop_dverk_with_ellipse_collapse_params(this%num_ode_vars, coop_ellipse_collapse_odes, this, a, y, a_end, tol, ind, c, this%num_ode_vars, w)
    class default
       stop "Evolve: Extended class this has not been implemented"
    end select
  end subroutine coop_ellipse_collapse_params_evolve

  subroutine coop_ellipse_collapse_params_get_solution(this, a_arr, x_arr)
    !!get multiple snapshots of of x
    !!input a_arr(1:n) is an array of scale factors in ascending order;
    !! return the result in x_arr(1:3, 1:n), where x_arr(1:3, i) is the solution of (x_1, x_2, x_3) at scale factor a_arr(i).
    class(coop_ellipse_collapse_params)::this
    COOP_REAL::a_arr(:), x_arr(:,:), y(6), a
    COOP_REAL,parameter::tol = max(1.d-10, coop_ellipse_collapse_accuracy*1.d-5)
    COOP_INT::ind, i, n, m
    COOP_REAL::c(24), w(this%num_ode_vars, 9)
    n = size(a_arr)
    if(size(x_arr, 2) .ne. n) stop "get_solution: input a_arr and x_arr have different sizes"
    m = min(size(x_arr, 1), 6)  !!normally m = 3 but not always
    select type(this)
    type is(coop_ellipse_collapse_params)
       ind = 1
       !modified by furen
       a = min(max(0.005d0, min(0.05d0, 100.d0*coop_ellipse_collapse_accuracy))/maxval(abs(this%lambda)), a_arr(1)*0.99d0, 0.03d0)
!       a=1.0/(20.0*sum(this%lambda))
       !end

       call this%set_initial_conditions(a, y)
       do i=1, n
          call coop_dverk_with_ellipse_collapse_params(this%num_ode_vars, coop_ellipse_collapse_odes, this, a, y, a_arr(i), tol, ind, c, this%num_ode_vars, w)
          x_arr(1:m, i) = y(1:m)
          if(size(x_arr,1).ge.7)then
             x_arr(7, i) = this%dadt(a_arr(i))*sqrt(a_arr(i))
          endif
          if(size(x_arr,1).ge.8)then
             x_arr(8, i) = this%Growth_D(a_arr(i))/a_arr(i)
          endif
       enddo
    class default
       stop "Evolve: Extended class this has not been implemented"
    end select
  end subroutine coop_ellipse_collapse_params_get_solution

!!return the redshift where the last virialized axis collapses==========================================================
  function coop_ellipse_collapse_params_zvir1(this) result(zvir1)
    class(coop_ellipse_collapse_params)::this
    COOP_REAL::a, a_last, a_next, y(this%num_ode_vars), Frho, ycopy(this%num_ode_vars), incr, zvir1
    COOP_REAL,parameter::tol = max(1.d-10, coop_ellipse_collapse_accuracy*1.d-5)
    COOP_INT::ind
    COOP_REAL::c(24), w(this%num_ode_vars, 9)
    COOP_INT::indcopy
    COOP_REAL::ccopy(24), wcopy(this%num_ode_vars, 9)
    type(coop_ellipse_collapse_params)::params_D_normal

    select type(this)
    type is(coop_ellipse_collapse_params)
       ind = 1
#ifdef normalize_zvir1
!       if(this%epsilon_s .ne. 0.d0 .or.this%w .ne. -1.d0 .or. this%wa .ne. 0.d0)then
!         call params_D_normal%pkp_init(Omega_m = this%Omega_m, Omega_k = this%Omega_k, h = this%h, w = -1.d0, epsilon_s = 0.d0, wa=0.d0, &
!                                       amin=this%amin, amax=this%amax, nmax=this%nmax, astart=this%astart)
!         this%lambda = this%lambda/params_D_normal%Growth_D(1.d0/1101.d0)*this%Growth_D(1.d0/1101.d0)
!       endif
       this%lambda = this%lambda*this%Growth_D(1.d0/(this%z_norm_zvir1 + 1.d0))/this%D_norm_zvir1

#endif
       Frho = sum(this%lambda)
       if(Frho .lt. 1.d0 .or. this%lambda(1) .gt. this%lambda(2) .or. this%lambda(2) .gt. this%lambda(3) .or. this%lambda(1) .lt. -50.d0 )then !!bad input or obviously won't collapse
          zvir1 = -1.d0
          return
       endif
       a = min(max(0.005d0, min(0.05d0, 50.d0*coop_ellipse_collapse_accuracy)/maxval(abs(this%lambda))), 0.03d0)
       call this%set_initial_conditions(a, y)
       a_next = 0.1d0/Frho
       call coop_dverk_with_ellipse_collapse_params(this%num_ode_vars, coop_ellipse_collapse_odes, this, a, y, a_next, tol/10.d0, ind, c, this%num_ode_vars, w)
       incr = 0.05d0
       do
          a_last = a
          indcopy = ind
          ccopy = c
          wcopy = w
          ycopy = y


#ifdef furen_version
          a_next = min(a*(1.d0+incr), this%amax)
#else
          a_next = min(a*(1.d0+incr), 1.d0)
#endif
          call coop_dverk_with_ellipse_collapse_params(this%num_ode_vars, coop_ellipse_collapse_odes, this, a, y, a_next, tol, ind, c, this%num_ode_vars, w)
          !use ccopy to make the value of all variables not change between to steps
          if(y(1)/a .lt.  this%collapse_a_ratio(1)  .and. y(4) .le. 0.d0)then
             if(incr .lt. coop_ellipse_collapse_accuracy)then
                zvir1 = 2.d0/(a+a_last) - 1.d0

                return
             endif
             a = a_last
             ind = indcopy
             c = ccopy
             w = wcopy
             y = ycopy
             incr = incr/2.d0
             cycle
          endif

#ifdef furen_version
          if(a_next .gt. this%amax-coop_ellipse_collapse_accuracy/2.d0)then
#else
          if(a_next .gt. 1.d0-coop_ellipse_collapse_accuracy/2.d0)then
#endif
             zvir1 =  coop_ellipse_collapse_bad_zvir
             return
          endif
       enddo
    class default
       stop "Evolve: Extended class this has not been implemented"
    end select
  end function coop_ellipse_collapse_params_zvir1

  !add by furen deng====================================================================================================
  subroutine coop_make_table_for_pkp(this, Omega_m, w, wa, epsilon_s, h, Omega_k, lambda, F_pk, e_nu, p_nu,amin, amax, astart, nmax, outfile)!element in a table is nmax+1

      !variable by furen
      COOP_INT::nmax
      COOP_REAL::amin,amax,dlnD_dlnatab(nmax+1),chitab(nmax+1),ttab(nmax+1),Dtab(nmax+1),omdetab(nmax+1),Htab(nmax+1),HBK_tab(nmax+1),astart
      COOP_REAL::y_chit(2),adynamic_t
      COOP_REAL::c_t(24), wspace_t(2, 9),anow,tnow,taunow,ainterp,Dnow,HD_Ha_now,aeq
      logical::now=.False.,full_evolution=.False.
      COOP_INT::ind_t,i0,istart
      COOP_STRING::outfile

      !defined by COOP a(na), Dbya(na), Dbyadot(na)
      !end

      class(coop_ellipse_collapse_params)::this
      type(coop_ellipse_collapse_params)::params_norm
      COOP_INT::na
      COOP_REAL::a(nmax+1), Dbya(nmax+1), Dbyadot(nmax+1), adynamic
      COOP_INT::i
      COOP_REAL,optional::Omega_m, w, wa, epsilon_s, Omega_k, h, lambda(3), F_pk, e_nu, p_nu
      logical cosmology_updated
      !!!!for D(a) solver
      COOP_REAL::y(2)
      COOP_REAL,parameter::tol = 1.d-8
      COOP_INT::ind
      COOP_REAL::c(24), wspace(2, 9)
      COOP_REAL::amin_D
      COOP_REAL::amax_D

      cosmology_updated = .false.
      !add by furen deng

      this%nmax=nmax
      na = nmax+1
      amin_D = amin
      amax_D = amax
      this%amax = amax
      this%amin = amin
      this%astart = astart
      !end

      if(present(Omega_m))then
          if(abs(this%Omega_m - Omega_m) .gt. 1.d-5)then
              this%Omega_m = Omega_m
              cosmology_updated = .true.
              write(*,*)'update Omega_m: ',Omega_m
          endif
      endif
      if(present(Omega_k))then
          if(abs(this%Omega_k - Omega_k) .gt. 1.d-5)then
              this%Omega_k = Omega_k
              cosmology_updated = .true.
              write(*,*)'update Omega_k: ',Omega_k
          endif
      endif

      if(present(w))then
          if(abs(this%w - w) .gt. 1.d-5)then
              this%w = w
              cosmology_updated = .true.
              write(*,*)'update w: ',w
          endif
      endif
      if(present(wa))then
          if(abs(this%wa - wa) .gt. 1.d-5)then
              this%wa = wa
              cosmology_updated = .true.
              write(*,*)'update wa: ',wa
          endif
      endif
      if(present(epsilon_s))then
          if(abs(this%epsilon_s - epsilon_s) .gt. 1.d-5)then
              this%epsilon_s = epsilon_s
              cosmology_updated = .true.
              write(*,*)'update epsilon_s: ',epsilon_s
          endif
      endif
      if(abs(this%epsilon_s) .gt. 1.d-5 .and. (abs(this%w+1.d0).gt.1.d-5 .or. abs(this%wa) .gt. 1.d-5))then
          write(*,*) "Warning: simultaneously setting w, wa, and epsilon_s. Use epsilon_s by default."
      endif
      if(present(h))then
        if(abs(this%h - h) .gt. 1.d-5)then
          this%h = h
          cosmology_updated = .true.
          write(*,*)'update h: ',h
        endif
      endif

#ifdef normalize_zvir1

      if(cosmology_updated.and.(abs(this%epsilon_s) .gt. 1.d-5 .or. abs(this%w+1.d0).gt.1.d-5 .or. abs(this%wa) .gt. 1.d-5))then

        call params_norm%init(Omega_m = this%Omega_m, Omega_k = this%Omega_k, h = this%h, w = -1.d0, epsilon_s = 0.d0, wa=0.d0)
        this%D_norm_zvir1 = params_norm%Growth_D(1.d0/(this%z_norm_zvir1 + 1.d0))
        write(*,*)'normalize to z = ',this%z_norm_zvir1

      end if
#endif

      this%Omega_r = 4.187d-5/this%h**2*(this%T_CMB/2.726)**4
      this%Omega_de = 1.d0 - this%Omega_m - this%Omega_r - this%Omega_k
      if(present(lambda))then
          this%lambda = lambda
          if(present(F_pk) .or. present(p_nu) .or. present(e_nu))then
              stop "You pass either lambda or (F_pk, e_nu, p_nu) to init, not both."
          endif
      elseif(present(F_pk) .and. present(p_nu) .and. present(e_nu))then
          if(F_pk .lt. 0.d0 .or. e_nu .lt. 0.d0 .or. abs(p_nu) .gt. e_nu)then
              write(*,*) "Invalid F_pk, e_nu, p_nu:"
              write(*,*) F_pk, e_nu, p_nu
              write(*,*) "check conditions: F_pk >=0, e_nu >= 0, and |p_nu| <= e_nu"
              stop
          endif
          this%lambda(3) = (F_pk/3.d0)*(1.d0 + 3.d0*e_nu + p_nu)
          this%lambda(2) = (F_pk/3.d0)*(1.d0 - 2.d0*p_nu)
          this%lambda(1) = (F_pk/3.d0)*(1.d0 - 3.d0*e_nu + p_nu)
      endif
      this%is_spherical = abs(this%lambda(1) - this%lambda(2)) .lt. 1.d-8 .and. abs(this%lambda(1) - this%lambda(3)) .lt. 1.d-8 .and. abs(this%collapse_a_ratio(1) - this%collapse_a_ratio(2)) .lt. 1.d-3 .and. abs(this%collapse_a_ratio(1)-this%collapse_a_ratio(3)) .lt. 1.d-3
      !!set b' function
      if(.not. this%bprime%initialized)call this%bprime%init_symmetric(f = coop_ellipse_collapse_bprime_reduced, nx = 501, xmin = 1.d-6, xmax = 1.d6, xlog = .true., name = "BPRIME")
      !set  D(a)/a
      !COOP_REAL,optional::dlnD_dlnatab(na),chitab(na),ttab(na),Dtab(na),omdetab(na),Htab(na),HBK_tab(na)
      !defined by COOP D_atab(na), atab(na)
      !type(coop_function)::chi_tab,omde_tab,t_tab,H_tab,HBK_w_tab

      if(cosmology_updated)then
        write(*,*)'============================'
        write(*,*)'update cosmology: '
        write(*,*)'Omega_m: ',this%Omega_m
        write(*,*)'Omega_k: ',this%Omega_k
        write(*,*)'w: ',this%w
        write(*,*)'wa: ',this%wa
        write(*,*)'epsilon_s: ',this%epsilon_s
        write(*,*)'h: ',this%h
        write(*,*)'Omega_de: ',this%Omega_de
        write(*,*)'Omega_r: ',this%Omega_r
        write(*,*)'============================'
      end if

      select type(this)
      type is(coop_ellipse_collapse_params)
          if( cosmology_updated .or. .not.this%Dbya%initialized)then
              !add by furen
              call this%chi_tab%free()
              call this%omde_tab%free()
              call this%t_tab%free()
              call this%H_tab%free()
              call this%HBK_w_tab%free()
              !end
              call this%Dbya%free()
              call this%Dbyadot%free()
              !d(D/a)/dt
              call coop_set_uniform_pkp(na, a, amin_D, amax_D, logscale=.true.)
!              call coop_set_uniform(na, a, amin_D, amax_D,logscale=.true.)
              !              call this%set_Growth_initial_conditions(a(1), y)
              call this%init_tab(a(1),y,y_chit)
              !dadt = this%dadt(a_ini)
              !H = dadt/a_ini
              Dbya(1) = y(1)
              Dbyadot(1) = y(2)
              !add by furen
              chitab(1)=y_chit(1)
              ttab(1)=y_chit(2)

              adynamic = a(1)
              adynamic_t = a(1)
              ind = 1
              ind_t = 1
              aeq = (this%Omega_m/this%omega_de)**(1.d0/(3.d0-1.08*(1.d0-this%omega_m)*this%epsilon_s))  !!see Huang, Bond, Kofman, 2011 ApJ
              omdetab(1)=fit_HBK_rho_DE_ratio(this%epsilon_s, aeq, a(1))
              HBK_tab(1)= HBK_w(this%epsilon_s, aeq, a(1))
              Htab(1) = (this%dadt(a(i)))/a(1)
              do i  = 2, na
                  !  subroutine coop_dverk_with_ellipse_collapse_params(n, fcn, params, x, y, xend, tol, ind, c, nw, w)
!                  write(*,*)'dverk for D'
                  call coop_dverk_with_ellipse_collapse_params(2, coop_ellipse_collapse_growth_ode, this, adynamic, y, a(i), tol, ind, c, 2, wspace)
                  Dbya(i) = y(1)
                  Dbyadot(i) = y(2)
                  !dadt = this%dadt(a_ini)
                  !add by furen deng
                  omdetab(i)=fit_HBK_rho_DE_ratio(this%epsilon_s, aeq, a(i))
                  HBK_tab(i)= HBK_w(this%epsilon_s, aeq, a(i))
                  Htab(i) = (this%dadt(a(i)))/a(i)
!                  write(*,*)'dverk for chi'
                  call coop_dverk_with_ellipse_collapse_params(2, ode_for_chi_t, this, adynamic_t, y_chit, a(i), tol, ind_t, c_t, 2, wspace_t)
                  chitab(i)=y_chit(1)
                  ttab(i)=y_chit(2)
                  if((.not.now) .and. (a(i).ge.1.0))then
                      i0=i
                      now=.true.
                  end if
                  !end
              enddo
              if(.not.now)stop 'amax < 1.0, cannot initialize universe!'

              ainterp=(a(i0)-1.0)/(a(i0)-a(i0-1))
!              write(*,*)ainterp
              Dtab=dbya*a
              Dnow=(Dtab(i0)*(1.0-ainterp)+Dtab(i0-1)*ainterp)

              Dtab=Dtab/Dnow
              dbyadot = dbyadot/Dnow
              dbya = dbya/Dnow
              Dnow=Dnow/Dnow

              Htab=Htab*this%h
              dlnD_dlnatab=1.0+a*dbyadot/Dtab/Htab*this%h

              tnow=ttab(i0)*(1.0-ainterp)+ttab(i0-1)*ainterp
              taunow=chitab(i0)*(1.0-ainterp)+chitab(i0-1)*ainterp
              anow=a(i0)*(1.0-ainterp)+a(i0-1)*ainterp
              HD_Ha_now=dlnD_dlnatab(i0)*(1.0-ainterp)+dlnD_dlnatab(i0-1)*ainterp

              chitab=3000.0*(taunow-chitab)
              taunow=3000.0*taunow
              !end
              !type(coop_function)::chi_tab,omde_tab,t_tab,H_tab,HBK_w_tab

              call this%Dbya%init(n = na, xmin = amin_D, xmax = amax_D, f = Dbya, method=COOP_INTERPOLATE_SPLINE, check_boundary = .false., name="D(a)/a")
              call this%Dbyadot%init(n = na, xmin = amin_D, xmax = amax_D, f = Dbyadot, method=COOP_INTERPOLATE_SPLINE, fleft = 0.d0, check_boundary = .false., name="dot(D(a)/a)")
              call this%chi_tab%init(n = na, xmin = amin_D, xmax = amax_D, f = chitab, method=COOP_INTERPOLATE_SPLINE, check_boundary = .false., name="chi(a)")
              call this%omde_tab%init(n = na, xmin = amin_D, xmax = amax_D, f = omdetab, method=COOP_INTERPOLATE_SPLINE,check_boundary = .false., name="omde(a)")
              call this%t_tab%init(n = na, xmin = amin_D, xmax = amax_D, f = ttab, method=COOP_INTERPOLATE_SPLINE,check_boundary = .false., name="t(a)")
              call this%H_tab%init(n = na, xmin = amin_D, xmax = amax_D, f = Htab, method=COOP_INTERPOLATE_SPLINE,check_boundary = .false., name="H(a)")
              call this%HBK_w_tab%init(n = na, xmin = amin_D, xmax = amax_D, f = HBK_tab, method=COOP_INTERPOLATE_SPLINE,check_boundary = .false., name="HBKw(a)")

              !write table by furen
              !for txt output, use read(,*)is ok
              !for binary output, the variable should be real(8) for COOP_REAL or integer(4) for COOP_INT

              if(full_evolution)then
                  open(unit=4,file=outfile,access='stream',status='replace')
                  write(4)nmax+1,anow,tnow,taunow,HD_Ha_now,Dnow,a,ttab,chitab,dlnD_dlnatab,Dtab,omdetab,HBK_tab,Htab
                  close(4)
              else
                  do i = 1, na
                      if(a(i).ge.astart)then
                          istart=i
                          exit
                      end if
                  end do
                  open(unit=4,file=outfile,access='stream',status='replace')
                  write(4)nmax+2-istart,anow,tnow,taunow,HD_Ha_now,Dnow,&
                          a(istart:),ttab(istart:),chitab(istart:),dlnD_dlnatab(istart:),Dtab(istart:),Htab(istart:),&
                          omdetab(istart:),HBK_tab(istart:)
                  close(4)
              end if

              !test, delete for use
              open(unit=4,file='./test_pkp.txt',status='replace')
              write(4,*)"parameters:"
              write(4,*)"     Omega_m      w      wa      epsilon_s      h      Omega_k      amin      amax      nmax"
              write(4,*)Omega_m, w, wa, epsilon_s, h, Omega_k, amin, amax, nmax
              write(4,*)"how to read binary table:"
              write(4,*)"COOP_INT=>integer(4)=>np.int16, COOP_REAL=>real(8)=>np.float64"
              write(4,*)"========================================================"
              write(4,*)"for linear table:"
              write(4,*)"size_of_table,anow,tnow,taunow,HD_Ha_now,Dnow,a,ttab,chitab,dlnD_dlnatab,Dtab,Htab,omde_omde0tab,HBK_tab"
              write(4,*)"int(4), real(8)*5, real(8,len=size_of_tab)*8"
              write(4,*)"========================================================"
              write(4,*)"for nonlinear table:"
              write(4,*)"nf, ne, np, fmin, fmax, emin, emax, p_emin, p_emax, zvir1"
              write(4,*)"int(4)*3,real(8)*6,real(8,len=product_of_first_three_int)*1"
              write(4,*)"========================================================"
              write(4,*)"size of table:",na
              write(4,*)"========================================================"
              write(4,*)'      anow      ','      tnow      ','      taunow      ','      HD_Ha_now      ','      Dnow      '
              write(4,*)anow,tnow,taunow,HD_Ha_now,Dnow
              write(4,*)"========================================================"
              write(4,*)"the full evolution table:"
              write(4,*)'      a      ','      ttab      ','      chitab      ','      dlnD_dlnatab      ','      Dtab      ','      Htab      ','      omde_omde0tab      ','      HBK_tab      '
              do i=1, na
                  write(4,*)a(i),ttab(i),chitab(i),dlnD_dlnatab(i),Dtab(i),Htab(i),omdetab(i),HBK_tab(i)
              end do
              close(4)
              !end

          endif
      class default
          stop "extended class for ode is not yet supported"
      end select
  end subroutine coop_make_table_for_pkp

  subroutine pkp_init_COOP(this, Omega_m, w, wa, epsilon_s, h, Omega_k, lambda, F_pk, e_nu, p_nu,amin, amax, astart,nmax)!element in a table is nmax+1

    !variable by furen
    COOP_INT,optional::nmax
    COOP_REAL,optional::amin,amax,astart

    COOP_REAL::y_chit(2),adynamic_t
    COOP_REAL::c_t(24), wspace_t(2, 9),anow,tnow,taunow,ainterp,Dnow,aeq
    logical::now=.False.
    COOP_INT::ind_t,i0,istart
!    type(coop_function)::chi_tab,omde_tab,t_tab,H_tab,HBK_w_tab
!    dlnD_dlnatab(nmax+1),chitab(nmax+1),ttab(nmax+1),Dtab(nmax+1),omdetab(nmax+1),Htab(nmax+1),HBK_tab(nmax+1),astart

    !defined by COOP a(na), Dbya(na), Dbyadot(na)
    !end

    class(coop_ellipse_collapse_params)::this
    type(coop_ellipse_collapse_params)::params_norm
    COOP_INT::na
    COOP_REAL,dimension(:),allocatable::a, Dbya, Dbyadot
    COOP_REAL,dimension(:),allocatable::chitab,omdetab,Htab,HBK_tab,ttab,Dtab
    COOP_REAL:: adynamic
    COOP_INT::i
    COOP_REAL,optional::Omega_m, w, wa, epsilon_s, Omega_k, h, lambda(3), F_pk, e_nu, p_nu
    logical cosmology_updated
    !!!!for D(a) solver
    COOP_REAL::y(2)
    COOP_REAL,parameter::tol = 1.d-8
    COOP_INT::ind
    COOP_REAL::c(24), wspace(2, 9)
    COOP_REAL::amin_D
    COOP_REAL::amax_D

    cosmology_updated = .false.
    !add by furen deng
    if(present(nmax))then
      this%nmax=nmax
    end if
    if(present(amin))then
      this%amin=amin
    end if
    if(present(amax))then
      this%amax=amax
    end if

    na = this%nmax+1
    amin_D = this%amin
    amax_D = this%amax

    allocate(a(na), Dbya(na), Dbyadot(na), chitab(na), omdetab(na), Htab(na), HBK_tab(na), ttab(na), Dtab(na))
    !end

    if(present(Omega_m))then
      if(abs(this%Omega_m - Omega_m) .gt. 1.d-5)then
        this%Omega_m = Omega_m
        cosmology_updated = .true.
        write(*,*)'update Omega_m: ',Omega_m
      endif
    endif
    if(present(Omega_k))then
      if(abs(this%Omega_k - Omega_k) .gt. 1.d-5)then
        this%Omega_k = Omega_k
        cosmology_updated = .true.
        write(*,*)'update Omega_k: ',Omega_k
      endif
    endif

    if(present(w))then
      if(abs(this%w - w) .gt. 1.d-5)then
        this%w = w
        cosmology_updated = .true.
        write(*,*)'update w: ',w
      endif
    endif
    if(present(wa))then
      if(abs(this%wa - wa) .gt. 1.d-5)then
        this%wa = wa
        cosmology_updated = .true.
        write(*,*)'update wa: ',wa
      endif
    endif
    if(present(epsilon_s))then
      if(abs(this%epsilon_s - epsilon_s) .gt. 1.d-5)then
        this%epsilon_s = epsilon_s
        cosmology_updated = .true.
        write(*,*)'update epsilon_s: ',epsilon_s
      endif
    endif
    if(abs(this%epsilon_s) .gt. 1.d-5 .and. (abs(this%w+1.d0).gt.1.d-5 .or. abs(this%wa) .gt. 1.d-5))then
      write(*,*) "Warning: simultaneously setting w, wa, and epsilon_s. Use epsilon_s by default."
    endif

    if(present(h))then
      if(abs(this%h - h) .gt. 1.d-5)then
        this%h = h
        cosmology_updated = .true.
        write(*,*)'update h: ',h
      endif
    endif

#ifdef normalize_zvir1

    if(cosmology_updated.and.(abs(this%epsilon_s) .gt. 1.d-5 .or. abs(this%w+1.d0).gt.1.d-5 .or. abs(this%wa) .gt. 1.d-5))then

      call params_norm%init(Omega_m = this%Omega_m, Omega_k = this%Omega_k, h = this%h, w = -1.d0, epsilon_s = 0.d0, wa=0.d0)

      this%D_norm_zvir1 = params_norm%Growth_D(1.d0/(this%z_norm_zvir1 + 1.d0))
      write(*,*)'normalize to z = ',this%z_norm_zvir1
    end if
#endif

    this%Omega_r = 4.187d-5/this%h**2*(this%T_CMB/2.726)**4
    this%Omega_de = 1.d0 - this%Omega_m - this%Omega_r - this%Omega_k
    if(present(lambda))then
      this%lambda = lambda
      if(present(F_pk) .or. present(p_nu) .or. present(e_nu))then
        stop "You pass either lambda or (F_pk, e_nu, p_nu) to init, not both."
      endif
    elseif(present(F_pk) .and. present(p_nu) .and. present(e_nu))then
      if(F_pk .lt. 0.d0 .or. e_nu .lt. 0.d0 .or. abs(p_nu) .gt. e_nu)then
        write(*,*) "Invalid F_pk, e_nu, p_nu:"
        write(*,*) F_pk, e_nu, p_nu
        write(*,*) "check conditions: F_pk >=0, e_nu >= 0, and |p_nu| <= e_nu"
        stop
      endif
      this%lambda(3) = (F_pk/3.d0)*(1.d0 + 3.d0*e_nu + p_nu)
      this%lambda(2) = (F_pk/3.d0)*(1.d0 - 2.d0*p_nu)
      this%lambda(1) = (F_pk/3.d0)*(1.d0 - 3.d0*e_nu + p_nu)
    endif
    this%is_spherical = abs(this%lambda(1) - this%lambda(2)) .lt. 1.d-8 .and. abs(this%lambda(1) - this%lambda(3)) .lt. 1.d-8 .and. abs(this%collapse_a_ratio(1) - this%collapse_a_ratio(2)) .lt. 1.d-3 .and. abs(this%collapse_a_ratio(1)-this%collapse_a_ratio(3)) .lt. 1.d-3
    !!set b' function
    if(.not. this%bprime%initialized)call this%bprime%init_symmetric(f = coop_ellipse_collapse_bprime_reduced, nx = 501, xmin = 1.d-6, xmax = 1.d6, xlog = .true., name = "BPRIME")
    !set  D(a)/a
    !COOP_REAL,optional::dlnD_dlnatab(na),chitab(na),ttab(na),Dtab(na),omdetab(na),Htab(na),HBK_tab(na)
    !defined by COOP D_atab(na), atab(na)
    !type(coop_function)::chi_tab,omde_tab,t_tab,H_tab,HBK_w_tab

    if(cosmology_updated)then
      write(*,*)'============================'
      write(*,*)'update cosmology: '
      write(*,*)'Omega_m: ',this%Omega_m
      write(*,*)'Omega_k: ',this%Omega_k
      write(*,*)'w: ',this%w
      write(*,*)'wa: ',this%wa
      write(*,*)'epsilon_s: ',this%epsilon_s
      write(*,*)'h: ',this%h
      write(*,*)'Omega_de: ',this%Omega_de
      write(*,*)'Omega_r: ',this%Omega_r
      write(*,*)'============================'
    end if

    select type(this)
    type is(coop_ellipse_collapse_params)
      if( cosmology_updated .or. .not.this%Dbya%initialized)then
        !add by furen
        call this%chi_tab%free()
        call this%omde_tab%free()
        call this%t_tab%free()
        call this%H_tab%free()
        call this%HBK_w_tab%free()
        !end
        call this%Dbya%free()
        call this%Dbyadot%free()
        !d(D/a)/dt
        call coop_set_uniform_pkp(na, a, amin_D, amax_D, logscale=.true.)
        !              call coop_set_uniform(na, a, amin_D, amax_D,logscale=.true.)
        !              call this%set_Growth_initial_conditions(a(1), y)
        call this%init_tab(a(1),y,y_chit)
        !dadt = this%dadt(a_ini)
        !H = dadt/a_ini
        Dbya(1) = y(1)
        Dbyadot(1) = y(2)
        !add by furen
        chitab(1)=y_chit(1)
        ttab(1)=y_chit(2)

        adynamic = a(1)
        adynamic_t = a(1)
        ind = 1
        ind_t = 1
        aeq = (this%Omega_m/this%omega_de)**(1.d0/(3.d0-1.08*(1.d0-this%omega_m)*this%epsilon_s))  !!see Huang, Bond, Kofman, 2011 ApJ
        omdetab(1)=fit_HBK_rho_DE_ratio(this%epsilon_s, aeq, a(1))
        HBK_tab(1)= HBK_w(this%epsilon_s, aeq, a(1))
        Htab(1) = (this%dadt(a(i)))/a(1)
        do i  = 2, na
          !  subroutine coop_dverk_with_ellipse_collapse_params(n, fcn, params, x, y, xend, tol, ind, c, nw, w)
          call coop_dverk_with_ellipse_collapse_params(2, coop_ellipse_collapse_growth_ode, this, adynamic, y, a(i), tol, ind, c, 2, wspace)
          Dbya(i) = y(1)
          Dbyadot(i) = y(2)
          !dadt = this%dadt(a_ini)
          !add by furen deng
          omdetab(i)=fit_HBK_rho_DE_ratio(this%epsilon_s, aeq, a(i))
          HBK_tab(i)= HBK_w(this%epsilon_s, aeq, a(i))
          Htab(i) = (this%dadt(a(i)))/a(i)
          call coop_dverk_with_ellipse_collapse_params(2, ode_for_chi_t, this, adynamic_t, y_chit, a(i), tol, ind_t, c_t, 2, wspace_t)
          chitab(i)=y_chit(1)
          ttab(i)=y_chit(2)
          if((.not.now) .and. (a(i).ge.1.0))then
            i0=i
            now=.true.
          end if
          !end
        enddo
        if(.not.now)stop 'amax < 1.0, cannot initialize universe!'

        ainterp=(a(i0)-1.0)/(a(i0)-a(i0-1))
!        write(*,*)ainterp
        Dtab=dbya*a
        Dnow=(Dtab(i0)*(1.0-ainterp)+Dtab(i0-1)*ainterp)

        Dtab=Dtab/Dnow
        dbyadot = dbyadot/Dnow
        dbya = dbya/Dnow
        Dnow=Dnow/Dnow

        Htab=Htab*this%h
!        dlnD_dlnatab=1.0+a*dbyadot/Dtab/Htab*this%h

        tnow=ttab(i0)*(1.0-ainterp)+ttab(i0-1)*ainterp
        taunow=chitab(i0)*(1.0-ainterp)+chitab(i0-1)*ainterp
        anow=a(i0)*(1.0-ainterp)+a(i0-1)*ainterp
!        HD_Ha_now=dlnD_dlnatab(i0)*(1.0-ainterp)+dlnD_dlnatab(i0-1)*ainterp

        chitab=3000.0*(taunow-chitab)
        taunow=3000.0*taunow
        !end
        !type(coop_function)::chi_tab,omde_tab,t_tab,H_tab,HBK_w_tab

        call this%Dbya%init(n = na, xmin = amin_D, xmax = amax_D, f = Dbya, method=COOP_INTERPOLATE_SPLINE, check_boundary = .false., name="D(a)/a")
        call this%Dbyadot%init(n = na, xmin = amin_D, xmax = amax_D, f = Dbyadot, method=COOP_INTERPOLATE_SPLINE, fleft = 0.d0, check_boundary = .false., name="dot(D(a)/a)")
        call this%chi_tab%init(n = na, xmin = amin_D, xmax = amax_D, f = chitab, method=COOP_INTERPOLATE_SPLINE, check_boundary = .false., name="chi(a)")
        call this%omde_tab%init(n = na, xmin = amin_D, xmax = amax_D, f = omdetab, method=COOP_INTERPOLATE_SPLINE,check_boundary = .false., name="omde(a)")
        call this%t_tab%init(n = na, xmin = amin_D, xmax = amax_D, f = ttab, method=COOP_INTERPOLATE_SPLINE,check_boundary = .false., name="t(a)")
        call this%H_tab%init(n = na, xmin = amin_D, xmax = amax_D, f = Htab, method=COOP_INTERPOLATE_SPLINE,check_boundary = .false., name="H(a)")
        call this%HBK_w_tab%init(n = na, xmin = amin_D, xmax = amax_D, f = HBK_tab, method=COOP_INTERPOLATE_SPLINE,check_boundary = .false., name="HBKw(a)")

        !write table by furen
        !for txt output, use read(,*)is ok
        !for binary output, the variable should be real(8) for COOP_REAL or integer(4) for COOP_INT

      endif
    class default
      stop "extended class for ode is not yet supported"
    end select

    deallocate(a, Dbya, Dbyadot,chitab,omdetab,Htab,HBK_tab,ttab,Dtab)

  end subroutine pkp_init_COOP


  subroutine ode_for_chi_t(n, a, y, dyda, params)
      type(coop_ellipse_collapse_params)::params
      COOP_INT::n
      COOP_REAL::a, dadt
      COOP_REAL::y(2), dyda(2)
      dadt=params%dadt(a)
      dyda(1)=1/(dadt)/a/params%h
      dyda(2)=1/(dadt)/params%h
  end subroutine ode_for_chi_t

!  subroutine ode_forw_chi_t(this,n, a, y, dyda, params)
!      class(coop_ellipse_collapse_params)::this
!      type(coop_ellipse_collapse_params)::params
!      COOP_INT::n
!      COOP_REAL::a, dadt
!      COOP_REAL::y(2), dyda(2)
!      dadt=params%dadt(a)
!      dyda(1)=1/(dadt)/a
!      dyda(2)=1/(dadt)
!  end subroutine ode_forw_chi_t

  subroutine set_initial_condition_for_table(this, a_ini, y, y_chit)
      class(coop_ellipse_collapse_params)::this
      COOP_REAL::a_ini, dadt, addbya, H
      COOP_REAL::y(2),y_chit(2)
      COOP_REAL::a_eq, t0, t1
      !aeq=(Omer+Omerd)/(Omnr+Omnrd)
      !t1=1./h/sqrt(Omnr)
      !t0=1./h/sqrt(Omnr+Omnrd)
      !y=sqrt(ai+aeq)
      !y0=sqrt(aeq)
      !t=t0*(0.66666667*(y**3-y0**3)-2.*aeq*(y-y0))
      !tau=t0*2.*(y-y0)
      !tin_Dlin=t
      !tauin_Dlin=tau
      dadt = this%dadt(a_ini)
      H = dadt/a_ini
      addbya = this%ddotabya(a_ini)

      a_eq = sqrt(this%Omega_r/this%Omega_m)
      t0 = sqrt(a_ini + a_eq)
      t1 = sqrt(a_eq)

      y(1) = 1.d0
      !when turning off radiation, y(2) -> 0 when a_ini -> 0
      y(2) =  - (2.d0*H**2+addbya-1.5d0*this%Omega_m/a_ini**3)*y(1)/4.d0/H

      !add influence of radiation

      y_chit(1) = t0 * 2.d0/sqrt(this%Omega_m)/this%h
      y_chit(2) = ((t0**3 - t1**3) * 2.d0/3.d0 - 2 * a_eq * (t0 - t1))/sqrt(this%Omega_m)/this%h
!      write(*,*)y_chit, t0, t1
!      y_chit(1)=sqrt(a_ini/this%Omega_m)*2.0/this%h
!      y_chit(2)=sqrt(a_ini**3/this%Omega_m)*2.0/3.0/this%h
  end subroutine set_initial_condition_for_table
  !end==================================================================================================================

  subroutine coop_ellipse_collapse_params_init(this, Omega_m, w, wa, epsilon_s, h, Omega_k, lambda, F_pk, e_nu, p_nu)
    class(coop_ellipse_collapse_params)::this
    COOP_INT,parameter::na = 256
    COOP_REAL::a(na), Dbya(na), Dbyadot(na), adynamic
    COOP_INT::i
    COOP_REAL,optional::Omega_m, w, wa, epsilon_s, Omega_k, h, lambda(3), F_pk, e_nu, p_nu
    logical cosmology_updated
    !!!!for D(a) solver
    COOP_REAL::y(2)
    COOP_REAL,parameter::tol = 1.d-8
    COOP_INT::ind
    COOP_REAL::c(24), wspace(2, 9)
    COOP_REAL,parameter::amin_D = 0.03d0

    cosmology_updated = .false.
    if(present(Omega_m))then
       if(abs(this%Omega_m - Omega_m) .gt. 1.d-5)then
          this%Omega_m = Omega_m
          cosmology_updated = .true.
       endif
    endif
    if(present(Omega_k))then
       if(abs(this%Omega_k - Omega_k) .gt. 1.d-5)then
          this%Omega_k = Omega_k
          cosmology_updated = .true.
       endif
    endif
    if(present(w))then
       if(abs(this%w - w) .gt. 1.d-5)then
          this%w = w
          cosmology_updated = .true.
       endif
    endif
    if(present(wa))then
       if(abs(this%wa - wa) .gt. 1.d-5)then
          this%wa = wa
          cosmology_updated = .true.
       endif
    endif
    if(present(epsilon_s))then
       if(abs(this%epsilon_s - epsilon_s) .gt. 1.d-5)then
          this%epsilon_s = epsilon_s
          cosmology_updated = .true.
       endif
    endif
    if(abs(this%epsilon_s) .gt. 1.d-5 .and. (abs(this%w+1.d0).gt.1.d-5 .or. abs(this%wa) .gt. 1.d-5))then
       write(*,*) "Warning: simultaneously setting w, wa, and epsilon_s. Use epsilon_s by default."
    endif
    if(present(h))then
       if(abs(this%h - h) .gt. 1.d-5)then
          this%h = h
          cosmology_updated = .true.
       endif
    endif
    this%Omega_r = 4.187d-5/this%h**2*(this%T_CMB/2.726)**4
    this%Omega_de = 1.d0 - this%Omega_m - this%Omega_r - this%Omega_k
    if(present(lambda))then
       this%lambda = lambda
       if(present(F_pk) .or. present(p_nu) .or. present(e_nu))then
          stop "You pass either lambda or (F_pk, e_nu, p_nu) to init, not both."
       endif
    elseif(present(F_pk) .and. present(p_nu) .and. present(e_nu))then
       if(F_pk .lt. 0.d0 .or. e_nu .lt. 0.d0 .or. abs(p_nu) .gt. e_nu)then
          write(*,*) "Invalid F_pk, e_nu, p_nu:"
          write(*,*) F_pk, e_nu, p_nu
          write(*,*) "check conditions: F_pk >=0, e_nu >= 0, and |p_nu| <= e_nu"
          stop
       endif
       this%lambda(3) = (F_pk/3.d0)*(1.d0 + 3.d0*e_nu + p_nu)
       this%lambda(2) = (F_pk/3.d0)*(1.d0 - 2.d0*p_nu)
       this%lambda(1) = (F_pk/3.d0)*(1.d0 - 3.d0*e_nu + p_nu)
    endif
    this%is_spherical = abs(this%lambda(1) - this%lambda(2)) .lt. 1.d-8 .and. abs(this%lambda(1) - this%lambda(3)) .lt. 1.d-8 .and. abs(this%collapse_a_ratio(1) - this%collapse_a_ratio(2)) .lt. 1.d-3 .and. abs(this%collapse_a_ratio(1)-this%collapse_a_ratio(3)) .lt. 1.d-3
    !!set b' function
    if(.not. this%bprime%initialized)call this%bprime%init_symmetric(f = coop_ellipse_collapse_bprime_reduced, nx = 501, xmin = 1.d-6, xmax = 1.d6, xlog = .true., name = "BPRIME")
    !!set  D(a)/a
    select type(this)
    type is(coop_ellipse_collapse_params)

       if( cosmology_updated .or. .not.this%Dbya%initialized)then
          call this%Dbya%free()
          call this%Dbyadot%free()
          call coop_set_uniform_pkp(na, a, amin_D, 1.d0)
!          call coop_set_uniform(na, a, amin_D, 1.d0)
          call this%set_Growth_initial_conditions(a(1), y)
          Dbya(1) = y(1)
          Dbyadot(1) = y(2)
          adynamic = a(1)
          ind = 1
          do i  = 2, na
             call coop_dverk_with_ellipse_collapse_params(2, coop_ellipse_collapse_growth_ode, this, adynamic, y, a(i), tol, ind, c, 2, wspace)
             Dbya(i) = y(1)
             Dbyadot(i) = y(2)
          enddo
          dbyadot = dbyadot/dbya(na)
          dbya = dbya/dbya(na)
          call this%Dbya%init(n = na, xmin = amin_D, xmax = 1.d0, f = Dbya, method=COOP_INTERPOLATE_SPLINE, check_boundary = .false., name="D(a)/a")
          call this%Dbyadot%init(n = na, xmin = amin_D, xmax = 1.d0, f = Dbyadot, method=COOP_INTERPOLATE_SPLINE, fleft = 0.d0, check_boundary = .false., name="dot(D(a)/a)")
       endif
    class default
       stop "extended class for ode is not yet supported"
    end select
  end subroutine coop_ellipse_collapse_params_init


  subroutine coop_ellipse_collapse_params_get_bprime(this, x, bprime)
    class(coop_ellipse_collapse_params)::this
    COOP_REAL,intent(IN)::x(3)
    COOP_REAL,intent(OUT)::bprime(3)
    bprime(1) = this%bprime%eval( (x(1)/x(2))**2,  (x(1)/x(3))**2 )
    bprime(2) = this%bprime%eval( (x(2)/x(1))**2,  (x(2)/x(3))**2 )
    bprime(3) = this%bprime%eval( (x(3)/x(2))**2,  (x(3)/x(1))**2 )
  end subroutine coop_ellipse_collapse_params_get_bprime

  !!D(a)
  function coop_ellipse_collapse_params_Growth_D(this, a) result(D)
    class(coop_ellipse_collapse_params)::this
    COOP_REAL::a, D
    D = this%Dbya%eval(a)*a
  end function coop_ellipse_collapse_params_Growth_D

  !!d ln D/dt 
  function coop_ellipse_collapse_params_Growth_H_D(this, a) result(H_D)
    class(coop_ellipse_collapse_params)::this
    COOP_REAL::a, H_D
    H_D = this%Dbyadot%eval(a)/this%Dbya%eval(a) + this%dadt(a)/a
  end function coop_ellipse_collapse_params_Growth_H_D
  
  !! H a / (H_0 a_0)
  function coop_ellipse_collapse_params_aH(this, a) result(aH)
    class(coop_ellipse_collapse_params)::this
    COOP_REAL::a, aH, aeq
    if(this%epsilon_s .ne. 0.d0)then
       aeq = (this%Omega_m/this%omega_de)**(1.d0/(3.d0-1.08*(1.d0-this%omega_m)*this%epsilon_s))  !!see Huang, Bond, Kofman, 2011 ApJ
       if(aeq .lt. 0.1 .or. aeq .gt. 0.99) stop "Bad arguments for HBK parametrization. Are you using a very large epsilon_s? "
       aH = sqrt(this%Omega_k + (this%Omega_m+this%Omega_r/a)/a + this%omega_de*fit_HBK_rho_DE_ratio(this%epsilon_s, aeq, a)*a**2)
    else
       aH = sqrt(this%Omega_k + (this%Omega_m+this%Omega_r/a)/a + this%omega_de*a**(-1.d0-3.d0*(this%w+this%wa))*exp(-3.d0*this%wa*(1.d0-a)))
    endif
  end function coop_ellipse_collapse_params_aH

!!return H_0^2  \ddot a / a
  function coop_ellipse_collapse_params_ddotabya(this, a) result(ddotabya)
    class(coop_ellipse_collapse_params)::this
    COOP_REAL::a, ddotabya, aeq
    if(this%epsilon_s .ne. 0.d0)then
       aeq = (this%Omega_m/this%omega_de)**(1.d0/(3.d0-1.08*(1.d0-this%omega_m)*this%epsilon_s))  !!see Huang, Bond, Kofman, 2011 ApJ
       ddotabya = -( (this%Omega_r/a*2.d0 + this%Omega_m)/a**3 + this%Omega_de*fit_HBK_rho_DE_ratio(this%epsilon_s, aeq, a)*(1.d0+3.d0*HBK_w(this%epsilon_s, aeq, a)) )/2.d0
    else
       ddotabya = -( (this%Omega_r/a*2.d0 + this%Omega_m)/a**3 + this%Omega_de*a**(-3.d0*(1.d0+this%w+this%wa))*exp(-3.d0*this%wa*(1.d0-a))*(1.d0+3.d0*(this%w+this%wa*(1.d0-a))) )/2.d0
    endif
  end function coop_ellipse_collapse_params_ddotabya

  subroutine coop_ellipse_collapse_params_free(this)
    class(coop_ellipse_collapse_params)::this
    call this%bprime%free()
    call this%Dbya%free()
    call this%Dbyadot%free()
  end subroutine coop_ellipse_collapse_params_free


!!=========================================================================
!!utilities; no need to understand or change anything below
  subroutine coop_dverk_with_ellipse_collapse_params(n, fcn, params, x, y, xend, tol, ind, c, nw, w)
    type(coop_ellipse_collapse_params) params
#define DVERK_ARGUMENTS ,params
#include "dverk.h"    
#undef DVERK_ARGUMENTS
  end subroutine coop_dverk_with_ellipse_collapse_params


  function coop_ellipse_collapse_bprime_reduced(lambda1, lambda2) result(bprime)
    COOP_REAL::lambda1, lambda2, bprime
    if(abs(lambda1-1.d0) .lt. 1.d-4 .and. abs(lambda2-1.d0) .lt. 1.d-4)then !!both close to 0
       bprime = (2.d0/15.d0)*(2.d0-lambda1 - lambda2) + (3.d0/28.d0-0.05d0)*((1.d0-lambda1)**2+(1.d0-lambda2)**2+(2.d0/3.d0)*(1.d0-lambda1)*(1.d0-lambda2))
       return
    endif
    bprime = (coop_elliptic_Rd(1.d0/lambda1, 1.d0/lambda2, 1.d0)/sqrt(lambda1*lambda2)-1.d0)*2.d0/3.d0
  end function coop_ellipse_collapse_bprime_reduced

  !!return rhoDE(a) / rhoDE_0
  function fit_HBK_rho_DE_ratio(epsilon_s, aeq, a) result(rat)
    COOP_REAL::a, aeq, epsilon_s, rat
    rat = exp(2.d0*epsilon_s*(fit(1.d0/aeq) - fit(a/aeq)))
  contains
    function fit(x)
      COOP_REAL::x, fit
      COOP_REAL,parameter::x1 = 1.49d0, x2 = 1.51d0
      if(x.lt. x1)then
         fit = x**3*(4.d0/27.d0)/(1.d0+0.3d0*x**3)*(1.d0+x**4.95*0.02)
      elseif(x .gt. x2)then
         fit = log(x+0.6d0/x**0.9) - 0.37123
      else
         fit = ((x1**3*(4.d0/27.d0)/(1.d0+0.3d0*x1**3)*(1.d0+x1**4.95*0.02))*(x2-x) + (log(x2+0.6d0/x2**0.9) - 0.37123)*(x-x1))/(x2-x1)
      endif
    end function fit
  end function fit_HBK_rho_DE_ratio

  function HBK_w(epsilon_s, aeq, a) result(w)
    COOP_REAL::epsilon_s, aeq, a, w
    w = -1.d0+ 2.d0/3.d0*epsilon_s* F(a/aeq)**2
  contains
    function F(x)
      COOP_REAL::x, F
      if(x.gt. 0.03d0)then
         F = sqrt(1.d0+ x**(-3)) - log(x**1.5d0 + sqrt(1.d0+x**3))/x**3
      else
         F = x**1.5*(2.d0/3.d0 - x**3/5.d0)
      endif
    end function F
  end function HBK_w

  !import by furen
  !mathematical methods
  FUNCTION coop_elliptic_Rd(x,y,z) result(rd)
    COOP_REAL, INTENT(IN) :: x,y,z
    COOP_REAL :: rd
    COOP_REAL, PARAMETER :: ERRTOL=0.05d0,TINY=1.d-25,BIG=4.5d21,&
            C1=3.d0/14.d0,C2=1.d0/6.d0,C3=9.d0/22.d0,&
            C4=3.d0/26.d0,C5=0.25d0*C3,C6=1.5d0*C4
    COOP_REAL :: alamb,ave,delx,dely,delz,ea,eb,ec,ed,&
            ee,fac,sqrtx,sqrty,sqrtz,sum,xt,yt,zt
    if(.not. (min(x,y) >= 0.0 .and. min(x+y,z) >= TINY .and. max(x,y,z) <= BIG))then
      write(*,*) x, y, z
      stop "elliptic_Rd parameters out of range"
    endif
    xt=x
    yt=y
    zt=z
    sum=0.0
    fac=1.0
    do
      sqrtx=sqrt(xt)
      sqrty=sqrt(yt)
      sqrtz=sqrt(zt)
      alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
      sum=sum+fac/(sqrtz*(zt+alamb))
      fac=0.25d0*fac
      xt=0.25d0*(xt+alamb)
      yt=0.25d0*(yt+alamb)
      zt=0.25d0*(zt+alamb)
      ave=0.2d0*(xt+yt+3.0d0*zt)
      delx=(ave-xt)/ave
      dely=(ave-yt)/ave
      delz=(ave-zt)/ave
      if (max(abs(delx),abs(dely),abs(delz)) <= ERRTOL) exit
    end do
    ea=delx*dely
    eb=delz*delz
    ec=ea-eb
    ed=ea-6.d0*eb
    ee=ed+ec+ec
    rd=3.d0*sum+fac*(1.0d0+ed*(-C1+C5*ed-C6*delz*ee)&
            +delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave))
  END FUNCTION coop_elliptic_Rd
!
  subroutine coop_set_uniform_pkp(n, x, lower, upper, logscale)
    COOP_INT n, i
    COOP_REAL x(n), lower, upper, rlow, dx
    logical,optional::logscale
    if(n.eq. 0)return
    x(1) = lower
    if(n.eq.1)return
    x(n) = upper
    if(present(logscale))then
      if(logscale)then
        dx = (log(upper) - log(lower))/(n-1)
        rlow = log(lower) - dx
        !$omp parallel do
        do i = 2, n-1
          x(i) = exp(rlow + dx*i)
        enddo
        !$omp end parallel do
        return
      endif
    endif
    dx = (upper-lower)/(n-1)
    rlow = lower-dx
    !$omp parallel do
    do i = 2, n-1
      x(i) =rlow + dx*i
    enddo
  end subroutine coop_set_uniform_pkp
!
  !methods for 2dfunction

  subroutine coop_2dfunction_init_symmetric(this, f, xmin, xmax, nx, xlog, zlog, name)
    !!save f(x, y) = f(y, x) into coop_2dfunction object this
    class(coop_2dfunction)::this
    external f
    COOP_REAL::f
    COOP_INT::nx
    COOP_REAL::xmin, xmax
    logical,optional::xlog, zlog
    COOP_UNKNOWN_STRING,optional::name
    COOP_INT::ix, iy
    COOP_REAL::dx, dy
    call this%free()
    if(present(name))this%name = trim(adjustl(name))
    if(nx .lt. 2 ) stop "2d function init assumes nx >=2"
    this%nx = nx
    this%ny = nx
    allocate(this%f(nx, nx), this%fxx(nx, nx), this%fyy(nx, nx))
    if(present(xlog))then
      this%xlog = xlog
    else
      this%xlog = .false.
    endif
    this%ylog = this%xlog
    if(present(zlog))then
      this%zlog = zlog
    else
      this%zlog = .false.
    endif
    if(this%xlog)then
      if(xmin .le. 0.d0 .or. xmax .le. 0.d0) stop "2dfunction_init_symmetric: xlog = T but xmin<=0 or xmax<=0"
      this%xmin = log(xmin)
      this%xmax = log(xmax)
    else
      this%xmin = xmin
      this%xmax = xmax
    endif

    this%dx = (this%xmax - this%xmin)/(this%nx-1)
    this%xmin = this%xmin + this%dx/1.d10
    this%xmax = this%xmax - this%dx/1.d10
    this%dx = (this%xmax - this%xmin)/(this%nx-1)
    if(abs(this%dx) .le. 0.d0) stop "2dfunction_init_symmetric: xmax == xmin?"

    this%dy = this%dx
    this%ymin = this%xmin
    this%ymax = this%xmax



    if(this%xlog)then
      !$omp parallel do private(ix, iy)
      do iy = 1, this%ny
        do ix = 1, iy
          this%f(ix, iy) = f(exp(this%xmin + this%dx*(ix-1)), exp(this%ymin + this%dy*(iy-1)))
        enddo
      enddo
      !$omp end parallel do
    else
      !$omp parallel do private(ix, iy)
      do iy = 1, this%ny
        do ix = 1, iy
          this%f(ix, iy) = f(this%xmin + this%dx*(ix-1), this%ymin + this%dy*(iy-1))
        enddo
      enddo
      !$omp end parallel do
    endif
    !$omp parallel do private(ix, iy)
    do iy = 1, this%ny
      do ix = 1, iy-1
        this%f(iy, ix) = this%f(ix, iy)
      enddo
    enddo
    !$omp end parallel do

    if(this%zlog)then
      if(any(this%f .le. 0.d0)) stop "2dfunction_init: negative values"
      this%f = log(this%f)
    endif
    !$omp parallel do private(ix, iy)
    do iy = 2, this%ny-1
      do ix = 1, this%nx
        this%fyy(ix, iy) = (this%f(ix, iy+1) + this%f(ix, iy-1) - 2.d0*this%f(ix, iy))
      enddo
    enddo
    !$omp end parallel do
    !$omp parallel do private(ix, iy)
    do iy = 1, this%ny
      do ix = 2, this%nx-1
        this%fxx(ix, iy) = (this%f(ix+1, iy) + this%f(ix-1, iy) - 2.d0*this%f(ix, iy))
      enddo
    enddo
    !$omp end parallel do
    this%fxx(1, :) = this%fyy(:, 1)
    this%fxx(this%nx, :) = this%fyy(:, this%ny)
    this%fyy(:,1) = this%fxx(1,:)
    this%fyy(:,this%ny) = this%fxx(this%nx, :)

    this%fxx(1, 1) = (this%fxx(2, 1)+this%fxx(1, 2))/2.d0
    this%fxx(1, this%nx) = (this%fxx(2, this%nx) + this%fxx(1, this%nx-1))/2.d0
    this%fxx(this%nx, 1) = (this%fxx(this%nx, 2)+this%fxx(this%nx-1, 1))/2.d0
    this%fxx(this%nx, this%nx) = (this%fxx(this%nx-1, this%nx) + this%fxx(this%nx, this%nx-1))/2.d0


    this%fyy(1, 1) = (this%fyy(2, 1)+this%fyy(1, 2))/2.d0
    this%fyy(1, this%nx) = (this%fyy(2, this%nx) + this%fyy(1, this%nx-1))/2.d0
    this%fyy(this%nx, 1) = (this%fyy(this%nx, 2)+this%fyy(this%nx-1, 1))/2.d0
    this%fyy(this%nx, this%nx) = (this%fyy(this%nx-1, this%nx) + this%fyy(this%nx, this%nx-1))/2.d0

    this%fxx = this%fxx/2.d0
    this%fyy = this%fyy/2.d0
    this%initialized = .true.
    return
  end subroutine coop_2dfunction_init_symmetric

  function coop_2dfunction_evaluate_bare(this, xbare, ybare) result(zbare)
    class(coop_2dfunction)::this
    COOP_REAL::xbare, ybare, zbare
    COOP_INT::ix, iy
    COOP_REAL::rx, ry
    rx = (xbare - this%xmin)/this%dx+1.d0
    ry = (ybare - this%ymin)/this%dy+1.d0
    ix = min(max(floor(rx), 1), this%nx-1)
    rx = max(min(rx - ix, 1.d0), 0.d0)
    iy = min(max(floor(ry), 1), this%ny-1)
    ry = max(min(ry - iy, 1.d0), 0.d0)
    zbare =(this%f(ix, iy)*(1.d0-ry) + this%f(ix, iy+1)*ry)*(1.d0-rx) &
            + (this%f(ix+1, iy)*(1.d0-ry) + this%f(ix+1, iy+1)*ry)*rx &
            + ((this%fxx(ix, iy)*(1.d0-ry) + this%fxx(ix, iy+1)*ry)*(1.d0-rx) &
                    + (this%fxx(ix+1, iy)*(1.d0-ry) + this%fxx(ix+1, iy+1)*ry)*rx)*rx*(rx-1.d0) &
            + ((this%fyy(ix, iy)*(1.d0-ry) + this%fyy(ix, iy+1)*ry)*(1.d0-rx) &
                    + (this%fyy(ix+1, iy)*(1.d0-ry) + this%fyy(ix+1, iy+1)*ry)*rx)*ry*(ry-1.d0)
  end function coop_2dfunction_evaluate_bare


  function coop_2dfunction_evaluate(this, x, y) result(z)
    class(coop_2dfunction)::this
    COOP_REAL::x, y, z
    if(this%xlog)then
      if(this%ylog)then
        z = this%eval_bare(log(x), log(y))
      else
        z = this%eval_bare(log(x), y)
      endif
    else
      if(this%ylog)then
        z = this%eval_bare(x, log(y))
      else
        z = this%eval_bare(x, y)
      endif
    endif
    if(this%zlog) z = exp(z)
  end function coop_2dfunction_evaluate

  subroutine coop_2dfunction_free(this)
    class(coop_2dfunction)::this
    COOP_DEALLOC(this%f)
    COOP_DEALLOC(this%fxx)
    COOP_DEALLOC(this%fyy)
    this%xlog = .false.
    this%ylog = .false.
    this%zlog = .false.
    this%nx  = 0
    this%ny = 0
    this%initialized = .false.
    this%name = "NoName2D"
  end subroutine coop_2dfunction_free

  !methods for function

  function coop_isnan_arrd(x) result(isnan)
    COOP_REAL,dimension(:):: x
    logical isnan
    isnan = .not. (all(x.gt.0.d0 .or. x.le.0.d0))
  end function coop_isnan_arrd

  subroutine coop_spline_uniform(n, y, y2)
    COOP_INT n, i
    COOP_REAL  y(n), y2(n)
    COOP_REAL yil, yir, bet, ypl, ypr
    COOP_REAL gam(n-1)
    if(n.le.2)then
      y2 = 0.
      return
    endif
    yir=(y(2)-y(1))
    ypl=(2.d0*y(2)-1.5d0*y(1)-0.5d0*y(3))
    ypr = -(2.d0*y(n-1)-1.5d0*y(n)-0.5d0*y(n-2))
    y2(1)=(yir-ypl)*3.
    gam(1)= 0.5
    do i=2, n-1
      bet=2.d0/3.d0-gam(i-1)/6.d0
      if(abs(bet) .lt. 1.d-30) stop 'Error in spline.'
      yil=yir
      yir=(y(i+1)-y(i))
      y2(i)=(yir-yil-y2(i-1)/6.d0)/bet
      gam(i)=(1.d0/6.d0)/bet
    enddo
    bet=1.d0/3.-gam(n-1)/6.d0
    if(abs(bet) .lt. 1.d-30) stop 'Error in spline.'
    y2(n)=(ypr-yir-y2(n-1)/6.d0)/bet
    do i=n-1, 1 , -1
      y2(i)=y2(i)-gam(i)*y2(i+1)
    enddo
    y2 =  y2/6.
  end subroutine coop_spline_uniform

  subroutine coop_smooth_data_d(n, y, sigma, logscale)
    COOP_INT::n
    COOP_REAL::y(n)
    COOP_INT::sigma
    COOP_REAL::w(-3*sigma:3*sigma), ycopy(1-3*sigma:n+3*sigma)
    COOP_INT i, m
    logical, optional::logscale
    w(0) = 1.d0
    m = 3*sigma
    !$omp parallel do
    do i = 1, m
      w(i) = exp(-(dble(i)/sigma)**2/2.d0)
      w(-i) = w(i)
    enddo
    !$omp end parallel do
    w = w/sum(w)
    ycopy(1:n) = y
    ycopy(1-m:0) = y(1)
    ycopy(n+1:n+m) = y(n)
    if(present(logscale))then
      if(logscale) ycopy = log(ycopy)
    endif
    !$omp parallel do
    do i=1, n
      y(i) = sum(w*ycopy(i-m:i+m))
    enddo
    !$omp end parallel do
    if(present(logscale))then
      if(logscale) y = exp(y)
    endif
  end subroutine coop_smooth_data_d


  subroutine coop_function_init_polynomial(this, p, xlog, ylog,  name)
    class(coop_function)::this
    logical,optional::xlog, ylog
    COOP_REAL,dimension(:)::p
    COOP_UNKNOWN_STRING,optional::name
    COOP_INT::i
    call this%free()
    this%method = COOP_INTERPOLATE_POLYNOMIAL
    if(all(p.eq.0.d0))then
      this%is_zero = .true.
      this%n = 1
      allocate(this%f(this%n), this%f1(this%n), this%f2(this%n))
      this%f =  0.d0
      this%f1(this%n) = 0.d0
      this%f2(this%n) = 0.d0
      this%xmax = 1.d99
      this%xmin = - this%xmax
      this%initialized = .true.
      return
    endif
    if(present(name))this%name = trim(adjustl(name))
    if(present(xlog))then
      this%xlog = xlog
    else
      this%xlog = .false.
    endif
    if(present(ylog))then
      this%ylog = ylog
    else
      this%ylog = .false.
    endif
    this%check_boundary = .false.
    this%n = size(p)
    allocate(this%f(this%n), this%f1(this%n), this%f2(this%n))
    this%f = p
    this%f1(this%n) = 0.d0
    this%f2(this%n) = 0.d0
    if(this%n .gt.1) &
            this%f2(this%n-1) = 0.d0
    do i = 1, this%n-1
      this%f1(i) = this%f(i+1)*i
    enddo
    do i=1, this%n - 2
      this%f2(i) = this%f1(i+1)*i
    enddo
    this%xmax = 1.d99**(1.d0/max(this%n-1,1))
    this%xmin = - this%xmax
    this%initialized = .true.
  end subroutine coop_function_init_polynomial

  subroutine coop_function_free(this)
    class(coop_function):: this
    COOP_DEALLOC(this%f)
    COOP_DEALLOC(this%f1)
    COOP_DEALLOC(this%f2)
    COOP_DEALLOC(this%f3)
    this%method = COOP_INTERPOLATE_LINEAR
    this%name = "NoName"
    this%check_boundary = .true.
    this%n = 0
    this%n_down  = 0
    this%scale = 1.d0
    this%shift = 0.d0
    this%xlog = .false.
    this%ylog = .false.
    this%initialized = .false.
    this%is_zero = .false.
  end subroutine coop_function_free

  function coop_polyvalue(n, p, x) result(y)
    COOP_INT n, j
    COOP_REAL p(n)
    COOP_REAL x, y
    if(n.le.0)then
      y = 0.d0
      return
    endif
    y = p(n)
    do j= n-1,1,-1
      y= y*x + p(j)
    enddo
  end function coop_polyvalue


  function coop_function_evaluate(this, x) result(f)
    class(coop_function):: this
    COOP_REAL x, f
    if(this%xlog)then
      f = coop_function_evaluate_bare(this, log(x))
    else
      f = coop_function_evaluate_bare(this, x)
    endif
    if(this%ylog)then
      f = exp(f)
    endif
    f = f*this%scale + this%shift
  end function coop_function_evaluate

  subroutine coop_function_init(this, n, xmin, xmax, f, method, fleft, fright, slopeleft, sloperight, chebyshev_order, xlog, ylog, check_boundary, smooth, name)
    class(coop_function):: this
    logical, optional::xlog, ylog
    COOP_INT,intent(IN):: n
    COOP_REAL,intent(IN):: xmin, xmax
    COOP_REAL,intent(IN)::f(n)
    COOP_REAL, optional::fleft, fright, slopeleft, sloperight
    COOP_INT, optional::method
    COOP_INT, optional::chebyshev_order
    logical,optional::check_boundary
    logical,optional::smooth
    COOP_UNKNOWN_STRING, optional::name
    COOP_INT i, count_tiny, count_small
    COOP_REAL::fmean, ftiny, curv, flarge, fsmall
    call this%free()
    if(all(f .eq. 0.d0))then
      call this%init_polynomial( (/ 0.d0 /) )
      return
    endif
    if(present(name))then
      this%name = trim(adjustl(name))
    endif
    if(coop_isnan_arrd(f))then
      write(*,*) "Cannot construct the function "//trim(this%name)//": found f = NAN within the specified range."
      stop
    endif
    if(n.eq.1)then
      call this%init_polynomial (f)
      return
    endif
    if(present(xlog))then
      this%xlog = xlog
    else
      this%xlog = .false.
    endif
    if(present(ylog))then
      this%ylog = ylog
    else
      this%ylog = .false.
    endif
    if(this%xlog .and. (xmin .le. 0.d0 .or. xmax .le. 0.d0))then
      write(*,*)"function:init"
      write(*,*)"cannot do xlog for x<=0"
      write(*,*)"stop"
    end if
    if(this%ylog)then
      if(any(f.le.0.d0))then
        write(*,*)"function:init"
        write(*,*)"cannot do ylog for y<=0"
        write(*,*)"stop"
        stop
      end if
    endif
    if(n .le. 1 .or. xmin .eq. xmax)then
      write(*,*) "coop function cannot be initialized for xmin = xmax"
      stop
    endif
    if(n.eq.2)then
      this%method = COOP_INTERPOLATE_LINEAR
    elseif(present(method))then
      this%method = method
    endif

    this%n = n
    if(.not. allocated(this%f))then
      allocate(this%f(this%n))
      !      if(this%method .ne. COOP_INTERPOLATE_CHEBYSHEV)then
      if(this%ylog)then
        this%f = log(f)
      else
        this%f = f
        !        endif
      endif
    endif
    allocate(this%f2(this%n))
    !    if(this%method .eq. COOP_INTERPOLATE_CHEBYSHEV) allocate(this%f1(this%n))

    if(this%xlog)then
      if(xmin .le. 0.d0 .or. xmax .le. 0.d0) stop "Error: cannot set xlog = .true. for xmin<0 or xmax<0"
      this%xmin = log(xmin)
      this%xmax = log(xmax)
    else
      this%xmin = xmin
      this%xmax = xmax
    endif
    this%dx = (this%xmax - this%xmin)/(this%n-1)
    if(present(fleft))then
      this%fleft = fleft
    else
      this%fleft = f(1)
    endif
    if(present(fright))then
      this%fright = fright
    else
      this%fright = f(n)
    endif
    if(this%ylog)then
      this%fleft = log(this%fleft)
      this%fright = log(this%fright)
    endif
    if(present(slopeleft))then
      this%slopeleft = slopeleft
    else
      this%slopeleft= 0.d0
    endif
    if(present(sloperight))then
      this%sloperight = sloperight
    else
      this%sloperight = 0.d0
    endif
    select case(this%method)
    case(COOP_INTERPOLATE_LINEAR)
      this%f2 = 0.

    case(COOP_INTERPOLATE_SPLINE)
      call coop_spline_uniform(this%n, this%f, this%f2)
      if(present(smooth))then
        if(smooth)then
          if(this%n .ge. 200)then  !!check f2 is smooth
            call coop_smooth_data_d(this%n, this%f2, min(this%n/200, 50))
          endif
        endif
      endif

    end select
    if(present(check_boundary))then
      this%check_boundary = check_boundary
    else
      this%check_boundary = .true.
    endif
    this%initialized = .true.
  end subroutine coop_function_init

  function coop_function_evaluate_bare(this, x) result(f)
    class(coop_function):: this
    COOP_REAL x, f, a, b, xdiff
    COOP_INT l, r
    select case(this%method)

    case(COOP_INTERPOLATE_POLYNOMIAL)
      f = coop_polyvalue(this%n, this%f, x)
      return

    end select
    xdiff = x - this%xmin
    b = xdiff/this%dx + 1.d0
    l = floor(b)
    if(l .lt. 1)then
      if(this%check_boundary)then
        if(b.gt. 0.9999d0)then
          f = this%fleft
          return
        endif
        write(*,*) this%xlog, this%ylog, this%n, b
        write(*,*) x, ":", this%xmin, " -- ", this%xmax, " ---", this%dx
        write(*,*) "coop_function cannot be evaluated out of its boundary"
        stop
      endif
      f = this%fleft + this%slopeleft*xdiff
    elseif(l.ge. this%n)then
      if(this%check_boundary)then
        if(b .le. dble(this%n)+1.d-4)then
          f = this%fright
          return
        endif
        write(*,*) "function "//trim(this%name)//": boundary overflow (from the left)"
        stop
      endif
      xdiff = x - this%xmax
      f = this%fright + this%sloperight*xdiff
    else
      select case(this%method)
      case(COOP_INTERPOLATE_LINEAR)
        b = b - l
        f = this%f(l) * (1.d0-b) + this%f(l+1) * b
      case(COOP_INTERPOLATE_QUADRATIC, COOP_INTERPOLATE_SPLINE)
        b = b - l
        r = l + 1
        a = 1. - b
        f = this%f(l) * a + this%f(r) * b + this%f2(l) * (a**2-1.)*a + this%f2(r)*(b**2-1.)*b

      case default
        stop "UNKNOWN interpolation method"
      end select
    endif
  end function coop_function_evaluate_bare

  !====================================================================================
!  function coop_function_evaluate_bare(this, x) result(f)
!    class(coop_function):: this
!    COOP_REAL x, f, a, b, xdiff
!    COOP_INT l, r
!    if(this%method==COOP_INTERPOLATE_POLYNOMIAL)then
!      f = coop_polyvalue(this%n, this%f, x)
!      return
!      end if
!      xdiff = x - this%xmin
!      b = xdiff/this%dx + 1.d0
!      l = floor(b)
!      if(l .lt. 1)then
!        if(this%check_boundary)then
!          if(b.gt. 0.9999d0)then
!            f = this%fleft
!            return
!          endif
!          write(*,*) this%xlog, this%ylog, this%n, b
!          write(*,*) x, ":", this%xmin, " -- ", this%xmax, " ---", this%dx
!          write(*,*) "coop_function cannot be evaluated out of its boundary"
!          stop
!        endif
!        f = this%fleft + this%slopeleft*xdiff
!      elseif(l.ge. this%n)then
!        if(this%check_boundary)then
!          if(b .le. dble(this%n)+1.d-4)then
!            f = this%fright
!            return
!          endif
!          write(*,*) "function "//trim(this%name)//": boundary overflow (from the left)"
!          stop
!        endif
!        xdiff = x - this%xmax
!        f = this%fright + this%sloperight*xdiff
!      else
!        select case(this%method)
!        case(COOP_INTERPOLATE_LINEAR)
!          b = b - l
!          f = this%f(l) * (1.d0-b) + this%f(l+1) * b
!        case(COOP_INTERPOLATE_SPLINE)!, COOP_INTERPOLATE_QUADRATIC)
!          b = b - l
!          r = l + 1
!          a = 1. - b
!          f = this%f(l) * a + this%f(r) * b + this%f2(l) * (a**2-1.)*a + this%f2(r)*(b**2-1.)*b
!          !       case(COOP_INTERPOLATE_CHEBYSHEV)
!          !          call coop_chebeval(this%n, this%xmin, this%xmax, this%f, x, f)
!        case default
!          stop "UNKNOWN interpolation method"
!        end select
!      endif
!  end function coop_function_evaluate_bare
!
!  subroutine coop_function_init(this, n, xmin, xmax, f, method, fleft, fright, slopeleft, sloperight, chebyshev_order, xlog, ylog, check_boundary, smooth, name)
!    class(coop_function):: this
!    logical, optional::xlog, ylog
!    COOP_INT,intent(IN):: n
!    COOP_REAL,intent(IN):: xmin, xmax
!    COOP_REAL,intent(IN)::f(n)
!    COOP_REAL, optional::fleft, fright, slopeleft, sloperight
!    COOP_INT, optional::method
!    COOP_INT, optional::chebyshev_order
!    logical,optional::check_boundary
!    logical,optional::smooth
!    COOP_UNKNOWN_STRING, optional::name
!    COOP_INT i, count_tiny, count_small
!    COOP_REAL::fmean, ftiny, curv, flarge, fsmall
!    call this%free()
!    if(all(f .eq. 0.d0))then
!      call this%init_polynomial( (/ 0.d0 /) )
!      return
!    endif
!    if(present(name))then
!      this%name = trim(adjustl(name))
!    endif
!    if(coop_isnan_arrd(f))then
!      write(*,*) "Cannot construct the function "//trim(this%name)//": found f = NAN within the specified range."
!      stop
!    endif
!    if(n.eq.1)then
!      call this%init_polynomial (f)
!      return
!    endif
!    if(present(xlog))then
!      this%xlog = xlog
!    else
!      this%xlog = .false.
!    endif
!    if(present(ylog))then
!      this%ylog = ylog
!    else
!      this%ylog = .false.
!    endif
!    if(this%xlog .and. (xmin .le. 0.d0 .or. xmax .le. 0.d0))then
!      write(*,*)"function:init"
!      write(*,*)"cannot do xlog for x<=0"
!      write(*,*)"stop"
!      stop
!    end if! call coop_return_error("function:init", "cannot do xlog for x<=0", "stop")
!    if(this%ylog)then
!      if(any(f.le.0.d0))then
!        write(*,*)"function:init"
!        write(*,*)"cannot do ylog for y<=0"
!        write(*,*)"stop"
!        stop
!      end if !call coop_return_error("function:init", "cannot do ylog for y<=0", "stop")
!    endif
!
!    if(n.eq.2)then
!      this%method = COOP_INTERPOLATE_LINEAR
!    elseif(present(method))then
!      this%method = method
!    endif
!
!    if(.not. allocated(this%f))then
!      allocate(this%f(this%n))
!      if(this%ylog)then
!        this%f = log(f)
!      else
!        this%f = f
!      endif
!    endif
!    allocate(this%f2(this%n))
!
!    if(this%xlog)then
!      if(xmin .le. 0.d0 .or. xmax .le. 0.d0) stop "Error: cannot set xlog = .true. for xmin<0 or xmax<0"
!      this%xmin = log(xmin)
!      this%xmax = log(xmax)
!    else
!      this%xmin = xmin
!      this%xmax = xmax
!    endif
!    this%dx = (this%xmax - this%xmin)/(this%n-1)
!    if(present(fleft))then
!      this%fleft = fleft
!    else
!      this%fleft = f(1)
!    endif
!    if(present(fright))then
!      this%fright = fright
!    else
!      this%fright = f(n)
!    endif
!    if(this%ylog)then
!      this%fleft = log(this%fleft)
!      this%fright = log(this%fright)
!    endif
!    if(present(slopeleft))then
!      this%slopeleft = slopeleft
!    else
!      this%slopeleft= 0.d0
!    endif
!    if(present(sloperight))then
!      this%sloperight = sloperight
!    else
!      this%sloperight = 0.d0
!    endif
!    select case(this%method)
!    case(COOP_INTERPOLATE_LINEAR)
!      this%f2 = 0.
!    case(COOP_INTERPOLATE_SPLINE)
!      call coop_spline_uniform(this%n, this%f, this%f2)
!      if(present(smooth))then
!        if(smooth)then
!          if(this%n .ge. 200)then  !!check f2 is smooth
!            call coop_smooth_data_d(this%n, this%f2, min(this%n/200, 50))
!          endif
!        endif
!      endif
!    end select
!    if(present(check_boundary))then
!      this%check_boundary = check_boundary
!    else
!      this%check_boundary = .true.
!    endif
!    this%initialized = .true.
!  end subroutine coop_function_init

!====================================================================================
  !end
end module coop_ellipse_collapse_mod
