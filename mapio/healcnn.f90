module coop_healcnn_mod
  use coop_wrapper_utils
  use coop_sphere_mod
  use coop_fitsio_mod
  use coop_healpix_mod
  use coop_stacking_mod
  use coop_fitswrap_mod
  use coop_gstack_mod
  implicit none
#include "constants.h"
#define ACTIVATE_FUNCTION 0
  !! 0 = sigmoid
  !! 1 = ReLU (Rectified Linear Unit)
  !! 2 = tanh
  !! 3 = modified ReLU

  type coop_healcnn_connection
     COOP_INT::imap_in = 0
     COOP_INT::imap_out = 0     
     COOP_INT::nb = 0
     COOP_INT::nx = 0     
     COOP_INT,dimension(:),allocatable::lmin
     COOP_INT,dimension(:),allocatable::lmax
     COOP_REAL,dimension(:),allocatable::w
     COOP_REAL,dimension(:),allocatable::dw     
     COOP_REAL,dimension(:),allocatable::x
     COOP_REAL,dimension(:,:),allocatable::P
   contains
     procedure::init => coop_healcnn_connection_init
     procedure::default => coop_healcnn_connection_default     
     procedure::free => coop_healcnn_connection_free
     procedure::kernel => coop_healcnn_connection_kernel
  end type coop_healcnn_connection


  type coop_healcnn_layer
     type(coop_healpix_maps)::h
     type(coop_healcnn_connection),dimension(:),allocatable::c
     COOP_REAL,dimension(:),allocatable::bias, db
     COOP_INT::nside_in, nside_pooling, nside_out, nmaps_in, nmaps_out, nc

   contains
     procedure::free => coop_healcnn_layer_free
     procedure::init => coop_healcnn_layer_init
  end type coop_healcnn_layer

  type coop_healcnn
     COOP_INT::nlayers = 0
     type(coop_healcnn_layer),dimension(:),allocatable::layers
     type(coop_healpix_maps)::mask
   contains
     procedure::free => coop_healcnn_free
     procedure::init => coop_healcnn_init     
  end type coop_healcnn

  

contains


  subroutine coop_healcnn_init(this, nlayers, map, mask)
    class(coop_healcnn)::this
    COOP_INT::nlayers
    COOP_UNKNOWN_STRING::map, mask
    this%nlayers = nlayers
    allocate(this%layers(nlayers))
    call this%layers(1)%h%read(filename=map)
    call this%mask%read(filename = mask)
  end subroutine coop_healcnn_init
  

  subroutine coop_healcnn_free(this)
    class(coop_healcnn)::this
    COOP_INT::i
    call this%mask%free()
    do i=1, this%nlayers
       call this%layers(i)%free()
    enddo
    this%nlayers = 0
  end subroutine coop_healcnn_free

  subroutine coop_healcnn_layer_free(this)
    class(coop_healcnn_layer)::this
    call this%h%free()
    COOP_DEALLOC(this%c)
    COOP_DEALLOC(this%bias)
    COOP_DEALLOC(this%db)    
  end subroutine coop_healcnn_layer_free

  subroutine coop_healcnn_layer_init(this, nside_in, nside_pooling, nside_out, nmaps_in, nmaps_out, nc)
    class(coop_healcnn_layer)::this
    COOP_INT::nside_in, nside_pooling, nside_out, nmaps_in, nmaps_out
    COOP_INT,optional::nc
    COOP_INT::i, j, k
    call this%free()
    this%nside_in = nside_in
    this%nside_pooling = nside_pooling    
    this%nside_out = nside_out    
    this%nmaps_in = nmaps_in
    this%nmaps_out = nmaps_out
    call this%h%init(nside = nside_in, nmaps = nmaps_in, genre="UNKNOWN")
    allocate(this%bias(nmaps_out), this%db(nmaps_out))
    if(present(nc))then
       this%nc = nc
       allocate(this%c(this%nc))
    else !!full connections
       this%nc = nmaps_in * nmaps_out
       allocate(this%c(this%nc))       
       k = 0
       do i=1, nmaps_in
          do j=1, nmaps_out
             k = k + 1
             call this%c(k)%default(this%nside_out)
             this%c(k)%imap_in = i
             this%c(k)%imap_out = j
          enddo
       enddo
    endif
  end subroutine coop_healcnn_layer_init

  subroutine coop_healcnn_connection_free(this)
    class(coop_healcnn_connection)::this
    COOP_DEALLOC(this%lmin)
    COOP_DEALLOC(this%lmax)    
    COOP_DEALLOC(this%w)
    COOP_DEALLOC(this%dw)    
    COOP_DEALLOC(this%x)
    COOP_DEALLOC(this%P)
    this%nb = 0
    this%nx = 0
    this%imap_in = 0
    this%imap_out = 0        
  end subroutine coop_healcnn_connection_free

  subroutine coop_healcnn_connection_default(this, nside)
    class(coop_healcnn_connection)::this
    COOP_REAL,dimension(:),allocatable::x
    COOP_INT,dimension(:),allocatable::ells    
    COOP_INT:: nside, nx, lmax, nl, i
    lmax = min(coop_healpix_default_lmax, floor(nside*coop_healpix_lmax_by_nside))
    nl = lmax/10
    nx = max(lmax * 5, 100)
    allocate(ells(nl))
    allocate(x(nx))
    call coop_set_uniform(nx, x, coop_pio2, 0.d0)
    do i=1, nl-1
       ells(i) = 10*i
    enddo
    ells(nl) = lmax
    call this%init(ells, cos(x))
    deallocate(x)
    deallocate(ells)
  end subroutine coop_healcnn_connection_default
  
  subroutine coop_healcnn_connection_init(this, lmax, x)
    class(coop_healcnn_connection)::this
    COOP_INT,dimension(:)::lmax
    COOP_REAL,dimension(:):: x
    COOP_INT::ix, ib, l, LL
    COOP_REAL,dimension(:),allocatable::Pls,cl1, cl2
    COOP_REAL::asymp, dis
    call this%free()
    this%nb = size(lmax)
    this%nx = size(x)
    allocate(this%x(this%nx), this%lmin(this%nb), this%lmax(this%nb), this%w(this%nb), this%dw(this%nb), this%P(this%nx, this%nb))
    this%x = x 
    this%lmax = lmax
    LL = lmax(this%nb)
    if(LL .lt. 2) stop "coop_healcnn_connection_init: does not support lmax < 2"        
    allocate(Pls(0: LL), cl1(0:LL), cl2(0:LL))
    do l=1,LL-1
       cl1(l) = l/(2.d0*l-1.d0)
       cl2(l) = (2.d0*l+3.d0)/(l+1)
    enddo
    this%lmin(1) = 0
    this%lmin(2:this%nb) = lmax(1:this%nb-1)+1
    this%w = 1.d0
    this%dw = 0.d0    
    asymp = 0.2d0 / LL ** 2
    do ix = 1, this%nx
       !! calculate Pls(l) = (2l+1) P_l(x)
       dis = (1.d0-x(ix))/2.d0
       Pls(0) = 1.d0
       Pls(1) = 3.d0*x(ix) 
       if(dis .gt. asymp)then  !!recurrence relation
          do l = 1, LL-1
             Pls(l+1) = (x(ix)*Pls(l) - cl1(l)*Pls(l-1))*cl2(l)
          enddo
       else  !!use approximation
          do l = 2, LL
             Pls(l) = (1.d0 - dis*l*(l+1)*(1.d0-dis*(l-1)*(l+2)/4.d0))*(2*l+1)
          enddo
       endif
       do ib = 1, this%nb
          this%p(ix, ib) = sum(Pls(this%lmin(ib):this%lmax(ib)))
       enddo
       
    enddo
    deallocate(Pls,cl1,cl2)
  end subroutine coop_healcnn_connection_init


  function coop_healcnn_connection_kernel(this, ib, x) result(k)
    class(coop_healcnn_connection)::this    
    COOP_REAL::k,x, r
    COOP_INT::ib, il, ir
    if(x.lt. this%x(1))then !!too far away, ignore correlation
       k = 0.d0
       return
    endif
    if(x.ge.this%x(this%nx))then
       k = this%p(this%nx, ib)
       return
    endif
    il = coop_left_index(this%nx, this%x, x)
    ir = il + 1
    r = (this%x(ir)-x)/(this%x(ir)-this%x(il))
    k = this%p(il, ib)*r + this%p(ir, ib)*(1.d0-r)
  end function coop_healcnn_connection_kernel


#if ACTIVATE_FUNCTION == 0 
  function coop_healcnn_activate(x) result(f)
    COOP_REAL::x, f
    f = 1.d0/(1.d0+exp(-x))
  end function coop_healcnn_activate

  function coop_healcnn_activate_derv(f) result(df)
    COOP_REAL::f,df
    df  = f*(1.d0-f)
  end function coop_healcnn_activate_derv

#elif ACTIVATE_FUNCTION == 1

  function coop_healcnn_activate(x) result(f)
    COOP_REAL::x, f
    f = max(0.d0, x)
  end function coop_healcnn_activate

  function coop_healcnn_activate_derv(f) result(df)
    COOP_REAL::f,df
    if(f.gt.0.d0)then
       df = 1.d0
    else
       df = 0.d0
    endif
  end function coop_healcnn_activate_derv

#elif ACTIVATE_FUNCTION == 2

  function coop_healcnn_activate(x) result(f)
    COOP_REAL::x, f
    f = tanh(x)
  end function coop_healcnn_activate

  function coop_healcnn_activate_derv(f) result(df)
    COOP_REAL::f,df
    df = 1.d0-f**2
  end function coop_healcnn_activate_derv

#elif ACTIVATE_FUNCTION == 3

  function coop_healcnn_activate(x) result(f)
    COOP_REAL::x, f
    if(x.gt.0.d0)then
       f = x
    else
       f = x/10.d0
    endif
  end function coop_healcnn_activate

  function coop_healcnn_activate_derv(f) result(df)
    COOP_REAL::f,df
    if(f.gt.0.d0)then
       df = 1.d0
    else
       df = 0.1d0
    endif
  end function coop_healcnn_activate_derv
  
#endif  

  !!==============================
  
  
end module coop_healcnn_mod


