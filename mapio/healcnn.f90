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

  COOP_INT,parameter::coop_healcnn_default_delta_ell = 10

  
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
     type(coop_healpix_maps)::dh
     type(coop_healpix_maps)::pooling     
     type(coop_healcnn_connection),dimension(:),allocatable::c
     COOP_REAL,dimension(:),allocatable::bias, db
     COOP_INT::nside_in,  nside_out, nmaps_in, nmaps_out, nc

   contains
     procedure::free => coop_healcnn_layer_free
     procedure::init => coop_healcnn_layer_init
     procedure::read => coop_healcnn_layer_read
     procedure::propogate => coop_healcnn_layer_propogate
     procedure::chain_derv => coop_healcnn_layer_chain_derv     
  end type coop_healcnn_layer

  type coop_healcnn
     COOP_INT::nlayers = 0
     type(coop_healcnn_layer),dimension(:),allocatable::layers
     type(coop_healpix_maps)::mask
   contains
     procedure::free => coop_healcnn_free
     procedure::init => coop_healcnn_init
     procedure::fp => coop_healcnn_fp
  end type coop_healcnn


  

contains


  subroutine  coop_healcnn_fp(this)
    class(coop_healcnn)::this
    COOP_INT::il
    do il=1, this%nlayers-1
       call this%layers(il)%propogate(next = this%layers(il+1))
    enddo
  end subroutine coop_healcnn_fp


  subroutine  coop_healcnn_bp(this)
    class(coop_healcnn)::this
    COOP_INT::il
    call this%layers(this%nlayers)%h%apply_mask(this%mask)
    this%layers(this%nlayers)%dh%map = this%layers(this%nlayers)%h%map
    do il = this%nlayers-1, 1, -1
       call this%layers(il)%chain_derv(next = this%layers(il+1))
    enddo
  end subroutine coop_healcnn_bp


  subroutine coop_healcnn_layer_chain_derv(this, next)
    class(coop_healcnn_layer)::this
    type(coop_healcnn_layer)::next
    type(coop_healpix_maps)::tmp
    COOP_REAL,parameter::radius = coop_SI_degree*15.d0
    COOP_INT::imap, i, ic, ib,  l, LL, pix
    COOP_INT::nlist
    COOP_REAL::vec1(3), vec2(3)
    COOP_INT,dimension(:),allocatable::listpix
    COOP_REAL::theta, phi
    if(next%h%nmaps .ne. this%nmaps_out) call coop_return_error("coop_healcnn_layer_chain_derv", "nmaps_out does not match", "stop")
    if(next%h%nside .eq. this%nside_out)then
       do i=1, this%nmaps_out
          this%db(i) = sum(next%dh%map(:, i))
       enddo
       this%dh%map = 0.
       do ic=1, this%nc
          this%c(ic)%dw = 0.d0
       enddo
       call tmp%init(nside = this%nside_out, nmaps = 1, genre="UNKNOWN")
       call tmp%allocate_alms()
       tmp%alm = 0.
       LL = min(tmp%lmax, this%h%lmax)
       allocate(listpix(0:this%h%npix-1))
       do ic=1, this%nc          
          do ib=1, this%c(ic)%nb
             do l = this%c(ic)%lmin(ib), this%c(ic)%lmax(ib)
                if(l .gt. LL)exit
                tmp%alm(l, 0:l, 1)=this%h%alm(l, 0:l, this%c(ic)%imap_in)
             enddo
             call tmp%alm2map()
             do i=0, tmp%npix-1
                this%c(ic)%dw(ib) = this%c(ic)%dw(ib) + next%dh%map(i, this%c(ic)%imap_out) * tmp%map(i,1) * coop_healcnn_activate_derv(next%h%map(i, this%c(ic)%imap_out))
             enddo
             do l = this%c(ic)%lmin(ib), this%c(ic)%lmax(ib)
                if(l .gt. LL)exit                
                tmp%alm(l, 0:l, 1)=0.
             enddo
          enddo
       enddo
       do pix = 0, next%h%npix-1
          call next%h%pix2ang(pix, theta, phi)
          call this%h%ang2pix(theta, phi, i)
          call this%h%query_disc(i, radius, listpix, nlist)
          call this%h%pix2vec(i, vec1)
          do ic=1, this%nc          
             do ib=1, this%c(ic)%nb
                do i=0, nlist-1
                   call this%h%pix2vec(listpix(i), vec2)
                   this%dh%map(listpix(i), this%c(ic)%imap_in) = this%dh%map(listpix(i), this%c(ic)%imap_in) + next%dh%map(pix, this%c(ic)%imap_out) * coop_healcnn_activate_derv(next%h%map(pix, this%c(ic)%imap_out))*this%c(ic)%kernel(ib, dot_product(vec1, vec2))*this%c(ic)%w(ib)
                enddo
             enddo
          enddo
       enddo
       this%dh%map = this%dh%map/ this%dh%npix
       deallocate(listpix)
       call tmp%free()
    elseif(next%h%nside .lt. this%nside_out)then !!do max pooling       
    endif
  end subroutine coop_healcnn_layer_chain_derv
  

  subroutine coop_healcnn_layer_propogate(this, next)
    class(coop_healcnn_layer)::this
    type(coop_healcnn_layer)::next
    COOP_INT::ic,ib, i, j, l, thislmax
    if(next%h%nmaps .ne. this%nmaps_out) call coop_return_error("coop_healcnn_layer_propogate", "nmaps_out does not match", "stop")
    call this%h%map2alm()
    call next%h%allocate_alms()
    if(next%h%nside .eq. this%nside_out)then
       thislmax = min(this%h%lmax, next%h%lmax)
       do ic = 1, this%nc
          do ib=1, this%c(ic)%nb
             if(this%c(ic)%lmin(ib) .gt. thislmax)exit
             do l=this%c(ic)%lmin(ib), this%c(ic)%lmax(ib)
                if(l .gt. thislmax) exit
                next%h%alm(l, 0:l,  this%c(ic)%imap_out) = this%h%alm(l, 0:l,  this%c(ic)%imap_in)*this%c(ic)%w(ib)
             enddo
          enddo
       enddo
       call next%h%alm2map()
       do j=1, next%h%nmaps
          do i=0, next%h%npix-1
             next%h%map(i,j) = coop_healcnn_activate(next%h%map(i,j))
          enddo
       enddo
    elseif(next%h%nside .lt. this%nside_out)then !!do max pooling
       call this%pooling%init(nside = this%nside_out, nmaps = this%nmaps_out, genre="UNKNOWN")
       call this%pooling%allocate_alms()
       thislmax = min(this%h%lmax, this%pooling%lmax)
       do ic = 1, this%nc
          do ib=1, this%c(ic)%nb
             if(this%c(ic)%lmin(ib) .gt. thislmax)exit
             do l=this%c(ic)%lmin(ib), this%c(ic)%lmax(ib)
                if(l .gt. thislmax) exit
                this%pooling%alm(l, 0:l,  this%c(ic)%imap_out) = this%h%alm(l, 0:l,  this%c(ic)%imap_in)*this%c(ic)%w(ib)
             enddo
          enddo
       enddo
       call this%pooling%alm2map()
       call this%pooling%write("this%pooling.fits")
       do j=1, this%pooling%nmaps
          do i=0, this%pooling%npix-1
             this%pooling%map(i,j) = coop_healcnn_activate(this%pooling%map(i,j))
          enddo
       enddo
       call coop_healpix_maps_max_udgrade(from = this%pooling, to = next%h)
    else
       call coop_return_error("coop_healcnn_layer_propogate", "pooling nside < output nside is not supported","stop") 
    endif
  end subroutine coop_healcnn_layer_propogate


  subroutine coop_healcnn_init(this, nlayers, map, mask, nmaps, nside, nside_pooling)
    class(coop_healcnn)::this
    COOP_INT::nlayers
    COOP_UNKNOWN_STRING::map, mask
    COOP_INT::nmaps(nlayers), nside(nlayers)
    COOP_INT,optional::nside_pooling(nlayers-1)
    COOP_INT::il
    this%nlayers = nlayers
    allocate(this%layers(nlayers))
    call this%mask%read(filename = mask)    
    if(present(nside_pooling))then
       call this%layers(1)%read(filename=map, nmaps_in=nmaps(1), nside_in = nside(1), nmaps_out = nmaps(2), nside_out = nside_pooling(1))
       do il = 2, this%nlayers-1
          call this%layers(il)%init(nside_in = nside(il), nside_out = nside_pooling(il), nmaps_in = nmaps(il), nmaps_out = nmaps(il+1))
       enddo       
    else
       call this%layers(1)%read(filename=map, nmaps_in=nmaps(1), nside_in = nside(1), nmaps_out = nmaps(2), nside_out = nside(2))       
       do il = 2, this%nlayers-1
          call this%layers(il)%init(nside_in = nside(il), nside_out = nside(il+1), nmaps_in = nmaps(il), nmaps_out = nmaps(il+1))
       enddo
    endif
    call this%layers(nlayers)%init(nside_in = nside(nlayers), nside_out = 4, nmaps_in = nmaps(nlayers), nmaps_out = 1) 
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
    call this%dh%free()
    call this%pooling%free()
    COOP_DEALLOC(this%c)
    COOP_DEALLOC(this%bias)
    COOP_DEALLOC(this%db)    
  end subroutine coop_healcnn_layer_free

  

  subroutine coop_healcnn_layer_init(this, nside_in,  nside_out, nmaps_in, nmaps_out, nc)
    class(coop_healcnn_layer)::this
    COOP_INT::nside_in,  nside_out, nmaps_in, nmaps_out
    COOP_INT,optional::nc
    COOP_INT::i, j, k
    call this%free()
    this%nside_in = nside_in
    this%nside_out = nside_out    
    this%nmaps_in = nmaps_in
    this%nmaps_out = nmaps_out
    
    call this%h%init(nside = nside_in, nmaps = nmaps_in, genre="UNKNOWN")
    this%dh = this%h
    allocate(this%bias(nmaps_out), this%db(nmaps_out))
    this%bias = 0.
    this%db=0.
    
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
             call this%c(k)%default(this%nside_out, n = nmaps_out, k = j, delta_ell = coop_healcnn_default_delta_ell)
             this%c(k)%imap_in = i
             this%c(k)%imap_out = j
          enddo
       enddo
    endif
  end subroutine coop_healcnn_layer_init


  subroutine coop_healcnn_layer_read(this, filename,  nmaps_in, nside_in, nside_out, nmaps_out, nc)
    class(coop_healcnn_layer)::this    
    COOP_UNKNOWN_STRING::filename
    COOP_INT::nside_out, nmaps_out, nmaps_in, nside_in
    COOP_INT,optional::nc
    COOP_INT::i, j, k
    call this%free()
    call this%h%read(filename = filename, nmaps_wanted = nmaps_in)
    this%dh = this%h
    if(this%h%nside .ne. nside_in) call this%h%udgrade(nside = nside_in)
    
    this%nside_in = this%h%nside
    this%nside_out = nside_out    
    this%nmaps_in = this%h%nmaps
    this%nmaps_out = nmaps_out

    allocate(this%bias(nmaps_out), this%db(nmaps_out))
    this%bias = 0.
    this%db=0.
    if(present(nc))then
       this%nc = nc
       allocate(this%c(this%nc))
    else !!full connections
       this%nc = this%nmaps_in * this%nmaps_out
       allocate(this%c(this%nc))       
       k = 0
       do i=1, this%nmaps_in
          do j=1, this%nmaps_out
             k = k + 1
             call this%c(k)%default(nside = this%nside_out, n = this%nmaps_out, k = j, delta_ell = coop_healcnn_default_delta_ell)
             this%c(k)%imap_in = i
             this%c(k)%imap_out = j
          enddo
       enddo
    endif
  end subroutine coop_healcnn_layer_read
  

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

  subroutine coop_healcnn_connection_default(this, nside, n, k, delta_ell)
    class(coop_healcnn_connection)::this
    COOP_REAL,dimension(:),allocatable::x
    COOP_INT,dimension(:),allocatable::ells    
    COOP_INT:: nside, nx, lmax, nl, i, lmin, dell
    COOP_INT,optional::n, k, delta_ell
    lmax = floor(nside*coop_healpix_lmax_by_nside)+1    
    if(present(n).and.present(k))then
       lmin = lmax*(k-1)/n       
       lmax = lmax*k/n 
    else
       lmin = 0
    endif
    if(present(delta_ell))then
       nl = max((lmax-lmin)/delta_ell, 1)       
    else
       nl = max((lmax-lmin)/15, 1)
    endif
    nx = max(lmax * 5, 100)
    dell = nint((lmax-lmin+1.d0)/nl)
    allocate(ells(0:nl))
    allocate(x(nx))
    call coop_set_uniform(nx, x, coop_pio2, 0.d0)
    ells(0) = lmin
    do i=1, nl-1
       ells(i) = ells(i-1)+dell
    enddo
    ells(nl) = lmax
    call this%init(ells, cos(x))
    deallocate(x)
    deallocate(ells)
  end subroutine coop_healcnn_connection_default
  
  subroutine coop_healcnn_connection_init(this, ells, x)
    class(coop_healcnn_connection)::this
    COOP_INT,dimension(:)::ells
    COOP_REAL,dimension(:):: x
    COOP_INT::ix, ib, l, LL
    COOP_REAL,dimension(:),allocatable::Pls,cl1, cl2
    COOP_REAL::asymp, dis
    call this%free()
    this%nb = size(ells)-1
    this%nx = size(x)
    allocate(this%x(this%nx), this%lmin(this%nb), this%lmax(this%nb), this%w(this%nb), this%dw(this%nb), this%P(this%nx, this%nb))
    this%x = x 
    this%lmax = ells(2:this%nb+1)-1
    this%lmin = ells(1:this%nb)
    LL = this%lmax(this%nb)
    if(LL .lt. 2) call coop_return_error("coop_healcnn_connection_init","does not support lmax < 2", "stop")
    allocate(Pls(0: LL), cl1(0:LL), cl2(0:LL))
    do l=1,LL-1
       cl1(l) = l/(2.d0*l-1.d0)
       cl2(l) = (2.d0*l+3.d0)/(l+1)
    enddo
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
    COOP_SINGLE::x, f
    f = 1.d0/(1.d0+exp(-x))
  end function coop_healcnn_activate

  function coop_healcnn_activate_derv(f) result(df)
    COOP_SINGLE::f,df
    df  = f*(1.d0-f)
  end function coop_healcnn_activate_derv

#elif ACTIVATE_FUNCTION == 1

  function coop_healcnn_activate(x) result(f)
    COOP_SINGLE::x, f
    f = max(0.d0, x)
  end function coop_healcnn_activate

  function coop_healcnn_activate_derv(f) result(df)
    COOP_SINGLE::f,df
    if(f.gt.0.d0)then
       df = 1.d0
    else
       df = 0.d0
    endif
  end function coop_healcnn_activate_derv

#elif ACTIVATE_FUNCTION == 2

  function coop_healcnn_activate(x) result(f)
    COOP_SINGLE::x, f
    f = tanh(x)
  end function coop_healcnn_activate

  function coop_healcnn_activate_derv(f) result(df)
    COOP_SINGLE::f,df
    df = 1.d0-f**2
  end function coop_healcnn_activate_derv

#elif ACTIVATE_FUNCTION == 3

  function coop_healcnn_activate(x) result(f)
    COOP_SINGLE::x, f
    if(x.gt.0.d0)then
       f = x
    else
       f = x/10.d0
    endif
  end function coop_healcnn_activate

  function coop_healcnn_activate_derv(f) result(df)
    COOP_SINGLE::f,df
    if(f.gt.0.d0)then
       df = 1.d0
    else
       df = 0.1d0
    endif
  end function coop_healcnn_activate_derv
  
#endif  

  !!==============================
  
  
end module coop_healcnn_mod


