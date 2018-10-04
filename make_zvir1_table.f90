program test
!  use coop_wrapper_firstorder
  use coop_ellipse_collapse_mod
  implicit none
!#include "constants.h"

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

#ifndef furen_version

#define furen_version

#endif


  !!table format:
  !!  (F_pk, e_nu, p_nu, zvir1)

  !variable defined by furen
  logical::if_pkp
  COOP_STRING::params_file
  !cuz the peak-patch read in real(4)
  real::pkp_fmin, pkp_fmax, pkp_p_emin, pkp_p_emax, pkp_emin, pkp_emax
  real,dimension(:,:,:),allocatable::pkp_zvir1
  integer::pkp_nf, pkp_np, pkp_ne
  !end

  COOP_INT, parameter::nthreads = 8
  type(coop_ellipse_collapse_params),dimension(nthreads)::params
  COOP_REAL::Omega_m,w, Omega_k, h, wa, epsilon_s
  COOP_REAL,dimension(:,:,:),allocatable::zvir1
  COOP_REAL,dimension(:),allocatable::f, p, e
  COOP_INT::nf, np, ne, if, ie, ip, ithread, nmax
  COOP_REAL::fmin, fmax, p_emin, p_emax, emin, emax, amax, amin, astart
  logical::logf, binary

  COOP_STRING::output, outfile
  write(*,*)"COOP read parameter from file!"
!  read(*,*)if_pkp
  if_pkp=.true.
  if(if_pkp)then
    binary=.true.
    logf=.true.
!    write(*,*)'Sequence of params:'
!    write(*,*) 'nf, ne, np, fmin, fmax, emin, emax, p_emin, p_emax, omm, omk, h, w, wa, eps, fr1, fr2, fr3, amin, amax, astart, nmax'
!    write(*,*)'input path of the parameters file , the linear file and the nonlinear file:'
    read(*,'(A128)')params_file,outfile,output
    params_file=trim(adjustl(params_file))
    outfile=trim(adjustl(outfile))
    output=trim(adjustl(output))
    write(*,*)'start to read!'
    open(unit=4,file=params_file)
    read(4,*)nf, ne, np, fmin, fmax, emin, emax, p_emin, p_emax, Omega_m, Omega_k, h, w, wa, epsilon_s, &
            params(1)%collapse_a_ratio(1), params(1)%collapse_a_ratio(2), params(1)%collapse_a_ratio(3), amin, amax, astart, nmax
    write(*,*)'finished reading parameters!'
    close(4)
    goto 6677
  end if
  write(*,*)'input parameters:'
  write(*,*)'table, file, nf, ne, np, fmin, fmax, emin, emax, p_emin, p_emax, Omega_m, Omega_k, h, w, wa, epsilon_s'
  write(*,*)'fr1, fr2, fr3, amin, amax, astart, nmax(nmax+1 is the number of element in the array),logf,binary'
  read(*,'(A128)')outfile,output
  read(*,*)nf, ne, np, fmin, fmax, emin, emax, p_emin, p_emax, Omega_m, Omega_k, h, w, wa, epsilon_s, &
          params(1)%collapse_a_ratio(1), params(1)%collapse_a_ratio(2), params(1)%collapse_a_ratio(3), amin, amax, astart, nmax ,logf ,binary


6677 continue

  if(fmin .gt. fmax .or. emin .gt. emax .or. p_emin.gt.p_emax .or. fmin .lt. 0.d0 .or. emin .lt. 0.d0 .or. abs(p_emin).gt. 1 .or. abs(p_emax) .gt. 1)then
     write(*,*) "check error: the input must satisfy "
     write(*,*) "numf>=1 and nume>=1 and nump>=1 and fmin<=fmax and emin<=emax and p_emin<=p_emax and fmin>0 and emin>=0 and abs(p_emin)<=1 and abs(p_emax)<=1"
     stop
  endif
  if(nf*ne*np .gt. 10000 .and. .not. binary)then
     write(*,*) "For large numf x nume x nump force binary format."
     binary = .true.
  endif
  if(nf*ne*np .gt. 100000000)then
     write(*,*) "check your numf x nume x nump; table size is too big"
     stop
  endif

  !subroutine coop_make_table_for_pkp(this, Omega_m, w, wa, epsilon_s, h, Omega_k, lambda, F_pk, e_nu, p_nu,amin, amax,nmax,outfile)!element in a table is nmax+1
  if(if_pkp)then
    call params(1)%make_tab(Omega_m = Omega_m, Omega_k = Omega_k, h = h, w = w, wa=wa, &
              epsilon_s = epsilon_s, amin=amin, amax=amax, nmax=nmax, astart=astart, outfile=outfile)

    write(*,*)"Finished linear Table!"
    write(*,*)"The linear table is dumped to "//trim(outfile)
    write(*,*)"Start to produce Zvir1 nonlinear Table!"

  else
    call params(1)%init(Omega_m = Omega_m, Omega_k = Omega_k, h = h, w = w, wa=wa, epsilon_s = epsilon_s)
  end if

  allocate(f(nf), p(np), e(ne), zvir1(nf, ne, np), pkp_zvir1(nf, ne, np))
  call coop_set_uniform_pkp(nf, f, fmin, fmax, logscale = logf)
  call coop_set_uniform_pkp(ne, e, emin, emax)
  call coop_set_uniform_pkp(np, p, p_emin, p_emax)
!  call coop_prtsystime(.true.)
  if( nf .gt. nthreads)then
     do ithread = 2, nthreads
        params(ithread) = params(1)
     enddo
     !$omp parallel do private(ithread, if, ie, ip)
     do ithread = 1, nthreads
       write(*,'(f6.2A1)')dble(ithread-1)/dble(nthreads)*100,'%'
        do if = ithread, nf, nthreads
           do ie = 1, ne
              do ip = 1, np
                 if(abs(p(ip)).gt. 1)then
                    zvir1(if, ie, ip) = coop_ellipse_collapse_bad_zvir
                    write(*,*)'Warning! |p/e| > 1!'
                 else
                    if(if_pkp)then
                      call params(ithread)%pkp_init(F_pk = f(if), e_nu = e(ie), p_nu=p(ip)*e(ie))
                    else
                      call params(ithread)%init(F_pk = f(if), e_nu = e(ie), p_nu=p(ip)*e(ie))
                    end if
		    !add for check
		    !if(if.eq.11.and.ie.eq.16.and.ip.eq.20)then
		    !	write(*,*)params(ithread)%zvir1(),f(if),e(ie),p(ip)*e(ie)
		    !	write(*,*)params(ithread)%Dbya%method
             !   write(*,*)params(ithread)%Dbyadot%eval(dble(1.3))
             !   write(*,*)params(ithread)%bprime%eval(dble(0.56),dble(0.23))
            !endif
		    !end
                    zvir1(if, ie, ip) = params(ithread)%zvir1()
		    if(zvir1(if, ie, ip).lt.-0.5)then
		        zvir1(if, ie, ip)=coop_ellipse_collapse_bad_zvir
		    else
		        zvir1(if, ie, ip)=zvir1(if, ie, ip)+1
		    endif
                 endif
              enddo
           enddo
        enddo
     enddo
     !$omp end parallel do
  else
     do if = 1, nf
       write(*,'(f6.2A1)')dble(if)/dble(nf)*100,'%'
        do ie = 1, ne
           do ip = 1, np
              if(abs(p(ip)).gt. 1)then
                 zvir1(if, ie, ip) = coop_ellipse_collapse_bad_zvir
              else
                if(if_pkp)then
                  call params(1)%pkp_init(F_pk = f(if), e_nu = e(ie), p_nu=p(ip)*e(ie))
                else
                  call params(1)%init(F_pk = f(if), e_nu = e(ie), p_nu=p(ip)*e(ie))

                end if
                 zvir1(if, ie, ip) = params(1)%zvir1()
		 if(zvir1(if, ie, ip).lt.0.0)then
		      zvir1(if, ie, ip)=coop_ellipse_collapse_bad_zvir
		 else
		      zvir1(if, ie, ip)=zvir1(if, ie, ip)+1
		 endif
              endif
           enddo
        enddo
     enddo
  endif
!  call coop_prtsystime()

  if(binary)then
!     call fp%open(output, "u")
     if(if_pkp)then
         pkp_nf=nf
         pkp_np=np
         pkp_ne=ne
         pkp_fmin=log10(fmin)
         pkp_fmax=log10(fmax)
         pkp_emin=emin
         pkp_emax=emax
         pkp_p_emin=p_emin
         pkp_p_emax=p_emax
         pkp_zvir1=zvir1
         open(unit=4,file=output,access='stream',status='replace')
         write(4) pkp_nf, pkp_ne, pkp_np, pkp_fmin, pkp_fmax, pkp_emin, pkp_emax, pkp_p_emin, pkp_p_emax, pkp_zvir1(1:nf, 1:ne, 1:np)
         deallocate(pkp_zvir1)
         close(4)
     else
         open(unit=4,file=output,access='stream',status='replace')
         write(4) nf, ne, np, log10(fmin), log10(fmax), emin, emax, p_emin, p_emax, zvir1(1:nf, 1:ne, 1:np)
         close(4)

     endif

  else
     open(unit=4,file=output,status='replace')
     write(4, *) "#  F_pk         ",  "    e_nu  ", "   p_nu   ", "  z_1  "
     do if = 1, nf
        do ie = 1, ne
           do ip = 1, np
              if(zvir1(if, ie, ip) .ge. 0.d0)then
                 write(4, *) f(if), e(ie), p(ip), zvir1(if, ie, ip)
              endif
           enddo
        enddo
     enddo
     close(4)
  endif
  write(*,*) "The nonlinear table is dumped to "//trim(output)

end program test
