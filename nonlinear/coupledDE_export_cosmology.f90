program test
  use coop_wrapper_firstorder
  use coop_coupledDE_collapse_mod
  implicit none
#include "constants.h"

#define furen_version

  type(coop_dictionary)::dict
  COOP_INT,parameter::n_threads = 8
  COOP_STRING::params_file, output_root, output_format
  type(coop_coupledDE_collapse_params)::params(n_threads)
  COOP_INT::n_fpk, n_e, n_pbye
  COOP_REAL::fpk_max, fpk_min, e_max, e_min, pbye_min, pbye_max, output_zmax
  COOP_REAL,dimension(:,:,:),allocatable::zcol
  COOP_REAL,dimension(:),allocatable::fpk, e, pbye, z, a
  COOP_INT::ifpk, ie, ipbye, output_nz, ithread, i
  COOP_REAL::H0, Omega_b, Omega_c, A_s, n_s, sigma_8, w, wa

#ifdef furen_version
  integer::pkp_n_fpk, pkp_n_e, pkp_n_pbye,h
  real::pkp_fpk_max, pkp_fpk_min, pkp_e_max, pkp_e_min, pkp_pbye_min, pkp_pbye_max
  COOP_REAL::output_amin,output_amax
  real,dimension(:,:,:),allocatable::pkp_zcol
  real,dimension(:),allocatable::pkp_a, pkp_ttab,pkp_chitab,pkp_dlnD_dlnatab,pkp_Dtab,pkp_Htab
#endif
  !end

  type(coop_file)::fp
  if(iargc().lt. 1)then
     write(*,*) "========================================================"
     write(*,*) "Syntax:"
     write(*,*) "./CDExport params.ini"
     write(*,*) "========================================================"
     stop
  endif

  !!=================import the cosmology from an ini file ===============
  call coop_get_Input(1, params_file)
!  write(*,*)params_file

  call coop_load_dictionary(params_file, dict)
  call coop_dictionary_lookup(dict, "output_root", output_root)

!  write(*,*)output_root

  call coop_dictionary_lookup(dict, "output_format", output_format, "BINARY")
#ifndef furen_version
  select case(trim(output_format))
  case("BINARY")
     call fp%open(trim(output_root)//"_zcol.dat", "s")
  case("FITS")
     stop "FITS format is to be implemented in future versions."
  case default
     stop "This format has not been implemented."
  end select
#endif
  write(*,*)'======================================'
  call params(1)%init(dict, update_cosmology = .true.)
  H0 = params(1)%cosmology%h()*100.
  Omega_b = params(1)%cosmology%Omega_b
  Omega_c = params(1)%cosmology%Omega_c
  A_s = params(1)%cosmology%As
  n_s = params(1)%cosmology%ns
  sigma_8 = params(1)%cosmology%sigma_8
  call coop_dictionary_lookup(dict, "de_w", w, -1.d0)
  call  coop_dictionary_lookup(dict, "de_wa", wa, 0.d0)

  do ithread = 2, n_threads
     params(ithread) = params(1)
  enddo
  !!==================import the halo settings ==========================
  call coop_dictionary_lookup(dict, "fpk_max", fpk_max, 30.d0)
  call coop_dictionary_lookup(dict, "fpk_min", fpk_min, 1.5d0)
  call coop_dictionary_lookup(dict, "e_min", e_min, 0.d0)
  call coop_dictionary_lookup(dict, "e_max", e_max, 0.5d0)
  call coop_dictionary_lookup(dict, "pbye_min", pbye_min, -1.d0)
  call coop_dictionary_lookup(dict, "pbye_max", pbye_max, 1.d0)

!  write(*,*)fpk_max

  call coop_dictionary_lookup(dict, "n_fpk", n_fpk, 50)
  call coop_dictionary_lookup(dict, "n_e", n_e, 30)
  call coop_dictionary_lookup(dict, "n_pbye", n_pbye, 30)
  allocate(zcol(n_fpk, n_e, n_pbye))
  allocate(fpk(n_fpk), e(n_e), pbye(n_pbye))
  call coop_set_uniform(n_fpk, fpk, fpk_min, fpk_max, logscale = .true.)
  call coop_set_uniform(n_e, e, e_min, e_max)
  call coop_set_uniform(n_pbye, pbye, pbye_min, pbye_max)

  write(*,*)"output_root: ", trim(output_root)
  write(*,*)"de_w: ", w
  write(*,*)"de_wa: ", wa
  write(*,*)"fpk_max: ", fpk_max
  write(*,*)"fpk_min: ", fpk_min
  write(*,*)"e_min: ", e_min
  write(*,*)"e_max: ", e_max
  write(*,*)"pbye_min: ", pbye_min
  write(*,*)"pbye_max: ", pbye_max
  write(*,*)"n_fpk: ", n_fpk
  write(*,*)"n_e: ", n_e
  write(*,*)"n_pbye: ", n_pbye

  !$omp parallel do
  do ithread = 1, n_threads
     !write(*,'(f6.2A1)')dble(ithread-1)/dble(n_threads)*100,'%'
     do ipbye = ithread, n_pbye, n_threads
        do ie = 1, n_e
           do ifpk = 1, n_fpk
              call params(ithread)%update_fep(fpk(ifpk), e(ie), pbye(ipbye)*e(ie))
              zcol(ifpk, ie, ipbye) = params(ithread)%zvir1()
#ifdef furen_version
              if(zcol(ifpk, ie, ipbye).gt.-0.5)then
                zcol(ifpk, ie, ipbye)=zcol(ifpk, ie, ipbye)+1.0
              end if
#endif
           enddo
        enddo
     enddo
  enddo
  !$omp end parallel do
  !!for other background functions

#ifdef furen_version
  allocate(pkp_zcol(n_fpk, n_e, n_pbye))
  pkp_zcol=zcol

  pkp_n_fpk = n_fpk
  pkp_n_e = n_e
  pkp_n_pbye = n_pbye

  pkp_fpk_max=fpk_max
  pkp_fpk_min=fpk_min
  pkp_e_max=e_max
  pkp_e_min=e_min
  pkp_pbye_min=pbye_min
  pkp_pbye_max=pbye_max
#endif

#ifdef furen_version

  call coop_dictionary_lookup(dict, "output_amin", output_amin, 0.0099d0)
  call coop_dictionary_lookup(dict, "output_amax", output_amax, 1.d0)
  if(output_amax.gt.1.d0)then
    write(*,*)'this version does not support a > 1!'
    write(*,*)'set amax = 1.0'
    output_amax = 1.d0
  end if
  write(*,*)"output_amin: ", output_amin
  write(*,*)"output_amax: ", output_amax
!  if(output_amin.lt.0.03d0)then
!    write(*,*)'this version does not support a < 0.03'
!    write(*,*)'set amin = 0.03'
!    output_amin = 0.03d0
!  end if

  call coop_dictionary_lookup(dict, "output_na", output_nz, 5001)
  allocate(z(output_nz), a(output_nz))
  call coop_set_uniform(output_nz, a, output_amin, output_amax, logscale = .true.)
  z=1/a

#else
  call coop_dictionary_lookup(dict, "output_zmax", output_zmax, 30.d0)
  call coop_dictionary_lookup(dict, "output_nz", output_nz, 5000)
  allocate(z(output_nz), a(output_nz))
  call coop_set_uniform(output_nz, z, 0.d0, output_zmax)
  a = 1.d0/(1.d0+z)
  write(*,*)"output_na: ", output_nz
#endif
#ifdef furen_version
  open(unit=4,file=trim(output_root)//'/HomelTab.dat',status='replace',access='stream')
  write(4)pkp_n_fpk, pkp_n_e, pkp_n_pbye,log10(pkp_fpk_min),log10(pkp_fpk_max),pkp_e_min,pkp_e_max,pkp_pbye_min,pkp_pbye_max
  write(4)pkp_zcol

  write(*,*)trim(output_root)//'/HomelTab.dat has been written in!'

  close(4)

!read(4)COOP_size
  !    read(4)COOP_anow,COOP_tnow,COOP_taunow,COOP_HD_Ha_now,COOP_Dnow
  !    read(4)COOP_a(1:COOP_size),COOP_ttab(1:COOP_size),COOP_chitab(1:COOP_size),COOP_dlnD_dlnatab(1:COOP_size),&
  !            COOP_Dtab(1:COOP_size),COOP_Htab(1:COOP_size)
  !  real,dimension(:),allocatable::pkp_a, pkp_ttab,pkp_chitab,pkp_dlnD_dlnatab,pkp_Dtab,pkp_Htab

  allocate(pkp_a(output_nz),pkp_ttab(output_nz),pkp_chitab(output_nz),&
          pkp_dlnD_dlnatab(output_nz),pkp_Dtab(output_nz),pkp_Htab(output_nz))

  do i=1,output_nz
    pkp_a(i)=a(i)
    pkp_ttab(i)=params(1)%cosmology%time(a(i))/H0*100.0
    pkp_chitab(i)=params(1)%cosmology%comoving_distance(a(i))/H0*3.0d5
    pkp_Htab(i)=params(1)%dadt(a(i))/a(i)*H0/100.0
    pkp_Dtab(i)=params(1)%growth_D(a(i))
    pkp_dlnD_dlnatab(i)=params(1)%Growth_H_D(a(i))/(params(1)%dadt(a(i))/a(i))
  end do

  open(unit=4,file=trim(output_root)//'/COOP_tables.dat',status='replace',access='stream')
  write(4)output_nz
!  write(4)pkp_a(output_nz),pkp_ttab(output_nz),pkp_chitab(output_nz),&
!          pkp_dlnD_dlnatab(output_nz),pkp_Dtab(output_nz)

  write(4)pkp_a,pkp_ttab,pkp_chitab,pkp_dlnD_dlnatab,pkp_Dtab,pkp_Htab

  write(*,*)trim(output_root)//'/COOP_tables.dat has been written in!'

  close(4)

#else
  select case(trim(output_format))
  case("BINARY")
     write(fp%unit) n_fpk, fpk_min, fpk_max
     write(fp%unit) n_e, e_min, e_max
     write(fp%unit) n_pbye, pbye_min, pbye_max
     write(fp%unit) zcol
     call fp%close()
     write(*,*) "====================================================================="
     write(*,*) "z_collapse(f_pk, e, p/e) is saved to "//trim(output_root)//"_zcol.dat"
     write(*,*) "sample code to read it:"
     write(*,*) "real*8 fpk_min, fpk_max, e_min, e_max, pbye_min, pbye_max"
     write(*,*) "integer n_fpk, n_e, n_pbye"
     write(*,*) "real*8,allocatable::zcol(:,:,:)"
     write(*,*) "open(..., access ='stream')"
     write(*,*) "read(...) n_fpk, fpk_min, fpk_max"
     write(*,*) "read(...) n_e, e_min, e_max"
     write(*,*) "read(...) n_pbye, pbye_min, pbye_max"
     write(*,*) "allocate(zcol(n_fpk, n_e, n_pbye))"
     write(*,*) "read(...) zcol"
     write(*,*) "****************************"
     write(*,*) "fpk is log uniform, e and p/e are both uniform"
     write(*,*) "====================================================================="

     call fp%open(trim(output_root)//"_background.dat", "s")
     write(fp%unit)H0, Omega_b, Omega_c, A_s, n_s, sigma_8, w, wa
     write(fp%unit) output_nz, output_zmax
     do i=1, output_nz
        write(fp%unit) z(i), params(1)%cosmology%time(a(i)),  params(1)%dadt(a(i))/a(i),&
                params(1)%cosmology%comoving_distance(a(i)),  params(1)%growth_D(a(i)), params(1)%Growth_H_D(a(i))
     enddo
     call fp%close()
     write(*,*) "====================================================================="     
     write(*,*) "background info is saved to "//trim(output_root)//"_background.dat"
     write(*,*) "to read it:"
     write(*,*) "real*8 zmax, H0, Omega_b, Omega_c, A_s, n_s, sigma_8, w, wa"
     write(*,*) "integer nz, iz"
     write(*,*) "real*8, dimension(:),allocatable::z, t, H, chi, D, H_D"
     write(*,*) "open(..., access ='stream')"
     write(*,*) "read(...) H0, Omega_b, Omega_c, A_s, n_s, sigma_8, w, wa"
     write(*,*) "read(...) nz, zmax"
     write(*,*) "allocate(z(nz), t(nz), H(nz), chi(nz), D(nz), H_D(nz)"
     write(*,*) "do iz = 1, nz"
     write(*,*) "  read(...) z(iz), t(iz), H(iz), chi(iz), D(iz), H_D(iz))"
     write(*,*) "enddo"
     write(*,*) "***************************************"
     write(*,*) "units: "
     write(*,*) "t (age of the universe):  1/H_0"
     write(*,*) "H (Hubble parameter):  H_0"
     write(*,*) "chi (comoving distance): c/H_0"
     write(*,*) "H_D (d ln D/ dt where D is the growth factor):  H_0  "
     write(*,*) "====================================================================="

  case("FITS")
     stop "FITS format is to be implemented in future versions."
  case default
     stop "This format has not been implemented."
  end select
#endif
#ifdef furen_version
  deallocate(pkp_zcol,pkp_a, pkp_ttab,pkp_chitab,pkp_dlnD_dlnatab,pkp_Dtab,pkp_Htab)
#endif

end program test
