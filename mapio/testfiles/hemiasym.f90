program test
  use coop_wrapper_utils
  use coop_healpix_mod
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::map, mask, m2
  type(coop_healpix_inpaint)::inp
  COOP_INT,parameter ::lmax = 200
  COOP_INT,parameter ::nrun = 300
  COOP_STRING::mask_spot = ""
  COOP_REAL::cls(0:lmax), Cls_sim(0:lmax, nrun), Cls_ave(0:lmax), delta_Cls(0:lmax), sigma
  COOP_INT::l, ell, i, irun
  COOP_REAL::theta, phi, theta1, theta2
  type(coop_file)::fp
  call coop_MPI_init()
  call coop_random_init()
  !!read in Cl's for fiducial LCDM model
  call fp%open("planck14best_lensedCls.dat", "r") 
  do l=2, lmax
     read(fp%unit, *) i, Cls(l)
     Cls(l) = Cls(l)*coop_2pi/l/(l+1.d0)*coop_gaussian_filter(60.d0, l)**2
  enddo  
  call fp%close()
  Cls(0:1) = 0.d0

  call map%read("lowl/commander_I_n0128_60a.fits")
  call mask%read("lowl/commander_mask_n0128_60a.fits")
  call coop_get_command_line_argument(key = "mask",  arg = mask_spot)
  call coop_get_command_line_argument(key = "theta1",  arg = theta1, default = 0.d0)
  call coop_get_command_line_argument(key = "theta2",  arg = theta2, default = 90.d0)  
  
  mask_spot = adjustl(mask_spot)

  !!mask out a disk
  select case(trim(mask_spot))
  case("COLDSPOT")
     call mask%mask_strip(l_deg = 207.8d0, b_deg = -56.3d0, r1_deg = theta1, r2_deg = theta2)
  case("NGP")
     call mask%mask_strip(l_deg = 0.d0, b_deg = 90.d0, r1_deg = theta1, r2_deg = theta2)
  case("SGP")
     call mask%mask_strip(l_deg = 0.d0, b_deg = -90.d0, r1_deg = theta1, r2_deg = theta2)
  case("NEP") 
     call mask%mask_strip(l_deg = 98.d0, b_deg = 31.d0, r1_deg = theta1, r2_deg = theta2)
  case("SEP")
     call mask%mask_strip(l_deg = 278.d0, b_deg = -31.d0, r1_deg = theta1, r2_deg = theta2)
  case("NCP")
     call mask%mask_strip(l_deg = 123.d0, b_deg = 28.d0, r1_deg = theta1, r2_deg = theta2)     
  case("SCP")
     call mask%mask_strip(l_deg = 303.d0, b_deg = -28.d0, r1_deg = theta1, r2_deg = theta2)
  case("NDP")
     call mask%mask_strip(l_deg = 263.85d0, b_deg = 48.25d0, r1_deg = theta1, r2_deg = theta2)
  case("SDP")
     call mask%mask_strip(l_deg = 83.85d0, b_deg = -48.25d0, r1_deg = theta1, r2_deg = theta2)               
  case("NASYM")
     call mask%mask_strip(l_deg = 212.d0, b_deg = -13.d0, r1_deg = theta1, r2_deg = theta2)
  case("SASYM")
     call mask%mask_strip(l_deg = 32.d0, b_deg = 13.d0,  r1_deg = theta1, r2_deg = theta2)
  case("NS")
     call mask%mask_strip(l_deg = 212.d0, b_deg = -13.d0,  r1_deg = 0.d0, r2_deg = theta1)
     call mask%mask_strip(l_deg = 32.d0, b_deg = 13.d0,  r1_deg = 0.d0, r2_deg = theta2)     
  case("NONE")
     !do nothing
  case default
     write(*,*) trim(mask_spot)//": unknown mask option"
     stop
  end select

  !!initialize
  !!Compute pixel-space covariance matrix from fiducial Cl's
  !!For more serious applications you should use a full covariance matrix from FFP9 simulations  
  call inp%init(map, mask, lmax, cls)

  Cls_sim = 0.d0
  Cls_ave = 0.d0
  do irun = 1, nrun
     write(*,*) "***** inpainting # ", irun, " ***********"
     call inp%upgrade(reset = .true., nside_want = inp%map%nside)
     inp%lMT%map = inp%lMT%map  + inp%lCT%map !!measured map + inpainted map
     call inp%lMT%map2alm(lmax = lmax)
     Cls_sim(2:lmax, irun) = inp%lMT%Cl(2:lmax, 1)
  enddo
  do l = 2, lmax
     Cls_ave(l) = sum(Cls_sim(l, 1:nrun))/nrun
     delta_Cls(l)  = sqrt(sum((Cls_sim(l, 1:nrun)-Cls_ave(l))**2)/nrun + Cls(l)**2*2./(2.d0*l+1.d0)) !!only compute the diagonal
  enddo

  call fp%open("clsout/clsout_"//trim(mask_spot)//"_"//COOP_STR_OF(theta1)//"to"//COOP_STR_OF(theta2)//".dat", "w")
  write(fp%unit, "(A8, 3A16)") "# ell ",  "  model C_l  ", " ave Cl  ", "  delta C_l "
  do l=2, lmax/2
     write(fp%unit, "(I8, 3E16.7)") l, Cls(l)*l*(l+1.d0)/coop_2pi, Cls_ave(l)*l*(l+1.d0)/coop_2pi,  delta_Cls(l)*l*(l+1.d0)/coop_2pi
  enddo
  call fp%close()
  call coop_MPI_finalize()  
end program test
