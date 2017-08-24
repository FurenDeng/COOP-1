program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_INT,parameter::n = 80
  COOP_REAL::f(-n:n, -n:n)
  COOP_INT,parameter::mmax = 4
  COOP_REAL::cr(1:n, 0:mmax), sr(1:n, 0:mmax), Ck(1:n, 0:mmax), Sk(1:n, 0:mmax)
  COOP_INT::i, j
  do i=-n, n
     do j=-n, n
        f(i, j) = map(i, j)
     end do
  enddo
  call coop_2D_radial_decompose(n, f, mmax, Cr, Sr)
  do i=1, n
     print*, cr(i, 0), c0(dble(i))
  enddo

contains

  function c0(r)
    COOP_REAL::r, c0
    c0 = 1.d0+sqrt(r)
  end function c0

  function c1(r)
    COOP_REAL::r, c1
    c1 = r
  end function c1
  
  function s1(r)
    COOP_REAL::r, s1
    s1 = r**2/(r**3+1)
  end function s1
  
  function c2(r)
    COOP_REAL::r, c2
    c2 = r**3/(r**4+1)
  end function c2
  
  function s2(r)
    COOP_REAL::r, s2
    s2 = r**4*exp(-r/5.d0)
  end function s2

  function c3(r)
    COOP_REAL::r, c3
    c3 = 0.d0
  end function c3
  
  function s3(r)
    COOP_REAL::r, s3
    s3 = sqrt(r)/(1.d0+r)
  end function s3

  function c4(r)
    COOP_REAL::r, c4
    c4 = r/(r**2+1.d0)
  end function c4
  
  function s4(r)
    COOP_REAL::r, s4
    s4 = r**2/(r**4+1.d0)
  end function s4
  
  

  function map(ix, iy)
    COOP_INT::ix, iy
    COOP_REAL::r, theta
    COOP_REAL::map
    r = sqrt(dble(ix**2+iy**2))
    theta = COOP_POLAR_ANGLE(dble(ix), dble(iy))
    map = c0(r) + c1(r)*cos(theta) + s1(r)*sin(theta) + c2(r)*cos(2.d0*theta) + s2(r)*sin(2.d0*theta) + c3(r)*cos(3.d0*theta) + s3(r)*sin(3.d0*theta) + c4(r)*cos(4.d0*theta) + s4(r) * sin(4.d0*theta)
  end function map
  
end program Test  
