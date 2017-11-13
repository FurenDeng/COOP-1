module qf
  
contains

  function qwp1f(a, lambda) result(f)
    !! calculate \int_0^a \sqrt{\frac{x^7}{1+x^3+\lambda x}} dx
    !! assuming that |lambda| is small and a <~ 1.5
    real*8,parameter::cutoff = 0.75d0    
    real*8 a, lambda, f, x, lam2, lam3, lam4
    lam2 = lambda**2
    lam3 = lam2*lambda
    lam4 = lam2*lam2
    if(a.lt. cutoff)then
       f = approx1(a)
    else
       f = approx1(cutoff) + approx2(a-1.d0) - approx2(cutoff-1.d0)
    endif

  contains
    
    function approx1(y)
      real*8 y, approx1
      approx1 = y**(9.d0/2.d0) * ( &
           2.d0/9.d0 + y*( &
           -lambda/11.d0 + y*( &
           3.d0/52.d0*lam2 + y*( &
           -1.d0/15.d0-lam3/24.d0 + y*( &
           3.d0/34.d0*lambda + 35.d0/1088.d0*lam4 + y*( &
           -15.d0/152.d0*lam2 + y*( &
           1.d0/28.d0+5.d0/48.d0*lam3 + y*( &
           -15.d0/184.d0*lambda-315.d0/2944.d0*lam4 + y*( &
           21.d0/160.d0*lam2 + y*(&
           -5.d0/216.d0 - 35.d0/192.d0*lam3 + y*( &
           35.d0/464.d0*lambda + 3465.d0/14848.d0*lam4)) )) ))))))) 
    end function approx1

    function approx2(y)
      real*8 y, approx2, t
      t = y/(2.d0+lambda)
      approx2 = y/sqrt(2.d0+lambda)*( &
           1.d0 + t* ( &
           (11.d0+6.d0*lambda)/4.d0  +  t*( &
           (24.d0*lam2 + 76.d0*lambda+59.d0)/24.d0 + t*( &
           (-49.d0-22.d0*lambda + 32.d0*lam2+16.d0*lam3)/64.d0 + t*( &
           19.d0-936.d0*lambda-1104.d0*lam2-320.d0*lam3)/640.d0 ))))
    end function approx2


    
  end function qwp1f
end module qf
  
program Test
  use qf
  implicit none
  real*8 a, lambda
  a=1.4d0
  lambda = 0.d0
  write(*,"(2F10.5)") qwp1f(a, lambda)*3.d0/a**3, sqrt((1+a**3)/a**3)-log(a**1.5+sqrt(1.d0+a**3))/a**3


end program Test



