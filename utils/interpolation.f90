module coop_interpolation_mod
  use coop_wrapper_typedef
  use coop_matrix_mod
  implicit none

  private

  integer,parameter::dl = kind(1.d0)
  integer,parameter::sp = kind(1.)


  public::coop_bilinear_interp, coop_bicubic_interp

  Interface coop_bilinear_interp
     module procedure bilinear_interp_s, bilinear_interp_v, bilinear_interp_s_sp, bilinear_interp_v_sp
  end Interface coop_bilinear_interp


  Interface coop_bicubic_interp
     module procedure bicubic_interp_s, bicubic_interp_v, bicubic_interp_s_sp, bicubic_interp_v_sp, bicubic_interp_v_sp2sp
  end Interface coop_bicubic_interp



contains

  subroutine bilinear_interp_s(f, rx, ry, z)
    real(dl) ,intent(IN)::f(2,2)
    real(dl) rx, ry
    real(dl) z
    z = (f(1,1)*(1.d0-rx)+f(2,1)*rx)*(1.d0-ry) + (f(1,2)*(1.d0-rx)+f(2,2)*rx)*ry
  end subroutine bilinear_interp_s


  subroutine bilinear_interp_v(n, f, rx, ry, z)
    integer n
    real(dl) ,intent(IN)::f(n, 2,2)
    real(dl) rx, ry
    real(dl) z(n)
    z = (f(:,1,1)*(1.d0-rx)+f(:,2,1)*rx)*(1.d0-ry) + (f(:,1,2)*(1.d0-rx)+f(:,2,2)*rx)*ry
  end subroutine bilinear_interp_v


  subroutine bicubic_interp_s(f, rx, ry, z)
    real(dl),intent(IN)::f(4,4)
    real(dl) rx, ry
    REAL(dl) :: x(16)
    REAL(dl), INTENT(OUT) :: z
    real(dl), DIMENSION(16,16),parameter:: wt = reshape( (/ &
         0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0, &
         0,-2,0,0,0,0,0,0,0,2,0,0,0,0,0,0, &
         0,4,0,0,0,-10,0,0,0,8,0,0,0,-2,0,0, &
         0,-2,0,0,0,6,0,0,0,-6,0,0,0,2,0,0, &
         0,0,0,0,-2,0,2,0,0,0,0,0,0,0,0,0, &
         1,0,-1,0,0,0,0,0,-1,0,1,0,0,0,0,0, &
         -2,0,2,0,5,0,-5,0,-4,0,4,0,1,0,-1,0, &
         1,0,-1,0,-3,0,3,0,3,0,-3,0,-1,0,1,0, &
         0,0,0,0,4,-10,8,-2,0,0,0,0,0,0,0,0, &
         -2,5,-4,1,0,0,0,0,2,-5,4,-1,0,0,0,0, &
         4,-10,8,-2,-10,25,-20,5,8,-20,16,-4,-2,5,-4,1, &
         -2,5,-4,1,6,-15,12,-3,-6,15,-12,3,2,-5,4,-1, &
         0,0,0,0,-2,6,-6,2,0,0,0,0,0,0,0,0, &
         1,-3,3,-1,0,0,0,0,-1,3,-3,1,0,0,0,0, &
         -2,6,-6,2,5,-15,15,-5,-4,12,-12,4,1,-3,3,-1, &
         1,-3,3,-1,-3,9,-9,3,3,-9,9,-3,-1,3,-3,1 &
         /), (/ 16,16 /)) / 4.d0
    x = reshape(matmul(reshape(f, (/1, 16/) ), wt), (/ 16 /))
    z=((x(16 )*ry+x(15 ))*ry+x(14 ))*ry+x(13 )
    z=rx*z+((x(12 )*ry+x(11 ))*ry+x(10 ))*ry+x(9 )
    z=rx*z+((x(8 )*ry+x(7 ))*ry+x(6 ))*ry+x(5 )
    z=rx*z+((x(4 )*ry+x(3 ))*ry+x(2 ))*ry+x(1 )
  End subroutine bicubic_interp_s

  Subroutine bicubic_interp_v(n, f, rx, ry, z)
    integer n
    real(dl),intent(IN)::f(n,4,4)
    real(dl) rx, ry
    REAL(dl) :: x(n, 16)
    REAL(dl), INTENT(OUT) :: z(n)
    real(dl), DIMENSION(16,16),parameter:: wt = reshape( (/ &
         0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0, &
         0,-2,0,0,0,0,0,0,0,2,0,0,0,0,0,0, &
         0,4,0,0,0,-10,0,0,0,8,0,0,0,-2,0,0, &
         0,-2,0,0,0,6,0,0,0,-6,0,0,0,2,0,0, &
         0,0,0,0,-2,0,2,0,0,0,0,0,0,0,0,0, &
         1,0,-1,0,0,0,0,0,-1,0,1,0,0,0,0,0, &
         -2,0,2,0,5,0,-5,0,-4,0,4,0,1,0,-1,0, &
         1,0,-1,0,-3,0,3,0,3,0,-3,0,-1,0,1,0, &
         0,0,0,0,4,-10,8,-2,0,0,0,0,0,0,0,0, &
         -2,5,-4,1,0,0,0,0,2,-5,4,-1,0,0,0,0, &
         4,-10,8,-2,-10,25,-20,5,8,-20,16,-4,-2,5,-4,1, &
         -2,5,-4,1,6,-15,12,-3,-6,15,-12,3,2,-5,4,-1, &
         0,0,0,0,-2,6,-6,2,0,0,0,0,0,0,0,0, &
         1,-3,3,-1,0,0,0,0,-1,3,-3,1,0,0,0,0, &
         -2,6,-6,2,5,-15,15,-5,-4,12,-12,4,1,-3,3,-1, &
         1,-3,3,-1,-3,9,-9,3,3,-9,9,-3,-1,3,-3,1 &
         /), (/ 16,16 /)) / 4.d0
    x = matmul(reshape(f, (/ n , 16 /)), wt) 
    z=((x(:,16)*ry+x(:,15))*ry+x(:,14))*ry+x(:,13)
    z=rx*z+((x(:,12)*ry+x(:,11))*ry+x(:,10))*ry+x(:,9)
    z=rx*z+((x(:,8)*ry+x(:,7))*ry+x(:,6))*ry+x(:,5)
    z=rx*z+((x(:,4)*ry+x(:,3))*ry+x(:,2))*ry+x(:,1)
  End subroutine bicubic_interp_v

  subroutine bilinear_interp_s_sp(f, rx, ry, z)
    real(sp) ,intent(IN)::f(2,2)
    real(dl) rx, ry
    real(dl) z
    z = (f(1,1)*(1.d0-rx)+f(2,1)*rx)*(1.d0-ry) + (f(1,2)*(1.d0-rx)+f(2,2)*rx)*ry
  end subroutine bilinear_interp_s_sp


  subroutine bilinear_interp_v_sp(n, f, rx, ry, z)
    integer n
    real(sp) ,intent(IN)::f(n, 2,2)
    real(dl) rx, ry
    real(dl) z(n)
    z = (f(:,1,1)*(1.d0-rx)+f(:,2,1)*rx)*(1.d0-ry) + (f(:,1,2)*(1.d0-rx)+f(:,2,2)*rx)*ry
  end subroutine bilinear_interp_v_sp


  subroutine bicubic_interp_s_sp(f, rx, ry, z)
    real(sp),intent(IN)::f(4,4)
    real(dl) rx, ry
    real(dl) :: x(16)
    real(dl), INTENT(OUT) :: z
    real(dl), DIMENSION(16,16),parameter:: wt = reshape( (/ &
         0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0, &
         0,-2,0,0,0,0,0,0,0,2,0,0,0,0,0,0, &
         0,4,0,0,0,-10,0,0,0,8,0,0,0,-2,0,0, &
         0,-2,0,0,0,6,0,0,0,-6,0,0,0,2,0,0, &
         0,0,0,0,-2,0,2,0,0,0,0,0,0,0,0,0, &
         1,0,-1,0,0,0,0,0,-1,0,1,0,0,0,0,0, &
         -2,0,2,0,5,0,-5,0,-4,0,4,0,1,0,-1,0, &
         1,0,-1,0,-3,0,3,0,3,0,-3,0,-1,0,1,0, &
         0,0,0,0,4,-10,8,-2,0,0,0,0,0,0,0,0, &
         -2,5,-4,1,0,0,0,0,2,-5,4,-1,0,0,0,0, &
         4,-10,8,-2,-10,25,-20,5,8,-20,16,-4,-2,5,-4,1, &
         -2,5,-4,1,6,-15,12,-3,-6,15,-12,3,2,-5,4,-1, &
         0,0,0,0,-2,6,-6,2,0,0,0,0,0,0,0,0, &
         1,-3,3,-1,0,0,0,0,-1,3,-3,1,0,0,0,0, &
         -2,6,-6,2,5,-15,15,-5,-4,12,-12,4,1,-3,3,-1, &
         1,-3,3,-1,-3,9,-9,3,3,-9,9,-3,-1,3,-3,1 &
         /), (/ 16,16 /)) / 4.d0
    x = reshape(matmul(reshape(f, (/1, 16/) ), wt), (/ 16 /))
    z=((x(16 )*ry+x(15 ))*ry+x(14 ))*ry+x(13 )
    z=rx*z+((x(12 )*ry+x(11 ))*ry+x(10 ))*ry+x(9 )
    z=rx*z+((x(8 )*ry+x(7 ))*ry+x(6 ))*ry+x(5 )
    z=rx*z+((x(4 )*ry+x(3 ))*ry+x(2 ))*ry+x(1 )
  End subroutine bicubic_interp_s_sp

  Subroutine bicubic_interp_v_sp(n, f, rx, ry, z)
    integer n
    real(sp),intent(IN)::f(n,4,4)
    real(dl) rx, ry
    real(dl) :: x(n, 16)
    real(dl), INTENT(OUT) :: z(n)
    real(dl), DIMENSION(16,16),parameter:: wt = reshape( (/ &
         0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0, &
         0,-2,0,0,0,0,0,0,0,2,0,0,0,0,0,0, &
         0,4,0,0,0,-10,0,0,0,8,0,0,0,-2,0,0, &
         0,-2,0,0,0,6,0,0,0,-6,0,0,0,2,0,0, &
         0,0,0,0,-2,0,2,0,0,0,0,0,0,0,0,0, &
         1,0,-1,0,0,0,0,0,-1,0,1,0,0,0,0,0, &
         -2,0,2,0,5,0,-5,0,-4,0,4,0,1,0,-1,0, &
         1,0,-1,0,-3,0,3,0,3,0,-3,0,-1,0,1,0, &
         0,0,0,0,4,-10,8,-2,0,0,0,0,0,0,0,0, &
         -2,5,-4,1,0,0,0,0,2,-5,4,-1,0,0,0,0, &
         4,-10,8,-2,-10,25,-20,5,8,-20,16,-4,-2,5,-4,1, &
         -2,5,-4,1,6,-15,12,-3,-6,15,-12,3,2,-5,4,-1, &
         0,0,0,0,-2,6,-6,2,0,0,0,0,0,0,0,0, &
         1,-3,3,-1,0,0,0,0,-1,3,-3,1,0,0,0,0, &
         -2,6,-6,2,5,-15,15,-5,-4,12,-12,4,1,-3,3,-1, &
         1,-3,3,-1,-3,9,-9,3,3,-9,9,-3,-1,3,-3,1 &
         /), (/ 16,16 /)) / 4.d0
    x = matmul(reshape(f, (/ n , 16 /)), wt) 
    z=((x(:,16)*ry+x(:,15))*ry+x(:,14))*ry+x(:,13)
    z=rx*z+((x(:,12)*ry+x(:,11))*ry+x(:,10))*ry+x(:,9)
    z=rx*z+((x(:,8)*ry+x(:,7))*ry+x(:,6))*ry+x(:,5)
    z=rx*z+((x(:,4)*ry+x(:,3))*ry+x(:,2))*ry+x(:,1)
  End subroutine bicubic_interp_v_sp


  Subroutine bicubic_interp_v_sp2sp(n, f, rx, ry, z)
    integer n
    real(sp),intent(IN)::f(n,4,4)
    real(dl) rx, ry
    real(dl) :: x(n, 16)
    real(sp), INTENT(OUT) :: z(n)
    real(dl), DIMENSION(16,16),parameter:: wt = reshape( (/ &
         0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0, &
         0,-2,0,0,0,0,0,0,0,2,0,0,0,0,0,0, &
         0,4,0,0,0,-10,0,0,0,8,0,0,0,-2,0,0, &
         0,-2,0,0,0,6,0,0,0,-6,0,0,0,2,0,0, &
         0,0,0,0,-2,0,2,0,0,0,0,0,0,0,0,0, &
         1,0,-1,0,0,0,0,0,-1,0,1,0,0,0,0,0, &
         -2,0,2,0,5,0,-5,0,-4,0,4,0,1,0,-1,0, &
         1,0,-1,0,-3,0,3,0,3,0,-3,0,-1,0,1,0, &
         0,0,0,0,4,-10,8,-2,0,0,0,0,0,0,0,0, &
         -2,5,-4,1,0,0,0,0,2,-5,4,-1,0,0,0,0, &
         4,-10,8,-2,-10,25,-20,5,8,-20,16,-4,-2,5,-4,1, &
         -2,5,-4,1,6,-15,12,-3,-6,15,-12,3,2,-5,4,-1, &
         0,0,0,0,-2,6,-6,2,0,0,0,0,0,0,0,0, &
         1,-3,3,-1,0,0,0,0,-1,3,-3,1,0,0,0,0, &
         -2,6,-6,2,5,-15,15,-5,-4,12,-12,4,1,-3,3,-1, &
         1,-3,3,-1,-3,9,-9,3,3,-9,9,-3,-1,3,-3,1 &
         /), (/ 16,16 /)) / 4.d0
    x = matmul(reshape(f, (/ n , 16 /)), wt) 
    z=real(((x(:,16)*ry+x(:,15))*ry+x(:,14))*ry+x(:,13))
    z=real(rx*z+((x(:,12)*ry+x(:,11))*ry+x(:,10))*ry+x(:,9))
    z=real(rx*z+((x(:,8)*ry+x(:,7))*ry+x(:,6))*ry+x(:,5))
    z=real(rx*z+((x(:,4)*ry+x(:,3))*ry+x(:,2))*ry+x(:,1))
  End subroutine bicubic_interp_v_sp2sp

end module coop_interpolation_mod
