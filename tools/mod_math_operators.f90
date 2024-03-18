!> This module contains multi-purposes mathematical operators
module mod_math_operators

  implicit none
  private
  public cross_product
  !> define operators for grad_R, grad_Z, lap_RZ, lap_RR, lap_ZZ, etc.
  public op_grad_R, op_grad_Z, op_lap_RZ, op_lap_RR, op_lap_ZZ, op_xjac, op_xjac_R, op_xjac_Z

contains
  
!> The vector cross product
pure function cross_product(a, b)
  real*8, dimension(3) :: cross_product
  real*8, dimension(3), intent(in) :: a, b

  cross_product(1) = a(2) * b(3) - a(3) * b(2)
  cross_product(2) = a(3) * b(1) - a(1) * b(3)
  cross_product(3) = a(1) * b(2) - a(2) * b(1)

end function cross_product

pure function op_grad_R(eq_s, eq_t, R_s, R_t, Z_s, Z_t, xjac)
real*8, intent(in)   :: eq_s, eq_t, R_s, Z_s, Z_t, R_t, xjac
real*8               :: op_grad_R
op_grad_R = (  Z_t * eq_s - Z_s * eq_t) / xjac
end function op_grad_R

pure function op_grad_Z(eq_s, eq_t, R_s, R_t, Z_s, Z_t, xjac)
real*8, intent(in)   :: eq_s, eq_t, R_s, Z_s, Z_t, R_t, xjac
real*8               :: op_grad_Z
op_grad_Z = (- R_t * eq_s + R_s * eq_t) / xjac
end function op_grad_Z

pure function op_lap_RZ(eq_s, eq_t, eq_st, eq_ss, eq_tt,              &
                        R_s, R_t, R_st, R_ss, R_tt,                   &
                        Z_s, Z_t, Z_st, Z_ss, Z_tt,                   &
                        xjac, xjac_R, xjac_Z)
real*8, intent(in)   :: eq_s, eq_t, eq_st, eq_ss, eq_tt,              &
                        R_s, R_t, R_st, R_ss, R_tt,                   &
                        Z_s, Z_t, Z_st, Z_ss, Z_tt,                   &
                        xjac, xjac_R, xjac_Z
real*8               :: op_lap_RZ
op_lap_RZ = (- eq_ss * Z_t*R_t - eq_tt * R_s*Z_s                      &
          + eq_st * (Z_s*R_t  + Z_t*R_s  )                            &
          - eq_s  * (R_st*Z_t - R_tt*Z_s )                            &
          - eq_t * (R_st*Z_s  - R_ss*Z_t ) ) / xjac**2                &
          - xjac_R * (- eq_s * R_t + eq_t * R_s )   / xjac**2  
end function op_lap_RZ

pure function op_lap_RR(eq_s, eq_t, eq_st, eq_ss, eq_tt,              &
                        R_s, R_t, R_st, R_ss, R_tt,                   &
                        Z_s, Z_t, Z_st, Z_ss, Z_tt,                   &
                        xjac, xjac_R, xjac_Z)
real*8, intent(in)   :: eq_s, eq_t, eq_st, eq_ss, eq_tt,              &
                        R_s, R_t, R_st, R_ss, R_tt,                   &
                        Z_s, Z_t, Z_st, Z_ss, Z_tt,                   &
                        xjac, xjac_R, xjac_Z
real*8               :: op_lap_RR
op_lap_RR = (eq_ss * Z_t**2 - 2.d0*eq_st * Z_s*Z_t + eq_tt * Z_s**2   &
          + eq_s * (Z_st*Z_t - Z_tt*Z_s )                             &
          + eq_t * (Z_st*Z_s - Z_ss*Z_t ) ) / xjac**2                 &
          - xjac_R * (eq_s* Z_t - eq_t * Z_s)  / xjac**2
end function op_lap_RR

pure function op_lap_ZZ(eq_s, eq_t, eq_st, eq_ss, eq_tt,              &
                        R_s, R_t, R_st, R_ss, R_tt,                   &
                        Z_s, Z_t, Z_st, Z_ss, Z_tt,                   &
                        xjac, xjac_R, xjac_Z)
real*8, intent(in)   :: eq_s, eq_t, eq_st, eq_ss, eq_tt,              &
                        R_s, R_t, R_st, R_ss, R_tt,                   &
                        Z_s, Z_t, Z_st, Z_ss, Z_tt,                   &
                        xjac, xjac_R, xjac_Z
real*8               :: op_lap_ZZ
op_lap_ZZ = (eq_ss * R_t**2 - 2.d0*eq_st * R_s*R_t + eq_tt * R_s**2   &
          + eq_s * (R_st*R_t - R_tt*R_s )                             &
          + eq_t * (R_st*R_s - R_ss*R_t ) ) / xjac**2                 &
          - xjac_Z * (- eq_s * R_t + eq_t * R_s )  / xjac**2
end function op_lap_ZZ

pure function op_xjac(R_s, R_t, Z_s, Z_t)
real*8, intent(in)   :: R_s, R_t, Z_s, Z_t
real*8               :: op_xjac
op_xjac = R_s * Z_t - R_t * Z_s
end function op_xjac

pure function op_xjac_R(R_s, R_t, R_st, R_ss, R_tt,                   &
                        Z_s, Z_t, Z_st, Z_ss, Z_tt,                   &
                        xjac)
real*8, intent(in)   :: R_s, R_t, R_st, R_ss, R_tt,                   &
                        Z_s, Z_t, Z_st, Z_ss, Z_tt,                   &
                        xjac
real*8               :: op_xjac_R
op_xjac_R = (R_ss*Z_t**2 - Z_ss*R_t*Z_t - 2.d0*R_st*Z_s*Z_t           &
                + Z_st*(R_s*Z_t + R_t*Z_s)                            &
                + R_tt*Z_s**2 - Z_tt*R_s*Z_s) / xjac
end function op_xjac_R

pure function op_xjac_Z(R_s, R_t, R_st, R_ss, R_tt,                   &
                        Z_s, Z_t, Z_st, Z_ss, Z_tt,                   &
                        xjac)
real*8, intent(in)   :: R_s, R_t, R_st, R_ss, R_tt,                   &
                        Z_s, Z_t, Z_st, Z_ss, Z_tt,                   &
                        xjac
real*8               :: op_xjac_Z
op_xjac_Z = (Z_tt*R_s**2 - R_tt*Z_s*R_s - 2.d0*Z_st*R_t*R_s           &
                + R_st*(Z_t*R_s + Z_s*R_t)                            &
                + Z_ss*R_t**2 - R_ss*Z_t*R_t) / xjac
end function op_xjac_Z

end module mod_math_operators
