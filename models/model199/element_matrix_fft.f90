module elm_fft
use parameters

  real*8 :: ELM_p(n_plane,n_vertex_max*n_var*(n_order+1),n_vertex_max*n_var*(n_order+1)*n_tor)
  real*8 :: ELM_n(n_plane,n_vertex_max*n_var*(n_order+1),n_vertex_max*n_var*(n_order+1)*n_tor)
  real*8 :: ELM_k(n_plane,n_vertex_max*n_var*(n_order+1),n_vertex_max*n_var*(n_order+1)*n_tor)
  real*8 :: ELM_kn(n_plane,n_vertex_max*n_var*(n_order+1),n_vertex_max*n_var*(n_order+1)*n_tor)
  real*8 :: RHS_p(n_plane,n_vertex_max*n_var*(n_order+1))
  real*8 :: RHS_k(n_plane,n_vertex_max*n_var*(n_order+1))

!$OMP THREADPRIVATE(ELM_p, ELM_n, ELM_k, ELM_kn, RHS_p, RHS_k)
end module

subroutine element_matrix_fft(element, nodes, xpoint2, psi_axis, psi_bnd, Z_xpoint, ELM,RHS)
!---------------------------------------------------------------
! calculates the matrix contribution of one element
!
! TO BE CHECKED (see element_matrix)
!
!---------------------------------------------------------------
use parameters
use data_structure
use gauss
use basis_at_gaussian
use phys_module
use elm_fft
implicit none
 
type (type_element)   :: element
type (type_node)      :: nodes(n_vertex_max)

real*8     :: ELM(n_vertex_max*n_var*(n_order+1)*n_tor,n_vertex_max*n_var*(n_order+1)*n_tor)
real*8     :: RHS(n_vertex_max*n_var*(n_order+1)*n_tor)

real*8     :: x_g(n_gauss,n_gauss),        x_s(n_gauss,n_gauss),        x_t(n_gauss,n_gauss)
real*8     :: x_ss(n_gauss,n_gauss),       x_st(n_gauss,n_gauss),       x_tt(n_gauss,n_gauss)
real*8     :: y_g(n_gauss,n_gauss),        y_s(n_gauss,n_gauss),        y_t(n_gauss,n_gauss)
real*8     :: y_ss(n_gauss,n_gauss),       y_st(n_gauss,n_gauss),       y_tt(n_gauss,n_gauss)

real*8     :: eq_g(n_plane,n_var,n_gauss,n_gauss), eq_s(n_plane,n_var,n_gauss,n_gauss)
real*8     :: eq_t(n_plane,n_var,n_gauss,n_gauss), eq_p(n_plane,n_var,n_gauss,n_gauss)
real*8     :: eq_ss(n_plane,n_var,n_gauss,n_gauss),eq_st(n_plane,n_var,n_gauss,n_gauss),eq_tt(n_plane,n_var,n_gauss,n_gauss)
real*8     :: delta_g(n_plane,n_var,n_gauss,n_gauss),delta_s(n_plane,n_var,n_gauss,n_gauss),delta_t(n_plane,n_var,n_gauss,n_gauss)

integer    :: i, j, ms, mt, mp, k, l, index_ij, index_kl, index, index_k, index_m, m, ik
integer    :: in, im, ij1, ij2, ij3, ij4, ij5, ij6, kl1, kl2, kl3, kl4, kl5, kl6
real*8     :: wst, xjac, xjac_s, xjac_t, BigR, r2, PI, phi, eps_cyl
real*8     :: current_source(n_gauss,n_gauss),particle_source(n_gauss,n_gauss),heat_source(n_gauss,n_gauss)
real*8     :: psi_axis, psi_bnd, Z_xpoint, dj_dpsi, dj_dz
real*8     :: Bgrad_rho_star,     Bgrad_rho,     Bgrad_T_star,  Bgrad_T, BB2
real*8     :: Bgrad_rho_star_psi, Bgrad_rho_psi, Bgrad_rho_rho, Bgrad_T_star_psi, Bgrad_T_psi, Bgrad_T_T, BB2_psi
real*8     :: Bgrad_rho_rho_n, Bgrad_T_T_n, Bgrad_rho_k_star, Bgrad_T_k_star
real*8     :: D_prof, ZK_prof, psi_norm
real*8     :: rhs_ij_1,   rhs_ij_2,   rhs_ij_3,   rhs_ij_4,   rhs_ij_5,   rhs_ij_6
real*8     :: rhs_ij_5_k, rhs_ij_6_k
real*8     :: rhs_stab_1, rhs_stab_2, rhs_stab_3, rhs_stab_4, rhs_stab_5, rhs_stab_6

real*8     :: v, v_x, v_y, v_s, v_t, v_p, v_ss, v_st, v_tt, v_xx, v_yy, v_xs, v_ys, v_xt, v_yt
real*8     :: ps0, ps0_x, ps0_y, ps0_p,ps0_s,ps0_t,  zj0, zj0_x, zj0_y, zj0_p, zj0_s, zj0_t
real*8     :: u0, u0_x, u0_y, u0_p, u0_s, u0_t,  w0, w0_x, w0_y, w0_p, w0_s, w0_t
real*8     :: r0, r0_x, r0_y, r0_p, r0_s, r0_t,  r0_hat, r0_x_hat, r0_y_hat, T0, T0_x, T0_y, T0_p, T0_s, T0_t
real*8     :: psi, psi_x, psi_y, psi_p, psi_s, psi_t, psi_ss, psi_st, psi_tt, psi_xs, psi_ys, psi_xt, psi_yt, psi_xx, psi_yy
real*8     :: zj, zj_x, zj_y, zj_p, zj_s, zj_t
real*8     :: u, u_x, u_y, u_p, u_s, u_t, w, w_x, w_y, w_p, w_s, w_t
real*8     :: rho, rho_x, rho_y, rho_s, rho_t, rho_p, rho_hat, rho_x_hat, rho_y_hat, T, T_x, T_y, T_s, T_t, T_p
real*8     :: w0_xs, w0_xt, w0_ys, w0_yt, w0_xx, w0_yy, P0, P0_s, P0_t, P0_x, P0_y
real*8     :: BigR_x, vv2, eta_T, visco_T, deta_dT, d2eta_d2T, dvisco_dT
real*8     :: amat_11, amat_12, amat_21, amat_22, amat_23, amat_24, amat_25, amat_26, amat_33, amat_31, amat_44, amat_42
real*8     :: amat_51, amat_52, amat_55, amat_61, amat_62, amat_66, amat_16, amat_13
real*8     :: amat_12_n, amat_23_n, amat_51_k, amat_55_kn, amat_55_k, amat_55_n, amat_61_k, amat_66_kn, amat_66_k, amat_66_n
real*8     :: amat_stab_11, amat_stab_12, amat_stab_13, amat_stab_14 ,amat_stab_21,amat_stab_22, amat_stab_23, amat_stab_24
real*8     :: amat_stab_31, amat_stab_32, amat_stab_33, amat_stab_34 ,amat_stab_41,amat_stab_42, amat_stab_43, amat_stab_44
real*8     :: theta, zeta, delta_u_x, delta_u_y
logical    :: xpoint2

integer*8  :: plan
real*8     :: in_fft(1:n_plane)
complex*16 :: out_fft(1:n_plane)
INTEGER    :: FFTW_FORWARD,  FFTW_BACKWARD, FFTW_ESTIMATE
PARAMETER (FFTW_FORWARD=-1,FFTW_BACKWARD=+1, FFTW_ESTIMATE=64)


ELM_p = 0.d0
ELM_n = 0.d0
ELM_k = 0.d0
ELM_kn = 0.d0
RHS_p = 0.d0
RHS_k = 0.d0
ELM   = 0.d0
RHS   = 0.d0

PI    = 2.d0*asin(1.d0)
GAMMA = 5.d0 / 3.d0

theta = 0.5d0  ; zeta = 0.0d0      ! Crank-Nicholson scheme
!theta = 1.0d0  ; zeta = 0.0d0      ! Euler scheme 
!theta = 1.0d0   ; zeta = 0.5d0      ! BDF2 (Gears) scheme

current_source  = 0.d0
particle_source = 0.d0
heat_source     = 0.d0

!---------------------------------------------------- value of (x,y) and derivatives on Gaussian points
x_g  = 0.d0; x_s  = 0.d0; x_t  = 0.d0; x_st  = 0.d0; x_ss  = 0.d0; x_tt  = 0.d0;
y_g  = 0.d0; y_s  = 0.d0; y_t  = 0.d0; y_st  = 0.d0; y_ss  = 0.d0; y_tt  = 0.d0;
eq_g = 0.d0; eq_s = 0.d0; eq_t = 0.d0; eq_st = 0.d0; eq_ss = 0.d0; eq_tt = 0.d0; eq_p = 0.d0;

delta_g = 0.d0; delta_s = 0.d0; delta_t = 0.d0

do i=1,n_vertex_max
 do j=1,n_order+1

   do ms=1, n_gauss
     do mt=1, n_gauss

       x_g(ms,mt)  = x_g(ms,mt)  + nodes(i)%x(j,1) * element%size(i,j) * H(i,j,ms,mt)
       x_s(ms,mt)  = x_s(ms,mt)  + nodes(i)%x(j,1) * element%size(i,j) * H_s(i,j,ms,mt)
       x_t(ms,mt)  = x_t(ms,mt)  + nodes(i)%x(j,1) * element%size(i,j) * H_t(i,j,ms,mt)

       y_g(ms,mt)  = y_g(ms,mt)  + nodes(i)%x(j,2) * element%size(i,j) * H(i,j,ms,mt)
       y_s(ms,mt)  = y_s(ms,mt)  + nodes(i)%x(j,2) * element%size(i,j) * H_s(i,j,ms,mt)
       y_t(ms,mt)  = y_t(ms,mt)  + nodes(i)%x(j,2) * element%size(i,j) * H_t(i,j,ms,mt)

       do mp=1,n_plane

         do k=1,n_var

           do in=1,n_tor

             eq_g(mp,k,ms,mt) = eq_g(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H(i,j,ms,mt)  * HZ(in,mp)

             eq_s(mp,k,ms,mt) = eq_s(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_s(i,j,ms,mt)* HZ(in,mp)

             eq_t(mp,k,ms,mt) = eq_t(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_t(i,j,ms,mt)* HZ(in,mp)

             eq_p(mp,k,ms,mt) = eq_p(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H(i,j,ms,mt)  * HZ_p(in,mp)

             delta_g(mp,k,ms,mt) = delta_g(mp,k,ms,mt) + nodes(i)%deltas(in,j,k) * element%size(i,j) * H(i,j,ms,mt)   * HZ(in,mp)
             delta_s(mp,k,ms,mt) = delta_s(mp,k,ms,mt) + nodes(i)%deltas(in,j,k) * element%size(i,j) * H_s(i,j,ms,mt) * HZ(in,mp)
             delta_t(mp,k,ms,mt) = delta_t(mp,k,ms,mt) + nodes(i)%deltas(in,j,k) * element%size(i,j) * H_t(i,j,ms,mt) * HZ(in,mp)

           enddo

         enddo

       enddo

       call current(xpoint2, x_g(ms,mt),y_g(ms,mt), Z_xpoint, eq_g(1,1,ms,mt),psi_axis,psi_bnd,current_source(ms,mt))
       call sources(xpoint2, y_g(ms,mt), Z_xpoint, eq_g(1,1,ms,mt),psi_axis,psi_bnd,particle_source(ms,mt),heat_source(ms,mt))

     enddo
   enddo
 enddo
enddo


!--------------------------------------------------- sum over the Gaussian integration points
do ms=1, n_gauss

 do mt=1, n_gauss

   wst = wgauss(ms)*wgauss(mt)

   xjac   = x_s(ms,mt)*y_t(ms,mt)  - x_t(ms,mt)*y_s(ms,mt)

   BigR    = x_g(ms,mt)
   BigR_x  = 1.d0

   eps_cyl = 1.d0          ! for cylinder geometry : epscyl = eps

   do mp = 1, n_plane

     ps0   = eq_g(mp,1,ms,mt)
     ps0_x = (   y_t(ms,mt) * eq_s(mp,1,ms,mt) - y_s(ms,mt) * eq_t(mp,1,ms,mt) ) / xjac
     ps0_y = ( - x_t(ms,mt) * eq_s(mp,1,ms,mt) + x_s(ms,mt) * eq_t(mp,1,ms,mt) ) / xjac
     ps0_p = eq_p(mp,1,ms,mt)
     ps0_s = eq_s(mp,1,ms,mt)
     ps0_t = eq_t(mp,1,ms,mt)

     u0    = eq_g(mp,2,ms,mt)
     u0_x  = (   y_t(ms,mt) * eq_s(mp,2,ms,mt) - y_s(ms,mt) * eq_t(mp,2,ms,mt) ) / xjac
     u0_y  = ( - x_t(ms,mt) * eq_s(mp,2,ms,mt) + x_s(ms,mt) * eq_t(mp,2,ms,mt) ) / xjac
     u0_p  = eq_p(mp,2,ms,mt)
     u0_s  = eq_s(mp,2,ms,mt)
     u0_t  = eq_t(mp,2,ms,mt)

     vv2   = BigR**2 *  ( u0_x * u0_x + u0_y *u0_y  )

     zj0   = eq_g(mp,3,ms,mt)
     zj0_x = (   y_t(ms,mt) * eq_s(mp,3,ms,mt) - y_s(ms,mt) * eq_t(mp,3,ms,mt) ) / xjac
     zj0_y = ( - x_t(ms,mt) * eq_s(mp,3,ms,mt) + x_s(ms,mt) * eq_t(mp,3,ms,mt) ) / xjac
     zj0_p = eq_p(mp,3,ms,mt)
     zj0_s = eq_s(mp,3,ms,mt)
     zj0_t = eq_t(mp,3,ms,mt)

     w0    = eq_g(mp,4,ms,mt)
     w0_x  = (   y_t(ms,mt) * eq_s(mp,4,ms,mt) - y_s(ms,mt) * eq_t(mp,4,ms,mt) ) / xjac
     w0_y  = ( - x_t(ms,mt) * eq_s(mp,4,ms,mt) + x_s(ms,mt) * eq_t(mp,4,ms,mt) ) / xjac
     w0_p  = eq_p(mp,4,ms,mt)
     w0_s  = eq_s(mp,4,ms,mt)
     w0_t  = eq_t(mp,4,ms,mt)

     r0    = abs(eq_g(mp,5,ms,mt))
     r0_x  = (   y_t(ms,mt) * eq_s(mp,5,ms,mt) - y_s(ms,mt) * eq_t(mp,5,ms,mt) ) / xjac
     r0_y  = ( - x_t(ms,mt) * eq_s(mp,5,ms,mt) + x_s(ms,mt) * eq_t(mp,5,ms,mt) ) / xjac
     r0_p  = eq_p(mp,5,ms,mt)
     r0_s  = eq_s(mp,5,ms,mt)
     r0_t  = eq_t(mp,5,ms,mt)

     r0_hat   = BigR**2 * r0
     r0_x_hat = 2.d0 * BigR * BigR_x  * r0 + BigR**2 * r0_x
     r0_y_hat = BigR**2 * r0_y

     T0    = abs(eq_g(mp,6,ms,mt))
     T0_x  = (   y_t(ms,mt) * eq_s(mp,6,ms,mt) - y_s(ms,mt) * eq_t(mp,6,ms,mt) ) / xjac
     T0_y  = ( - x_t(ms,mt) * eq_s(mp,6,ms,mt) + x_s(ms,mt) * eq_t(mp,6,ms,mt) ) / xjac
     T0_p  = eq_p(mp,6,ms,mt)
     T0_s  = eq_s(mp,6,ms,mt)
     T0_t  = eq_t(mp,6,ms,mt)

     P0    = r0 * T0
     P0_x  = r0_x * T0 + r0 * T0_x
     P0_y  = r0_y * T0 + r0 * T0_y
     P0_s  = r0_s * T0 + r0 * T0_s
     P0_t  = r0_t * T0 + r0 * T0_t

     delta_u_x = (   y_t(ms,mt) * delta_s(mp,2,ms,mt) - y_s(ms,mt) * delta_t(mp,2,ms,mt) ) / xjac
     delta_u_y = ( - x_t(ms,mt) * delta_s(mp,2,ms,mt) + x_s(ms,mt) * delta_t(mp,2,ms,mt) ) / xjac

     eta_T   = eta   * (T0/T_0)**(-1.5d0)                   ! temperature dependent resistivity
     visco_T = visco * (T0/T_0)**(-1.5d0)                   ! temperature dependent viscosity

     deta_dT   = - eta   * (1.5d0)  * T0**(-2.5d0) * T_0**(1.5d0)
     d2eta_d2T =   eta   * (3.75d0) * T0**(-3.5d0) * T_0**(1.5d0)
     dvisco_dT = - visco * (1.5d0)  * T0**(-2.5d0) * T_0**(1.5d0)


     psi_norm = (ps0 - psi_axis)/(psi_bnd - psi_axis)
     if (xpoint2) then
       if ((psi_norm .lt. 1.d0) .and. (y_g(ms,mt) .lt. Z_xpoint)) then
         psi_norm = 2.d0 - psi_norm
       endif
     endif
     D_prof  = D_perp(1)  * ((1.d0-D_perp(2))  + D_perp(2)  *(0.5d0 - 0.5d0*tanh((psi_norm-D_perp(5)) /D_perp(4))))
     ZK_prof = ZK_perp(1) * ((1.d0-ZK_perp(2)) + ZK_perp(2) *(0.5d0 - 0.5d0*tanh((psi_norm-ZK_perp(5))/ZK_perp(4))))

     do i=1,n_vertex_max

       do j=1,n_order+1

         index_ij = n_var*(n_order+1)*(i-1) + n_var * (j-1) + 1   ! index in the ELM matrix

         v   =  H(i,j,ms,mt) * element%size(i,j)
         v_x = (  y_t(ms,mt) * h_s(i,j,ms,mt) - y_s(ms,mt) * h_t(i,j,ms,mt) ) * element%size(i,j) / xjac
         v_y = (- x_t(ms,mt) * h_s(i,j,ms,mt) + x_s(ms,mt) * h_t(i,j,ms,mt) ) * element%size(i,j) / xjac

         v_s = h_s(i,j,ms,mt) * element%size(i,j)
         v_t = h_t(i,j,ms,mt) * element%size(i,j)
         v_p = H(i,j,ms,mt)   * element%size(i,j)

         Bgrad_rho_star   = ( v_x  * ps0_y - v_y  * ps0_x ) / BigR                         ! only the part without    d/d(phi)
         Bgrad_rho_k_star = ( F0 / BigR * v_p )           / BigR                           ! only the part containing d/d(phi)
         Bgrad_rho        = ( F0 / BigR * r0_p +  r0_x * ps0_y - r0_y * ps0_x ) / BigR
         Bgrad_T_star     = ( v_x  * ps0_y - v_y  * ps0_x ) / BigR                         ! only the part without    d/d(phi)
         Bgrad_T_k_star   = ( F0 / BigR * v_p           ) / BigR                         ! only the part containing d/d(phi)
         Bgrad_T          = ( F0 / BigR * T0_p +  T0_x * ps0_y - T0_y * ps0_x ) / BigR

         BB2 = (F0*F0 + ps0_x * ps0_x + ps0_y * ps0_y )/BigR**2

         rhs_ij_1 =   v * eta_T  * (zj0 - current_source(ms,mt))/ BigR    * xjac * tstep &
                      + v * (ps0_s * u0_t - ps0_t * u0_s)                        * tstep &
                      - v * eps_cyl * F0 / BigR  * u0_p                   * xjac * tstep &          ! F0 due to absence of normalisation
                      + eta_num * (v_x * zj0_x + v_y * zj0_y)             * xjac * tstep &
                      + zeta * v * delta_g(mp,1,ms,mt) / BigR             * xjac 

         rhs_ij_2 = - 0.5d0 * vv2 * (v_x * r0_y_hat - v_y * r0_x_hat)     * xjac * tstep &
                      - r0_hat * BigR**2 * w0 * (v_s * u0_t - v_t * u0_s)        * tstep &
                      + v * (ps0_s * zj0_t - ps0_t * zj0_s )                     * tstep &
                      - visco_T * BigR * (v_x * w0_x + v_y * w0_y)        * xjac * tstep &
                      - v * eps_cyl * F0 / BigR * zj0_p                   * xjac * tstep &         ! F0 due to absence of normalisation
                      + BigR**2 * (v_s * p0_t - v_t * p0_s)                      * tstep &
                      + visco_num * (v_s * w0_s + v_t * w0_t)                    * tstep &
                      - zeta * BigR * r0_hat * (v_x * delta_u_x + v_y * delta_u_y) * xjac  

         rhs_ij_3 = 0.d0 !- ( v_x * ps0_x  + v_y * ps0_y + v*zj0) / BigR * xjac * tstep
         rhs_ij_4 = 0.d0 !- ( v_x * u0_x   + v_y * u0_y  + v*w0)  * BigR * xjac * tstep

         rhs_ij_5   = v * BigR * particle_source(ms,mt)                                        * xjac * tstep &
                    + v * BigR**2 * ( r0_s * u0_t - r0_t * u0_s)                                      * tstep &
                    + v * 2.d0 * BigR * r0 * u0_y                                              * xjac * tstep &
                    - (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star * Bgrad_rho                 * xjac * tstep &
                    - D_prof * BigR  * (v_x*r0_x + v_y*r0_y                                  ) * xjac * tstep &
                    + zeta * v * delta_g(mp,5,ms,mt) * BigR                                    * xjac 

         rhs_ij_5_k =  - (D_par-D_prof) * BigR / BB2 * Bgrad_rho_k_star * Bgrad_rho            * xjac * tstep &
                       - D_prof * BigR  * (                  v_p*r0_p * eps_cyl**2 /BigR**2 )  * xjac * tstep

         rhs_ij_6 =   v * BigR * heat_source(ms,mt)                               * xjac * tstep &
                    + v * BigR**2 * ( T0_s * u0_t - T0_t * u0_s)                         * tstep &
                    + v * 2.d0* (GAMMA-1.d0) * BigR * T0 * u0_y                   * xjac * tstep &
                    - (ZK_par-ZK_prof) * BigR / BB2 * Bgrad_T_star * Bgrad_T      * xjac * tstep &
                    - ZK_prof * BigR * (v_x*T0_x + v_y*T0_y                     ) * xjac * tstep &
                    + zeta * v * delta_g(mp,6,ms,mt) * BigR                       * xjac 

         rhs_ij_6_k = - (ZK_par-ZK_prof) * BigR / BB2 * Bgrad_T_k_star * Bgrad_T * xjac * tstep &
                      - ZK_prof * BigR * (                + v_p*T0_p /BigR**2 )  * xjac * tstep

         ij1 = index_ij
         ij2 = index_ij + 1
         ij3 = index_ij + 2
         ij4 = index_ij + 3
         ij5 = index_ij + 4
         ij6 = index_ij + 5

         RHS_p(mp,ij1) = RHS_p(mp,ij1) + rhs_ij_1 * wst
         RHS_p(mp,ij2) = RHS_p(mp,ij2) + rhs_ij_2 * wst
         RHS_p(mp,ij3) = RHS_p(mp,ij3) + rhs_ij_3 * wst
         RHS_p(mp,ij4) = RHS_p(mp,ij4) + rhs_ij_4 * wst
         RHS_p(mp,ij5) = RHS_p(mp,ij5) + rhs_ij_5 * wst
         RHS_p(mp,ij6) = RHS_p(mp,ij6) + rhs_ij_6 * wst

         RHS_k(mp,ij5) = RHS_k(mp,ij5) + rhs_ij_5_k * wst
         RHS_k(mp,ij6) = RHS_k(mp,ij6) + rhs_ij_6_k * wst

         do k=1,n_vertex_max

           do l=1,n_order+1

             psi   = H(k,l,ms,mt) * element%size(k,l)

             psi_x = (   y_t(ms,mt) * h_s(k,l,ms,mt) - y_s(ms,mt) * h_t(k,l,ms,mt) ) / xjac * element%size(k,l)
             psi_y = ( - x_t(ms,mt) * h_s(k,l,ms,mt) + x_s(ms,mt) * h_t(k,l,ms,mt) ) / xjac * element%size(k,l)

             psi_p = H(k,l,ms,mt)   * element%size(k,l)
             psi_s = h_s(k,l,ms,mt) * element%size(k,l)
             psi_t = h_t(k,l,ms,mt) * element%size(k,l)

             u   = psi   ;    zj   = psi   ;    w   = psi    ;    rho   = psi    ;    T   = psi
             u_x = psi_x ;    zj_x = psi_x ;    w_x = psi_x  ;    rho_x = psi_x  ;    T_x = psi_x
             u_y = psi_y ;    zj_y = psi_y ;    w_y = psi_y  ;    rho_y = psi_y  ;    T_y = psi_y
             u_p = psi_p ;    zj_p = psi_p ;    w_p = psi_p  ;    rho_p = psi_p  ;    T_p = psi_p
             u_s = psi_s ;    zj_s = psi_s ;    w_s = psi_s  ;    rho_s = psi_s  ;    T_s = psi_s
             u_t = psi_t ;    zj_t = psi_t ;    w_t = psi_t  ;    rho_t = psi_t  ;    T_t = psi_t

             rho_hat   = BigR**2 * rho
             rho_x_hat = 2.d0 * BigR * BigR_x  * rho + BigR**2 * rho_x
             rho_y_hat = BigR**2 * rho_y

             index_kl = n_var*(n_order+1)*(k-1) + n_var * (l-1) + 1   ! index in the ELM matrix

!------------------------------------------------------------ equation 1
             amat_11 = v * psi / BigR * xjac * (1.d0+zeta)                                             &
                     - v * (psi_s * u0_t - psi_t * u0_s)                              * theta * tstep

             amat_12 = -  v * (ps0_s * u_t - ps0_t * u_s)                             * theta * tstep

             amat_12_n = +  eps_cyl * F0 / BigR * v * u_p * xjac                      * theta * tstep

             amat_13 = - eta_num * (v_x * zj_x + v_y * zj_y)                   * xjac * theta * tstep  &
                       - eta_T * v * zj / BigR                                 * xjac * theta * tstep

             amat_16 = - deta_dT * v * T * (zj0 - current_source(ms,mt))/ BigR * xjac * theta * tstep


!------------------------------------------------------------ equation 2
             amat_22 = - BigR * r0_hat * (v_x * u_x + v_y * u_y) * xjac * (1.d0+zeta)                                   &
                       + r0_hat * BigR**2 * w0 * (v_s * u_t  - v_t  * u_s)                              * theta * tstep &
                       + BigR**2 * (u_x * u0_x + u_y * u0_y) * (v_x * r0_y_hat - v_y * r0_x_hat) * xjac * theta * tstep

             amat_21 = - v * (psi_s * zj0_t - psi_t * zj0_s )              * theta * tstep

             amat_23 = - v * (ps0_s * zj_t  - ps0_t * zj_s)                * theta * tstep

             amat_23_n = + eps_cyl * F0 / BigR * v * zj_p  * xjac               * theta * tstep

             amat_24 = r0_hat * BigR**2 * w  * ( v_s * u0_t - v_t * u0_s)  * theta * tstep  &
                     + BigR * ( v_x * w_x + v_y * w_y) * visco_T  * xjac   * theta * tstep  &
                     - visco_num * (v_s * w_s + v_t * w_t)                         * tstep

             amat_25 = + 0.5d0 * vv2 * (v_x * rho_y_hat - v_y * rho_x_hat)   * xjac * theta * tstep &
                       + rho_hat * BigR**2 * w0 * (v_s * u0_t - v_t * u0_s)         * theta * tstep &
                       - BigR**2 * (v_s * rho_t * T0   - v_t * rho_s * T0  )        * theta * tstep &
                       - BigR**2 * (v_s * rho   * T0_t - v_t * rho   * T0_s)        * theta * tstep

             amat_26 = - BigR**2 * (v_s * r0_t * T   - v_t * r0_s * T)      * theta * tstep  &
                       - BigR**2 * (v_s * r0   * T_t - v_t * r0   * T_s)    * theta * tstep  &
                       + dvisco_dT * T * ( v_x * w0_x + v_y * w0_y ) * BigR * xjac * theta * tstep

!------------------------------------------------------------ equation 3
             amat_33 = v * zj / BigR * xjac                                * tstep
             amat_31 = (v_x * psi_x + v_y * psi_y ) / BigR * xjac          * tstep

!------------------------------------------------------------ equation 4
             amat_44 =  v * w * BigR * xjac                                * tstep
             amat_42 = (v_x * u_x + v_y * u_y) * BigR * xjac               * tstep

!------------------------------------------------------------ equation 5
             Bgrad_rho_star_psi = ( v_x  * psi_y - v_y  * psi_x ) / BigR
             Bgrad_rho_psi      = ( r0_x * psi_y - r0_y * psi_x ) / BigR
             Bgrad_rho_rho      = ( rho_x * ps0_y - rho_y * ps0_x ) / BigR
             Bgrad_rho_rho_n    = ( F0 / BigR * rho_p ) / BigR
             BB2_psi            = 2.d0 * (psi_x * ps0_x + psi_y * ps0_y ) /BigR**2


             amat_51 = - (D_par-D_prof) * BigR * BB2_psi/ BB2**2 * Bgrad_rho_star     * Bgrad_rho     * xjac * theta * tstep &
                       + (D_par-D_prof) * BigR / BB2             * Bgrad_rho_star_psi * Bgrad_rho     * xjac * theta * tstep &
                       + (D_par-D_prof) * BigR / BB2             * Bgrad_rho_star     * Bgrad_rho_psi * xjac * theta * tstep

             amat_51_k = - (D_par-D_prof) * BigR * BB2_psi/ BB2**2 * Bgrad_rho_k_star * Bgrad_rho     * xjac * theta * tstep &
                         + (D_par-D_prof) * BigR / BB2             * Bgrad_rho_k_star * Bgrad_rho_psi * xjac * theta * tstep

             amat_52 =  - v * BigR**2 * ( r0_s * u_t - r0_t * u_s)                                        * theta * tstep &
                        - v * 2.d0 * BigR * r0 * u_y                                               * xjac * theta * tstep

             amat_55 = v * rho * BigR * (1.d0 + zeta)                                              * xjac   &
                     - v * BigR**2 * ( rho_s * u0_t - rho_t * u0_s)                                       * theta * tstep &
                     - v * 2.d0 * BigR * rho * u0_y                                                * xjac * theta * tstep &
                     + (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star * Bgrad_rho_rho                * xjac * theta * tstep &
                     + D_prof * BigR  * (v_x*rho_x + v_y*rho_y )                                   * xjac * theta * tstep

             amat_55_k = + (D_par-D_prof) * BigR / BB2 * Bgrad_rho_k_star * Bgrad_rho_rho          * xjac * theta * tstep
             amat_55_n = + (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star   * Bgrad_rho_rho_n        * xjac * theta * tstep

             amat_55_kn = + (D_par-D_prof) * BigR / BB2 * Bgrad_rho_k_star * Bgrad_rho_rho_n       * xjac * theta * tstep &
                          + D_prof * BigR  * ( v_p*rho_p * eps_cyl**2 /BigR**2 )                   * xjac * theta * tstep


!------------------------------------------------------------ equation 6
             Bgrad_T_star_psi = ( v_x  * psi_y - v_y  * psi_x  ) / BigR
             Bgrad_T_psi      = ( T0_x * psi_y - T0_y * psi_x )  / BigR
             Bgrad_T_T        = ( T_x * ps0_y - T_y * ps0_x ) / BigR
             Bgrad_T_T_n      = ( F0 / BigR * T_p) / BigR

             amat_61 = - (ZK_par-ZK_prof) * BigR * BB2_psi / BB2**2 * Bgrad_T_star     * Bgrad_T     * xjac * theta * tstep &
                       + (ZK_par-ZK_prof) * BigR / BB2              * Bgrad_T_star_psi * Bgrad_T     * xjac * theta * tstep &
                       + (ZK_par-ZK_prof) * BigR / BB2              * Bgrad_T_star     * Bgrad_T_psi * xjac * theta * tstep

             amat_61_k = - (ZK_par-ZK_prof) * BigR * BB2_psi / BB2**2 * Bgrad_T_k_star * Bgrad_T     * xjac * theta * tstep &
                         + (ZK_par-ZK_prof) * BigR / BB2              * Bgrad_T_k_star * Bgrad_T_psi * xjac * theta * tstep

             amat_62 = - v * BigR**2 * ( T0_x * u_y - T0_y * u_x)             * xjac * theta * tstep &
                       - v * 2.d0* (GAMMA-1.d0) * BigR * T0 * u_y             * xjac * theta * tstep

             amat_66 =   v * T   * BigR * xjac * (1.d0 + zeta)                                           &
                       - v * BigR**2 * ( T_s * u0_t - T_t * u0_s)                        * theta * tstep &
                       - v * 2.d0* (GAMMA-1.d0) * BigR * T * u0_y                 * xjac * theta * tstep &
                       + (ZK_par-ZK_prof) * BigR / BB2 * Bgrad_T_star * Bgrad_T_T * xjac * theta * tstep &
                       + ZK_prof * BigR * ( v_x*T_x + v_y*T_y )                   * xjac * theta * tstep

             amat_66_k = + (ZK_par-ZK_prof) * BigR / BB2 * Bgrad_T_k_star * Bgrad_T_T    * xjac * theta * tstep
             amat_66_n = + (ZK_par-ZK_prof) * BigR / BB2 * Bgrad_T_star   * Bgrad_T_T_n  * xjac * theta * tstep

             amat_66_kn = + (ZK_par-ZK_prof) * BigR / BB2 * Bgrad_T_k_star * Bgrad_T_T_n * xjac * theta * tstep &
                          + ZK_prof * BigR   * (v_p*T_p /BigR**2 )                       * xjac * theta * tstep


             kl1 = index_kl
             kl2 = index_kl + 1
             kl3 = index_kl + 2
             kl4 = index_kl + 3
             kl5 = index_kl + 4
             kl6 = index_kl + 5

             ELM_p(mp,ij1,kl1)  =  ELM_p(mp,ij1,kl1) + wst * amat_11
             ELM_n(mp,ij1,kl2)  =  ELM_n(mp,ij1,kl2) + wst * amat_12_n

             ELM_p(mp,ij1,kl2)  =  ELM_p(mp,ij1,kl2) + wst * amat_12
             ELM_p(mp,ij1,kl3)  =  ELM_p(mp,ij1,kl3) + wst * amat_13
             ELM_p(mp,ij1,kl6)  =  ELM_p(mp,ij1,kl6) + wst * amat_16


             ELM_p(mp,ij2,kl1)  =  ELM_p(mp,ij2,kl1) + wst * amat_21
             ELM_p(mp,ij2,kl2)  =  ELM_p(mp,ij2,kl2) + wst * amat_22
             ELM_p(mp,ij2,kl3)  =  ELM_p(mp,ij2,kl3) + wst * amat_23
             ELM_n(mp,ij2,kl3)  =  ELM_n(mp,ij2,kl3) + wst * amat_23_n

             ELM_p(mp,ij2,kl4)  =  ELM_p(mp,ij2,kl4) + wst * amat_24
             ELM_p(mp,ij2,kl5)  =  ELM_p(mp,ij2,kl5) + wst * amat_25
             ELM_p(mp,ij2,kl6)  =  ELM_p(mp,ij2,kl6) + wst * amat_26


             ELM_p(mp,ij3,kl1)  =  ELM_p(mp,ij3,kl1) + wst * amat_31
             ELM_p(mp,ij3,kl3)  =  ELM_p(mp,ij3,kl3) + wst * amat_33

             ELM_p(mp,ij4,kl2)  =  ELM_p(mp,ij4,kl2) + wst * amat_42
             ELM_p(mp,ij4,kl4)  =  ELM_p(mp,ij4,kl4) + wst * amat_44

 
             ELM_p(mp,ij5,kl1)  =  ELM_p(mp,ij5,kl1)  + wst * amat_51
             ELM_k(mp,ij5,kl1)  =  ELM_k(mp,ij5,kl1)  + wst * amat_51_k

             ELM_p(mp,ij5,kl2)  =  ELM_p(mp,ij5,kl2)  + wst * amat_52

             ELM_p(mp,ij5,kl5)  =  ELM_p(mp,ij5,kl5)  + wst * amat_55
             ELM_k(mp,ij5,kl5)  =  ELM_k(mp,ij5,kl5)  + wst * amat_55_k
             ELM_n(mp,ij5,kl5)  =  ELM_n(mp,ij5,kl5)  + wst * amat_55_n
             ELM_kn(mp,ij5,kl5) =  ELM_kn(mp,ij5,kl5) + wst * amat_55_kn


             ELM_p(mp,ij6,kl1)  =  ELM_p(mp,ij6,kl1)  + wst * amat_61
             ELM_k(mp,ij6,kl1)  =  ELM_k(mp,ij6,kl1)  + wst * amat_61_k

             ELM_p(mp,ij6,kl2)  =  ELM_p(mp,ij6,kl2)  + wst * amat_62

             ELM_p(mp,ij6,kl6)  =  ELM_p(mp,ij6,kl6)  + wst * amat_66
             ELM_k(mp,ij6,kl6)  =  ELM_k(mp,ij6,kl6)  + wst * amat_66_k
             ELM_n(mp,ij6,kl6)  =  ELM_n(mp,ij6,kl6)  + wst * amat_66_n
             ELM_kn(mp,ij6,kl6) =  ELM_kn(mp,ij6,kl6) + wst * amat_66_kn

           enddo

         enddo
       enddo

     enddo
   enddo

 enddo
enddo


do i=1,n_vertex_max*n_var*(n_order+1)

  do j=1, n_vertex_max*n_var*(n_order+1)

!    call dfftw_plan_dft_r2c_1d(plan,n_plane,in_fft,out_fft,FFTW_FORWARD,FFTW_ESTIMATE)

    in_fft =  ELM_p(1:n_plane,i,j)

!    call dfftw_execute(plan)
    
    call my_fft(in_fft,out_fft,n_plane)
    
    do k=1,(n_tor+1)/2

      index_k = n_tor*(i-1) + max(2*(k-1),1)

      do m=1,(n_tor+1)/2

        index_m = n_tor*(j-1) + max(2*(m-1),1)

        l = (k-1) + (m-1)

        if ( (l .ge. 0) .and. (l .lt. n_plane/2) ) then

          ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + real(out_fft(l+1))
          ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   - imag(out_fft(l+1))
          ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) - imag(out_fft(l+1))
          ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) - real(out_fft(l+1))

        elseif ( (l .lt. 0) .and. (abs(l) .lt. n_plane/2) ) then

          ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + real(out_fft(abs(l)+1))
          ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + imag(out_fft(abs(l)+1))
          ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + imag(out_fft(abs(l)+1))
          ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) - real(out_fft(abs(l)+1))

        endif

        l = (k-1) - (m-1)

        if ( (l .ge. 0) .and. (l .lt. n_plane/2) ) then

          ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + real(out_fft(l+1))
          ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   - imag(out_fft(l+1))
          ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + imag(out_fft(l+1))
          ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + real(out_fft(l+1))

        elseif ( (l .lt. 0) .and. (abs(l) .lt. n_plane/2) ) then

          ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + real(out_fft(abs(l)+1))
          ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + imag(out_fft(abs(l)+1))
          ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) - imag(out_fft(abs(l)+1))
          ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + real(out_fft(abs(l)+1))

        endif

      enddo

    enddo

!    call dfftw_destroy_plan(plan)

  enddo

enddo

do i=1,n_vertex_max*n_var*(n_order+1)

  do j=1, n_vertex_max*n_var*(n_order+1)

  if (maxval(abs(ELM_n(1:n_plane,i,j))) .ne. 0.d0) then

!    call dfftw_plan_dft_r2c_1d(plan,n_plane,in_fft,out_fft,FFTW_FORWARD,FFTW_ESTIMATE)

    in_fft =  ELM_n(1:n_plane,i,j)

!    call dfftw_execute(plan)

    call my_fft(in_fft,out_fft,n_plane)

    do k=1,(n_tor+1)/2

      index_k = n_tor*(i-1) + max(2*(k-1),1)

      do m=1,(n_tor+1)/2

        im = max(2*(m-1),1)
        index_m = n_tor*(j-1) + max(2*(m-1),1)

        l = (k-1) + (m-1)

        if ( (l .ge. 0) .and. (l .lt. n_plane/2) ) then

          ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + imag(out_fft(l+1)) * float(mode(im))
          ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + real(out_fft(l+1)) * float(mode(im))
          ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + real(out_fft(l+1)) * float(mode(im))
          ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) - imag(out_fft(l+1)) * float(mode(im))

        elseif ( (l .lt. 0) .and. (abs(l) .lt. n_plane/2) ) then

          ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   - imag(out_fft(abs(l)+1)) * float(mode(im))
          ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + real(out_fft(abs(l)+1)) * float(mode(im))
          ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + real(out_fft(abs(l)+1)) * float(mode(im))
          ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + imag(out_fft(abs(l)+1)) * float(mode(im))

        endif

        l = (k-1) - (m-1)

        if ( (l .ge. 0) .and. (l .lt. n_plane/2) ) then

          ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   - imag(out_fft(l+1)) * float(mode(im))
          ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   - real(out_fft(l+1)) * float(mode(im))
          ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + real(out_fft(l+1)) * float(mode(im))
          ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) - imag(out_fft(l+1)) * float(mode(im))

        elseif ( (l .lt. 0) .and. (abs(l) .lt. n_plane/2) ) then

          ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + imag(out_fft(abs(l)+1)) * float(mode(im))
          ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   - real(out_fft(abs(l)+1)) * float(mode(im))
          ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + real(out_fft(abs(l)+1)) * float(mode(im))
          ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + imag(out_fft(abs(l)+1)) * float(mode(im))

        endif

      enddo

    enddo

!    call dfftw_destroy_plan(plan)

  endif
  enddo

enddo

do i=1,n_vertex_max*n_var*(n_order+1)

  do j=1, n_vertex_max*n_var*(n_order+1)

  if (maxval(abs(ELM_k(1:n_plane,i,j))) .ne. 0.d0) then

!    call dfftw_plan_dft_r2c_1d(plan,n_plane,in_fft,out_fft,FFTW_FORWARD,FFTW_ESTIMATE)

    in_fft =  ELM_k(1:n_plane,i,j)

!    call dfftw_execute(plan)

    call my_fft(in_fft,out_fft,n_plane)

    do k=1,(n_tor+1)/2

      ik      = max(2*(k-1),1)
      index_k = n_tor*(i-1) + max(2*(k-1),1)

      do m=1,(n_tor+1)/2

        index_m = n_tor*(j-1) + max(2*(m-1),1)

        l = (k-1) + (m-1)

        if ( (l .ge. 0) .and. (l .lt. n_plane/2) ) then

          ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + imag(out_fft(l+1)) * float(mode(ik))
          ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + real(out_fft(l+1)) * float(mode(ik))
          ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + real(out_fft(l+1)) * float(mode(ik))
          ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) - imag(out_fft(l+1)) * float(mode(ik))

        elseif ( (l .lt. 0) .and. (abs(l) .lt. n_plane/2) ) then

          ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   - imag(out_fft(abs(l)+1)) * float(mode(ik))
          ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + real(out_fft(abs(l)+1)) * float(mode(ik))
          ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + real(out_fft(abs(l)+1)) * float(mode(ik))
          ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + imag(out_fft(abs(l)+1)) * float(mode(ik))

        endif

        l = (k-1) - (m-1)

        if ( (l .ge. 0) .and. (l .lt. n_plane/2) ) then

          ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + imag(out_fft(l+1)) * float(mode(ik))
          ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + real(out_fft(l+1)) * float(mode(ik))
          ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) - real(out_fft(l+1)) * float(mode(ik))
          ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + imag(out_fft(l+1)) * float(mode(ik))

        elseif ( (l .lt. 0) .and. (abs(l) .lt. n_plane/2) ) then

          ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   - imag(out_fft(abs(l)+1)) * float(mode(ik))
          ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + real(out_fft(abs(l)+1)) * float(mode(ik))
          ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) - real(out_fft(abs(l)+1)) * float(mode(ik))
          ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) - imag(out_fft(abs(l)+1)) * float(mode(ik))

        endif

      enddo

    enddo

!    call dfftw_destroy_plan(plan)

  endif
  enddo

enddo


do i=1,n_vertex_max*n_var*(n_order+1)

  do j=1, n_vertex_max*n_var*(n_order+1)

  if (maxval(abs(ELM_kn(1:n_plane,i,j))) .ne. 0.d0) then

!    call dfftw_plan_dft_r2c_1d(plan,n_plane,in_fft,out_fft,FFTW_FORWARD,FFTW_ESTIMATE)

    in_fft =  ELM_kn(1:n_plane,i,j)

!    call dfftw_execute(plan)

    call my_fft(in_fft,out_fft,n_plane)

    do k=1,(n_tor+1)/2

      ik      = max(2*(k-1),1)
      index_k = n_tor*(i-1) + max(2*(k-1),1)

      do m=1,(n_tor+1)/2

        im      = max(2*(m-1),1)
        index_m = n_tor*(j-1) + max(2*(m-1),1)

        l = (k-1) + (m-1)

        if ( (l .ge. 0) .and. (l .lt. n_plane/2) ) then

           ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   - real(out_fft(l+1)) * float(mode(im)) * float(mode(ik))
           ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + imag(out_fft(l+1)) * float(mode(im)) * float(mode(ik))
           ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + imag(out_fft(l+1)) * float(mode(im)) * float(mode(ik))
           ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + real(out_fft(l+1)) * float(mode(im)) * float(mode(ik))

        elseif ( (l .lt. 0) .and. (abs(l) .lt. n_plane/2) ) then

           ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   - real(out_fft(abs(l)+1)) * float(mode(im)) * float(mode(ik))
           ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   - imag(out_fft(abs(l)+1)) * float(mode(im)) * float(mode(ik))
           ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) - imag(out_fft(abs(l)+1)) * float(mode(im)) * float(mode(ik))
           ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + real(out_fft(abs(l)+1)) * float(mode(im)) * float(mode(ik))

        endif

        l = (k-1) - (m-1)

        if ( (l .ge. 0) .and. (l .lt. n_plane/2) ) then

          ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + real(out_fft(l+1)) * float(mode(im)) * float(mode(ik))
          ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   - imag(out_fft(l+1)) * float(mode(im)) * float(mode(ik))
          ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + imag(out_fft(l+1)) * float(mode(im)) * float(mode(ik))
          ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + real(out_fft(l+1)) * float(mode(im)) * float(mode(ik))

        elseif ( (l .lt. 0) .and. (abs(l) .lt. n_plane/2) ) then

          ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + real(out_fft(abs(l)+1)) * float(mode(im)) * float(mode(ik))
          ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + imag(out_fft(abs(l)+1)) * float(mode(im)) * float(mode(ik))
          ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) - imag(out_fft(abs(l)+1)) * float(mode(im)) * float(mode(ik))
          ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + real(out_fft(abs(l)+1)) * float(mode(im)) * float(mode(ik))

        endif

      enddo

    enddo

!    call dfftw_destroy_plan(plan)

  endif
  enddo

enddo

ELM = 0.5d0 * ELM

do j=1, n_vertex_max*n_var*(n_order+1)

!  call dfftw_plan_dft_r2c_1d(plan,n_plane,in_fft,out_fft,FFTW_FORWARD,FFTW_ESTIMATE)

  in_fft = RHS_p(1:n_plane,j)

!  call dfftw_execute(plan)

  call my_fft(in_fft,out_fft,n_plane)

  index = n_tor*(j-1) + 1

  RHS(index) = real(out_fft(1))

  do k=2,(n_tor+1)/2

    index = n_tor*(j-1) + 2*(k-1)

    RHS(index)   =   real(out_fft(k))
    RHS(index+1) = - imag(out_fft(k))

  enddo

!  call dfftw_destroy_plan(plan)
 
enddo

do j=1, n_vertex_max*n_var*(n_order+1)

!  call dfftw_plan_dft_r2c_1d(plan,n_plane,in_fft,out_fft,FFTW_FORWARD,FFTW_ESTIMATE)

  in_fft = RHS_k(1:n_plane,j)

!  call dfftw_execute(plan)
  call my_fft(in_fft,out_fft,n_plane)

  index = n_tor*(j-1) + 1
  ik    = 1

  RHS(index) = RHS(index) + imag(out_fft(1)) * float(mode(ik))

  do k=2,(n_tor+1)/2

    ik    = max(2*(k-1),1)
    index = n_tor*(j-1) + 2*(k-1)

    RHS(index)   = RHS(index)   + imag(out_fft(k)) * float(mode(ik))
    RHS(index+1) = RHS(index+1) + real(out_fft(k)) * float(mode(ik))

  enddo

!  call dfftw_destroy_plan(plan)

enddo

return
end

subroutine my_fft(in_fft,out_fft,n)
      
implicit none
      
real*8     :: in_fft(*)
complex*16 :: out_fft(*)
integer    :: n
      
real*8     :: tmp_fft(2*n+2)
integer    :: i
      
tmp_fft(1:n) = in_fft(1:n)
      
call RFT2(tmp_fft,n,1)
      
do i=1,n
  out_fft(i) = cmplx(tmp_fft(2*i-1),tmp_fft(2*i))
enddo
      
return
end
      
