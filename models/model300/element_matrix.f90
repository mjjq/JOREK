module mod_elt_matrix
  implicit none
contains

subroutine element_matrix(element,nodes, xpoint2, xcase2, psi_axis, psi_bnd, Z_xpoint, ELM, RHS, tid)
!---------------------------------------------------------------
! calculates the matrix contribution of one element
!---------------------------------------------------------------
use parameters
use data_structure
use gauss
use basis_at_gaussian
use phys_module
use profiles, only: interpolProf

implicit none

type (type_element)   :: element
type (type_node)      :: nodes(n_vertex_max)

real*8, dimension (:,:), pointer  :: ELM
real*8, dimension (:)  , pointer  :: RHS
integer, intent(in) :: tid

integer    :: i, j, ms, mt, mp, k, l, index_ij, index_kl, index, xcase2
integer    :: in, im, ij1, ij2, ij3, ij4, ij5, ij6, ij7, kl1, kl2, kl3, kl4, kl5, kl6, kl7
real*8     :: wst, xjac, xjac_s, xjac_t, BigR, r2, PI, phi, eps_cyl
real*8     :: current_source(n_gauss,n_gauss),particle_source(n_gauss,n_gauss),heat_source(n_gauss,n_gauss)
real*8     :: psi_axis, psi_bnd, Z_xpoint(2), dj_dpsi, dj_dz
real*8     :: Bgrad_rho_star,     Bgrad_rho,     Bgrad_T_star,  Bgrad_T, BB2
real*8     :: Bgrad_rho_star_psi, Bgrad_rho_psi, Bgrad_rho_rho, Bgrad_T_star_psi, Bgrad_T_psi, Bgrad_T_T, BB2_psi
real*8     :: rhs_ij_1,   rhs_ij_2,   rhs_ij_3,   rhs_ij_4,   rhs_ij_5,   rhs_ij_6, rhs_ij_7
real*8     :: rhs_stab_1, rhs_stab_2, rhs_stab_3, rhs_stab_4, rhs_stab_5, rhs_stab_6
real*8     :: ZK_prof, D_prof, psi_norm, theta, zeta, delta_u_x, delta_u_y

real*8     :: v, v_x, v_y, v_s, v_t, v_p, v_ss, v_st, v_tt, v_xx, v_yy, v_xs, v_ys, v_xt, v_yt
real*8     :: ps0, ps0_x, ps0_y, ps0_p,ps0_s,ps0_t,  zj0, zj0_x, zj0_y, zj0_p, zj0_s, zj0_t
real*8     :: u0, u0_x, u0_y, u0_p, u0_s, u0_t,  w0, w0_x, w0_y, w0_p, w0_s, w0_t
real*8     :: Vpar0, Vpar0_s, Vpar0_t, Vpar0_p, Vpar0_x, Vpar0_y, Vpar0_ss, Vpar0_st, Vpar0_tt
real*8     :: r0, r0_x, r0_y, r0_p, r0_s, r0_t, r0_ss, r0_st, r0_tt, r0_hat, r0_x_hat, r0_y_hat
real*8     :: T0, T0_x, T0_y, T0_p, T0_s, T0_t
real*8     :: psi, psi_x, psi_y, psi_p, psi_s, psi_t, psi_ss, psi_st, psi_tt, psi_xs, psi_ys, psi_xt, psi_yt, psi_xx, psi_yy
real*8     :: zj, zj_x, zj_y, zj_p, zj_s, zj_t
real*8     :: Vpar, Vpar_x, Vpar_y, Vpar_p, Vpar_s, Vpar_t, Vpar_ss, Vpar_st, Vpar_tt
real*8     :: u, u_x, u_y, u_p, u_s, u_t, w, w_x, w_y, w_p, w_s, w_t
real*8     :: rho, rho_x, rho_y, rho_s, rho_t, rho_p, rho_ss, rho_st, rho_tt, rho_hat, rho_x_hat, rho_y_hat
real*8     :: T, T_x, T_y, T_s, T_t, T_p
real*8     :: P0, P0_x, P0_y, P0_s, P0_t, P0_p
real*8     :: w0_ss, w0_st, w0_tt, w_ss, w_st, w_tt
real*8     :: BigR_x, vv2, eta_T, visco_T, deta_dT, d2eta_d2T, dvisco_dT, visco_num_T, eta_num_T
real*8     :: amat_11, amat_12, amat_21, amat_22, amat_23, amat_24, amat_25, amat_26, amat_33, amat_31, amat_44, amat_42
real*8     :: amat_51, amat_52, amat_55, amat_57, amat_61, amat_62, amat_66, amat_67, amat_16, amat_13
real*8     :: amat_71, amat_72, amat_75, amat_76, amat_77
real*8     :: amat_stab_11, amat_stab_12, amat_stab_13, amat_stab_14 ,amat_stab_21,amat_stab_22, amat_stab_23, amat_stab_24
real*8     :: amat_stab_31, amat_stab_32, amat_stab_33, amat_stab_34 ,amat_stab_41,amat_stab_42, amat_stab_43, amat_stab_44
logical    :: xpoint2

real*8, dimension(n_gauss,n_gauss)    :: x_g, x_s, x_t
real*8, dimension(n_gauss,n_gauss)    :: x_ss, x_st, x_tt
real*8, dimension(n_gauss,n_gauss)    :: y_g, y_s, y_t
real*8, dimension(n_gauss,n_gauss)    :: y_ss, y_st, y_tt

real*8, dimension(:,:,:,:) , pointer :: eq_g, eq_s, eq_t
real*8, dimension(:,:,:,:) , pointer :: eq_p
real*8, dimension(:,:,:,:) , pointer :: eq_ss, eq_st, eq_tt   
real*8, dimension(:,:,:,:) , pointer :: delta_g, delta_s, delta_t

eq_g    => thread_struct(tid)%eq_g   
eq_s    => thread_struct(tid)%eq_s   
eq_t    => thread_struct(tid)%eq_t   
eq_p    => thread_struct(tid)%eq_p   
eq_ss   => thread_struct(tid)%eq_ss  
eq_st   => thread_struct(tid)%eq_st  
eq_tt   => thread_struct(tid)%eq_tt  
delta_g => thread_struct(tid)%delta_g
delta_s => thread_struct(tid)%delta_s
delta_t => thread_struct(tid)%delta_t

ELM = 0.d0
RHS = 0.d0

PI    = 2.d0*asin(1.d0)

!theta = 0.5d0 ! Crank-Nicholson parameter
!theta = 1.0d0  ; zeta = 0.0d0      ! Euler scheme
theta = 1.0d0   ; zeta = 0.5d0      ! BDF2 (Gears) scheme

!---------------------------------------------------- value of (x,y) and derivatives on Gaussian points
x_g  = 0.d0; x_s  = 0.d0; x_t  = 0.d0; x_st  = 0.d0; x_ss  = 0.d0; x_tt  = 0.d0;
y_g  = 0.d0; y_s  = 0.d0; y_t  = 0.d0; y_st  = 0.d0; y_ss  = 0.d0; y_tt  = 0.d0;
eq_g = 0.d0; eq_s = 0.d0; eq_t = 0.d0; eq_st = 0.d0; eq_ss = 0.d0; eq_tt = 0.d0; eq_p = 0.d0;

delta_g = 0.d0; delta_s = 0.d0; delta_t = 0.d0

current_source  = 0.d0
particle_source = 0.d0
heat_source     = 0.d0

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

             eq_ss(mp,k,ms,mt) = eq_ss(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_ss(i,j,ms,mt)* HZ(in,mp)
             eq_st(mp,k,ms,mt) = eq_st(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_st(i,j,ms,mt)* HZ(in,mp)
             eq_tt(mp,k,ms,mt) = eq_tt(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_tt(i,j,ms,mt)* HZ(in,mp)

             delta_g(mp,k,ms,mt) = delta_g(mp,k,ms,mt) + nodes(i)%deltas(in,j,k) * element%size(i,j) * H(i,j,ms,mt)   * HZ(in,mp)
             delta_s(mp,k,ms,mt) = delta_s(mp,k,ms,mt) + nodes(i)%deltas(in,j,k) * element%size(i,j) * H_s(i,j,ms,mt) * HZ(in,mp)
             delta_t(mp,k,ms,mt) = delta_t(mp,k,ms,mt) + nodes(i)%deltas(in,j,k) * element%size(i,j) * H_t(i,j,ms,mt) * HZ(in,mp)

           enddo

         enddo

       enddo

       call current(xpoint2, xcase2, x_g(ms,mt),y_g(ms,mt), Z_xpoint, eq_g(1,1,ms,mt),psi_axis,psi_bnd,current_source(ms,mt))
       call sources(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, eq_g(1,1,ms,mt),psi_axis,psi_bnd,particle_source(ms,mt),heat_source(ms,mt))

     enddo
   enddo
 enddo
enddo


!--------------------------------------------------- sum over the Gaussian integration points
do ms=1, n_gauss

 do mt=1, n_gauss

   wst = wgauss(ms)*wgauss(mt)

   xjac    = x_s(ms,mt)*y_t(ms,mt)  - x_t(ms,mt)*y_s(ms,mt)

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
     w0_ss = eq_ss(mp,4,ms,mt)
     w0_tt = eq_tt(mp,4,ms,mt)
     w0_st = eq_st(mp,4,ms,mt)

     r0    = eq_g(mp,5,ms,mt)
     r0_x  = (   y_t(ms,mt) * eq_s(mp,5,ms,mt) - y_s(ms,mt) * eq_t(mp,5,ms,mt) ) / xjac
     r0_y  = ( - x_t(ms,mt) * eq_s(mp,5,ms,mt) + x_s(ms,mt) * eq_t(mp,5,ms,mt) ) / xjac
     r0_p  = eq_p(mp,5,ms,mt)
     r0_s  = eq_s(mp,5,ms,mt)
     r0_t  = eq_t(mp,5,ms,mt)
     r0_ss = eq_ss(mp,5,ms,mt)
     r0_st = eq_st(mp,5,ms,mt)
     r0_tt = eq_tt(mp,5,ms,mt)

     r0_hat   = BigR**2 * abs(r0)
     r0_x_hat = 2.d0 * BigR * BigR_x  * r0 + BigR**2 * r0_x
     r0_y_hat = BigR**2 * r0_y

     T0    = eq_g(mp,6,ms,mt)
     T0_x  = (   y_t(ms,mt) * eq_s(mp,6,ms,mt) - y_s(ms,mt) * eq_t(mp,6,ms,mt) ) / xjac
     T0_y  = ( - x_t(ms,mt) * eq_s(mp,6,ms,mt) + x_s(ms,mt) * eq_t(mp,6,ms,mt) ) / xjac
     T0_p  = eq_p(mp,6,ms,mt)
     T0_s  = eq_s(mp,6,ms,mt)
     T0_t  = eq_t(mp,6,ms,mt)

     Vpar0    = eq_g(mp,7,ms,mt)
     Vpar0_x  = (   y_t(ms,mt) * eq_s(mp,7,ms,mt) - y_s(ms,mt) * eq_t(mp,7,ms,mt) ) / xjac
     Vpar0_y  = ( - x_t(ms,mt) * eq_s(mp,7,ms,mt) + x_s(ms,mt) * eq_t(mp,7,ms,mt) ) / xjac
     Vpar0_p  = eq_p(mp,7,ms,mt)
     Vpar0_s  = eq_s(mp,7,ms,mt)
     Vpar0_t  = eq_t(mp,7,ms,mt)
     Vpar0_ss = eq_ss(mp,7,ms,mt)
     Vpar0_st = eq_st(mp,7,ms,mt)
     Vpar0_tt = eq_tt(mp,7,ms,mt)

     P0    = abs(r0 * T0)
     P0_x  = r0_x * T0 + r0 * T0_x
     P0_y  = r0_y * T0 + r0 * T0_y
     P0_s  = r0_s * T0 + r0 * T0_s
     P0_t  = r0_t * T0 + r0 * T0_t
     P0_p  = r0_p * T0 + r0 * T0_p

     delta_u_x = (   y_t(ms,mt) * delta_s(mp,2,ms,mt) - y_s(ms,mt) * delta_t(mp,2,ms,mt) ) / xjac
     delta_u_y = ( - x_t(ms,mt) * delta_s(mp,2,ms,mt) + x_s(ms,mt) * delta_t(mp,2,ms,mt) ) / xjac

     ! --- Temperature dependent resistivity
     if ( eta_T_dependent ) then
       eta_T     = eta   * (abs(T0)/T_0)**(-1.5d0)
       deta_dT   = - eta   * (1.5d0)  * abs(T0)**(-2.5d0) * T_0**(1.5d0)
       d2eta_d2T =   eta   * (3.75d0) * abs(T0)**(-3.5d0) * T_0**(1.5d0)
     else
       eta_T     = eta
       deta_dT   = 0.d0
       d2eta_d2T = 0.d0
     end if
     
     ! --- Temperature dependent viscosity
     if ( visco_T_dependent ) then
       visco_T   = visco * (abs(T0)/T_0)**(-1.5d0)
       dvisco_dT = - visco * (1.5d0)  * abs(T0)**(-2.5d0) * T_0**(1.5d0)
     else
       visco_T   = visco
       dvisco_dT = 0.d0
     end if
     
     eta_num_T   = eta_num                         ! hyperresistivity
     visco_num_T = visco_num                       ! hyperviscosity

     psi_norm = (ps0 - psi_axis)/(psi_bnd - psi_axis)
     if (xpoint2) then
       if ((psi_norm .lt. 1.d0) .and. (y_g(ms,mt) .lt. Z_xpoint(1)) .and. (xcase2 .ne. 2)) then
         psi_norm = 2.d0 - psi_norm
       endif
       if ((psi_norm .lt. 1.d0) .and. (y_g(ms,mt) .gt. Z_xpoint(2)) .and. (xcase2 .ne. 1)) then
         psi_norm = 2.d0 - psi_norm
       endif
     endif

     if ( num_d_perp ) then
       D_prof  = interpolProf(num_d_perp_x, num_d_perp_y, num_d_perp_len, psi_norm)
     else
       D_prof  = D_perp(1)  * ((1.d0-D_perp(2))  + D_perp(2)  *(0.5d0 - 0.5d0*tanh((psi_norm-D_perp(5)) /D_perp(4))))
     end if
     
     if ( num_zk_perp ) then
       ZK_prof = interpolProf(num_zk_perp_x, num_zk_perp_y, num_zk_perp_len, psi_norm)
     else
       ZK_prof = ZK_perp(1) * ((1.d0-ZK_perp(2)) + ZK_perp(2) *(0.5d0 - 0.5d0*tanh((psi_norm-ZK_perp(5))/ZK_perp(4))))
     end if

     do i=1,n_vertex_max

       do j=1,n_order+1

         do im=1,n_tor

           index_ij = n_tor*n_var*(n_order+1)*(i-1) + n_tor * n_var * (j-1) + im   ! index in the ELM matrix

           v   =  H(i,j,ms,mt) * element%size(i,j) * HZ(im,mp)
           v_x = (  y_t(ms,mt) * h_s(i,j,ms,mt) - y_s(ms,mt) * h_t(i,j,ms,mt) ) * element%size(i,j) / xjac * HZ(im,mp)
           v_y = (- x_t(ms,mt) * h_s(i,j,ms,mt) + x_s(ms,mt) * h_t(i,j,ms,mt) ) * element%size(i,j) / xjac * HZ(im,mp)

           v_s = h_s(i,j,ms,mt) * element%size(i,j) * HZ(im,mp)
           v_t = h_t(i,j,ms,mt) * element%size(i,j) * HZ(im,mp)
           v_p = H(i,j,ms,mt)   * element%size(i,j) * HZ_p(im,mp)

           v_ss = h_ss(i,j,ms,mt) * element%size(i,j) * HZ(im,mp)
           v_tt = h_tt(i,j,ms,mt) * element%size(i,j) * HZ(im,mp)
           v_st = h_st(i,j,ms,mt) * element%size(i,j) * HZ(im,mp)

           Bgrad_rho_star = ( F0 / BigR * v_p  +  v_x  * ps0_y - v_y  * ps0_x ) / BigR    ! F0 due to absence of normalisation
           Bgrad_rho      = ( F0 / BigR * r0_p +  r0_x * ps0_y - r0_y * ps0_x ) / BigR    ! F0 due to absence of normalisation
           Bgrad_T_star   = ( F0 / BigR * v_p  +  v_x  * ps0_y - v_y  * ps0_x ) / BigR    ! F0 due to absence of normalisation
           Bgrad_T        = ( F0 / BigR * T0_p +  T0_x * ps0_y - T0_y * ps0_x ) / BigR    ! F0 due to absence of normalisation

           BB2 = (F0*F0 + ps0_x * ps0_x + ps0_y * ps0_y )/BigR**2                         ! F0 due to absence of normalisation

           rhs_ij_1 =   v * eta_T  * (zj0 - current_source(ms,mt))/ BigR  * xjac * tstep &
                      + v * (ps0_s * u0_t - ps0_t * u0_s)                        * tstep &
                      - v * eps_cyl * F0 / BigR  * u0_p                   * xjac * tstep &          ! F0 due to absence of normalisation
                      + eta_num_T * (v_x * zj0_x + v_y * zj0_y)           * xjac * tstep &
                      + zeta * v * delta_g(mp,1,ms,mt) / BigR             * xjac

           rhs_ij_2 = - 0.5d0 * vv2 * (v_x * r0_y_hat - v_y * r0_x_hat)   * xjac * tstep &
                      - r0_hat * BigR**2 * w0 * (v_s * u0_t - v_t * u0_s)        * tstep &
                      + v * (ps0_s * zj0_t - ps0_t * zj0_s )                     * tstep &
                      - visco_T * BigR * (v_x * w0_x + v_y * w0_y)        * xjac * tstep &
                      - v * eps_cyl * F0 / BigR * zj0_p                   * xjac * tstep &         ! F0 due to absence of normalisation
                      + BigR**2 * (v_s * p0_t - v_t * p0_s)                      * tstep &
                      - visco_num_T  * ( (v_ss * (x_t(ms,mt)**2+y_t(ms,mt)**2) + v_tt*(x_s(ms,mt)**2+y_s(ms,mt)**2)   &
                               -2.d0*  v_st * (x_s(ms,mt)*x_t(ms,mt) + y_s(ms,mt)*y_t(ms,mt)) )                       &
                                    * (w0_ss * (x_t(ms,mt)**2+y_t(ms,mt)**2) + w0_tt*(x_s(ms,mt)**2+y_s(ms,mt)**2)    &
                               -2.d0*  w0_st * (x_s(ms,mt)*x_t(ms,mt) + y_s(ms,mt)*y_t(ms,mt)) ) ) / xjac**3  * tstep &
                      - zeta * BigR * r0_hat * (v_x * delta_u_x + v_y * delta_u_y) * xjac

           rhs_ij_3 = 0.d0 !- ( v_x * ps0_x  + v_y * ps0_y + v*zj0) / BigR * xjac * tstep
           rhs_ij_4 = 0.d0 !- ( v_x * u0_x   + v_y * u0_y  + v*w0)  * BigR * xjac * tstep

           rhs_ij_5 = v * BigR * particle_source(ms,mt)                                        * xjac * tstep &
                    + v * BigR**2 * ( r0_s * u0_t - r0_t * u0_s)                                      * tstep &
                    + v * 2.d0 * BigR * r0 * u0_y                                              * xjac * tstep &
                    - (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star * Bgrad_rho                 * xjac * tstep &
                    - D_prof * BigR  * (v_x*r0_x + v_y*r0_y + v_p*r0_p * eps_cyl**2 /BigR**2 ) * xjac * tstep &
                    - v * F0 / BigR * Vpar0 * r0_p                                             * xjac * tstep &
                    - v * Vpar0 * (r0_s * ps0_t - r0_t * ps0_s)                                       * tstep &
                    - v * F0 / BigR * r0 * vpar0_p                                             * xjac * tstep &
                    - v * r0 * (vpar0_s * ps0_t - vpar0_t * ps0_s)                                    * tstep &
                    + zeta * v * delta_g(mp,5,ms,mt) * BigR                                    * xjac &

                    - D_perp_num  *                                                                   &
                           ( (v_ss * (x_t(ms,mt)**2 + y_t(ms,mt)**2)                                  &
                            + v_tt * (x_s(ms,mt)**2 + y_s(ms,mt)**2)                                  &
                     -2.d0 *  v_st * (x_s(ms,mt)*x_t(ms,mt) + y_s(ms,mt)*y_t(ms,mt)) )                &
                           * (r0_ss * (x_t(ms,mt)**2 + y_t(ms,mt)**2)                                 &
                            + r0_tt * (x_s(ms,mt)**2 + y_s(ms,mt)**2)                                 &
                     -2.d0 *  r0_st * (x_s(ms,mt)*x_t(ms,mt) + y_s(ms,mt)*y_t(ms,mt)) ) ) / xjac**3  * theta*tstep


           rhs_ij_6 = v * BigR * heat_source(ms,mt)                               * xjac * tstep &
                    + v * BigR**2 * ( T0_s * u0_t - T0_t * u0_s)                         * tstep &
                    + v * 2.d0* (GAMMA-1.d0) * BigR * T0 * u0_y                   * xjac * tstep &
                    - (ZK_par-ZK_prof) * BigR / BB2 * Bgrad_T_star * Bgrad_T      * xjac * tstep &
                    - ZK_prof * BigR * (v_x*T0_x + v_y*T0_y + v_p*T0_p /BigR**2 ) * xjac * tstep &
                    - v * F0 / BigR * Vpar0 * T0_p                                * xjac * tstep &
                    - v * Vpar0 * (T0_s * ps0_t - T0_t * ps0_s)                          * tstep &
                    - v * (GAMMA-1.d0) * T0 * (vpar0_s * ps0_t - vpar0_t * ps0_s)        * tstep &
                    - v * (GAMMA-1.d0) * T0 * F0 / BigR * vpar0_p                 * xjac * tstep &
                    + zeta * v * delta_g(mp,6,ms,mt) * BigR                       * xjac

           rhs_ij_7 = - v * F0 / BigR * P0_p                                      * xjac * tstep &
                      - v * (P0_s * ps0_t - P0_t * ps0_s)                                * tstep &
                      - visco_par * (v_x * vpar0_x + v_y * vpar0_y) * BigR        * xjac * tstep &
                      + zeta * v * delta_g(mp,7,ms,mt) * R0 * F0**2 / BigR        * xjac &
!                      - F0 * (u0_s * vpar0_t - u0_t * vpar0_s)                           * tstep &
!                      - vpar0 * (F0/BigR)**2 * (vpar0_s * ps0_t - vpar0_t * ps0_s)       * tstep &
!                      - (F0/BigR)**3 * vpar0 * vpar0_p                            * xjac * tstep

                    - visco_par_num  *                                                                &
                           ( (v_ss * (x_t(ms,mt)**2 + y_t(ms,mt)**2)                                  &
                            + v_tt * (x_s(ms,mt)**2 + y_s(ms,mt)**2)                                  &
                     -2.d0 *  v_st * (x_s(ms,mt)*x_t(ms,mt) + y_s(ms,mt)*y_t(ms,mt)) )                &
                           * (Vpar0_ss * (x_t(ms,mt)**2 + y_t(ms,mt)**2)                              &
                            + Vpar0_tt * (x_s(ms,mt)**2 + y_s(ms,mt)**2)                              &
                     -2.d0 *  Vpar0_st * (x_s(ms,mt)*x_t(ms,mt) + y_s(ms,mt)*y_t(ms,mt)) ) ) / xjac**3  * theta*tstep

           ij1 = index_ij
           ij2 = index_ij + 1*n_tor
           ij3 = index_ij + 2*n_tor
           ij4 = index_ij + 3*n_tor
           ij5 = index_ij + 4*n_tor
           ij6 = index_ij + 5*n_tor
           ij7 = index_ij + 6*n_tor

           RHS(ij1) = RHS(ij1) + rhs_ij_1 * wst
           RHS(ij2) = RHS(ij2) + rhs_ij_2 * wst
           RHS(ij3) = RHS(ij3) + rhs_ij_3 * wst
           RHS(ij4) = RHS(ij4) + rhs_ij_4 * wst
           RHS(ij5) = RHS(ij5) + rhs_ij_5 * wst
           RHS(ij6) = RHS(ij6) + rhs_ij_6 * wst
           RHS(ij7) = RHS(ij7) + rhs_ij_7 * wst

           do k=1,n_vertex_max

             do l=1,n_order+1

               do in = 1, n_tor

                 psi   = H(k,l,ms,mt) * element%size(k,l) * HZ(in,mp)

                 psi_x = (   y_t(ms,mt) * h_s(k,l,ms,mt) - y_s(ms,mt) * h_t(k,l,ms,mt) ) / xjac    &
                              * element%size(k,l) * HZ(in,mp)

                 psi_y = ( - x_t(ms,mt) * h_s(k,l,ms,mt) + x_s(ms,mt) * h_t(k,l,ms,mt) )  / xjac   &
                              * element%size(k,l) * HZ(in,mp)

                 psi_p = H(k,l,ms,mt)   * element%size(k,l) * HZ_p(in,mp)
                 psi_s = h_s(k,l,ms,mt) * element%size(k,l) * HZ(in,mp)
                 psi_t = h_t(k,l,ms,mt) * element%size(k,l) * HZ(in,mp)

                 u   = psi   ;    zj   = psi   ;    w   = psi    ;    rho   = psi    ;    T   = psi   ;    vpar   = psi
                 u_x = psi_x ;    zj_x = psi_x ;    w_x = psi_x  ;    rho_x = psi_x  ;    T_x = psi_x ;    vpar_x = psi_x
                 u_y = psi_y ;    zj_y = psi_y ;    w_y = psi_y  ;    rho_y = psi_y  ;    T_y = psi_y ;    vpar_y = psi_y
                 u_p = psi_p ;    zj_p = psi_p ;    w_p = psi_p  ;    rho_p = psi_p  ;    T_p = psi_p ;    vpar_p = psi_p
                 u_s = psi_s ;    zj_s = psi_s ;    w_s = psi_s  ;    rho_s = psi_s  ;    T_s = psi_s ;    vpar_s = psi_s
                 u_t = psi_t ;    zj_t = psi_t ;    w_t = psi_t  ;    rho_t = psi_t  ;    T_t = psi_t ;    vpar_t = psi_t

                 w_ss = h_ss(k,l,ms,mt) * element%size(k,l) * HZ(in,mp)
                 w_tt = h_tt(k,l,ms,mt) * element%size(k,l) * HZ(in,mp)
                 w_st = h_st(k,l,ms,mt) * element%size(k,l) * HZ(in,mp)

                 Vpar_ss = h_ss(k,l,ms,mt) * element%size(k,l)
                 Vpar_tt = h_tt(k,l,ms,mt) * element%size(k,l)
                 Vpar_st = h_st(k,l,ms,mt) * element%size(k,l)

                 rho_ss = h_ss(k,l,ms,mt) * element%size(k,l) * HZ(in,mp)
                 rho_tt = h_tt(k,l,ms,mt) * element%size(k,l) * HZ(in,mp)
                 rho_st = h_st(k,l,ms,mt) * element%size(k,l) * HZ(in,mp)

                 rho_hat   = BigR**2 * rho
                 rho_x_hat = 2.d0 * BigR * BigR_x  * rho + BigR**2 * rho_x
                 rho_y_hat = BigR**2 * rho_y


                 index_kl = n_tor*n_var*(n_order+1)*(k-1) + n_tor * n_var * (l-1) + in   ! index in the ELM matrix

!---------------------------------------------------------------- equation 1
                 amat_11 = v * psi / BigR * xjac * (1.d0 + zeta)                                     &
                         - v * (psi_s * u0_t - psi_t * u0_s)                        * theta * tstep

                 amat_12 = - v * (ps0_s * u_t - ps0_t * u_s)                        * theta * tstep  &
                           + eps_cyl * F0 / BigR * v * u_p * xjac                   * theta * tstep

                 amat_13 = - eta_num_T * (v_x * zj_x + v_y * zj_y)           * xjac * theta * tstep  &
                           - eta_T * v * zj / BigR                           * xjac * theta * tstep

                 amat_16 = - deta_dT * v * T * (zj0 - current_source(ms,mt)) / BigR * xjac * theta * tstep

!---------------------------------------------------------------- equation 1
!                 amat_11 = v * psi / BigR * xjac                                                 &
!                         + eta_T * (psi_x * v_x + psi_y * v_y) / BigR * xjac           * theta * tstep  &
!                         - v * (psi_s * u0_t - psi_t * u0_s)                           * theta * tstep  &
!                         + v * deta_dT * (T0_x * psi_x  + T0_y * psi_y ) / BigR * xjac * theta * tstep
!
!                 amat_12 = -  v * (ps0_s * u_t - ps0_t * u_s)                 * theta * tstep  &
!                           +  eps_cyl * F0 / BigR * v * u_p * xjac            * theta * tstep
!
!                 amat_13 = - eta_num * (v_s * zj_t + v_t * zj_s)              * theta * tstep
!
!                 amat_16 = + deta_dT * T * ( v_x * ps0_x + v_y * ps0_y ) / BigR         * xjac * theta * tstep &
!                           + v * deta_dT * ( T_x * ps0_x + T_y * ps0_y ) / BigR         * xjac * theta * tstep &
!                           + v * d2eta_d2T * T * ( T0_x * ps0_x + T0_y * ps0_y ) / BigR * xjac * theta * tstep &
!                           + v * deta_dT   * T * current_source(ms,mt)   / BigR         * xjac * theta * tstep
!
!---------------------------------------------------------------- equation 2
                 amat_22 = - BigR * r0_hat * (v_x * u_x + v_y * u_y) * xjac * (1.d0 + zeta)                                 &
                           + r0_hat * BigR**2 * w0 * (v_s * u_t  - v_t  * u_s)                              * theta * tstep &
                           + BigR**2 * (u_x * u0_x + u_y * u0_y) * (v_x * r0_y_hat - v_y * r0_x_hat) * xjac * theta * tstep

                 amat_21 = - v * (psi_s * zj0_t - psi_t * zj0_s)               * theta * tstep
                 amat_23 = - v * (ps0_s * zj_t  - ps0_t * zj_s)                * theta * tstep  &
                           + eps_cyl * F0 / BigR * v * zj_p  * xjac            * theta * tstep      ! F0 due to absence of normalisation

                 amat_24 = r0_hat * BigR**2 * w  * ( v_s * u0_t - v_t * u0_s)  * theta * tstep  &
                         + BigR * ( v_x * w_x + v_y * w_y) * visco_T  * xjac   * theta * tstep  &
                         + visco_num_T  * ( (v_ss * (x_t(ms,mt)**2+y_t(ms,mt)**2) + v_tt*(x_s(ms,mt)**2+y_s(ms,mt)**2) &
                               -2.d0*  v_st * (x_s(ms,mt)*x_t(ms,mt) + y_s(ms,mt)*y_t(ms,mt)) )                  &
                                    * (w_ss * (x_t(ms,mt)**2+y_t(ms,mt)**2) + w_tt*(x_s(ms,mt)**2+y_s(ms,mt)**2) &
                               -2.d0*  w_st * (x_s(ms,mt)*x_t(ms,mt) + y_s(ms,mt)*y_t(ms,mt)) ) ) / xjac**3  * theta*tstep

                 amat_25 = + 0.5d0 * vv2 * (v_x * rho_y_hat - v_y * rho_x_hat)   * xjac * theta * tstep &
                           + rho_hat * BigR**2 * w0 * (v_s * u0_t - v_t * u0_s)         * theta * tstep &
                           - BigR**2 * (v_s * rho_t * T0   - v_t * rho_s * T0  )        * theta * tstep &
                           - BigR**2 * (v_s * rho   * T0_t - v_t * rho   * T0_s)        * theta * tstep

                 amat_26 = - BigR**2 * (v_s * r0_t * T   - v_t * r0_s * T)      * theta * tstep  &
                           - BigR**2 * (v_s * r0   * T_t - v_t * r0   * T_s)    * theta * tstep  &
                           + dvisco_dT * T * ( v_x * w0_x + v_y * w0_y ) * BigR * xjac * theta * tstep

!---------------------------------------------------------------- equation 3
                 amat_33 = v * zj / BigR * xjac                               * tstep
                 amat_31 = (v_x * psi_x + v_y * psi_y ) / BigR * xjac         * tstep

!---------------------------------------------------------------- equation 4
                 amat_44 =  v * w * BigR * xjac                                * tstep
                 amat_42 = (v_x * u_x + v_y * u_y) * BigR * xjac               * tstep

!---------------------------------------------------------------- equation 5
                 Bgrad_rho_star_psi = ( v_x  * psi_y - v_y  * psi_x ) / BigR
                 Bgrad_rho_psi      = ( r0_x * psi_y - r0_y * psi_x ) / BigR
                 Bgrad_rho_rho      = ( F0 / BigR * rho_p +  rho_x * ps0_y - rho_y * ps0_x ) / BigR   ! F0 due to absence of normalisation
                 BB2_psi            = 2.d0 * (psi_x * ps0_x + psi_y * ps0_y ) /BigR**2

                 amat_51 = - (D_par-D_prof) * BigR * BB2_psi / BB2**2 * Bgrad_rho_star     * Bgrad_rho     * xjac * theta * tstep &
                           + (D_par-D_prof) * BigR / BB2              * Bgrad_rho_star_psi * Bgrad_rho     * xjac * theta * tstep &
                           + (D_par-D_prof) * BigR / BB2              * Bgrad_rho_star     * Bgrad_rho_psi * xjac * theta * tstep &
                           + v * Vpar0 * (r0_s * psi_t - r0_t * psi_s)                                            * theta * tstep &
                           + v * r0 * (vpar0_s * psi_t - vpar0_t * psi_s)                                         * theta * tstep

                 amat_52 =  - v * BigR**2 * ( r0_s * u_t - r0_t * u_s)                                     * theta * tstep &
                            - v * 2.d0 * BigR * r0 * u_y                                            * xjac * theta * tstep

                 amat_55 = v * rho * BigR * xjac * (1.d0 + zeta)  &
                         - v * BigR**2 * ( rho_s * u0_t - rho_t * u0_s)                                       * theta * tstep &
                         - v * 2.d0 * BigR * rho * u0_y                                                * xjac * theta * tstep &
                         + (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star * Bgrad_rho_rho                * xjac * theta * tstep &
                         + D_prof * BigR  * (v_x*rho_x + v_y*rho_y + v_p*rho_p * eps_cyl**2 /BigR**2 ) * xjac * theta * tstep &
                         + v * F0 / BigR * Vpar0 * rho_p                                               * xjac * theta * tstep &
                         + v * Vpar0 * (rho_s * ps0_t - rho_t * ps0_s)                                        * theta * tstep &
                         + v * rho * (vpar0_s * ps0_t - vpar0_t * ps0_s)                                      * theta * tstep &
                         + v * rho * F0 / BigR * vpar0_p                                               * xjac * theta * tstep &

                     + D_perp_num  *                                                                   &
                           ( (v_ss * (x_t(ms,mt)**2 + y_t(ms,mt)**2)                                   &
                            + v_tt * (x_s(ms,mt)**2 + y_s(ms,mt)**2)                                   &
                     -2.d0 *  v_st * (x_s(ms,mt)*x_t(ms,mt) + y_s(ms,mt)*y_t(ms,mt)) )                 &
                           * (rho_ss * (x_t(ms,mt)**2 + y_t(ms,mt)**2)                                 &
                            + rho_tt * (x_s(ms,mt)**2 + y_s(ms,mt)**2)                                 &
                     -2.d0 *  rho_st * (x_s(ms,mt)*x_t(ms,mt) + y_s(ms,mt)*y_t(ms,mt)) ) ) / xjac**3  * theta*tstep


                 amat_57 = + v * F0 / BigR * Vpar * r0_p                                             * xjac * theta * tstep &
                           + v * Vpar * (r0_s * ps0_t - r0_t * ps0_s)                                       * theta * tstep &
                           + v * r0 * (vpar_s * ps0_t - vpar_t * ps0_s)                                     * theta * tstep &
                           + v * r0 * F0 / BigR * vpar_p                                             * xjac * theta * tstep

!---------------------------------------------------------------- equation 6
                 Bgrad_T_star_psi = ( v_x  * psi_y - v_y  * psi_x  ) / BigR
                 Bgrad_T_psi      = ( T0_x * psi_y - T0_y * psi_x )  / BigR
                 Bgrad_T_T        = ( F0 / BigR * T_p +  T_x * ps0_y - T_y * ps0_x ) / BigR          ! F0 due to absence of normalisation

                 amat_61 = - (ZK_par-ZK_prof) * BigR * BB2_psi / BB2**2 * Bgrad_T_star * Bgrad_T     * xjac * theta * tstep &
                           + (ZK_par-ZK_prof) * BigR / BB2     * Bgrad_T_star_psi      * Bgrad_T     * xjac * theta * tstep &
                           + (ZK_par-ZK_prof) * BigR / BB2     * Bgrad_T_star          * Bgrad_T_psi * xjac * theta * tstep &
                           + v * Vpar0 * (T0_s * psi_t - T0_t * psi_s)                                      * theta * tstep &
                           + v * (GAMMA-1.d0) * T0 * (vpar0_s * psi_t - vpar0_t * psi_s)                    * theta * tstep


                 amat_62 = - v * BigR**2 * ( T0_s * u_t - T0_t * u_s)                    * theta * tstep &
                           - v * 2.d0* (GAMMA-1.d0) * BigR * T0 * u_y             * xjac * theta * tstep

                 amat_66 =   v * T   * BigR * xjac * (1.d0 + zeta)    &
                           - v * BigR**2 * ( T_s * u0_t - T_t * u0_s)                        * theta * tstep &
                           - v * 2.d0* (GAMMA-1.d0) * BigR * T * u0_y                 * xjac * theta * tstep &
                           + (ZK_par-ZK_prof) * BigR / BB2 * Bgrad_T_star * Bgrad_T_T * xjac * theta * tstep &
                           + ZK_prof * BigR * (v_x*T_x + v_y*T_y + v_p*T_p /BigR**2 ) * xjac * theta * tstep &
                           + v * F0 / BigR * Vpar0 * T_p                              * xjac * theta * tstep &
                           + v * Vpar0 * (T_s * ps0_t - T_t * ps0_s)                         * theta * tstep &
                           + v * (GAMMA-1.d0) * T * (vpar0_s * ps0_t - vpar0_t * ps0_s)      * theta * tstep &
                           + v * (GAMMA-1.d0) * T * F0 / BigR * vpar0_p               * xjac * theta * tstep

                 amat_67 = + v * F0 / BigR * Vpar * T0_p                              * xjac * theta * tstep &
                           + v * Vpar * (T0_s * ps0_t - T0_t * ps0_s)                        * theta * tstep &
                           + v * (GAMMA-1.d0) * T0 * (vpar_s * ps0_t - vpar_t * ps0_s)       * theta * tstep &
                           + v * (GAMMA-1.d0) * T0 * F0 / BigR * vpar_p               * xjac * theta * tstep

!---------------------------------------------------------------- equation 7


                 amat_71 = + v * (P0_s * psi_t - P0_t * psi_s)                             * theta * tstep ! &
!                           + vpar0 * (F0/BigR)**2 * (vpar0_s * ps_t - vpar0_t * ps_s)      * theta * tstep

                 amat_72 = 0.d0 ! &
!                         + F0 * (u_s * vpar0_t - u_t * vpar0_s)                            * theta * tstep &

                 amat_75 = + v * (rho_s * T0 * ps0_t - rho_t * T0 * ps0_s)                 * theta * tstep &
                           + v * (rho * T0_s * ps0_t - rho * T0_t * ps0_s)                 * theta * tstep &
                           + v * F0 / BigR * rho_p * T0                             * xjac * theta * tstep

                 amat_76 = + v * (T_s * R0 * ps0_t - T_t * R0 * ps0_s)                     * theta * tstep &
                           + v * (T * R0_s * ps0_t - T * R0_t * ps0_s)                     * theta * tstep &
                           + v * F0 / BigR * T_p * R0                               * xjac * theta * tstep

                 amat_77 = v * Vpar * R0 * F0**2 / BigR * xjac * (1.d0 + zeta) &
                         + visco_par * (v_x * Vpar_x + v_y * Vpar_y) * BigR         * xjac  * theta * tstep &
!                         + F0 * (u0_s * vpar_t - u0_t * vpar_s)                             * theta * tstep &
!                         + vpar * (F0/BigR)**2  * (vpar0_s * ps0_t - vpar0_t * ps0_s)       * theta * tstep &
!                         + vpar0 * (F0/BigR)**2 * (vpar_s  * ps0_t - vpar_t * ps0_s)        * theta * tstep &
!                         + (F0/BigR)**3 * vpar  * vpar0_p                            * xjac * theta * tstep &
!                         + (F0/BigR)**3 * vpar0 * vpar_p                             * xjac * theta * tstep

                     + visco_par_num  *                                                                &
                           ( (v_ss * (x_t(ms,mt)**2 + y_t(ms,mt)**2)                                   &
                            + v_tt * (x_s(ms,mt)**2 + y_s(ms,mt)**2)                                   &
                     -2.d0 *  v_st * (x_s(ms,mt)*x_t(ms,mt) + y_s(ms,mt)*y_t(ms,mt)) )                 &
                           * (Vpar_ss * (x_t(ms,mt)**2 + y_t(ms,mt)**2)                                &
                            + Vpar_tt * (x_s(ms,mt)**2 + y_s(ms,mt)**2)                                &
                     -2.d0 *  Vpar_st * (x_s(ms,mt)*x_t(ms,mt) + y_s(ms,mt)*y_t(ms,mt)) ) ) / xjac**3  * theta*tstep



                 kl1 = index_kl
                 kl2 = index_kl + 1*n_tor
                 kl3 = index_kl + 2*n_tor
                 kl4 = index_kl + 3*n_tor
                 kl5 = index_kl + 4*n_tor
                 kl6 = index_kl + 5*n_tor
                 kl7 = index_kl + 6*n_tor

                 ELM(ij1,kl1) =  ELM(ij1,kl1) + wst * amat_11
                 ELM(ij1,kl2) =  ELM(ij1,kl2) + wst * amat_12
                 ELM(ij1,kl3) =  ELM(ij1,kl3) + wst * amat_13
                 ELM(ij1,kl6) =  ELM(ij1,kl6) + wst * amat_16

                 ELM(ij2,kl1) =  ELM(ij2,kl1) + wst * amat_21
                 ELM(ij2,kl2) =  ELM(ij2,kl2) + wst * amat_22
                 ELM(ij2,kl3) =  ELM(ij2,kl3) + wst * amat_23
                 ELM(ij2,kl4) =  ELM(ij2,kl4) + wst * amat_24
                 ELM(ij2,kl5) =  ELM(ij2,kl5) + wst * amat_25
                 ELM(ij2,kl6) =  ELM(ij2,kl6) + wst * amat_26

                 ELM(ij3,kl1) =  ELM(ij3,kl1) + wst * amat_31
                 ELM(ij3,kl3) =  ELM(ij3,kl3) + wst * amat_33

                 ELM(ij4,kl2) =  ELM(ij4,kl2) + wst * amat_42
                 ELM(ij4,kl4) =  ELM(ij4,kl4) + wst * amat_44

                 ELM(ij5,kl1) =  ELM(ij5,kl1) + wst * amat_51
                 ELM(ij5,kl2) =  ELM(ij5,kl2) + wst * amat_52
                 ELM(ij5,kl5) =  ELM(ij5,kl5) + wst * amat_55
                 ELM(ij5,kl7) =  ELM(ij5,kl7) + wst * amat_57

                 ELM(ij6,kl1) =  ELM(ij6,kl1) + wst * amat_61
                 ELM(ij6,kl2) =  ELM(ij6,kl2) + wst * amat_62
                 ELM(ij6,kl6) =  ELM(ij6,kl6) + wst * amat_66
                 ELM(ij6,kl7) =  ELM(ij6,kl7) + wst * amat_67

                 ELM(ij7,kl1) =  ELM(ij7,kl1) + wst * amat_71
                 ELM(ij7,kl2) =  ELM(ij7,kl2) + wst * amat_72
                 ELM(ij7,kl5) =  ELM(ij7,kl5) + wst * amat_75
                 ELM(ij7,kl6) =  ELM(ij7,kl6) + wst * amat_76
                 ELM(ij7,kl7) =  ELM(ij7,kl7) + wst * amat_77

               enddo
             enddo
           enddo

         enddo
       enddo

     enddo
   enddo

 enddo
enddo

return
end subroutine element_matrix
end module mod_elt_matrix
