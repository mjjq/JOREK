module mod_elt_matrix
  implicit none
contains

subroutine element_matrix(element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, ELM, RHS, tid)
!---------------------------------------------------------------
! calculates the matrix contribution of one element
!---------------------------------------------------------------
use constants
use mod_parameters
use data_structure
use gauss
use basis_at_gaussian
use phys_module

implicit none

type (type_element)   :: element
type (type_node)      :: nodes(n_vertex_max)

real*8, dimension (:,:), allocatable  :: ELM
real*8, dimension (:)  , allocatable  :: RHS
integer, intent(in) :: tid

integer    :: i, j, ms, mt, mp, k, l, index_ij, index_kl, index, xcase2
integer    :: in, im, ij1, ij2, ij3, ij4, ij5, ij6, ij7, kl1, kl2, kl3, kl4, kl5, kl6, kl7
real*8     :: wst, xjac, xjac_s, xjac_t, xjac_x, xjac_y, xjac3, BigR, phi
real*8     :: current_source(n_gauss,n_gauss),particle_source(n_gauss,n_gauss),heat_source(n_gauss,n_gauss), source_pellet
real*8     :: R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint(2), Z_xpoint(2), dj_dpsi, dj_dz
real*8     :: UR0, UR0_x, UR0_y, UR0_p, UR0_s, UR0_t, UR0_ss, UR0_st, UR0_tt
real*8     :: UZ0, UZ0_x, UZ0_y, UZ0_p, UZ0_s, UZ0_t, UZ0_ss, UZ0_st, UZ0_tt
real*8     :: UP0, UP0_x, UP0_y, UP0_p, UP0_s, UP0_t, UP0_ss, UP0_st, UP0_tt
real*8     :: AR0, AR0_x, AR0_y, AR0_p, AR0_s, AR0_t, AR0_ss, AR0_st, AR0_tt
real*8     :: AZ0, AZ0_x, AZ0_y, AZ0_p, AZ0_s, AZ0_t, AZ0_ss, AZ0_st, AZ0_tt
real*8     :: AP0, AP0_x, AP0_y, AP0_p, AP0_s, AP0_t, AP0_ss, AP0_st, AP0_tt
real*8     :: r0
real*8     :: T0, T0_x, T0_y, T0_p, T0_s, T0_t, T0_ss, T0_st, T0_tt
real*8     :: P0, P0_x, P0_y, P0_p, P0_s, P0_t, P0_ss, P0_st, P0_tt
real*8     :: UR, UR_x, UR_y, UR_p, UR_s, UR_t, UR_ss, UR_st, UR_tt
real*8     :: UZ, UZ_x, UZ_y, UZ_p, UZ_s, UZ_t, UZ_ss, UZ_st, UZ_tt
real*8     :: UP, UP_x, UP_y, UP_p, UP_s, UP_t, UP_ss, UP_st, UP_tt
real*8     :: AR, AR_x, AR_y, AR_p, AR_s, AR_t, AR_ss, AR_st, AR_tt
real*8     :: AZ, AZ_x, AZ_y, AZ_p, AZ_s, AZ_t, AZ_ss, AZ_st, AZ_tt
real*8     :: AP, AP_x, AP_y, AP_p, AP_s, AP_t, AP_ss, AP_st, AP_tt
real*8     :: T, T_x, T_y, T_p, T_s, T_t, T_ss, T_st, T_tt

real*8     :: v, v_x, v_y, v_s, v_t, v_p, v_ss, v_st, v_tt, v_xx, v_yy, v_xy
real*8     :: bf, bf_x, bf_y, bf_s, bf_t, bf_p, bf_ss, bf_st, bf_tt, bf_xx, bf_yy, bf_xy
   
real*8     :: Bgrad_T_star, Bgrad_T_star_AR, Bgrad_T_star_AZ, Bgrad_T_star_AP
real*8     :: Bgrad_T0,     Bgrad_T0_AR,     Bgrad_T0_AZ,     Bgrad_T0_AP,      Bgrad_T0_T
real*8     :: BB2,          BB2_AR,          BB2_AZ,          BB2_AP
real*8     :: BR0,          BR0_AR,          BR0_AZ,          BR0_AP 
real*8     :: BZ0,          BZ0_AR,          BZ0_AZ,          BZ0_AP 
real*8     :: BP0,          BP0_AR,          BP0_AZ,          BP0_AP 
real*8     :: Bgrad_vstar,  Bgrad_vstar_AR,  Bgrad_vstar_AZ,  Bgrad_vstar_AP
real*8     :: Ugrad_vstar,  Ugrad_vstar_AR,  Ugrad_vstar_AZ,  Ugrad_vstar_AP
real*8     :: BR0_A1star,  BR0_A2star  ,    BR0_A3star
real*8     :: BZ0_A1star,  BZ0_A2star  ,    BZ0_A3star
real*8     :: BP0_A1star,  BP0_A2star  ,    BP0_A3star

real*8     :: ZK_prof, D_prof, psi_norm, theta, zeta, delta_u_x, delta_u_y, delta_ps_x, delta_ps_y


real*8     :: eta_T, visco_T, deta_dT, d2eta_d2T, dvisco_dT, visco_num_T, eta_num_T
real*8     :: amat(n_var,n_var), rhs_ij(n_var)
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

! --- Take time evolution parameters from phys_module
theta = time_evol_theta
zeta  = time_evol_zeta

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

       if (keep_current_prof) &
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
   
   xjac_x  = (x_ss(ms,mt)*y_t(ms,mt)**2 - y_ss(ms,mt)*x_t(ms,mt)*y_t(ms,mt) - 2.d0*x_st(ms,mt)*y_s(ms,mt)*y_t(ms,mt)   &           
	   + y_st(ms,mt)*(x_s(ms,mt)*y_t(ms,mt) + x_t(ms,mt)*y_s(ms,mt)) + x_tt(ms,mt)*y_s(ms,mt)**2 - y_tt(ms,mt)*x_s(ms,mt)*y_s(ms,mt)) / xjac
	   
   xjac_y  = (y_tt(ms,mt)*x_s(ms,mt)**2 - x_tt(ms,mt)*y_s(ms,mt)*x_s(ms,mt) - 2.d0*y_st(ms,mt)*x_t(ms,mt)*x_s(ms,mt)   &           
	   + x_st(ms,mt)*(y_t(ms,mt)*x_s(ms,mt) + y_s(ms,mt)*x_t(ms,mt)) + y_ss(ms,mt)*x_t(ms,mt)**2 - x_ss(ms,mt)*y_t(ms,mt)*x_t(ms,mt)) / xjac

   BigR    = x_g(ms,mt)

   xjac3 = BigR * xjac 

   do mp = 1, n_plane

     uR0   = eq_g(mp,var_u1,ms,mt)
     uR0_x = (   y_t(ms,mt) * eq_s(mp,var_u1,ms,mt) - y_s(ms,mt) * eq_t(mp,var_u1,ms,mt) ) / xjac
     uR0_y = ( - x_t(ms,mt) * eq_s(mp,var_u1,ms,mt) + x_s(ms,mt) * eq_t(mp,var_u1,ms,mt) ) / xjac
     uR0_p = eq_p(mp,var_u1,ms,mt)
     uR0_s = eq_s(mp,var_u1,ms,mt)
     uR0_t = eq_t(mp,var_u1,ms,mt)

     uZ0   = eq_g(mp,var_u2,ms,mt)
     uZ0_x = (   y_t(ms,mt) * eq_s(mp,var_u2,ms,mt) - y_s(ms,mt) * eq_t(mp,var_u2,ms,mt) ) / xjac
     uZ0_y = ( - x_t(ms,mt) * eq_s(mp,var_u2,ms,mt) + x_s(ms,mt) * eq_t(mp,var_u2,ms,mt) ) / xjac
     uZ0_p = eq_p(mp,var_u2,ms,mt)
     uZ0_s = eq_s(mp,var_u2,ms,mt)
     uZ0_t = eq_t(mp,var_u2,ms,mt)

!-------------------------------------------UP0 is defined U0_phi : V = .. + UP0 * e_phi (physical component)
     uP0   = eq_g(mp,var_u3,ms,mt)
     uP0_x = (   y_t(ms,mt) * eq_s(mp,var_u3,ms,mt) - y_s(ms,mt) * eq_t(mp,var_u3,ms,mt) ) / xjac
     uP0_y = ( - x_t(ms,mt) * eq_s(mp,var_u3,ms,mt) + x_s(ms,mt) * eq_t(mp,var_u3,ms,mt) ) / xjac
     uP0_p = eq_p(mp,var_u3,ms,mt)
     uP0_s = eq_s(mp,var_u3,ms,mt)
     uP0_t = eq_t(mp,var_u3,ms,mt)

     AR0   = eq_g(mp,var_A1,ms,mt)
     AR0_x = (   y_t(ms,mt) * eq_s(mp,var_A1,ms,mt) - y_s(ms,mt) * eq_t(mp,var_A1,ms,mt) ) / xjac
     AR0_y = ( - x_t(ms,mt) * eq_s(mp,var_A1,ms,mt) + x_s(ms,mt) * eq_t(mp,var_A1,ms,mt) ) / xjac
     AR0_p = eq_p(mp,var_A1,ms,mt)
     AR0_s = eq_s(mp,var_A1,ms,mt)
     AR0_t = eq_t(mp,var_A1,ms,mt)

     AZ0   = eq_g(mp,var_A2,ms,mt)
     AZ0_x = (   y_t(ms,mt) * eq_s(mp,var_A2,ms,mt) - y_s(ms,mt) * eq_t(mp,var_A2,ms,mt) ) / xjac
     AZ0_y = ( - x_t(ms,mt) * eq_s(mp,var_A2,ms,mt) + x_s(ms,mt) * eq_t(mp,var_A2,ms,mt) ) / xjac
     AZ0_p = eq_p(mp,var_A2,ms,mt)
     AZ0_s = eq_s(mp,var_A2,ms,mt)
     AZ0_t = eq_t(mp,var_A2,ms,mt)
!---------------------------------------------- AP0 == psi
!-------------------------------------------AP0 is defined as A0_3 : A = .. + A0P * grad(phi)
     AP0   = eq_g(mp,var_A3,ms,mt)
     AP0_x = (   y_t(ms,mt) * eq_s(mp,var_A3,ms,mt) - y_s(ms,mt) * eq_t(mp,var_A3,ms,mt) ) / xjac
     AP0_y = ( - x_t(ms,mt) * eq_s(mp,var_A3,ms,mt) + x_s(ms,mt) * eq_t(mp,var_A3,ms,mt) ) / xjac
     AP0_p = eq_p(mp,var_A3,ms,mt)
     AP0_s = eq_s(mp,var_A3,ms,mt)
     AP0_t = eq_t(mp,var_A3,ms,mt)
    
     r0    = 1.d0
     
     T0    = eq_g(mp,var_T,ms,mt)
     T0_x  = (   y_t(ms,mt) * eq_s(mp,var_T,ms,mt) - y_s(ms,mt) * eq_t(mp,var_T,ms,mt) ) / xjac
     T0_y  = ( - x_t(ms,mt) * eq_s(mp,var_T,ms,mt) + x_s(ms,mt) * eq_t(mp,var_T,ms,mt) ) / xjac
     T0_p  = eq_p(mp,var_T,ms,mt)
     T0_s  = eq_s(mp,var_T,ms,mt)
     T0_t  = eq_t(mp,var_T,ms,mt)
     T0_ss = eq_ss(mp,var_T,ms,mt)
     T0_tt = eq_tt(mp,var_T,ms,mt)
     T0_st = eq_st(mp,var_T,ms,mt)

     P0    = abs(r0 * T0)
     P0_x  = r0_x * T0 + r0 * T0_x
     P0_y  = r0_y * T0 + r0 * T0_y
     P0_s  = r0_s * T0 + r0 * T0_s
     P0_t  = r0_t * T0 + r0 * T0_t
     P0_p  = r0_p * T0 + r0 * T0_p
     P0_ss = r0_ss * T0 + 2.d0 * r0_s * T0_s + r0 * T0_ss
     P0_tt = r0_tt * T0 + 2.d0 * r0_t * T0_t + r0 * T0_tt
     P0_st = r0_st * T0 + r0_s * T0_t + r0_t * T0_s + r0 * T0_st

     D_prof  = D_perp(1)  
     ZK_prof = ZK_perp(1)

     phi = 2.d0*PI*float(mp-1)/float(n_plane) / float(n_period)
     
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
	   
! to be done : verify second derivatives           
	   v_xx = (v_ss * y_t(ms,mt)**2 - 2.d0*v_st * y_s(ms,mt)*y_t(ms,mt) + v_tt * y_s(ms,mt)**2  &	        
		+ v_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                          &	   
	        + v_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )  / xjac**2             &		
		- xjac_x * (v_s * y_t(ms,mt) - v_t * y_s(ms,mt)) / xjac**2

	   v_yy = (v_ss * x_t(ms,mt)**2 - 2.d0*v_st * x_s(ms,mt)*x_t(ms,mt) + v_tt * x_s(ms,mt)**2  &	        
		+ v_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                          &	   
	        + v_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )     / xjac**2          &		
		- xjac_y * (- v_s * x_t(ms,mt) + v_t * x_s(ms,mt) ) / xjac**2

           v_xy = (- v_ss * y_t(ms,mt)*x_t(ms,mt) - v_tt * x_s(ms,mt)*y_s(ms,mt)                   &
     	        + v_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                        &        
                - v_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                        &	   
	        - v_t  * (x_st(ms,mt)*y_s(ms,mt) - x_ss(ms,mt)*y_t(ms,mt) )  )  / xjac**2          &		
                - xjac_x * (- v_s * x_t(ms,mt) + v_t * x_s(ms,mt) )   / xjac**2

! to be done replace F0 by F(psi_eq)
           Bgrad_T_star   = ( F0 / BigR * v_p  +  v_x   * (AP0_y - AZ0_p) + v_y  * (AR0_p- AP0_x) &
                                               +  v_p   * (AZ0_x - AR0_y) ) / BigR   
           Bgrad_T0       = ( F0 / BigR * T0_p +  T0_x  * (AP0_y - AZ0_p) + T0_y  * (AR0_p- AP0_x) &
                                               +  T0_p  * (AZ0_x - AR0_y) ) / BigR   

           BR0 = ( AZ0_p - AP0_y )/ BigR
           BZ0 = ( AP0_x - AR0_p )/ BigR
           BP0 = ( AR0_y - AZ0_x )/ BigR

           BB2 = (F0+BP0)**2/BigR**2 + BR0 * BR0 + BZ0 * BZ0 

           Bgrad_vstar = BR0 * v_x + BZ0 * v_y + BP0 * v_p / BigR
           Ugrad_vstar = UR0 * v_x + UZ0 * v_y + UP0 * v_p / BigR

!###################################################################################################
!#  equation 1   (R component momentum equation)                                                   #
!###################################################################################################
 
           rhs_ij(var_u1) = (v_x + v/BigR) * (r0*T0 + 0.5d0*BB2)         * xjac3 * tstep &
	                  + (BR0 * Bgrad_vstar + r0 * UR0 * Ugrad_vstar) * xjac3 * tstep

!###################################################################################################
!#  equation 2   (Z component momentum equation)                                                   #
!###################################################################################################

           rhs_ij(var_u2) =  v_y  * (r0*T0 + 0.5d0*BB2)                   * xjac3 * tstep &
	                  + (BZ0 * Bgrad_vstar + r0 * UZ0 * Ugrad_vstar) * xjac3 * tstep
          
!###################################################################################################
!#  equation 3 (Phi component momentum equation)                                                   #
!###################################################################################################

            rhs_ij(var_u3) = v_p / BigR  * (r0*T0 + 0.5d0*BB2)            * xjac3 * tstep &
	                   + (BP0 * Bgrad_vstar + r0 * UP0 * Ugrad_vstar) * xjac3 * tstep

!###################################################################################################
!#  equation 4 (R component induction equation)                                                    #
!###################################################################################################
           BR0_A1star = 0.d0                      ! rot(A_star) component AR
           BZ0_A1star = - v_p / BigR
           BP0_A1star =   v_y / BigR

           rhs_ij(var_A1) = - eta * (BR0_A1star * BR0 + BZ0_A1star * BZ0 + BP0_A1star * BP0) * xjac3 * tstep &
	                    + v * ( UZ0 * BP0 - UP0 * BZ0) * xjac3 * tstep

!###################################################################################################
!#  equation 5 (Z component induction equation)                                                    #
!###################################################################################################
           BR0_A2star =   v_p / BigR
           BZ0_A2star = 0.d0
           BP0_A2star = - v_x / BigR

           rhs_ij(var_A2) = - eta * (BR0_A2star * BR0 + BZ0_A2star * BZ0 + BP0_A2star * BP0) * xjac3 * tstep &
	                    + v * ( UP0 * BR0 - UR0 * BP0) * xjac3 * tstep
         
!###################################################################################################
!#  equation 6 (PHI component induction equation)                                                  #
!###################################################################################################
           BR0_A3star  =  - v_y / BigR
           BZ0_A3star  =    v_x / BigR
           BP0_A3star  = 0.d0

           rhs_ij(var_A3) = - eta * (BR0_A3star * BR0 + BZ0_A3star * BZ0 + BP0_A3star * BP0) * xjac3 * tstep &
	                    + v * ( UR0 * BZ0 - UZ0 * BR0) * xjac3 * tstep

!###################################################################################################
!#  equation 7 (parallel velocity  equation)                                                       #
!###################################################################################################

           rhs_ij(var_T) = (ZK_par-ZK_perp) * B_grad_T_star * B_grad_T0            * xjac3 * tstep &
                         + ZK_perp * ( v_x * T0_x + v_y * T0_y + v_p * T0_p/ BigR) * xjac3 * tstep &
                         + v * (UR0 * T0_x + UZ0 * T0_y + UP0 * T0_p/BigR )        * xjac3 * tstep &
			 + v * gamma * r0 * T0 * (VR0_x + VR0/BigR + VZ0_y + VP0_p / BigR) * xjac3 * tstep

!###################################################################################################
!#  equations end                                                                                  #
!###################################################################################################

           do ivar=1,n_var
             ij = index_ij + (ivar-1)*n_tor
             RHS(ij) =  RHS(ij) + wst * rhs_ij(ivar)
           enddo

           do k=1,n_vertex_max

             do l=1,n_order+1

               do in = 1, n_tor

                 bf   = H(k,l,ms,mt) * element%size(k,l) * HZ(in,mp)

                 bf_x = (   y_t(ms,mt) * h_s(k,l,ms,mt) - y_s(ms,mt) * h_t(k,l,ms,mt) ) / xjac    &
                              * element%size(k,l) * HZ(in,mp)

                 bf_y = ( - x_t(ms,mt) * h_s(k,l,ms,mt) + x_s(ms,mt) * h_t(k,l,ms,mt) )  / xjac   &
                              * element%size(k,l) * HZ(in,mp)

                 bf_p  = H(k,l,ms,mt)    * element%size(k,l) * HZ_p(in,mp)
                 bf_s  = h_s(k,l,ms,mt)  * element%size(k,l) * HZ(in,mp)
                 bf_t  = h_t(k,l,ms,mt)  * element%size(k,l) * HZ(in,mp)
                 bf_ss = h_ss(k,l,ms,mt) * element%size(k,l) * HZ(in,mp)
                 bf_tt = h_tt(k,l,ms,mt) * element%size(k,l) * HZ(in,mp)
                 bf_st = h_st(k,l,ms,mt) * element%size(k,l) * HZ(in,mp)

                 bf_xx = (bf_ss * y_t(ms,mt)**2 - 2.d0 * bf_st * y_s(ms,mt)*y_t(ms,mt) + bf_tt * y_s(ms,mt)**2  &	        
		        + bf_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                         &	   
	                + bf_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )    / xjac**2             &		
		        - xjac_x * (bf_s * y_t(ms,mt) - bf_t * y_s(ms,mt)) / xjac**2

	         bf_yy = (bf_ss * x_t(ms,mt)**2 - 2.d0 * bf_st * x_s(ms,mt)*x_t(ms,mt) + bf_tt * x_s(ms,mt)**2   &	        
		        + bf_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                          &	   
	                + bf_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )    / xjac**2           &		
		        - xjac_y * (- bf_s * x_t(ms,mt) + bf_t * x_s(ms,mt) ) / xjac**2
           
	         bf_xy = (- bf_ss * y_t(ms,mt)*x_t(ms,mt) - bf_tt * x_s(ms,mt)*y_s(ms,mt)                 &
     	                + bf_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                        &        
                        - bf_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                        &	   
	                - bf_t  * (x_st(ms,mt)*y_s(ms,mt) - x_ss(ms,mt)*y_t(ms,mt) )  )  / xjac**2             &		
                        - xjac_x * (- bf_s * x_t(ms,mt) + bf_t * x_s(ms,mt) )   / xjac**2

                 uR    = bf    ;  uZ    = bf    ;  uP    = bf    ; 
                 uR_x  = bf_x  ;  uZ_x  = bf_x  ;  uP_x  = bf_x  ;
                 uR_y  = bf_y  ;  uZ_y  = bf_y  ;  uP_y  = bf_y  ; 
                 uR_p  = bf_p  ;  uZ_p  = bf_p  ;  uP_p  = bf_p  ; 
                 uR_s  = bf_s  ;  uZ_s  = bf_s  ;  uP_s  = bf_s  ;
                 uR_t  = bf_t  ;  uZ_t  = bf_t  ;  uP_t  = bf_t  ; 
                 uR_ss = bf_ss ;  uZ_ss = bf_ss ;  uP_ss = bf_ss ; 
                 uR_tt = bf_tt ;  uZ_tt = bf_tt ;  uP_tt = bf_tt ; 
                 uR_st = bf_st ;  uZ_st = bf_st ;  uP_st = bf_st ;

                 AR    = bf    ;  AZ    = bf    ;  AP    = bf    ; T    = bf
                 AR_x  = bf_x  ;  AZ_x  = bf_x  ;  AP_x  = bf_x  ; T_x  = bf_x
                 AR_y  = bf_y  ;  AZ_y  = bf_y  ;  AP_y  = bf_y  ; T_y  = bf_y
                 AR_p  = bf_p  ;  AZ_p  = bf_p  ;  AP_p  = bf_p  ; T_p  = bf_p
                 AR_s  = bf_s  ;  AZ_s  = bf_s  ;  AP_s  = bf_s  ; T_s  = bf_s
                 AR_t  = bf_t  ;  AZ_t  = bf_t  ;  AP_t  = bf_t  ; T_t  = bf_T
                 AR_ss = bf_ss ;  AZ_ss = bf_ss ;  AP_ss = bf_ss ; T_ss = bf_ss
                 AR_tt = bf_tt ;  AZ_tt = bf_tt ;  AP_tt = bf_tt ; T_tt = bf_tt
                 AR_st = bf_st ;  AZ_st = bf_st ;  AP_st = bf_st ; T_st = bf_st

		 BR0_AR = 0.d0
                 BZ0_AR = - AR_p / BigR
                 BP0_AR =   AR_y / BigR

	         BR0_AZ =  AZ_p  / BigR
                 BZ0_AZ = 0.d0
                 BP0_AZ = - AZ_x / BigR

                 BR0_AP = - AP_y / BigR
                 BZ0_AP =   AP_x / BigR
                 BP0_AP = 0.d0

		 BB2_AR = 2.d0*(BR0_AR * BR0 + BZ0_AR * BZ0 + BP0_AR * (F0+BP0)/BigR**2 )
                 BB2_AZ = 2.d0*(BR0_AZ * BR0 + BZ0_AZ * BZ0 + BP0_AZ * (F0+BP0)/BigR**2 )
                 BB2_AP = 2.d0*(BR0_AP * BR0 + BZ0_AP * BZ0 + BP0_AP * (F0+BP0)/BigR**2 )
       
		 Bgrad_vstar_AR = BR0_AR * v_x + BZ0_AR * v_y + BP0_AR * v_p / BigR
		 Bgrad_vstar_AZ = BR0_AZ * v_x + BZ0_AZ * v_y + BP0_AZ * v_p / BigR
		 Bgrad_vstar_AP = BR0_AP * v_x + BZ0_AP * v_y + BP0_AP * v_p / BigR
 
                 Ugrad_vstar_UR = UR * v_x 
                 Ugrad_vstar_UZ = UZ * v_y 
                 Ugrad_vstar_UP = UP * v_p / BigR

                 index_kl = n_tor*n_var*(n_order+1)*(k-1) + n_tor * n_var * (l-1) + in   ! index in the ELM matrix

!###################################################################################################
!#  equation 1   (R component momentum equation)                                                   #
!###################################################################################################
   	     
		 amat(var_u1,var_u1) = v * r0 * uR                 * xjac3                 &		                    
	                             - r0 * UR * Ugrad_vstar       * xjac3 * theta * tstep &
				     - r0 * UR0 * Ugrad_vstar_UR   * xjac3 * theta * tstep
		  
                 amat(var_u1,var_u2) = - r0 * UR0 * Ugrad_vstar_UZ * xjac3 * theta * tstep 

                 amat(var_u1,var_u3) = - r0 * UR0 * Ugrad_vstar_UP * xjac3 * theta * tstep 	
		 
                 amat(var_u1,var_a1) = -(v_x + v/BigR) * ( 0.5d0*BB2_AR)              * xjac3 * theta * tstep &
		                     - (BR0_AR * Bgrad_vstar + BR0 * Bgrad_vstar_AR ) * xjac3 * theta * tstep 

		 amat(var_u1,var_a2) = -(v_x + v/BigR) * ( 0.5d0*BB2_AZ)              * xjac3 * theta * tstep &
		                     - (BR0_AZ * Bgrad_vstar + BR0 * Bgrad_vstar_AZ ) * xjac3 * theta * tstep 
                 
		 amat(var_u1,var_a3) = -(v_x + v/BigR) * ( 0.5d0*BB2_AP)              * xjac3 * theta * tstep &
		                     - (BR0_AP * Bgrad_vstar + BR0 * Bgrad_vstar_AP ) * xjac3 * theta * tstep 

                 amat(var_u1,var_T)  = - (v_x + v/BigR) * r0*T                        * xjac3 * theta * tstep 

!###################################################################################################
!#  equation 2   (Z component momentum equation)                                                   #
!###################################################################################################
 
                 amat(var_u2,var_u1) = - r0 * UZ0 * Ugrad_vstar_UR * xjac3 * theta * tstep

                 amat(var_u2,var_u2) = v * r0 * uZ               * xjac3                 &	
		                     - r0 * UZ  * Ugrad_vstar    * xjac3 * theta * tstep &
				     - r0 * UZ0 * Ugrad_vstar_UZ * xjac3 * theta * tstep
 
                 amat(var_u2,var_u3) = - r0 * UZ0 * Ugrad_vstar_UP * xjac3 * theta * tstep

                 amat(var_u2,var_a1) = - v_y  * 0.5d0*BB2_AR                         * xjac3 * theta * tstep &
		                     - (BZ0_AR * Bgrad_vstar + BZ0 * Bgrad_vstar_AR) * xjac3 * theta *tstep

		 amat(var_u2,var_a2) = - v_y  * 0.5d0*BB2_AZ                         * xjac3 * theta * tstep &
		                     - (BZ0_AZ * Bgrad_vstar + BZ0 * Bgrad_vstar_AZ) * xjac3 * theta *tstep

		 amat(var_u2,var_a3) = - v_y  * 0.5d0*BB2_AP                         * xjac3 * theta * tstep &
		                     - (BZ0_AP * Bgrad_vstar + BZ0 * Bgrad_vstar_AP) * xjac3 * theta *tstep

                 amat(var_u2,var_T) = - v_y  * r0*T                                  * xjac3 * theta * tstep 

!###################################################################################################
!#  equation 3   (Phi component momentum equation)                                                 #
!###################################################################################################
           rhs_ij(var_u3) = v_p / BigR  * (r0*T0 + 0.5d0*BB2)            * xjac3 * tstep &
	                   + (BP0 * Bgrad_vstar + r0 * UP0 * Ugrad_vstar) * xjac3 * tstep

                 amat(var_u3,var_u1) = - r0 * UP0 * Ugrad_vstar_UR * xjac3 * theta * tstep

		 amat(var_u3,var_u2) = - r0 * UP0 * Ugrad_vstar_UZ * xjac3 * theta * tstep

		 amat(var_u3,var_u3) = v * r0 * uP / BigR        * xjac3                 &
		                     - r0 * UP  * Ugrad_vstar    * xjac3 * theta * tstep &
		                     - r0 * UP0 * Ugrad_vstar_UP * xjac3 * theta * tstep

                 amat(var_u3,var_a1) = - v_p / BigR  * 0.5d0*BB2_AR                   * xjac3 * theta * tstep &
		                     - ( BP0_AR * Bgrad_vstar + BP0 * Bgrad_vstar_AR) * xjac3 * theta * tstep

                amat(var_u3,var_a2) = - v_p / BigR  * 0.5d0*BB2_AZ                    * xjac3 * theta * tstep &
		                     - ( BP0_AZ * Bgrad_vstar + BP0 * Bgrad_vstar_AZ) * xjac3 * theta * tstep

                amat(var_u3,var_a3) = - v_p / BigR  * 0.5d0*BB2_AP                    * xjac3 * theta * tstep &
		                     - ( BP0_AP * Bgrad_vstar + BP0 * Bgrad_vstar_AP) * xjac3 * theta * tstep

                amat(var_u3,var_T) = - v_p / BigR  * r0*T                             * xjac3 * theta * tstep 
		 
!###################################################################################################
!#  equation 4   (R component induction equation)                                                  #
!###################################################################################################
            

                amat(var_A1,var_u1) = 0.d0
	        amat(var_A1,var_u2) = - v * ( UZ * BP0 ) * xjac3 * theta * tstep
	        amat(var_A1,var_u3) = + v * ( UP * BZ0 ) * xjac3 * theta * tstep                    
	    
           amat(var_A1,var_A1) = v * AR * xjac3  &
	                       + eta * (BR0_AR_star * BR0_AR + BZ0_AR_star * BZ0_AR + BP0_AR_star * BP0_AR) * xjac3 * theta * tstep &
	                       - v * ( UZ0 * BP0_AR - UP0 * BZ0_AR)                                         * xjac3 * theta * tstep

           amat(var_A1,var_A2) = eta * (BR0_AR_star * BR0_AZ + BZ0_AR_star * BZ0_AZ + BP0_AR_star * BP0_AZ) * xjac3 * theta * tstep &
	                       - v * ( UZ0 * BP0_AZ - UP0 * BZ0_AZ)                                         * xjac3 * theta *tstep
           
	   amat(var_A1,var_A3) = eta * (BR0_AR_star * BR0_AP + BZ0_AR_star * BZ0_AP + BP0_AR_star * BP0_AP) * xjac3 * theta * tstep &
	                       - v * ( UZ0 * BP0_AP - UP0 * BZ0_AP)                                         * xjac3 * theta * tstep

 	   amat(var_A1,var_T) = 0.d0

!###################################################################################################
!#  equation 5   (Z component induction equation)                                                  #
!###################################################################################################
           
           amat(var_A2,var_u1) =  + v * UR * BP0 * xjac3 * theta * tstep

	   amat(var_A2,var_u2) = 0.d0

           amat(var_A2,var_u3) = - v *  UP * BR0 * xjac3 * theta * tstep

	   amat(var_A2,var_A1) = eta * (BR0_AZ_star * BR0_AR + BZ0_AZ_star * BZ0_AR + BP0_AZ_star * BP0_AR) * xjac3 * theta * tstep &
	                       - v * ( UP0 * BR0_AR - UR0 * BP0_AR)                                         * xjac3 * theta * tstep

	   amat(var_A2,var_A2) = v * AZ * xjac3  &
                               + eta * (BR0_AZ_star * BR0_AZ + BZ0_AZ_star * BZ0_AZ + BP0_AZ_star * BP0_AZ) * xjac3 * theta *tstep &
	                       - v * ( UP0 * BR0_AZ - UR0 * BP0_AZ)                                         * xjac3 * theta *tstep
 
	   amat(var_A2,var_A3) = eta * (BR0_AZ_star * BR0_AP + BZ0_AZ_star * BZ0_AP + BP0_AZ_star * BP0_AP) * xjac3 * theta *tstep &
	                       - v * ( UP0 * BR0_AP - UR0 * BP0_AP)                                         * xjac3 * theta * tstep

           amat(var_A2,var_T) = 0.d0

!###################################################################################################
!#  equation 6   (Phi component induction equation)                                                #
!###################################################################################################

            amat(var_A3,var_u1) = - v * UR * BZ0 * xjac3 * theta * tstep

            amat(var_A3,var_u2) = + v * UZ * BR0 * xjac3 * theta * tstep

	    amat(var_A3,var_u3) = 0.d0

            amat(var_A3,var_A1) = eta * (BR0_AP_star * BR0_AR + BZ0_AP_star * BZ0_AR + BP0_AP_star * BP0_AR) * xjac3 * theta * tstep &
	                        - v * ( UR0 * BZ0_AR - UZ0 * BR0_AR)                                         * xjac3 * theta * tstep

	    amat(var_A3,var_A2) = eta * (BR0_AP_star * BR0_AZ + BZ0_AP_star * BZ0_AZ + BP0_AP_star * BP0_AZ) * xjac3 * theta * tstep &
	                        - v * ( UR0 * BZ0_AZ - UZ0 * BR0_AZ)                                         * xjac3 * theta * tstep

	    amat(var_A3,var_A3) =  v * AP * xjac3  &
	                        + eta * (BR0_AP_star * BR0_AP + BZ0_AP_star * BZ0_AP + BP0_AP_star * BP0_AP) * xjac3 * theta * tstep &
	                        - v * ( UR0 * BZ0_AP - UZ0 * BR0_AP)                                         * xjac3 * theta * tstep

!###################################################################################################
!#  equation 7   parallel velocity equation                                                        #
!###################################################################################################

	   Bgrad_T_star   = ( F0 / BigR * v_p  +  v_x   * (AP0_y - AZ0_p) + v_y  * (AR0_p- AP0_x) &
                                               +  v_p   * (AZ0_x - AR0_y) ) / BigR   
           Bgrad_T0        = ( F0 / BigR * T0_p +  T0_x  * (AP0_y - AZ0_p) + T0_y  * (AR0_p- AP0_x) &
                                               +  T0_p  * (AZ0_x - AR0_y) ) / BigR   

	   Bgrad_T_star_AR = (  v_y  * AR_p - v_p * AR_y ) / BigR   
	   Bgrad_T_star_AZ = (- v_x  * AZ_p + v_p * AZ_x ) / BigR
	   Bgrad_T_star_AP = (  v_x  * AP_y - v_y * AP_x ) / BigR   

	   Bgrad_T0_AR     = (  T0_y * AR_p - T0_p * AR_y ) / BigR   
	   Bgrad_T0_AZ     = (- T0_x * AZ_p + T0_p * AZ_x ) / BigR
	   Bgrad_T0_AP     = (  T0_x * AP_y - T0_y * AP_x ) / BigR   

           Bgrad_T0_T     = ( F0 / BigR * T_p +  T_x  * (AP0_y - AZ0_p) + T_y  * (AR0_p- AP0_x) &
                                             +  T_p  * (AZ0_x - AR0_y) ) / BigR   

           amat(var_T,var_u1) =  - v * UR * T0_x                          * xjac3 * theta * tstep &
			         - v * gamma * r0 * T0 * (VR_x + VR/BigR) * xjac3 * theta * tstep

           amat(var_T,var_u2) = - v * UZ * T0_y                           * xjac3 * theta * tstep &
			        - v * gamma * r0 * T0 * VZ_y              * xjac3 * theta * tstep

           amat(var_T,var_u3) = - v * UP * T0_p/BigR                      * xjac3 * theta * tstep &
			        - v * gamma * r0 * T0 * VP_p / BigR       * xjac3 * theta * tstep

           amat(var_T,var_a1) = - (ZK_par-ZK_perp) * B_grad_T_star_AR * B_grad_T0    * xjac3 * theta * tstep &
	                        - (ZK_par-ZK_perp) * B_grad_T_star    * B_grad_T0_AR * xjac3 * theta * tstep 
                            
           amat(var_T,var_a2) = - (ZK_par-ZK_perp) * B_grad_T_star_AZ * B_grad_T0    * xjac3 * theta * tstep &
	                        - (ZK_par-ZK_perp) * B_grad_T_star    * B_grad_T0_AZ * xjac3 * theta * tstep 

           amat(var_T,var_a3) = - (ZK_par-ZK_perp) * B_grad_T_star_AP * B_grad_T0    * xjac3 * theta * tstep &
	                        - (ZK_par-ZK_perp) * B_grad_T_star    * B_grad_T0_AP * xjac3 * theta * tstep 

           amat(var_T,var_T)  =  v * r0 * T *  xjac3 * (1.d0 + zeta)    &
                              - (ZK_par-ZK_perp) * B_grad_T_star * B_grad_T0_T                 * xjac3 * theta * tstep &
                              - ZK_perp * ( v_x * T_x + v_y * T_y + v_p * T_p/ BigR)           * xjac3 * theta * tstep &
                              - v * (UR0 * T_x + UZ0 * T_y + UP0 * T_p/BigR )                  * xjac3 * theta * tstep &
			      - v * gamma * r0 * T * (VR0_x + VR0/BigR + VZ0_y + VP0_p / BigR) * xjac3 * theta * tstep

!###################################################################################################



                 kl1 = index_kl
                 kl2 = index_kl + 1*n_tor
                 kl3 = index_kl + 2*n_tor
                 kl4 = index_kl + 3*n_tor
                 kl5 = index_kl + 4*n_tor
                 kl6 = index_kl + 5*n_tor
                 kl7 = index_kl + 6*n_tor

                 do ivar=1,n_var
	           do kvar=1, n_var

	             ij = index_ij + (ivar-1)*n_tor
		     kl = index_kl + (kvar-1)*n_tor

                     ELM(ij,kl) =  ELM(ij,kl) + wst * amat(ivar,kvar)
		   enddo
		 enddo

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
