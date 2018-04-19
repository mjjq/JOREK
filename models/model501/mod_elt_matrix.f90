module mod_elt_matrix

  implicit none

contains

subroutine element_matrix(element, nodes, xpoint2, xcase2, minRad, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, ELM,RHS, tid)
!---------------------------------------------------------------
! calculates the matrix contribution of one element
!---------------------------------------------------------------
use constants
use mod_parameters
use data_structure
use gauss
use basis_at_gaussian
use phys_module
use pellet_module
use diffusivities, only: get_dperp, get_zkperp
use corr_neg
use mgi_module
use mod_coronal
use vacuum, only: freeb_fact

implicit none

type (type_element)   :: element
type (type_node)      :: nodes(n_vertex_max)

real*8, dimension (:,:), pointer  :: ELM
real*8, dimension (:)  , pointer  :: RHS
integer, intent(in) :: tid

integer    :: i, j, ms, mt, mp, k, l, index_ij, index_kl, index, xcase2
integer    :: in, im, ij1, ij2, ij3, ij4, ij5, ij6, ij7, kl1, kl2, kl3, kl4, kl5, kl6, kl7
real*8     :: wst, xjac, xjac_s, xjac_t, xjac_x, xjac_y, BigR, r2, phi, delta_phi, eps_cyl
real*8     :: current_source(n_gauss,n_gauss), particle_source(n_gauss,n_gauss), heat_source(n_gauss,n_gauss)
real*8     :: source_volume, source_pellet, source_pellet2
real*8     :: minRad, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint(2), Z_xpoint(2), dj_dpsi, dj_dz
real*8     :: Bgrad_rho_star,     Bgrad_rho,     Bgrad_rhon, Bgrad_T_star, Bgrad_T, BB2
real*8     :: Bgrad_rho_star_psi, Bgrad_rho_psi, Bgrad_rhon_psi, Bgrad_rho_rho, Bgrad_rho_rhon, Bgrad_T_star_psi, Bgrad_T_psi, Bgrad_T_T, BB2_psi
real*8     :: rhs_ij_1,   rhs_ij_2,   rhs_ij_3,   rhs_ij_4,   rhs_ij_5,   rhs_ij_6, rhs_ij_7
real*8     :: rhs_stab_1, rhs_stab_2, rhs_stab_3, rhs_stab_4, rhs_stab_5, rhs_stab_6    
real*8     :: ZK_prof, D_prof, psi_norm, theta, zeta, delta_u_x, delta_u_y, delta_ps_x, delta_ps_y, ZKpar_T, dZKpar_dT
real*8     :: v, v_x, v_y, v_s, v_t, v_p, v_ss, v_st, v_tt, v_xx, v_yy, v_xy
real*8     :: ps0, ps0_x, ps0_y, ps0_p,ps0_s,ps0_t, ps0_ss, ps0_st, ps0_tt, ps0_xx, ps0_xy, ps0_yy
real*8     :: zj0, zj0_x, zj0_y, zj0_p, zj0_s, zj0_t
real*8     :: u0, u0_x, u0_y, u0_p, u0_s, u0_t, u0_ss, u0_st, u0_tt, u0_xx, u0_xy, u0_yy 
real*8     :: w0, w0_x, w0_y, w0_p, w0_s, w0_t, w0_ss, w0_st, w0_tt, w0_xx, w0_xy, w0_yy
real*8     :: Vpar0, Vpar0_x, Vpar0_y, Vpar0_p, Vpar0_s, Vpar0_t, Vpar0_ss, Vpar0_st, Vpar0_tt, vpar0_xx, vpar0_xy, vpar0_yy
real*8     :: r0, r0_x, r0_y, r0_p, r0_s, r0_t, r0_ss, r0_st, r0_tt, r0_xx, r0_xy, r0_yy
real*8     :: r0_corr, rn0_corr, r0_hat, r0_x_hat, r0_y_hat, T0_corr, dT0_corr_dT, d2T0_corr_dT2
real*8     :: T0, T0_x, T0_y, T0_p, T0_s, T0_t, T0_ss, T0_st, T0_tt, T0_xx, T0_xy, T0_yy, T_corr
real*8     :: psi, psi_x, psi_y, psi_p, psi_s, psi_t, psi_ss, psi_st, psi_tt, psi_xx, psi_yy, psi_xy
real*8     :: zj, zj_x, zj_y, zj_p, zj_s, zj_t, zj_ss, zj_st, zj_tt
real*8     :: vpar, vpar_x, vpar_y, vpar_s, vpar_t, vpar_p, vpar_ss, vpar_st, vpar_tt, vpar_xx, vpar_xy, vpar_yy
real*8     :: u, u_x, u_y, u_p, u_s, u_t, u_ss, u_st, u_tt, u_xx, u_xy, u_yy
real*8     :: w, w_x, w_y, w_p, w_s, w_t, w_ss, w_st, w_tt, w_xx, w_xy, w_yy
real*8     :: rho, rho_x, rho_y, rho_s, rho_t, rho_p, rho_ss, rho_st, rho_tt, rho_xx, rho_xy, rho_yy, rho_hat, rho_x_hat, rho_y_hat
real*8     :: T, T_x, T_y, T_s, T_t, T_p, T_ss, T_st, T_tt, T_xx, T_xy, T_yy
real*8     :: P0, P0_x, P0_y, P0_s, P0_t, P0_ss, P0_st, P0_tt, P0_p, P0_xx, P0_xy, P0_yy
real*8     :: BigR_x, vv2, eta_T, visco_T, deta_dT, d2eta_d2T, dvisco_dT, d2visco_d2T
real*8     :: visco_num_T, eta_num_T, eta_Sp, detaSp_dT
real*8     :: amat_11, amat_12, amat_21, amat_22, amat_23, amat_24, amat_25, amat_26, amat_33, amat_31, amat_44, amat_42
real*8     :: amat_51, amat_52, amat_55, amat_56, amat_57, amat_61, amat_62, amat_63, amat_65, amat_66, amat_67, amat_16, amat_13
real*8     :: amat_71, amat_72, amat_75, amat_76, amat_77, amat_15
real*8     :: ZK_par_num, T0_ps0_x, T_ps0_x, T0_psi_x, T0_ps0_y, T_ps0_y, T0_psi_y, v_ps0_x, v_psi_x, v_ps0_y, v_psi_y
real*8     :: TG_num1, TG_num2, TG_num5, TG_num6, TG_num7, TG_num8
logical    :: xpoint2

!==================MB: velocity profile is kept by a source which compensating diffusion
real*8     :: Vt0,Vt0_x,Vt0_y
real*8     :: V_source(n_gauss,n_gauss)
real*8     :: dV_dpsi_source(n_gauss,n_gauss),dV_dz_source(n_gauss,n_gauss)
real*8     :: dV_dpsi2,dV_dz2,dV_dpsi_dz,dV_dpsi3,dV_dpsi_dz2,dV_dpsi2_dz
!=======================================
real*8     :: eq_zne(n_gauss,n_gauss), eq_zTe(n_gauss,n_gauss)
real*8     :: dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2,dn_dpsi2_dz
real*8     :: dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2,dT_dpsi2_dz
real*8     :: w00_xx, w00_yy 
!======================================= NEO
real*8     :: amat_27, Btheta2
real*8     :: epsil, Btheta2_psi
real*8, dimension(n_gauss,n_gauss)    :: amu_neo_prof, aki_neo_prof
!======================================= NEO

!================== Parameters specific to model5XX
! Matrix, RHS and neutrals-related variables
real*8     :: amat_28, amat_58, amat_68, amat_78
real*8     :: amat_81, amat_82, amat_85, amat_86, amat_87, amat_88
real*8     :: rhs_ij_8
real*8     :: ij8, kl8
real*8     :: rn0, rn0_x, rn0_y, rn0_p, rn0_s, rn0_t, rn0_ss, rn0_st, rn0_tt, rn0_hat, rn0_x_hat, rn0_y_hat
real*8     :: rhon, rhon_x, rhon_y, rhon_s, rhon_t, rhon_p, rhon_ss, rhon_st, rhon_tt, rhon_hat, rhon_x_hat, rhon_y_hat
real*8     :: rn0_xx, rn0_yy, rhon_xx, rhon_yy

! Neutral source
real*8     :: source_mgi
real*8     :: source_mgi_tmp


! time normalization
real*8     :: t_norm

! Temporary variables serving the SPI module
integer    :: spi_i

real*8     :: spi_R_tmp
real*8     :: spi_Z_tmp
real*8     :: spi_phi_tmp
real*8     :: spi_abl_tmp
real*8     :: ng_radius !< Radius of neutral gas cloud as a result of the ablation
! Additional variables reserved for future implementation
!real*8     :: spi_Vel_R_tmp
!real*8     :: spi_Vel_Z_tmp
!real*8     :: spi_Vel_phi_tmp

! Neutral diffusion coefficients
real*8     :: Dn0x, Dn0y, Dn0p

! Atomic physics coefficients:
!   -Mass ratio between main ions and impurites (m_i/m_imp)
real*8     :: m_i_over_m_imp
!   -Mean impurity ionization state
real*8     :: Z_imp, dZ_imp_dT, T0_Zimp, alpha_Zimp
!   -Coefficients related to Z_imp
real*8     :: alpha_imp, dalpha_imp_dT, alpha_imp_bis
real*8     :: beta_imp, dbeta_imp_dT

!   -Radiation from injected impurities
real*8     :: Lrad, dLrad_dT                                  ! Radiation rate and its derivative wrt. temperature
real*8     :: T_rad, dT_rad_dT                                ! Temperature used in radiation rate
real*8     :: ne_rad                                          ! Electron density used in radiation rate
real*8     :: coef_rad_1, A0_rad, A1_rad, T1_rad, sig1_rad    ! Radiation rate parameters
real*8     :: A2_rad, T2_rad, sig2_rad

!   -Radiation from background impurities
real*8     :: Arad_bg, Brad_bg, Crad_bg, frad_bg, dfrad_bg_dT

!   -Temporary variable for charge state distribution
real*8, allocatable :: dP_imp_dT(:), P_imp(:)
real*8     :: E_ion, dE_ion_dT
integer*8  :: ion_i, ion_k

! Parameters related to negative temperature handling
real*8     :: T_neg, delta_neg
!===============================

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

zk_par_num = 0.d0

! --- Taylor-Galerkin Stabilisation coefficients
TG_num1    = TGNUM(1); TG_num2    = TGNUM(2); TG_num5    = TGNUM(5); TG_num6    = TGNUM(6); TG_num7    = TGNUM(7);

TG_num8    = TGNUM(8)

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
V_source=0.d0
dV_dpsi_source=0.d0
dV_dz_source=0.d0
eq_zne          = 0.d0
eq_zTe          = 0.d0         



do i=1,n_vertex_max
 do j=1,n_order+1

   do ms=1, n_gauss
     do mt=1, n_gauss

       x_g(ms,mt)  = x_g(ms,mt)  + nodes(i)%x(j,1) * element%size(i,j) * H(i,j,ms,mt)
       x_s(ms,mt)  = x_s(ms,mt)  + nodes(i)%x(j,1) * element%size(i,j) * H_s(i,j,ms,mt)
       x_t(ms,mt)  = x_t(ms,mt)  + nodes(i)%x(j,1) * element%size(i,j) * H_t(i,j,ms,mt)

       x_ss(ms,mt) = x_ss(ms,mt) + nodes(i)%x(j,1) * element%size(i,j) * H_ss(i,j,ms,mt)
       x_st(ms,mt) = x_st(ms,mt) + nodes(i)%x(j,1) * element%size(i,j) * H_st(i,j,ms,mt)
       x_tt(ms,mt) = x_tt(ms,mt) + nodes(i)%x(j,1) * element%size(i,j) * H_tt(i,j,ms,mt)

       y_g(ms,mt)  = y_g(ms,mt)  + nodes(i)%x(j,2) * element%size(i,j) * H(i,j,ms,mt)
       y_s(ms,mt)  = y_s(ms,mt)  + nodes(i)%x(j,2) * element%size(i,j) * H_s(i,j,ms,mt)
       y_t(ms,mt)  = y_t(ms,mt)  + nodes(i)%x(j,2) * element%size(i,j) * H_t(i,j,ms,mt)

       y_ss(ms,mt) = y_ss(ms,mt) + nodes(i)%x(j,2) * element%size(i,j) * H_ss(i,j,ms,mt)
       y_st(ms,mt) = y_st(ms,mt) + nodes(i)%x(j,2) * element%size(i,j) * H_st(i,j,ms,mt)
       y_tt(ms,mt) = y_tt(ms,mt) + nodes(i)%x(j,2) * element%size(i,j) * H_tt(i,j,ms,mt)

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
     enddo
   enddo
 enddo
enddo

do ms=1, n_gauss
  do mt=1, n_gauss

       call current(xpoint2, xcase2, x_g(ms,mt),y_g(ms,mt), Z_xpoint, eq_g(1,1,ms,mt),psi_axis,psi_bnd,current_source(ms,mt))

       call sources(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, eq_g(1,1,ms,mt),psi_axis,psi_bnd,particle_source(ms,mt),heat_source(ms,mt))
      
       current_source(ms,mt) = 0.d0

       call density(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, eq_g(1,1,ms,mt),psi_axis,psi_bnd,eq_zne(ms,mt), &
                    dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz)

       call temperature(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, eq_g(1,1,ms,mt),psi_axis,psi_bnd,eq_zTe(ms,mt), &
                        dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2, dT_dpsi2_dz)
  enddo
enddo

eq_zTe = eq_zTe / 2.d0  ! electron temperature

!--------------------------------------------------- sum over the Gaussian integration points
do ms=1, n_gauss

 do mt=1, n_gauss

   wst = wgauss(ms)*wgauss(mt)

   xjac    = x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)

   xjac_x  = (x_ss(ms,mt)*y_t(ms,mt)**2 - y_ss(ms,mt)*x_t(ms,mt)*y_t(ms,mt) - 2.d0*x_st(ms,mt)*y_s(ms,mt)*y_t(ms,mt)   &
           + y_st(ms,mt)*(x_s(ms,mt)*y_t(ms,mt) + x_t(ms,mt)*y_s(ms,mt))                                               &
           + x_tt(ms,mt)*y_s(ms,mt)**2 - y_tt(ms,mt)*x_s(ms,mt)*y_s(ms,mt)) / xjac

   xjac_y  = (y_tt(ms,mt)*x_s(ms,mt)**2 - x_tt(ms,mt)*y_s(ms,mt)*x_s(ms,mt) - 2.d0*y_st(ms,mt)*x_t(ms,mt)*x_s(ms,mt)   &
	   + x_st(ms,mt)*(y_t(ms,mt)*x_s(ms,mt) + y_s(ms,mt)*x_t(ms,mt))                                               &
           + y_ss(ms,mt)*x_t(ms,mt)**2 - x_ss(ms,mt)*y_t(ms,mt)*x_t(ms,mt)) / xjac

   BigR    = x_g(ms,mt)
   BigR_x  = 1.d0

   eps_cyl = 1.d0          ! for cylinder geometry : epscyl = eps

   do mp = 1, n_plane

     ps0    = eq_g(mp,1,ms,mt)
     ps0_x  = (   y_t(ms,mt) * eq_s(mp,1,ms,mt) - y_s(ms,mt) * eq_t(mp,1,ms,mt) ) / xjac
     ps0_y  = ( - x_t(ms,mt) * eq_s(mp,1,ms,mt) + x_s(ms,mt) * eq_t(mp,1,ms,mt) ) / xjac
     ps0_p  = eq_p(mp,1,ms,mt)
     ps0_s  = eq_s(mp,1,ms,mt)
     ps0_t  = eq_t(mp,1,ms,mt)
     ps0_ss = eq_ss(mp,1,ms,mt)
     ps0_tt = eq_tt(mp,1,ms,mt)
     ps0_st = eq_st(mp,1,ms,mt)

     u0    = eq_g(mp,2,ms,mt)
     u0_x  = (   y_t(ms,mt) * eq_s(mp,2,ms,mt) - y_s(ms,mt) * eq_t(mp,2,ms,mt) ) / xjac
     u0_y  = ( - x_t(ms,mt) * eq_s(mp,2,ms,mt) + x_s(ms,mt) * eq_t(mp,2,ms,mt) ) / xjac
     u0_p  = eq_p(mp,2,ms,mt)
     u0_s  = eq_s(mp,2,ms,mt)
     u0_t  = eq_t(mp,2,ms,mt)
     u0_ss = eq_ss(mp,2,ms,mt)
     u0_tt = eq_tt(mp,2,ms,mt)
     u0_st = eq_st(mp,2,ms,mt)

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

     r0_corr = corr_neg_dens(r0,(/1.d-8,1.d-5/)) ! Correction for negative r0 ...

     r0_hat   = BigR**2 * r0
     r0_x_hat = 2.d0 * BigR * BigR_x  * r0 + BigR**2 * r0_x
     r0_y_hat = BigR**2 * r0_y

     rn0    = eq_g(mp,8,ms,mt)
     rn0_x  = (   y_t(ms,mt) * eq_s(mp,8,ms,mt) - y_s(ms,mt) * eq_t(mp,8,ms,mt) ) / xjac    
     rn0_y  = ( - x_t(ms,mt) * eq_s(mp,8,ms,mt) + x_s(ms,mt) * eq_t(mp,8,ms,mt) ) / xjac   
     rn0_p  = eq_p(mp,8,ms,mt)                                                             
     rn0_s  = eq_s(mp,8,ms,mt)                                                             
     rn0_t  = eq_t(mp,8,ms,mt)                                                             
     rn0_ss = eq_ss(mp,8,ms,mt)                                                            
     rn0_st = eq_st(mp,8,ms,mt)                                                            
     rn0_tt = eq_tt(mp,8,ms,mt)                                                            

     rn0_corr = corr_neg_dens(rn0, (/ 1.d-12, 1.d-5 /)) ! Correction for negative rn0 ...

     rn0_xx = (rn0_ss * y_t(ms,mt)**2 - 2.d0*rn0_st * y_s(ms,mt)*y_t(ms,mt) + rn0_tt * y_s(ms,mt)**2     &
            + rn0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                                 &
            + rn0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )    /    xjac**2               &
            - xjac_x * (rn0_s* y_t(ms,mt) - rn0_t * y_s(ms,mt))  / xjac**2

     rn0_yy = (rn0_ss * x_t(ms,mt)**2 - 2.d0*rn0_st * x_s(ms,mt)*x_t(ms,mt) + rn0_tt * x_s(ms,mt)**2     &
            + rn0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                                 &
            + rn0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )          / xjac**2            &
            - xjac_y * (- rn0_s * x_t(ms,mt) + rn0_t * x_s(ms,mt) )  / xjac**2



     rn0_hat   = BigR**2 * rn0                                                        
     rn0_x_hat = 2.d0 * BigR * BigR_x  * rn0 + BigR**2 * rn0_x                             
     rn0_y_hat = BigR**2 * rn0_y                                                            

     T0    = eq_g(mp,6,ms,mt)
     T0_x  = (   y_t(ms,mt) * eq_s(mp,6,ms,mt) - y_s(ms,mt) * eq_t(mp,6,ms,mt) ) / xjac
     T0_y  = ( - x_t(ms,mt) * eq_s(mp,6,ms,mt) + x_s(ms,mt) * eq_t(mp,6,ms,mt) ) / xjac
     T0_p  = eq_p(mp,6,ms,mt)
     T0_s  = eq_s(mp,6,ms,mt)
     T0_t  = eq_t(mp,6,ms,mt)
     T0_ss = eq_ss(mp,6,ms,mt)
     T0_tt = eq_tt(mp,6,ms,mt)
     T0_st = eq_st(mp,6,ms,mt)

     T0_corr = corr_neg_temp(T0,(/1.d-2,1.d-1/)) ! For use in eta(T), visco(T), ...
     dT0_corr_dT = dcorr_neg_temp_dT(T0,(/1.d-2,1.d-1/)) ! Improve the correction
     d2T0_corr_dT2 = d2corr_neg_temp_dT2(T0,(/1.d-2,1.d-1/))

     Vpar0    = eq_g(mp,7,ms,mt)
     Vpar0_x  = (   y_t(ms,mt) * eq_s(mp,7,ms,mt) - y_s(ms,mt) * eq_t(mp,7,ms,mt) ) / xjac
     Vpar0_y  = ( - x_t(ms,mt) * eq_s(mp,7,ms,mt) + x_s(ms,mt) * eq_t(mp,7,ms,mt) ) / xjac
     Vpar0_p  = eq_p(mp,7,ms,mt)
     Vpar0_s  = eq_s(mp,7,ms,mt)
     Vpar0_t  = eq_t(mp,7,ms,mt)
     Vpar0_ss = eq_ss(mp,7,ms,mt)
     Vpar0_st = eq_st(mp,7,ms,mt)
     Vpar0_tt = eq_tt(mp,7,ms,mt)

     P0    = r0 * T0
     P0_x  = r0_x * T0 + r0 * T0_x
     P0_y  = r0_y * T0 + r0 * T0_y
     P0_s  = r0_s * T0 + r0 * T0_s
     P0_t  = r0_t * T0 + r0 * T0_t
     P0_p  = r0_p * T0 + r0 * T0_p
     P0_ss = r0_ss * T0 + 2.d0 * r0_s * T0_s + r0 * T0_ss
     P0_tt = r0_tt * T0 + 2.d0 * r0_t * T0_t + r0 * T0_tt
     P0_st = r0_st * T0 + r0_s * T0_t + r0_t * T0_s + r0 * T0_st

     ps0_xx = (ps0_ss * y_t(ms,mt)**2 - 2.d0*ps0_st * y_s(ms,mt)*y_t(ms,mt) + ps0_tt * y_s(ms,mt)**2 &
             + ps0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                            &
             + ps0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )    / xjac**2             &
             - xjac_x * (ps0_s* y_t(ms,mt) - ps0_t * y_s(ms,mt))  / xjac**2

     ps0_yy = (ps0_ss * x_t(ms,mt)**2 - 2.d0*ps0_st * x_s(ms,mt)*x_t(ms,mt) + ps0_tt * x_s(ms,mt)**2 &
             + ps0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                            &
             + ps0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )    / xjac**2             &
             - xjac_y * (- ps0_s * x_t(ms,mt) + ps0_t * x_s(ms,mt) )  / xjac**2

     ps0_xy = (- ps0_ss * y_t(ms,mt)*x_t(ms,mt) - ps0_tt * x_s(ms,mt)*y_s(ms,mt)                     &
              + ps0_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                          &
              - ps0_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                          &
              - ps0_t * (x_st(ms,mt)*y_s(ms,mt)  - x_ss(ms,mt)*y_t(ms,mt) ) )  / xjac**2             &
              - xjac_x * (- ps0_s * x_t(ms,mt) + ps0_t * x_s(ms,mt) )   / xjac**2

     u0_xx = (u0_ss * y_t(ms,mt)**2 - 2.d0*u0_st * y_s(ms,mt)*y_t(ms,mt) + u0_tt * y_s(ms,mt)**2     &
           + u0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                               &
           + u0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )      / xjac**2              &
           - xjac_x * (u0_s * y_t(ms,mt) - u0_t * y_s(ms,mt)) / xjac**2

     u0_yy = (u0_ss * x_t(ms,mt)**2 - 2.d0*u0_st * x_s(ms,mt)*x_t(ms,mt) + u0_tt * x_s(ms,mt)**2     &
           + u0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                               &
           + u0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )      / xjac**2              &
           - xjac_y * (- u0_s * x_t(ms,mt) + u0_t * x_s(ms,mt) ) / xjac**2

     u0_xy = (- u0_ss * y_t(ms,mt)*x_t(ms,mt) - u0_tt * x_s(ms,mt)*y_s(ms,mt)                        &
              + u0_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                           &
              - u0_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                           &
              - u0_t * (x_st(ms,mt)*y_s(ms,mt)  - x_ss(ms,mt)*y_t(ms,mt) )  )  / xjac**2             &
              - xjac_x * (- u0_s * x_t(ms,mt) + u0_t * x_s(ms,mt) )   / xjac**2

!     w0_xx = (w0_ss * y_t(ms,mt)**2 - 2.d0*w0_st * y_s(ms,mt)*y_t(ms,mt) + w0_tt * y_s(ms,mt)**2     &
!            + w0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                              &
!            + w0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )     / xjac**2              &
!            - xjac_x * (w0_s* y_t(ms,mt) - w0_t * y_s(ms,mt))  / xjac**2

!     w0_yy = (w0_ss * x_t(ms,mt)**2 - 2.d0*w0_st * x_s(ms,mt)*x_t(ms,mt) + w0_tt * x_s(ms,mt)**2     &
!            + w0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                              &
!            + w0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )     / xjac**2              &
!            - xjac_y * (- w0_s * x_t(ms,mt) + w0_t * x_s(ms,mt) )  / xjac**2

     w0_xy = (- w0_ss * y_t(ms,mt)*x_t(ms,mt) - w0_tt * x_s(ms,mt)*y_s(ms,mt)                        &
              + w0_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                           &
              - w0_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                           &
              - w0_t  * (x_st(ms,mt)*y_s(ms,mt) - x_ss(ms,mt)*y_t(ms,mt) ) )  / xjac**2              &
              - xjac_x * (- w0_s * x_t(ms,mt) + w0_t * x_s(ms,mt) )   / xjac**2

!----------------- simplified version of 2nd derivatives (for some unknown reason this is more stable!)
     w0_xx = (w0_ss * y_t(ms,mt)**2 - 2.d0*w0_st * y_s(ms,mt)*y_t(ms,mt) + w0_tt * y_s(ms,mt)**2 )     / xjac**2
     w0_yy = (w0_ss * x_t(ms,mt)**2 - 2.d0*w0_st * x_s(ms,mt)*x_t(ms,mt) + w0_tt * x_s(ms,mt)**2 )     / xjac**2

     r0_xx = (r0_ss * y_t(ms,mt)**2 - 2.d0*r0_st * y_s(ms,mt)*y_t(ms,mt) + r0_tt * y_s(ms,mt)**2     &
	    + r0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                              &
	    + r0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )    / xjac**2               &
            - xjac_x * (r0_s* y_t(ms,mt) - r0_t * y_s(ms,mt))  / xjac**2

     r0_yy = (r0_ss * x_t(ms,mt)**2 - 2.d0*r0_st * x_s(ms,mt)*x_t(ms,mt) + r0_tt * x_s(ms,mt)**2     &
            + r0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                              &
	    + r0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )       / xjac**2            &	
            - xjac_y * (- r0_s * x_t(ms,mt) + r0_t * x_s(ms,mt) )  / xjac**2

     r0_xy = (- r0_ss * y_t(ms,mt)*x_t(ms,mt) - r0_tt * x_s(ms,mt)*y_s(ms,mt)                        &
     	      + r0_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                           &
              - r0_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                           &
	      - r0_t * (x_st(ms,mt)*y_s(ms,mt)  - x_ss(ms,mt)*y_t(ms,mt) ) )  / xjac**2              &	
              - xjac_x * (- r0_s * x_t(ms,mt) + r0_t * x_s(ms,mt) )   / xjac**2

     T0_xx = (T0_ss * y_t(ms,mt)**2 - 2.d0*T0_st * y_s(ms,mt)*y_t(ms,mt) + T0_tt * y_s(ms,mt)**2     &
	    + T0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                              &
	    + T0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )    / xjac**2               &	
            - xjac_x * (T0_s * y_t(ms,mt) - T0_t * y_s(ms,mt))  / xjac**2

     T0_yy = (T0_ss * x_t(ms,mt)**2 - 2.d0*T0_st * x_s(ms,mt)*x_t(ms,mt) + T0_tt * x_s(ms,mt)**2     &
            + T0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                              &
	    + T0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )       / xjac**2            &	
            - xjac_y * (- T0_s * x_t(ms,mt) + T0_t * x_s(ms,mt) )  / xjac**2

     T0_xy = (- T0_ss * y_t(ms,mt)*x_t(ms,mt) - T0_tt * x_s(ms,mt)*y_s(ms,mt)                        &
     	      + T0_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                           &
              - T0_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                           &
	      - T0_t  * (x_st(ms,mt)*y_s(ms,mt) - x_ss(ms,mt)*y_t(ms,mt) )  )  / xjac**2             &
              - xjac_x * (- T0_s * x_t(ms,mt) + T0_t * x_s(ms,mt) )   / xjac**2

     vpar0_xx = (vpar0_ss * y_t(ms,mt)**2 - 2.d0*vpar0_st * y_s(ms,mt)*y_t(ms,mt) + vpar0_tt * y_s(ms,mt)**2 &
               + vpar0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                                &
               + vpar0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )   / xjac**2                  &
               - xjac_x * (vpar0_s * y_t(ms,mt) - vpar0_t * y_s(ms,mt)) / xjac**2

     vpar0_yy = (vpar0_ss * x_t(ms,mt)**2 - 2.d0*vpar0_st * x_s(ms,mt)*x_t(ms,mt) + vpar0_tt * x_s(ms,mt)**2 &
               + vpar0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                                &
               + vpar0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )      / xjac**2               &
               - xjac_y * (- vpar0_s * x_t(ms,mt) + vpar0_t * x_s(ms,mt) ) / xjac**2

     vpar0_xy = (- vpar0_ss * y_t(ms,mt)*x_t(ms,mt) - vpar0_tt * x_s(ms,mt)*y_s(ms,mt)                       &
                 + vpar0_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                             &
                 - vpar0_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                             &
                 - vpar0_t * (x_st(ms,mt)*y_s(ms,mt)  - x_ss(ms,mt)*y_t(ms,mt) )  )  / xjac**2               &
                 - xjac_x * (- vpar0_s * x_t(ms,mt) + vpar0_t * x_s(ms,mt) )   / xjac**2


     P0_xx = r0_xx * T0 + 2.d0 * r0_x * T0_x + r0 * T0_xx
     P0_yy = r0_yy * T0 + 2.d0 * r0_y * T0_y + r0 * T0_yy
     P0_xy = r0_xy * T0 + r0_x * T0_y + r0_y * T0_x + r0 * T0_xy
 
     T0_ps0_x = T0_xx * ps0_y - T0_xy * ps0_x + T0_x * ps0_xy - T0_y * ps0_xx
     T0_ps0_y = T0_xy * ps0_y - T0_yy * ps0_x + T0_x * ps0_yy - T0_y * ps0_xy

     delta_u_x = (   y_t(ms,mt) * delta_s(mp,2,ms,mt) - y_s(ms,mt) * delta_t(mp,2,ms,mt) ) / xjac
     delta_u_y = ( - x_t(ms,mt) * delta_s(mp,2,ms,mt) + x_s(ms,mt) * delta_t(mp,2,ms,mt) ) / xjac

     delta_ps_x = (   y_t(ms,mt) * delta_s(mp,1,ms,mt) - y_s(ms,mt) * delta_t(mp,1,ms,mt) ) / xjac
     delta_ps_y = ( - x_t(ms,mt) * delta_s(mp,1,ms,mt) + x_s(ms,mt) * delta_t(mp,1,ms,mt) ) / xjac
     
     ! --- Temperature dependent resistivity
     if ( eta_T_dependent ) then
       eta_T     = eta   * (T0_corr/T_0)**(-1.5d0)
       deta_dT   = ( - eta   * (1.5d0)  * T0_corr**(-2.5d0) * T_0**(1.5d0) ) * dT0_corr_dT
       d2eta_d2T = (  eta   * (3.75d0) * T0_corr**(-3.5d0) * T_0**(1.5d0) ) * d2T0_corr_dT2
     else
       eta_T     = eta
       deta_dT   = 0.d0
       d2eta_d2T = 0.d0
     end if

     !eta_Sp = 1.65d-9*17*(1.d-3*T0_corr/(2*EL_CHG*MU_ZERO*central_density*1.d20))**(-1.5d0) &
     !                          *sqrt(central_mass*MASS_PROTON*1.d20*central_density/MU_ZERO) 

     eta_Sp = 0. ! For intear benchmark
     
     detaSp_dT = -1.65d-9*17 * (1.5d0) * T0_corr**(-2.5d0) * (1.d-3/(2*EL_CHG*MU_ZERO*central_density*1.d20))**(-1.5d0) &
                      * sqrt(central_mass*MASS_PROTON*1.d20*central_density/MU_ZERO) * dT0_corr_dT

     ! --- Temperature dependent viscosity
     if ( visco_T_dependent ) then       
       visco_T   = visco * (T0_corr/T_0)**(-1.5d0)
       dvisco_dT = - visco * (1.5d0)  * T0_corr**(-2.5d0) * T_0**(1.5d0) * dT0_corr_dT
     else
       visco_T   = visco
       dvisco_dT = 0.d0
     end if
     
     ! --- Temperature dependent parallel heat diffusivity
     if ( ZKpar_T_dependent ) then
       ZKpar_T   = ZK_par * (T0_corr/T_0)**(+2.5d0)              ! temperature dependent parallel conductivity
       dZKpar_dT = ZK_par * (2.5d0)  * T0_corr**(+1.5d0) * T_0**(-2.5d0) * dT0_corr_dT
       if (ZKpar_T .gt. ZK_par_max) then
         ZKpar_T   = Zk_par_max
         dZKpar_dT = 0.d0
       endif
     else
       ZKpar_T   = ZK_par                                            ! parallel conductivity
       dZKpar_dT = 0.d0
     endif

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

     D_prof  = get_dperp (psi_norm)
     ZK_prof = get_zkperp(psi_norm)

!========= For intear benchmark=====
     !T_neg = 1.d-5
     !delta_neg = 1.d-4

     !ZK_prof = ZK_prof * (1+100*(0.5+0.5*tanh(-(T0-T_neg)/delta_neg))) 
!============End==============

     !if (T0 .lt. T_neg) then
     !heat_source = 1.d-4
     !endif    

     ! --- Increase diffusivity if very small density/temperature
     if (xpoint2) then
       if (r0 .lt. D_prof_neg_thresh)  then
         D_prof  = D_prof_neg
       endif
       if (T0 .lt. ZK_prof_neg_thresh) then
         ZK_prof = ZK_prof_neg
       endif
     endif

     phi       = 2.d0*PI*float(mp-1)/float(n_plane) / float(n_period)
     delta_phi = 2.d0*PI/float(n_plane) / float(n_period)

     source_pellet = 0.d0
     source_volume = 0.d0
   
     Dn0x = D_neutral_x      
     Dn0y = D_neutral_y      
     Dn0p = D_neutral_p      

  !-------------------------------------------
  ! Atomic physics parameters for Argon
  !-------------------------------------------

     select case ( trim(gas_type) )
       case('D2')
         m_i_over_m_imp = 1.
       case('Ar')
         m_i_over_m_imp = 1./20. ! Argon mass = 40 u and main ion (D) mass = 2 u
       case('Ne')
         m_i_over_m_imp = 1./10. ! Neon mass = 20 u and main ion (D) mass = 2 u
       case default
         write(*,*) '!! Gas type "', trim(gas_type), '" unknown (in mgi_source.f90) !!'
         write(*,*) '=> We assume the gas is D2.'
         m_i_over_m_imp = 1.
     end select


     ! Te in eV:
     T_rad = T0_corr/(2.d0*EL_CHG*MU_ZERO*central_density*1.d20)
     dT_rad_dT = dT0_corr_dT/(2.d0*EL_CHG*MU_ZERO*central_density*1.d20)

     ! We estimate the effective charge by a test density 10^20/m^3
     ! Later maybe we should implement a iterative method
     if (flag_adas) then

       if (allocated(imp_adas(1)%ionisation_energy)) then

         if (allocated(P_imp)) deallocate(P_imp)
         if (allocated(dP_imp_dT)) deallocate(dP_imp_dT)

         allocate(P_imp(1:imp_adas(1)%n_Z))
         allocate(dP_imp_dT(1:imp_adas(1)%n_Z))

         call imp_cor(1)%interp(density=20.,temperature=log10(T_rad*EL_CHG/K_BOLTZ),&
                                p_out=P_imp,z_eff=Z_imp)
         call imp_cor(1)%interp_gradients(density=20.,temperature=log10(T_rad*EL_CHG/K_BOLTZ),&
                                          p_Te_out=dP_imp_dT,z_eff_Te=dZ_imp_dT)

         ! Calculate the ionization potential energy and it's time gradient
         E_ion     = 0.
         dE_ion_dT = 0.

         do ion_i=1, imp_adas(1)%n_Z
           do ion_k=1, ion_i
             E_ion     = E_ion + P_imp(ion_i)*imp_adas(1)%ionisation_energy(ion_k)
             dE_ion_dT = dE_ion_dT + dP_imp_dT(ion_i)*imp_adas(1)%ionisation_energy(ion_k)
           end do
         end do
         ! Convert from eV to JOREK unit
         E_ion     = E_ion * EL_CHG*MU_ZERO*central_density*1.d20
         dE_ion_dT = dE_ion_dT * EL_CHG*MU_ZERO*central_density*1.d20 
         ! Convert the gradient in K to gradient in JOREK unit
         dE_ion_dT = dE_ion_dT * dT_rad_dT * EL_CHG / K_BOLTZ

       else
         call imp_cor(1)%interp(density=20.,temperature=log10(T_rad*EL_CHG/K_BOLTZ),z_eff=Z_imp)
         call imp_cor(1)%interp_gradients(density=20.,temperature=log10(T_rad*EL_CHG/K_BOLTZ),z_eff_Te=dZ_imp_dT)
         E_ion     = 0.
         dE_ion_dT = 0.
       end if

       ! Convert gradient in T(K) in to gradient in T (eV)
       dZ_imp_dT = dZ_imp_dT *EL_CHG / K_BOLTZ
       ! Derivative wrt to T, with T in JOREK units
       dZ_imp_dT = dZ_imp_dT / (2.d0*EL_CHG*MU_ZERO*central_density*1.d20)
       dZ_imp_dT = dZ_imp_dT * dT0_corr_dT
     else

       T0_Zimp        = 437.  ! eV
       alpha_Zimp     = 0.415

       Z_imp     = 10. !18.*tanh((T_rad/T0_Zimp)**alpha_Zimp)
     ! Derivative wrt to Te, with Te in eV
       dZ_imp_dT = 0. !(18./T0_Zimp)*alpha_Zimp*((T_rad/T0_Zimp)**(alpha_Zimp-1))*(1.-(tanh(T_rad/T0_Zimp))**2.) * dT_rad_dT
     ! Derivative wrt to T, with T in JOREK units
       dZ_imp_dT = dZ_imp_dT / (2.d0*EL_CHG*MU_ZERO*central_density*1.d20)

       E_ion     = 0.
       dE_ion_dT = 0.

     end if

     alpha_imp     = 0.5*m_i_over_m_imp*(Z_imp+1.) - 1.
     dalpha_imp_dT = 0.5*m_i_over_m_imp*dZ_imp_dT
     alpha_imp_bis = alpha_imp + dalpha_imp_dT*T0

     beta_imp     = m_i_over_m_imp*Z_imp - 1.
     dbeta_imp_dT = m_i_over_m_imp*dZ_imp_dT

     ne_rad       = (r0_corr + beta_imp * rn0_corr) * 1.d20 * central_density ! electron density (SI)

  !-------------------------------------------
  ! --- Radiative function, if flag_adas is enabled use interpolation, if not use simple model
  ! ------------------------------------------

     ! Normalization coefficient for radiation rate from SI units (W.m^3) to JOREK units:
     coef_rad_1 = 2.d0/3.d0*MU_ZERO**1.5d0*(central_mass*MASS_PROTON)**0.5d0&
                  *(central_density*1.d20)**2.5d0*m_i_over_m_imp


     if (flag_adas .and. ne_rad > 1.d16 .and. T_rad > 1.) then

       Lrad = 0.0
       dLrad_dT = 0.0

       ! Here we are temperarily only considering one impurity species, in the
       ! future maybe a do loop will is needed
       call radiation_function(imp_adas(1),imp_cor(1),log10(ne_rad),log10(T_rad*EL_CHG/K_BOLTZ),Lrad,dLrad_dT)

       Lrad = Lrad * coef_rad_1

       ! Convert gradient in T(K) in to gradient in T (eV)
       dLrad_dT = dLrad_dT * coef_rad_1 * EL_CHG / K_BOLTZ
       ! Derivative wrt to T, with T in JOREK units
       dLrad_dT = dLrad_dT / (2.d0*EL_CHG*MU_ZERO*central_density*1.d20)
       dLrad_dT = dLrad_dT * dT0_corr_dT

     else if (flag_adas) then
       Lrad = 0.
       dLrad_dT = 0.
       E_ion = 0.
       dE_ion_dT = 0.
     else

  !-------------------------------------------
  ! --- Radiative cooling rate for Argon (approximate fit of cooling rate at coronal equilibrium)
  ! ------------------------------------------
    
!   if (T_rad .gt. 5.) then

     A0_rad   = 2.8*1.d-33    ! W.m^3
     A1_rad   = 2.335*1.d-31  ! W.m^3
     T1_rad   = 23.           ! eV
     sig1_rad = 14.           ! eV
     A2_rad   = 3.846*1.d-32  ! W.m^3
     T2_rad   = 236.          ! eV
     sig2_rad = 150.          ! eV

!     Lrad     = coef_rad_1*(A0_rad + A1_rad*exp(-((T_rad-T1_rad)/sig1_rad)**4.) + A2_rad*exp(-((T_rad-T2_rad)/sig2_rad)**2))
     !Lrad     = (1./2.)*coef_rad_1*5.d-32 * (tanh((T_rad-20.)/10.)-tanh(-20./10.))
     Lrad      = 0. ! For Test

     ! Derivative wrt to T, with T in JOREK units
     !dLrad_dT = (1./2.)*coef_rad_1*5.d-32 * (1./10.) * dT_rad_dT * (1-tanh((T_rad-20.)/10.)**2) * dT0_corr_dT
     dLrad_dT = 0. ! For Test

!   else
!     Lrad = 0.d0
!     dLrad_dT = 0.d0
!   endif
     end if

     ! This is to detect N/A
     if (Lrad/=Lrad .or. dLrad_dT/=dLrad_dT .or. E_ion/=E_ion .or. dE_ion_dT/=dE_ion_dT) then
       write(*,*) "WARNING: Lrad, dLrad_dT, E_ion/=E_ion, dE_ion_dT/=dE_ion_dT = ",&
                            Lrad, dLrad_dT, E_ion, dE_ion_dT
       stop
     end if


   !--------------------------------------------------------
   ! --- Source of neutrals from Massive Gas Injection (MGI)
   !--------------------------------------------------------

     source_mgi = 0.d0                    

!============================================================!
! Important note: in order to implementing more complicated  !
!    model, we should add more arguments to mgi_source       !
!============================================================!

     if (using_spi) then

       if (JET_MGI .or. ASDEX_MGI) then
         write(*,*) "WARNING: Using SPI, disabling MGI settings"
         JET_MGI = .false.
         ASDEX_MGI = .false.
       end if

       do spi_i=1, n_spi

         source_mgi_tmp = 0.d0

         if (pellets(spi_i)%spi_radius > 0.0) then
           spi_R_tmp   = pellets(spi_i)%spi_R
           spi_Z_tmp   = pellets(spi_i)%spi_Z
           spi_phi_tmp = pellets(spi_i)%spi_phi
           spi_abl_tmp = pellets(spi_i)%spi_abl

           ng_radius   = pellets(spi_i)%spi_radius * ng_radius_ratio

           if (ng_radius < ng_radius_min) then
             ng_radius = ng_radius_min
           end if

           call mgi_source(spi_abl_tmp,spi_R_tmp,spi_Z_tmp,spi_phi_tmp,ng_radius,mgi_sig,mgi_deltaphi,&
                         mgi_tor_norm, A_Dmv,K_Dmv,V_Dmv,P_Dmv,t_mgi,L_tube,x_g(ms,mt),y_g(ms,mt),    &
                         phi,source_mgi_tmp,t_now,JET_MGI,ASDEX_MGI,central_density,central_mass)
         end if

         source_mgi = source_mgi + source_mgi_tmp

       end do

     else

       call mgi_source(mgi_amplitude,mgi_R,mgi_Z,mgi_phi,mgi_radius,mgi_sig,mgi_deltaphi,mgi_tor_norm, &
                     A_Dmv,K_Dmv,V_Dmv,P_Dmv,t_mgi,L_tube,x_g(ms,mt),y_g(ms,mt),phi,source_mgi,t_now,  &
                     JET_MGI,ASDEX_MGI,central_density,central_mass)

     end if

     ! This is to detect N/A
     if (source_mgi /= source_mgi) then
       write(*,*) "WARNING: source_mgi = ", source_mgi
       stop
     end if     

     if (source_mgi .lt. 0.d0) then
      source_mgi = 0.d0
     endif


   !--------------------------------------------------------
   ! --- Radiation from background impurity
   !--------------------------------------------------------

    Arad_bg = 2.4d-31
    Brad_bg = 20.
    Crad_bg = 0.8

    frad_bg     = (2./3.)*(1./(central_mass*MASS_PROTON))*((MU_ZERO*central_mass*MASS_PROTON*central_density*1.d20)**(1.5d0))                &
                  *nimp_bg*Arad_bg*exp(-((log(T_rad)-log(Brad_bg))**2.)/Crad_bg**2.)

    dfrad_bg_dT = -(1./3.)*((MU_ZERO*central_mass*MASS_PROTON*central_density*1.d20)**(0.5d0))*(1./EL_CHG)                                   &
                  *2.*(nimp_bg*Arad_bg/Crad_bg**2.)*(log(T_rad)-log(Brad_bg))*(1./T_rad)*exp(-((log(T_rad)-log(Brad_bg))**2.)/Crad_bg**2.)

!--------------------------------------------------------

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
	   
	   v_xx = (v_ss * y_t(ms,mt)**2 - 2.d0*v_st * y_s(ms,mt)*y_t(ms,mt) + v_tt * y_s(ms,mt)**2  &	        
		+ v_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                          &	   
	        + v_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )  / xjac**2             &		
		- xjac_x * (v_s * y_t(ms,mt) - v_t * y_s(ms,mt)) / xjac**2

	   v_yy = (v_ss * x_t(ms,mt)**2 - 2.d0*v_st * x_s(ms,mt)*x_t(ms,mt) + v_tt * x_s(ms,mt)**2  &	        
		+ v_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                          &	   
	        + v_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )     / xjac**2          &		
		- xjac_y * (- v_s * x_t(ms,mt) + v_t * x_s(ms,mt) ) / xjac**2

           v_xy = (- v_ss * y_t(ms,mt)*x_t(ms,mt) - v_tt * x_s(ms,mt)*y_s(ms,mt)                    &
     	        + v_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                         &        
                - v_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                         &   
	        - v_t  * (x_st(ms,mt)*y_s(ms,mt) - x_ss(ms,mt)*y_t(ms,mt) )  )  / xjac**2           &
                - xjac_x * (- v_s * x_t(ms,mt) + v_t * x_s(ms,mt) )   / xjac**2

           Bgrad_rho_star = ( F0 / BigR * v_p  +  v_x  * ps0_y - v_y  * ps0_x ) / BigR    ! F0 due to absence of normalisation
           Bgrad_rho      = ( F0 / BigR * r0_p +  r0_x * ps0_y - r0_y * ps0_x ) / BigR    ! F0 due to absence of normalisation
           Bgrad_rhon     = ( F0 / BigR * rn0_p +  rn0_x * ps0_y - rn0_y * ps0_x ) / BigR    ! F0 due to absence of normalisation
           Bgrad_T_star   = ( F0 / BigR * v_p  +  v_x  * ps0_y - v_y  * ps0_x ) / BigR    ! F0 due to absence of normalisation
           Bgrad_T        = ( F0 / BigR * T0_p +  T0_x * ps0_y - T0_y * ps0_x ) / BigR    ! F0 due to absence of normalisation

           BB2            = (F0*F0 + ps0_x * ps0_x + ps0_y * ps0_y )/BigR**2
           Btheta2        = (ps0_x * ps0_x + ps0_y * ps0_y )/BigR
   !        Btheta2_psi  = 2.d0 * (psi_x * ps0_x + psi_y * ps0_y ) /BigR**2
           
           v_ps0_x  = v_xx  * ps0_y - v_xy  * ps0_x + v_x  * ps0_xy - v_y * ps0_xx
           v_ps0_y  = v_xy  * ps0_y - v_yy  * ps0_x + v_x  * ps0_yy - v_y * ps0_xy

!###################################################################################################
!#  equation 1   (induction equation)                                                              #
!###################################################################################################

              rhs_ij_1 =   v * eta_T  * (zj0 - current_source(ms,mt))/ BigR  * xjac * tstep &
                      + v * (ps0_s * u0_t - ps0_t * u0_s)                        * tstep &
                      - v * eps_cyl * F0 / BigR  * u0_p                   * xjac * tstep &
                      + eta_num_T * (v_x * zj0_x + v_y * zj0_y)           * xjac * tstep &

                      - v * tauIC/(r0_corr*BB2) * F0**2/BigR**2 * (ps0_s * p0_t - ps0_t * p0_s) * tstep &
                      + v * tauIC/(r0_corr*BB2) * F0**3/BigR**3 * eps_cyl * p0_p * xjac * tstep &

                      + zeta * v * delta_g(mp,1,ms,mt) / BigR             * xjac


!###################################################################################################
!#  equation 2   (perpendicular momentum equation)                                                 #
!###################################################################################################

         rhs_ij_2 = - 0.5d0 * vv2 * (v_x * r0_y_hat - v_y * r0_x_hat)     * xjac * tstep &
                      - r0_hat * BigR**2 * w0 * (v_s * u0_t - v_t * u0_s)        * tstep &
                      + v * (ps0_s * zj0_t - ps0_t * zj0_s )                     * tstep &
                      - visco_T * BigR * (v_x * w0_x + v_y * w0_y)        * xjac * tstep &
                      - v * eps_cyl * F0 / BigR * zj0_p                   * xjac * tstep &
                      + BigR**2 * (v_s * p0_t - v_t * p0_s)                      * tstep &

                      - visco_num_T * (v_xx + v_x/Bigr + v_yy)*(w0_xx + w0_x/Bigr + w0_yy) * xjac * tstep &

                      - TG_num2 * 0.25d0 * r0_hat * BigR**3 * (w0_x * u0_y - w0_y * u0_x) &
                                * ( v_x * u0_y - v_y * u0_x) * xjac * tstep * tstep &

!====================================New TG_num terms=================================
                      - TG_num2 * 0.25d0 * w0 * BigR**3 * BigR**2 * (r0_x * u0_y - r0_y * u0_x) &
                                * ( v_x * u0_y - v_y * u0_x) * xjac * tstep * tstep &
!===============================End of NewTG_num terms==============================


                      - v * tauIC * BigR**4 * (p0_s * w0_t - p0_t * w0_s)        * tstep &

                      - tauIC * BigR**3 * p0_y * (v_x* u0_x + v_y * u0_y) * xjac * tstep &

                      - v * tauIC * BigR**4 * (u0_xy * (p0_xx - p0_yy) - p0_xy * (u0_xx - u0_yy) ) * xjac * tstep &

                      !+ BigR**3 * (particle_source(ms,mt) + source_pellet) * (v_x * u0_x + v_y * u0_y) * xjac* tstep &

                      - zeta * BigR * r0_hat * (v_x * delta_u_x + v_y * delta_u_y) * xjac &

                      - zeta * BigR * BigR**2 * delta_g(mp,5,ms,mt) * (v_x * u0_x + v_y * u0_y) * xjac &
                      ! Additional zeta term for timesteping here (new momentum)
                      - BigR**2 * (r0_x_hat * u0_y - r0_y_hat * u0_x) * (v_x * u0_x + v_y * u0_y) * xjac* tstep &
                      ! Perp component of the third term of Eq.20 here (new momentum)
                      + BigR**2 * r0 * F0 / BigR * vpar0_p * (v_x * u0_x + v_y * u0_y) * xjac* tstep &
                      + BigR**2 * vpar0 * F0 / BigR * r0_p * (v_x * u0_x + v_y * u0_y) * xjac* tstep &
                      ! Toroidal para component of the third term of Eq.20 here (new momentum)
                      + BigR**2 * r0 * (vpar0_x * ps0_y - vpar0_y * ps0_x) &
                        * (v_x * u0_x + v_y * u0_y) * xjac* tstep  &
                      + BigR**2 * vpar0 * (r0_x * ps0_y - r0_y * ps0_x) &
                        * (v_x * u0_x + v_y * u0_y) * xjac* tstep
                      ! Poloidal para component of the third term of Eq.20 here (new momentum)


!###################################################################################################
!#  equation 3   (current definition)                                                              #
!###################################################################################################

         rhs_ij_3 = - ( v_x * ps0_x  + v_y * ps0_y + v*zj0) / BigR * xjac * freeb_fact

!###################################################################################################
!#  equation 4   (vorticity definition)                                                            #
!###################################################################################################

         rhs_ij_4 = 0.d0 !- ( v_x * u0_x   + v_y * u0_y  + v*w0)  * BigR * xjac 

!###################################################################################################
!#  equation 5 (total density equation)                                                                  #
!###################################################################################################

         rhs_ij_5   = v * BigR * (particle_source(ms,mt) + source_pellet + source_mgi)                                 * xjac * tstep &
                    + v * BigR**2 * ( r0_s * u0_t - r0_t * u0_s)                                                              * tstep &
                    + v * 2.d0 * BigR * r0 * u0_y                                                                      * xjac * tstep &
                    !- (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star * (Bgrad_rho-Bgrad_rhon)                            * xjac * tstep &
                    ! The new diffusion scheme for the impurities
                    - (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star * (Bgrad_rho)                                       * xjac * tstep &
                    !- D_prof * BigR  * (v_x*(r0_x-rn0_x) + v_y*(r0_y-rn0_y) + v_p*(r0_p-rn0_p) * eps_cyl**2 /BigR**2 ) * xjac * tstep &
                    ! The new diffusion scheme for the impurities
                    - D_prof * BigR  * (v_x*(r0_x) + v_y*(r0_y) + v_p*(r0_p) * eps_cyl**2 /BigR**2 )                   * xjac * tstep &
                    - v * F0 / BigR * Vpar0 * r0_p                                                                     * xjac * tstep &
                    - v * Vpar0 * (r0_s * ps0_t - r0_t * ps0_s)                                                               * tstep &
                    - v * F0 / BigR * r0 * vpar0_p                                                                     * xjac * tstep &
                    - v * r0 * (vpar0_s * ps0_t - vpar0_t * ps0_s)                                                            * tstep &

                    + v * 2.d0 * tauIC * p0_y * BigR                                                                   * xjac * tstep &

                    + zeta * v * delta_g(mp,5,ms,mt) * BigR                                                                   * xjac  &

	            + BigR* (-Dn0x*rn0_x*v_x - Dn0y*rn0_y*v_y - Dn0p*rn0_p*v_p*eps_cyl**2/BigR**2)                     * xjac * tstep &  


                    - D_perp_num * (v_xx + v_x/Bigr + v_yy)*(r0_xx + r0_x/Bigr + r0_yy) * BigR * xjac * tstep                         &

                    - TG_num5 * 0.25d0 * BigR**3 * (r0_x * u0_y - r0_y * u0_x)                                                        &
                                                 * ( v_x * u0_y - v_y * u0_x) * xjac * tstep * tstep                                  &
                    - TG_num5 * 0.25d0 / BigR * vpar0**2                                                                              &
                              * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                                                      &
                              * ( v_x * ps0_y -  v_y * ps0_x + F0 / BigR * v_p) * xjac * tstep * tstep        

!###################################################################################################
!#  equation 6 (energy  equation)                                                                  #
!###################################################################################################

         rhs_ij_6 =   v * BigR * heat_source(ms,mt)                                    * xjac * tstep &
 
                    + v * (r0 + rn0*alpha_imp_bis) * BigR**2 * ( T0_s * u0_t - T0_t * u0_s)   * tstep &
                    + v * T0 * BigR**2 * ( r0_s * u0_t - r0_t * u0_s)                         * tstep &
                    + v * alpha_imp * T0 * BigR**2 * (rn0_s * u0_t - rn0_t * u0_s)            * tstep &

                    + v * (r0 + rn0*alpha_imp) * T0 * 2.d0* GAMMA * BigR * u0_y        * xjac * tstep &

                    - v * (r0 + rn0*alpha_imp_bis) * F0 / BigR * Vpar0 * T0_p          * xjac * tstep &
                    - v * T0 * F0 / BigR * Vpar0 * (r0_p + alpha_imp * rn0_p)          * xjac * tstep &

                    - v * (r0 + rn0*alpha_imp_bis) * Vpar0 * (T0_s * ps0_t - T0_t * ps0_s)    * tstep &
                    - v * T0 * Vpar0 * (r0_s * ps0_t - r0_t * ps0_s)                          * tstep &
                    - v * T0 * Vpar0 * alpha_imp * (rn0_s * ps0_t - rn0_t * ps0_s)            * tstep &

                    - v * (r0 + rn0*alpha_imp) * T0 * GAMMA * (vpar0_s * ps0_t - vpar0_t * ps0_s)        * tstep &
                    - v * (r0 + rn0*alpha_imp) * T0 * GAMMA * F0 / BigR * vpar0_p                 * xjac * tstep &

                    - (ZKpar_T-ZK_prof) * BigR / BB2 * Bgrad_T_star * Bgrad_T          * xjac * tstep &
                    - ZK_prof * BigR * (v_x*T0_x + v_y*T0_y + v_p*T0_p /BigR**2 )      * xjac * tstep &

                    - ZK_perp_num  *  (v_xx + v_x/Bigr + v_yy)*(T0_xx + T0_x/Bigr + T0_yy) * BigR * xjac * tstep &
                     
                    - ZK_par_num * (v_ps0_x  * ps0_y - v_ps0_y  * ps0_x) &
                                 * (T0_ps0_x * ps0_y - T0_ps0_y * ps0_x)               * xjac * tstep &

                    - TG_num6 * 0.25d0 * BigR**3 * T0 * (r0_x * u0_y - r0_y * u0_x) &
                                       * ( v_x * u0_y - v_y * u0_x) * xjac * tstep * tstep  &
                    - TG_num6 * 0.25d0 * BigR**3 * r0 * (T0_x * u0_y - T0_y * u0_x) &
                                       * ( v_x * u0_y - v_y * u0_x) * xjac * tstep * tstep  &

                    - TG_num6 * 0.25d0 / BigR * vpar0**2 &
                              * T0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                         &
                              * ( v_x * ps0_y -  v_y * ps0_x + F0 / BigR * v_p) * xjac * tstep * tstep        &

                    - TG_num6 * 0.25d0 / BigR * vpar0**2 &
                              * r0 * (T0_x * ps0_y - T0_y * ps0_x + F0 / BigR * T0_p)                         &
                              * ( v_x * ps0_y -  v_y * ps0_x + F0 / BigR * v_p) * xjac * tstep * tstep        &

                    + zeta * v * (r0 + rn0 * alpha_imp_bis) * delta_g(mp,6,ms,mt) * BigR               * xjac &
                    + zeta * v * T0 * delta_g(mp,5,ms,mt) * BigR                                       * xjac &
                    + zeta * v * alpha_imp * T0 * delta_g(mp,8,ms,mt) * BigR                           * xjac &   

                    + v * BigR * (2/(3 * BigR**2)) * eta_Sp * zj0**2                    * xjac * tstep  &
                    - v * BigR * (r0_corr+beta_imp*rn0_corr) * rn0_corr * Lrad          * xjac * tstep  &
                    - v * BigR * r0_corr * frad_bg                                      * xjac * tstep

!###################################################################################################
!#  equation 7 (parallel velocity  equation)                                                       #
!###################################################################################################

         rhs_ij_7 = - v * F0 / BigR * P0_p                                              * xjac * tstep &
                    - v * (P0_s * ps0_t - P0_t * ps0_s)                                        * tstep &

                      - visco_par * (v_x * vpar0_x + v_y * vpar0_y) * BigR                * xjac * tstep &


                    !- v*(particle_source(ms,mt) + source_pellet) * vpar0 * BB2   * BigR * xjac * tstep &
		       
		      - 0.5d0 * r0 * vpar0**2 * BB2 * (ps0_s * v_t - ps0_t * v_s)                * tstep &
                      + 0.5d0 * r0 * vpar0**2 * BB2 * F0 / BigR * v_p                     * xjac * tstep &

                      - 0.5d0 * v * vpar0**2 * BB2 * (ps0_s * r0_t - ps0_t * r0_s)                * tstep &
                      + 0.5d0 * v * vpar0**2 * BB2 * F0 / BigR * r0_p                      * xjac * tstep &

                    - visco_par_num * (v_xx + v_x/Bigr + v_yy)*(vpar0_xx + vpar0_x/Bigr + vpar0_yy) * BigR * xjac * tstep &
                                  
                      + zeta * v * delta_g(mp,7,ms,mt) * R0 * F0**2 / BigR                       * xjac  &
                      ! Why it is F0**2 rather than BB2 here?
		      + zeta * v * r0 * vpar0 * (ps0_x * delta_ps_x + ps0_y * delta_ps_y) / BigR * xjac  &

                      + zeta * v * delta_g(mp,5,ms,mt) * vpar0 * F0**2 / BigR                    * xjac  &
                      ! Additional zeta term for timesteping here (new momentum)
                      + (r0_x_hat * u0_y - r0_y_hat * u0_x) * vpar0 * BB2 * v * xjac* tstep &
                      ! Perp component of the fifth term of Eq.29 here (new momentum)
                      - r0 * F0 / BigR * vpar0_p * vpar0 * BB2 * v * xjac* tstep &
                      - vpar0 * F0 / BigR * r0_p * vpar0 * BB2 * v * xjac* tstep &
                      ! Toroidal para component of the fifth term of Eq.29 here (new momentum)
                      - r0 * (vpar0_x * ps0_y - vpar0_y * ps0_x) &
                        * vpar0 * BB2 * v * xjac* tstep  &
                      - vpar0 * (r0_x * ps0_y - r0_y * ps0_x) &
                        * vpar0 * BB2 * v * xjac* tstep  &
                      ! Poloidal para component of the fifth term of Eq.29 here (new momentum)


                      - TG_NUM7 * 0.25d0 * r0 * Vpar0**2 * BB2 &
                             * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR  &
                             * (-(ps0_s * v_t     - ps0_t * v_s)    /xjac + F0 / BigR * v_p)  * xjac * tstep * tstep &

!                      - TG_NUM7 * 0.25d0 * v  * Vpar0**2 * BB2 &
!                             * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR  &
!                             * (-(ps0_s * r0_t    - ps0_t * r0_s)   /xjac + F0 / BigR * r0_p)  * xjac * tstep * tstep &

!=============================== New TG_num terms==================================

                      - TG_NUM7 * 0.25d0 * vpar0 * Vpar0**2 * BB2 &
                             * (-(ps0_s * r0_t - ps0_t * r0_s)/xjac + F0 / BigR * r0_p) / BigR  &
                             * (-(ps0_s * v_t     - ps0_t * v_s)    /xjac + F0 / BigR * v_p)  * xjac * tstep * tstep 

!===============================End of new TG_num terms============================


!=============================== New TG_num !terms==================================
!
!            - TG_NUM7 * 0.25d0 * vpar0 * Vpar0**2 * BB2 &
!                      * (-(ps0_s * r0_t - ps0_t * r0_s)/xjac + F0 / BigR * r0_p) / BigR  &
!                      * (-(ps0_s * v_t  - ps0_t * v_s) /xjac + F0 / BigR * v_p )  * xjac * tstep * tstep &
!            - TG_NUM7 * 0.25d0 * v  * Vpar0**2 * BB2 &
!                      * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR  &
!                      * (-(ps0_s * r0_t    - ps0_t * r0_s)   /xjac + F0 / BigR * r0_p)  * xjac * tstep * tstep
!
!===============================End of new TG_num terms============================


!################################################################################################### 
!#  equation 8 (impurity density equation)                                                          # 
!################################################################################################### 


	   rhs_ij_8 = BigR* (- Dn0x * rn0_x * v_x - Dn0y * rn0_y * v_y - Dn0p * rn0_p * v_p*eps_cyl**2/BigR**2)   * xjac * tstep &        
                      ! New diffusion scheme for impurities
                      - (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star * (Bgrad_rhon)                                      * xjac * tstep & 
                      - D_prof * BigR  * (v_x*(rn0_x) + v_y*(rn0_y) + v_p*(rn0_p) * eps_cyl**2 /BigR**2 )                * xjac * tstep & 
                      + v * BigR**2 * ( rn0_s * u0_t - rn0_t * u0_s)                                                     * tstep &
                      + v * 2.d0 * BigR * rn0 * u0_y                                                              * xjac * tstep &
                      - v * F0 / BigR * Vpar0 * rn0_p                                                             * xjac * tstep &
                      - v * Vpar0 * (rn0_s * ps0_t - rn0_t * ps0_s)                                                      * tstep &
                      - v * F0 / BigR * rn0 * vpar0_p                                                             * xjac * tstep &
                      - v * rn0 * (vpar0_s * ps0_t - vpar0_t * ps0_s)                                                    * tstep &

                    - TG_num8 * 0.25d0 * BigR**3 * (rn0_x * u0_y - rn0_y * u0_x)                              &
                                                 * ( v_x * u0_y - v_y * u0_x) * xjac * tstep * tstep          &
                    - TG_num8 * 0.25d0 / BigR * vpar0**2                                                      &
                              * (rn0_x * ps0_y - rn0_y * ps0_x + F0 / BigR * rn0_p)                           &
                              * ( v_x * ps0_y -  v_y * ps0_x + F0 / BigR * v_p) * xjac * tstep * tstep        &

                      + BigR * v * source_mgi                                                                     * xjac * tstep &

                      + v * delta_g(mp,8,ms,mt) * BigR * xjac * zeta                                          &
                      - Dn_perp_num * (v_xx + v_x/Bigr + v_yy)*(rn0_xx + rn0_x/Bigr + rn0_yy) * BigR * xjac * tstep

!###################################################################################################
!#  RHS equations end                                                                                  #
!###################################################################################################

         ij1 = index_ij
         ij2 = index_ij + 1*n_tor
         ij3 = index_ij + 2*n_tor
         ij4 = index_ij + 3*n_tor
         ij5 = index_ij + 4*n_tor
         ij6 = index_ij + 5*n_tor
         ij7 = index_ij + 6*n_tor
         ij8 = index_ij + 7*n_tor

           RHS(ij1) = RHS(ij1) + rhs_ij_1 * wst
           RHS(ij2) = RHS(ij2) + rhs_ij_2 * wst
           RHS(ij3) = RHS(ij3) + rhs_ij_3 * wst
           RHS(ij4) = RHS(ij4) + rhs_ij_4 * wst
           RHS(ij5) = RHS(ij5) + rhs_ij_5 * wst
           RHS(ij6) = RHS(ij6) + rhs_ij_6 * wst
           RHS(ij7) = RHS(ij7) + rhs_ij_7 * wst
           RHS(ij8) = RHS(ij8) + rhs_ij_8 * wst
	   
           do k=1,n_vertex_max

             do l=1,n_order+1

               do in = 1, n_tor

                 psi   = H(k,l,ms,mt) * element%size(k,l) * HZ(in,mp)

                 psi_x = (   y_t(ms,mt) * h_s(k,l,ms,mt) - y_s(ms,mt) * h_t(k,l,ms,mt) ) / xjac    &
                              * element%size(k,l) * HZ(in,mp)

                 psi_y = ( - x_t(ms,mt) * h_s(k,l,ms,mt) + x_s(ms,mt) * h_t(k,l,ms,mt) )  / xjac   &
                              * element%size(k,l) * HZ(in,mp)

                 psi_p  = H(k,l,ms,mt)    * element%size(k,l) * HZ_p(in,mp)
                 psi_s  = h_s(k,l,ms,mt)  * element%size(k,l) * HZ(in,mp)
                 psi_t  = h_t(k,l,ms,mt)  * element%size(k,l) * HZ(in,mp)
                 psi_ss = h_ss(k,l,ms,mt) * element%size(k,l) * HZ(in,mp)
                 psi_tt = h_tt(k,l,ms,mt) * element%size(k,l) * HZ(in,mp)
                 psi_st = h_st(k,l,ms,mt) * element%size(k,l) * HZ(in,mp)

                 psi_xx = (psi_ss * y_t(ms,mt)**2 - 2.d0*psi_st * y_s(ms,mt)*y_t(ms,mt) + psi_tt * y_s(ms,mt)**2  &	        
		        + psi_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                              &	   
	                + psi_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )    / xjac**2               &		
		        - xjac_x * (psi_s * y_t(ms,mt) - psi_t * y_s(ms,mt)) / xjac**2

	         psi_yy = (psi_ss * x_t(ms,mt)**2 - 2.d0*psi_st * x_s(ms,mt)*x_t(ms,mt) + psi_tt * x_s(ms,mt)**2  &	        
		        + psi_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                              &	   
	                + psi_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )    / xjac**2               &		
		        - xjac_y * (- psi_s * x_t(ms,mt) + psi_t * x_s(ms,mt) ) / xjac**2
		
	         psi_xy = (- psi_ss * y_t(ms,mt)*x_t(ms,mt) - psi_tt * x_s(ms,mt)*y_s(ms,mt)                      &
     	                + psi_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                             &        
                        - psi_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                             &	   
	                - psi_t  * (x_st(ms,mt)*y_s(ms,mt) - x_ss(ms,mt)*y_t(ms,mt) )  )  / xjac**2               &		
                        - xjac_x * (- psi_s * x_t(ms,mt) + psi_t * x_s(ms,mt) )   / xjac**2
 
                 u    = psi    ;  zj    = psi    ;  w    = psi    ; rho    = psi    ;  T    = psi    ; vpar    = psi    ; rhon   = psi  
                 u_x  = psi_x  ;  zj_x  = psi_x  ;  w_x  = psi_x  ; rho_x  = psi_x  ;  T_x  = psi_x  ; vpar_x  = psi_x  ; rhon_x = psi_x
                 u_y  = psi_y  ;  zj_y  = psi_y  ;  w_y  = psi_y  ; rho_y  = psi_y  ;  T_y  = psi_y  ; vpar_y  = psi_y  ; rhon_y = psi_y
                 u_p  = psi_p  ;  zj_p  = psi_p  ;  w_p  = psi_p  ; rho_p  = psi_p  ;  T_p  = psi_p  ; vpar_p  = psi_p  ; rhon_p = psi_p
                 u_s  = psi_s  ;  zj_s  = psi_s  ;  w_s  = psi_s  ; rho_s  = psi_s  ;  T_s  = psi_s  ; vpar_s  = psi_s  ; rhon_s = psi_s
                 u_t  = psi_t  ;  zj_t  = psi_t  ;  w_t  = psi_t  ; rho_t  = psi_t  ;  T_t  = psi_t  ; vpar_t  = psi_t  ; rhon_t = psi_t
                 u_ss = psi_ss ;  zj_ss = psi_ss ;  w_ss = psi_ss ; rho_ss = psi_ss ;  T_ss = psi_ss ; vpar_ss = psi_ss ; rhon_ss = psi_ss
                 u_tt = psi_tt ;  zj_tt = psi_tt ;  w_tt = psi_tt ; rho_tt = psi_tt ;  T_tt = psi_tt ; vpar_tt = psi_tt ; rhon_tt = psi_tt
                 u_st = psi_st ;  zj_st = psi_st ;  w_st = psi_st ; rho_st = psi_st ;  T_st = psi_st ; vpar_st = psi_st ; rhon_st = psi_st

                 u_xx = psi_xx ;                                    rho_xx = psi_xx ;  T_xx = psi_xx ; vpar_xx = psi_xx ; rhon_xx = psi_xx
                 u_yy = psi_yy ;                                    rho_yy = psi_yy ;  T_yy = psi_yy ; vpar_yy = psi_yy ; rhon_yy = psi_yy
                 u_xy = psi_xy ;                                    rho_xy = psi_xy ;  T_xy = psi_xy ; vpar_xy = psi_xy

                 w_xx = (psi_ss * y_t(ms,mt)**2 - 2.d0*psi_st * y_s(ms,mt)*y_t(ms,mt) + psi_tt * y_s(ms,mt)**2 ) / xjac**2       
                 w_yy = (psi_ss * x_t(ms,mt)**2 - 2.d0*psi_st * x_s(ms,mt)*x_t(ms,mt) + psi_tt * x_s(ms,mt)**2 ) / xjac**2
                 w_xy = psi_xy          

                 rho_hat   = BigR**2 * rho
                 rho_x_hat = 2.d0 * BigR * BigR_x  * rho + BigR**2 * rho_x
                 rho_y_hat = BigR**2 * rho_y
                  Btheta2_psi  = 2.d0 * (psi_x * ps0_x + psi_y * ps0_y ) /BigR**2
                 rhon_hat   = BigR**2 * rhon                                     
                 rhon_x_hat = 2.d0 * BigR * BigR_x  * rhon + BigR**2 * rhon_x    
                 rhon_y_hat = BigR**2 * rhon_y                                   

                 index_kl = n_tor*n_var*(n_order+1)*(k-1) + n_tor * n_var * (l-1) + in   ! index in the ELM matrix

!###################################################################################################
!#  equation 1   (induction equation)                                                              #
!###################################################################################################

                 amat_11 = v * psi / BigR * xjac * (1.d0 + zeta)                                     &
                         - v * (psi_s * u0_t - psi_t * u0_s)                        * theta * tstep  &

                          + v * tauIC/(r0_corr*BB2)*F0**2/BigR**2 * (psi_s * p0_t - psi_t * p0_s) * theta * tstep 
 
                 ! term with BB2 still missing

                 amat_12 = - v * (ps0_s * u_t - ps0_t * u_s)                        * theta * tstep  &
                           + eps_cyl * F0 / BigR * v * u_p * xjac                   * theta * tstep

                 amat_13 = - eta_num_T * (v_x * zj_x + v_y * zj_y)           * xjac * theta * tstep  &
                           - eta_T * v * zj / BigR                           * xjac * theta * tstep
   
                  amat_15 =  v * tauIC/(r0_corr*BB2)*F0**2/BigR**2* T0  * (ps0_s * rho_t - ps0_t * rho_s) * theta * tstep &
                           + v * tauIC/(r0_corr*BB2)*F0**2/BigR**2* rho * (ps0_s * T0_t  - ps0_t * T0_s)  * theta * tstep & 
                           - v * tauIC/(r0_corr*BB2)*F0**3/BigR**3* eps_cyl * T0  * rho_p * xjac * theta * tstep &
                           - v * tauIC/(r0_corr*BB2)*F0**3/BigR**3* eps_cyl * rho * T0_p  * xjac * theta * tstep &
   
                           - v * tauIC * rho /(r0_corr**2 * BB2) * F0**2/BigR**2 * (ps0_s * p0_t - ps0_t * p0_s) * theta * tstep &      
                           + v * tauIC * rho /(r0_corr**2 * BB2) * F0**3/BigR**3 * eps_cyl * p0_p * xjac         * theta * tstep 


                 amat_16 = - deta_dT * v * T * (zj0-current_source(ms,mt)) / BigR * xjac * theta * tstep &
 
                        + v * tauIC/(r0_corr*BB2) * F0**2/BigR**2 * r0 * (ps0_s * T_t  - ps0_t * T_s) * theta * tstep &
                        + v * tauIC/(r0_corr*BB2) * F0**2/BigR**2 * T  * (ps0_s * r0_t - ps0_t * r0_s)* theta * tstep &   
                        - v * tauIC/(r0_corr*BB2) * F0**3/BigR**3 * eps_cyl * r0 * T_p  * xjac * theta * tstep &
                        - v * tauIC/(r0_corr*BB2) * F0**3/BigR**3 * eps_cyl * T  * r0_p * xjac * theta * tstep 


!---------------------------------------------------------------- equation 1
!                 amat_11 = v * psi / BigR * xjac * (1.d0 + zeta)                                        &
!                         + eta_T * (psi_x * v_x + psi_y * v_y) / BigR * xjac           * theta * tstep  &
!                         - v * (psi_s * u0_t - psi_t * u0_s)                           * theta * tstep  &
!                         + v * deta_dT * (T0_x * psi_x  + T0_y * psi_y ) / BigR * xjac * theta * tstep

!###################################################################################################
!#  equation 2   (perpendicular momentum equation)                                                 #
!###################################################################################################

                 amat_21 = - v * (psi_s * zj0_t - psi_t * zj0_s)                          * theta * tstep &

                      - BigR**2 * r0 * (vpar0_x * psi_y - vpar0_y * psi_x) &
                        * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep &
                      - BigR**2 * vpar0 * (r0_x * psi_y - r0_y * psi_x) &
                        * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep
                      ! Poloidal para component of the third term of Eq.20 here (new momentum)


                 amat_22 = - BigR * r0_hat * (v_x * u_x + v_y * u_y) * xjac * (1.d0 + zeta)                                 &
                           + r0_hat * BigR**2 * w0 * (v_s * u_t  - v_t  * u_s)                              * theta * tstep &
                           + BigR**2 * (u_x * u0_x + u_y * u0_y) * (v_x * r0_y_hat - v_y * r0_x_hat) * xjac * theta * tstep &
			   
		           + tauIC * BigR**3 * p0_y * (v_x* u_x + v_y * u_y)                         * xjac * theta * tstep &
		           + v * tauIC * BigR**4 * (u_xy * (p0_xx - p0_yy) - p0_xy * (u_xx - u_yy))  * xjac * theta * tstep &
                           !- BigR**3 * (particle_source(ms,mt)+source_pellet) * (v_x * u_x + v_y * u_y) * xjac * theta * tstep &

                       + BigR**2 * (r0_x_hat * u0_y - r0_y_hat * u0_x) * (v_x * u_x + v_y * u_y) &
                         * xjac * theta * tstep &
                       + BigR**2 * (r0_x_hat * u_y - r0_y_hat * u_x) * (v_x * u0_x + v_y * u0_y) &
                         * xjac * theta * tstep &
                      ! Perp component of the third term of Eq.20 here (new momentum)
                      - BigR**2 * r0 * F0 / BigR * vpar0_p * (v_x * u_x + v_y * u_y) * xjac * theta * tstep &
                      - BigR**2 * vpar0 * F0 / BigR * r0_p * (v_x * u_x + v_y * u_y) * xjac * theta * tstep &
                      ! Toroidal para component of the third term of Eq.20 here (new momentum)
                      - BigR**2 * r0 * (vpar0_x * ps0_y - vpar0_y * ps0_x) &
                        * (v_x * u_x + v_y * u_y) * xjac * theta * tstep  &
                      - BigR**2 * vpar0 * (r0_x * ps0_y - r0_y * ps0_x) &
                        * (v_x * u_x + v_y * u_y) * xjac * theta * tstep  &
                      ! Poloidal para component of the third term of Eq.20 here (new momentum)


                           + TG_num2 * 0.25d0 * r0_hat * BigR**3 * (w0_x * u_y - w0_y * u_x)       &
                                     * ( v_x * u0_y - v_y * u0_x) * xjac * theta * tstep * tstep   &

                           + TG_num2 * 0.25d0 * r0_hat * BigR**3 * (w0_x * u0_y - w0_y * u0_x)     &
                                     * ( v_x * u_y - v_y * u_x)   * xjac * theta * tstep * tstep   & 

!====================================New TG_num terms=================================
                      + TG_num2 * 0.25d0 * w0 * BigR**3 * BigR**2 * (r0_x * u_y - r0_y * u_x) &
                                * ( v_x * u0_y - v_y * u0_x) * theta * xjac * tstep * tstep &

                      + TG_num2 * 0.25d0 * w0 * BigR**3 * BigR**2 * (r0_x * u0_y - r0_y * u0_x) &
                                * ( v_x * u_y - v_y * u_x) * theta * xjac * tstep * tstep

!===============================End of NewTG_num terms==============================

   
                 amat_23 = - v * (ps0_s * zj_t  - ps0_t * zj_s)                           * theta * tstep  &
                           + eps_cyl * F0 / BigR * v * zj_p  * xjac                       * theta * tstep

                 amat_24 = r0_hat * BigR**2 * w  * ( v_s * u0_t - v_t * u0_s)  * theta * tstep  &
                         + BigR * ( v_x * w_x + v_y * w_y) * visco_T  * xjac   * theta * tstep  &

                         + visco_num_T * (v_xx + v_x/BigR + v_yy)*(w_xx + w_x/BigR + w_yy) * xjac * theta * tstep    &
     
                        + v * tauIC * BigR**4 * (p0_s * w_t - p0_t * w_s)              * theta * tstep &

                         + TG_num2 * 0.25d0 * r0_hat * BigR**3 * (w_x * u0_y - w_y * u0_x)     &
                                   * ( v_x * u0_y - v_y * u0_x) * xjac * theta * tstep * tstep &

!====================================New TG_num terms=================================
                      + TG_num2 * 0.25d0 * w * BigR**3 * BigR**2 * (r0_x * u0_y - r0_y * u0_x) &
                                * ( v_x * u0_y - v_y * u0_x) * theta * xjac * tstep * tstep

!===============================End of NewTG_num terms==============================


                 amat_25 = - BigR**3 * rho * (v_x * u0_x + v_y * u0_y) * xjac  * (1.d0 + zeta)          &
                      ! Additional zeta term for timesteping here (new momentum)

                      + BigR**2 * (rho_x_hat * u0_y - rho_y_hat * u0_x) &
                        * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep &
                      ! Perp component of the third term of Eq.20 here (new momentum)
                      - BigR**2 * rho * F0 / BigR * vpar0_p * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep &
                      - BigR**2 * vpar0 * F0 / BigR * rho_p * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep &
                      ! Toroidal para component of the third term of Eq.20 here (new momentum)
                      - BigR**2 * rho * (vpar0_x * ps0_y - vpar0_y * ps0_x) &
                        * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep  &
                      - BigR**2 * vpar0 * (rho_x * ps0_y - rho_y * ps0_x) &
                        * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep  &
                      ! Poloidal para component of the third term of Eq.20 here (new momentum)

                           + 0.5d0 * vv2 * (v_x * rho_y_hat - v_y * rho_x_hat)   * xjac * theta * tstep &
                           + rho_hat * BigR**2 * w0 * (v_s * u0_t - v_t * u0_s)         * theta * tstep &
                           - BigR**2 * (v_s * rho_t * T0   - v_t * rho_s * T0  )        * theta * tstep &
                           - BigR**2 * (v_s * rho   * T0_t - v_t * rho   * T0_s)        * theta * tstep &

                           + v * tauIC * BigR**4 * T0  * (rho_s * w0_t - rho_t * w0_s)  * theta * tstep &
                           + v * tauIC * BigR**4 * rho * (T0_s  * w0_t - T0_t  * w0_s)  * theta * tstep &

                           + tauIC * BigR**3 * (T0_y * rho + T0 * rho_y) * (v_x* u0_x + v_y * u0_y) * xjac * theta * tstep &

                           + v * tauIC * BigR**4 * ( (u0_xy * (rho_xx*T0 + 2.d0*rho_x*T0_x + rho*T0_xx           &
                                                            -  rho_yy*T0 - 2.d0*rho_y*T0_y - rho*T0_yy))         &
                                                   - (rho_xy * T0 + rho_x*T0_y + rho_y*T0_x + rho*T0_xy) * (u0_xx - u0_yy)  )   &
                                       * xjac * theta * tstep                                          &

                         + TG_num2 * 0.25d0 * rho_hat * BigR**3 * (w0_x * u0_y - w0_y * u0_x) &
                                  * ( v_x * u0_y - v_y * u0_x) * xjac * theta * tstep * tstep &

!====================================New TG_num terms=================================
                      + TG_num2 * 0.25d0 * w * BigR**3 * BigR**2 * (r0_x * u0_y - r0_y * u0_x) &
                                * ( v_x * u0_y - v_y * u0_x) * theta * xjac * tstep * tstep

!===============================End of NewTG_num terms==============================


                 amat_26 = - BigR**2 * (v_s * r0_t * T   - v_t * r0_s * T)      * theta * tstep  &
                           - BigR**2 * (v_s * r0   * T_t - v_t * r0   * T_s)    * theta * tstep  &
                           + dvisco_dT * T * ( v_x * w0_x + v_y * w0_y ) * BigR * xjac  * theta * tstep &

                           + v * tauIC * BigR**4 * r0 * (T_s  * w0_t - T_t  * w0_s)  * theta * tstep &
                           + v * tauIC * BigR**4 * T  * (r0_s * w0_t - r0_t * w0_s)  * theta * tstep &

		            + tauIC * BigR**3 * (r0_y * T + r0 * T_y) * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep &
		            + v * tauIC * BigR**4 * ( (u0_xy * (T_xx * r0 + 2.d0 * T_x * r0_x + T * r0_xx         &
			                                     - T_yy*r0 - 2.d0 * T_y * r0_y - T * r0_yy))       &			 
			                           - (T_xy * r0 + T_x*r0_y + T_y*r0_x + T*r0_xy) * (u0_xx - u0_yy)  )         &
						 * xjac * theta * tstep 

             amat_27 = - BigR**2 * r0 * F0 / BigR * vpar_p * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep &
                      - BigR**2 * vpar * F0 / BigR * r0_p * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep &
                      ! Toroidal para component of the third term of Eq.20 here (new momentum)
                      - BigR**2 * r0 * (vpar_x * ps0_y - vpar_y * ps0_x) &
                        * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep  &
                      - BigR**2 * vpar * (r0_x * ps0_y - r0_y * ps0_x) &
                        * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep
                      ! Poloidal para component of the third term of Eq.20 here (new momentum)

                amat_28 = 0     ! Place holder

!###################################################################################################
!#  equation 3                                                                                     #
!###################################################################################################

                 amat_33 = v * zj / BigR * xjac                              
                 amat_31 = (v_x * psi_x + v_y * psi_y ) / BigR * xjac         

!###################################################################################################
!#  equation 4                                                                                     #
!###################################################################################################

                 amat_44 =  v * w * BigR * xjac                               
                 amat_42 = (v_x * u_x + v_y * u_y) * BigR * xjac              

!###################################################################################################
!#  equation 5    continuity equation (total density)                                              #
!###################################################################################################

                 Bgrad_rho_star_psi = ( v_x  * psi_y - v_y  * psi_x ) / BigR
                 Bgrad_rho_psi      = ( r0_x * psi_y - r0_y * psi_x ) / BigR
                 Bgrad_rhon_psi     = ( rn0_x * psi_y - rn0_y * psi_x ) / BigR
                 Bgrad_rho_rho      = ( F0 / BigR * rho_p +  rho_x * ps0_y - rho_y * ps0_x ) / BigR
                 Bgrad_rho_rhon     = ( F0 / BigR * rhon_p +  rhon_x * ps0_y - rhon_y * ps0_x ) / BigR
                 BB2_psi            = 2.d0 * (psi_x * ps0_x + psi_y * ps0_y ) /BigR**2

                 ! New impurity diffusion scheme
                 amat_51 = - (D_par-D_prof) * BigR * BB2_psi / BB2**2 * Bgrad_rho_star * (Bgrad_rho)                     * xjac * theta * tstep &
                           !- (D_par-D_prof) * BigR * BB2_psi / BB2**2 * Bgrad_rho_star * (Bgrad_rho-Bgrad_rhon)         * xjac * theta * tstep &
                           !+ (D_par-D_prof) * BigR / BB2              * Bgrad_rho_star_psi * (Bgrad_rho-Bgrad_rhon)     * xjac * theta * tstep &
                           + (D_par-D_prof) * BigR / BB2              * Bgrad_rho_star_psi * (Bgrad_rho)                 * xjac * theta * tstep &
                           !+ (D_par-D_prof) * BigR / BB2              * Bgrad_rho_star * (Bgrad_rho_psi-Bgrad_rhon_psi) * xjac * theta * tstep &
                           + (D_par-D_prof) * BigR / BB2              * Bgrad_rho_star * (Bgrad_rho_psi)                 * xjac * theta * tstep &
                           + v * Vpar0 * (r0_s * psi_t - r0_t * psi_s)                                            * theta * tstep &
                           + v * r0 * (vpar0_s * psi_t - vpar0_t * psi_s)                                         * theta * tstep &

                           + TG_num5 * 0.25d0 / BigR * vpar0**2                                                              &
                                     * (r0_x * psi_y - r0_y * psi_x)                                                         &
                                     * ( v_x * ps0_y -  v_y * ps0_x + F0 / BigR * v_p) * xjac * theta * tstep * tstep        &

                           + TG_num5 * 0.25d0 / BigR * vpar0**2                                                              &
                                     * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                                      &
                                     * ( v_x * psi_y -  v_y * psi_x                   ) * xjac * theta * tstep * tstep


                 amat_52 = - v * BigR**2 * ( r0_s * u_t - r0_t * u_s)                                     * theta * tstep &
                           - v * 2.d0 * BigR * r0 * u_y                                            * xjac * theta * tstep &

                           + TG_num5 * 0.25d0 * BigR**3 * (r0_x * u_y  - r0_y * u_x)                                      &
                                                        * ( v_x * u0_y - v_y  * u0_x) * xjac * theta * tstep * tstep      &

                           + TG_num5 * 0.25d0 * BigR**3 * (r0_x * u0_y - r0_y * u0_x)                                     &
                                                        * ( v_x * u_y  - v_y  * u_x)  * xjac * theta * tstep * tstep 
			   
                 amat_55 = v * rho * BigR * xjac * (1.d0 + zeta)  &
                         - v * BigR**2 * ( rho_s * u0_t - rho_t * u0_s)                                       * theta * tstep &
                         - v * 2.d0 * BigR * rho * u0_y                                                * xjac * theta * tstep &
                         + (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star * Bgrad_rho_rho                * xjac * theta * tstep &
                         + D_prof * BigR  * (v_x*rho_x + v_y*rho_y + v_p*rho_p * eps_cyl**2 /BigR**2 ) * xjac * theta * tstep &
                         + v * F0 / BigR * Vpar0 * rho_p                                               * xjac * theta * tstep &
                         + v * Vpar0 * (rho_s * ps0_t - rho_t * ps0_s)                                        * theta * tstep &
                         + v * rho * (vpar0_s * ps0_t - vpar0_t * ps0_s)                                      * theta * tstep &
                         + v * rho * F0 / BigR * vpar0_p                                               * xjac * theta * tstep &

                        - v * 2.d0 * tauIC * (rho_y * T0 + rho*T0_y) * BigR                           * xjac * theta * tstep &
                                                 
                         + D_perp_num     * (v_xx + v_x/BigR + v_yy)*(rho_xx + rho_x/BigR + rho_yy)   * BigR * xjac * theta * tstep &

                         + TG_num5 * 0.25d0 * BigR**3 * (rho_x * u0_y - rho_y * u0_x)                                &
                                                      * ( v_x  * u0_y - v_y   * u0_x) * xjac * theta * tstep * tstep &

                         + TG_num5 * 0.25d0 / BigR * vpar0**2                                                &
                              * (rho_x * ps0_y - rho_y * ps0_x + F0 / BigR * rho_p)                          &
                              * ( v_x * ps0_y -  v_y * ps0_x   + F0 / BigR * v_p) * xjac * theta * tstep * tstep 

		 amat_56 = - v * 2.d0 * tauIC * (T_y * r0 + T*r0_y) * BigR                               * xjac * theta * tstep 

                 amat_57 = + v * F0 / BigR * Vpar * r0_p                                             * xjac * theta * tstep &
                           + v * Vpar * (r0_s * ps0_t - r0_t * ps0_s)                                       * theta * tstep &
                           + v * r0 * (vpar_s * ps0_t - vpar_t * ps0_s)                                     * theta * tstep &
                           + v * r0 * F0 / BigR * vpar_p                                             * xjac * theta * tstep &

                           + TG_num5 * 0.25d0 / BigR * 2.d0*vpar0*vpar                                 &
                              * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                       &
                              * ( v_x * ps0_y -  v_y * ps0_x + F0 / BigR * v_p) * xjac * theta * tstep * tstep 

                 amat_58 = &
                           ! New diffusion scheme for impurities
                           !- (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star * Bgrad_rho_rhon                  * xjac * theta * tstep &
                           !- D_prof * BigR  * (v_x*rhon_x + v_y*rhon_y + v_p*rhon_p * eps_cyl**2 /BigR**2 ) * xjac * theta * tstep &

                           + BigR * (Dn0x * rhon_x * v_x + Dn0y * rhon_y * v_y + Dn0p * rhon_p * v_p*eps_cyl**2/BigR**2)    &
                                                                                                     * xjac * theta * tstep

!###################################################################################################
!#  equation 6   energy equation                                                                   #
!###################################################################################################

                 Bgrad_T_star_psi = ( v_x  * psi_y - v_y  * psi_x  ) / BigR
                 Bgrad_T_psi      = ( T0_x * psi_y - T0_y * psi_x )  / BigR
                 Bgrad_T_T        = ( F0 / BigR * T_p +  T_x * ps0_y - T_y * ps0_x ) / BigR          ! F0 due to absence of normalisation

                 T_ps0_x = T_xx * ps0_y - T_xy * ps0_x + T_x * ps0_xy - T_y * ps0_xx
                 T_ps0_y = T_xy * ps0_y - T_yy * ps0_x + T_x * ps0_yy - T_y * ps0_xy
 
                 T0_psi_x = T0_xx * psi_y - T0_xy * psi_x + T0_x * psi_xy - T0_y * psi_xx
                 T0_psi_y = T0_xy * psi_y - T0_yy * psi_x + T0_x * psi_yy - T0_y * psi_xy
                 
                 v_psi_x = v_xx * psi_y - v_xy * psi_x + v_x * psi_xy - v_y * psi_xx
                 v_psi_y = v_xy * psi_y - v_yy * psi_x + v_x * psi_yy - v_y * psi_xy


                 amat_61 = - (ZKpar_T-ZK_prof) * BigR * BB2_psi / BB2**2 * Bgrad_T_star * Bgrad_T     * xjac * theta * tstep &
                           + (ZKpar_T-ZK_prof) * BigR / BB2     * Bgrad_T_star_psi      * Bgrad_T     * xjac * theta * tstep &
                           + (ZKpar_T-ZK_prof) * BigR / BB2     * Bgrad_T_star          * Bgrad_T_psi * xjac * theta * tstep &
                           + v * (r0 + rn0 * alpha_imp) * Vpar0 * (T0_s * psi_t - T0_t * psi_s)              * theta * tstep &
                           + v * T0 * Vpar0 * ((r0_s+rn0_s*alpha_imp)*psi_t - (r0_t+rn0_t*alpha_imp)*psi_s)  * theta * tstep &
                           + v * (r0 + rn0 * alpha_imp) * GAMMA * T0 * (vpar0_s * psi_t - vpar0_t * psi_s)   * theta * tstep &

                           + ZK_par_num * (v_psi_x  * ps0_y - v_psi_y  * ps0_x + v_ps0_x * psi_y - v_ps0_y * psi_x)          &
                                        * (T0_ps0_x * ps0_y - T0_ps0_y * ps0_x)                       * xjac * theta * tstep &
                           + ZK_par_num * (T0_psi_x * ps0_y - T0_psi_y * ps0_x + T0_ps0_x * psi_y - T0_ps0_y * psi_x)        &
                                        * (v_ps0_x  * ps0_y - v_ps0_y  * ps0_x)                       * xjac * theta * tstep &

                           + TG_num6 * 0.25d0 / BigR * vpar0**2                                                       &
                                     * T0 * (r0_x * psi_y - r0_y * psi_x)                                             &
                                     * ( v_x * ps0_y -  v_y * ps0_x + F0 / BigR * v_p) * xjac * theta * tstep * tstep &
                     
                           + TG_num6 * 0.25d0 / BigR * vpar0**2                                                       &
                                     * r0 * (T0_x * psi_y - T0_y * psi_x)                                             &
                                     * ( v_x * ps0_y -  v_y * ps0_x + F0 / BigR * v_p) * xjac * theta * tstep * tstep &
                     
                           + TG_num6 * 0.25d0 / BigR * vpar0**2                                                       &
                                     * T0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                          &
                                     * ( v_x * psi_y -  v_y * psi_x ) * xjac * theta * tstep * tstep                  &
                     
                           + TG_num6 * 0.25d0 / BigR * vpar0**2                                                       &
                                     * r0 * (T0_x * ps0_y - T0_y * ps0_x + F0 / BigR * T0_p)                          &
                                     * ( v_x * psi_y -  v_y * psi_x ) * xjac * theta * tstep * tstep


                 amat_62 = - v * (r0 + rn0 * alpha_imp) * BigR**2 * ( T0_s * u_t - T0_t * u_s)             * theta * tstep &
		           - v * T0 * BigR**2 * ((r0_s+rn0_s*alpha_imp)*u_t - (r0_t+rn0_t*alpha_imp)*u_s)  * theta * tstep &
                           - v * (r0 + rn0 * alpha_imp) * 2.d0* GAMMA * BigR * T0 * u_y             * xjac * theta * tstep &

                           + TG_num6 * 0.25d0 * BigR**2 * T0* (r0_x * u_y - r0_y * u_x)       &
                                     * ( v_x * u0_y - v_y * u0_x) * xjac * theta*tstep*tstep  &

                           + TG_num6 * 0.25d0 * BigR**2 * r0* (T0_x * u_y - T0_y * u_x)       &
                                     * ( v_x * u0_y - v_y * u0_x) * xjac * theta*tstep*tstep  &

                           + TG_num6 * 0.25d0 * BigR**2 * T0* (r0_x * u0_y - r0_y * u0_x)     &
                                     * ( v_x * u_y - v_y * u_x) * xjac * theta*tstep*tstep    &

                           + TG_num6 * 0.25d0 * BigR**2 * r0* (T0_x * u0_y - T0_y * u0_x)     &
                                     * ( v_x * u_y - v_y * u_x) * xjac * theta*tstep*tstep 

                amat_63 = - v * BigR * zj * (4/(3 * BigR**2)) * eta_Sp * zj0                      * xjac * theta * tstep


                amat_65 =   v * rho * T0   * BigR * xjac * (1.d0 + zeta)    &

		           - v * rho * BigR**2 * ( T0_s  * u0_t - T0_t  * u0_s)                        * theta * tstep &
		           - v * T0  * BigR**2 * ( rho_s * u0_t - rho_t * u0_s)                        * theta * tstep &

                           - v * rho * 2.d0* GAMMA * BigR * T0 * u0_y                           * xjac * theta * tstep &

                           + v * rho * F0 / BigR * Vpar0 * T0_p                                 * xjac * theta * tstep &
                           + v * T0  * F0 / BigR * Vpar0 * rho_p                                * xjac * theta * tstep &

                           + v * rho * Vpar0 * (T0_s  * ps0_t - T0_t  * ps0_s)                         * theta * tstep &
                           + v * T0  * Vpar0 * (rho_s * ps0_t - rho_t * ps0_s)                         * theta * tstep &

                           + v * rho * GAMMA * T0 * (vpar0_s * ps0_t - vpar0_t * ps0_s)                * theta * tstep &
                           + v * rho * GAMMA * T0 * F0 / BigR * vpar0_p                         * xjac * theta * tstep &

                           + TG_num6 * 0.25d0 * BigR**2 * T0* (rho_x * u0_y - rho_y * u0_x)      &
                                     * ( v_x * u0_y - v_y * u0_x) * xjac * theta*tstep*tstep     &

                           + TG_num6 * 0.25d0 * BigR**2 * rho * (T0_x * u0_y - T0_y * u0_x)      &
                                     * ( v_x * u0_y - v_y * u0_x) * xjac* theta*tstep*tstep      &

                           + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                              * T0 * (rho_x * ps0_y - rho_y * ps0_x + F0 / BigR * rho_p)                      &
                              * ( v_x * ps0_y -  v_y * ps0_x + F0 / BigR * v_p) * xjac * theta * tstep * tstep&

                           + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                              * rho * (T0_x * ps0_y - T0_y * ps0_x + F0 / BigR * T0_p)                        &
                              * ( v_x * ps0_y -  v_y * ps0_x + F0 / BigR * v_p) * xjac * theta * tstep * tstep &

                           + v * BigR * rho * rn0_corr * Lrad                                   * xjac * theta * tstep &
                           + v * BigR * rho * frad_bg                                           * xjac * theta * tstep


                 amat_66 =   v * (r0 + rn0 * alpha_imp_bis) * T * BigR * xjac * (1.d0 + zeta)                     &
!=============== The ionization potential energy term=========================
                           + v * rn0 * dE_ion_dT            * T * BigR * xjac * (1.d0 + zeta)                     &
!================= End ionization potential energy ===========================
                           - v * (r0 + rn0 * alpha_imp_bis) * BigR**2 * ( T_s  * u0_t - T_t  * u0_s) * theta * tstep &
                           - v * T  * BigR**2 * ( r0_s * u0_t - r0_t * u0_s)                       * theta * tstep &
                           - v * alpha_imp_bis * T * BigR**2 * (rn0_s * u0_t - rn0_t * u0_s)       * theta * tstep &

                           - v * (r0 + rn0 * alpha_imp_bis) * 2.d0* GAMMA * BigR * T * u0_y * xjac * theta * tstep &


                           + v * (r0 + rn0 * alpha_imp_bis) * F0 / BigR * Vpar0 * T_p       * xjac * theta * tstep &
                           + v * T * F0  / BigR * Vpar0 * (r0_p + rn0_p * alpha_imp_bis)    * xjac * theta * tstep &

                           + v * (r0 + rn0 * alpha_imp_bis) * Vpar0 * (T_s  * ps0_t - T_t  * ps0_s)* theta * tstep &
                           + v * T  * Vpar0 * (r0_s * ps0_t - r0_t * ps0_s)                        * theta * tstep &
                           + v * alpha_imp_bis * T * Vpar0 * (rn0_s * ps0_t - rn0_t * ps0_s)       * theta * tstep &

                           + v * (r0 + rn0 * alpha_imp_bis) * GAMMA * T * (vpar0_s * ps0_t - vpar0_t * ps0_s) * theta * tstep &
                           + v * (r0 + rn0 * alpha_imp_bis) * GAMMA * T * F0 / BigR * vpar0_p          * xjac * theta * tstep &

                           + (ZKpar_T-ZK_prof) * BigR / BB2 * Bgrad_T_star * Bgrad_T_T      * xjac * theta * tstep &
                           + ZK_prof * BigR * (v_x*T_x + v_y*T_y + v_p*T_p /BigR**2 )       * xjac * theta * tstep &

                           + dZKpar_dT * T * BigR / BB2 * Bgrad_T_star * Bgrad_T            * xjac * theta * tstep &

                           + ZK_perp_num * (v_xx + v_x/BigR + v_yy)*(T_xx + T_x/BigR + T_yy) * BigR * xjac * theta * tstep &

                           + TG_num6 * 0.25d0 * BigR**2 * T* (r0_x * u0_y - r0_y * u0_x)         &
                                     * ( v_x * u0_y - v_y * u0_x) * xjac * theta * tstep * tstep &

                           + TG_num6 * 0.25d0 * BigR**2 * r0* (T_x * u0_y - T_y * u0_x)          &
                                     * ( v_x * u0_y - v_y * u0_x) * xjac * theta * tstep * tstep &

                           + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                              * T * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                          &
                              * ( v_x * ps0_y -  v_y * ps0_x + F0 / BigR * v_p) * xjac * theta * tstep * tstep&

                           + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                              * r0 * (T_x * ps0_y - T_y * ps0_x + F0 / BigR * T_p)                            &
                              * ( v_x * ps0_y -  v_y * ps0_x + F0 / BigR * v_p) * xjac * theta * tstep * tstep &

                           - v * BigR * T * ((2d0)/(3*BigR**2)) * detaSp_dT * zj0**2        * xjac * theta * tstep  &
                           + v * BigR * T * (r0_corr + beta_imp*rn0_corr) * rn0_corr * dLrad_dT  * xjac * theta * tstep  &
                           + v * BigR * T * dbeta_imp_dT * rn0_corr**2 * Lrad                    * xjac * theta * tstep  &
                           + v * BigR * T * r0_corr * dfrad_bg_dT                                * xjac * theta * tstep


                 amat_67 = + v * (r0 + rn0 * alpha_imp) * F0 / BigR * Vpar * T0_p           * xjac * theta * tstep &
		           + v * T0 * F0 / BigR * Vpar * (r0_p + rn0_p * alpha_imp)         * xjac * theta * tstep &

                           + v * (r0 + rn0 * alpha_imp) * Vpar * (T0_s * ps0_t - T0_t * ps0_s)             * theta * tstep &
                           + v * T0 * Vpar * ((r0_s+rn0_s*alpha_imp)*ps0_t - (r0_t+rn0_t*alpha_imp)*ps0_s) * theta * tstep &

                           + v * (r0 + rn0 * alpha_imp) * GAMMA * T0 * (vpar_s * ps0_t - vpar_t * ps0_s)   * theta * tstep &
                           + v * (r0 + rn0 * alpha_imp) * GAMMA * T0 * F0 / BigR * vpar_p           * xjac * theta * tstep &

                           + TG_num6 * 0.25d0 / BigR * 2.d0 * vpar0*vpar &
                              * T0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                          &
                              * ( v_x * ps0_y -  v_y * ps0_x + F0 / BigR * v_p) * xjac * theta * tstep * tstep &

                           + TG_num6 * 0.25d0 / BigR * 2.d0 * vpar0*vpar &
                              * r0 * (T0_x * ps0_y - T0_y * ps0_x + F0 / BigR * T0_p)                          &
                              * ( v_x * ps0_y -  v_y * ps0_x + F0 / BigR * v_p) * xjac * theta * tstep * tstep

                 amat_68 =   v * rhon * alpha_imp * T0 * BigR * xjac * (1.d0 + zeta)                              &
!=============== The ionization potential energy term=========================
                           + v * rhon * E_ion          * BigR * xjac * (1.d0 + zeta)                              &
!================= End ionization potential energy ===========================

                           - v * rhon * BigR**2 * alpha_imp_bis * (T0_s * u0_t - T0_t * u0_s)     * theta * tstep &
                           - v * alpha_imp * T0 * BigR**2 * (rhon_s * u0_t - rhon_t * u0_s)       * theta * tstep &
                           + v * rhon * F0 / BigR * Vpar0 * alpha_imp_bis * T0_p           * xjac * theta * tstep &
                           + v * alpha_imp * T0 * F0 / BigR * Vpar0 * rhon_p               * xjac * theta * tstep &
                           + v * rhon * Vpar0 * alpha_imp_bis * (T0_s * ps0_t - T0_t * ps0_s)     * theta * tstep &
                           + v * alpha_imp * T0 * Vpar0 * (rhon_s * ps0_t - rhon_t * ps0_s)       * theta * tstep &

                           - v * alpha_imp * rhon * 2.d0* GAMMA * BigR * T0 * u0_y                   * xjac * theta * tstep &
                           + v * alpha_imp * rhon * GAMMA * T0 * (vpar0_s * ps0_t - vpar0_t * ps0_s)        * theta * tstep &
                           + v * alpha_imp * rhon * GAMMA * T0 * F0 / BigR * vpar0_p                 * xjac * theta * tstep &

                           + v * BigR * rhon * (r0_corr + 2*beta_imp*rn0_corr) * Lrad                          * xjac * theta * tstep			   

!###################################################################################################
!#  equation 7   parallel velocity equation                                                        #
!###################################################################################################

             amat_71 = v * r0 * vpar0 / BigR * (ps0_x * psi_x + ps0_y * psi_y) * xjac * (1.d0 + zeta) &

			 + v * (P0_s * psi_t - P0_t * psi_s)                                            * theta * tstep &
                         + 0.5d0 * r0 * vpar0**2 * BB2 * (psi_s * v_t - psi_t * v_s)                    * theta * tstep &
                         + 0.5d0 * v  * vpar0**2 * BB2 * (psi_s * r0_t - psi_t * r0_s)                  * theta * tstep &

                         !+ v * (particle_source(ms,mt) + source_pellet) * vpar0 * BB2_psi * BigR * xjac * theta * tstep &

                      + r0 * (vpar0_x * psi_y - vpar0_y * psi_x) &
                        * vpar0 * BB2 * v * xjac * theta * tstep  &
                      + vpar0 * (r0_x * psi_y - r0_y * psi_x) &
                        * vpar0 * BB2 * v * xjac * theta * tstep  &
                      ! Poloidal para component of the fifth term of Eq.29 here (new momentum)


                         + TG_NUM7 * 0.25d0 * r0 * Vpar0**2 * BB2 &
                                   * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac) / BigR  &
                                   * (-(psi_s * v_t     - psi_t * v_s)    /xjac)  * xjac * theta * tstep*tstep &
            
                         + TG_NUM7 * 0.25d0 * r0 * Vpar0**2 * BB2 &
                                   * (-(psi_s * vpar0_t - psi_t * vpar0_s)/xjac) / BigR  &
                                   * (-(ps0_s * v_t     - ps0_t * v_s)    /xjac)  * xjac * theta * tstep*tstep &
            
!                         + TG_NUM7 * 0.25d0 * v  * Vpar0**2 * BB2 &
!                                   * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac) / BigR  &
!                                   * (-(psi_s * r0_t    - psi_t * r0_s)   /xjac)  * xjac * theta * tstep*tstep &
            
!                         + TG_NUM7 * 0.25d0 * v  * Vpar0**2 * BB2 &
!                                   * (-(psi_s * vpar0_t - psi_t * vpar0_s)/xjac) / BigR  &
!                                   * (-(ps0_s * r0_t    - ps0_t * r0_s)   /xjac)  * xjac * theta * tstep*tstep

!=============================== New TG_num terms==================================
                       + TG_NUM7 * 0.25d0 * vpar0 * Vpar0**2 * BB2 &
                                 * (-(ps0_s * r0_t - ps0_t * r0_s)/xjac) / BigR &
                                 * (-(psi_s * v_t     - psi_t * v_s)    /xjac) * xjac * theta * tstep*tstep &

                       + TG_NUM7 * 0.25d0 * vpar0 * Vpar0**2 * BB2 &
                                 * (-(psi_s * r0_t - psi_t * r0_s)/xjac) / BigR &
                                 * (-(ps0_s * v_t     - ps0_t * v_s)    /xjac) * xjac * theta * tstep*tstep
!===============================End of new TG_num terms============================


                 amat_72 = - (r0_x_hat * u_y - r0_y_hat * u_x) * vpar0 * BB2 * v * theta * xjac* tstep
                      ! Perp component of the third term of Eq.20 here (new momentum)
 

                 amat_75 = + v * (rho_s * T0 * ps0_t - rho_t * T0 * ps0_s)                 * theta * tstep &
                           + v * (rho * T0_s * ps0_t - rho * T0_t * ps0_s)                 * theta * tstep &
                           + v * F0 / BigR * (rho_p * T0 + rho * T0_p)              * xjac * theta * tstep &

		           + 0.5d0 * rho * vpar0**2 * BB2 * (ps0_s * v_t - ps0_t * v_s)    * theta * tstep &
		           - 0.5d0 * rho * vpar0**2 * BB2 * F0 / BigR * v_p         * xjac * theta * tstep &
                           + 0.5d0 * v   * vpar0**2 * BB2 * (ps0_s * rho_t - ps0_t * rho_s)* theta * tstep &
                           - 0.5d0 * v   * vpar0**2 * BB2 * F0 / BigR * rho_p       * xjac * theta * tstep &

                       + v * rho * vpar0 * F0**2 / BigR * xjac * (1.d0 + zeta)  &
                      ! Additional zeta term for timesteping here (new momentum)

                      - (rho_x_hat * u0_y - rho_y_hat * u0_x) * vpar0 * BB2 * v * theta * xjac* tstep &
                      ! Perp component of the fifth term of Eq.29 here (new momentum)
                      + rho * F0 / BigR * vpar0_p * vpar0 * BB2 * v * theta * xjac* tstep &
                      + vpar0 * F0 / BigR * rho_p * vpar0 * BB2 * v * theta * xjac* tstep &
                      ! Toroidal para component of the fifth term of Eq.29 here (new momentum)
                      + rho * (vpar0_x * ps0_y - vpar0_y * ps0_x) &
                        * vpar0 * BB2 * v * theta * xjac* tstep  &
                      + vpar0 * (rho_x * ps0_y - rho_y * ps0_x) &
                        * vpar0 * BB2 * v * theta * xjac* tstep  &
                      ! Poloidal para component of the fifth term of Eq.29 here (new momentum)


                           + TG_NUM7 * 0.25d0 * rho * Vpar0**2 * BB2 &
                                     * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR  &
                                     * (-(ps0_s * v_t     - ps0_t * v_s)    /xjac + F0 / BigR * v_p)  * xjac * theta * tstep*tstep &
               
!                           + TG_NUM7 * 0.25d0 * v * Vpar0**2 * BB2 &
!                                     * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR  &
!                                     * (-(ps0_s * rho_t   - ps0_t * rho_s)  /xjac + F0 / BigR * rho_p)* xjac * theta * tstep*tstep 
   
!=============================== New TG_num terms==================================

                           + TG_NUM7 * 0.25d0 * vpar0 * Vpar0**2 * BB2 &
                                     * (-(ps0_s * rho_t - ps0_t * rho_s)/xjac + F0 / BigR * rho_p) / BigR  &
                                     * (-(ps0_s * v_t     - ps0_t * v_s) /xjac + F0 / BigR * v_p)  * xjac * theta * tstep*tstep

!===============================End of new TG_num terms============================



                 amat_76 = + v * (T_s * r0 * ps0_t - T_t * r0 * ps0_s)                     * theta * tstep &
                           + v * (T * r0_s * ps0_t - T * r0_t * ps0_s)                     * theta * tstep &
                           + v * F0 / BigR * (T_p * r0 + T * r0_p)                  * xjac * theta * tstep 

                 amat_77 = v * Vpar * r0 * F0**2 / BigR * xjac * (1.d0 + zeta) &
                         + visco_par * (v_x * Vpar_x + v_y * Vpar_y) * BigR        * xjac  * theta * tstep &

                         !+ v * (particle_source(ms,mt) + source_pellet)*vpar*BB2 * BigR  * xjac * theta * tstep &

                      - (r0_x_hat * u0_y - r0_y_hat * u0_x) * vpar * BB2 * v * xjac * theta * tstep &
                      ! Perp component of the fifth term of Eq.29 here (new momentum)
                      + r0 * F0 /BigR * vpar0_p * vpar * BB2 * v * xjac * theta * tstep &
                      + r0 * F0 /BigR * vpar_p * vpar0 * BB2 * v * xjac * theta * tstep &
                      + 2.0 * vpar0 * F0 / BigR * r0_p * vpar * BB2 * v * xjac * theta * tstep &
                      ! Toroidal para component of the fifth term of Eq.29 here (new momentum)
                      + r0 * (vpar_x * ps0_y - vpar_y * ps0_x) &
                        * vpar0 * BB2 * v * xjac * theta * tstep  &
                      + r0 * (vpar0_x * ps0_y - vpar0_y * ps0_x) &
                        * vpar * BB2 * v * xjac * theta * tstep  &
                      + 2.0 * vpar0 * (r0_x * ps0_y - r0_y * ps0_x) &
                        * vpar * BB2 * v * xjac * theta * tstep  &
                      ! Poloidal para component of the fifth term of Eq.29 here (new momentum)


                         + r0 * vpar0 * vpar * BB2 * (ps0_s * v_t - ps0_t * v_s)                * theta * tstep &
                         - r0 * vpar0 * vpar * BB2 * F0 / BigR * v_p                     * xjac * theta * tstep &
                         + v  * vpar0 * vpar * BB2 * (ps0_s * r0_t - ps0_t * r0_s)              * theta * tstep &
                         - v  * vpar0 * vpar * BB2 * F0 / BigR * r0_p                    * xjac * theta * tstep &
                     
                         + visco_par_num * (v_xx + v_x/BigR + v_yy)*(vpar_xx + vpar_x/BigR + vpar_yy) * BigR * xjac * theta * tstep &

                         + TG_NUM7 * 0.5d0 * r0 * Vpar * Vpar0 * BB2 &
                                   * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR                     &
                                   * (-(ps0_s * v_t     - ps0_t * v_s)    /xjac + F0 / BigR * v_p)  * xjac * theta * tstep*tstep  &
             
!                         + TG_NUM7 * 0.5d0 * v * Vpar * Vpar0 * BB2 &
!                                   * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR                     &
!                                   * (-(ps0_s * r0_t    - ps0_t * r0_s)   /xjac + F0 / BigR * r0_p)  * xjac * theta * tstep*tstep &
             
                         + TG_NUM7 * 0.25d0 * r0 * Vpar0**2 * BB2 &
                                   * (-(ps0_s * vpar_t - ps0_t * vpar_s)/xjac + F0 / BigR * vpar_p) / BigR                        &
                                   * (-(ps0_s * v_t    - ps0_t * v_s)   /xjac + F0 / BigR * v_p)  * xjac * theta * tstep*tstep    &
             
!                         + TG_NUM7 * 0.25d0 * v * Vpar0**2 * BB2 &
!                                   * (-(ps0_s * vpar_t - ps0_t * vpar_s)/xjac + F0 / BigR * vpar_p) / BigR                        &
!                                   * (-(ps0_s * r0_t   - ps0_t * r0_s)  /xjac + F0 / BigR * r0_p)  * xjac * theta * tstep*tstep   

!=============================== New TG_num terms==================================

            + TG_NUM7 * 0.75d0 * Vpar * Vpar0**2 * BB2 &
                      * (-(ps0_s * r0_t - ps0_t * r0_s)/xjac + F0 / BigR * r0_p) / BigR             &
                      * (-(ps0_s * v_t  - ps0_t * v_s) /xjac + F0 / BigR * v_p)  * xjac * theta * tstep*tstep

!===============================End of new TG_num terms============================


                amat_78 = 0 ! Place holder

!################################################################################################### 
!#  equation 8   impurity density equation                                                          # 
!################################################################################################### 

                amat_81 =  + v * rn0 * (vpar0_s * psi_t - vpar0_t * psi_s)                                * theta * tstep &
                           + v * Vpar0 * (rn0_s * psi_t - rn0_t * psi_s)                                  * theta * tstep &
                           ! New diffusion scheme for impurities
                           - (D_par-D_prof) * BigR * BB2_psi / BB2**2 * Bgrad_rho_star * (Bgrad_rhon)                     * xjac * theta * tstep &
                           + (D_par-D_prof) * BigR / BB2              * Bgrad_rho_star_psi * (Bgrad_rhon)                 * xjac * theta * tstep &
                           + (D_par-D_prof) * BigR / BB2              * Bgrad_rho_star * (Bgrad_rhon_psi)                 * xjac * theta * tstep &
                           + TG_num8 * 0.25d0 / BigR * vpar0**2                                                           &
                                     * (rn0_x * psi_y - rn0_y * psi_x)                                                    &
                                     * ( v_x * ps0_y -  v_y * ps0_x + F0 / BigR * v_p) * xjac * theta * tstep * tstep     &

                           + TG_num8 * 0.25d0 / BigR * vpar0**2                                                           &
                                     * (rn0_x * ps0_y - rn0_y * ps0_x + F0 / BigR * rn0_p)                                &
                                     * ( v_x * psi_y -  v_y * psi_x                   ) * xjac * theta * tstep * tstep

                amat_82 =  + v * BigR**2 * ( rn0_s * u_t - rn0_t * u_s)                                   * theta * tstep &
                                              + v * 2.d0 * BigR * rn0 * u_y                        * xjac * theta * tstep &

                           + TG_num8 * 0.25d0 * BigR**3 * (rn0_x * u_y  - rn0_y * u_x)                                    &
                                                        * ( v_x * u0_y - v_y  * u0_x) * xjac * theta * tstep * tstep      &

                           + TG_num8 * 0.25d0 * BigR**3 * (rn0_x * u0_y - rn0_y * u0_x)                                   &
                                                        * ( v_x * u_y  - v_y  * u_x)  * xjac * theta * tstep * tstep 

               ! We do not include the term coming from div(rhon * v_star_i) because they are prop. to rho_n/rho, because they may cause problems
               ! in areas where rho is small.                   

                amat_85 = 0 ! Place holder
                amat_86 = 0 ! Place holder                 

                amat_87 = + v * F0 / BigR * Vpar * rn0_p                                           * xjac * theta * tstep &
                          + v * Vpar * (rn0_s * ps0_t - rn0_t * ps0_s)                                    * theta * tstep &
                          + v * rn0 * (vpar_s * ps0_t - vpar_t * ps0_s)                                   * theta * tstep &
                          + v * rn0 * F0 / BigR * vpar_p                                           * xjac * theta * tstep &                  

                           + TG_num8 * 0.25d0 / BigR * 2.d0*vpar0*vpar                                                    &
                              * (rn0_x * ps0_y - rn0_y * ps0_x + F0 / BigR * rn0_p)                                       &
                              * ( v_x * ps0_y -  v_y * ps0_x + F0 / BigR * v_p) * xjac * theta * tstep * tstep 

		amat_88 = + v * rhon * BigR * xjac * (1.d0 + zeta)   &
                          - v * BigR**2 * ( rhon_s * u0_t - rhon_t * u0_s)                                * theta * tstep &
                          - v * 2.d0 * BigR * rhon * u0_y                                          * xjac * theta * tstep &
                          + v * F0 / BigR * Vpar0 * rhon_p                                         * xjac * theta * tstep &
                          + v * Vpar0 * (rhon_s * ps0_t - rhon_t * ps0_s)                                 * theta * tstep &
                          + v * rhon * (vpar0_s * ps0_t - vpar0_t * ps0_s)                                * theta * tstep &
                          + v * rhon * F0 / BigR  * vpar0_p                                        * xjac * theta * tstep &

                         + TG_num8 * 0.25d0 * BigR**3 * (rhon_x * u0_y - rhon_y * u0_x)                                   &
                                                      * ( v_x  * u0_y - v_y   * u0_x) * xjac * theta * tstep * tstep      &

                         + TG_num8 * 0.25d0 / BigR * vpar0**2                                                             &
                              * (rhon_x * ps0_y - rhon_y * ps0_x + F0 / BigR * rhon_p)                                    &
                              * ( v_x * ps0_y -  v_y * ps0_x   + F0 / BigR * v_p) * xjac * theta * tstep * tstep          &
                          !New diffusion scheme for impurities
                         + (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star * Bgrad_rho_rhon                  * xjac * theta * tstep &
                         + D_prof * BigR  * (v_x*rhon_x + v_y*rhon_y + v_p*rhon_p * eps_cyl**2 /BigR**2 ) * xjac * theta * tstep &
         
		          + BigR * (Dn0x * rhon_x * v_x + Dn0y * rhon_y * v_y + Dn0p * rhon_p * v_p*eps_cyl**2/BigR**2) * xjac * theta * tstep  & 

                          + Dn_perp_num     * (v_xx + v_x/BigR + v_yy)*(rhon_xx + rhon_x/BigR + rhon_yy)   * BigR * xjac * theta * tstep 


                 kl1 = index_kl
                 kl2 = index_kl + 1*n_tor
                 kl3 = index_kl + 2*n_tor
                 kl4 = index_kl + 3*n_tor
                 kl5 = index_kl + 4*n_tor
                 kl6 = index_kl + 5*n_tor
                 kl7 = index_kl + 6*n_tor
                 kl8 = index_kl + 7*n_tor
		 
                 ELM(ij1,kl1) =  ELM(ij1,kl1) + wst * amat_11
                 ELM(ij1,kl2) =  ELM(ij1,kl2) + wst * amat_12
                 ELM(ij1,kl3) =  ELM(ij1,kl3) + wst * amat_13
                 ELM(ij1,kl5) =  ELM(ij1,kl5) + wst * amat_15
                 ELM(ij1,kl6) =  ELM(ij1,kl6) + wst * amat_16

                 ELM(ij2,kl1) =  ELM(ij2,kl1) + wst * amat_21
                 ELM(ij2,kl2) =  ELM(ij2,kl2) + wst * amat_22
                 ELM(ij2,kl3) =  ELM(ij2,kl3) + wst * amat_23
                 ELM(ij2,kl4) =  ELM(ij2,kl4) + wst * amat_24
                 ELM(ij2,kl5) =  ELM(ij2,kl5) + wst * amat_25
                 ELM(ij2,kl6) =  ELM(ij2,kl6) + wst * amat_26
                 ELM(ij2,kl7) =  ELM(ij2,kl7) + wst * amat_27 ! New term due to the new momentum eq.
		 ELM(ij2,kl8) =  ELM(ij2,kl8) + wst * amat_28

                 ELM(ij3,kl1) =  ELM(ij3,kl1) + wst * amat_31
                 ELM(ij3,kl3) =  ELM(ij3,kl3) + wst * amat_33

                 ELM(ij4,kl2) =  ELM(ij4,kl2) + wst * amat_42
                 ELM(ij4,kl4) =  ELM(ij4,kl4) + wst * amat_44

                 ELM(ij5,kl1) =  ELM(ij5,kl1) + wst * amat_51
                 ELM(ij5,kl2) =  ELM(ij5,kl2) + wst * amat_52
                 ELM(ij5,kl5) =  ELM(ij5,kl5) + wst * amat_55
                 ELM(ij5,kl6) =  ELM(ij5,kl6) + wst * amat_56
                 ELM(ij5,kl7) =  ELM(ij5,kl7) + wst * amat_57
		 ELM(ij5,kl8) =  ELM(ij5,kl8) + wst * amat_58
                 
		 ELM(ij6,kl1) =  ELM(ij6,kl1) + wst * amat_61
                 ELM(ij6,kl2) =  ELM(ij6,kl2) + wst * amat_62
                 ELM(ij6,kl5) =  ELM(ij6,kl5) + wst * amat_65
                 ELM(ij6,kl6) =  ELM(ij6,kl6) + wst * amat_66
                 ELM(ij6,kl7) =  ELM(ij6,kl7) + wst * amat_67
		 ELM(ij6,kl8) =  ELM(ij6,kl8) + wst * amat_68

                 ELM(ij7,kl1) =  ELM(ij7,kl1) + wst * amat_71
                 ELM(ij7,kl2) =  ELM(ij7,kl2) + wst * amat_72
                 ELM(ij7,kl5) =  ELM(ij7,kl5) + wst * amat_75
                 ELM(ij7,kl6) =  ELM(ij7,kl6) + wst * amat_76
                 ELM(ij7,kl7) =  ELM(ij7,kl7) + wst * amat_77
		 ELM(ij7,kl8) =  ELM(ij7,kl8) + wst * amat_78

                 ELM(ij8,kl1) =  ELM(ij8,kl1) + wst * amat_81
                 ELM(ij8,kl2) =  ELM(ij8,kl2) + wst * amat_82		 
		 ELM(ij8,kl5) =  ELM(ij8,kl5) + wst * amat_85
                 ELM(ij8,kl6) =  ELM(ij8,kl6) + wst * amat_86
                 ELM(ij8,kl7) =  ELM(ij8,kl7) + wst * amat_87
                 ELM(ij8,kl8) =  ELM(ij8,kl8) + wst * amat_88

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
