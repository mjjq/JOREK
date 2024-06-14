module mod_elt_matrix_fft
  implicit none
contains

#include "corr_neg_include.f90"

subroutine element_matrix_fft(element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, ELM, RHS, tid, &
                              ELM_p, ELM_n, ELM_k, ELM_kn, RHS_p, RHS_k,  eq_g, eq_s, eq_t, eq_p, eq_ss, eq_st, eq_tt, delta_g, delta_s, delta_t, & 
                              i_tor_min, i_tor_max, aux_nodes, ELM_pnn, ELM_knn)
! NOT YET IMPLEMENTED

use mod_parameters
use data_structure
use gauss
use basis_at_gaussian
use phys_module
use diffusivities, only: get_dperp, get_zk_iperp, get_zk_eperp
use equil_info, only : get_psi_n, ES
use mod_F_profile
use mod_bootstrap_functions
use mod_sources

implicit none

! --- Input Variables
type (type_element)        :: element
type (type_node)           :: nodes(n_vertex_max)
type (type_node), optional :: aux_nodes(n_vertex_max)

logical, intent(in)    :: xpoint2
integer, intent(in)    :: xcase2
real*8,  intent(in)    :: R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint(2), Z_xpoint(2)

#define DIM0 n_tor*n_vertex_max*n_degrees*n_var
#define DIM1 n_plane
#define DIM2 1:n_vertex_max*n_var*n_degrees

#define ENABLE_VISCO_PAR_AND_PARPAR
#define ENABLE_DIA_DRIFT_TEMPERATURE
#define MOMENTUM_UPHI_IN_PARALLEL
#define TEMPERATURE_DEPENDENT_UPDATE
#define VISCOSITY_TENSOR_FORM

#ifndef VISCOSITY_TENSOR_FORM
#define VISCOSITY_INSIDE_LAPLACIAN
#define VISCOSITY_COMPLETE_TERMS
#endif

real*8, dimension (DIM0,DIM0)  :: ELM
real*8, dimension (DIM0)       :: RHS
integer, intent(in)            :: tid, i_tor_min, i_tor_max
real*8, dimension(DIM1, DIM2, DIM2), intent(inout) :: ELM_p
real*8, dimension(DIM1, DIM2, DIM2), intent(inout) :: ELM_n
real*8, dimension(DIM1, DIM2, DIM2), intent(inout) :: ELM_k
real*8, dimension(DIM1, DIM2, DIM2), intent(inout) :: ELM_kn
real*8, dimension(DIM1, DIM2),       intent(inout) :: RHS_p
real*8, dimension(DIM1, DIM2),       intent(inout) :: RHS_k

real*8, dimension(n_plane,n_var,n_gauss,n_gauss), intent(inout) :: eq_g, eq_s, eq_t
real*8, dimension(n_plane,n_var,n_gauss,n_gauss), intent(inout) :: eq_p
real*8, dimension(n_plane,n_var,n_gauss,n_gauss), intent(inout) :: eq_ss, eq_st, eq_tt
real*8, dimension(n_plane,n_var,n_gauss,n_gauss), intent(inout) :: delta_g, delta_s, delta_t

! --- Variables outside the OMP loop
integer    :: n_tor_start, n_tor_end, n_tor_local
integer    :: i, j, index, index_k, index_m, i_v, j_loc, i_loc, m, ik
integer    :: in, im, ivar, kvar, ms, mt, mp
real*8     :: wst, xjac, xjac_R, xjac_Z, R, Z, theta, zeta
real*8     :: current_source_JR(n_gauss,n_gauss), current_source_JZ(n_gauss,n_gauss), current_source_Jp(n_gauss,n_gauss)
real*8     :: particle_source(n_gauss,n_gauss),heat_source_i(n_gauss,n_gauss),heat_source_e(n_gauss,n_gauss)
real*8     :: in_fft(1:n_plane)
complex*16 :: out_fft(1:n_plane)
real*8     :: rho_axisym(n_gauss,n_gauss), Ti_axisym(n_gauss,n_gauss), Te_axisym(n_gauss,n_gauss)
real*8     :: rho_GS,Ti_GS,Te_GS, dd1,dd2,dd3,dd4,dd5,dd6,dd7,dd8
real*8     :: psi_axisym(n_gauss,n_gauss), psi_axisym_s(n_gauss,n_gauss), psi_axisym_t(n_gauss,n_gauss)
real*8     ::                              psi_axisym_R(n_gauss,n_gauss), psi_axisym_Z(n_gauss,n_gauss)
real*8     :: Fprof_time_dep,dF_dpsi(n_gauss,n_gauss)      ,dF_dz      ,dF_dpsi2      ,dF_dz2      ,dF_dpsi_dz
real*8     :: zFFprime      ,dFFprime_dpsi,dFFprime_dz,dFFprime_dpsi2,dFFprime_dz2,dFFprime_dpsi_dz
real*8     :: rho_initial(n_gauss,n_gauss),dn_dpsi, dn_dz, dn_dpsi2, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi_dz2,  dn_dpsi2_dz
real*8     :: Ti_initial (n_gauss,n_gauss),dTi_dpsi,dTi_dz,dTi_dpsi2,dTi_dz2,dTi_dpsi_dz,dTi_dpsi3,dTi_dpsi_dz2, dTi_dpsi2_dz
real*8     :: Te_initial (n_gauss,n_gauss),dTe_dpsi,dTe_dz,dTe_dpsi2,dTe_dz2,dTe_dpsi_dz,dTe_dpsi3,dTe_dpsi_dz2, dTe_dpsi2_dz
real*8     :: Jb, Jb_0
real*8     :: amu_neo_prof(n_gauss,n_gauss), aki_neo_prof(n_gauss,n_gauss)
real*8     :: V_source(n_gauss,n_gauss), dV_dpsi_source(n_gauss,n_gauss),dV_dz_source(n_gauss,n_gauss)
real*8     :: dV_dpsi2,dV_dz2,dV_dpsi_dz,dV_dpsi3,dV_dpsi_dz2, dV_dpsi2_dz
real*8     :: eta_ARAZ, tauIC_ARAZ, tauIC_AP, tauIC_rho, tauIC_tem
logical    :: use_fft


real*8, dimension(n_gauss,n_gauss)    :: x_g, x_s, x_t, x_ss, x_st, x_tt
real*8, dimension(n_gauss,n_gauss)    :: y_g, y_s, y_t, y_ss, y_st, y_tt
real*8, dimension(n_gauss,n_gauss)    :: Fprofile
real*8, dimension(n_var)              :: TG_NUM
real*8, dimension(n_tor,n_plane) :: HHZ, HHZ_p, HHZ_pp

#ifndef TEMPERATURE_DEPENDENT_UPDATE
real*8     :: Te_GS_initial(n_gauss,n_gauss), dTe_GS_dR_initial(n_gauss,n_gauss), dTe_GS_dZ_initial(n_gauss,n_gauss)
real*8     :: Ti_GS_initial(n_gauss,n_gauss), dTi_GS_dR_initial(n_gauss,n_gauss), dTi_GS_dZ_initial(n_gauss,n_gauss)
#endif

! --- Variables inside the OMP loop
integer    :: k, l, index_ij, index_kl, ij, kl, sgn_l

real*8     :: AR0,  AR0_R,  AR0_Z,  AR0_p,  AR0_s,  AR0_t, AR0_ss, AR0_tt, AR0_st, AR0_RR, AR0_ZZ, AR0_RZ, AR0_pp
real*8     :: AZ0,  AZ0_R,  AZ0_Z,  AZ0_p,  AZ0_s,  AZ0_t, AZ0_ss, AZ0_tt, AZ0_st, AZ0_RR, AZ0_ZZ, AZ0_RZ, AZ0_pp
real*8     :: A30,  A30_R,  A30_Z,  A30_p,  A30_s,  A30_t, A30_ss, A30_tt, A30_st, A30_RR, A30_ZZ, A30_RZ, A30_pp
real*8     :: UR0,  UR0_R,  UR0_Z,  UR0_p,  UR0_s,  UR0_t, UR0_ss, UR0_st, UR0_tt, UR0_RR, UR0_ZZ, UR0_RZ, UR0_pp
real*8     :: UZ0,  UZ0_R,  UZ0_Z,  UZ0_p,  UZ0_s,  UZ0_t, UZ0_ss, UZ0_st, UZ0_tt, UZ0_RR, UZ0_ZZ, UZ0_RZ, UZ0_pp
real*8     :: Up0,  Up0_R,  Up0_Z,  Up0_p,  Up0_s,  Up0_t, Up0_ss, Up0_st, Up0_tt, Up0_RR, Up0_ZZ, Up0_RZ, Up0_pp
real*8     :: rho0, rho0_R, rho0_Z, rho0_p, rho0_s, rho0_t,rho0_ss,rho0_st,rho0_tt,rho0_RR,rho0_ZZ,rho0_RZ,rho0_pp, rho0_corr
real*8     :: Ti0,  Ti0_R,  Ti0_Z,  Ti0_p,  Ti0_s,  Ti0_t, Ti0_ss, Ti0_st, Ti0_tt, Ti0_RR, Ti0_ZZ, Ti0_RZ, Ti0_pp,  Ti0_corr
real*8     :: Te0,  Te0_R,  Te0_Z,  Te0_p,  Te0_s,  Te0_t, Te0_ss, Te0_st, Te0_tt, Te0_RR, Te0_ZZ, Te0_RZ, Te0_pp,  Te0_corr
real*8     :: pi0,  pi0_R,  pi0_Z,  pi0_p,  pi0_s,  pi0_t,  pi0_corr
real*8     :: pe0,  pe0_R,  pe0_Z,  pe0_p,  pe0_s,  pe0_t,  pe0_corr
real*8     :: p0,   p0_R,   p0_Z,   p0_p,   p0_s,   p0_t,   p0_corr

real*8     :: AR,  AR_R,  AR_Z,  AR_p__n,  AR_s,  AR_t
real*8     :: AZ,  AZ_R,  AZ_Z,  AZ_p__n,  AZ_s,  AZ_t
real*8     :: A3,  A3_R,  A3_Z,  A3_p__n,  A3_s,  A3_t
real*8     :: UR,  UR_R,  UR_Z,  UR_p__n,  UR_s,  UR_t, UR_RR, UR_ZZ, UR_pp
real*8     :: UZ,  UZ_R,  UZ_Z,  UZ_p__n,  UZ_s,  UZ_t, UZ_RR, UZ_ZZ, UZ_pp
real*8     :: Up,  Up_R,  Up_Z,  Up_p__n,  Up_s,  Up_t, Up_RR, Up_ZZ, Up_pp
real*8     :: Ti,  Ti_R,  Ti_Z,  Ti_p__n,  Ti_s,  Ti_t, Ti_RR, Ti_ZZ, Ti_pp
real*8     :: Te,  Te_R,  Te_Z,  Te_p__n,  Te_s,  Te_t, Te_RR, Te_ZZ, Te_pp
real*8     :: rho, rho_R, rho_Z, rho_p__n, rho_s, rho_t, rho_RR, rho_ZZ, rho_pp

real*8     :: v,  v_R,  v_Z,  v_s,  v_t,  v_p__k,  v_ss,  v_st,  v_tt,  v_RR,  v_ZZ
real*8     :: bf, bf_R, bf_Z, bf_s, bf_t, bf_p, bf_ss, bf_st, bf_tt, bf_RR, bf_ZZ

real*8     :: Fprof
real*8     :: BR0, BR0_AR,    BR0_AZ__n, BR0_A3
real*8     :: BZ0, BZ0_AR__n, BZ0_AZ,    BZ0_A3
real*8     :: Bp0, Bp0_AR,    Bp0_AZ,    Bp0_A3
real*8     :: BB2, BB2_AR__p, BB2_AR__n, BB2_AZ__p, BB2_AZ__n, BB2_A3
real*8     :: Bp00

real*8     :: BgradTi, BgradTi_AR__p, BgradTi_AR__n, BgradTi_AZ__p, BgradTi_AZ__n, BgradTi_A3, BgradTi_Ti__p, BgradTi_Ti__n
real*8     :: BgradTe, BgradTe_AR__p, BgradTe_AR__n, BgradTe_AZ__p, BgradTe_AZ__n, BgradTe_A3, BgradTe_Te__p, BgradTe_Te__n

real*8     :: BgradRho, BgradRho_AR__p, BgradRho_AR__n, BgradRho_AZ__p, BgradRho_AZ__n, BgradRho_A3, BgradRho_rho__p, BgradRho_rho__n

real*8     :: BgradPe, BgradPe_AR__p, BgradPe_AR__n, BgradPe_AZ__p, BgradPe_AZ__n, BgradPe_A3, BgradPe_Te__p, BgradPe_Te__n, BgradPe_rho__p, BgradPe_rho__n

real*8     :: BgradVstar__p, BgradVstar__k
real*8     :: BgradVstar_AR__p, BgradVstar_AR__k, BgradVstar_AR__n
real*8     :: BgradVstar_AZ__p, BgradVstar_AZ__k, BgradVstar_AZ__n
real*8     :: BgradVstar_A3__p, BgradVstar_A3__k

real*8     :: UgradRho, UgradRho_UR, UgradRho_UZ, UgradRho_Up, UgradRho_rho__p, UgradRho_rho__n
real*8     :: UgradTi,  UgradTi_UR,  UgradTi_UZ,  UgradTi_Up,  UgradTi_Ti__p,   UgradTi_Ti__n
real*8     :: UgradTe,  UgradTe_UR,  UgradTe_UZ,  UgradTe_Up,  UgradTe_Te__p,   UgradTe_Te__n

real*8     :: UgradVstar__p, UgradVstar__k, UgradVstar_UR, UgradVstar_UZ, UgradVstar_Up__k

real*8     :: gradRho_gradVstar__p, gradRho_gradVstar__k, gradRho_gradVstar_rho__p, gradRho_gradVstar_rho__kn
real*8     :: gradTi_gradVstar__p,  gradTi_gradVstar__k,  gradTi_gradVstar_Ti__p,   gradTi_gradVstar_Ti__kn
real*8     :: gradTe_gradVstar__p,  gradTe_gradVstar__k,  gradTe_gradVstar_Te__p,   gradTe_gradVstar_Te__kn

real*8     :: UgradUR, UgradUR_UR__p, UgradUR_UR__n, UgradUR_UZ,    UgradUR_Up
real*8     :: UgradUZ, UgradUZ_UR,    UgradUZ_UZ__p, UgradUZ_UZ__n, UgradUZ_Up
real*8     :: UgradUp, UgradUp_UR,    UgradUp_UZ   , UgradUp_Up__p, UgradUp_Up__n

real*8     :: divU, divU_UR, divU_UZ, divU_Up__n

real*8     :: divRhoU, divRhoU_UR, divRhoU_UZ, divRhoU_Up__p, divRhoU_Up__n, divRhoU_rho__p, divRhoU_rho__n

real*8     :: Vt0, Vt0_R, Vt0_Z

real*8     :: tau_IC, distance_bnd

real*8     :: VdiaR0
real*8     :: VdiaR0_AR__p,  VdiaR0_AR__n
real*8     :: VdiaR0_AZ__p,  VdiaR0_AZ__n
real*8     :: VdiaR0_A3__p,  VdiaR0_A3__n
real*8     :: VdiaR0_rho__p, VdiaR0_rho__n
real*8     :: VdiaR0_Ti__p,  VdiaR0_Ti__n

real*8     :: VdiaZ0
real*8     :: VdiaZ0_AR__p,  VdiaZ0_AR__n
real*8     :: VdiaZ0_AZ__p,  VdiaZ0_AZ__n
real*8     :: VdiaZ0_A3__p,  VdiaZ0_A3__n
real*8     :: VdiaZ0_rho__p, VdiaZ0_rho__n
real*8     :: VdiaZ0_Ti__p,  VdiaZ0_Ti__n

real*8     :: VdiaP0
real*8     :: VdiaP0_AR__p,  VdiaP0_AR__n
real*8     :: VdiaP0_AZ__p,  VdiaP0_AZ__n
real*8     :: VdiaP0_A3__p,  VdiaP0_A3__n
real*8     :: VdiaP0_rho__p, VdiaP0_rho__n
real*8     :: VdiaP0_Ti__p,  VdiaP0_Ti__n

real*8     :: VdiaGradUR
real*8     :: VdiaGradUR_AR__p,  VdiaGradUR_AR__n
real*8     :: VdiaGradUR_AZ__p,  VdiaGradUR_AZ__n
real*8     :: VdiaGradUR_A3__p,  VdiaGradUR_A3__n
real*8     :: VdiaGradUR_rho__p, VdiaGradUR_rho__n
real*8     :: VdiaGradUR_Ti__p,  VdiaGradUR_Ti__n
real*8     :: VdiaGradUR_UR__p,  VdiaGradUR_UR__n

real*8     :: VdiaGradUZ
real*8     :: VdiaGradUZ_AR__p,  VdiaGradUZ_AR__n
real*8     :: VdiaGradUZ_AZ__p,  VdiaGradUZ_AZ__n
real*8     :: VdiaGradUZ_A3__p,  VdiaGradUZ_A3__n
real*8     :: VdiaGradUZ_rho__p, VdiaGradUZ_rho__n
real*8     :: VdiaGradUZ_Ti__p,  VdiaGradUZ_Ti__n
real*8     :: VdiaGradUZ_UZ__p,  VdiaGradUZ_UZ__n

real*8     :: VdiaGradUp
real*8     :: VdiaGradUp_AR__p,  VdiaGradUp_AR__n
real*8     :: VdiaGradUp_AZ__p,  VdiaGradUp_AZ__n
real*8     :: VdiaGradUp_A3__p,  VdiaGradUp_A3__n
real*8     :: VdiaGradUp_rho__p, VdiaGradUp_rho__n
real*8     :: VdiaGradUp_Ti__p,  VdiaGradUp_Ti__n
real*8     :: VdiaGradUp_Up__p,  VdiaGradUp_Up__n

real*8     :: VdiaGradVstar__p, VdiaGradVstar__k
real*8     :: VdiaGradVstar_AR__p,  VdiaGradVstar_AR__n,  VdiaGradVstar_AR__k,  VdiaGradVstar_AR__kn
real*8     :: VdiaGradVstar_AZ__p,  VdiaGradVstar_AZ__n,  VdiaGradVstar_AZ__k,  VdiaGradVstar_AZ__kn
real*8     :: VdiaGradVstar_A3__p,  VdiaGradVstar_A3__n,  VdiaGradVstar_A3__k,  VdiaGradVstar_A3__kn
real*8     :: VdiaGradVstar_rho__p, VdiaGradVstar_rho__n, VdiaGradVstar_rho__k, VdiaGradVstar_rho__kn
real*8     :: VdiaGradVstar_Ti__p,  VdiaGradVstar_Ti__n,  VdiaGradVstar_Ti__k,  VdiaGradVstar_Ti__kn

real*8     :: Btht, Btht_AR__p, Btht_AR__n, Btht_AZ__p, Btht_AZ__n, Btht_A3
real*8     :: Vtht
real*8     :: Vtht_AR__p,  Vtht_AR__n
real*8     :: Vtht_AZ__p,  Vtht_AZ__n
real*8     :: Vtht_A3__p,  Vtht_A3__n
real*8     :: Vtht_rho__p, Vtht_rho__n
real*8     :: Vtht_Ti__p,  Vtht_Ti__n
real*8     :: Vtht_UR, Vtht_UZ
real*8     :: Vneo
real*8     :: Vneo_AR__p,  Vneo_AR__n
real*8     :: Vneo_AZ__p,  Vneo_AZ__n
real*8     :: Vneo_A3__p,  Vneo_A3__n
real*8     :: Vneo_Ti__p,  Vneo_Ti__n
real*8     :: PneoR
real*8     :: PneoR_AR__p,  PneoR_AR__n
real*8     :: PneoR_AZ__p,  PneoR_AZ__n
real*8     :: PneoR_A3__p,  PneoR_A3__n
real*8     :: PneoR_rho__p, PneoR_rho__n
real*8     :: PneoR_Ti__p,  PneoR_Ti__n
real*8     :: PneoR_UR,     PneoR_UZ
real*8     :: PneoZ
real*8     :: PneoZ_AR__p,  PneoZ_AR__n
real*8     :: PneoZ_AZ__p,  PneoZ_AZ__n
real*8     :: PneoZ_A3__p,  PneoZ_A3__n
real*8     :: PneoZ_rho__p, PneoZ_rho__n
real*8     :: PneoZ_Ti__p,  PneoZ_Ti__n
real*8     :: PneoZ_UR,     PneoZ_UZ

real*8     :: ZKi_prof, ZKe_prof, D_par_prof, D_prof, psi_norm

real*8     :: eta_T, visco_T, deta_dT, d2eta_d2T, dvisco_dT, visco_num_T
real*8     :: visco_R, visco_Z, visco_p, d2visco_d2T
real*8     :: visco_T_T, visco_R_T, visco_Z_T, visco_p_T__p, visco_p_T__n ! basis functions
real*8     :: eta_num_T, eta_R, eta_Z, eta_p, Zki_par_T, dZKi_par_dT, Zke_par_T, dZKe_par_dT
real*8     :: eta_T_T, eta_R_T, eta_Z_T, eta_p_T__p, eta_p_T__n
real*8     :: eta_T_ohm, deta_dT_ohm

real*8     :: Qconv_UR
real*8     :: Qconv_UR_AR__p,  Qconv_UR_AR__n
real*8     :: Qconv_UR_AZ__p,  Qconv_UR_AZ__n
real*8     :: Qconv_UR_A3__p,  Qconv_UR_A3__n
real*8     :: Qconv_UR_UR__p,  Qconv_UR_UR__n
real*8     :: Qconv_UR_UZ__p,  Qconv_UR_UZ__n
real*8     :: Qconv_UR_Up__p,  Qconv_UR_Up__n
real*8     :: Qconv_UR_rho__p, Qconv_UR_rho__n
real*8     :: Qconv_UR_Ti__p,  Qconv_UR_Ti__n

real*8     :: Qconv_UZ
real*8     :: Qconv_UZ_AR__p,  Qconv_UZ_AR__n
real*8     :: Qconv_UZ_AZ__p,  Qconv_UZ_AZ__n
real*8     :: Qconv_UZ_A3__p,  Qconv_UZ_A3__n
real*8     :: Qconv_UZ_UR__p,  Qconv_UZ_UR__n
real*8     :: Qconv_UZ_UZ__p,  Qconv_UZ_UZ__n
real*8     :: Qconv_UZ_Up__p,  Qconv_UZ_Up__n
real*8     :: Qconv_UZ_rho__p, Qconv_UZ_rho__n
real*8     :: Qconv_UZ_Ti__p,  Qconv_UZ_Ti__n

real*8     :: Qconv_Up
real*8     :: Qconv_Up_AR__p,  Qconv_Up_AR__n
real*8     :: Qconv_Up_AZ__p,  Qconv_Up_AZ__n
real*8     :: Qconv_Up_A3__p,  Qconv_Up_A3__n
real*8     :: Qconv_Up_UR__p,  Qconv_Up_UR__n
real*8     :: Qconv_Up_UZ__p,  Qconv_Up_UZ__n
real*8     :: Qconv_Up_Up__p,  Qconv_Up_Up__n
real*8     :: Qconv_Up_rho__p, Qconv_Up_rho__n
real*8     :: Qconv_Up_Ti__p,  Qconv_Up_Ti__n

real*8     :: JxB_UR__p, JxB_UR__k
real*8     :: JxB_UZ__p, JxB_UZ__k
real*8     :: JxB_Up__p, JxB_Up__k
real*8     :: JxB_UR_AR__p, JxB_UR_AR__k, JxB_UR_AR__n, JxB_UR_AR__kn
real*8     :: JxB_UZ_AR__p, JxB_UZ_AR__k, JxB_UZ_AR__n, JxB_UZ_AR__kn
real*8     :: JxB_Up_AR__p, JxB_Up_AR__k, JxB_Up_AR__n, JxB_Up_AR__kn
real*8     :: JxB_UR_AZ__p, JxB_UR_AZ__k, JxB_UR_AZ__n, JxB_UR_AZ__kn
real*8     :: JxB_UZ_AZ__p, JxB_UZ_AZ__k, JxB_UZ_AZ__n, JxB_UZ_AZ__kn
real*8     :: JxB_Up_AZ__p, JxB_Up_AZ__k, JxB_Up_AZ__n, JxB_Up_AZ__kn
real*8     :: JxB_UR_A3__p, JxB_UR_A3__k, JxB_UR_A3__n, JxB_UR_A3__kn
real*8     :: JxB_UZ_A3__p, JxB_UZ_A3__k, JxB_UZ_A3__n, JxB_UZ_A3__kn
real*8     :: JxB_Up_A3__p, JxB_Up_A3__k, JxB_Up_A3__n, JxB_Up_A3__kn

real*8     :: Qvisc_UR__p, Qvisc_UR__k
real*8     :: Qvisc_UR_UR__p, Qvisc_UR_UR__k, Qvisc_UR_UR__n, Qvisc_UR_UR__kn
real*8     :: Qvisc_UR_UZ__p, Qvisc_UR_UZ__k, Qvisc_UR_UZ__n, Qvisc_UR_UZ__kn
real*8     :: Qvisc_UR_Up__p, Qvisc_UR_Up__k, Qvisc_UR_Up__n, Qvisc_UR_Up__kn

real*8     :: Qvisc_UZ__p, Qvisc_UZ__k
real*8     :: Qvisc_UZ_UR__p, Qvisc_UZ_UR__k, Qvisc_UZ_UR__n, Qvisc_UZ_UR__kn
real*8     :: Qvisc_UZ_UZ__p, Qvisc_UZ_UZ__k, Qvisc_UZ_UZ__n, Qvisc_UZ_UZ__kn
real*8     :: Qvisc_UZ_Up__p, Qvisc_UZ_Up__k, Qvisc_UZ_Up__n, Qvisc_UZ_Up__kn

real*8     :: Qvisc_Up__p, Qvisc_Up__k
real*8     :: Qvisc_Up_UR__p, Qvisc_Up_UR__k, Qvisc_Up_UR__n, Qvisc_Up_UR__kn
real*8     :: Qvisc_Up_UZ__p, Qvisc_Up_UZ__k, Qvisc_Up_UZ__n, Qvisc_Up_UZ__kn
real*8     :: Qvisc_Up_Up__p, Qvisc_Up_Up__k, Qvisc_Up_Up__n, Qvisc_Up_Up__kn

real*8     :: Qgrad_visc_dot_grad_UR__UR__p, Qgrad_visc_dot_grad_UR__UR__n
real*8     :: Qgrad_visc_dot_grad_UZ__UZ__p, Qgrad_visc_dot_grad_UZ__UZ__n
real*8     :: Qgrad_visc_dot_grad_Up__Up__p, Qgrad_visc_dot_grad_Up__Up__n
real*8     :: Qgrad_visc_dot_grad_UR__Te__p, Qgrad_visc_dot_grad_UR__Te__n
real*8     :: Qgrad_visc_dot_grad_UZ__Te__p, Qgrad_visc_dot_grad_UZ__Te__n
real*8     :: Qgrad_visc_dot_grad_Up__Te__p, Qgrad_visc_dot_grad_Up__Te__n

real*8     :: Qgrad_visc_dot_grad_UR__UZ__p, Qgrad_visc_dot_grad_UR__Up__p, Qgrad_visc_dot_grad_UR__Up__n
real*8     :: Qgrad_visc_dot_grad_UZ__UR__p, Qgrad_visc_dot_grad_UZ__Up__p, Qgrad_visc_dot_grad_UZ__Up__n
real*8     :: Qgrad_visc_dot_grad_Up__UR__p, Qgrad_visc_dot_grad_Up__UR__n, Qgrad_visc_dot_grad_Up__UZ__p, Qgrad_visc_dot_grad_Up__UZ__n

real*8     :: QdivU_BdotGradVisc__UR__p, QdivU_BdotGradVisc__UZ__p, QdivU_BdotGradVisc__Up__n
real*8     :: QdivU_BdotGradVisc__AR__p
real*8     :: QdivU_BdotGradVisc__AR__n
real*8     :: QdivU_BdotGradVisc__AZ__p
real*8     :: QdivU_BdotGradVisc__AZ__n
real*8     :: QdivU_BdotGradVisc__A3__p
real*8     :: QdivU_BdotGradVisc__Te__p
real*8     :: QdivU_BdotGradVisc__Te__n

real*8     :: curl_U0_dot_gradmu_x_B0__UR__p, curl_U0_dot_gradmu_x_B0__UR__n, curl_U0_dot_gradmu_x_B0__UZ__p, curl_U0_dot_gradmu_x_B0__UZ__n, curl_U0_dot_gradmu_x_B0__Up__p
real*8     :: curl_U0_dot_gradmu_x_B0__AR__p, curl_U0_dot_gradmu_x_B0__AR__n, curl_U0_dot_gradmu_x_B0__AZ__p, curl_U0_dot_gradmu_x_B0__AZ__n, curl_U0_dot_gradmu_x_B0__A3__p
real*8     :: curl_U0_dot_gradmu_x_B0__Te__p, curl_U0_dot_gradmu_x_B0__Te__n

real*8     :: curl_U0_dot_curl_B0__UR__p, curl_U0_dot_curl_B0__UR__n
real*8     :: curl_U0_dot_curl_B0__UZ__p, curl_U0_dot_curl_B0__UZ__n
real*8     :: curl_U0_dot_curl_B0__Up__p

real*8     :: curl_U0_dot_curl_B0__AR__p, curl_U0_dot_curl_B0__AR__n, curl_U0_dot_curl_B0__AR__nn
real*8     :: curl_U0_dot_curl_B0__AZ__p, curl_U0_dot_curl_B0__AZ__n, curl_U0_dot_curl_B0__AZ__nn
real*8     :: curl_U0_dot_curl_B0__A3__p, curl_U0_dot_curl_B0__A3__n

real*8     :: Qvisc_grad_BR_dot_grad_UR__UR__p, Qvisc_grad_BR_dot_grad_UR__UR__n
real*8     :: Qvisc_grad_BZ_dot_grad_UZ__UZ__p, Qvisc_grad_BZ_dot_grad_UZ__UZ__n
real*8     :: Qvisc_grad_Bp_dot_grad_Up__Up__p, Qvisc_grad_Bp_dot_grad_Up__Up__n
real*8     :: Qvisc_grad_BR_dot_grad_UR__AR__p
real*8     :: Qvisc_grad_BR_dot_grad_UR__AZ__n, Qvisc_grad_BR_dot_grad_UR__AZ__nn
real*8     :: Qvisc_grad_BR_dot_grad_UR__A3__p, Qvisc_grad_BR_dot_grad_UR__A3__n
real*8     :: Qvisc_grad_BZ_dot_grad_UZ__AR__p, Qvisc_grad_BZ_dot_grad_UZ__AR__nn
real*8     :: Qvisc_grad_BZ_dot_grad_UZ__AZ__p
real*8     :: Qvisc_grad_BZ_dot_grad_UZ__A3__p, Qvisc_grad_BZ_dot_grad_UZ__A3__n
real*8     :: Qvisc_grad_Bp_dot_grad_Up__AR__p, Qvisc_grad_Bp_dot_grad_Up__AR__n
real*8     :: Qvisc_grad_Bp_dot_grad_Up__AZ__p, Qvisc_grad_Bp_dot_grad_Up__AZ__n
real*8     :: Qvisc_grad_Bp_dot_grad_Up__A3__p

real*8     :: Qvisc_gradmu_dot_grad_vecV_UR__UR__p
real*8     :: Qvisc_gradmu_dot_grad_vecV_UR__UR__n
real*8     :: Qvisc_gradmu_dot_grad_vecV_UR__Up__p
real*8     :: Qvisc_gradmu_dot_grad_vecV_UR__Te__p
real*8     :: Qvisc_gradmu_dot_grad_vecV_UR__Te__n

real*8     :: Qvisc_gradmu_dot_grad_vecV_UZ__UZ__p
real*8     :: Qvisc_gradmu_dot_grad_vecV_UZ__UZ__n
real*8     :: Qvisc_gradmu_dot_grad_vecV_UZ__Te__p
real*8     :: Qvisc_gradmu_dot_grad_vecV_UZ__Te__n

real*8     :: Qvisc_gradmu_dot_grad_vecV_Up__UR__p
real*8     :: Qvisc_gradmu_dot_grad_vecV_Up__Up__p
real*8     :: Qvisc_gradmu_dot_grad_vecV_Up__Up__n
real*8     :: Qvisc_gradmu_dot_grad_vecV_Up__Te__p
real*8     :: Qvisc_gradmu_dot_grad_vecV_Up__Te__n

real*8     :: Qvisc_T
real*8     :: Qvisc_T_UR__p, Qvisc_T_UR__n
real*8     :: Qvisc_T_UZ__p, Qvisc_T_UZ__n
real*8     :: Qvisc_T_Up__p, Qvisc_T_Up__n
real*8     :: Qvisc_T_T__p,  Qvisc_T_T__n

real*8     :: Qgrad_visc_dot_grad_UR_p
real*8     :: Qgrad_visc_dot_grad_UZ_p
real*8     :: Qgrad_visc_dot_grad_Up_p
real*8     :: Qvisc_grad_BR_dot_grad_UR_p
real*8     :: Qvisc_grad_BZ_dot_grad_UZ_p
real*8     :: Qvisc_grad_Bp_dot_grad_Up_p

real*8     :: QdivU_BdotGradVisc, curl_U0_dot_gradmu_x_B0, curl_U0_dot_curl_B0

real*8     :: Qvisc_gradmu_dot_grad_vecV_UR_p
real*8     :: Qvisc_gradmu_dot_grad_vecV_UZ_p
real*8     :: Qvisc_gradmu_dot_grad_vecV_Up_p

real*8     :: Qgrad_visc_dot_gradU_T_UR__p
real*8     :: Qgrad_visc_dot_gradU_T_UZ__p
real*8     :: Qgrad_visc_dot_gradU_T_Up__p

real*8     :: Q_minus_2_over_3_divU_grad_visc_UR__p
real*8     :: Q_minus_2_over_3_divU_grad_visc_UZ__p
real*8     :: Q_minus_2_over_3_divU_grad_visc_Up__p

real*8     :: Q_1_over_3_visc_grad_divU_UR__p
real*8     :: Q_1_over_3_visc_grad_divU_UZ__p
real*8     :: Q_1_over_3_visc_grad_divU_Up__p
real*8     :: Q_1_over_3_visc_grad_divU_Up__k

real*8     :: Q_1_over_3_B_dot_visc_grad_divU__p
real*8     :: Q_1_over_3_B_dot_visc_grad_divU__k

real*8     :: Qgrad_visc_dot_gradU_T_UR__UR__p
real*8     :: Qgrad_visc_dot_gradU_T_UR__UZ__p
real*8     :: Qgrad_visc_dot_gradU_T_UR__Up__p
real*8     :: Qgrad_visc_dot_gradU_T_UR__Te__p
real*8     :: Qgrad_visc_dot_gradU_T_UR__Te__n

real*8     :: Qgrad_visc_dot_gradU_T_UZ__UR__p
real*8     :: Qgrad_visc_dot_gradU_T_UZ__UZ__p
real*8     :: Qgrad_visc_dot_gradU_T_UZ__Up__p
real*8     :: Qgrad_visc_dot_gradU_T_UZ__Te__p
real*8     :: Qgrad_visc_dot_gradU_T_UZ__Te__n

real*8     :: Qgrad_visc_dot_gradU_T_Up__UR__p
real*8     :: Qgrad_visc_dot_gradU_T_Up__UR__n
real*8     :: Qgrad_visc_dot_gradU_T_Up__UZ__n
real*8     :: Qgrad_visc_dot_gradU_T_Up__Up__p
real*8     :: Qgrad_visc_dot_gradU_T_Up__Up__n
real*8     :: Qgrad_visc_dot_gradU_T_Up__Te__p
real*8     :: Qgrad_visc_dot_gradU_T_Up__Te__n

real*8     :: Q_minus_2_over_3_divU_grad_visc_UR__UR__p
real*8     :: Q_minus_2_over_3_divU_grad_visc_UR__UZ__p
real*8     :: Q_minus_2_over_3_divU_grad_visc_UR__Up__n
real*8     :: Q_minus_2_over_3_divU_grad_visc_UR__Te__p

real*8     :: Q_minus_2_over_3_divU_grad_visc_UZ__UR__p
real*8     :: Q_minus_2_over_3_divU_grad_visc_UZ__UZ__p
real*8     :: Q_minus_2_over_3_divU_grad_visc_UZ__Up__n
real*8     :: Q_minus_2_over_3_divU_grad_visc_UZ__Te__p

real*8     :: Q_minus_2_over_3_divU_grad_visc_Up__UR__p
real*8     :: Q_minus_2_over_3_divU_grad_visc_Up__UZ__p
real*8     :: Q_minus_2_over_3_divU_grad_visc_Up__Up__n
real*8     :: Q_minus_2_over_3_divU_grad_visc_Up__Te__p
real*8     :: Q_minus_2_over_3_divU_grad_visc_Up__Te__n

real*8     :: Q_1_over_3_visc_grad_divU_UR__UR__p
real*8     :: Q_1_over_3_visc_grad_divU_UR__UZ__p
real*8     :: Q_1_over_3_visc_grad_divU_UR__Up__n
real*8     :: Q_1_over_3_visc_grad_divU_UR__Te__p

real*8     :: Q_1_over_3_visc_grad_divU_UZ__UR__p
real*8     :: Q_1_over_3_visc_grad_divU_UZ__UZ__p
real*8     :: Q_1_over_3_visc_grad_divU_UZ__Up__n
real*8     :: Q_1_over_3_visc_grad_divU_UZ__Te__p

real*8     :: Q_1_over_3_visc_grad_divU_Up__UR__p 
real*8     :: Q_1_over_3_visc_grad_divU_Up__UZ__p 
real*8     :: Q_1_over_3_visc_grad_divU_Up__Up__n 
real*8     :: Q_1_over_3_visc_grad_divU_Up__Te__p 
real*8     :: Q_1_over_3_visc_grad_divU_Up__Te__n 
real*8     :: Q_1_over_3_visc_grad_divU_Up__UR__k 
real*8     :: Q_1_over_3_visc_grad_divU_Up__UZ__k 
real*8     :: Q_1_over_3_visc_grad_divU_Up__Up__kn
real*8     :: Q_1_over_3_visc_grad_divU_Up__Te__k 

real*8     :: Q_1_over_3_B_dot_visc_grad_divU__AR__p
real*8     :: Q_1_over_3_B_dot_visc_grad_divU__AR__n
real*8     :: Q_1_over_3_B_dot_visc_grad_divU__AZ__p
real*8     :: Q_1_over_3_B_dot_visc_grad_divU__AZ__n
real*8     :: Q_1_over_3_B_dot_visc_grad_divU__A3__p
real*8     :: Q_1_over_3_B_dot_visc_grad_divU__UR__p
real*8     :: Q_1_over_3_B_dot_visc_grad_divU__UZ__p
real*8     :: Q_1_over_3_B_dot_visc_grad_divU__Up__n
real*8     :: Q_1_over_3_B_dot_visc_grad_divU__Te__p
real*8     :: Q_1_over_3_B_dot_visc_grad_divU__Te__n

real*8     :: Q_1_over_3_B_dot_visc_grad_divU__AR__k 
real*8     :: Q_1_over_3_B_dot_visc_grad_divU__AZ__k 
real*8     :: Q_1_over_3_B_dot_visc_grad_divU__A3__k 
real*8     :: Q_1_over_3_B_dot_visc_grad_divU__UR__k 
real*8     :: Q_1_over_3_B_dot_visc_grad_divU__UZ__k 
real*8     :: Q_1_over_3_B_dot_visc_grad_divU__Up__kn
real*8     :: Q_1_over_3_B_dot_visc_grad_divU__Te__k 

! --- fourth order stabilization (numerical diffusion terms)
real*8     :: lap_Vstar, lap_bf
real*8     :: lap_AR, lap_AZ, lap_A3
real*8     :: lap_UR, lap_UZ, lap_Up
real*8     :: lap_rho,lap_Ti, lap_Te

! --- Ion-electron energy transfer
real*8     :: t_norm, nu_e_bg, lambda_e_bg, dTi_e, dTe_i
real*8     :: dnu_e_bg_dTi, dnu_e_bg_dTe
real*8     :: dnu_e_bg_drho, dnu_e_bg_drhon
real*8     :: ddTi_e_dTi, ddTi_e_dTe, ddTi_e_drho, ddTi_e_drhon
real*8     :: ddTe_i_dTi, ddTe_i_dTe, ddTe_i_drho, ddTe_i_drhon
real*8     :: Te_corr_eV, dTe_corr_eV_dT                      ! Electron temperature in eV
real*8     :: ne_SI                                          ! Electron density in SI unit
real*8     :: drho0_corr_dn, dTi0_corr_dT, dTe0_corr_dT

! --- Matrix
real*8, dimension(n_var      )   :: rhs_p_ij, rhs_k_ij, Pvec_prev, Qvec_p, Qvec_k
real*8, dimension(n_var,n_var)   :: amat, Pjac, Qjac_p, Qjac_k, Qjac_n, Qjac_kn

! --- Ohmic heating, for details please see:
! https://www.jorek.eu/wiki/doku.php?id=ohmic_heating
real*8, dimension(DIM1, DIM2, DIM2), intent(inout) :: ELM_pnn
real*8, dimension(n_plane,n_var,n_gauss,n_gauss) :: eq_pp, eq_sp, eq_tp
real*8, dimension(n_var,n_var)   :: Qjac_pnn
real*8     :: AR0_Rp, AR0_Zp, AR0_sp, AR0_tp
real*8     :: AZ0_Rp, AZ0_Zp, AZ0_sp, AZ0_tp
real*8     :: A30_Rp, A30_Zp, A30_sp, A30_tp
real*8     :: bf_RZ, bf_pp, bf_sp, bf_tp, bf_Rp, bf_Zp
real*8     :: AR_RR, AR_RZ, AR_ZZ, AR_pp, AR_Rp, AR_Zp
real*8     :: AZ_RR, AZ_RZ, AZ_ZZ, AZ_pp, AZ_Rp, AZ_Zp
real*8     :: A3_RR, A3_RZ, A3_ZZ, A3_pp, A3_Rp, A3_Zp
real*8     :: psieq_R,  psieq_Z
real*8     :: BR0_R, BR0_Z, BR0_p, BZ0_R, BZ0_Z, BZ0_p, Bp0_R, Bp0_Z, Bp0_p
real*8     :: JR0 , JR0_AR__p, JR0_AR__n, JR0_AR__nn, JR0_AZ__p, JR0_AZ__n, JR0_AZ__nn, JR0_A3__p, JR0_A3__n, JR0_A3__nn
real*8     :: JZ0 , JZ0_AR__p, JZ0_AR__n, JZ0_AR__nn, JZ0_AZ__p, JZ0_AZ__n, JZ0_AZ__nn, JZ0_A3__p, JZ0_A3__n, JZ0_A3__nn
real*8     :: Jp0 , Jp0_AR__p, Jp0_AR__n, Jp0_AR__nn, Jp0_AZ__p, Jp0_AZ__n, Jp0_AZ__nn, Jp0_A3__p, Jp0_A3__n, Jp0_A3__nn
real*8     :: JJ2,  JJ2_AR__p, JJ2_AR__n, JJ2_AR__nn, JJ2_AZ__p, JJ2_AZ__n, JJ2_AZ__nn, JJ2_A3__p, JJ2_A3__n, JJ2_A3__nn

!  --- For shock capturing stabilization
real*8     :: midp_edge1(1:2), midp_edge2(1:2), midp_edge3(1:2), midp_edge4(1:2)
real*8     :: len1, len2, h_e
real*8     :: Ptot_R,  Ptot_Z,  Ptot_p, Ptot_corr
real*8     :: f_p, d_p, tau_sc, R_rho, R_Ti, R_Te
real*8     :: s_p, src_rho, src_pi, src_pe
real*8     :: rhoi_eff, rhoe_eff

!  --- For VMS stabilization
integer    :: jj
real*8     :: tscale
real*8     :: speed(n_plane,n_gauss,n_gauss)
real*8     :: vms_AR__p(n_var), vms_AR__k(n_var)
real*8     :: vms_AZ__p(n_var), vms_AZ__k(n_var)
real*8     :: vms_A3__p(n_var), vms_A3__k(n_var)
real*8     :: vms_UR__p(n_var), vms_UR__k(n_var)
real*8     :: vms_UZ__p(n_var), vms_UZ__k(n_var)
real*8     :: vms_Up__p(n_var), vms_Up__k(n_var)
real*8     :: vms_rho__p(n_var), vms_rho__k(n_var)
real*8     :: vms_Ti__p(n_var), vms_Ti__k(n_var)
real*8     :: vms_Te__p(n_var), vms_Te__k(n_var)
real*8     :: Pvec_prev_k(n_var), Pjac_k(n_var,n_var)
real*8     :: res(n_var), res_jac__p(n_var, n_var), res_jac__n(n_var, n_var), res_jac__nn(n_var,n_var)
real*8     :: vsR, vsR_UR__p, vsR_UR__n, vsR_UR__nn, vsR_UZ__p, vsR_UZ__n, vsR_UZ__nn, vsR_Up__p, vsR_Up__n, vsR_Up__nn
real*8     :: vsZ, vsZ_UR__p, vsZ_UR__n, vsZ_UR__nn, vsZ_UZ__p, vsZ_UZ__n, vsZ_UZ__nn, vsZ_Up__p, vsZ_Up__n, vsZ_Up__nn
real*8     :: vsp, vsp_UR__p, vsp_UR__n, vsp_UR__nn, vsp_UZ__p, vsp_UZ__n, vsp_UZ__nn, vsp_Up__p, vsp_Up__n, vsp_Up__nn

real*8, dimension(DIM1, DIM2, DIM2), intent(inout) :: ELM_knn
#ifdef ENABLE_VISCO_PAR_AND_PARPAR
!  --- For visco_par and visco_par_par
real*8, dimension(n_var,n_var)                     :: Qjac_knn
real*8     :: qparR, qparZ, qparPhi, qparpar, qparparR, qparparZ, qparparPhi
real*8     :: grad_BRVstar_R, grad_BRVstar_Z, grad_BRVstar_phi, grad_BRVstar_phi_k
real*8     :: grad_BpVstar_R, grad_BpVstar_Z, grad_BpVstar_phi, grad_BpVstar_phi_k
real*8     :: grad_BZVstar_R, grad_BZVstar_Z, grad_BZVstar_phi, grad_BZVstar_phi_k
real*8     :: BB2_R, BB2_Z, BB2_phi
real*8     :: B_dot_grad_BB2, B_dot_U0, B_dot_grad_BdotU0
real*8     :: grad_BdotU0_R, grad_BdotU0_Z, grad_BdotU0_phi

real*8     :: grad_BR_AR__R, grad_BZ_AR__R__n, grad_Bp_AR__R
real*8     :: grad_BR_AZ__R__n, grad_BZ_AZ__R, grad_Bp_AZ__R
real*8     :: grad_BR_A3__R, grad_BZ_A3__R, grad_Bp_A3__R

real*8     :: grad_BR_AR__Z, grad_BZ_AR__Z__n, grad_Bp_AR__Z
real*8     :: grad_BR_AZ__Z__n, grad_BZ_AZ__Z, grad_Bp_AZ__Z
real*8     :: grad_BR_A3__Z, grad_BZ_A3__Z, grad_Bp_A3__Z

real*8     :: grad_BR_AR__phi, grad_BZ_AR__phi__nn, grad_Bp_AR__phi__n
real*8     :: grad_BR_AZ__phi__nn, grad_BZ_AZ__phi, grad_Bp_AZ__phi__n
real*8     :: grad_BR_A3__phi__n, grad_BZ_A3__phi__n, grad_Bp_A3__phi

real*8     :: B_dot_grad_BR, B_dot_grad_BZ, B_dot_grad_Bp
real*8     :: B_AR_dot_U__p, B_AR_dot_U__n, B_AZ_dot_U__p, B_AZ_dot_U__n, B_A3_dot_U__p, B_A3_dot_U__n
real*8     :: grad_B_AR_dot_U__R__p, grad_B_AR_dot_U__R__n, grad_B_AZ_dot_U__R__p, grad_B_AZ_dot_U__R__n, grad_B_A3_dot_U__R__p, grad_B_A3_dot_U__R__n
real*8     :: grad_B_AR_dot_U__Z__p, grad_B_AR_dot_U__Z__n, grad_B_AZ_dot_U__Z__p, grad_B_AZ_dot_U__Z__n, grad_B_A3_dot_U__Z__p, grad_B_A3_dot_U__Z__n
real*8     :: grad_B_AR_dot_U__phi__p, grad_B_AR_dot_U__phi__n, grad_B_AR_dot_U__phi__nn, grad_B_AZ_dot_U__phi__p, grad_B_AZ_dot_U__phi__n, grad_B_AZ_dot_U__phi__nn, grad_B_A3_dot_U__phi__p, grad_B_A3_dot_U__phi__n, grad_B_A3_dot_U__phi__nn
real*8     :: grad_BB2_AR__R__p, grad_BB2_AR__R__n, grad_BB2_AR__R__nn, grad_BB2_AZ__R__p, grad_BB2_AZ__R__n, grad_BB2_AZ__R__nn, grad_BB2_A3__R__p, grad_BB2_A3__R__n, grad_BB2_A3__R__nn
real*8     :: grad_BB2_AR__Z__p, grad_BB2_AR__Z__n, grad_BB2_AR__Z__nn, grad_BB2_AZ__Z__p, grad_BB2_AZ__Z__n, grad_BB2_AZ__Z__nn, grad_BB2_A3__Z__p, grad_BB2_A3__Z__n, grad_BB2_A3__Z__nn
real*8     :: grad_BB2_AR__phi__p, grad_BB2_AR__phi__n, grad_BB2_AR__phi__nn, grad_BB2_AZ__phi__p, grad_BB2_AZ__phi__n, grad_BB2_AZ__phi__nn, grad_BB2_A3__phi__p, grad_BB2_A3__phi__n, grad_BB2_A3__phi__nn
real*8     :: B_dot_grad_BRVstar__p, B_dot_grad_BRVstar__k, B_dot_grad_BZVstar__p, B_dot_grad_BZVstar__k, B_dot_grad_BpVstar__p, B_dot_grad_BpVstar__k

real*8     :: B_dot_grad_BRVstar__AR__p, B_dot_grad_BRVstar__AR__n, B_dot_grad_BRVstar__AR__k, B_dot_grad_BRVstar__AR__kn, B_dot_grad_BRVstar__AR__nn
real*8     :: B_dot_grad_BZVstar__AR__p, B_dot_grad_BZVstar__AR__n, B_dot_grad_BZVstar__AR__k, B_dot_grad_BZVstar__AR__kn, B_dot_grad_BZVstar__AR__nn
real*8     :: B_dot_grad_BpVstar__AR__p, B_dot_grad_BpVstar__AR__n, B_dot_grad_BpVstar__AR__k, B_dot_grad_BpVstar__AR__kn, B_dot_grad_BpVstar__AR__nn

real*8     :: B_dot_grad_BRVstar__AZ__p, B_dot_grad_BRVstar__AZ__n, B_dot_grad_BRVstar__AZ__k, B_dot_grad_BRVstar__AZ__kn, B_dot_grad_BRVstar__AZ__nn
real*8     :: B_dot_grad_BZVstar__AZ__p, B_dot_grad_BZVstar__AZ__n, B_dot_grad_BZVstar__AZ__k, B_dot_grad_BZVstar__AZ__kn, B_dot_grad_BZVstar__AZ__nn
real*8     :: B_dot_grad_BpVstar__AZ__p, B_dot_grad_BpVstar__AZ__n, B_dot_grad_BpVstar__AZ__k, B_dot_grad_BpVstar__AZ__kn, B_dot_grad_BpVstar__AZ__nn

real*8     :: B_dot_grad_BRVstar__A3__p, B_dot_grad_BRVstar__A3__n, B_dot_grad_BRVstar__A3__k, B_dot_grad_BRVstar__A3__kn, B_dot_grad_BRVstar__A3__nn
real*8     :: B_dot_grad_BZVstar__A3__p, B_dot_grad_BZVstar__A3__n, B_dot_grad_BZVstar__A3__k, B_dot_grad_BZVstar__A3__kn, B_dot_grad_BZVstar__A3__nn
real*8     :: B_dot_grad_BpVstar__A3__p, B_dot_grad_BpVstar__A3__n, B_dot_grad_BpVstar__A3__k, B_dot_grad_BpVstar__A3__kn, B_dot_grad_BpVstar__A3__nn

real*8     :: B_dot_grad_BB2__AR__p, B_dot_grad_BB2__AR__n, B_dot_grad_BB2__AR__nn
real*8     :: B_dot_grad_BB2__AZ__p, B_dot_grad_BB2__AZ__n, B_dot_grad_BB2__AZ__nn
real*8     :: B_dot_grad_BB2__A3__p, B_dot_grad_BB2__A3__n, B_dot_grad_BB2__A3__nn
real*8     :: B_AR_dot_grad_BdotU0__p, B_AR_dot_grad_BdotU0__n, B_AZ_dot_grad_BdotU0__p, B_AZ_dot_grad_BdotU0__n, B_A3_dot_grad_BdotU0__p, B_A3_dot_grad_BdotU0__n

real*8     :: B_dot_grad__BdotU0_AR__p, B_dot_grad__BdotU0_AR__n, B_dot_grad__BdotU0_AR__nn
real*8     :: B_dot_grad__BdotU0_AZ__p, B_dot_grad__BdotU0_AZ__n, B_dot_grad__BdotU0_AZ__nn
real*8     :: B_dot_grad__BdotU0_A3__p, B_dot_grad__BdotU0_A3__n, B_dot_grad__BdotU0_A3__nn
real*8     :: B_dot_grad__UR__p, B_dot_grad__UR__n, B_dot_grad__UZ__p, B_dot_grad__UZ__n, B_dot_grad__Up__p, B_dot_grad__Up__n

real*8     :: Qvec_UR_visco_par__p, Qvec_UR_visco_par__k, Qvec_UR_visco_parpar__p, Qvec_UR_visco_parpar__k
real*8     :: Qvec_UZ_visco_par__p, Qvec_UZ_visco_par__k, Qvec_UZ_visco_parpar__p, Qvec_UZ_visco_parpar__k
real*8     :: Qvec_Up_visco_par__p, Qvec_Up_visco_par__k, Qvec_Up_visco_parpar__p, Qvec_Up_visco_parpar__k

real*8     :: Qjac_UR_visco_par_AR_p, Qjac_UR_visco_par_AR_n, Qjac_UR_visco_par_AR_k, Qjac_UR_visco_par_AR_kn, Qjac_UR_visco_par_AR_pnn, Qjac_UR_visco_par_AR_knn
real*8     :: Qjac_UR_visco_parpar_AR_p, Qjac_UR_visco_parpar_AR_n, Qjac_UR_visco_parpar_AR_k, Qjac_UR_visco_parpar_AR_kn, Qjac_UR_visco_parpar_AR_pnn, Qjac_UR_visco_parpar_AR_knn
real*8     :: Qjac_UR_visco_par_AZ_p, Qjac_UR_visco_par_AZ_n, Qjac_UR_visco_par_AZ_k, Qjac_UR_visco_par_AZ_kn, Qjac_UR_visco_par_AZ_pnn, Qjac_UR_visco_par_AZ_knn
real*8     :: Qjac_UR_visco_parpar_AZ_p, Qjac_UR_visco_parpar_AZ_n, Qjac_UR_visco_parpar_AZ_k, Qjac_UR_visco_parpar_AZ_kn, Qjac_UR_visco_parpar_AZ_pnn, Qjac_UR_visco_parpar_AZ_knn
real*8     :: Qjac_UR_visco_par_A3_p, Qjac_UR_visco_par_A3_n, Qjac_UR_visco_par_A3_k, Qjac_UR_visco_par_A3_kn, Qjac_UR_visco_par_A3_pnn, Qjac_UR_visco_par_A3_knn
real*8     :: Qjac_UR_visco_parpar_A3_p, Qjac_UR_visco_parpar_A3_n, Qjac_UR_visco_parpar_A3_k, Qjac_UR_visco_parpar_A3_kn, Qjac_UR_visco_parpar_A3_pnn, Qjac_UR_visco_parpar_A3_knn

real*8     :: Qjac_UR_visco_par_UR_p, Qjac_UR_visco_par_UR_n, Qjac_UR_visco_par_UR_k, Qjac_UR_visco_par_UR_kn, Qjac_UR_visco_par_UR_pnn, Qjac_UR_visco_par_UR_knn
real*8     :: Qjac_UR_visco_parpar_UR_p, Qjac_UR_visco_parpar_UR_n, Qjac_UR_visco_parpar_UR_k, Qjac_UR_visco_parpar_UR_kn, Qjac_UR_visco_parpar_UR_pnn, Qjac_UR_visco_parpar_UR_knn
real*8     :: Qjac_UR_visco_par_UZ_p, Qjac_UR_visco_par_UZ_n, Qjac_UR_visco_par_UZ_k, Qjac_UR_visco_par_UZ_kn, Qjac_UR_visco_par_UZ_pnn, Qjac_UR_visco_par_UZ_knn
real*8     :: Qjac_UR_visco_parpar_UZ_p, Qjac_UR_visco_parpar_UZ_n, Qjac_UR_visco_parpar_UZ_k, Qjac_UR_visco_parpar_UZ_kn, Qjac_UR_visco_parpar_UZ_pnn, Qjac_UR_visco_parpar_UZ_knn
real*8     :: Qjac_UR_visco_par_Up_p, Qjac_UR_visco_par_Up_n, Qjac_UR_visco_par_Up_k, Qjac_UR_visco_par_Up_kn, Qjac_UR_visco_par_Up_pnn, Qjac_UR_visco_par_Up_knn
real*8     :: Qjac_UR_visco_parpar_Up_p, Qjac_UR_visco_parpar_Up_n, Qjac_UR_visco_parpar_Up_k, Qjac_UR_visco_parpar_Up_kn, Qjac_UR_visco_parpar_Up_pnn, Qjac_UR_visco_parpar_Up_knn

real*8     :: Qjac_UZ_visco_par_AR_p, Qjac_UZ_visco_par_AR_n, Qjac_UZ_visco_par_AR_k, Qjac_UZ_visco_par_AR_kn, Qjac_UZ_visco_par_AR_pnn, Qjac_UZ_visco_par_AR_knn
real*8     :: Qjac_UZ_visco_parpar_AR_p, Qjac_UZ_visco_parpar_AR_n, Qjac_UZ_visco_parpar_AR_k, Qjac_UZ_visco_parpar_AR_kn, Qjac_UZ_visco_parpar_AR_pnn, Qjac_UZ_visco_parpar_AR_knn
real*8     :: Qjac_UZ_visco_par_AZ_p, Qjac_UZ_visco_par_AZ_n, Qjac_UZ_visco_par_AZ_k, Qjac_UZ_visco_par_AZ_kn, Qjac_UZ_visco_par_AZ_pnn, Qjac_UZ_visco_par_AZ_knn
real*8     :: Qjac_UZ_visco_parpar_AZ_p, Qjac_UZ_visco_parpar_AZ_n, Qjac_UZ_visco_parpar_AZ_k, Qjac_UZ_visco_parpar_AZ_kn, Qjac_UZ_visco_parpar_AZ_pnn, Qjac_UZ_visco_parpar_AZ_knn
real*8     :: Qjac_UZ_visco_par_A3_p, Qjac_UZ_visco_par_A3_n, Qjac_UZ_visco_par_A3_k, Qjac_UZ_visco_par_A3_kn, Qjac_UZ_visco_par_A3_pnn, Qjac_UZ_visco_par_A3_knn
real*8     :: Qjac_UZ_visco_parpar_A3_p, Qjac_UZ_visco_parpar_A3_n, Qjac_UZ_visco_parpar_A3_k, Qjac_UZ_visco_parpar_A3_kn, Qjac_UZ_visco_parpar_A3_pnn, Qjac_UZ_visco_parpar_A3_knn

real*8     :: Qjac_UZ_visco_par_UR_p, Qjac_UZ_visco_par_UR_n, Qjac_UZ_visco_par_UR_k, Qjac_UZ_visco_par_UR_kn, Qjac_UZ_visco_par_UR_pnn, Qjac_UZ_visco_par_UR_knn
real*8     :: Qjac_UZ_visco_parpar_UR_p, Qjac_UZ_visco_parpar_UR_n, Qjac_UZ_visco_parpar_UR_k, Qjac_UZ_visco_parpar_UR_kn, Qjac_UZ_visco_parpar_UR_pnn, Qjac_UZ_visco_parpar_UR_knn
real*8     :: Qjac_UZ_visco_par_UZ_p, Qjac_UZ_visco_par_UZ_n, Qjac_UZ_visco_par_UZ_k, Qjac_UZ_visco_par_UZ_kn, Qjac_UZ_visco_par_UZ_pnn, Qjac_UZ_visco_par_UZ_knn
real*8     :: Qjac_UZ_visco_parpar_UZ_p, Qjac_UZ_visco_parpar_UZ_n, Qjac_UZ_visco_parpar_UZ_k, Qjac_UZ_visco_parpar_UZ_kn, Qjac_UZ_visco_parpar_UZ_pnn, Qjac_UZ_visco_parpar_UZ_knn
real*8     :: Qjac_UZ_visco_par_Up_p, Qjac_UZ_visco_par_Up_n, Qjac_UZ_visco_par_Up_k, Qjac_UZ_visco_par_Up_kn, Qjac_UZ_visco_par_Up_pnn, Qjac_UZ_visco_par_Up_knn
real*8     :: Qjac_UZ_visco_parpar_Up_p, Qjac_UZ_visco_parpar_Up_n, Qjac_UZ_visco_parpar_Up_k, Qjac_UZ_visco_parpar_Up_kn, Qjac_UZ_visco_parpar_Up_pnn, Qjac_UZ_visco_parpar_Up_knn

real*8     :: Qjac_Up_visco_par_AR_p, Qjac_Up_visco_par_AR_n, Qjac_Up_visco_par_AR_k, Qjac_Up_visco_par_AR_kn, Qjac_Up_visco_par_AR_pnn, Qjac_Up_visco_par_AR_knn
real*8     :: Qjac_Up_visco_parpar_AR_p, Qjac_Up_visco_parpar_AR_n, Qjac_Up_visco_parpar_AR_k, Qjac_Up_visco_parpar_AR_kn, Qjac_Up_visco_parpar_AR_pnn, Qjac_Up_visco_parpar_AR_knn
real*8     :: Qjac_Up_visco_par_AZ_p, Qjac_Up_visco_par_AZ_n, Qjac_Up_visco_par_AZ_k, Qjac_Up_visco_par_AZ_kn, Qjac_Up_visco_par_AZ_pnn, Qjac_Up_visco_par_AZ_knn
real*8     :: Qjac_Up_visco_parpar_AZ_p, Qjac_Up_visco_parpar_AZ_n, Qjac_Up_visco_parpar_AZ_k, Qjac_Up_visco_parpar_AZ_kn, Qjac_Up_visco_parpar_AZ_pnn, Qjac_Up_visco_parpar_AZ_knn
real*8     :: Qjac_Up_visco_par_A3_p, Qjac_Up_visco_par_A3_n, Qjac_Up_visco_par_A3_k, Qjac_Up_visco_par_A3_kn, Qjac_Up_visco_par_A3_pnn, Qjac_Up_visco_par_A3_knn
real*8     :: Qjac_Up_visco_parpar_A3_p, Qjac_Up_visco_parpar_A3_n, Qjac_Up_visco_parpar_A3_k, Qjac_Up_visco_parpar_A3_kn, Qjac_Up_visco_parpar_A3_pnn, Qjac_Up_visco_parpar_A3_knn

real*8     :: Qjac_Up_visco_par_UR_p, Qjac_Up_visco_par_UR_n, Qjac_Up_visco_par_UR_k, Qjac_Up_visco_par_UR_kn, Qjac_Up_visco_par_UR_pnn, Qjac_Up_visco_par_UR_knn
real*8     :: Qjac_Up_visco_parpar_UR_p, Qjac_Up_visco_parpar_UR_n, Qjac_Up_visco_parpar_UR_k, Qjac_Up_visco_parpar_UR_kn, Qjac_Up_visco_parpar_UR_pnn, Qjac_Up_visco_parpar_UR_knn
real*8     :: Qjac_Up_visco_par_UZ_p, Qjac_Up_visco_par_UZ_n, Qjac_Up_visco_par_UZ_k, Qjac_Up_visco_par_UZ_kn, Qjac_Up_visco_par_UZ_pnn, Qjac_Up_visco_par_UZ_knn
real*8     :: Qjac_Up_visco_parpar_UZ_p, Qjac_Up_visco_parpar_UZ_n, Qjac_Up_visco_parpar_UZ_k, Qjac_Up_visco_parpar_UZ_kn, Qjac_Up_visco_parpar_UZ_pnn, Qjac_Up_visco_parpar_UZ_knn
real*8     :: Qjac_Up_visco_par_Up_p, Qjac_Up_visco_par_Up_n, Qjac_Up_visco_par_Up_k, Qjac_Up_visco_par_Up_kn, Qjac_Up_visco_par_Up_pnn, Qjac_Up_visco_par_Up_knn
real*8     :: Qjac_Up_visco_parpar_Up_p, Qjac_Up_visco_parpar_Up_n, Qjac_Up_visco_parpar_Up_k, Qjac_Up_visco_parpar_Up_kn, Qjac_Up_visco_parpar_Up_pnn, Qjac_Up_visco_parpar_Up_knn


real*8     :: B_dot_grad_BR__AR__p, B_dot_grad_BR__AR__n
real*8     :: B_dot_grad_BZ__AR__p, B_dot_grad_BZ__AR__n, B_dot_grad_BZ__AR__nn
real*8     :: B_dot_grad_Bp__AR__p, B_dot_grad_Bp__AR__n
real*8     :: B_dot_grad_BR__AZ__p, B_dot_grad_BR__AZ__n, B_dot_grad_BR__AZ__nn
real*8     :: B_dot_grad_BZ__AZ__p, B_dot_grad_BZ__AZ__n
real*8     :: B_dot_grad_Bp__AZ__p, B_dot_grad_Bp__AZ__n
real*8     :: B_dot_grad_BR__A3__p, B_dot_grad_BR__A3__n
real*8     :: B_dot_grad_BZ__A3__p, B_dot_grad_BZ__A3__n
real*8     :: B_dot_grad_Bp__A3__p

real*8     :: Qvec_UR_visco_par_for_Uphi__p, Qvec_UR_visco_parpar_for_Uphi__p
real*8     :: Qvec_UZ_visco_par_for_Uphi__p, Qvec_UZ_visco_parpar_for_Uphi__p
real*8     :: Qvec_Up_visco_par_for_Uphi__p, Qvec_Up_visco_parpar_for_Uphi__p

! Expanded into _p, _n, _nn, and apply on Uphi equation
real*8     :: Qjac_UR_visco_par_AR_for_Uphi_p, Qjac_UR_visco_par_AR_for_Uphi_n, Qjac_UR_visco_par_AR_for_Uphi_pnn
real*8     :: Qjac_UR_visco_parpar_AR_for_Uphi_p, Qjac_UR_visco_parpar_AR_for_Uphi_n, Qjac_UR_visco_parpar_AR_for_Uphi_pnn
real*8     :: Qjac_UR_visco_par_AZ_for_Uphi_p, Qjac_UR_visco_par_AZ_for_Uphi_n, Qjac_UR_visco_par_AZ_for_Uphi_pnn
real*8     :: Qjac_UR_visco_parpar_AZ_for_Uphi_p, Qjac_UR_visco_parpar_AZ_for_Uphi_n, Qjac_UR_visco_parpar_AZ_for_Uphi_pnn
real*8     :: Qjac_UR_visco_par_A3_for_Uphi_p, Qjac_UR_visco_par_A3_for_Uphi_n, Qjac_UR_visco_par_A3_for_Uphi_pnn
real*8     :: Qjac_UR_visco_parpar_A3_for_Uphi_p, Qjac_UR_visco_parpar_A3_for_Uphi_n, Qjac_UR_visco_parpar_A3_for_Uphi_pnn
real*8     :: Qjac_UR_visco_par_UR_for_Uphi_p, Qjac_UR_visco_par_UR_for_Uphi_n, Qjac_UR_visco_par_UR_for_Uphi_pnn
real*8     :: Qjac_UR_visco_parpar_UR_for_Uphi_p, Qjac_UR_visco_parpar_UR_for_Uphi_n, Qjac_UR_visco_parpar_UR_for_Uphi_pnn
real*8     :: Qjac_UR_visco_par_UZ_for_Uphi_p, Qjac_UR_visco_par_UZ_for_Uphi_n, Qjac_UR_visco_par_UZ_for_Uphi_pnn
real*8     :: Qjac_UR_visco_parpar_UZ_for_Uphi_p, Qjac_UR_visco_parpar_UZ_for_Uphi_n, Qjac_UR_visco_parpar_UZ_for_Uphi_pnn
real*8     :: Qjac_UR_visco_par_Up_for_Uphi_p, Qjac_UR_visco_par_Up_for_Uphi_n, Qjac_UR_visco_par_Up_for_Uphi_pnn
real*8     :: Qjac_UR_visco_parpar_Up_for_Uphi_p, Qjac_UR_visco_parpar_Up_for_Uphi_n, Qjac_UR_visco_parpar_Up_for_Uphi_pnn
real*8     :: Qjac_UZ_visco_par_AR_for_Uphi_p, Qjac_UZ_visco_par_AR_for_Uphi_n, Qjac_UZ_visco_par_AR_for_Uphi_pnn
real*8     :: Qjac_UZ_visco_parpar_AR_for_Uphi_p, Qjac_UZ_visco_parpar_AR_for_Uphi_n, Qjac_UZ_visco_parpar_AR_for_Uphi_pnn
real*8     :: Qjac_UZ_visco_par_AZ_for_Uphi_p, Qjac_UZ_visco_par_AZ_for_Uphi_n, Qjac_UZ_visco_par_AZ_for_Uphi_pnn
real*8     :: Qjac_UZ_visco_parpar_AZ_for_Uphi_p, Qjac_UZ_visco_parpar_AZ_for_Uphi_n, Qjac_UZ_visco_parpar_AZ_for_Uphi_pnn
real*8     :: Qjac_UZ_visco_par_A3_for_Uphi_p, Qjac_UZ_visco_par_A3_for_Uphi_n, Qjac_UZ_visco_par_A3_for_Uphi_pnn
real*8     :: Qjac_UZ_visco_parpar_A3_for_Uphi_p, Qjac_UZ_visco_parpar_A3_for_Uphi_n, Qjac_UZ_visco_parpar_A3_for_Uphi_pnn
real*8     :: Qjac_UZ_visco_par_UR_for_Uphi_p, Qjac_UZ_visco_par_UR_for_Uphi_n, Qjac_UZ_visco_par_UR_for_Uphi_pnn
real*8     :: Qjac_UZ_visco_parpar_UR_for_Uphi_p, Qjac_UZ_visco_parpar_UR_for_Uphi_n, Qjac_UZ_visco_parpar_UR_for_Uphi_pnn
real*8     :: Qjac_UZ_visco_par_UZ_for_Uphi_p, Qjac_UZ_visco_par_UZ_for_Uphi_n, Qjac_UZ_visco_par_UZ_for_Uphi_pnn
real*8     :: Qjac_UZ_visco_parpar_UZ_for_Uphi_p, Qjac_UZ_visco_parpar_UZ_for_Uphi_n, Qjac_UZ_visco_parpar_UZ_for_Uphi_pnn
real*8     :: Qjac_UZ_visco_par_Up_for_Uphi_p, Qjac_UZ_visco_par_Up_for_Uphi_n, Qjac_UZ_visco_par_Up_for_Uphi_pnn
real*8     :: Qjac_UZ_visco_parpar_Up_for_Uphi_p, Qjac_UZ_visco_parpar_Up_for_Uphi_n, Qjac_UZ_visco_parpar_Up_for_Uphi_pnn
real*8     :: Qjac_Up_visco_par_AR_for_Uphi_p, Qjac_Up_visco_par_AR_for_Uphi_n, Qjac_Up_visco_par_AR_for_Uphi_pnn
real*8     :: Qjac_Up_visco_parpar_AR_for_Uphi_p, Qjac_Up_visco_parpar_AR_for_Uphi_n, Qjac_Up_visco_parpar_AR_for_Uphi_pnn
real*8     :: Qjac_Up_visco_par_AZ_for_Uphi_p, Qjac_Up_visco_par_AZ_for_Uphi_n, Qjac_Up_visco_par_AZ_for_Uphi_pnn
real*8     :: Qjac_Up_visco_parpar_AZ_for_Uphi_p, Qjac_Up_visco_parpar_AZ_for_Uphi_n, Qjac_Up_visco_parpar_AZ_for_Uphi_pnn
real*8     :: Qjac_Up_visco_par_A3_for_Uphi_p, Qjac_Up_visco_par_A3_for_Uphi_n, Qjac_Up_visco_par_A3_for_Uphi_pnn
real*8     :: Qjac_Up_visco_parpar_A3_for_Uphi_p, Qjac_Up_visco_parpar_A3_for_Uphi_n, Qjac_Up_visco_parpar_A3_for_Uphi_pnn
real*8     :: Qjac_Up_visco_par_UR_for_Uphi_p, Qjac_Up_visco_par_UR_for_Uphi_n, Qjac_Up_visco_par_UR_for_Uphi_pnn
real*8     :: Qjac_Up_visco_parpar_UR_for_Uphi_p, Qjac_Up_visco_parpar_UR_for_Uphi_n, Qjac_Up_visco_parpar_UR_for_Uphi_pnn
real*8     :: Qjac_Up_visco_par_UZ_for_Uphi_p, Qjac_Up_visco_par_UZ_for_Uphi_n, Qjac_Up_visco_par_UZ_for_Uphi_pnn
real*8     :: Qjac_Up_visco_parpar_UZ_for_Uphi_p, Qjac_Up_visco_parpar_UZ_for_Uphi_n, Qjac_Up_visco_parpar_UZ_for_Uphi_pnn
real*8     :: Qjac_Up_visco_par_Up_for_Uphi_p, Qjac_Up_visco_par_Up_for_Uphi_n, Qjac_Up_visco_par_Up_for_Uphi_pnn
real*8     :: Qjac_Up_visco_parpar_Up_for_Uphi_p, Qjac_Up_visco_parpar_Up_for_Uphi_n, Qjac_Up_visco_parpar_Up_for_Uphi_pnn

real*8     :: visco_pp_diff ! = visco_par_par - visco_par
real*8     :: visco_par_prof, visco_pp_diff_prof
visco_pp_diff = visco_par_par - visco_par
#endif

! --- Switches for numerical stability of resistive and diamagnetic terms in AR and AZ equations
eta_ARAZ  = 0.d0  ! =0.0 to switch off resistive   terms for AR and AZ equations
tauIC_ARAZ= 0.d0  ! =0.0 to switch off diamagnetic terms for AR and AZ equations
tauIC_AP  = 0.d0  ! =0.0 to switch off diamagnetic terms for A_phi equation
tauIC_rho = 0.d0
tauIC_tem = 0.d0
if (eta_ARAZ_on  ) eta_ARAZ   = 1.d0 ! switched on by default
if (tauIC_ARAZ_on) tauIC_ARAZ = 1.d0 ! switched on by default
if (tauIC_AP_on)   tauIC_AP   = 1.d0 ! switched on by default
if (tauIC_rho_on)  tauIC_rho  = 1.d0 ! switched on by default
if (tauIC_tem_on)  tauIC_tem  = 1.d0 ! switched on by default

! --- Initialise
ELM_p  = 0.d0
ELM_n  = 0.d0
ELM_k  = 0.d0
ELM_kn = 0.d0
ELM_pnn = 0.d0
#ifdef ENABLE_VISCO_PAR_AND_PARPAR
ELM_knn  = 0.d0
#endif
RHS_p  = 0.d0
RHS_k  = 0.d0
ELM    = 0.d0
RHS    = 0.d0

rhs_p_ij  = 0.d0
rhs_k_ij  = 0.d0
amat      = 0.d0
Pjac      = 0.d0
Qjac_p    = 0.d0
Qjac_k    = 0.d0
Qjac_n    = 0.d0
Qjac_kn   = 0.d0
Qjac_pnn  = 0.d0
#ifdef ENABLE_VISCO_PAR_AND_PARPAR
Qjac_knn  = 0.d0
#endif
Pvec_prev = 0.d0
Qvec_p    = 0.d0
Qvec_k    = 0.d0

! --- Implicit scheme
theta = time_evol_theta
!zeta  = time_evol_zeta
! change zeta for variable dt
zeta  = time_evol_zeta * 2.0d0 * tstep / (tstep + tstep_prev)

! --- Do we need to use the FFT or non-FFT version?
if ( (i_tor_min == 1) .and. (i_tor_max == n_tor) ) then
  ! In case of global matrix construction:
  use_fft = n_tor > n_tor_fft_thresh
else
  ! In case of "direct construction" of harmonic matrix never FFT:
  use_fft = .false.
end if

if ( use_fft ) then
  ! In case of FFT, don't loop over toroidal harmonics:
  n_tor_start = 1
  n_tor_end   = 1
else
  n_tor_start = i_tor_min
  n_tor_end   = i_tor_max
end if

n_tor_local = n_tor_end - n_tor_start +1

! --- Toroidal basis functions
if (use_fft) then
  ! --- Not needed in case of FFT
  HHZ    = 1.d0
  HHZ_p  = 1.d0
  HHZ_pp = 1.d0
else
  ! --- Needed in case of non-FFT
  do mp=1,n_plane
    do in = 1,n_tor
      HHZ   (in,mp) = HZ   (in,mp)
      HHZ_p (in,mp) = HZ_p (in,mp)
      HHZ_pp(in,mp) = HZ_pp(in,mp)
    enddo
  enddo
endif

! --- RZ variables, Equations variables, and GS-Equilibrium variables
x_g  = 0.d0; x_s  = 0.d0; x_t  = 0.d0; x_ss = 0.d0; x_st = 0.d0; x_tt = 0.d0 
y_g  = 0.d0; y_s  = 0.d0; y_t  = 0.d0; y_ss = 0.d0; y_st = 0.d0; y_tt = 0.d0
rho_axisym  = 0.d0 ; Ti_axisym  = 0.d0; Te_axisym  = 0.d0

eq_g = 0.d0; eq_s = 0.d0; eq_t = 0.d0; eq_p = 0.d0; eq_ss = 0.d0; eq_st = 0.d0; eq_tt = 0.d0
psi_axisym = 0.d0 ; psi_axisym_s = 0.d0 ; psi_axisym_t = 0.d0
eq_pp    = 0.d0   ; eq_sp    = 0.d0     ; eq_tp=0.d0
delta_g = 0.d0; delta_s = 0.d0; delta_t = 0.d0
Fprofile= 0.d0
do i=1,n_vertex_max
  do j=1,n_degrees
    do mt=1, n_gauss
      do ms=1, n_gauss

        x_g(ms,mt)  = x_g(ms,mt)  + nodes(i)%x(1,j,1) * element%size(i,j) * H(i,j,ms,mt)
        x_s(ms,mt)  = x_s(ms,mt)  + nodes(i)%x(1,j,1) * element%size(i,j) * H_s(i,j,ms,mt)
        x_t(ms,mt)  = x_t(ms,mt)  + nodes(i)%x(1,j,1) * element%size(i,j) * H_t(i,j,ms,mt)

        x_ss(ms,mt) = x_ss(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_ss(i,j,ms,mt)
        x_st(ms,mt) = x_st(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_st(i,j,ms,mt)
        x_tt(ms,mt) = x_tt(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_tt(i,j,ms,mt)

        y_g(ms,mt)  = y_g(ms,mt)  + nodes(i)%x(1,j,2) * element%size(i,j) * H(i,j,ms,mt)
        y_s(ms,mt)  = y_s(ms,mt)  + nodes(i)%x(1,j,2) * element%size(i,j) * H_s(i,j,ms,mt)
        y_t(ms,mt)  = y_t(ms,mt)  + nodes(i)%x(1,j,2) * element%size(i,j) * H_t(i,j,ms,mt)

        y_ss(ms,mt) = y_ss(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_ss(i,j,ms,mt)
        y_st(ms,mt) = y_st(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_st(i,j,ms,mt)
        y_tt(ms,mt) = y_tt(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_tt(i,j,ms,mt)

        Fprofile  (ms,mt) = Fprofile  (ms,mt) + nodes(i)%Fprof_eq(j) * element%size(i,j) * H  (i,j,ms,mt)

        ! --- Equilibrium psi (n=0 only) for sources
        rho_axisym(ms,mt) = rho_axisym(ms,mt) + nodes(i)%values(1,j,var_rho) * element%size(i,j) * H(i,j,ms,mt)  * HZ(1,1)
        Ti_axisym(ms,mt)  = Ti_axisym(ms,mt)  + nodes(i)%values(1,j,var_Ti ) * element%size(i,j) * H(i,j,ms,mt)  * HZ(1,1)
        Te_axisym(ms,mt)  = Te_axisym(ms,mt)  + nodes(i)%values(1,j,var_Te ) * element%size(i,j) * H(i,j,ms,mt)  * HZ(1,1)
        psi_axisym(ms,mt) = psi_axisym(ms,mt) + nodes(i)%values(1,j,var_A3) * element%size(i,j) * H(i,j,ms,mt)  * HZ(1,1)
        psi_axisym_s(ms,mt) = psi_axisym_s(ms,mt) + nodes(i)%values(1,j,var_A3) * element%size(i,j) * H_s(i,j,ms,mt)
        psi_axisym_t(ms,mt) = psi_axisym_t(ms,mt) + nodes(i)%values(1,j,var_A3) * element%size(i,j) * H_t(i,j,ms,mt)
      end do
    end do

    do mt=1, n_gauss
      do ms=1, n_gauss

        do k=1,n_var

          do mp=1,n_plane
            do in=1,n_tor
              ! --- store variables
              eq_g(mp,k,ms,mt) = eq_g(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H(i,j,ms,mt)  * HZ(in,mp)
              eq_s(mp,k,ms,mt) = eq_s(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_s(i,j,ms,mt)* HZ(in,mp)
              eq_t(mp,k,ms,mt) = eq_t(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_t(i,j,ms,mt)* HZ(in,mp)
              eq_p(mp,k,ms,mt) = eq_p(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H(i,j,ms,mt)  * HZ_p(in,mp)

              eq_ss(mp,k,ms,mt) = eq_ss(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_ss(i,j,ms,mt)* HZ(in,mp)
              eq_st(mp,k,ms,mt) = eq_st(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_st(i,j,ms,mt)* HZ(in,mp)
              eq_tt(mp,k,ms,mt) = eq_tt(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_tt(i,j,ms,mt)* HZ(in,mp)


              eq_pp(mp,k,ms,mt) = eq_pp(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H(i,j,ms,mt)  * HZ_pp(in,mp)
              eq_sp(mp,k,ms,mt) = eq_sp(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_s(i,j,ms,mt)* HZ_p(in,mp)
              eq_tp(mp,k,ms,mt) = eq_tp(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_t(i,j,ms,mt)* HZ_p(in,mp)

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

! approximate estimate of the element length h_e
! needed for Shock capturing stabilzation
midp_edge1(:) =  0.5d0 * ( nodes(1)%x(1,1,:) + nodes(2)%x(1,1,:) )
midp_edge2(:) =  0.5d0 * ( nodes(2)%x(1,1,:) + nodes(3)%x(1,1,:) )
midp_edge3(:) =  0.5d0 * ( nodes(3)%x(1,1,:) + nodes(4)%x(1,1,:) )
midp_edge4(:) =  0.5d0 * ( nodes(4)%x(1,1,:) + nodes(1)%x(1,1,:) )

len1 = sqrt( (midp_edge1(1)-midp_edge3(1))**2 + (midp_edge1(2)-midp_edge3(2))**2)
len2 = sqrt( (midp_edge2(1)-midp_edge4(1))**2 + (midp_edge2(2)-midp_edge4(2))**2)
h_e = dmin1(len1, len2)

! --- Sources
! --- Note about the current sources:
! --- The toroidal current source can be taken from the routine current.f90, as usual.
! --- The JR and JZ current sources, however, need to be calculated using the initial Grad-Shafranov equilibrium
! --- At t=0, with GS equilibrium, we have:
! --- JR = + psi_Z * dF_dpsi / R = + (psi_bnd_init - psi_axis_init) * psi_norm_Z * dF_dpsi / R
! --- JZ = - psi_R * dF_dpsi / R = - (psi_bnd_init - psi_axis_init) * psi_norm_R * dF_dpsi / R
! --- Thus, these are our current sources as a funcion of psi_norm. Now, at any time, the psi-map will change
! --- And so denormalising psi_norm again, we get:
! --- JR = + (psi_bnd_init - psi_axis_init) / (psi_bnd - psi_axis) * psi_Z * dF_dpsi / R
! --- JZ = - (psi_bnd_init - psi_axis_init) / (psi_bnd - psi_axis) * psi_R * dF_dpsi / R
current_source_JR = 0.d0
current_source_JZ = 0.d0
current_source_Jp = 0.d0
particle_source   = 0.d0
heat_source_i     = 0.d0
heat_source_e     = 0.d0
V_source          = 0.d0
dV_dpsi_source    = 0.d0
dV_dz_source      = 0.d0
dF_dpsi           = 0.d0
do mt=1, n_gauss
  do ms=1, n_gauss
    ! --- These are just temporary, only for the sources (they are recalculated afterwards at each n_plane)
    R     = x_g(ms,mt)
    Z     = y_g(ms,mt)
    xjac  = x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)
    ! --- This will be used in the main loop!
    psi_axisym_R(ms,mt) = (   y_t(ms,mt) * psi_axisym_s(ms,mt)  - y_s(ms,mt) * psi_axisym_t(ms,mt) ) / xjac
    psi_axisym_Z(ms,mt) = ( - x_t(ms,mt) * psi_axisym_s(ms,mt)  + x_s(ms,mt) * psi_axisym_t(ms,mt) ) / xjac
    ! --- These are just temporary, only for the sources (they are recalculated afterwards at each n_plane)
    A30   = eq_g(1,var_A3,ms,mt)
    A30_s = eq_s(1,var_A3,ms,mt)
    A30_t = eq_t(1,var_A3,ms,mt)
    A30_R = (   y_t(ms,mt) * A30_s  - y_s(ms,mt) * A30_t ) / xjac
    A30_Z = ( - x_t(ms,mt) * A30_s  + x_s(ms,mt) * A30_t ) / xjac
    ! --- The dF_dpsi function calculated on time-dependent psi_norm
    ! --- Note: Fprof is be taken from the node values (cleaner)
    call F_profile(xpoint2, xcase2, Z, Z_xpoint, psi_axisym(ms,mt), psi_axis, psi_bnd, &
                   Fprof_time_dep,dF_dpsi(ms,mt) ,dF_dz      ,dF_dpsi2 ,dF_dz2      ,dF_dpsi_dz , &
                   zFFprime ,dFFprime_dpsi,dFFprime_dz,dFFprime_dpsi2,dFFprime_dz2,dFFprime_dpsi_dz)
    if (keep_current_prof) then
      ! --- Toroidal current source. Historically JOREK uses a negative current, so we need to reverse it.
      call current(xpoint2, xcase2, R,Z, Z_xpoint, psi_axisym(ms,mt),psi_axis,psi_bnd,current_source_Jp(ms,mt))
      current_source_Jp(ms,mt) = - current_source_Jp(ms,mt)
      ! --- Poloidal current sources
      current_source_JR(ms,mt) = + (ES%psi_bnd_init - ES%psi_axis_init) / (psi_bnd - psi_axis) * psi_axisym_Z(ms,mt) * dF_dpsi(ms,mt) / R
      current_source_JZ(ms,mt) = - (ES%psi_bnd_init - ES%psi_axis_init) / (psi_bnd - psi_axis) * psi_axisym_R(ms,mt) * dF_dpsi(ms,mt) / R
      if (eta_ARAZ_simple) then
        current_source_JR(ms,mt) = 0.d0
        current_source_JZ(ms,mt) = 0.d0
      endif
    endif
    call sources(xpoint2, xcase2, Z, Z_xpoint, psi_axisym(ms,mt),psi_axis,psi_bnd,particle_source(ms,mt),heat_source_i(ms,mt),heat_source_e(ms,mt))
    ! --- Adaptive source to make sure density/temperature stays exactly at input profiles
    psi_norm = get_psi_n(psi_axisym(ms,mt), y_g(ms,mt))
    if (particlesource_fix_profile) then
      if (psi_norm .lt. particlesource_fix_profile_psi_max) then
        call density    (xpoint2, xcase2, Z, Z_xpoint, psi_axisym(ms,mt),psi_axis,psi_bnd, rho_GS,dd1,dd2,dd3,dd4,dd5,dd6,dd7,dd8)
        if ( (rho_axisym(ms,mt) .lt. rho_GS) .or. particlesource_fix_profile_allow_negative ) &
          particle_source(ms,mt) = (1.0/particlesource_fix_profile_tsteps) * (rho_GS-rho_axisym(ms,mt)) ! / tstep
      endif
    endif
    if (heatsource_i_fix_profile) then
      if (psi_norm .lt. heatsource_i_fix_profile_psi_max) then
        call temperature_i(xpoint2, xcase2, Z, Z_xpoint, psi_axisym(ms,mt),psi_axis,psi_bnd, Ti_GS  ,dd1,dd2,dd3,dd4,dd5,dd6,dd7,dd8)
        if ( (Ti_axisym(ms,mt)   .lt. Ti_GS  ) .or. heatsource_i_fix_profile_allow_negative ) &
          heat_source_i(ms,mt)     = (1.0/heatsource_i_fix_profile_tsteps) * (Ti_GS  -Ti_axisym(ms,mt)  ) ! / tstep
      endif
    endif
    if (heatsource_e_fix_profile) then
      if (psi_norm .lt. heatsource_e_fix_profile_psi_max) then
        call temperature_e(xpoint2, xcase2, Z, Z_xpoint, psi_axisym(ms,mt),psi_axis,psi_bnd, Te_GS  ,dd1,dd2,dd3,dd4,dd5,dd6,dd7,dd8)
        if ( (Te_axisym(ms,mt)   .lt. Te_GS  ) .or. heatsource_e_fix_profile_allow_negative ) &
          heat_source_e(ms,mt)     = (1.0/heatsource_e_fix_profile_tsteps) * (Te_GS  -Te_axisym(ms,mt)  ) ! / tstep
      endif
    endif

#ifndef TEMPERATURE_DEPENDENT_UPDATE
    call temperature_i(xpoint2, xcase2, Z, Z_xpoint, psi_axisym(ms,mt),psi_axis,psi_bnd, Ti_GS  ,dd1,dd2,dd3,dd4,dd5,dd6,dd7,dd8)
    Ti_GS_initial(ms,mt) = max(Ti_GS, 1.d-12)
    dTi_GS_dR_initial(ms,mt) = dd1 * psi_axisym_R(ms,mt)
    dTi_GS_dZ_initial(ms,mt) = dd1 * psi_axisym_Z(ms,mt)
    call temperature_e(xpoint2, xcase2, Z, Z_xpoint, psi_axisym(ms,mt),psi_axis,psi_bnd, Te_GS  ,dd1,dd2,dd3,dd4,dd5,dd6,dd7,dd8)
    Te_GS_initial(ms,mt) = max(Te_GS, 1.d-12)
    dTe_GS_dR_initial(ms,mt) = dd1 * psi_axisym_R(ms,mt)
    dTe_GS_dZ_initial(ms,mt) = dd1 * psi_axisym_Z(ms,mt)
#endif

    ! --- Bootstrap current 
    if (bootstrap) then
      Ti0    = eq_g(1,var_Ti,ms,mt)
      Ti0_s  = eq_s(1,var_Ti,ms,mt)
      Ti0_t  = eq_t(1,var_Ti,ms,mt)
      Ti0_R  = (   y_t(ms,mt) * Ti0_s  - y_s(ms,mt) * Ti0_t ) / xjac
      Ti0_Z  = ( - x_t(ms,mt) * Ti0_s  + x_s(ms,mt) * Ti0_t ) / xjac
      Te0    = eq_g(1,var_Te,ms,mt)
      Te0_s  = eq_s(1,var_Te,ms,mt)
      Te0_t  = eq_t(1,var_Te,ms,mt)
      Te0_R  = (   y_t(ms,mt) * Te0_s  - y_s(ms,mt) * Te0_t ) / xjac
      Te0_Z  = ( - x_t(ms,mt) * Te0_s  + x_s(ms,mt) * Te0_t ) / xjac
      rho0   = eq_g(1,var_rho,ms,mt)
      rho0_s = eq_s(1,var_rho,ms,mt)
      rho0_t = eq_t(1,var_rho,ms,mt)
      rho0_R = (   y_t(ms,mt) * rho0_s  - y_s(ms,mt) * rho0_t ) / xjac
      rho0_Z = ( - x_t(ms,mt) * rho0_s  + x_s(ms,mt) * rho0_t ) / xjac
      call bootstrap_current(R, Z,                                                        &
                             R_axis,   Z_axis,   psi_axis,                                &
                             R_xpoint, Z_xpoint, psi_bnd, psi_norm,                       &
                             psi_axisym(ms,mt), psi_axisym_R(ms,mt), psi_axisym_Z(ms,mt), &
                             rho0, rho0_R, rho0_Z,                                        &
                             Ti0,  Ti0_R,  Ti0_Z,                                         &
                             Te0,  Te0_R,  Te0_Z,                                         &
                             Jb)
      ! --- Full Sauter formula for initial profiles
      call density      (xpoint2, xcase2, Z, Z_xpoint, psi_axisym(ms,mt),psi_axis,psi_bnd, &
                         rho_initial(ms,mt),dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz)
      call temperature_i(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, eq_g(1,1,ms,mt),psi_axis,psi_bnd, &
                         Ti_initial (ms,mt),dTi_dpsi,dTi_dz,dTi_dpsi2,dTi_dz2,dTi_dpsi_dz,dTi_dpsi3,dTi_dpsi_dz2, dTi_dpsi2_dz)
      call temperature_e(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, eq_g(1,1,ms,mt),psi_axis,psi_bnd, &
                         Te_initial (ms,mt),dTe_dpsi,dTe_dz,dTe_dpsi2,dTe_dz2,dTe_dpsi_dz,dTe_dpsi3,dTe_dpsi_dz2, dTe_dpsi2_dz)
      Ti0   = Ti_initial(ms,mt)
      Ti0_R = dTi_dpsi * psi_axisym_R(ms,mt)
      Ti0_Z = dTi_dpsi * psi_axisym_Z(ms,mt)
      Te0   = Te_initial(ms,mt)
      Te0_R = dTe_dpsi * psi_axisym_R(ms,mt)
      Te0_Z = dTe_dpsi * psi_axisym_Z(ms,mt)
      rho0   = rho_initial(ms,mt)
      rho0_R = dn_dpsi * psi_axisym_R(ms,mt)
      rho0_Z = dn_dpsi * psi_axisym_Z(ms,mt)
      call bootstrap_current(R, Z,                                                        &
                             R_axis,   Z_axis,   psi_axis,                                &
                             R_xpoint, Z_xpoint, psi_bnd, psi_norm,                       &
                             psi_axisym(ms,mt), psi_axisym_R(ms,mt), psi_axisym_Z(ms,mt), &
                             rho0, rho0_R, rho0_Z,                                        &
                             Ti0,  Ti0_R,  Ti0_Z,                                         &
                             Te0,  Te0_R,  Te0_Z,                                         &
                             Jb_0)
      ! --- Subtract the initial equilibrium part
      Jb = Jb - Jb_0
      ! --- Historically JOREK uses a negative current, so we need to reverse it.
      Jb = -Jb
      ! --- Add bootstrap to main current source
      current_source_Jp(ms,mt) = current_source_Jp(ms,mt) + Jb
    endif
    ! --- Neoclassical friction
    if ( NEO ) then 
      if (num_neo_file) then
        call neo_coef(xpoint2, xcase2, Z, Z_xpoint, psi_axisym(ms,mt),psi_axis,psi_bnd, amu_neo_prof(ms,mt), aki_neo_prof(ms,mt))
      else
        amu_neo_prof(ms,mt) = amu_neo_const
        aki_neo_prof(ms,mt) = aki_neo_const
      endif
      psi_norm = get_psi_n(psi_axisym(ms,mt), Z)
      if (psi_norm .gt. 0.995) amu_neo_prof(ms,mt) = 0.d0
      if (psi_norm .lt. 0.100) amu_neo_prof(ms,mt) = 0.d0
    endif
    ! --- Source of toroidal velocity
    if ( ( abs(V_0) .ge. 1.e-12 ) .or. ( num_rot ) ) then
      call velocity(xpoint2, xcase2, Z, Z_xpoint, psi_axisym(ms,mt), psi_axis, psi_bnd, &
                    V_source(ms,mt), dV_dpsi_source(ms,mt),dV_dz_source(ms,mt),         &
                    dV_dpsi2,dV_dz2,dV_dpsi_dz,dV_dpsi3,dV_dpsi_dz2, dV_dpsi2_dz)
    else
      V_source(ms,mt)       = 0.d0
      dV_dpsi_source(ms,mt) = 0.d0
      dV_dz_source(ms,mt)   = 0.d0
    endif

  enddo
enddo


! --- Main loops
do j=1,n_degrees
  do i=1,n_vertex_max
    ELM_p(:,:,1:n_var)  = 0
    ELM_n(:,:,1:n_var)  = 0
    ELM_k(:,:,1:n_var)  = 0
    ELM_kn(:,:,1:n_var) = 0

    do mt=1, n_gauss
      do ms=1, n_gauss

        wst = wgauss(ms)*wgauss(mt)

        xjac    = x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)
        xjac_R  = (x_ss(ms,mt)*y_t(ms,mt)**2 - y_ss(ms,mt)*x_t(ms,mt)*y_t(ms,mt) - 2.d0*x_st(ms,mt)*y_s(ms,mt)*y_t(ms,mt)   &
                + y_st(ms,mt)*(x_s(ms,mt)*y_t(ms,mt) + x_t(ms,mt)*y_s(ms,mt))                                               &
                + x_tt(ms,mt)*y_s(ms,mt)**2 - y_tt(ms,mt)*x_s(ms,mt)*y_s(ms,mt)) / xjac
        xjac_Z  = (y_tt(ms,mt)*x_s(ms,mt)**2 - x_tt(ms,mt)*y_s(ms,mt)*x_s(ms,mt) - 2.d0*y_st(ms,mt)*x_t(ms,mt)*x_s(ms,mt)   &
                + x_st(ms,mt)*(y_t(ms,mt)*x_s(ms,mt) + y_s(ms,mt)*x_t(ms,mt))                                               &
                + y_ss(ms,mt)*x_t(ms,mt)**2 - x_ss(ms,mt)*y_t(ms,mt)*x_t(ms,mt)) / xjac

        R = x_g(ms,mt)
        Z = y_g(ms,mt)

        do mp = 1, n_plane

          ! --- AR
          AR0   = eq_g(mp,var_AR,ms,mt)
          AR0_p = eq_p(mp,var_AR,ms,mt)
          AR0_s = eq_s(mp,var_AR,ms,mt)
          AR0_t = eq_t(mp,var_AR,ms,mt)
          AR0_R = (   y_t(ms,mt) * AR0_s  - y_s(ms,mt) * AR0_t ) / xjac
          AR0_Z = ( - x_t(ms,mt) * AR0_s  + x_s(ms,mt) * AR0_t ) / xjac
          AR0_ss = eq_ss(mp,var_AR,ms,mt)
          AR0_st = eq_st(mp,var_AR,ms,mt)
          AR0_tt = eq_tt(mp,var_AR,ms,mt)
          AR0_RR = (AR0_ss * y_t(ms,mt)**2 - 2.d0*AR0_st * y_s(ms,mt)*y_t(ms,mt) + AR0_tt * y_s(ms,mt)**2 &
                 + AR0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) ) &
                 + AR0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) ) / xjac**2             &
                 - xjac_R * (AR0_s* y_t(ms,mt) - AR0_t * y_s(ms,mt))  / xjac**2
          AR0_ZZ = (AR0_ss * x_t(ms,mt)**2 - 2.d0*AR0_st * x_s(ms,mt)*x_t(ms,mt) + AR0_tt * x_s(ms,mt)**2 &
                 + AR0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) ) &
                 + AR0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) ) / xjac**2             &
                 - xjac_Z * (- AR0_s * x_t(ms,mt) + AR0_t * x_s(ms,mt) )  / xjac**2

          AR0_sp = eq_sp(mp,var_AR,ms,mt)
          AR0_tp = eq_tp(mp,var_AR,ms,mt)
          AR0_pp = eq_pp(mp,var_AR,ms,mt)
          AR0_Rp = (   y_t(ms,mt) * AR0_sp  - y_s(ms,mt) * AR0_tp ) / xjac
          AR0_Zp = ( - x_t(ms,mt) * AR0_sp  + x_s(ms,mt) * AR0_tp ) / xjac
          AR0_RZ = (- AR0_ss * y_t(ms,mt)*x_t(ms,mt) - AR0_tt * x_s(ms,mt)*y_s(ms,mt) &
                 + AR0_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  ) &
                 - AR0_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) ) &
                 - AR0_t * (x_st(ms,mt)*y_s(ms,mt)  - x_ss(ms,mt)*y_t(ms,mt) ) ) / xjac**2             &
                 - xjac_R * (- AR0_s * x_t(ms,mt) + AR0_t * x_s(ms,mt) )   / xjac**2

          ! --- AZ
          AZ0   = eq_g(mp,var_AZ,ms,mt)
          AZ0_p = eq_p(mp,var_AZ,ms,mt)
          AZ0_s = eq_s(mp,var_AZ,ms,mt)
          AZ0_t = eq_t(mp,var_AZ,ms,mt)
          AZ0_R = (   y_t(ms,mt) * AZ0_s  - y_s(ms,mt) * AZ0_t ) / xjac
          AZ0_Z = ( - x_t(ms,mt) * AZ0_s  + x_s(ms,mt) * AZ0_t ) / xjac
          AZ0_ss = eq_ss(mp,var_AZ,ms,mt)
          AZ0_st = eq_st(mp,var_AZ,ms,mt)
          AZ0_tt = eq_tt(mp,var_AZ,ms,mt)
          AZ0_RR = (AZ0_ss * y_t(ms,mt)**2 - 2.d0*AZ0_st * y_s(ms,mt)*y_t(ms,mt) + AZ0_tt * y_s(ms,mt)**2 &
                 + AZ0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) ) &
                 + AZ0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) ) / xjac**2             &
                 - xjac_R * (AZ0_s* y_t(ms,mt) - AZ0_t * y_s(ms,mt))  / xjac**2
          AZ0_ZZ = (AZ0_ss * x_t(ms,mt)**2 - 2.d0*AZ0_st * x_s(ms,mt)*x_t(ms,mt) + AZ0_tt * x_s(ms,mt)**2 &
                 + AZ0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) ) &
                 + AZ0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) ) / xjac**2             &
                 - xjac_Z * (- AZ0_s * x_t(ms,mt) + AZ0_t * x_s(ms,mt) )  / xjac**2

          AZ0_sp = eq_sp(mp,var_AZ,ms,mt)
          AZ0_tp = eq_tp(mp,var_AZ,ms,mt)
          AZ0_pp = eq_pp(mp,var_AZ,ms,mt)
          AZ0_Rp = (   y_t(ms,mt) * AZ0_sp  - y_s(ms,mt) * AZ0_tp ) / xjac
          AZ0_Zp = ( - x_t(ms,mt) * AZ0_sp  + x_s(ms,mt) * AZ0_tp ) / xjac
          AZ0_RZ = (- AZ0_ss * y_t(ms,mt)*x_t(ms,mt) - AZ0_tt * x_s(ms,mt)*y_s(ms,mt) &
                 + AZ0_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  ) &
                 - AZ0_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) ) &
                 - AZ0_t * (x_st(ms,mt)*y_s(ms,mt)  - x_ss(ms,mt)*y_t(ms,mt) ) ) / xjac**2             &
                 - xjac_R * (- AZ0_s * x_t(ms,mt) + AZ0_t * x_s(ms,mt) )   / xjac**2

          ! --- A3:=psi is defined as: A = ... + A03 * grad(phi)
          ! --- as opposed to the magnetic and velocity fieds that are defined as
          ! --- B = ... + Bp * e_phi             and     V = ... + Vp * e_phi       ie.
          ! --- B = ... + Bp * R * grad(phi)     and     V = ... + Vp * R * grad(phi)
          A30   = eq_g(mp,var_A3,ms,mt)
          A30_p = eq_p(mp,var_A3,ms,mt)
          A30_s = eq_s(mp,var_A3,ms,mt)
          A30_t = eq_t(mp,var_A3,ms,mt)
          A30_R = (   y_t(ms,mt) * A30_s  - y_s(ms,mt) * A30_t ) / xjac
          A30_Z = ( - x_t(ms,mt) * A30_s  + x_s(ms,mt) * A30_t ) / xjac
          A30_ss = eq_ss(mp,var_A3,ms,mt)
          A30_st = eq_st(mp,var_A3,ms,mt)
          A30_tt = eq_tt(mp,var_A3,ms,mt)
          A30_RR = (A30_ss * y_t(ms,mt)**2 - 2.d0*A30_st * y_s(ms,mt)*y_t(ms,mt) + A30_tt * y_s(ms,mt)**2 &
                 + A30_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) ) &
                 + A30_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) ) / xjac**2             &
                 - xjac_R * (A30_s* y_t(ms,mt) - A30_t * y_s(ms,mt))  / xjac**2
          A30_ZZ = (A30_ss * x_t(ms,mt)**2 - 2.d0*A30_st * x_s(ms,mt)*x_t(ms,mt) + A30_tt * x_s(ms,mt)**2 &
                 + A30_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) ) &
                 + A30_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) ) / xjac**2             &
                 - xjac_Z * (- A30_s * x_t(ms,mt) + A30_t * x_s(ms,mt) )  / xjac**2

          A30_sp = eq_sp(mp,var_A3,ms,mt)
          A30_tp = eq_tp(mp,var_A3,ms,mt)
          A30_pp = eq_pp(mp,var_A3,ms,mt)
          A30_Rp = (   y_t(ms,mt) * A30_sp  - y_s(ms,mt) * A30_tp ) / xjac
          A30_Zp = ( - x_t(ms,mt) * A30_sp  + x_s(ms,mt) * A30_tp ) / xjac
          A30_RZ = (- A30_ss * y_t(ms,mt)*x_t(ms,mt) - A30_tt * x_s(ms,mt)*y_s(ms,mt)                     &
                 + A30_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  ) &
                 - A30_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) ) &
                 - A30_t * (x_st(ms,mt)*y_s(ms,mt)  - x_ss(ms,mt)*y_t(ms,mt) ) ) / xjac**2             &
                 - xjac_R * (- A30_s * x_t(ms,mt) + A30_t * x_s(ms,mt) )   / xjac**2

          ! --- UR
          UR0   = eq_g(mp,var_UR,ms,mt)
          UR0_p = eq_p(mp,var_UR,ms,mt)
          UR0_s = eq_s(mp,var_UR,ms,mt)
          UR0_t = eq_t(mp,var_UR,ms,mt)
          UR0_R = (   y_t(ms,mt) * UR0_s  - y_s(ms,mt) * UR0_t ) / xjac
          UR0_Z = ( - x_t(ms,mt) * UR0_s  + x_s(ms,mt) * UR0_t ) / xjac
          UR0_ss = eq_ss(mp,var_UR,ms,mt)
          UR0_st = eq_st(mp,var_UR,ms,mt)
          UR0_tt = eq_tt(mp,var_UR,ms,mt)
          UR0_RR =  (   UR0_ss * y_t(ms,mt)**2 - 2.d0*UR0_st * y_s(ms,mt)*y_t(ms,mt)              &
                      + UR0_tt * y_s(ms,mt)**2                                                    &
                      + UR0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                &
                      + UR0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) )   )/ xjac**2   &
                      - UR0_R * xjac_R / xjac
          UR0_ZZ =  (   UR0_ss * x_t(ms,mt)**2 - 2.d0*UR0_st * x_s(ms,mt)*x_t(ms,mt)              &
                      + UR0_tt * x_s(ms,mt)**2                                                    &
                      + UR0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                &
                      + UR0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) )   )/ xjac**2   &
                      - UR0_Z * xjac_Z / xjac

          ! --- UZ
          UZ0   = eq_g(mp,var_UZ,ms,mt)
          UZ0_p = eq_p(mp,var_UZ,ms,mt)
          UZ0_s = eq_s(mp,var_UZ,ms,mt)
          UZ0_t = eq_t(mp,var_UZ,ms,mt)
          UZ0_R = (   y_t(ms,mt) * UZ0_s  - y_s(ms,mt) * UZ0_t ) / xjac
          UZ0_Z = ( - x_t(ms,mt) * UZ0_s  + x_s(ms,mt) * UZ0_t ) / xjac
          UZ0_ss = eq_ss(mp,var_UZ,ms,mt)
          UZ0_st = eq_st(mp,var_UZ,ms,mt)
          UZ0_tt = eq_tt(mp,var_UZ,ms,mt)
          UZ0_RR =  (   UZ0_ss * y_t(ms,mt)**2 - 2.d0*UZ0_st * y_s(ms,mt)*y_t(ms,mt)              &
                      + UZ0_tt * y_s(ms,mt)**2                                                    &
                      + UZ0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                &
                      + UZ0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) )   )/ xjac**2   &
                      - UZ0_R * xjac_R / xjac
          UZ0_ZZ =  (   UZ0_ss * x_t(ms,mt)**2 - 2.d0*UZ0_st * x_s(ms,mt)*x_t(ms,mt)              &
                      + UZ0_tt * x_s(ms,mt)**2                                                    &
                      + UZ0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                &
                      + UZ0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) )   )/ xjac**2   &
                      - UZ0_Z * xjac_Z / xjac

          ! --- Up is defined u0_phi : V = .. + Up0 * e_phi (physical component)
          Up0   = eq_g(mp,var_Up,ms,mt)
          Up0_p = eq_p(mp,var_Up,ms,mt)
          Up0_s = eq_s(mp,var_Up,ms,mt)
          Up0_t = eq_t(mp,var_Up,ms,mt)
          Up0_R = (   y_t(ms,mt) * Up0_s  - y_s(ms,mt) * Up0_t ) / xjac
          Up0_Z = ( - x_t(ms,mt) * Up0_s  + x_s(ms,mt) * Up0_t ) / xjac
          Up0_ss = eq_ss(mp,var_Up,ms,mt)
          Up0_st = eq_st(mp,var_Up,ms,mt)
          Up0_tt = eq_tt(mp,var_Up,ms,mt)
          Up0_RR =  (   Up0_ss * y_t(ms,mt)**2 - 2.d0*Up0_st * y_s(ms,mt)*y_t(ms,mt)              &
                      + Up0_tt * y_s(ms,mt)**2                                                    &
                      + Up0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                &
                      + Up0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) )   )/ xjac**2   &
                      - Up0_R * xjac_R / xjac
          Up0_ZZ =  (   Up0_ss * x_t(ms,mt)**2 - 2.d0*Up0_st * x_s(ms,mt)*x_t(ms,mt)              &
                      + Up0_tt * x_s(ms,mt)**2                                                    &
                      + Up0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                &
                      + Up0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) )   )/ xjac**2   &
                      - Up0_Z * xjac_Z / xjac

          ! --- rho
          rho0      = eq_g(mp,var_rho,ms,mt)
          rho0_corr = max(rho0,1.d-12)!corr_neg_dens1(rho0) ! CAREFUL! FULL-MHD DOESN'T LIKE THE CORR FUNCTIONS AT ALL
          drho0_corr_dn = 0.d0!dcorr_neg_dens_drho(rho0)
          rho0_p    = eq_p(mp,var_rho,ms,mt)
          rho0_s    = eq_s(mp,var_rho,ms,mt)
          rho0_t    = eq_t(mp,var_rho,ms,mt)
          rho0_R    = (   y_t(ms,mt) * rho0_s  - y_s(ms,mt) * rho0_t ) / xjac
          rho0_Z    = ( - x_t(ms,mt) * rho0_s  + x_s(ms,mt) * rho0_t ) / xjac
          rho0_ss   = eq_ss(mp,var_rho,ms,mt)
          rho0_st   = eq_st(mp,var_rho,ms,mt)
          rho0_tt   = eq_tt(mp,var_rho,ms,mt)
          rho0_RR   = ( rho0_ss * y_t(ms,mt)**2 - 2.d0*rho0_st * y_s(ms,mt)*y_t(ms,mt)              &
                      + rho0_tt * y_s(ms,mt)**2 &
                      + rho0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                &
                      + rho0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) )   )/ xjac**2   &
                      - rho0_R * xjac_R / xjac
          rho0_ZZ   = ( rho0_ss * x_t(ms,mt)**2 - 2.d0*rho0_st * x_s(ms,mt)*x_t(ms,mt)              &
                      + rho0_tt * x_s(ms,mt)**2 &
                      + rho0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                &
                      + rho0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) )   )/ xjac**2   &
                      - rho0_Z * xjac_Z / xjac

          ! --- Ti
          Ti0      = eq_g(mp,var_Ti,ms,mt)
          Ti0_corr = max(Ti0,1.d-12)!corr_neg_temp1(Ti0) ! CAREFUL! FULL-MHD DOESN'T LIKE THE CORR FUNCTIONS AT ALL
          dTi0_corr_dT = 0.d0 !dcorr_neg_temp_dT(Ti0) ! Improve the correction
          Ti0_p    = eq_p(mp,var_Ti,ms,mt)
          Ti0_s    = eq_s(mp,var_Ti,ms,mt)
          Ti0_t    = eq_t(mp,var_Ti,ms,mt)
          Ti0_R    = (   y_t(ms,mt) * Ti0_s  - y_s(ms,mt) * Ti0_t ) / xjac
          Ti0_Z    = ( - x_t(ms,mt) * Ti0_s  + x_s(ms,mt) * Ti0_t ) / xjac
          Ti0_ss   = eq_ss(mp,var_Ti,ms,mt)
          Ti0_st   = eq_st(mp,var_Ti,ms,mt)
          Ti0_tt   = eq_tt(mp,var_Ti,ms,mt)
          Ti0_RR   = ( Ti0_ss * y_t(ms,mt)**2 - 2.d0*Ti0_st * y_s(ms,mt)*y_t(ms,mt)              &
                       + Ti0_tt * y_s(ms,mt)**2 &
                       + Ti0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                &
                       + Ti0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) )   )/ xjac**2   &
                       - Ti0_R * xjac_R / xjac
          Ti0_ZZ   = ( Ti0_ss * x_t(ms,mt)**2 - 2.d0*Ti0_st * x_s(ms,mt)*x_t(ms,mt)              &
                       + Ti0_tt * x_s(ms,mt)**2 &
                       + Ti0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                &
                       + Ti0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) )   )/ xjac**2   &
                       - Ti0_Z * xjac_Z / xjac

          ! --- Te
          Te0      = eq_g(mp,var_Te,ms,mt)
          Te0_corr = max(Te0,1.d-12)!corr_neg_temp1(Te0) ! CAREFUL! FULL-MHD DOESN'T LIKE THE CORR FUNCTIONS AT ALL
          dTe0_corr_dT = 0.d0 !dcorr_neg_temp_dT(Te0) ! Improve the correction
          Te0_p    = eq_p(mp,var_Te,ms,mt)
          Te0_s    = eq_s(mp,var_Te,ms,mt)
          Te0_t    = eq_t(mp,var_Te,ms,mt)
          Te0_R    = (   y_t(ms,mt) * Te0_s  - y_s(ms,mt) * Te0_t ) / xjac
          Te0_Z    = ( - x_t(ms,mt) * Te0_s  + x_s(ms,mt) * Te0_t ) / xjac
          Te0_ss   = eq_ss(mp,var_Te,ms,mt)
          Te0_st   = eq_st(mp,var_Te,ms,mt)
          Te0_tt   = eq_tt(mp,var_Te,ms,mt)
          Te0_RR   = ( Te0_ss * y_t(ms,mt)**2 - 2.d0*Te0_st * y_s(ms,mt)*y_t(ms,mt)              &
                       + Te0_tt * y_s(ms,mt)**2 &
                       + Te0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                &
                       + Te0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) )   )/ xjac**2   &
                       - Te0_R * xjac_R / xjac
          Te0_ZZ   = ( Te0_ss * x_t(ms,mt)**2 - 2.d0*Te0_st * x_s(ms,mt)*x_t(ms,mt)              &
                       + Te0_tt * x_s(ms,mt)**2 &
                       + Te0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                &
                       + Te0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) )   )/ xjac**2   &
                       - Te0_Z * xjac_Z / xjac

          ! --- P
          pi0      = rho0 * Ti0
          pi0_corr = rho0_corr * Ti0_corr
          pi0_R    = rho0_R * Ti0 + rho0 * Ti0_R
          pi0_Z    = rho0_Z * Ti0 + rho0 * Ti0_Z
          pi0_s    = rho0_s * Ti0 + rho0 * Ti0_s
          pi0_t    = rho0_t * Ti0 + rho0 * Ti0_t
          pi0_p    = rho0_p * Ti0 + rho0 * Ti0_p

          ! --- P
          pe0      = rho0 * Te0
          pe0_corr = rho0_corr * Te0_corr
          pe0_R    = rho0_R * Te0 + rho0 * Te0_R
          pe0_Z    = rho0_Z * Te0 + rho0 * Te0_Z
          pe0_s    = rho0_s * Te0 + rho0 * Te0_s
          pe0_t    = rho0_t * Te0 + rho0 * Te0_t
          pe0_p    = rho0_p * Te0 + rho0 * Te0_p

          ! --- P
          p0      = rho0 * (Ti0+Te0)
          p0_corr = rho0_corr * (Ti0_corr+Te0_corr)
          p0_R    = rho0_R * (Ti0+Te0) + rho0 * (Ti0_R+Te0_R)
          p0_Z    = rho0_Z * (Ti0+Te0) + rho0 * (Ti0_Z+Te0_Z)
          p0_s    = rho0_s * (Ti0+Te0) + rho0 * (Ti0_s+Te0_s)
          p0_t    = rho0_t * (Ti0+Te0) + rho0 * (Ti0_t+Te0_t)
          p0_p    = rho0_p * (Ti0+Te0) + rho0 * (Ti0_p+Te0_p)

          ! --- psi_norm
          psi_norm = get_psi_n(psi_axisym(ms,mt), y_g(ms,mt))

          ! --- Diffusions
          D_par_prof = D_par
          D_prof   = get_dperp (psi_norm)
          ZKi_prof = get_zk_iperp(psi_norm)
          ZKe_prof = get_zk_eperp(psi_norm)
          visco_par_prof = visco_par
          visco_pp_diff_prof = visco_pp_diff
#ifdef TEMPERATURE_DEPENDENT_UPDATE
          ! --- Resistivity
          if ( eta_T_dependent .and. Te0_corr <= T_max_eta) then
            eta_T     = eta   * (Te0_corr/Te_0)**(-1.5d0)
            deta_dT   = - eta   * (1.5d0)  * Te0_corr**(-2.5d0) * Te_0**(1.5d0)
            d2eta_d2T =   eta   * (3.75d0) * Te0_corr**(-3.5d0) * Te_0**(1.5d0)
          else if ( eta_T_dependent .and. Te0_corr > T_max_eta) then
            eta_T     = eta   * (T_max_eta/Te_0)**(-1.5d0)
            deta_dT   = 0.
            d2eta_d2T = 0.     
          else
            eta_T     = eta
            deta_dT   = 0.d0
            d2eta_d2T = 0.d0
          end if
          if ( eta_T_dependent .and.  xpoint2 .and. (Te0 .lt. T_min) ) then
              eta_T     = eta    * (T_min/Te_0)**(-1.5d0)
              deta_dT   = 0.d0
              d2eta_d2T = 0.d0
          end if
          eta_R = deta_dT * Te0_R
          eta_Z = deta_dT * Te0_Z
          eta_p = deta_dT * Te0_p

          ! --- Eta for ohmic heating ! NOT YET IMPLEMENTED !!! NEED THE CURRENT !!!
          if ( eta_T_dependent .and. Te0_corr <= T_max_eta_ohm) then
            eta_T_ohm     = eta_ohmic   * (Te0_corr/Te_0)**(-1.5d0)
            deta_dT_ohm   = - eta_ohmic   * (1.5d0)  * Te0_corr**(-2.5d0) * Te_0**(1.5d0)
          else if ( eta_T_dependent .and. Te0_corr > T_max_eta_ohm) then
            eta_T_ohm     = eta_ohmic   * (T_max_eta_ohm/Te_0)**(-1.5d0)
            deta_dT_ohm   = 0.    
          else
            eta_T_ohm     = eta_ohmic
            deta_dT_ohm   = 0.d0
          end if

          ! --- Viscosity
          if ( visco_T_dependent ) then
            visco_T     =   visco * (Te0_corr/Te_0)**(-1.5d0)
            dvisco_dT   = - visco * (1.5d0)  * Te0_corr**(-2.5d0) * Te_0**(1.5d0)
            d2visco_d2T =   visco * (3.75d0) * Te0_corr**(-3.5d0) * Te_0**(1.5d0)
            if ( xpoint2 .and. (Te0 .lt. T_min) ) then
              visco_T     = visco  * (T_min/Te_0)**(-1.5d0)
              dvisco_dT   = 0.d0
              d2visco_d2T = 0.d0
            endif
          else
            visco_T     = visco
            dvisco_dT   = 0.d0
            d2visco_d2T = 0.d0
          end if
          visco_R = dvisco_dT * Te0_R
          visco_Z = dvisco_dT * Te0_Z
          visco_p = dvisco_dT * Te0_p

          ! --- Temperature dependent parallel heat diffusivity
          if ( ZKpar_T_dependent ) then
            ZKi_par_T   = ZK_i_par * (Ti0_corr/Ti_0)**(+2.5d0)              ! temperature dependent parallel conductivity
            dZKi_par_dT = ZK_i_par * (2.5d0)  * Ti0_corr**(+1.5d0) * Ti_0**(-2.5d0) * dTi0_corr_dT
            if (ZKi_par_T .gt. ZK_par_max) then
              ZKi_par_T   = Zk_par_max
              dZKi_par_dT = 0.d0
            endif
            if ( xpoint2 .and. (Ti0 .lt. T_min) ) then
              ZKi_par_T   = ZK_i_par * (T_min/Ti_0)**(+2.5d0)
              dZKi_par_dT = 0.d0
            endif

            ZKe_par_T   = ZK_e_par * (Te0_corr/Te_0)**(+2.5d0)              ! temperature dependent parallel conductivity
            dZKe_par_dT = ZK_e_par * (2.5d0)  * Te0_corr**(+1.5d0) * Te_0**(-2.5d0) * dTe0_corr_dT
            if (ZKe_par_T .gt. ZK_par_max) then
              ZKe_par_T   = Zk_par_max
              dZKe_par_dT = 0.d0
            endif
            if ( xpoint2 .and. (Te0 .lt. T_min) ) then
              ZKe_par_T   = ZK_e_par * (T_min/Te_0)**(+2.5d0)
              dZKe_par_dT = 0.d0
            endif
          else
            ZKi_par_T   = ZK_i_par                                            ! parallel conductivity
            dZKi_par_dT = 0.d0
            ZKe_par_T   = ZK_e_par                                            ! parallel conductivity
            dZKe_par_dT = 0.d0
          endif
#else
          Te_GS = Te_GS_initial(ms,mt)
          Ti_GS = Ti_GS_initial(ms,mt)
          ! --- Resistivity
          if ( eta_T_dependent .and. Te_GS <= T_max_eta) then
            eta_T     = eta   * (Te_GS/Te_0)**(-1.5d0)
            deta_dT   = - eta   * (1.5d0)  * Te_GS**(-2.5d0) * Te_0**(1.5d0)
            d2eta_d2T =   eta   * (3.75d0) * Te_GS**(-3.5d0) * Te_0**(1.5d0)
          else if ( eta_T_dependent .and. Te_GS > T_max_eta) then
            eta_T     = eta   * (T_max_eta/Te_0)**(-1.5d0)
            deta_dT   = 0.
            d2eta_d2T = 0.     
          else
            eta_T     = eta
            deta_dT   = 0.d0
            d2eta_d2T = 0.d0
          end if
          if ( eta_T_dependent .and.  xpoint2 .and. (Te_GS .lt. T_min) ) then
              eta_T     = eta    * (T_min/Te_0)**(-1.5d0)
              deta_dT   = 0.d0
              d2eta_d2T = 0.d0
          end if
          eta_R = deta_dT * dTe_GS_dR_initial(ms,mt)
          eta_Z = deta_dT * dTe_GS_dZ_initial(ms,mt)
          eta_p = deta_dT * 0.d0

          ! --- Eta for ohmic heating ! NOT YET IMPLEMENTED !!! NEED THE CURRENT !!!
          if ( eta_T_dependent .and. Te_GS <= T_max_eta_ohm) then
            eta_T_ohm     = eta_ohmic   * (Te_GS/Te_0)**(-1.5d0)
            deta_dT_ohm   = - eta_ohmic   * (1.5d0)  * Te_GS**(-2.5d0) * Te_0**(1.5d0)
          else if ( eta_T_dependent .and. Te_GS > T_max_eta_ohm) then
            eta_T_ohm     = eta_ohmic   * (T_max_eta_ohm/Te_0)**(-1.5d0)
            deta_dT_ohm   = 0.    
          else
            eta_T_ohm     = eta_ohmic
            deta_dT_ohm   = 0.d0
          end if

          ! --- Viscosity
          if ( visco_T_dependent ) then
            visco_T     =   visco * (Te_GS/Te_0)**(-1.5d0)
            dvisco_dT   = - visco * (1.5d0)  * Te_GS**(-2.5d0) * Te_0**(1.5d0)
            d2visco_d2T =   visco * (3.75d0) * Te_GS**(-3.5d0) * Te_0**(1.5d0)
            if ( xpoint2 .and. (Te_GS .lt. T_min) ) then
              visco_T     = visco  * (T_min/Te_0)**(-1.5d0)
              dvisco_dT   = 0.d0
              d2visco_d2T = 0.d0
            endif
          else
            visco_T     = visco
            dvisco_dT   = 0.d0
            d2visco_d2T = 0.d0
          end if
          visco_R = dvisco_dT * dTe_GS_dR_initial(ms,mt)
          visco_Z = dvisco_dT * dTe_GS_dZ_initial(ms,mt)
          visco_p = dvisco_dT * 0.d0

          ! --- Temperature dependent parallel heat diffusivity
          if ( ZKpar_T_dependent ) then
            ZKi_par_T   = ZK_i_par * (Ti_GS/Ti_0)**(+2.5d0)              ! temperature dependent parallel conductivity
            dZKi_par_dT = ZK_i_par * (2.5d0)  * Ti_GS**(+1.5d0) * Ti_0**(-2.5d0) * dTi0_corr_dT
            if (ZKi_par_T .gt. ZK_par_max) then
              ZKi_par_T   = Zk_par_max
              dZKi_par_dT = 0.d0
            endif
            if ( xpoint2 .and. (Ti_GS .lt. T_min) ) then
              ZKi_par_T   = ZK_i_par * (T_min/Ti_0)**(+2.5d0)
              dZKi_par_dT = 0.d0
            endif

            ZKe_par_T   = ZK_e_par * (Te_GS/Te_0)**(+2.5d0)              ! temperature dependent parallel conductivity
            dZKe_par_dT = ZK_e_par * (2.5d0)  * Te_GS**(+1.5d0) * Te_0**(-2.5d0) * dTe0_corr_dT
            if (ZKe_par_T .gt. ZK_par_max) then
              ZKe_par_T   = Zk_par_max
              dZKe_par_dT = 0.d0
            endif
            if ( xpoint2 .and. (Te_GS .lt. T_min) ) then
              ZKe_par_T   = ZK_e_par * (T_min/Te_0)**(+2.5d0)
              dZKe_par_dT = 0.d0
            endif
          else
            ZKi_par_T   = ZK_i_par                                            ! parallel conductivity
            dZKi_par_dT = 0.d0
            ZKe_par_T   = ZK_e_par                                            ! parallel conductivity
            dZKe_par_dT = 0.d0
          endif
#endif
          ! --- Ion-electron energy transfer
          if (thermalization) then
            ! Te in eV:
            Te_corr_eV     = Te0_corr/(EL_CHG*MU_ZERO*central_density*1.d20)
            dTe_corr_eV_dT = dTe0_corr_dT/(EL_CHG*MU_ZERO*central_density*1.d20)

            ne_SI          = rho0_corr * 1.d20 * central_density ! electron density (SI)
            if (ne_SI < 1.d16) ne_SI = 1.d16 ! To prevent absurd number in the coulomb lambda

            lambda_e_bg  = 23. - log((ne_SI*1.d-6)**0.5*Te_corr_eV**(-1.5)) ! Assuming bg_charge is 1! 
            nu_e_bg      = 1.8d-19*(1.d6*MASS_ELECTRON*MASS_PROTON*central_mass) ** 0.5&
                           * (1.d14*central_density*rho0_corr) * lambda_e_bg &
                           / (1.d3*(MASS_ELECTRON*Ti0_corr+Te0_corr*MASS_PROTON*central_mass)&
                           / (EL_CHG * MU_ZERO * central_density * 1.d20)) ** 1.5 ! Assuming bg_charge is 1!

            if (nu_e_bg < 0.)  nu_e_bg  = 0.

            !Converting the energy transfer rate from s^-1 to JOREK unit
            t_norm   = sqrt(MU_ZERO * central_mass * MASS_PROTON * central_density * 1.d20)
            nu_e_bg  = nu_e_bg * t_norm    

            dTe_i    = nu_e_bg * (Ti0_corr - Te0_corr) * rho0_corr
            dTi_e    = -dTe_i

            !Calculating the density and temperature derivative for amats
            !We negelect the coulomb log's dericatives due to their smallness
            ! IMPORTANT NOTE: in full-MHD these derivatives are very unstable for some reason
            !                 this may be for the same reason that correction functions cannot be used for
            !                 the density and temperatures. Will need to be investigated in the future

            dnu_e_bg_dTi    = 0.d0!-1.5*MASS_ELECTRON*nu_e_bg*dTi0_corr_dT / (MASS_ELECTRON*Ti0_corr + MASS_PROTON*central_mass*Te0_corr)
            dnu_e_bg_dTe    = 0.d0!-1.5*MASS_PROTON*central_mass*nu_e_bg*dTe0_corr_dT &
                                  !/ (MASS_ELECTRON*Ti0_corr + MASS_PROTON*central_mass*Te0_corr)

            dnu_e_bg_drho   = 0.d0!nu_e_bg * drho0_corr_dn / rho0_corr

            ddTe_i_dTi      = 0.d0!dnu_e_bg_dTi * (Ti0_corr - Te0_corr) * rho0_corr + nu_e_bg * dTi0_corr_dT * rho0_corr
            ddTe_i_dTe      = 0.d0!dnu_e_bg_dTe * (Ti0_corr - Te0_corr) * rho0_corr - nu_e_bg * dTe0_corr_dT * rho0_corr
            ddTe_i_drho     = 0.d0!dnu_e_bg_drho * (Ti0_corr - Te0_corr) * rho0_corr &
                              !+ nu_e_bg * (Ti0_corr - Te0_corr) * drho0_corr_dn
        
            ddTi_e_dTi      = -ddTe_i_dTi
            ddTi_e_dTe      = -ddTe_i_dTe
            ddTi_e_drho     = -ddTe_i_drho

          else
            dTe_i       = 0.d0
            dTi_e       = 0.d0
            ddTe_i_dTi  = 0.d0
            ddTe_i_dTe  = 0.d0
            ddTe_i_drho = 0.d0
            ddTi_e_dTi  = 0.d0
            ddTi_e_dTe  = 0.d0
            ddTi_e_drho = 0.d0
          endif


          ! --- Magnetic field
          Fprof = Fprofile(ms,mt)
          BR0 = ( A30_Z - AZ0_p )/ R
          BZ0 = ( AR0_p - A30_R )/ R
          Bp0 = ( AZ0_R - AR0_Z )    + Fprof / R
          BB2 = Bp0**2 + BR0**2 + BZ0**2

          ! --- RESISTIVITY SWITCHES FOR AR AND AZ EQUATIONS
          ! --- 1.
          ! --- Default set-up is eta_ARAZ_on = .true.
          ! --- In this case, the same resistivity is used for all AR,AZ,A3 components
          ! --- 2.
          ! --- If (eta_ARAZ_on = .true.) and (eta_ARAZ_simple = .true.), then the Fprof component of Bphi is removed from the resistive term.
          ! --- Note that in a stationary equilibrium, the Fprof component of Bphi from the current and the current source should exactly cancel anyway.
          ! --- However, this Fprof component can lead to strong noise at high resistivity, so it is often better to remove it.
          ! --- 3.
          ! --- If (eta_ARAZ_on = .false.), then resistivity is switched off by default for the AR and AZ equations.
          ! --- However, regardless of which eta-model is used for A3, (ie. eta_T_dependent or not), you can still set eta_ARAZ_const 
          ! --- which will use a constant resistivity at the level eta_ARAZ_const.
          ! --- 4.
          ! --- Here again, setting (eta_ARAZ_simple = .true.) with eta_ARAZ_const removes the Fprof dependence of Bphi in the resistive term and the current term for AR&AZ

          ! --- Toroidal magnetic field without F contribution (for eta_ARAZ_simple)
          Bp00 = Bp0
          if (eta_ARAZ_simple) Bp00 = ( AZ0_R - AR0_Z )

          ! --- B.grad and V.grad
          BgradTi   = BR0 * Ti0_R   + BZ0 * Ti0_Z   + Bp0 * Ti0_p   / R
          BgradTe   = BR0 * Te0_R   + BZ0 * Te0_Z   + Bp0 * Te0_p   / R
          UgradTi   = UR0 * Ti0_R   + UZ0 * Ti0_Z   + Up0 * Ti0_p   / R
          UgradTe   = UR0 * Te0_R   + UZ0 * Te0_Z   + Up0 * Te0_p   / R
          BgradRho = BR0 * rho0_R + BZ0 * rho0_Z + Bp0 * rho0_p / R
          UgradRho = UR0 * rho0_R + UZ0 * rho0_Z + Up0 * rho0_p / R
          BgradPe  = rho0*BgradTe + Te0*BgradRho

          ! --- V.grad.V
          UgradUR  = UR0 * UR0_R + UZ0 * UR0_Z + Up0 * UR0_p / R
          UgradUZ  = UR0 * UZ0_R + UZ0 * UZ0_Z + Up0 * UZ0_p / R
          UgradUp  = UR0 * Up0_R + UZ0 * Up0_Z + Up0 * Up0_p / R

          ! --- div.V
          divU     = UR0_R + UR0/R + UZ0_Z + Up0_p/R
          divRhoU  = rho0 * divU + UgradRho

          ! --- Diamagnetic velocity Vdia = tau / (rho*BB2) * B x grad(p)
          ! --- Note-1: Phi component defined as physical component VdiaP*e_phi, like V and B
          ! --- Note-2: Factor of F0 is here so that we have the same definition of tau_IC in RMHD and FMHD
          tau_IC = tauIC
          if (tauIC_tanh) then ! suppress the diamagnetic term near the boundary
            ! distance_bnd = 1.d10
            ! do im=1,n_vertex_max
            !   if (nodes(im)%boundary .eq. 1) then
            !     distance_bnd = min(distance_bnd, sqrt((R-nodes(im)%x(1,1,1))**2 + (Z-nodes(im)%x(1,1,2))**2) )
            !   endif
            ! enddo
            distance_bnd = 1.d0 - dsqrt(psi_norm)
            tau_IC = tauIC * (0.5d0 - 0.5d0 * tanh(-(distance_bnd - tauIC_shift)/tauIC_steep) )
          endif
          VdiaR0 = tau_IC*F0 / (R * rho0_corr * BB2) * (  BZ0*pi0_p - R*Bp0*pi0_Z)
          VdiaZ0 = tau_IC*F0 / (R * rho0_corr * BB2) * (R*BP0*pi0_R -   BR0*pi0_p)
          VdiaP0 = tau_IC*F0 / (    rho0_corr * BB2) * (  BR0*pi0_Z -   BZ0*pi0_R)

          ! --- Vdia.grad.V
          VdiaGradUR  = VdiaR0 * UR0_R + VdiaZ0 * UR0_Z + VdiaP0 * UR0_p / R
          VdiaGradUZ  = VdiaR0 * UZ0_R + VdiaZ0 * UZ0_Z + VdiaP0 * UZ0_p / R
          VdiaGradUp  = VdiaR0 * Up0_R + VdiaZ0 * Up0_Z + VdiaP0 * Up0_p / R
          
          ! --- Toroidal velocity
          Vt0   = V_source(ms,mt)
          Vt0_R = dV_dpsi_source(ms,mt) * psi_axisym_R(ms,mt)
          Vt0_Z = dV_dz_source(ms,mt) + dV_dpsi_source(ms,mt) * psi_axisym_Z(ms,mt)

          ! --- Neoclassical friction, we assume that the magnetic field is from the axisymmetric component only
          if (NEO) then
            Btht  = sqrt(BR0**2 + BZ0**2)
            Btht  = max(1.d-10,Btht)
            Vtht  = ( BR0*(UR0+VdiaR0) + BZ0*(UZ0+VdiaZ0) ) / Btht
            Vneo  = +aki_neo_prof(ms,mt) * tau_IC*F0 / BB2 / Btht * ( BR0*Bp0*Ti0_Z - BZ0*Bp0*Ti0_R)
            PneoR = +amu_neo_prof(ms,mt) * rho0 * BB2/Btht**3 * BR0 * (Vtht - Vneo)
            PneoZ = +amu_neo_prof(ms,mt) * rho0 * BB2/Btht**3 * BZ0 * (Vtht - Vneo)
          else
            PneoR = 0.d0
            PneoZ = 0.d0
          endif
          
          psieq_R = (   y_t(ms,mt) * psi_axisym_s(ms,mt)  - y_s(ms,mt) * psi_axisym_t(ms,mt) ) / xjac
          psieq_Z = ( - x_t(ms,mt) * psi_axisym_s(ms,mt)  + x_s(ms,mt) * psi_axisym_t(ms,mt) ) / xjac

          BR0_R = -1.0d0/R**2 * ( A30_Z - AZ0_p ) + ( A30_RZ - AZ0_Rp )/R
          BR0_Z = ( A30_ZZ - AZ0_Zp ) / R
          BR0_p = ( A30_Zp - AZ0_pp ) / R
          BZ0_R = -1.0d0/R**2 * ( AR0_p - A30_R ) + ( AR0_Rp - A30_RR )/R
          BZ0_Z = ( AR0_Zp - A30_RZ ) / R
          BZ0_p = ( AR0_pp - A30_Rp ) / R
          Bp0_R = ( AZ0_RR - AR0_RZ )  &
                + (psi_bnd - psi_axis) / (ES%psi_bnd_init - ES%psi_axis_init) * psieq_R * dF_dpsi(ms,mt) / R  &
                - Fprofile(ms,mt)/R**2
          Bp0_Z = ( AZ0_RZ - AR0_ZZ ) &
                + (psi_bnd - psi_axis) / (ES%psi_bnd_init - ES%psi_axis_init) * psieq_Z * dF_dpsi(ms,mt) / R
          Bp0_p = ( AZ0_Rp - AR0_Zp )

          JR0 = Bp0_Z - BZ0_p / R
          JZ0 = (BR0_p - R*Bp0_R - Bp0) / R
          Jp0 = BZ0_R - BR0_Z

          JJ2 = JR0*JR0 + JZ0*JZ0 + Jp0*Jp0

          ! For shock capturing stabilization
          tau_sc = 0.d0
          if (use_sc) call calculate_sc_quantities()

          do im=n_tor_start, n_tor_end

            ! --- test functions (V*)
            v   = H(i,j,ms,mt)   * element%size(i,j) * HHZ(im,mp)
            v_s = H_s(i,j,ms,mt) * element%size(i,j) * HHZ(im,mp)
            v_t = H_t(i,j,ms,mt) * element%size(i,j) * HHZ(im,mp)
            v_p__k = H(i,j,ms,mt)   * element%size(i,j) * HHZ_p(im,mp)
            v_R = (  y_t(ms,mt) * v_s - y_s(ms,mt) * v_t ) / xjac
            v_Z = (- x_t(ms,mt) * v_s + x_s(ms,mt) * v_t ) / xjac 
            v_ss= H_ss(i,j,ms,mt) * element%size(i,j) * HHZ(im,mp)
            v_st= H_st(i,j,ms,mt) * element%size(i,j) * HHZ(im,mp)
            v_tt= H_tt(i,j,ms,mt) * element%size(i,j) * HHZ(im,mp)
            v_RR= ( v_ss * y_t(ms,mt)**2 - 2.d0*v_st * y_s(ms,mt)*y_t(ms,mt) &
                  + v_tt * y_s(ms,mt)**2 &
                  + v_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                &
                  + v_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) )   )/ xjac**2   &
                  - v_R * xjac_R / xjac
            v_ZZ= ( v_ss * x_t(ms,mt)**2 - 2.d0*v_st * x_s(ms,mt)*x_t(ms,mt) &
                  + v_tt * x_s(ms,mt)**2 &
                  + v_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                &
                  + v_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) )   )/ xjac**2   &
                  - v_Z * xjac_Z / xjac

            ! --- Terms for Integration-by-parts
            BgradVstar__p = BR0 * v_R + BZ0 * v_Z
            BgradVstar__k = Bp0 * v_p__k / R

            UgradVstar__p = UR0 * v_R + UZ0 * v_Z
            UgradVstar__k = Up0 * v_p__k / R

            gradRho_gradVstar__p = rho0_R * v_R  + rho0_Z * v_Z
            gradRho_gradVstar__k = (rho0_p / R) * (v_p__k  / R)

            gradTi_gradVstar__p = Ti0_R * v_R  + Ti0_Z * v_Z
            gradTi_gradVstar__k = (Ti0_p / R) * (v_p__k  / R)

            gradTe_gradVstar__p = Te0_R * v_R  + Te0_Z * v_Z
            gradTe_gradVstar__k = (Te0_p / R) * (v_p__k  / R)

            ! --- (V.grad)V        normal part                       diamagnetic part
            Qconv_UR = - v * rho0 * ( UgradUR - Up0**2  / R )    - v * rho0 * ( VdiaGradUR - VdiaP0*Up0 / R )
            Qconv_UZ = - v * rho0 * UgradUZ                      - v * rho0 * VdiaGradUZ
            Qconv_Up = - v * rho0 * ( UgradUp + UR0*Up0 / R )    - v * rho0 * ( VdiaGradUp + VdiaP0*UR0 / R )

            JxB_UR__p = - BR0 * BgradVstar__p - Bp0**2 * v / R + ( v_R + v / R ) * ( 0.5d0 * BB2 )
            JxB_UR__k = - BR0 * BgradVstar__k
            JxB_UZ__p = - BZ0 * BgradVstar__p + v_Z * ( 0.5d0 * BB2)
            JxB_UZ__k = - BZ0 * BgradVstar__k
            JxB_Up__p = - Bp0 * BgradVstar__p + Bp0*BR0* v / R
            JxB_Up__k = - Bp0 * BgradVstar__k + ( v_p__k/ R ) * ( 0.5d0 * BB2 )

            ! --- Diamagnetic terms
            VdiaGradVstar__p = VdiaR0 * v_R + VdiaZ0 * v_Z
            VdiaGradVstar__k = VdiaP0 * v_p__k / R

            ! --- Viscous terms (The correct viscosity (there is only one...))
            Qvisc_UR__p = 0.d0 ; Qvisc_UR__k = 0.d0
            Qvisc_UZ__p = 0.d0 ; Qvisc_UZ__k = 0.d0
            Qvisc_Up__p = 0.d0 ; Qvisc_Up__k = 0.d0
            Qvisc_T     = 0.d0

            Qgrad_visc_dot_grad_UR_p = 0.d0
            Qgrad_visc_dot_grad_UZ_p = 0.d0
            Qgrad_visc_dot_grad_Up_p = 0.d0
            Qvisc_grad_BR_dot_grad_UR_p = 0.d0
            Qvisc_grad_BZ_dot_grad_UZ_p = 0.d0
            Qvisc_grad_Bp_dot_grad_Up_p = 0.d0
            QdivU_BdotGradVisc = 0.d0; curl_U0_dot_gradmu_x_B0 = 0.d0; curl_U0_dot_curl_B0 = 0.d0

            Qvisc_gradmu_dot_grad_vecV_UR_p = 0.d0
            Qvisc_gradmu_dot_grad_vecV_UZ_p = 0.d0
            Qvisc_gradmu_dot_grad_vecV_Up_p = 0.d0

            Qgrad_visc_dot_gradU_T_UR__p = 0.d0
            Qgrad_visc_dot_gradU_T_UZ__p = 0.d0
            Qgrad_visc_dot_gradU_T_Up__p = 0.d0

            Q_minus_2_over_3_divU_grad_visc_UR__p = 0.d0
            Q_minus_2_over_3_divU_grad_visc_UZ__p = 0.d0
            Q_minus_2_over_3_divU_grad_visc_Up__p = 0.d0

            Q_1_over_3_visc_grad_divU_UR__p = 0.d0
            Q_1_over_3_visc_grad_divU_UZ__p = 0.d0
            Q_1_over_3_visc_grad_divU_Up__p = 0.d0
            Q_1_over_3_visc_grad_divU_Up__k = 0.d0

            Q_1_over_3_B_dot_visc_grad_divU__p = 0.d0
            Q_1_over_3_B_dot_visc_grad_divU__k = 0.d0

#ifdef VISCOSITY_TENSOR_FORM
            Qvisc_UR__p = - (v_R + v/R) * divU + v_Z * (UZ0_R - UR0_Z)
            Qvisc_UR__k = - v_p__k/R**2 * (UR0_p - Up0 - R*Up0_R)

            Qvisc_UZ__p = - v_Z * divU - v_R * (UZ0_R - UR0_Z)
            Qvisc_UZ__k = + v_p__k/R**2 * (R*Up0_Z - UZ0_p)

            Qvisc_Up__p = - v_Z/R * (R*Up0_Z - UZ0_p) + (v_R + v/R)/R * (UR0_p - Up0 - R*Up0_R)
            Qvisc_Up__k = - v_p__k/R * divU

            Qgrad_visc_dot_grad_UR_p = - (visco_R + visco_T / R) * divU + visco_Z * (UZ0_R-UR0_Z) - 1.d0 / R**2 * visco_p * (UR0_p - Up0 - R * Up0_R)
            Qgrad_visc_dot_grad_UZ_p = - visco_Z * divU + 1.d0 / R**2 * visco_p * (R*Up0_Z-UZ0_p) - visco_R * (UZ0_R - UR0_Z)
            Qgrad_visc_dot_grad_Up_p = - 1.d0 / R * visco_p * divU + 1.d0 / R * (visco_R + visco_T / R) * (UR0_p - Up0 - R*Up0_R) - 1.d0 / R * visco_Z * (R * Up0_Z - UZ0_p)
            
            QdivU_BdotGradVisc = divU * (BR0 * visco_R + BZ0 * visco_Z + Bp0 * visco_p / R)

            curl_U0_dot_gradmu_x_B0 = (Up0_Z - UZ0_p / R) * (visco_Z * Bp0 - 1.d0 / R * visco_p * BZ0) &
                                    + 1.d0 / R**2 * (UR0_p - Up0 - R*Up0_R) * (visco_p * BR0 - R * visco_R * Bp0) &
                                    + (UZ0_R - UR0_Z) * (visco_R * BZ0 - visco_Z * BR0)

            curl_U0_dot_curl_B0 = (Up0_Z - UZ0_p / R) * (Bp0_Z - 1.d0 / R * BZ0_p) &
                                + 1.d0 / R**2 * (UR0_p - Up0 - R*Up0_R) * (BR0_p - Bp0 - R * Bp0_R) &
                                + (UZ0_R - UR0_Z) * (BZ0_R - BR0_Z)
            
            ! --- The toroidal momentum source (same as above replacing all Up0 by -Vt0)
            Qvisc_UR__k = Qvisc_UR__k + v_p__k/R**2 * (- Vt0 - R*Vt0_R)
            Qvisc_UZ__k = Qvisc_UZ__k - v_p__k/R**2 * (R*Vt0_Z)
            Qvisc_Up__p = Qvisc_Up__p + v_Z/R * (R*Vt0_Z) - (v_R + v/R)/R * (- Vt0 - R*Vt0_R)

            ! TODO: implement the toroidal momentum source (all related terms above and below need to be modified)

#else
            ! --- full MHD viscosity
            Qvisc_UR__p = - (v_R * UR0_R + v_Z * UR0_Z) - v * ( UR0 / R**2 + 2.d0 * Up0_p / R**2 ) 
            Qvisc_UR__k = - (v_p__k * UR0_p / R**2)
            Qvisc_UZ__p = - (v_R * UZ0_R + v_Z * UZ0_Z)
            Qvisc_UZ__k = - (v_p__k * UZ0_p / R**2)
            Qvisc_Up__p = - (v_R * Up0_R + v_Z * Up0_Z) - v * ( Up0 / R**2 - 2.d0 * UR0_p / R**2 )
            Qvisc_Up__k = - (v_p__k * Up0_p / R**2)

            Qgrad_visc_dot_grad_UR_p = visco_R * UR0_R + visco_Z * UR0_Z + visco_p * UR0_p / R**2
            Qgrad_visc_dot_grad_UZ_p = visco_R * UZ0_R + visco_Z * UZ0_Z + visco_p * UZ0_p / R**2
            Qgrad_visc_dot_grad_Up_p = visco_R * Up0_R + visco_Z * Up0_Z + visco_p * Up0_p / R**2

            Qvisc_grad_BR_dot_grad_UR_p = BR0_R * UR0_R + BR0_Z * UR0_Z + BR0_p * UR0_p / R**2
            Qvisc_grad_BZ_dot_grad_UZ_p = BZ0_R * UZ0_R + BZ0_Z * UZ0_Z + BZ0_p * UZ0_p / R**2
            Qvisc_grad_Bp_dot_grad_Up_p = Bp0_R * Up0_R + Bp0_Z * Up0_Z + Bp0_p * Up0_p / R**2

#ifdef VISCOSITY_INSIDE_LAPLACIAN
            Qvisc_gradmu_dot_grad_vecV_UR_p = visco_R * UR0_R + visco_Z * UR0_Z + visco_p * UR0_p / R**2 - visco_p * Up0 / R**2
            Qvisc_gradmu_dot_grad_vecV_UZ_p = visco_R * UZ0_R + visco_Z * UZ0_Z + visco_p * UZ0_p / R**2
            Qvisc_gradmu_dot_grad_vecV_Up_p = visco_R * Up0_R + visco_Z * Up0_Z + visco_p * Up0_p / R**2 + visco_p * UR0 / R**2
#endif

#ifdef VISCOSITY_COMPLETE_TERMS
             Qgrad_visc_dot_gradU_T_UR__p = visco_R * UR0_R + visco_Z * UZ0_R + visco_p / R * Up0_R
             Qgrad_visc_dot_gradU_T_UZ__p = visco_R * UR0_Z + visco_Z * UZ0_Z + visco_p / R * Up0_Z
             Qgrad_visc_dot_gradU_T_Up__p = visco_R / R * (UR0_p - Up0) + visco_Z / R * UZ0_p + visco_p / R**2 * (Up0_p + UR0)

             Q_minus_2_over_3_divU_grad_visc_UR__p = - 2.d0/3.d0 * visco_R * divU
             Q_minus_2_over_3_divU_grad_visc_UZ__p = - 2.d0/3.d0 * visco_Z * divU
             Q_minus_2_over_3_divU_grad_visc_Up__p = - 2.d0/3.d0 * visco_p * divU / R

             Q_1_over_3_visc_grad_divU_UR__p = -1.d0/3.d0 * divU * (visco * (v_R + v/R) + v * (visco_R + visco/R))
             Q_1_over_3_visc_grad_divU_UZ__p = -1.d0/3.d0 * divU * (visco * v_Z         + v * visco_Z)
             Q_1_over_3_visc_grad_divU_Up__p = -1.d0/3.d0 * divU * (                    + v * visco_p) / R
             Q_1_over_3_visc_grad_divU_Up__k = -1.d0/3.d0 * divU * (visco * v_p__k                   ) / R

             Q_1_over_3_B_dot_visc_grad_divU__p = -1.d0/3.d0 * divU * (visco * (BR0 * v_R + BZ0 * v_Z) + v * (BR0 * visco_R + BZ0 * visco_Z + Bp0 * visco_p / R))
             Q_1_over_3_B_dot_visc_grad_divU__k = -1.d0/3.d0 * divU * (visco * (Bp0 * v_p__k / R))
#endif

            ! TODO: implement the toroidal momentum source (all related terms above and below need to be modified)
#endif

            ! --- fourth order diffusion terms
            lap_Vstar = v_R / R + v_RR + v_ZZ

            lap_AR  = AR0_R  / R  + AR0_RR   + AR0_ZZ
            lap_AZ  = AZ0_R  / R  + AZ0_RR   + AZ0_ZZ 
            lap_A3  = A30_R  / R  + A30_RR   + A30_ZZ 
            lap_UR  = UR0_R  / R  + UR0_RR   + UR0_ZZ 
            lap_UZ  = UZ0_R  / R  + UZ0_RR   + UZ0_ZZ 
            lap_Up  = Up0_R  / R  + Up0_RR   + Up0_ZZ 
            lap_rho = rho0_R / R  + rho0_RR  + rho0_ZZ
            lap_Ti  = Ti0_R  / R  + Ti0_RR   + Ti0_ZZ  
            lap_Te  = Te0_R  / R  + Te0_RR   + Te0_ZZ  

#ifdef ENABLE_VISCO_PAR_AND_PARPAR
            ! --- For visco_par and visco_par_par
            grad_BRVstar_R =  v * BR0_R + v_R * BR0
            grad_BRVstar_Z =  v * BR0_Z + v_Z * BR0
            grad_BRVstar_phi = (v * BR0_p            )
            grad_BRVstar_phi_k = (            v_p__k * BR0)

            grad_BZVstar_R =  v * BZ0_R + v_R * BZ0
            grad_BZVstar_Z =  v * BZ0_Z + v_Z * BZ0
            grad_BZVstar_phi = (v * BZ0_p            )
            grad_BZVstar_phi_k = (            v_p__k * BZ0)

            grad_BpVstar_R =  v * Bp0_R + v_R * Bp0
            grad_BpVstar_Z =  v * Bp0_Z + v_Z * Bp0
            grad_BpVstar_phi = (v * Bp0_p            ) 
            grad_BpVstar_phi_k = (            v_p__k * Bp0) 

            BB2_R   = 2 * (BR0 * BR0_R + Bp0 * Bp0_R + BZ0 * BZ0_R)
            BB2_Z   = 2 * (BR0 * BR0_Z + Bp0 * Bp0_Z + BZ0 * BZ0_Z)
            BB2_phi = 2 * (BR0 * BR0_p + Bp0 * Bp0_p + BZ0 * BZ0_p)

            B_dot_grad_BB2  = BR0 * BB2_R + BZ0 * BB2_Z + Bp0 * BB2_phi / R
            B_dot_U0        = BR0 * UR0 + BZ0 * UZ0 + Bp0 * Up0
            grad_BdotU0_R   = (BR0_R * UR0 + BZ0_R * UZ0 + Bp0_R * Up0 + BR0 * UR0_R + BZ0 * UZ0_R + Bp0 * Up0_R)
            grad_BdotU0_Z   = (BR0_Z * UR0 + BZ0_Z * UZ0 + Bp0_Z * Up0 + BR0 * UR0_Z + BZ0 * UZ0_Z + Bp0 * Up0_Z)
            grad_BdotU0_phi = (BR0_p * UR0 + BZ0_p * UZ0 + Bp0_p * Up0 + BR0 * UR0_p + BZ0 * UZ0_p + Bp0 * Up0_p)
            B_dot_grad_BdotU0 = BR0 * grad_BdotU0_R &
                              + BZ0 * grad_BdotU0_Z &
                              + Bp0 * grad_BdotU0_phi / R

            B_dot_grad_BR   = BR0 * BR0_R + BZ0 * BR0_Z + Bp0 * BR0_p / R
            B_dot_grad_BZ   = BR0 * BZ0_R + BZ0 * BZ0_Z + Bp0 * BZ0_p / R
            B_dot_grad_Bp   = BR0 * Bp0_R + BZ0 * Bp0_Z + Bp0 * Bp0_p / R

            qparR   = 1.d0 / BB2 * grad_BdotU0_R   - B_dot_U0 / BB2 ** 2 * BB2_R
            qparZ   = 1.d0 / BB2 * grad_BdotU0_Z   - B_dot_U0 / BB2 ** 2 * BB2_Z
            qparPhi = 1.d0 / BB2 * grad_BdotU0_phi / R - B_dot_U0 / BB2 ** 2 * BB2_phi / R

            qparpar    = 1.d0 / BB2 ** 2 * (-1.d0 / BB2 * B_dot_grad_BB2 * B_dot_U0 + B_dot_grad_BdotU0)
            qparparR   = BR0 * qparpar
            qparparZ   = BZ0 * qparpar
            qparparPhi = Bp0 * qparpar

            Qvec_UR_visco_par__p    = qparR    * grad_BRVstar_R + qparZ    * grad_BRVstar_Z + qparPhi    * grad_BRVstar_phi / R 
            Qvec_UR_visco_par__k    = qparPhi    * grad_BRVstar_phi_k / R
            Qvec_UR_visco_parpar__p = qparparR * grad_BRVstar_R + qparparZ * grad_BRVstar_Z + qparparPhi * grad_BRVstar_phi / R
            Qvec_UR_visco_parpar__k = qparparPhi * grad_BRVstar_phi_k / R
            ! define v * q \cdot \nabla B_R, specialized due to (B \cdot U) in the Uphi equation
            Qvec_UR_visco_par_for_Uphi__p    = v * (qparR    * BR0_R + qparZ    * BR0_Z + qparPhi    * BR0_p / R)
            Qvec_UR_visco_parpar_for_Uphi__p = v * (qparparR * BR0_R + qparparZ * BR0_Z + qparparPhi * BR0_p / R)

            Qvec_UZ_visco_par__p    = qparR    * grad_BZVstar_R + qparZ    * grad_BZVstar_Z + qparPhi    * grad_BZVstar_phi / R 
            Qvec_UZ_visco_par__k    = qparPhi    * grad_BZVstar_phi_k / R
            Qvec_UZ_visco_parpar__p = qparparR * grad_BZVstar_R + qparparZ * grad_BZVstar_Z + qparparPhi * grad_BZVstar_phi / R 
            Qvec_UZ_visco_parpar__k = qparparPhi * grad_BZVstar_phi_k / R
            ! define v * q \cdot \nabla B_Z, specialized due to (B \cdot U) in the Uphi equation
            Qvec_UZ_visco_par_for_Uphi__p    = v * (qparR    * BZ0_R + qparZ    * BZ0_Z + qparPhi    * BZ0_p / R)
            Qvec_UZ_visco_parpar_for_Uphi__p = v * (qparparR * BZ0_R + qparparZ * BZ0_Z + qparparPhi * BZ0_p / R)

            Qvec_Up_visco_par__p    = qparR    * grad_BpVstar_R + qparZ    * grad_BpVstar_Z + qparPhi    * grad_BpVstar_phi / R 
            Qvec_Up_visco_par__k    = qparPhi    * grad_BpVstar_phi_k / R
            Qvec_Up_visco_parpar__p = qparparR * grad_BpVstar_R + qparparZ * grad_BpVstar_Z + qparparPhi * grad_BpVstar_phi / R 
            Qvec_Up_visco_parpar__k = qparparPhi * grad_BpVstar_phi_k / R
            ! define v * q \cdot \nabla B_p, specialized due to (B \cdot U) in the Uphi equation
            Qvec_Up_visco_par_for_Uphi__p    = v * (qparR    * Bp0_R + qparZ    * Bp0_Z + qparPhi    * Bp0_p / R)
            Qvec_Up_visco_parpar_for_Uphi__p = v * (qparparR * Bp0_R + qparparZ * Bp0_Z + qparparPhi * Bp0_p / R)

#endif

            rhs_p_ij  = 0.d0
            rhs_k_ij  = 0.d0
            Pvec_prev = 0.d0 ! The time derivative part
            Qvec_p    = 0.d0 ! The rest of the RHS (poloidal part)
            Qvec_k    = 0.d0 ! The rest of the RHS (toroidal part that has phi-derivatives of the test-function)

            !###################################################################################################
            !#  equation 1 (R component induction equation)                                                    #
            !###################################################################################################
            Pvec_prev(var_AR) =   v * delta_g(mp,var_AR,ms,mt)

            Qvec_p(var_AR) = + v * (UZ0 * Bp0 - Up0 * BZ0)                             &
                             + tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2 * BR0 * BgradPe&
                             + eta_ARAZ * v * (eta_Z * Bp00 - eta_p * BZ0 / R )        &
                             - eta_ARAZ * eta_T * ( - v_Z * Bp00)                      &
                             - eta_ARAZ_const   * ( - v_Z * Bp00)                      &
                             + eta_ARAZ * eta_T * v * current_source_JR(ms,mt)         &
                             - eta_ARAZ * eta_num * lap_Vstar * lap_AR
            Qvec_k(var_AR) = - eta_ARAZ * eta_T * ( + v_p__k * BZ0 / R)

            !###################################################################################################
            !#  equation 2 (Z component induction equation)                                                    #
            !###################################################################################################
            Pvec_prev(var_AZ) =   v * delta_g(mp,var_AZ,ms,mt) 

            Qvec_p(var_AZ) = + v * (Up0 * BR0 - UR0 * Bp0)                             &
                             + tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2 * BZ0 * BgradPe&
                             + eta_ARAZ * v * (eta_p / R * BR0 - eta_R * Bp00)         &
                             - eta_ARAZ * eta_T * ( + v_R * Bp00)                      &
                             - eta_ARAZ_const   * ( + v_R * Bp00)                      &
                             + eta_ARAZ * eta_T * v * current_source_JZ(ms,mt)         &
                             - eta_ARAZ * eta_num * lap_Vstar * lap_AZ
            Qvec_k(var_AZ) = - eta_ARAZ * eta_T * ( - v_p__k * BR0 / R)

            !###################################################################################################
            !#  equation 3 (PHI component induction equation)                                                  #
            !###################################################################################################
            Pvec_prev(var_A3) =   v * delta_g(mp,var_A3,ms,mt)

            Qvec_p(var_A3) = + R * v * (UR0 * BZ0 - UZ0 * BR0)              &
                             + tauIC_AP * v * tau_IC*F0/rho0_corr/BB2 * R*Bp0 * BgradPe&
                             - eta_T * v_Z * R                * BR0         &
                             + eta_T * ( 2.d0 * v + R * v_R ) * BZ0         &
                             + R * v * (eta_R * BZ0 - eta_Z * BR0)          &
                             + eta_T * v * current_source_Jp(ms,mt)         &
                             - eta_num * lap_Vstar * lap_A3

            !###################################################################################################
            !#  equation 4   (R component momentum equation)                                                   #
            !###################################################################################################
            Pvec_prev(var_UR) =   v * rho0_corr * delta_g(mp,var_UR, ms,mt)
                      
            Qvec_p(var_UR)    = + Qconv_UR              &
                                + (v_R+v/R) * p0        &
                                + JxB_UR__p             &
                                + visco_T * Qvisc_UR__p &
#ifdef VISCOSITY_TENSOR_FORM
            ! --- the missed term in tensor form viscosity
                                + v * Qgrad_visc_dot_grad_UR_p &
#else
            ! --- the missed term in full MHD viscosity
                                - v * Qgrad_visc_dot_grad_UR_p &
#ifdef VISCOSITY_INSIDE_LAPLACIAN
                                + v * Qvisc_gradmu_dot_grad_vecV_UR_p &
#endif
#ifdef VISCOSITY_COMPLETE_TERMS
                                + v * Qgrad_visc_dot_gradU_T_UR__p          &
                                + v * Q_minus_2_over_3_divU_grad_visc_UR__p &
                                +     Q_1_over_3_visc_grad_divU_UR__p       &
#endif
#endif
                                - v * PneoR             &
                                - v * particle_source(ms,mt) * UR0 &
                                - visco_num * lap_Vstar * lap_UR
            Qvec_k(var_UR)    = + JxB_UR__k             &
                                + visco_T * Qvisc_UR__k

#ifdef ENABLE_VISCO_PAR_AND_PARPAR
            Qvec_p(var_UR)     = Qvec_p(var_UR) + &
                                - visco_par_prof     * Qvec_UR_visco_par__p &
                                - visco_pp_diff_prof * Qvec_UR_visco_parpar__p
            Qvec_k(var_UR)     = Qvec_k(var_UR) + &
                                - visco_par_prof     * Qvec_UR_visco_par__k &
                                - visco_pp_diff_prof * Qvec_UR_visco_parpar__k
#endif

            !###################################################################################################
            !#  equation 5   (Z component momentum equation)                                                   #
            !###################################################################################################
            Pvec_prev(var_UZ) =   v * rho0_corr * delta_g(mp,var_UZ, ms,mt)

            Qvec_p(var_UZ)  = + Qconv_UZ              &
                              + v_Z * p0              &
                              + JxB_UZ__p             &
                              + visco_T * Qvisc_UZ__p &
#ifdef VISCOSITY_TENSOR_FORM
            ! --- the missed term in tensor form viscosity
                              + v * Qgrad_visc_dot_grad_UZ_p &
#else
            ! --- the missed term in full MHD viscosity
                              - v * Qgrad_visc_dot_grad_UZ_p &
#ifdef VISCOSITY_INSIDE_LAPLACIAN
                              + v * Qvisc_gradmu_dot_grad_vecV_UZ_p &
#endif
#ifdef VISCOSITY_COMPLETE_TERMS
                              + v * Qgrad_visc_dot_gradU_T_UZ__p          &
                              + v * Q_minus_2_over_3_divU_grad_visc_UZ__p &
                              +     Q_1_over_3_visc_grad_divU_UZ__p       &
#endif
#endif
                              - v * PneoZ             &
                              - v * particle_source(ms,mt) * UZ0 &
                              - visco_num * lap_Vstar * lap_UZ
            Qvec_k(var_UZ)  = + JxB_UZ__k             &
                              + visco_T * Qvisc_UZ__k

#ifdef ENABLE_VISCO_PAR_AND_PARPAR
            Qvec_p(var_UZ)  = Qvec_p(var_UZ) + &
                              - visco_par_prof     * Qvec_UZ_visco_par__p &
                              - visco_pp_diff_prof * Qvec_UZ_visco_parpar__p
            Qvec_k(var_UZ)  = Qvec_k(var_UZ) + &
                              - visco_par_prof     * Qvec_UZ_visco_par__k &
                              - visco_pp_diff_prof * Qvec_UZ_visco_parpar__k
#endif
                           
            !###################################################################################################
            !#  equation 6 (Phi component momentum equation)                                                   #
            !###################################################################################################
#ifdef MOMENTUM_UPHI_IN_PARALLEL
            Pvec_prev(var_Up) = + v * rho0_corr * BR0 * delta_g(mp,var_UR, ms,mt) &
                                + v * rho0_corr * BZ0 * delta_g(mp,var_UZ, ms,mt) &
                                + v * rho0_corr * BP0 * delta_g(mp,var_Up, ms,mt)

            Qvec_p(var_Up)    = + BR0 * Qconv_UR              &
                                + BZ0 * Qconv_UZ              &
                                + Bp0 * Qconv_Up              &
                                + BR0 * visco_T * Qvisc_UR__p &
                                + BZ0 * visco_T * Qvisc_UZ__p &
                                + Bp0 * visco_T * Qvisc_Up__p &
#ifdef VISCOSITY_TENSOR_FORM
            ! --- the missed term in tensor form viscosity, the B\cdot operation should be considerd when perform integration-by-parts
                                - v * QdivU_BdotGradVisc      &
                                - v * curl_U0_dot_gradmu_x_B0 &
                                - v * visco_T * curl_U0_dot_curl_B0 &
#else
            ! --- the missed term in full MHD viscosity, the B\cdot operation should be considerd when perform integration-by-parts
                                - BR0 * v * Qgrad_visc_dot_grad_UR_p &
                                - BZ0 * v * Qgrad_visc_dot_grad_UZ_p &
                                - Bp0 * v * Qgrad_visc_dot_grad_Up_p &
                                - visco_T * v * Qvisc_grad_BR_dot_grad_UR_p   &
                                - visco_T * v * Qvisc_grad_BZ_dot_grad_UZ_p   &
                                - visco_T * v * Qvisc_grad_Bp_dot_grad_Up_p   &
#ifdef VISCOSITY_INSIDE_LAPLACIAN
                                + BR0 * v * Qvisc_gradmu_dot_grad_vecV_UR_p &
                                + BZ0 * v * Qvisc_gradmu_dot_grad_vecV_UZ_p &
                                + Bp0 * v * Qvisc_gradmu_dot_grad_vecV_Up_p &
#endif
#ifdef VISCOSITY_COMPLETE_TERMS
                                + BR0 * v * Qgrad_visc_dot_gradU_T_UR__p &
                                + BZ0 * v * Qgrad_visc_dot_gradU_T_UZ__p &
                                + Bp0 * v * Qgrad_visc_dot_gradU_T_Up__p &
                                + BR0 * v * Q_minus_2_over_3_divU_grad_visc_UR__p &
                                + BZ0 * v * Q_minus_2_over_3_divU_grad_visc_UZ__p &
                                + Bp0 * v * Q_minus_2_over_3_divU_grad_visc_Up__p &
                                + Q_1_over_3_B_dot_visc_grad_divU__p &
#endif
#endif
                                + p0 * BgradVstar__p          &
                                - v * (BR0*PneoR + BZ0*PneoZ) &
                                - v * particle_source(ms,mt) * (BR0*UR0 + BZ0*UZ0 + Bp0*Up0) &
                                - visco_num * lap_Vstar * (BR0 * lap_UR + BZ0 * lap_UZ + Bp0 * lap_Up)
            Qvec_k(var_Up)    = + BR0 * visco_T * Qvisc_UR__k &
                                + BZ0 * visco_T * Qvisc_UZ__k &
                                + Bp0 * visco_T * Qvisc_Up__k &
#ifndef VISCOSITY_TENSOR_FORM
#ifdef VISCOSITY_COMPLETE_TERMS
                                + Q_1_over_3_B_dot_visc_grad_divU__k &
#endif
#endif
                                + p0 * BgradVstar__k

#ifdef ENABLE_VISCO_PAR_AND_PARPAR
            Qvec_p(var_Up)    = Qvec_p(var_Up) + &
                                - visco_par_prof     * (BR0 * Qvec_UR_visco_par__p    + BZ0 * Qvec_UZ_visco_par__p    + Bp0 * Qvec_Up_visco_par__p   ) &
                                - visco_pp_diff_prof * (BR0 * Qvec_UR_visco_parpar__p + BZ0 * Qvec_UZ_visco_parpar__p + Bp0 * Qvec_Up_visco_parpar__p) &
                                - visco_par_prof     * (BR0 * Qvec_UR_visco_par_for_Uphi__p    + BZ0 * Qvec_UZ_visco_par_for_Uphi__p    + Bp0 * Qvec_Up_visco_par_for_Uphi__p) &
                                - visco_pp_diff_prof * (BR0 * Qvec_UR_visco_parpar_for_Uphi__p + BZ0 * Qvec_UZ_visco_parpar_for_Uphi__p + Bp0 * Qvec_Up_visco_parpar_for_Uphi__p)
            Qvec_k(var_Up)    = Qvec_k(var_Up) + &
                                - visco_par_prof     * (BR0 * Qvec_UR_visco_par__k    + BZ0 * Qvec_UZ_visco_par__k    + Bp0 * Qvec_Up_visco_par__k  ) &
                                - visco_pp_diff_prof * (BR0 * Qvec_UR_visco_parpar__k + BZ0 * Qvec_UZ_visco_parpar__k + Bp0 * Qvec_Up_visco_parpar__k)
#endif

#else
            Pvec_prev(var_Up) = v * rho0_corr * delta_g(mp,var_Up, ms,mt)

            Qvec_p(var_Up)    = + Qconv_Up                    &
                                + JxB_Up__p                   &
                                + visco_T * Qvisc_Up__p       &
#ifdef VISCOSITY_TENSOR_FORM
            ! --- the missed term in tensor form viscosity
                                + v * Qgrad_visc_dot_grad_Up_p &
#else
                                - v * Qgrad_visc_dot_grad_Up_p &
#ifdef VISCOSITY_INSIDE_LAPLACIAN
                                + v * Qvisc_gradmu_dot_grad_vecV_Up_p &
#endif
#ifdef VISCOSITY_COMPLETE_TERMS
                                + v * Qgrad_visc_dot_gradU_T_Up__p          &
                                + v * Q_minus_2_over_3_divU_grad_visc_Up__p &
                                +     Q_1_over_3_visc_grad_divU_Up__p       &
#endif
#endif
                                - v * particle_source(ms,mt) * Up0 &
                                - visco_num * lap_Vstar * lap_Up
            Qvec_k(var_Up)    = + v_p__k * p0 / R             &
                                + JxB_Up__k                   &
#ifndef VISCOSITY_TENSOR_FORM
#ifdef VISCOSITY_COMPLETE_TERMS
                                +     Q_1_over_3_visc_grad_divU_Up__k       &
#endif
#endif
                                + visco_T * Qvisc_Up__k

#ifdef ENABLE_VISCO_PAR_AND_PARPAR
            Qvec_p(var_Up)    = Qvec_p(var_Up) + &
                                - visco_par_prof     * (Qvec_Up_visco_par__p   ) &
                                - visco_pp_diff_prof * (Qvec_Up_visco_parpar__p) 
            Qvec_k(var_Up)    = Qvec_k(var_Up) + &
                                - visco_par_prof     * (Qvec_Up_visco_par__k  ) &
                                - visco_pp_diff_prof * (Qvec_Up_visco_parpar__k)
#endif
#endif
            
            !###################################################################################################
            !#  equation 7 (Density equation)                                                                  #
            !###################################################################################################
            Pvec_prev(var_rho) =   v * delta_g(mp,var_rho,ms,mt)

            Qvec_p(var_rho) = - v * ( rho0 * divU + UgradRho )                  &
                              + tauIC_rho * rho0 * VdiaGradVstar__p             &
                              - D_prof * gradRho_gradVstar__p                   &
                              - (D_par_prof-D_prof) * BgradVstar__p * BgradRho / BB2 &
                              + v * particle_source(ms,mt)                      &
                              - D_perp_num * lap_Vstar * lap_Rho
            Qvec_k(var_rho) = + tauIC_rho * rho0 * VdiaGradVstar__k             &
                              - D_prof * gradRho_gradVstar__k                   &
                              - (D_par_prof-D_prof) * BgradVstar__k * BgradRho / BB2

            !###################################################################################################
            !#  equation 8 (Ion Pressure equation)                                                             #
            !###################################################################################################
            Pvec_prev(var_Ti) =   v * rho0_corr * delta_g(mp,var_Ti, ms,mt) &
                                + v * Ti0_corr  * delta_g(mp,var_rho,ms,mt)

            Qvec_p(var_Ti) = + v * ( - rho0 * UgradTi -  Ti0 * UgradRho - gamma * pi0 * divU ) &
                             + v * heat_source_i(ms,mt)                                        &
                             + v * (gamma-1.d0) * Qvisc_T                                      &
                             - ZKi_prof * gradTi_gradVstar__p                                  &
                             - (ZKi_par_T-ZKi_prof) * BgradVstar__p * BgradTi / BB2            &
                             + v * dTi_e                                                       &
                             - ZK_i_perp_num * lap_Vstar * lap_Ti
            Qvec_k(var_Ti) = - ZKi_prof * gradTi_gradVstar__k                                  &
                             - (ZKi_par_T-ZKi_prof) * BgradVstar__k * BgradTi / BB2
#ifdef ENABLE_DIA_DRIFT_TEMPERATURE
            Qvec_p(var_Ti) = Qvec_p(var_Ti)                                                    &
                             + tauIC_tem * (gamma-1.d0) * v                                    & 
                             * (Ti0 * (VdiaR0 * rho0_R + VdiaZ0 * rho0_Z + VdiaP0 * rho0_p / R)&
                             + rho0 * (VdiaR0 * Ti0_R + VdiaZ0 * Ti0_Z + VdiaP0 * Ti0_p / R))  &
                             + tauIC_tem * gamma * rho0 * Ti0 * VdiaGradVstar__p
            Qvec_k(var_Ti) = Qvec_k(var_Ti)                                                    &
                             + tauIC_tem * gamma * rho0 * Ti0 * VdiaGradVstar__k
#endif


            !###################################################################################################
            !#  equation 9 (Electron Pressure equation)                                                        #
            !###################################################################################################
            Pvec_prev(var_Te) =   v * rho0_corr * delta_g(mp,var_Te, ms,mt) &
                                + v * Te0_corr  * delta_g(mp,var_rho,ms,mt)

            Qvec_p(var_Te) = + v * ( - rho0 * UgradTe - Te0 * UgradRho -  gamma * pe0 * divU ) &
                             + v * heat_source_e(ms,mt)                                        &
                             + v * (gamma-1.d0) * Qvisc_T                                      &
                             - ZKe_prof * gradTe_gradVstar__p                                  &
                             - (ZKe_par_T-ZKe_prof) * BgradVstar__p * BgradTe / BB2            &
                             + v * dTe_i                                                       &
                             + v * (gamma-1.0d0) * eta_T_ohm * JJ2                             &
                             - ZK_e_perp_num * lap_Vstar * lap_Te
            Qvec_k(var_Te) = - ZKe_prof * gradTe_gradVstar__k                                  &
                             - (ZKe_par_T-ZKe_prof) * BgradVstar__k * BgradTe / BB2
#ifdef ENABLE_DIA_DRIFT_TEMPERATURE
            Qvec_p(var_Te) = Qvec_p(var_Te)                                                    &
                             + tauIC_tem * (gamma-1.d0) * v                                    & 
                             * (Te0 * (VdiaR0 * rho0_R + VdiaZ0 * rho0_Z + VdiaP0 * rho0_p / R)&
                             + rho0 * (VdiaR0 * Te0_R + VdiaZ0 * Te0_Z + VdiaP0 * Te0_p / R))  &
                             + tauIC_tem * gamma * rho0 * Te0 * VdiaGradVstar__p
            Qvec_k(var_Te) = Qvec_k(var_Te)                                                    &
                             + tauIC_tem * gamma * rho0 * Te0 * VdiaGradVstar__k
#endif

            if (use_vms) call add_vms_to_rhs()
            
            ! --- Fill Up the RHS
            if (use_fft) then
              index_ij =       n_var*n_degrees*(i-1) +       n_var*(j-1) + 1
              do ivar= 1,n_var
                RHS_p_ij(ivar) = tstep * Qvec_p(ivar) + zeta * Pvec_prev(ivar) * tstep / tstep_prev
                RHS_k_ij(ivar) = tstep * Qvec_k(ivar)

                ij = index_ij + (ivar-1)
                RHS_p(mp,ij)   =  RHS_p(mp,ij) + RHS_p_ij(ivar) * wst * R * xjac
                RHS_k(mp,ij)   =  RHS_k(mp,ij) + RHS_k_ij(ivar) * wst * R * xjac
              enddo
            else
              index_ij = n_tor_local*n_var*n_degrees*(i-1) + n_tor_local*n_var*(j-1) + im - n_tor_start +1
              do ivar= 1,n_var
                RHS_p_ij(ivar) = tstep * Qvec_p(ivar) + zeta * Pvec_prev(ivar) * tstep / tstep_prev
                RHS_k_ij(ivar) = tstep * Qvec_k(ivar)

                ij = index_ij + (ivar-1)*n_tor_local
                RHS(ij)        =  RHS(ij) + (RHS_p_ij(ivar) + RHS_k_ij(ivar)) * wst * R * xjac
              enddo
            endif

            do l=1,n_degrees

              do k=1,n_vertex_max

                do in =  n_tor_start, n_tor_end

                  ! --- Basis functions
                  bf    = H(k,l,ms,mt)   * element%size(k,l) * HHZ(in,mp)
                  bf_p  = H(k,l,ms,mt)   * element%size(k,l) * HHZ_p(in,mp)
                  bf_s  = H_s(k,l,ms,mt) * element%size(k,l) * HHZ(in,mp)
                  bf_t  = H_t(k,l,ms,mt) * element%size(k,l) * HHZ(in,mp)
                  bf_R = (   y_t(ms,mt) * bf_s - y_s(ms,mt) * bf_t ) / xjac
                  bf_Z = ( - x_t(ms,mt) * bf_s + x_s(ms,mt) * bf_t ) / xjac
                  bf_ss = H_ss(k,l,ms,mt) * element%size(k,l) * HHZ(in,mp)
                  bf_tt = H_tt(k,l,ms,mt) * element%size(k,l) * HHZ(in,mp)
                  bf_st = H_st(k,l,ms,mt) * element%size(k,l) * HHZ(in,mp)
                  bf_RR = (bf_ss * y_t(ms,mt)**2 - 2.d0*bf_st * y_s(ms,mt)*y_t(ms,mt) + bf_tt * y_s(ms,mt)**2  &
                         + bf_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                              &
                         + bf_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )    / xjac**2               & 
                         - xjac_R * (bf_s * y_t(ms,mt) - bf_t * y_s(ms,mt)) / xjac**2
                  bf_ZZ = (bf_ss * x_t(ms,mt)**2 - 2.d0*bf_st * x_s(ms,mt)*x_t(ms,mt) + bf_tt * x_s(ms,mt)**2  &
                         + bf_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                              &
                         + bf_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )       / xjac**2            & 
                         - xjac_Z * (- bf_s * x_t(ms,mt) + bf_t * x_s(ms,mt) ) / xjac**2

                  bf_pp  = H(k,l,ms,mt)   * element%size(k,l) * HHZ_pp(in,mp)
                  bf_sp  = H_s(k,l,ms,mt) * element%size(k,l) * HHZ_p(in,mp)
                  bf_tp  = H_t(k,l,ms,mt) * element%size(k,l) * HHZ_p(in,mp)
                  bf_Rp = (   y_t(ms,mt) * bf_sp - y_s(ms,mt) * bf_tp ) / xjac
                  bf_Zp = ( - x_t(ms,mt) * bf_sp + x_s(ms,mt) * bf_tp ) / xjac
                  bf_RZ = (- bf_ss * y_t(ms,mt)*x_t(ms,mt) - bf_tt * x_s(ms,mt)*y_s(ms,mt) &
                         + bf_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  ) &
                         - bf_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) ) &
                         - bf_t * (x_st(ms,mt)*y_s(ms,mt)  - x_ss(ms,mt)*y_t(ms,mt) ) ) / xjac**2             &
                         - xjac_R * (- bf_s * x_t(ms,mt) + bf_t * x_s(ms,mt) ) / xjac**2

                  UR    = bf    ;  UZ    = bf    ;  Up    = bf
                  UR_R  = bf_R  ;  UZ_R  = bf_R  ;  Up_R  = bf_R
                  UR_Z  = bf_Z  ;  UZ_Z  = bf_Z  ;  Up_Z  = bf_Z
                  UR_p__n = bf_p;  UZ_p__n = bf_p;  Up_p__n = bf_p
                  UR_s  = bf_s  ;  UZ_s  = bf_s  ;  Up_s  = bf_s
                  UR_t  = bf_t  ;  UZ_t  = bf_t  ;  Up_t  = bf_t

                  AR    = bf    ;  AZ    = bf    ;  A3    = bf    ; Ti    = bf    ; Te    = bf    ; rho    = bf
                  AR_R  = bf_R  ;  AZ_R  = bf_R  ;  A3_R  = bf_R  ; Ti_R  = bf_R  ; Te_R  = bf_R  ; rho_R  = bf_R 
                  AR_Z  = bf_Z  ;  AZ_Z  = bf_Z  ;  A3_Z  = bf_Z  ; Ti_Z  = bf_Z  ; Te_Z  = bf_Z  ; rho_Z  = bf_Z
                  AR_p__n = bf_p;  AZ_p__n = bf_p;  A3_p__n = bf_p; Ti_p__n = bf_p; Te_p__n = bf_p; rho_p__n = bf_p
                  AR_s  = bf_s  ;  AZ_s  = bf_s  ;  A3_s  = bf_s  ; Ti_s  = bf_s  ; Te_s  = bf_s  ; rho_s  = bf_s 
                  AR_t  = bf_t  ;  AZ_t  = bf_t  ;  A3_t  = bf_t  ; Ti_t  = bf_T  ; Te_t  = bf_T  ; rho_t  = bf_T 

                  AR_RR = bf_RR ; AR_RZ = bf_RZ ; AR_ZZ = bf_ZZ
                  AZ_RR = bf_RR ; AZ_RZ = bf_RZ ; AZ_ZZ = bf_ZZ
                  A3_RR = bf_RR ; A3_RZ = bf_RZ ; A3_ZZ = bf_ZZ

                  AR_Rp = bf_Rp ; AR_Zp  = bf_Zp ; AR_pp = bf_pp
                  AZ_Rp = bf_Rp ; AZ_Zp  = bf_Zp ; AZ_pp = bf_pp
                  A3_Rp = bf_Rp ; A3_Zp  = bf_Zp ; A3_pp = bf_pp

                  UR_RR = bf_RR ; UR_ZZ  = bf_ZZ ; UR_pp = bf_pp
                  UZ_RR = bf_RR ; UZ_ZZ  = bf_ZZ ; UZ_pp = bf_pp
                  Up_RR = bf_RR ; Up_ZZ  = bf_ZZ ; Up_pp = bf_pp

                  rho_RR= bf_RR ; rho_ZZ = bf_ZZ ; rho_pp= bf_pp
                  Ti_RR = bf_RR ; Ti_ZZ  = bf_ZZ ; Ti_pp = bf_pp
                  Te_RR = bf_RR ; Te_ZZ  = bf_ZZ ; Te_pp = bf_pp

                  ! --- Linearised quantities
                  BR0_AR    =   0.d0        ; BR0_AZ__n = - AZ_p__n / R ; BR0_A3 =   A3_Z / R
                  BZ0_AR__n =   AR_p__n / R ; BZ0_AZ    =   0.d0        ; BZ0_A3 = - A3_R / R
                  Bp0_AR    = - AR_Z        ; Bp0_AZ    =   AZ_R        ; Bp0_A3 =   0.d0

                  BB2_AR__p = 2.d0*(BR0_AR    * BR0                   + Bp0_AR * Bp0 )
                  BB2_AR__n = 2.d0*(                + BZ0_AR__n * BZ0                )
                  BB2_AZ__p = 2.d0*(                + BZ0_AZ    * BZ0 + Bp0_AZ * Bp0 )
                  BB2_AZ__n = 2.d0*(BR0_AZ__n * BR0                                  )
                  BB2_A3    = 2.d0*(BR0_A3    * BR0 + BZ0_A3    * BZ0 + Bp0_A3 * Bp0 )

                  BgradTi_AR__p = BR0_AR    * Ti0_R                     + Bp0_AR * Ti0_p / R
                  BgradTi_AR__n =                   + BZ0_AR__n * Ti0_Z
                  BgradTi_AZ__p =                   + BZ0_AZ    * Ti0_Z + Bp0_AZ * Ti0_p / R
                  BgradTi_AZ__n = BR0_AZ__n * Ti0_R
                  BgradTi_A3    = BR0_A3    * Ti0_R + BZ0_A3    * Ti0_Z + Bp0_A3 * Ti0_p / R
                  BgradTi_Ti__p = BR0 * Ti_R + BZ0 * Ti_Z
                  BgradTi_Ti__n = Bp0 * Ti_p__n / R

                  BgradTe_AR__p = BR0_AR    * Te0_R                     + Bp0_AR * Te0_p / R
                  BgradTe_AR__n =                   + BZ0_AR__n * Te0_Z
                  BgradTe_AZ__p =                   + BZ0_AZ    * Te0_Z + Bp0_AZ * Te0_p / R
                  BgradTe_AZ__n = BR0_AZ__n * Te0_R
                  BgradTe_A3    = BR0_A3    * Te0_R + BZ0_A3    * Te0_Z + Bp0_A3 * Te0_p / R
                  BgradTe_Te__p = BR0 * Te_R + BZ0 * Te_Z
                  BgradTe_Te__n = Bp0 * Te_p__n / R

                  BgradRho_AR__p  = BR0_AR    * rho0_R                      + Bp0_AR * rho0_p / R
                  BgradRho_AR__n  =                    + BZ0_AR__n * rho0_Z
                  BgradRho_AZ__p  =                    + BZ0_AZ    * rho0_Z + Bp0_AZ * rho0_p / R
                  BgradRho_AZ__n  = BR0_AZ__n * rho0_R
                  BgradRho_A3     = BR0_A3    * rho0_R + BZ0_A3    * rho0_Z + Bp0_A3 * rho0_p / R
                  BgradRho_rho__p = BR0 * rho_R + BZ0 * rho_Z
                  BgradRho_rho__n = Bp0 * rho_p__n / R

                  BgradPe_AR__p  = rho0*BgradTe_AR__p + Te0*BgradRho_AR__p
                  BgradPe_AR__n  = rho0*BgradTe_AR__n + Te0*BgradRho_AR__n
                  BgradPe_AZ__p  = rho0*BgradTe_AZ__p + Te0*BgradRho_AZ__p
                  BgradPe_AZ__n  = rho0*BgradTe_AZ__n + Te0*BgradRho_AZ__n
                  BgradPe_A3     = rho0*BgradTe_A3    + Te0*BgradRho_A3   
                  BgradPe_rho__p = rho*BgradTe        + Te0*BgradRho_rho__p
                  BgradPe_rho__n =                    + Te0*BgradRho_rho__n
                  BgradPe_Te__p  = rho0*BgradTe_Te__p + Te *BgradRho
                  BgradPe_Te__n  = rho0*BgradTe_Te__n

                  UgradUR_UR__p = UR  * UR0_R + UR0 * UR_R + UZ0 * UR_Z
                  UgradUR_UR__n = Up0 * UR_p__n  / R
                  UgradUR_UZ    = UZ  * UR0_Z
                  UgradUR_Up    = Up  * UR0_p / R

                  UgradUZ_UR    = UR  * UZ0_R
                  UgradUZ_UZ__p = UR0 * UZ_R + UZ * UZ0_Z + UZ0 * UZ_Z
                  UgradUZ_UZ__n = Up0 * UZ_p__n  / R
                  UgradUZ_Up    = Up  * UZ0_p / R

                  UgradUp_UR    = UR  * Up0_R
                  UgradUp_UZ    = UZ  * Up0_Z
                  UgradUp_Up__p = UR0 * Up_R + UZ0 * Up_Z + Up  * Up0_p / R
                  UgradUp_Up__n = Up0 * Up_p__n  / R

                  UgradTi_UR    = UR * Ti0_R
                  UgradTi_UZ    = UZ * Ti0_Z
                  UgradTi_Up    = Up * Ti0_p / R
                  UgradTi_Ti__p = UR0 * Ti_R + UZ0 * Ti_Z
                  UgradTi_Ti__n = Up0 * Ti_p__n / R

                  UgradTe_UR    = UR * Te0_R
                  UgradTe_UZ    = UZ * Te0_Z
                  UgradTe_Up    = Up * Te0_p / R
                  UgradTe_Te__p = UR0 * Te_R + UZ0 * Te_Z
                  UgradTe_Te__n = Up0 * Te_p__n / R

                  UgradRho_UR     = UR * rho0_R
                  UgradRho_UZ     = UZ * rho0_Z
                  UgradRho_Up     = Up * rho0_p / R
                  UgradRho_rho__p = UR0 * rho_R + UZ0 * rho_Z
                  UgradRho_rho__n = Up0 * rho_p__n / R

                  divU_UR    = UR_R + UR / R
                  divU_UZ    = UZ_Z
                  divU_Up__n = Up_p__n / R

                  divRhoU_UR    = rho0 * divU_UR + UgradRho_UR
                  divRhoU_UZ    = rho0 * divU_UZ + UgradRho_UZ
                  divRhoU_Up__p = UgradRho_Up
                  divRhoU_Up__n = rho0 * divU_Up__n
                  divRhoU_rho__p = rho * divU + UgradRho_rho__p
                  divRhoU_rho__n = UgradRho_rho__n

                  VdiaR0_AR__p  = + tau_IC*F0 / (R * rho0_corr * BB2   ) * (                  - R*Bp0_AR*pi0_Z) &
                                  - tau_IC*F0 / (R * rho0_corr * BB2**2) * (  BZ0      *pi0_p - R*Bp0   *pi0_Z) * BB2_AR__p
                  VdiaR0_AR__n  = + tau_IC*F0 / (R * rho0_corr * BB2   ) * (  BZ0_AR__n*pi0_p                 ) &
                                  - tau_IC*F0 / (R * rho0_corr * BB2**2) * (  BZ0      *pi0_p - R*Bp0   *pi0_Z) * BB2_AR__n
                  VdiaR0_AZ__p  = + tau_IC*F0 / (R * rho0_corr * BB2   ) * (  BZ0_AZ   *pi0_p - R*Bp0_AZ*pi0_Z) &
                                  - tau_IC*F0 / (R * rho0_corr * BB2**2) * (  BZ0      *pi0_p - R*Bp0   *pi0_Z) * BB2_AZ__p
                  VdiaR0_AZ__n  = - tau_IC*F0 / (R * rho0_corr * BB2**2) * (  BZ0      *pi0_p - R*Bp0   *pi0_Z) * BB2_AZ__n
                  VdiaR0_A3__p  = + tau_IC*F0 / (R * rho0_corr * BB2   ) * (  BZ0_A3   *pi0_p - R*Bp0_A3*pi0_Z) &
                                  - tau_IC*F0 / (R * rho0_corr * BB2**2) * (  BZ0      *pi0_p - R*Bp0   *pi0_Z) * BB2_A3
                  VdiaR0_A3__n  = 0.d0
                  VdiaR0_rho__p = + tau_IC*F0 / (R * rho0_corr * BB2   ) * (  BZ0*(rho*Ti0_p) - R*Bp0*(rho*Ti0_Z+rho_Z*Ti0)) &
                                  - tau_IC*F0 / (R * rho0_corr**2 * BB2) * (  BZ0*pi0_p       - R*Bp0*pi0_Z                ) * rho
                  VdiaR0_rho__n = + tau_IC*F0 / (R * rho0_corr * BB2   ) * (  BZ0*(rho_p__n*Ti0)                              )
                  VdiaR0_Ti__p  = + tau_IC*F0 / (R * rho0_corr * BB2   ) * (  BZ0*(rho0_p*Ti) - R*Bp0*(rho0*Ti_Z+rho0_Z*Ti))
                  VdiaR0_Ti__n  = + tau_IC*F0 / (R * rho0_corr * BB2   ) * (  BZ0*(rho0*Ti_p__n)                              )

                  VdiaZ0_AR__p  = + tau_IC*F0 / (R * rho0_corr * BB2   ) * (R*BP0_AR*pi0_R -   BR0_AR   *pi0_p) &
                                  - tau_IC*F0 / (R * rho0_corr * BB2**2) * (R*BP0   *pi0_R -   BR0      *pi0_p) * BB2_AR__p
                  VdiaZ0_AR__n  = - tau_IC*F0 / (R * rho0_corr * BB2**2) * (R*BP0   *pi0_R -   BR0      *pi0_p) * BB2_AR__n
                  VdiaZ0_AZ__p  = + tau_IC*F0 / (R * rho0_corr * BB2   ) * (R*BP0_AZ*pi0_R                    ) &
                                  - tau_IC*F0 / (R * rho0_corr * BB2**2) * (R*BP0   *pi0_R -   BR0      *pi0_p) * BB2_AZ__p
                  VdiaZ0_AZ__n  = + tau_IC*F0 / (R * rho0_corr * BB2   ) * (               -   BR0_AZ__n*pi0_p) &
                                  - tau_IC*F0 / (R * rho0_corr * BB2**2) * (R*BP0   *pi0_R -   BR0      *pi0_p) * BB2_AZ__n
                  VdiaZ0_A3__p  = + tau_IC*F0 / (R * rho0_corr * BB2   ) * (R*BP0_A3*pi0_R -   BR0_A3   *pi0_p) &
                                  - tau_IC*F0 / (R * rho0_corr * BB2**2) * (R*BP0   *pi0_R -   BR0      *pi0_p) * BB2_A3
                  VdiaZ0_A3__n  = 0.d0
                  VdiaZ0_rho__p = + tau_IC*F0 / (R * rho0_corr * BB2   ) * (R*BP0*(rho*Ti0_R+rho_R*Ti0) -   BR0*(rho*Ti0_p)) &
                                  - tau_IC*F0 / (R * rho0_corr**2 * BB2) * (R*BP0*pi0_R                 -   BR0*pi0_p      ) * rho
                  VdiaZ0_rho__n = + tau_IC*F0 / (R * rho0_corr * BB2   ) * (                            -   BR0*(rho_p__n*Ti0))
                  VdiaZ0_Ti__p  = + tau_IC*F0 / (R * rho0_corr * BB2   ) * (R*BP0*(Ti*rho0_R+Ti_R*rho0) -   BR0*(Ti*rho0_p))
                  VdiaZ0_Ti__n  = + tau_IC*F0 / (R * rho0_corr * BB2   ) * (                            -   BR0*(Ti_p__n*rho0))

                  VdiaP0_AR__p  = + tau_IC*F0 / (    rho0_corr * BB2   ) * (  BR0_AR   *pi0_Z                    ) &
                                  - tau_IC*F0 / (    rho0_corr * BB2**2) * (  BR0      *pi0_Z -   BZ0      *pi0_R) * BB2_AR__p
                  VdiaP0_AR__n  = + tau_IC*F0 / (    rho0_corr * BB2   ) * (                  -   BZ0_AR__n*pi0_R) &
                                  - tau_IC*F0 / (    rho0_corr * BB2**2) * (  BR0      *pi0_Z -   BZ0      *pi0_R) * BB2_AR__n
                  VdiaP0_AZ__p  = + tau_IC*F0 / (    rho0_corr * BB2   ) * (                  -   BZ0_AZ   *pi0_R) &
                                  - tau_IC*F0 / (    rho0_corr * BB2**2) * (  BR0      *pi0_Z -   BZ0      *pi0_R) * BB2_AZ__p
                  VdiaP0_AZ__n  = + tau_IC*F0 / (    rho0_corr * BB2   ) * (  BR0_AZ__n*pi0_Z                    ) &
                                  - tau_IC*F0 / (    rho0_corr * BB2**2) * (  BR0      *pi0_Z -   BZ0      *pi0_R) * BB2_AZ__n
                  VdiaP0_A3__p  = + tau_IC*F0 / (    rho0_corr * BB2   ) * (  BR0_A3   *pi0_Z -   BZ0_A3   *pi0_R) &
                                  - tau_IC*F0 / (    rho0_corr * BB2**2) * (  BR0      *pi0_Z -   BZ0      *pi0_R) * BB2_A3
                  VdiaP0_A3__n  = 0.d0
                  VdiaP0_rho__p = + tau_IC*F0 / (    rho0_corr * BB2   ) * (  BR0*(rho*Ti0_Z+rho_Z*Ti0) -   BZ0*(rho*Ti0_R+rho_R*Ti0)) &
                                  - tau_IC*F0 / (    rho0_corr**2 * BB2) * (  BR0*pi0_Z                 -   BZ0*pi0_R                ) * rho
                  VdiaP0_rho__n = 0.d0
                  VdiaP0_Ti__p  = + tau_IC*F0 / (    rho0_corr * BB2   ) * (  BR0*(rho0*Ti_Z+rho0_Z*Ti) -   BZ0*(rho0*Ti_R+rho0_R*Ti))
                  VdiaP0_Ti__n  = 0.d0

                  VdiaGradUR_AR__p  = VdiaR0_AR__p  * UR0_R + VdiaZ0_AR__p  * UR0_Z + VdiaP0_AR__p  * UR0_p / R
                  VdiaGradUR_AR__n  = VdiaR0_AR__n  * UR0_R + VdiaZ0_AR__n  * UR0_Z + VdiaP0_AR__n  * UR0_p / R
                  VdiaGradUR_AZ__p  = VdiaR0_AZ__p  * UR0_R + VdiaZ0_AZ__p  * UR0_Z + VdiaP0_AZ__p  * UR0_p / R
                  VdiaGradUR_AZ__n  = VdiaR0_AZ__n  * UR0_R + VdiaZ0_AZ__n  * UR0_Z + VdiaP0_AZ__n  * UR0_p / R
                  VdiaGradUR_A3__p  = VdiaR0_A3__p  * UR0_R + VdiaZ0_A3__p  * UR0_Z + VdiaP0_A3__p  * UR0_p / R
                  VdiaGradUR_A3__n  = VdiaR0_A3__n  * UR0_R + VdiaZ0_A3__n  * UR0_Z + VdiaP0_A3__n  * UR0_p / R
                  VdiaGradUR_rho__p = VdiaR0_rho__p * UR0_R + VdiaZ0_rho__p * UR0_Z + VdiaP0_rho__p * UR0_p / R
                  VdiaGradUR_rho__n = VdiaR0_rho__n * UR0_R + VdiaZ0_rho__n * UR0_Z + VdiaP0_rho__n * UR0_p / R
                  VdiaGradUR_Ti__p  = VdiaR0_Ti__p  * UR0_R + VdiaZ0_Ti__p  * UR0_Z + VdiaP0_Ti__p  * UR0_p / R
                  VdiaGradUR_Ti__n  = VdiaR0_Ti__n  * UR0_R + VdiaZ0_Ti__n  * UR0_Z + VdiaP0_Ti__n  * UR0_p / R
                  VdiaGradUR_UR__p  = VdiaR0        * UR_R  + VdiaZ0        * UR_Z
                  VdiaGradUR_UR__n  =                                               + VdiaP0        * UR_p__n  / R

                  VdiaGradUZ_AR__p  = VdiaR0_AR__p  * UZ0_R + VdiaZ0_AR__p  * UZ0_Z + VdiaP0_AR__p  * UZ0_p / R
                  VdiaGradUZ_AR__n  = VdiaR0_AR__n  * UZ0_R + VdiaZ0_AR__n  * UZ0_Z + VdiaP0_AR__n  * UZ0_p / R
                  VdiaGradUZ_AZ__p  = VdiaR0_AZ__p  * UZ0_R + VdiaZ0_AZ__p  * UZ0_Z + VdiaP0_AZ__p  * UZ0_p / R
                  VdiaGradUZ_AZ__n  = VdiaR0_AZ__n  * UZ0_R + VdiaZ0_AZ__n  * UZ0_Z + VdiaP0_AZ__n  * UZ0_p / R
                  VdiaGradUZ_A3__p  = VdiaR0_A3__p  * UZ0_R + VdiaZ0_A3__p  * UZ0_Z + VdiaP0_A3__p  * UZ0_p / R
                  VdiaGradUZ_A3__n  = VdiaR0_A3__n  * UZ0_R + VdiaZ0_A3__n  * UZ0_Z + VdiaP0_A3__n  * UZ0_p / R
                  VdiaGradUZ_rho__p = VdiaR0_rho__p * UZ0_R + VdiaZ0_rho__p * UZ0_Z + VdiaP0_rho__p * UZ0_p / R
                  VdiaGradUZ_rho__n = VdiaR0_rho__n * UZ0_R + VdiaZ0_rho__n * UZ0_Z + VdiaP0_rho__n * UZ0_p / R
                  VdiaGradUZ_Ti__p  = VdiaR0_Ti__p  * UZ0_R + VdiaZ0_Ti__p  * UZ0_Z + VdiaP0_Ti__p  * UZ0_p / R
                  VdiaGradUZ_Ti__n  = VdiaR0_Ti__n  * UZ0_R + VdiaZ0_Ti__n  * UZ0_Z + VdiaP0_Ti__n  * UZ0_p / R
                  VdiaGradUZ_UZ__p  = VdiaR0        * UZ_R  + VdiaZ0        * UZ_Z
                  VdiaGradUZ_UZ__n  =                                               + VdiaP0        * UZ_p__n  / R

                  VdiaGradUp_AR__p  = VdiaR0_AR__p  * Up0_R + VdiaZ0_AR__p  * Up0_Z + VdiaP0_AR__p  * Up0_p / R
                  VdiaGradUp_AR__n  = VdiaR0_AR__n  * Up0_R + VdiaZ0_AR__n  * Up0_Z + VdiaP0_AR__n  * Up0_p / R
                  VdiaGradUp_AZ__p  = VdiaR0_AZ__p  * Up0_R + VdiaZ0_AZ__p  * Up0_Z + VdiaP0_AZ__p  * Up0_p / R
                  VdiaGradUp_AZ__n  = VdiaR0_AZ__n  * Up0_R + VdiaZ0_AZ__n  * Up0_Z + VdiaP0_AZ__n  * Up0_p / R
                  VdiaGradUp_A3__p  = VdiaR0_A3__p  * Up0_R + VdiaZ0_A3__p  * Up0_Z + VdiaP0_A3__p  * Up0_p / R
                  VdiaGradUp_A3__n  = VdiaR0_A3__n  * Up0_R + VdiaZ0_A3__n  * Up0_Z + VdiaP0_A3__n  * Up0_p / R
                  VdiaGradUp_rho__p = VdiaR0_rho__p * Up0_R + VdiaZ0_rho__p * Up0_Z + VdiaP0_rho__p * Up0_p / R
                  VdiaGradUp_rho__n = VdiaR0_rho__n * Up0_R + VdiaZ0_rho__n * Up0_Z + VdiaP0_rho__n * Up0_p / R
                  VdiaGradUp_Ti__p  = VdiaR0_Ti__p  * Up0_R + VdiaZ0_Ti__p  * Up0_Z + VdiaP0_Ti__p  * Up0_p / R
                  VdiaGradUp_Ti__n  = VdiaR0_Ti__n  * Up0_R + VdiaZ0_Ti__n  * Up0_Z + VdiaP0_Ti__n  * Up0_p / R
                  VdiaGradUp_Up__p  = VdiaR0        * Up_R  + VdiaZ0        * Up_Z
                  VdiaGradUp_Up__n  =                                               + VdiaP0        * Up_p__n  / R

                  !if (NEO) then ! can't remember, I think only the RHS is needed because it's small. maybe all this can be removed...
                  if (.false.) then
                    if (Btht.le. 1.d-10) then
                      Btht_AR__p = 0.d0
                      Btht_AR__n = 0.d0
                      Btht_AZ__p = 0.d0
                      Btht_AZ__n = 0.d0
                      Btht_A3    = 0.d0
                    else
                      Btht_AR__p = 0.d0
                      Btht_AR__n = 0.5*(BR0**2 + BZ0**2)**(-0.5) * 2.0*(BZ0*BZ0_AR__n)
                      Btht_AZ__p = 0.d0
                      Btht_AZ__n = 0.5*(BR0**2 + BZ0**2)**(-0.5) * 2.0*(BR0*BR0_AZ__n)
                      Btht_A3    = 0.5*(BR0**2 + BZ0**2)**(-0.5) * 2.0*(BR0*BR0_A3 + BZ0*BZ0_A3)
                    endif
                    
                    Vtht_AR__p  = + ( BR0      *(VdiaR0_AR__p) + BZ0      *(VdiaZ0_AR__p) ) / Btht &
                                  - ( BR0      *(UR0+VdiaR0)   + BZ0      *(UZ0+VdiaZ0)   ) / Btht**2 * Btht_AR__p
                    Vtht_AR__n  = + ( BR0      *(VdiaR0_AR__n) + BZ0      *(VdiaZ0_AR__n) ) / Btht &
                                  + (                          + BZ0_AR__n*(UZ0+VdiaZ0)   ) / Btht &
                                  - ( BR0      *(UR0+VdiaR0)   + BZ0      *(UZ0+VdiaZ0)   ) / Btht**2 * Btht_AR__n
                    Vtht_AZ__p  = + ( BR0      *(VdiaR0_AZ__p) + BZ0      *(VdiaZ0_AZ__p) ) / Btht &
                                  - ( BR0      *(UR0+VdiaR0)   + BZ0      *(UZ0+VdiaZ0)   ) / Btht**2 * Btht_AZ__p
                    Vtht_AZ__n  = + ( BR0      *(VdiaR0_AZ__n) + BZ0      *(VdiaZ0_AZ__n) ) / Btht &
                                  + ( BR0_AZ__n*(UR0+VdiaR0)                              ) / Btht &
                                  - ( BR0      *(UR0+VdiaR0)   + BZ0      *(UZ0+VdiaZ0)   ) / Btht**2 * Btht_AZ__n
                    Vtht_A3__p  = + ( BR0      *(VdiaR0_A3__p) + BZ0      *(VdiaZ0_A3__p) ) / Btht &
                                  + ( BR0_A3   *(UR0+VdiaR0)   + BZ0_A3   *(UZ0+VdiaZ0)   ) / Btht &
                                  - ( BR0      *(UR0+VdiaR0)   + BZ0      *(UZ0+VdiaZ0)   ) / Btht**2 * Btht_A3
                    Vtht_A3__n  = + ( BR0      *(VdiaR0_A3__n) + BZ0      *(VdiaZ0_A3__n) ) / Btht
                    Vtht_rho__p = + ( BR0*(UR0+VdiaR0_rho__p) + BZ0*(UZ0+VdiaZ0_rho__p) ) / Btht
                    Vtht_rho__n = + ( BR0*(UR0+VdiaR0_rho__n) + BZ0*(UZ0+VdiaZ0_rho__n) ) / Btht
                    Vtht_Ti__p  = + ( BR0*(UR0+VdiaR0_Ti__p ) + BZ0*(UZ0+VdiaZ0_Ti__p ) ) / Btht
                    Vtht_Ti__n  = + ( BR0*(UR0+VdiaR0_Ti__n ) + BZ0*(UZ0+VdiaZ0_Ti__n ) ) / Btht
                    Vtht_UR     = + ( BR0*UR ) / Btht
                    Vtht_UZ     = + ( BZ0*UZ ) / Btht
                    
                    Vneo_AR__p = -aki_neo_prof(ms,mt) * tau_IC*F0 / BB2    / Btht    * ( BR0*Bp0_AR*Ti0_Z    - BZ0*Bp0_AR*Ti0_R   ) &
                                 +aki_neo_prof(ms,mt) * tau_IC*F0 / BB2**2 / Btht    * ( BR0*Bp0*Ti0_Z       - BZ0*Bp0*Ti0_R      ) * BB2_AR__p &
                                 +aki_neo_prof(ms,mt) * tau_IC*F0 / BB2    / Btht**2 * ( BR0*Bp0*Ti0_Z       - BZ0*Bp0*Ti0_R      ) * Btht_AR__p
                    Vneo_AR__n = -aki_neo_prof(ms,mt) * tau_IC*F0 / BB2    / Btht    * (                     - BZ0_AR__n*Bp0*Ti0_R) &
                                 +aki_neo_prof(ms,mt) * tau_IC*F0 / BB2**2 / Btht    * ( BR0*Bp0*Ti0_Z       - BZ0*Bp0*Ti0_R      ) * BB2_AR__n &
                                 +aki_neo_prof(ms,mt) * tau_IC*F0 / BB2    / Btht**2 * ( BR0*Bp0*Ti0_Z       - BZ0*Bp0*Ti0_R      ) * Btht_AR__n
                    Vneo_AZ__p = -aki_neo_prof(ms,mt) * tau_IC*F0 / BB2    / Btht    * ( BR0*Bp0_AZ*Ti0_Z    - BZ0*Bp0_AZ*Ti0_R   ) &
                                 +aki_neo_prof(ms,mt) * tau_IC*F0 / BB2**2 / Btht    * ( BR0*Bp0*Ti0_Z       - BZ0*Bp0*Ti0_R      ) * BB2_AZ__p &
                                 +aki_neo_prof(ms,mt) * tau_IC*F0 / BB2    / Btht**2 * ( BR0*Bp0*Ti0_Z       - BZ0*Bp0*Ti0_R      ) * Btht_AZ__p
                    Vneo_AZ__n = -aki_neo_prof(ms,mt) * tau_IC*F0 / BB2    / Btht    * ( BR0_AZ__n*Bp0*Ti0_Z                      ) &
                                 +aki_neo_prof(ms,mt) * tau_IC*F0 / BB2**2 / Btht    * ( BR0*Bp0*Ti0_Z       - BZ0*Bp0*Ti0_R      ) * BB2_AZ__n &
                                 +aki_neo_prof(ms,mt) * tau_IC*F0 / BB2    / Btht**2 * ( BR0*Bp0*Ti0_Z       - BZ0*Bp0*Ti0_R      ) * Btht_AZ__n
                    Vneo_A3__p = -aki_neo_prof(ms,mt) * tau_IC*F0 / BB2    / Btht    * ( BR0_A3*Bp0*Ti0_Z    - BZ0_A3*Bp0*Ti0_R   ) &
                                 +aki_neo_prof(ms,mt) * tau_IC*F0 / BB2**2 / Btht    * ( BR0*Bp0*Ti0_Z       - BZ0*Bp0*Ti0_R      ) * BB2_A3 &
                                 +aki_neo_prof(ms,mt) * tau_IC*F0 / BB2    / Btht**2 * ( BR0*Bp0*Ti0_Z       - BZ0*Bp0*Ti0_R      ) * Btht_A3
                    Vneo_A3__n = 0.d0
                    Vneo_Ti__p = -aki_neo_prof(ms,mt) * tau_IC*F0 / BB2    / Btht    * ( BR0*Bp0*Ti_Z        - BZ0*Bp0*Ti_R       )
                    Vneo_Ti__n = 0.d0
                    
                    PneoR_AR__p  = +amu_neo_prof(ms,mt) * rho0 * BB2_AR__p/Btht**3 * BR0       * (Vtht        - Vneo      ) &
                                   +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BR0       * (Vtht_AR__p  - Vneo_AR__p) &
                                   -amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**4 * BR0       * (Vtht        - Vneo      ) * 3.0 * Btht_AR__p
                    PneoR_AR__n  = +amu_neo_prof(ms,mt) * rho0 * BB2_AR__n/Btht**3 * BR0       * (Vtht        - Vneo      ) &
                                   +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BR0       * (Vtht_AR__n  - Vneo_AR__n) &
                                   -amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**4 * BR0       * (Vtht        - Vneo      ) * 3.0 * Btht_AR__n
                    PneoR_AZ__p  = +amu_neo_prof(ms,mt) * rho0 * BB2_AZ__p/Btht**3 * BR0       * (Vtht        - Vneo      ) &
                                   +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BR0       * (Vtht_AZ__p  - Vneo_AZ__p) &
                                   -amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**4 * BR0       * (Vtht        - Vneo      ) * 3.0 * Btht_AZ__p
                    PneoR_AZ__n  = +amu_neo_prof(ms,mt) * rho0 * BB2_AZ__n/Btht**3 * BR0       * (Vtht        - Vneo      ) &
                                   +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BR0_AZ__n * (Vtht        - Vneo      ) &
                                   +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BR0       * (Vtht_AZ__n  - Vneo_AZ__n) &
                                   -amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**4 * BR0       * (Vtht        - Vneo      ) * 3.0 * Btht_AZ__n
                    PneoR_A3__p  = +amu_neo_prof(ms,mt) * rho0 * BB2_A3   /Btht**3 * BR0       * (Vtht        - Vneo      ) &
                                   +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BR0_A3    * (Vtht        - Vneo      ) &
                                   +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BR0       * (Vtht_A3__p  - Vneo_A3__p) &
                                   -amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**4 * BR0       * (Vtht        - Vneo      ) * 3.0 * Btht_A3
                    PneoR_A3__n  = +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BR0       * (Vtht_A3__n  - Vneo_A3__n)
                    PneoR_rho__p = +amu_neo_prof(ms,mt) * rho  * BB2      /Btht**3 * BR0       * (Vtht        - Vneo      ) &
                                   +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BR0       * (Vtht_rho__p             )
                    PneoR_rho__n = +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BR0       * (Vtht_rho__n             )
                    PneoR_Ti__p  = +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BR0       * (Vtht_Ti__p  - Vneo_Ti__p)
                    PneoR_Ti__n  = +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BR0       * (Vtht_Ti__n  - Vneo_Ti__n)
                    PneoR_UR     = +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BR0       * (Vtht_UR                 )
                    PneoR_UZ     = +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BR0       * (Vtht_UZ                 )
                    
                    PneoZ_AR__p  = +amu_neo_prof(ms,mt) * rho0 * BB2_AR__p/Btht**3 * BZ0       * (Vtht        - Vneo      ) &
                                   +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BZ0       * (Vtht_AR__p  - Vneo_AR__p) &
                                   -amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**4 * BZ0       * (Vtht        - Vneo      ) * 3.0 * Btht_AR__p
                    PneoZ_AR__n  = +amu_neo_prof(ms,mt) * rho0 * BB2_AR__n/Btht**3 * BZ0       * (Vtht        - Vneo      ) &
                                   +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BZ0_AR__n * (Vtht        - Vneo      ) &
                                   +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BZ0       * (Vtht_AR__n  - Vneo_AR__n) &
                                   -amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**4 * BZ0       * (Vtht        - Vneo      ) * 3.0 * Btht_AR__n
                    PneoZ_AZ__p  = +amu_neo_prof(ms,mt) * rho0 * BB2_AZ__p/Btht**3 * BZ0       * (Vtht        - Vneo      ) &
                                   +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BZ0       * (Vtht_AZ__p  - Vneo_AZ__p) &
                                   -amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**4 * BZ0       * (Vtht        - Vneo      ) * 3.0 * Btht_AZ__p
                    PneoZ_AZ__n  = +amu_neo_prof(ms,mt) * rho0 * BB2_AZ__n/Btht**3 * BZ0       * (Vtht        - Vneo      ) &
                                   +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BZ0       * (Vtht_AZ__n  - Vneo_AZ__n) &
                                   -amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**4 * BZ0       * (Vtht        - Vneo      ) * 3.0 * Btht_AZ__n
                    PneoZ_A3__p  = +amu_neo_prof(ms,mt) * rho0 * BB2_A3   /Btht**3 * BZ0       * (Vtht        - Vneo      ) &
                                   +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BZ0_A3    * (Vtht        - Vneo      ) &
                                   +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BZ0       * (Vtht_A3__p  - Vneo_A3__p) &
                                   -amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**4 * BZ0       * (Vtht        - Vneo      ) * 3.0 * Btht_A3
                    PneoZ_A3__n  = +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BZ0       * (Vtht_A3__n  - Vneo_A3__n)
                    PneoZ_rho__p = +amu_neo_prof(ms,mt) * rho  * BB2      /Btht**3 * BZ0       * (Vtht        - Vneo      ) &
                                   +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BZ0       * (Vtht_rho__p             )
                    PneoZ_rho__n = +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BZ0       * (Vtht_rho__n             )
                    PneoZ_Ti__p  = +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BZ0       * (Vtht_Ti__p  - Vneo_Ti__p)
                    PneoZ_Ti__n  = +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BZ0       * (Vtht_Ti__n  - Vneo_Ti__n)
                    PneoZ_UR     = +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BZ0       * (Vtht_UR                 )
                    PneoZ_UZ     = +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BZ0       * (Vtht_UZ                 )
                  else
                    PneoR_AR__p  = 0.d0 ; PneoZ_AR__p  = 0.d0
                    PneoR_AR__n  = 0.d0 ; PneoZ_AR__n  = 0.d0
                    PneoR_AZ__p  = 0.d0 ; PneoZ_AZ__p  = 0.d0
                    PneoR_AZ__n  = 0.d0 ; PneoZ_AZ__n  = 0.d0
                    PneoR_A3__p  = 0.d0 ; PneoZ_A3__p  = 0.d0
                    PneoR_A3__n  = 0.d0 ; PneoZ_A3__n  = 0.d0
                    PneoR_rho__p = 0.d0 ; PneoZ_rho__p = 0.d0
                    PneoR_rho__n = 0.d0 ; PneoZ_rho__n = 0.d0
                    PneoR_Ti__p  = 0.d0 ; PneoZ_Ti__p  = 0.d0
                    PneoR_Ti__n  = 0.d0 ; PneoZ_Ti__n  = 0.d0
                    PneoR_UR     = 0.d0 ; PneoZ_UR     = 0.d0
                    PneoR_UZ     = 0.d0 ; PneoZ_UZ     = 0.d0
                  endif

                  Qconv_UR_AR__p  =                                           - v * rho0 * ( VdiaGradUR_AR__p  - VdiaP0_AR__p *Up0 / R )
                  Qconv_UR_AR__n  =                                           - v * rho0 * ( VdiaGradUR_AR__n  - VdiaP0_AR__n *Up0 / R )
                  Qconv_UR_AZ__p  =                                           - v * rho0 * ( VdiaGradUR_AZ__p  - VdiaP0_AZ__p *Up0 / R )
                  Qconv_UR_AZ__n  =                                           - v * rho0 * ( VdiaGradUR_AZ__n  - VdiaP0_AZ__n *Up0 / R )
                  Qconv_UR_A3__p  =                                           - v * rho0 * ( VdiaGradUR_A3__p  - VdiaP0_A3__p *Up0 / R )
                  Qconv_UR_A3__n  =                                           - v * rho0 * ( VdiaGradUR_A3__n  - VdiaP0_A3__n *Up0 / R )
                  Qconv_UR_rho__p = - v * rho  * ( UgradUR - Up0**2  / R )    - v * rho0 * ( VdiaGradUR_rho__p - VdiaP0_rho__p*Up0 / R ) &
                                                                              - v * rho  * ( VdiaGradUR        - VdiaP0       *Up0 / R )
                  Qconv_UR_rho__n =                                           - v * rho0 * ( VdiaGradUR_rho__n - VdiaP0_rho__n*Up0 / R )
                  Qconv_UR_Ti__p  =                                           - v * rho0 * ( VdiaGradUR_Ti__p  - VdiaP0_Ti__p *Up0 / R )
                  Qconv_UR_Ti__n  =                                           - v * rho0 * ( VdiaGradUR_Ti__n  - VdiaP0_Ti__n *Up0 / R )
                  Qconv_UR_UR__p  = - v * rho0 * ( UgradUR_UR__p         )    - v * rho0 * ( VdiaGradUR_UR__p                          )
                  Qconv_UR_UR__n  = - v * rho0 * ( UgradUR_UR__n         )    - v * rho0 * ( VdiaGradUR_UR__n                          )
                  Qconv_UR_UZ__p  = - v * rho0 * ( UgradUR_UZ            )
                  Qconv_UR_UZ__n  = 0.d0
                  Qconv_UR_Up__p  = - v * rho0 * ( UgradUR_Up-2.0*Up*Up0**2/R)- v * rho0 * (                   - VdiaP0       *Up  / R )
                  Qconv_UR_Up__n  = 0.d0

                  Qconv_UZ_AR__p  =                                           - v * rho0 * VdiaGradUZ_AR__p 
                  Qconv_UZ_AR__n  =                                           - v * rho0 * VdiaGradUZ_AR__n 
                  Qconv_UZ_AZ__p  =                                           - v * rho0 * VdiaGradUZ_AZ__p 
                  Qconv_UZ_AZ__n  =                                           - v * rho0 * VdiaGradUZ_AZ__n 
                  Qconv_UZ_A3__p  =                                           - v * rho0 * VdiaGradUZ_A3__p 
                  Qconv_UZ_A3__n  =                                           - v * rho0 * VdiaGradUZ_A3__n 
                  Qconv_UZ_rho__p = - v * rho  * UgradUZ                      - v * rho0 * VdiaGradUZ_rho__p - v * rho  * VdiaGradUZ
                  Qconv_UZ_rho__n =                                           - v * rho0 * VdiaGradUZ_rho__n
                  Qconv_UZ_Ti__p  =                                           - v * rho0 * VdiaGradUZ_Ti__p
                  Qconv_UZ_Ti__n  =                                           - v * rho0 * VdiaGradUZ_Ti__n
                  Qconv_UZ_UR__p  = - v * rho0 * UgradUZ_UR
                  Qconv_UZ_UR__n  = 0.d0
                  Qconv_UZ_UZ__p  = - v * rho0 * UgradUZ_UZ__p                - v * rho0 * VdiaGradUZ_UZ__p
                  Qconv_UZ_UZ__n  = - v * rho0 * UgradUZ_UZ__n                - v * rho0 * VdiaGradUZ_UZ__n
                  Qconv_UZ_Up__p  = - v * rho0 * UgradUZ_Up
                  Qconv_UZ_Up__n  = 0.d0

                  Qconv_Up_AR__p  =                                           - v * rho0 * ( VdiaGradUp_AR__p  + VdiaP0_AR__p *UR0 / R )
                  Qconv_Up_AR__n  =                                           - v * rho0 * ( VdiaGradUp_AR__n  + VdiaP0_AR__n *UR0 / R )
                  Qconv_Up_AZ__p  =                                           - v * rho0 * ( VdiaGradUp_AZ__p  + VdiaP0_AZ__p *UR0 / R )
                  Qconv_Up_AZ__n  =                                           - v * rho0 * ( VdiaGradUp_AZ__n  + VdiaP0_AZ__n *UR0 / R )
                  Qconv_Up_A3__p  =                                           - v * rho0 * ( VdiaGradUp_A3__p  + VdiaP0_A3__p *UR0 / R )
                  Qconv_Up_A3__n  =                                           - v * rho0 * ( VdiaGradUp_A3__n  + VdiaP0_A3__n *UR0 / R )
                  Qconv_Up_rho__p = - v * rho  * ( UgradUp + UR0*Up0 / R )    - v * rho0 * ( VdiaGradUp_rho__p + VdiaP0_rho__p*UR0 / R ) &
                                                                              - v * rho  * ( VdiaGradUp        + VdiaP0       *UR0 / R )
                  Qconv_Up_rho__n =                                           - v * rho0 * ( VdiaGradUp_rho__n + VdiaP0_rho__n*UR0 / R )
                  Qconv_Up_Ti__p  =                                           - v * rho0 * ( VdiaGradUp_Ti__p  + VdiaP0_Ti__p *UR0 / R )
                  Qconv_Up_Ti__n  =                                           - v * rho0 * ( VdiaGradUp_Ti__n  + VdiaP0_Ti__n *UR0 / R )
                  Qconv_Up_UR__p  = - v * rho0 * ( UgradUp_UR    + UR*Up0 / R)- v * rho0 * (                   + VdiaP0       *UR  / R )
                  Qconv_Up_UR__n  = 0.d0
                  Qconv_Up_UZ__p  = - v * rho0 * ( UgradUp_UZ                )
                  Qconv_Up_UZ__n  = 0.d0
                  Qconv_Up_Up__p  = - v * rho0 * ( UgradUp_Up__p + UR0*Up / R)- v * rho0 * ( VdiaGradUp_Up__p                          )
                  Qconv_Up_Up__n  = - v * rho0 * ( UgradUp_Up__n             )- v * rho0 * ( VdiaGradUp_Up__n                          )

                  BgradVstar_AR__p = BR0_AR    * v_R
                  BgradVstar_AR__n = BZ0_AR__n * v_Z
                  BgradVstar_AR__k = Bp0_AR    * v_p__k / R

                  BgradVstar_AZ__n = BR0_AZ__n * v_R
                  BgradVstar_AZ__p = BZ0_AZ    * v_Z
                  BgradVstar_AZ__k = Bp0_AZ    * v_p__k / R

                  BgradVstar_A3__p = BR0_A3 * v_R  + BZ0_A3 * v_Z
                  BgradVstar_A3__k = Bp0_A3 * v_p__k / R

                  VdiaGradVstar_AR__p   = VdiaR0_AR__p * v_R + VdiaZ0_AR__p * v_Z
                  VdiaGradVstar_AR__n   = VdiaR0_AR__n * v_R + VdiaZ0_AR__n * v_Z
                  VdiaGradVstar_AR__k   = VdiaP0_AR__p * v_p__k / R
                  VdiaGradVstar_AR__kn  = VdiaP0_AR__n * v_p__k / R

                  VdiaGradVstar_AZ__p   = VdiaR0_AZ__p * v_R + VdiaZ0_AZ__p * v_Z
                  VdiaGradVstar_AZ__n   = VdiaR0_AZ__n * v_R + VdiaZ0_AZ__n * v_Z
                  VdiaGradVstar_AZ__k   = VdiaP0_AZ__p * v_p__k / R
                  VdiaGradVstar_AZ__kn  = VdiaP0_AZ__n * v_p__k / R

                  VdiaGradVstar_A3__p   = VdiaR0_A3__p * v_R + VdiaZ0_A3__p * v_Z
                  VdiaGradVstar_A3__n   = VdiaR0_A3__n * v_R + VdiaZ0_A3__n * v_Z
                  VdiaGradVstar_A3__k   = VdiaP0_A3__p * v_p__k / R
                  VdiaGradVstar_A3__kn  = VdiaP0_A3__n * v_p__k / R

                  VdiaGradVstar_rho__p  = VdiaR0_rho__p * v_R + VdiaZ0_rho__p * v_Z
                  VdiaGradVstar_rho__n  = VdiaR0_rho__n * v_R + VdiaZ0_rho__n * v_Z
                  VdiaGradVstar_rho__k  = VdiaP0_rho__p * v_p__k / R
                  VdiaGradVstar_rho__kn = VdiaP0_rho__n * v_p__k / R

                  VdiaGradVstar_Ti__p   = VdiaR0_Ti__p * v_R + VdiaZ0_Ti__p * v_Z
                  VdiaGradVstar_Ti__n   = VdiaR0_Ti__n * v_R + VdiaZ0_Ti__n * v_Z
                  VdiaGradVstar_Ti__k   = VdiaP0_Ti__p * v_p__k / R
                  VdiaGradVstar_Ti__kn  = VdiaP0_Ti__n * v_p__k / R

                  JxB_UR_AR__p  = - BR0_AR * BgradVstar__p - BR0 * BgradVstar_AR__p &
                                  - 2.0*Bp0*Bp0_AR * v / R + ( v_R + v / R ) * ( 0.5d0 * BB2_AR__p )
                  JxB_UR_AR__n  = - BR0 * BgradVstar_AR__n + ( v_R + v / R ) * ( 0.5d0 * BB2_AR__n )
                  JxB_UR_AR__k  = - BR0_AR * BgradVstar__k - BR0 * BgradVstar_AR__k
                  JxB_UR_AR__kn = 0.d0

                  JxB_UR_AZ__p  = - BR0       * BgradVstar_AZ__p - 2.0*Bp0*Bp0_AZ * v / R + ( v_R + v / R ) * ( 0.5d0 * BB2_AZ__p )
                  JxB_UR_AZ__n  = - BR0_AZ__n * BgradVstar__p    - BR0 * BgradVstar_AZ__n + ( v_R + v / R ) * ( 0.5d0 * BB2_AZ__n )
                  JxB_UR_AZ__k  = - BR0       * BgradVstar_AZ__k
                  JxB_UR_AZ__kn = - BR0_AZ__n * BgradVstar__k

                  JxB_UR_A3__p  = - BR0_A3 * BgradVstar__p - BR0 * BgradVstar_A3__p &
                                  - 2.0*Bp0*Bp0_A3 * v / R + ( v_R + v / R ) * ( 0.5d0 * BB2_A3 )
                  JxB_UR_A3__n  = 0.d0
                  JxB_UR_A3__k  = - BR0_A3 * BgradVstar__k - BR0 * BgradVstar_A3__k
                  JxB_UR_A3__kn = 0.d0

                  JxB_UZ_AR__p  = - BZ0       * BgradVstar_AR__p                          + v_Z * ( 0.5d0 * BB2_AR__p)
                  JxB_UZ_AR__n  = - BZ0_AR__n * BgradVstar__p - BZ0 * BgradVstar_AR__n    + v_Z * ( 0.5d0 * BB2_AR__n)
                  JxB_UZ_AR__k  = - BZ0       * BgradVstar_AR__k
                  JxB_UZ_AR__kn = - BZ0_AR__n * BgradVstar__k

                  JxB_UZ_AZ__p  = - BZ0_AZ * BgradVstar__p - BZ0 * BgradVstar_AZ__p + v_Z * ( 0.5d0 * BB2_AZ__p)
                  JxB_UZ_AZ__n  = - BZ0    * BgradVstar_AZ__n                       + v_Z * ( 0.5d0 * BB2_AZ__n)
                  JxB_UZ_AZ__k  = - BZ0_AZ * BgradVstar__k - BZ0 * BgradVstar_AZ__k
                  JxB_UZ_AZ__kn = 0.d0

                  JxB_UZ_A3__p  = - BZ0_A3 * BgradVstar__p - BZ0 * BgradVstar_A3__p + v_Z * ( 0.5d0 * BB2_A3)
                  JxB_UZ_A3__n  = 0.d0
                  JxB_UZ_A3__k  = - BZ0_A3 * BgradVstar__k - BZ0 * BgradVstar_A3__k
                  JxB_UZ_A3__kn = 0.d0

                  JxB_Up_AR__p  = - Bp0_AR * BgradVstar__p - Bp0 * BgradVstar_AR__p + Bp0_AR*BR0* v / R + Bp0*BR0_AR* v / R
                  JxB_Up_AR__n  = - Bp0    * BgradVstar_AR__n
                  JxB_Up_AR__k  = - Bp0_AR * BgradVstar__k - Bp0 * BgradVstar_AR__k + ( v_p__k/ R ) * ( 0.5d0 * BB2_AR__p )
                  JxB_Up_AR__kn = + ( v_p__k/ R ) * ( 0.5d0 * BB2_AR__n )

                  JxB_Up_AZ__p  = - Bp0_AZ * BgradVstar__p - Bp0 * BgradVstar_AZ__p + Bp0_AZ*BR0* v / R
                  JxB_Up_AZ__n  = - Bp0    * BgradVstar_AZ__n                       + Bp0*BR0_AZ__n* v / R
                  JxB_Up_AZ__k  = - Bp0_AZ * BgradVstar__k - Bp0 * BgradVstar_AZ__k + ( v_p__k/ R ) * ( 0.5d0 * BB2_AZ__p )
                  JxB_Up_AZ__kn = + ( v_p__k/ R ) * ( 0.5d0 * BB2_AZ__n )

                  JxB_Up_A3__p  = - Bp0_A3 * BgradVstar__p - Bp0 * BgradVstar_A3__p + Bp0_A3*BR0* v / R + Bp0*BR0_A3* v / R
                  JxB_Up_A3__n  = 0.d0
                  JxB_Up_A3__k  = - Bp0_A3 * BgradVstar__k - Bp0 * BgradVstar_A3__k + ( v_p__k/ R ) * ( 0.5d0 * BB2_A3 )
                  JxB_Up_A3__kn = 0.d0

                  gradRho_gradVstar_rho__p  = rho_R * v_R  + rho_Z * v_Z
                  gradRho_gradVstar_rho__kn = (rho_p__n / R) * (v_p__k  / R)

                  gradTi_gradVstar_Ti__p  = Ti_R * v_R  + Ti_Z * v_Z
                  gradTi_gradVstar_Ti__kn = (Ti_p__n / R) * (v_p__k  / R)

                  gradTe_gradVstar_Te__p  = Te_R * v_R  + Te_Z * v_Z
                  gradTe_gradVstar_Te__kn = (Te_p__n / R) * (v_p__k  / R)
#ifdef TEMPERATURE_DEPENDENT_UPDATE
                  eta_T_T    = deta_dT * Te
                  eta_R_T    = d2eta_d2T * Te * Te0_R + deta_dT * Te_R
                  eta_Z_T    = d2eta_d2T * Te * Te0_Z + deta_dT * Te_Z
                  eta_p_T__p = d2eta_d2T * Te * Te0_p
                  eta_p_T__n = deta_dT * Te_p__n

                  visco_T_T  = dvisco_dT * Te
                  visco_R_T  = d2visco_d2T * Te * Te0_R + dvisco_dT * Te_R
                  visco_Z_T  = d2visco_d2T * Te * Te0_Z + dvisco_dT * Te_Z
                  visco_p_T__p = d2visco_d2T * Te * Te0_p
                  visco_p_T__n = dvisco_dT * Te_p__n
#else
                  eta_T_T    = 0.d0
                  eta_R_T    = 0.d0
                  eta_Z_T    = 0.d0
                  eta_p_T__p = 0.d0
                  eta_p_T__n = 0.d0

                  deta_dT_ohm = 0.d0

                  visco_T_T    = 0.d0
                  visco_R_T    = 0.d0
                  visco_Z_T    = 0.d0
                  visco_p_T__p = 0.d0
                  visco_p_T__n = 0.d0

                  dZKi_par_dT = 0.d0
                  dZKe_par_dT = 0.d0
#endif
                  JR0_AR__p  = 0.d0  ; JR0_AR__n  = 0.d0  ; JR0_AR__nn = 0.d0
                  JR0_AZ__p  = 0.d0  ; JR0_AZ__n  = 0.d0  ; JR0_AZ__nn = 0.d0
                  JR0_A3__p  = 0.d0  ; JR0_A3__n  = 0.d0  ; JR0_A3__nn = 0.d0

                  JZ0_AR__p  = 0.d0  ; JZ0_AR__n  = 0.d0  ; JZ0_AR__nn = 0.d0
                  JZ0_AZ__p  = 0.d0  ; JZ0_AZ__n  = 0.d0  ; JZ0_AZ__nn = 0.d0
                  JZ0_A3__p  = 0.d0  ; JZ0_A3__n  = 0.d0  ; JZ0_A3__nn = 0.d0

                  Jp0_AR__p  = 0.d0  ; Jp0_AR__n  = 0.d0  ; Jp0_AR__nn = 0.d0
                  Jp0_AZ__p  = 0.d0  ; Jp0_AZ__n  = 0.d0  ; Jp0_AZ__nn = 0.d0
                  Jp0_A3__p  = 0.d0  ; Jp0_A3__n  = 0.d0  ; Jp0_A3__nn = 0.d0

                  JJ2_AR__p  = 0.d0  ; JJ2_AR__n  = 0.d0  ; JJ2_AR__nn = 0.d0
                  JJ2_AZ__p  = 0.d0  ; JJ2_AZ__n  = 0.d0  ; JJ2_AZ__nn = 0.d0
                  JJ2_A3__p  = 0.d0  ; JJ2_A3__n  = 0.d0  ; JJ2_A3__nn = 0.d0

                  if(eta_ohmic .gt. 1.d-12)then

                    JR0_AR__p  = - AR_ZZ
                    JR0_AR__n  = 0.d0
                    JR0_AR__nn = - AR_pp / R**2
                    JR0_AZ__p  = AZ_RZ
                    JR0_AZ__n  = 0.0d0
                    JR0_AZ__nn = 0.0d0
                    JR0_A3__p  = dF_dpsi(ms,mt) * A3_Z / R
                    JR0_A3__n  = A3_Rp / R**2
                    JR0_A3__nn = 0.d0

                    JZ0_AR__p  = AR_RZ / R**2 + AR_Z / R
                    JZ0_AR__n  = 0.0d0
                    JZ0_AR__nn = 0.0d0
                    JZ0_AZ__p  = - AZ_RR - AZ_R / R
                    JZ0_AZ__n  = 0.d0
                    JZ0_AZ__nn = - AZ_pp / R**2
                    JZ0_A3__p  = - (dF_dpsi(ms,mt)*A3_R) / R ! + Fprofile(ms,mt)/R**2 - Fprofile(ms,mt)/R**2
                    JZ0_A3__n  = A3_Zp / R**2
                    JZ0_A3__nn = 0.d0

                    Jp0_AR__p  = 0.0d0
                    Jp0_AR__n  = AR_Rp / R - AR_p__n / R**2
                    Jp0_AR__nn = 0.d0
                    Jp0_AZ__p  = 0.0d0
                    Jp0_AZ__n  = AZ_Zp / R
                    Jp0_AZ__nn = 0.d0
                    Jp0_A3__p  = - A3_RR / R + A3_R / R**2 - A3_ZZ / R
                    Jp0_A3__n  = 0.0d0
                    Jp0_A3__nn = 0.0d0

                    JJ2_AR__p  = 2.0d0 * (JR0*JR0_AR__p  + JZ0*JZ0_AR__p  + Jp0*Jp0_AR__p)
                    JJ2_AR__n  = 2.0d0 * (JR0*JR0_AR__n  + JZ0*JZ0_AR__n  + Jp0*Jp0_AR__n)
                    JJ2_AR__nn = 2.0d0 * (JR0*JR0_AR__nn + JZ0*JZ0_AR__nn + Jp0*Jp0_AR__nn)

                    JJ2_AZ__p  = 2.0d0 * (JR0*JR0_AZ__p  + JZ0*JZ0_AZ__p  + Jp0*Jp0_AZ__p)
                    JJ2_AZ__n  = 2.0d0 * (JR0*JR0_AZ__n  + JZ0*JZ0_AZ__n  + Jp0*Jp0_AZ__n)
                    JJ2_AZ__nn = 2.0d0 * (JR0*JR0_AZ__nn + JZ0*JZ0_AZ__nn + Jp0*Jp0_AZ__nn)

                    JJ2_A3__p  = 2.0d0 * (JR0*JR0_A3__p  + JZ0*JZ0_A3__p  + Jp0*Jp0_A3__p)
                    JJ2_A3__n  = 2.0d0 * (JR0*JR0_A3__n  + JZ0*JZ0_A3__n  + Jp0*Jp0_A3__n)
                    JJ2_A3__nn = 2.0d0 * (JR0*JR0_A3__nn + JZ0*JZ0_A3__nn + Jp0*Jp0_A3__nn)

                  endif

#ifdef ENABLE_VISCO_PAR_AND_PARPAR
                  !  --- For visco_par and visco_par_par
                  grad_BR_AR__R    = 0.d0
                  grad_BZ_AR__R__n = -1.d0 / R**2 * AR_p__n + 1.d0 / R * AR_Rp
                  grad_Bp_AR__R    = -AR_RZ

                  grad_BR_AR__Z    = 0.d0
                  grad_BZ_AR__Z__n = 1.d0 / R * AR_Zp
                  grad_Bp_AR__Z    = -AR_ZZ

                  grad_BR_AR__phi = 0.d0
                  grad_BZ_AR__phi__nn = 1.d0 / R * AR_pp
                  grad_Bp_AR__phi__n = -AR_Zp

                  grad_BR_AZ__R__n = 1.d0 / R**2 * AZ_p__n - 1.d0 / R * AZ_Rp
                  grad_BZ_AZ__R    = 0.d0
                  grad_Bp_AZ__R    = AZ_RR

                  grad_BR_AZ__Z__n = - 1.d0 / R * AZ_Zp
                  grad_BZ_AZ__Z    = 0.d0
                  grad_Bp_AZ__Z    = AZ_RZ

                  grad_BR_AZ__phi__nn = - 1.d0 / R * AZ_pp
                  grad_BZ_AZ__phi    = 0.d0
                  grad_Bp_AZ__phi__n    = AZ_Rp

                  grad_BR_A3__R    = -1.d0 / R**2 * A3_Z + 1.d0 / R * A3_RZ
                  grad_BZ_A3__R    = 1.d0 / R**2 * A3_R - 1.d0 / R * A3_RR
                  grad_Bp_A3__R    = 0.d0

                  grad_BR_A3__Z    =  1.d0 / R * A3_ZZ
                  grad_BZ_A3__Z    = -1.d0 / R * A3_RZ
                  grad_Bp_A3__Z    = 0.d0

                  grad_BR_A3__phi__n    =  1.d0 / R * A3_Zp
                  grad_BZ_A3__phi__n    = -1.d0 / R * A3_Rp
                  grad_Bp_A3__phi    = 0.d0

                  B_AR_dot_U__p = BR0_AR * UR0 + Bp0_AR * Up0
                  B_AR_dot_U__n = BZ0_AR__n * UZ0 
                  B_AZ_dot_U__p = BZ0_AZ * UZ0 + Bp0_AZ * Up0
                  B_AZ_dot_U__n = BR0_AZ__n * UR0
                  B_A3_dot_U__p = BR0_A3 * UR0 + BZ0_A3 * UZ0 + Bp0_A3 * Up0
                  B_A3_dot_U__n = 0.d0

                  grad_B_AR_dot_U__R__p = grad_BR_AR__R * UR0 + grad_Bp_AR__R * Up0 &
                                         + BR0_AR * UR0_R + Bp0_AR * Up0_R
                  grad_B_AR_dot_U__R__n = grad_BZ_AR__R__n * UZ0  &
                                         + BZ0_AR__n * UZ0_R
                  grad_B_AZ_dot_U__R__p = grad_BZ_AZ__R * UZ0 + grad_Bp_AZ__R * Up0 &
                                         + BZ0_AZ * UZ0_R + Bp0_AZ * Up0_R
                  grad_B_AZ_dot_U__R__n = grad_BR_AZ__R__n * UR0 &
                                         + BR0_AZ__n * UR0_R
                  grad_B_A3_dot_U__R__p = grad_BR_A3__R * UR0 + grad_BZ_A3__R * UZ0 + grad_Bp_A3__R * Up0    &
                                         + BR0_A3 * UR0_R + BZ0_A3 * UZ0_R + Bp0_A3 * Up0_R
                  grad_B_A3_dot_U__R__n = 0.d0

                  grad_B_AR_dot_U__Z__p = grad_BR_AR__Z * UR0 + grad_Bp_AR__Z * Up0 &
                                         + BR0_AR * UR0_Z + Bp0_AR * Up0_Z
                  grad_B_AR_dot_U__Z__n = grad_BZ_AR__Z__n * UZ0 &
                                         + BZ0_AR__n * UZ0_Z
                  grad_B_AZ_dot_U__Z__p = grad_BZ_AZ__Z * UZ0 + grad_Bp_AZ__Z * Up0 &
                                         + BZ0_AZ * UZ0_Z + Bp0_AZ * Up0_Z
                  grad_B_AZ_dot_U__Z__n = grad_BR_AZ__Z__n * UR0 &
                                         + BR0_AZ__n * UR0_Z
                  grad_B_A3_dot_U__Z__p = grad_BR_A3__Z * UR0 + grad_BZ_A3__Z * UZ0 + grad_Bp_A3__Z * Up0    &
                                         + BR0_A3 * UR0_Z + BZ0_A3 * UZ0_Z + Bp0_A3 * Up0_Z
                  grad_B_A3_dot_U__Z__n = 0.d0

                  grad_B_AR_dot_U__phi__p = grad_BR_AR__phi * UR0 &
                                            + BR0_AR * UR0_p + Bp0_AR * Up0_p
                  grad_B_AR_dot_U__phi__n = grad_Bp_AR__phi__n * Up0 &
                                            + BZ0_AR__n * UZ0_p 
                  grad_B_AR_dot_U__phi__nn = grad_BZ_AR__phi__nn * UZ0 
                                            
                  grad_B_AZ_dot_U__phi__p = grad_BZ_AZ__phi * UZ0  &
                                            + BZ0_AZ * UZ0_p + Bp0_AZ* Up0_p
                  grad_B_AZ_dot_U__phi__n = grad_Bp_AZ__phi__n * Up0 &
                                            + BR0_AZ__n * UR0_p
                  grad_B_AZ_dot_U__phi__nn = grad_BR_AZ__phi__nn * UR0 

                  grad_B_A3_dot_U__phi__p = grad_Bp_A3__phi * Up0 &
                                            + BR0_A3 * UR0_p + BZ0_A3 * UZ0_p + Bp0_A3 * Up0_p
                  grad_B_A3_dot_U__phi__n = grad_BR_A3__phi__n * UR0 + grad_BZ_A3__phi__n * UZ0
                  grad_B_A3_dot_U__phi__nn = 0.d0

                  grad_BB2_AR__R__p = 2.d0 * (grad_BR_AR__R * BR0 + BR0_AR * BR0_R + &
                                                grad_Bp_AR__R * Bp0 + Bp0_AR * Bp0_R )
                  grad_BB2_AR__R__n = 2.d0 * (grad_BZ_AR__R__n * BZ0 + BZ0_AR__n * BZ0_R)
                  grad_BB2_AR__R__nn = 0.d0
                  grad_BB2_AZ__R__p = 2.d0 * (grad_BZ_AZ__R * BZ0 + BZ0_AZ * BZ0_R + &
                                                grad_Bp_AZ__R * Bp0 + Bp0_AZ * Bp0_R)
                  grad_BB2_AZ__R__n = 2.d0 * (grad_BR_AZ__R__n * BR0 + BR0_AZ__n * BR0_R)
                  grad_BB2_AZ__R__nn = 0.d0
                  grad_BB2_A3__R__p = 2.d0 * (grad_BR_A3__R * BR0 + BR0_A3 * BR0_R + &
                                                grad_BZ_A3__R * BZ0 + BZ0_A3 * BZ0_R + &
                                                grad_Bp_A3__R * Bp0 + Bp0_A3 * Bp0_R)
                  grad_BB2_A3__R__n = 0.d0
                  grad_BB2_A3__R__nn = 0.d0

                  grad_BB2_AR__Z__p = 2.d0 * (grad_BR_AR__Z * BR0 + BR0_AR * BR0_Z + &
                                                grad_Bp_AR__Z * Bp0 + Bp0_AR * Bp0_Z )
                  grad_BB2_AR__Z__n = 2.d0 * (grad_BZ_AR__Z__n * BZ0 + BZ0_AR__n * BZ0_Z)
                  grad_BB2_AR__Z__nn = 0.d0
                  grad_BB2_AZ__Z__p = 2.d0 * (grad_BZ_AZ__Z * BZ0 + BZ0_AZ * BZ0_Z + &
                                                grad_Bp_AZ__Z * Bp0 + Bp0_AZ * Bp0_Z)
                  grad_BB2_AZ__Z__n = 2.d0 * (grad_BR_AZ__Z__n * BR0 + BR0_AZ__n * BR0_Z)
                  grad_BB2_AZ__Z__nn = 0.d0
                  grad_BB2_A3__Z__p = 2.d0 * (grad_BR_A3__Z * BR0 + BR0_A3 * BR0_Z + &
                                                grad_BZ_A3__Z * BZ0 + BZ0_A3 * BZ0_Z + &
                                                grad_Bp_A3__Z * Bp0 + Bp0_A3 * Bp0_Z)
                  grad_BB2_A3__Z__n = 0.d0
                  grad_BB2_A3__Z__nn = 0.d0

                  grad_BB2_AR__phi__p = 2.d0 * (grad_BR_AR__phi * BR0 + BR0_AR * BR0_p + &
                                                Bp0_AR * Bp0_p )
                  grad_BB2_AR__phi__n = 2.d0 * (BZ0_AR__n * BZ0_p + &
                                                grad_Bp_AR__phi__n * Bp0)
                  grad_BB2_AR__phi__nn = 2.d0 * (grad_BZ_AR__phi__nn * BZ0)
                  grad_BB2_AZ__phi__p = 2.d0 * (grad_BZ_AZ__phi * BZ0 + BZ0_AZ * BZ0_p + &
                                                Bp0_AZ * Bp0_p)
                  grad_BB2_AZ__phi__n = 2.d0 * (BR0_AZ__n * BR0_p + &
                                                grad_Bp_AZ__phi__n * Bp0)
                  grad_BB2_AZ__phi__nn = 2.d0 * (grad_BR_AZ__phi__nn * BR0)
                  grad_BB2_A3__phi__p = 2.d0 * (BR0_A3 * BR0_p + &
                                                BZ0_A3 * BZ0_p + &
                                                grad_Bp_A3__phi * Bp0 + Bp0_A3 * Bp0_p)
                  grad_BB2_A3__phi__n = 2.d0 * (grad_BR_A3__phi__n * BR0 + &
                                                grad_BZ_A3__phi__n * BZ0)
                  grad_BB2_A3__phi__nn = 0.d0

                  B_dot_grad_BRVstar__p = BR0 * (v * BR0_R + BR0 * v_R) + &
                                          BZ0 * (v * BR0_Z + BR0 * v_Z) + &
                                          Bp0 * (v * BR0_p) / R
                  B_dot_grad_BRVstar__k = Bp0 * (BR0 * v_p__k) / R

                  B_dot_grad_BZVstar__p = BR0 * (v * BZ0_R + BZ0 * v_R) + &
                                          BZ0 * (v * BZ0_Z + BZ0 * v_Z) + &
                                          Bp0 * (v * BZ0_p) / R
                  B_dot_grad_BZVstar__k = Bp0 * (BZ0 * v_p__k) / R

                  B_dot_grad_BpVstar__p = BR0 * (v * Bp0_R + Bp0 * v_R) + &
                                          BZ0 * (v * Bp0_Z + Bp0 * v_Z) + &
                                          Bp0 * (v * Bp0_p) / R
                  B_dot_grad_BpVstar__k = Bp0 * (Bp0 * v_p__k) / R

                  B_dot_grad_BRVstar__AR__p =  BR0_AR * (v * BR0_R + BR0 * v_R) + &
                                               Bp0_AR * (v * BR0_p) / R + &
                                               BR0 * (v * grad_BR_AR__R + BR0_AR * v_R) + &
                                               BZ0 * (v * grad_BR_AR__Z + BR0_AR * v_Z) + &
                                               Bp0 * (v * grad_BR_AR__phi) / R
                  B_dot_grad_BRVstar__AR__n =  BZ0_AR__n * (v * BR0_Z + BR0 * v_Z) 
                  B_dot_grad_BRVstar__AR__k =  Bp0_AR * (BR0 * v_p__k) / R + &
                                               Bp0 * (BR0_AR * v_p__k) / R
                  B_dot_grad_BRVstar__AR__kn = 0.d0
                  B_dot_grad_BRVstar__AR__nn = 0.d0

                  B_dot_grad_BZVstar__AR__p =  BR0_AR * (v * BZ0_R + BZ0 * v_R) + &
                                               Bp0_AR * (v * BZ0_p) / R 
                  B_dot_grad_BZVstar__AR__n =  BZ0_AR__n * (v * BZ0_Z + BZ0 * v_Z) + &
                                               BR0 * (v * grad_BZ_AR__R__n + BZ0_AR__n * v_R) + &
                                               BZ0 * (v * grad_BZ_AR__Z__n + BZ0_AR__n * v_Z) 
                  B_dot_grad_BZVstar__AR__k =  Bp0_AR * (BZ0 * v_p__k) / R 
                  B_dot_grad_BZVstar__AR__kn = Bp0 * (BZ0_AR__n * v_p__k) / R
                  B_dot_grad_BZVstar__AR__nn = Bp0 * (v * grad_BZ_AR__phi__nn) / R

                  B_dot_grad_BpVstar__AR__p =  BR0_AR * (v * Bp0_R + Bp0 * v_R) + &
                                               Bp0_AR * (v * Bp0_p) / R + &
                                               BR0 * (v * grad_Bp_AR__R + Bp0_AR * v_R) + &
                                               BZ0 * (v * grad_Bp_AR__Z + Bp0_AR * v_Z)
                  B_dot_grad_BpVstar__AR__n =  BZ0_AR__n * (v * Bp0_Z + Bp0 * v_Z) + &
                                               Bp0 * (v * grad_Bp_AR__phi__n) / R
                  B_dot_grad_BpVstar__AR__k =  Bp0_AR * (Bp0 * v_p__k) / R + &
                                               Bp0 * (Bp0_AR * v_p__k) / R
                  B_dot_grad_BpVstar__AR__kn = 0.d0
                  B_dot_grad_BpVstar__AR__nn = 0.d0

                  B_dot_grad_BRVstar__AZ__p =  BZ0_AZ * (v * BR0_Z + BR0 * v_Z) + &
                                               Bp0_AZ * (v * BR0_p) / R
                  B_dot_grad_BRVstar__AZ__n =  BR0_AZ__n * (v * BR0_R + BR0 * v_R) + &
                                               BR0 * (v * grad_BR_AZ__R__n + BR0_AZ__n * v_R) + &
                                               BZ0 * (v * grad_BR_AZ__Z__n + BR0_AZ__n * v_Z)
                  B_dot_grad_BRVstar__AZ__k =  Bp0_AZ * (BR0 * v_p__k) / R 
                  B_dot_grad_BRVstar__AZ__kn = Bp0 * (BR0_AZ__n * v_p__k) / R
                  B_dot_grad_BRVstar__AZ__nn = Bp0 * (v * grad_BR_AZ__phi__nn) / R

                  B_dot_grad_BZVstar__AZ__p =  BZ0_AZ * (v * BZ0_Z + BZ0 * v_Z) + &
                                               Bp0_AZ * (v * BZ0_p) / R + &
                                               BR0 * (v * grad_BZ_AZ__R + BZ0_AZ * v_R) + &
                                               BZ0 * (v * grad_BZ_AZ__Z + BZ0_AZ * v_Z) + &
                                               Bp0 * (v * grad_BZ_AZ__phi) / R
                  B_dot_grad_BZVstar__AZ__n =  BR0_AZ__n * (v * BZ0_R + BZ0 * v_R)
                  B_dot_grad_BZVstar__AZ__k =  Bp0_AZ * (BZ0 * v_p__k) / R + &
                                               Bp0 * (BZ0_AZ * v_p__k) / R
                  B_dot_grad_BZVstar__AZ__kn = 0.d0
                  B_dot_grad_BZVstar__AZ__nn = 0.d0

                  B_dot_grad_BpVstar__AZ__p =  BZ0_AZ * (v * Bp0_Z + Bp0 * v_Z) + &
                                               Bp0_AZ * (v * Bp0_p) / R + &
                                               BR0 * (v * grad_Bp_AZ__R + Bp0_AZ * v_R) + &
                                               BZ0 * (v * grad_Bp_AZ__Z + Bp0_AZ * v_Z)
                  B_dot_grad_BpVstar__AZ__n =  BR0_AZ__n * (v * Bp0_R + Bp0 * v_R) + &
                                               Bp0 * (v * grad_Bp_AZ__phi__n) / R
                  B_dot_grad_BpVstar__AZ__k =  Bp0_AZ * (Bp0 * v_p__k) / R + &
                                               Bp0 * (Bp0_AZ * v_p__k) / R
                  B_dot_grad_BpVstar__AZ__kn = 0.d0
                  B_dot_grad_BpVstar__AZ__nn = 0.d0

                  B_dot_grad_BRVstar__A3__p =  BR0_A3 * (v * BR0_R + BR0 * v_R) + &
                                               BZ0_A3 * (v * BR0_Z + BR0 * v_Z) + &
                                               Bp0_A3 * (v * BR0_p) / R + &
                                               BR0 * (v * grad_BR_A3__R + BR0_A3 * v_R) + &
                                               BZ0 * (v * grad_BR_A3__Z + BR0_A3 * v_Z)
                  B_dot_grad_BRVstar__A3__n =  Bp0 * (v * grad_BR_A3__phi__n) / R
                  B_dot_grad_BRVstar__A3__k =  Bp0_A3 * (BR0 * v_p__k) / R + &
                                               Bp0 * (BR0_A3 * v_p__k) / R
                  B_dot_grad_BRVstar__A3__kn = 0.d0
                  B_dot_grad_BRVstar__A3__nn = 0.d0

                  B_dot_grad_BZVstar__A3__p =  BR0_A3 * (v * BZ0_R + BZ0 * v_R) + &
                                               BZ0_A3 * (v * BZ0_Z + BZ0 * v_Z) + &
                                               Bp0_A3 * (v * BZ0_p) / R + &
                                               BR0 * (v * grad_BZ_A3__R + BZ0_A3 * v_R) + &
                                               BZ0 * (v * grad_BZ_A3__Z + BZ0_A3 * v_Z) 
                  B_dot_grad_BZVstar__A3__n =  Bp0 * (v * grad_BZ_A3__phi__n) / R
                  B_dot_grad_BZVstar__A3__k =  Bp0_A3 * (BZ0 * v_p__k) / R + &
                                               Bp0 * (BZ0_A3 * v_p__k) / R
                  B_dot_grad_BZVstar__A3__kn = 0.d0
                  B_dot_grad_BZVstar__A3__nn = 0.d0

                  B_dot_grad_BpVstar__A3__p =  BR0_A3 * (v * Bp0_R + Bp0 * v_R) + &
                                               BZ0_A3 * (v * Bp0_Z + Bp0 * v_Z) + &
                                               Bp0_A3 * (v * Bp0_p) / R + &
                                               BR0 * (v * grad_Bp_A3__R + Bp0_A3 * v_R) + &
                                               BZ0 * (v * grad_Bp_A3__Z + Bp0_A3 * v_Z) + &
                                               Bp0 * (v * grad_Bp_A3__phi) / R
                  B_dot_grad_BpVstar__A3__n =  0.d0
                  B_dot_grad_BpVstar__A3__k =  Bp0_A3 * (Bp0 * v_p__k) / R + &
                                               Bp0 * (Bp0_A3 * v_p__k) / R
                  B_dot_grad_BpVstar__A3__kn = 0.d0
                  B_dot_grad_BpVstar__A3__nn = 0.d0

                  B_dot_grad_BB2__AR__p  = BR0_AR * BB2_R + Bp0_AR * BB2_phi / R + &
                                           BR0 * (grad_BB2_AR__R__p) + &
                                           BZ0 * (grad_BB2_AR__Z__p) + &
                                           Bp0 * (grad_BB2_AR__phi__p) / R
                  B_dot_grad_BB2__AR__n  = BZ0_AR__n * BB2_Z + &
                                           BR0 * (grad_BB2_AR__R__n) + &
                                           BZ0 * (grad_BB2_AR__Z__n) + &
                                           Bp0 * (grad_BB2_AR__phi__n) / R
                  B_dot_grad_BB2__AR__nn = BR0 * (grad_BB2_AR__R__nn) + &
                                           BZ0 * (grad_BB2_AR__Z__nn) + &
                                           Bp0 * (grad_BB2_AR__phi__nn) / R

                  B_dot_grad_BB2__AZ__p  = BZ0_AZ * BB2_Z + Bp0_AZ * BB2_phi / R + &
                                           BR0 * (grad_BB2_AZ__R__p) + &
                                           BZ0 * (grad_BB2_AZ__Z__p) + &
                                           Bp0 * (grad_BB2_AZ__phi__p) / R
                  B_dot_grad_BB2__AZ__n  = BR0_AZ__n * BB2_R + &
                                           BR0 * (grad_BB2_AZ__R__n) + &
                                           BZ0 * (grad_BB2_AZ__Z__n) + &
                                           Bp0 * (grad_BB2_AZ__phi__n) / R
                  B_dot_grad_BB2__AZ__nn = BR0 * (grad_BB2_AZ__R__nn) + &
                                           BZ0 * (grad_BB2_AZ__Z__nn) + &
                                           Bp0 * (grad_BB2_AZ__phi__nn) / R

                  B_dot_grad_BB2__A3__p  = BR0_A3 * BB2_R + BZ0_A3 * BB2_Z + Bp0_A3 * BB2_phi / R + &
                                           BR0 * (grad_BB2_A3__R__p) + &
                                           BZ0 * (grad_BB2_A3__Z__p) + &
                                           Bp0 * (grad_BB2_A3__phi__p) / R
                  B_dot_grad_BB2__A3__n  = BR0 * (grad_BB2_A3__R__n) + &
                                           BZ0 * (grad_BB2_A3__Z__n) + &
                                           Bp0 * (grad_BB2_A3__phi__n) / R
                  B_dot_grad_BB2__A3__nn = BR0 * (grad_BB2_A3__R__nn) + &
                                           BZ0 * (grad_BB2_A3__Z__nn) + &
                                           Bp0 * (grad_BB2_A3__phi__nn) / R

                  B_AR_dot_grad_BdotU0__p = BR0_AR * grad_BdotU0_R + Bp0_AR * grad_BdotU0_phi / R
                  B_AR_dot_grad_BdotU0__n = BZ0_AR__n * grad_BdotU0_Z
                  B_AZ_dot_grad_BdotU0__p = BZ0_AZ * grad_BdotU0_Z + Bp0_AZ * grad_BdotU0_phi / R
                  B_AZ_dot_grad_BdotU0__n = BR0_AZ__n * grad_BdotU0_R
                  B_A3_dot_grad_BdotU0__p = BR0_A3 * grad_BdotU0_R + BZ0_A3 * grad_BdotU0_Z + Bp0_A3 * grad_BdotU0_phi / R
                  B_A3_dot_grad_BdotU0__n = 0.d0

                  B_dot_grad__BdotU0_AR__p  = BR0 * (UR0_R * BR0_AR + UR0 * grad_BR_AR__R   + Up0_R * Bp0_AR + Up0 * grad_Bp_AR__R) + &
                                              BZ0 * (UR0_Z * BR0_AR + UR0 * grad_BR_AR__Z   + Up0_Z * Bp0_AR + Up0 * grad_Bp_AR__Z) + &
                                              Bp0 * (UR0_p * BR0_AR + UR0 * grad_BR_AR__phi + Up0_p * Bp0_AR) / R
                  B_dot_grad__BdotU0_AR__n  = BR0 * (UZ0_R * BZ0_AR__n + UZ0 * grad_BZ_AR__R__n) + &
                                              BZ0 * (UZ0_Z * BZ0_AR__n + UZ0 * grad_BZ_AR__Z__n) + &
                                              Bp0 * (UZ0_p * BZ0_AR__n + Up0 * grad_Bp_AR__phi__n) / R
                  B_dot_grad__BdotU0_AR__nn = Bp0 * (UZ0 * grad_BZ_AR__phi__nn) / R

                  B_dot_grad__BdotU0_AZ__p  = BR0 * (UZ0_R * BZ0_AZ + UZ0 * grad_BZ_AZ__R    + Up0_R * Bp0_AZ + Up0 * grad_Bp_AZ__R) + &
                                              BZ0 * (UZ0_Z * BZ0_AZ + UZ0 * grad_BZ_AZ__Z    + Up0_Z * Bp0_AZ + Up0 * grad_Bp_AZ__Z) + &
                                              Bp0 * (UZ0_p * BZ0_AZ + UZ0 * grad_BZ_AZ__phi  + Up0_p * Bp0_AZ) / R
                  B_dot_grad__BdotU0_AZ__n  = BR0 * (UR0_R * BR0_AZ__n + UR0 * grad_BR_AZ__R__n) + &
                                              BZ0 * (UR0_Z * BR0_AZ__n + UR0 * grad_BR_AZ__Z__n) + &
                                              Bp0 * (UR0_p * BR0_AZ__n + Up0 * grad_Bp_AZ__phi__n) / R
                  B_dot_grad__BdotU0_AZ__nn = Bp0 * (UR0 * grad_BR_AZ__phi__nn) / R

                  B_dot_grad__BdotU0_A3__p  = BR0 * (UR0_R * BR0_A3 + UR0 * grad_BR_A3__R      + UZ0_R * BZ0_A3 + UZ0 * grad_BZ_A3__R      + Up0_R * Bp0_A3 + Up0 * grad_Bp_A3__R) + &
                                              BZ0 * (UR0_Z * BR0_A3 + UR0 * grad_BR_A3__Z      + UZ0_Z * BZ0_A3 + UZ0 * grad_BZ_A3__Z      + Up0_Z * Bp0_A3 + Up0 * grad_Bp_A3__Z) + &
                                              Bp0 * (UR0_p * BR0_A3 + UZ0_p * BZ0_A3 + Up0_p * Bp0_A3 + Up0 * grad_Bp_A3__phi) / R
                  B_dot_grad__BdotU0_A3__n  = Bp0 * (UR0 * grad_BR_A3__phi__n + UZ0 * grad_BZ_A3__phi__n) / R
                  B_dot_grad__BdotU0_A3__nn = 0.d0

                  B_dot_grad__UR__p = BR0 *  UR_R + BZ0 * UR_Z; B_dot_grad__UR__n = Bp0 * UR_p__n / R
                  B_dot_grad__UZ__p = BR0 *  UZ_R + BZ0 * UZ_Z; B_dot_grad__UZ__n = Bp0 * UZ_p__n / R
                  B_dot_grad__Up__p = BR0 *  Up_R + BZ0 * Up_Z; B_dot_grad__Up__n = Bp0 * Up_p__n / R

                  ! B_dot_grad_BR[Z/p]__AR[Z/3]
                  B_dot_grad_BR__AR__p =  BR0_AR    * BR0_R + BR0 * grad_BR_AR__R &
                                       +                      BZ0 * grad_BR_AR__Z &
                                       + (Bp0_AR    * BR0_p + Bp0 * grad_BR_AR__phi) / R
                  B_dot_grad_BR__AR__n = BZ0_AR__n * BR0_Z

                  B_dot_grad_BZ__AR__p =  BR0_AR    * BZ0_R                          &
                                       + (Bp0_AR    * BZ0_p                            ) / R
                  B_dot_grad_BZ__AR__n =                      BR0 * grad_BZ_AR__R__n &
                                       +  BZ0_AR__n * BZ0_Z + BZ0 * grad_BZ_AR__Z__n
                  B_dot_grad_BZ__AR__nn= (                    Bp0 * grad_BZ_AR__phi__nn) / R

                  B_dot_grad_Bp__AR__p =  BR0_AR    * Bp0_R + BR0 * grad_Bp_AR__R &
                                       +                      BZ0 * grad_Bp_AR__Z &
                                       + (Bp0_AR    * Bp0_p                           ) / R
                  B_dot_grad_Bp__AR__n =  BZ0_AR__n * Bp0_Z                       &
                                       + (                    Bp0 * grad_Bp_AR__phi__n) / R

                  B_dot_grad_BR__AZ__p =  BZ0_AZ    * BR0_Z                          &
                                       + (Bp0_AZ    * BR0_p                            ) / R
                  B_dot_grad_BR__AZ__n =  BR0_AZ__n * BR0_R + BR0 * grad_BR_AZ__R__n &
                                       +                      BZ0 * grad_BR_AZ__Z__n 

                  B_dot_grad_BR__AZ__nn= (                    Bp0 * grad_BR_AZ__phi__nn) / R

                  B_dot_grad_BZ__AZ__p =                      BR0 * grad_BZ_AZ__R &
                                       +  BZ0_AZ    * BZ0_Z + BZ0 * grad_BZ_AZ__Z &
                                       + (Bp0_AZ    * BZ0_p + Bp0 * grad_BZ_AZ__phi) / R
                  B_dot_grad_BZ__AZ__n =  BR0_AZ__n * BZ0_R

                  B_dot_grad_Bp__AZ__p =                      BR0 * grad_Bp_AZ__R &
                                       +  BZ0_AZ    * Bp0_Z + BZ0 * grad_Bp_AZ__Z &
                                       + (Bp0_AZ    * Bp0_p                           ) / R
                  B_dot_grad_Bp__AZ__n =  BR0_AZ__n * Bp0_R                       &
                                       + (                    Bp0 * grad_Bp_AZ__phi__n) / R


                  B_dot_grad_BR__A3__p =  BR0_A3    * BR0_R + BR0 * grad_BR_A3__R &
                                       +  BZ0_A3    * BR0_Z + BZ0 * grad_BR_A3__Z &
                                       + (Bp0_A3    * BR0_p                           ) / R
                  B_dot_grad_BR__A3__n = (                    Bp0 * grad_BR_A3__phi__n) / R

                  B_dot_grad_BZ__A3__p =  BR0_A3    * BZ0_R + BR0 * grad_BZ_A3__R &
                                       +  BZ0_A3    * BZ0_Z + BZ0 * grad_BZ_A3__Z &
                                       + (Bp0_A3    * BZ0_p                           ) / R
                  B_dot_grad_BZ__A3__n = (                    Bp0 * grad_BZ_A3__phi__n) / R
                  B_dot_grad_Bp__A3__p =  BR0_A3    * Bp0_R + BR0 * grad_Bp_A3__R &
                                       +  BZ0_A3    * Bp0_Z + BZ0 * grad_Bp_A3__Z &
                                       + (Bp0_A3    * Bp0_p + Bp0 * grad_Bp_A3__phi) / R

                  ! --- var_UR,var_AR
                  Qjac_UR_visco_par_AR_p = ( &
                    -1.d0/BB2**2*(BB2_AR__p)*( &
                      (grad_BdotU0_R - 1.d0/BB2*B_dot_U0*BB2_R)*grad_BRVstar_R  &
                    + (grad_BdotU0_Z - 1.d0/BB2*B_dot_U0*BB2_Z)*grad_BRVstar_Z  &
                    + (grad_BdotU0_phi - 1.d0/BB2*B_dot_U0*BB2_phi)*(grad_BRVstar_phi)/R**2)  &
                    + 1.d0/BB2*( &
                      ((grad_B_AR_dot_U__R__p)+1.d0/BB2**2*(BB2_AR__p)*B_dot_U0*BB2_R &
                    - 1.d0/BB2*(B_AR_dot_U__p)*BB2_R-1.d0/BB2*B_dot_U0*(grad_BB2_AR__R__p))*grad_BRVstar_R  &
                    + (grad_BdotU0_R-1.d0/BB2*B_dot_U0*BB2_R)*(v*grad_BR_AR__R+BR0_AR*v_R) &
                    + ((grad_B_AR_dot_U__Z__p)+1.d0/BB2**2*(BB2_AR__p)*B_dot_U0*BB2_Z &
                    - 1.d0/BB2*(B_AR_dot_U__p)*BB2_Z-1.d0/BB2*B_dot_U0*(grad_BB2_AR__Z__p))*grad_BRVstar_Z  &
                    + (grad_BdotU0_Z-1.d0/BB2*B_dot_U0*BB2_Z)*(v*grad_BR_AR__Z+BR0_AR*v_Z) &
                    + ((grad_B_AR_dot_U__phi__p)+1.d0/BB2**2*(BB2_AR__p)*B_dot_U0*BB2_phi  &
                    - 1.d0/BB2*(B_AR_dot_U__p)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_AR__phi__p))*(grad_BRVstar_phi)/R**2  &
                    + (grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*(v*grad_BR_AR__phi)/R**2 ) &
                    )
                  Qjac_UR_visco_parpar_AR_p = ( &
                    -2.d0/BB2**3*(BB2_AR__p)*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BRVstar__p) &
                    +1.d0/BB2**2*(( &
                    1.d0/BB2**2*(BB2_AR__p)*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__AR__p)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_AR_dot_U__p)+(B_AR_dot_grad_BdotU0__p)+(B_dot_grad__BdotU0_AR__p) &
                    )*(B_dot_grad_BRVstar__p) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BRVstar__AR__p)) &
                    )
                  Qjac_UR_visco_par_AR_n = ( &
                    -1.d0/BB2**2*(BB2_AR__n)*( &
                      (grad_BdotU0_R - 1.d0/BB2*B_dot_U0*BB2_R)*grad_BRVstar_R  &
                    + (grad_BdotU0_Z - 1.d0/BB2*B_dot_U0*BB2_Z)*grad_BRVstar_Z  &
                    + (grad_BdotU0_phi - 1.d0/BB2*B_dot_U0*BB2_phi)*(grad_BRVstar_phi)/R**2)  &
                    + 1.d0/BB2*( &
                      ((grad_B_AR_dot_U__R__n)+1.d0/BB2**2*(BB2_AR__n)*B_dot_U0*BB2_R &
                    - 1.d0/BB2*(B_AR_dot_U__n)*BB2_R-1.d0/BB2*B_dot_U0*(grad_BB2_AR__R__n))*grad_BRVstar_R  &
                    + ((grad_B_AR_dot_U__Z__n)+1.d0/BB2**2*(BB2_AR__n)*B_dot_U0*BB2_Z &
                    - 1.d0/BB2*(B_AR_dot_U__n)*BB2_Z-1.d0/BB2*B_dot_U0*(grad_BB2_AR__Z__n))*grad_BRVstar_Z  &
                    + ((grad_B_AR_dot_U__phi__n)+1.d0/BB2**2*(BB2_AR__n)*B_dot_U0*BB2_phi  &
                    - 1.d0/BB2*(B_AR_dot_U__n)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_AR__phi__n))*(grad_BRVstar_phi)/R**2  &
                    ) &
                    )
                  Qjac_UR_visco_parpar_AR_n = ( &
                    -2.d0/BB2**3*(BB2_AR__n)*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BRVstar__p) &
                    +1.d0/BB2**2*(( &
                    1.d0/BB2**2*(BB2_AR__n)*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__AR__n)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_AR_dot_U__n)+(B_AR_dot_grad_BdotU0__n)+(B_dot_grad__BdotU0_AR__n) &
                    )*(B_dot_grad_BRVstar__p) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BRVstar__AR__n)) &
                    )
                  Qjac_UR_visco_par_AR_k = ( &
                    -1.d0/BB2**2*(BB2_AR__p)*( &
                    + (grad_BdotU0_phi - 1.d0/BB2*B_dot_U0*BB2_phi)*(grad_BRVstar_phi_k)/R**2)  &
                    + 1.d0/BB2*( &
                    + ((grad_B_AR_dot_U__phi__p)+1.d0/BB2**2*(BB2_AR__p)*B_dot_U0*BB2_phi  &
                    - 1.d0/BB2*(B_AR_dot_U__p)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_AR__phi__p))*(grad_BRVstar_phi_k)/R**2  &
                    + (grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*(BR0_AR*v_p__k)/R**2 ) &
                    )
                  Qjac_UR_visco_parpar_AR_k = ( &
                    -2.d0/BB2**3*(BB2_AR__p)*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BRVstar__k) &
                    +1.d0/BB2**2*(( &
                    1.d0/BB2**2*(BB2_AR__p)*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__AR__p)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_AR_dot_U__p)+(B_AR_dot_grad_BdotU0__p)+(B_dot_grad__BdotU0_AR__p) &
                    )*(B_dot_grad_BRVstar__k) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BRVstar__AR__k)) &
                    )
                  Qjac_UR_visco_par_AR_kn = ( &
                    -1.d0/BB2**2*(BB2_AR__n)*( &
                    + (grad_BdotU0_phi - 1.d0/BB2*B_dot_U0*BB2_phi)*(grad_BRVstar_phi_k)/R**2)  &
                    + 1.d0/BB2*( &
                    + ((grad_B_AR_dot_U__phi__n)+1.d0/BB2**2*(BB2_AR__n)*B_dot_U0*BB2_phi  &
                    - 1.d0/BB2*(B_AR_dot_U__n)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_AR__phi__n))*(grad_BRVstar_phi_k)/R**2  &
                    ) &
                    )
                  Qjac_UR_visco_parpar_AR_kn = ( &
                    -2.d0/BB2**3*(BB2_AR__n)*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BRVstar__k) &
                    +1.d0/BB2**2*(( &
                    1.d0/BB2**2*(BB2_AR__n)*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__AR__n)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_AR_dot_U__n)+(B_AR_dot_grad_BdotU0__n)+(B_dot_grad__BdotU0_AR__n) &
                    )*(B_dot_grad_BRVstar__k) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BRVstar__AR__kn)) &
                    )
                  Qjac_UR_visco_par_AR_pnn = ( &
                    + 1.d0/BB2*( &
                      ( &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AR__R__nn))*grad_BRVstar_R  &
                    + ( &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AR__Z__nn))*grad_BRVstar_Z  &
                    + ((grad_B_AR_dot_U__phi__nn)  &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AR__phi__nn))*(grad_BRVstar_phi)/R**2  &
                    ) &
                    )
                  Qjac_UR_visco_parpar_AR_pnn = ( &
                    +1.d0/BB2**2*(( &
                    -1.d0/BB2*(B_dot_grad_BB2__AR__nn)*B_dot_U0 &
                    +(B_dot_grad__BdotU0_AR__nn) &
                    )*(B_dot_grad_BRVstar__p) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BRVstar__AR__nn)) &
                    )
                  Qjac_UR_visco_par_AR_knn = ( &
                    + 1.d0/BB2*( &
                    + ((grad_B_AR_dot_U__phi__nn)  &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AR__phi__nn))*(grad_BRVstar_phi_k)/R**2  &
                    ) &
                    )
                  Qjac_UR_visco_parpar_AR_knn = ( &
                    +1.d0/BB2**2*(( &
                    -1.d0/BB2*(B_dot_grad_BB2__AR__nn)*B_dot_U0 &
                    +(B_dot_grad__BdotU0_AR__nn) &
                    )*(B_dot_grad_BRVstar__k) &
                    ) &
                    )
                  
                  Qjac_UR_visco_par_AR_for_Uphi_p = &
                    - v * 1.d0/BB2**2*(BB2_AR__p)*( &
                      (grad_BdotU0_R - 1.d0/BB2*B_dot_U0*BB2_R)*BR0_R  &
                    + (grad_BdotU0_Z - 1.d0/BB2*B_dot_U0*BB2_Z)*BR0_Z  &
                    + (grad_BdotU0_phi - 1.d0/BB2*B_dot_U0*BB2_phi)*BR0_p/R**2)  &
                    + v * 1.d0/BB2*( &
                      ((grad_B_AR_dot_U__R__p)+1.d0/BB2**2*(BB2_AR__p)*B_dot_U0*BB2_R &
                    - 1.d0/BB2*(B_AR_dot_U__p)*BB2_R-1.d0/BB2*B_dot_U0*(grad_BB2_AR__R__p))*BR0_R  &
                    + (grad_BdotU0_R-1.d0/BB2*B_dot_U0*BB2_R)*grad_BR_AR__R &
                    + ((grad_B_AR_dot_U__Z__p)+1.d0/BB2**2*(BB2_AR__p)*B_dot_U0*BB2_Z &
                    - 1.d0/BB2*(B_AR_dot_U__p)*BB2_Z-1.d0/BB2*B_dot_U0*(grad_BB2_AR__Z__p))*BR0_Z  &
                    + (grad_BdotU0_Z-1.d0/BB2*B_dot_U0*BB2_Z)*grad_BR_AR__Z &
                    + ((grad_B_AR_dot_U__phi__p)+1.d0/BB2**2*(BB2_AR__p)*B_dot_U0*BB2_phi  &
                    - 1.d0/BB2*(B_AR_dot_U__p)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_AR__phi__p))*BR0_p/R**2  &
                    + (grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*grad_BR_AR__phi/R**2 )

                  Qjac_UR_visco_par_AR_for_Uphi_n = &
                    - v * 1.d0/BB2**2*(BB2_AR__n)*( &
                      (grad_BdotU0_R - 1.d0/BB2*B_dot_U0*BB2_R)*BR0_R  &
                    + (grad_BdotU0_Z - 1.d0/BB2*B_dot_U0*BB2_Z)*BR0_Z  &
                    + (grad_BdotU0_phi - 1.d0/BB2*B_dot_U0*BB2_phi)*BR0_p/R**2)  &
                    + v * 1.d0/BB2*( &
                      ((grad_B_AR_dot_U__R__n)+1.d0/BB2**2*(BB2_AR__n)*B_dot_U0*BB2_R &
                    - 1.d0/BB2*(B_AR_dot_U__n)*BB2_R-1.d0/BB2*B_dot_U0*(grad_BB2_AR__R__n))*BR0_R  &
                    + ((grad_B_AR_dot_U__Z__n)+1.d0/BB2**2*(BB2_AR__n)*B_dot_U0*BB2_Z &
                    - 1.d0/BB2*(B_AR_dot_U__n)*BB2_Z-1.d0/BB2*B_dot_U0*(grad_BB2_AR__Z__n))*BR0_Z  &
                    + ((grad_B_AR_dot_U__phi__n)+1.d0/BB2**2*(BB2_AR__n)*B_dot_U0*BB2_phi  &
                    - 1.d0/BB2*(B_AR_dot_U__n)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_AR__phi__n))*BR0_p/R**2  &
                    )
                  Qjac_UR_visco_par_AR_for_Uphi_pnn = &
                    + v * 1.d0/BB2*( &
                      ( &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AR__R__nn))*BR0_R  &
                    + ( &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AR__Z__nn))*BR0_Z  &
                    + ((grad_B_AR_dot_U__phi__nn)  &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AR__phi__nn))*BR0_p/R**2  &
                    )
                  
                  Qjac_UR_visco_parpar_AR_for_Uphi_p = &
                    - v * 2.d0/BB2**3*(BB2_AR__p)*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*B_dot_grad_BR &
                    + v * 1.d0/BB2**2*(( &
                    1.d0/BB2**2*(BB2_AR__p)*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__AR__p)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_AR_dot_U__p)+(B_AR_dot_grad_BdotU0__p)+(B_dot_grad__BdotU0_AR__p) &
                    )*B_dot_grad_BR &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BR__AR__p))

                  Qjac_UR_visco_parpar_AR_for_Uphi_n = &
                    - v * 2.d0/BB2**3*(BB2_AR__n)*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*B_dot_grad_BR &
                    + v * 1.d0/BB2**2*(( &
                    1.d0/BB2**2*(BB2_AR__n)*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__AR__n)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_AR_dot_U__n)+(B_AR_dot_grad_BdotU0__n)+(B_dot_grad__BdotU0_AR__n) &
                    )*B_dot_grad_BR &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BR__AR__n))

                  Qjac_UR_visco_parpar_AR_for_Uphi_pnn = &
                    + v * 1.d0/BB2**2*(( &
                    -1.d0/BB2*(B_dot_grad_BB2__AR__nn)*B_dot_U0 &
                    +(B_dot_grad__BdotU0_AR__nn) &
                    )*B_dot_grad_BR &
                    )

                  ! --- var_UR,var_AZ
                  Qjac_UR_visco_par_AZ_p = ( &
                    -1.d0/BB2**2*(BB2_AZ__p)*( &
                    (grad_BdotU0_R-1.d0/BB2*B_dot_U0*BB2_R)*grad_BRVstar_R &
                    +(grad_BdotU0_Z-1.d0/BB2*B_dot_U0*BB2_Z)*grad_BRVstar_Z &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*(grad_BRVstar_phi)/R**2) &
                    +1.d0/BB2*( &
                    ((grad_B_AZ_dot_U__R__p) + 1.d0/BB2**2*(BB2_AZ__p)*B_dot_U0*BB2_R &
                    -1.d0/BB2*(B_AZ_dot_U__p)*BB2_R-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__R__p))*grad_BRVstar_R &
                    +((grad_B_AZ_dot_U__Z__p)+1.d0/BB2**2*(BB2_AZ__p)*B_dot_U0*BB2_Z &
                    -1.d0/BB2*(B_AZ_dot_U__p)*BB2_Z-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__Z__p))*grad_BRVstar_Z &
                    +((grad_B_AZ_dot_U__phi__p)+1.d0/BB2**2*(BB2_AZ__p)*B_dot_U0*BB2_phi &
                    -1.d0/BB2*(B_AZ_dot_U__p)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__phi__p))*(grad_BRVstar_phi)/R**2 &
                    ) &
                    )
                  Qjac_UR_visco_parpar_AZ_p = ( &
                    -2.d0/BB2**3*(BB2_AZ__p)*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BRVstar__p) &
                    +1.d0/BB2**2*(( &
                    1.d0/BB2**2*(BB2_AZ__p)*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__AZ__p)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_AZ_dot_U__p)+(B_AZ_dot_grad_BdotU0__p)+(B_dot_grad__BdotU0_AZ__p) &
                    )*(B_dot_grad_BRVstar__p) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BRVstar__AZ__p)) &
                    )
                  Qjac_UR_visco_par_AZ_n = ( &
                    -1.d0/BB2**2*(BB2_AZ__n)*( &
                    (grad_BdotU0_R-1.d0/BB2*B_dot_U0*BB2_R)*grad_BRVstar_R &
                    +(grad_BdotU0_Z-1.d0/BB2*B_dot_U0*BB2_Z)*grad_BRVstar_Z &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*(grad_BRVstar_phi)/R**2) &
                    +1.d0/BB2*( &
                    ((grad_B_AZ_dot_U__R__n) + 1.d0/BB2**2*(BB2_AZ__n)*B_dot_U0*BB2_R &
                    -1.d0/BB2*(B_AZ_dot_U__n)*BB2_R-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__R__n))*grad_BRVstar_R &
                    +(grad_BdotU0_R-1.d0/BB2*B_dot_U0*BB2_R)*(v*grad_BR_AZ__R__n+BR0_AZ__n*v_R) &
                    +((grad_B_AZ_dot_U__Z__n)+1.d0/BB2**2*(BB2_AZ__n)*B_dot_U0*BB2_Z &
                    -1.d0/BB2*(B_AZ_dot_U__n)*BB2_Z-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__Z__n))*grad_BRVstar_Z &
                    +(grad_BdotU0_Z-1.d0/BB2*B_dot_U0*BB2_Z)*(v*grad_BR_AZ__Z__n+BR0_AZ__n*v_Z) &
                    +((grad_B_AZ_dot_U__phi__n)+1.d0/BB2**2*(BB2_AZ__n)*B_dot_U0*BB2_phi &
                    -1.d0/BB2*(B_AZ_dot_U__n)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__phi__n))*(grad_BRVstar_phi)/R**2 &
                    ) &
                    )
                  Qjac_UR_visco_parpar_AZ_n = ( &
                    -2.d0/BB2**3*(BB2_AZ__n)*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BRVstar__p) &
                    +1.d0/BB2**2*(( &
                    1.d0/BB2**2*(BB2_AZ__n)*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__AZ__n)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_AZ_dot_U__n)+(B_AZ_dot_grad_BdotU0__n)+(B_dot_grad__BdotU0_AZ__n) &
                    )*(B_dot_grad_BRVstar__p) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BRVstar__AZ__n)) &
                    )
                    Qjac_UR_visco_par_AZ_k = ( &
                    -1.d0/BB2**2*(BB2_AZ__p)*( &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*(grad_BRVstar_phi_k)/R**2) &
                    +1.d0/BB2*( &
                    +((grad_B_AZ_dot_U__phi__p)+1.d0/BB2**2*(BB2_AZ__p)*B_dot_U0*BB2_phi &
                    -1.d0/BB2*(B_AZ_dot_U__p)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__phi__p))*(grad_BRVstar_phi_k)/R**2 &
                    ) &
                    )
                  Qjac_UR_visco_parpar_AZ_k = ( &
                    -2.d0/BB2**3*(BB2_AZ__p)*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BRVstar__k) &
                    +1.d0/BB2**2*(( &
                    1.d0/BB2**2*(BB2_AZ__p)*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__AZ__p)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_AZ_dot_U__p)+(B_AZ_dot_grad_BdotU0__p)+(B_dot_grad__BdotU0_AZ__p) &
                    )*(B_dot_grad_BRVstar__k) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BRVstar__AZ__k)) &
                    )
                  Qjac_UR_visco_par_AZ_kn = ( &
                    -1.d0/BB2**2*(BB2_AZ__n)*( &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*(grad_BRVstar_phi_k)/R**2) &
                    +1.d0/BB2*( &
                    +((grad_B_AZ_dot_U__phi__n)+1.d0/BB2**2*(BB2_AZ__n)*B_dot_U0*BB2_phi &
                    -1.d0/BB2*(B_AZ_dot_U__n)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__phi__n))*(grad_BRVstar_phi_k)/R**2 &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*(BR0_AZ__n*v_p__k)/R**2) &
                    )
                  Qjac_UR_visco_parpar_AZ_kn = ( &
                    -2.d0/BB2**3*(BB2_AZ__n)*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BRVstar__k) &
                    +1.d0/BB2**2*(( &
                    1.d0/BB2**2*(BB2_AZ__n)*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__AZ__n)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_AZ_dot_U__n)+(B_AZ_dot_grad_BdotU0__n)+(B_dot_grad__BdotU0_AZ__n) &
                    )*(B_dot_grad_BRVstar__k) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BRVstar__AZ__kn)) &
                    )
                  Qjac_UR_visco_par_AZ_pnn = ( &
                    +1.d0/BB2*( &
                    ( &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AZ__R__nn))*grad_BRVstar_R &
                    +( &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AZ__Z__nn))*grad_BRVstar_Z &
                    +((grad_B_AZ_dot_U__phi__nn) &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AZ__phi__nn))*(grad_BRVstar_phi)/R**2 &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*(v*grad_BR_AZ__phi__nn)/R**2) &
                    )
                  Qjac_UR_visco_parpar_AZ_pnn = ( &
                    +1.d0/BB2**2*(( &
                    -1.d0/BB2*(B_dot_grad_BB2__AZ__nn)*B_dot_U0 &
                    +(B_dot_grad__BdotU0_AZ__nn) &
                    )*(B_dot_grad_BRVstar__p) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BRVstar__AZ__nn)) &
                    )
                  Qjac_UR_visco_par_AZ_knn = ( &
                    +1.d0/BB2*( &
                    +((grad_B_AZ_dot_U__phi__nn) &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AZ__phi__nn))*(grad_BRVstar_phi_k)/R**2 &
                    ) &
                    )
                  Qjac_UR_visco_parpar_AZ_knn = ( &
                    +1.d0/BB2**2*(( &
                    -1.d0/BB2*(B_dot_grad_BB2__AZ__nn)*B_dot_U0 &
                    +(B_dot_grad__BdotU0_AZ__nn) &
                    )*(B_dot_grad_BRVstar__k) &
                    ) &
                    )

                  Qjac_UR_visco_par_AZ_for_Uphi_p = &
                    - v * 1.d0/BB2**2*(BB2_AZ__p)*( &
                     (grad_BdotU0_R-1.d0/BB2*B_dot_U0*BB2_R)*BR0_R &
                    +(grad_BdotU0_Z-1.d0/BB2*B_dot_U0*BB2_Z)*BR0_Z &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*BR0_p/R**2) &
                    + v * 1.d0/BB2*( &
                    ((grad_B_AZ_dot_U__R__p) + 1.d0/BB2**2*(BB2_AZ__p)*B_dot_U0*BB2_R &
                    -1.d0/BB2*(B_AZ_dot_U__p)*BB2_R-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__R__p))*BR0_R &
                    +((grad_B_AZ_dot_U__Z__p)+1.d0/BB2**2*(BB2_AZ__p)*B_dot_U0*BB2_Z &
                    -1.d0/BB2*(B_AZ_dot_U__p)*BB2_Z-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__Z__p))*BR0_Z &
                    +((grad_B_AZ_dot_U__phi__p)+1.d0/BB2**2*(BB2_AZ__p)*B_dot_U0*BB2_phi &
                    -1.d0/BB2*(B_AZ_dot_U__p)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__phi__p))*BR0_p/R**2 &
                    ) 

                  Qjac_UR_visco_par_AZ_for_Uphi_n = &
                    - v * 1.d0/BB2**2*(BB2_AZ__n)*( &
                     (grad_BdotU0_R-1.d0/BB2*B_dot_U0*BB2_R)*BR0_R &
                    +(grad_BdotU0_Z-1.d0/BB2*B_dot_U0*BB2_Z)*BR0_Z &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*BR0_p/R**2) &
                    + v * 1.d0/BB2*( &
                    ((grad_B_AZ_dot_U__R__n) + 1.d0/BB2**2*(BB2_AZ__n)*B_dot_U0*BB2_R &
                    -1.d0/BB2*(B_AZ_dot_U__n)*BB2_R-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__R__n))*BR0_R &
                    +(grad_BdotU0_R-1.d0/BB2*B_dot_U0*BB2_R)*grad_BR_AZ__R__n &
                    +((grad_B_AZ_dot_U__Z__n)+1.d0/BB2**2*(BB2_AZ__n)*B_dot_U0*BB2_Z &
                    -1.d0/BB2*(B_AZ_dot_U__n)*BB2_Z-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__Z__n))*BR0_Z &
                    +(grad_BdotU0_Z-1.d0/BB2*B_dot_U0*BB2_Z)*grad_BR_AZ__Z__n &
                    +((grad_B_AZ_dot_U__phi__n)+1.d0/BB2**2*(BB2_AZ__n)*B_dot_U0*BB2_phi &
                    -1.d0/BB2*(B_AZ_dot_U__n)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__phi__n))*BR0_p/R**2 &
                    ) 

                  Qjac_UR_visco_par_AZ_for_Uphi_pnn = &
                    + v * 1.d0/BB2*( &
                    ( &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AZ__R__nn))*BR0_R &
                    +( &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AZ__Z__nn))*BR0_Z &
                    +((grad_B_AZ_dot_U__phi__nn) &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AZ__phi__nn))*BR0_p/R**2 &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*grad_BR_AZ__phi__nn/R**2) 

                  Qjac_UR_visco_parpar_AZ_for_Uphi_p = &
                    - v * 2.d0/BB2**3*(BB2_AZ__p)*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*B_dot_grad_BR &
                    + v * 1.d0/BB2**2*(( &
                    1.d0/BB2**2*(BB2_AZ__p)*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__AZ__p)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_AZ_dot_U__p)+(B_AZ_dot_grad_BdotU0__p)+(B_dot_grad__BdotU0_AZ__p) &
                    )*B_dot_grad_BR &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BR__AZ__p)) 

                  Qjac_UR_visco_parpar_AZ_for_Uphi_n = &
                    - v * 2.d0/BB2**3*(BB2_AZ__n)*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*B_dot_grad_BR &
                    + v * 1.d0/BB2**2*(( &
                    1.d0/BB2**2*(BB2_AZ__n)*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__AZ__n)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_AZ_dot_U__n)+(B_AZ_dot_grad_BdotU0__n)+(B_dot_grad__BdotU0_AZ__n) &
                    )*B_dot_grad_BR &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BR__AZ__n)) 

                  Qjac_UR_visco_parpar_AZ_for_Uphi_pnn = &
                    + v * 1.d0/BB2**2*(( &
                    -1.d0/BB2*(B_dot_grad_BB2__AZ__nn)*B_dot_U0 &
                    +(B_dot_grad__BdotU0_AZ__nn) &
                    )*B_dot_grad_BR &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BR__AZ__nn)) 

                  ! --- var_UR,var_A3
                  Qjac_UR_visco_par_A3_p = ( &
                    -1.d0/BB2**2*BB2_A3*( &
                    (grad_BdotU0_R-1.d0/BB2*B_dot_U0*BB2_R) * grad_BRVstar_R &
                    +(grad_BdotU0_Z-1.d0/BB2*B_dot_U0*BB2_Z) * grad_BRVstar_Z &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi) * (grad_BRVstar_phi)/R**2) &
                    +1.d0/BB2*( &
                    ((grad_B_A3_dot_U__R__p)+1.d0/BB2**2*BB2_A3*B_dot_U0*BB2_R &
                    -1.d0/BB2*(B_A3_dot_U__p)*BB2_R-1.d0/BB2*B_dot_U0*(grad_BB2_A3__R__p))*grad_BRVstar_R &
                    +(grad_BdotU0_R-1.d0/BB2*B_dot_U0*BB2_R)*(v*grad_BR_A3__R+BR0_A3*v_R) &
                    +((grad_B_A3_dot_U__Z__p)+1.d0/BB2**2*BB2_A3*B_dot_U0*BB2_Z &
                    -1.d0/BB2*(B_A3_dot_U__p)*BB2_Z-1.d0/BB2*B_dot_U0*(grad_BB2_A3__Z__p))*grad_BRVstar_Z &
                    +(grad_BdotU0_Z-1.d0/BB2*B_dot_U0*BB2_Z)*(v*grad_BR_A3__Z+BR0_A3*v_Z) &
                    +((grad_B_A3_dot_U__phi__p)+1.d0/BB2**2*BB2_A3*B_dot_U0*BB2_phi &
                    -1.d0/BB2*(B_A3_dot_U__p)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_A3__phi__p))*(grad_BRVstar_phi)/R**2 &
                    ) &
                    )
                  Qjac_UR_visco_parpar_A3_p = ( &
                    -2.d0/BB2**3*BB2_A3*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BRVstar__p) &
                    +1.d0/BB2**2*(( &
                    1.d0/BB2**2*BB2_A3*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__A3__p)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_A3_dot_U__p)+(B_A3_dot_grad_BdotU0__p)+(B_dot_grad__BdotU0_A3__p) &
                    )*(B_dot_grad_BRVstar__p) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BRVstar__A3__p)) &
                    )
                  Qjac_UR_visco_par_A3_n = ( &
                    +1.d0/BB2*( &
                    ((grad_B_A3_dot_U__R__n) &
                    -1.d0/BB2*(B_A3_dot_U__n)*BB2_R-1.d0/BB2*B_dot_U0*(grad_BB2_A3__R__n))*grad_BRVstar_R &
                    +((grad_B_A3_dot_U__Z__n) &
                    -1.d0/BB2*(B_A3_dot_U__n)*BB2_Z-1.d0/BB2*B_dot_U0*(grad_BB2_A3__Z__n))*grad_BRVstar_Z &
                    +((grad_B_A3_dot_U__phi__n) &
                    -1.d0/BB2*(B_A3_dot_U__n)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_A3__phi__n))*(grad_BRVstar_phi)/R**2 &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*(v*grad_BR_A3__phi__n)/R**2) &
                    )
                  Qjac_UR_visco_parpar_A3_n = ( &
                    +1.d0/BB2**2*(( &
                    -1.d0/BB2*(B_dot_grad_BB2__A3__n)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_A3_dot_U__n)+(B_A3_dot_grad_BdotU0__n)+(B_dot_grad__BdotU0_A3__n) &
                    )*(B_dot_grad_BRVstar__p) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BRVstar__A3__n)) &
                    )
                  Qjac_UR_visco_par_A3_k = ( &
                    -1.d0/BB2**2*BB2_A3*( &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi) * (grad_BRVstar_phi_k)/R**2) &
                    +1.d0/BB2*( &
                    +((grad_B_A3_dot_U__phi__p)+1.d0/BB2**2*BB2_A3*B_dot_U0*BB2_phi &
                    -1.d0/BB2*(B_A3_dot_U__p)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_A3__phi__p))*(grad_BRVstar_phi_k)/R**2 &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*(BR0_A3*v_p__k)/R**2) &
                    )
                  Qjac_UR_visco_parpar_A3_k = ( &
                    -2.d0/BB2**3*BB2_A3*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BRVstar__k) &
                    +1.d0/BB2**2*(( &
                    1.d0/BB2**2*BB2_A3*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__A3__p)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_A3_dot_U__p)+(B_A3_dot_grad_BdotU0__p)+(B_dot_grad__BdotU0_A3__p) &
                    )*(B_dot_grad_BRVstar__k) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BRVstar__A3__k)) &
                    )
                  Qjac_UR_visco_par_A3_kn = ( &
                    +1.d0/BB2*( &
                    +((grad_B_A3_dot_U__phi__n) &
                    -1.d0/BB2*(B_A3_dot_U__n)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_A3__phi__n))*(grad_BRVstar_phi_k)/R**2 &
                    ) &
                    )
                  Qjac_UR_visco_parpar_A3_kn = ( &
                    +1.d0/BB2**2*(( &
                    -1.d0/BB2*(B_dot_grad_BB2__A3__n)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_A3_dot_U__n)+(B_A3_dot_grad_BdotU0__n)+(B_dot_grad__BdotU0_A3__n) &
                    )*(B_dot_grad_BRVstar__k) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BRVstar__A3__kn)) &
                    )
                  Qjac_UR_visco_par_A3_pnn = ( &
                    +1.d0/BB2*( &
                    ( &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_A3__R__nn))*grad_BRVstar_R &
                    +( &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_A3__Z__nn))*grad_BRVstar_Z &
                    +((grad_B_A3_dot_U__phi__nn) &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_A3__phi__nn))*(grad_BRVstar_phi)/R**2 &
                    ) &
                    )
                  Qjac_UR_visco_parpar_A3_pnn = ( &
                    +1.d0/BB2**2*(( &
                    -1.d0/BB2*(B_dot_grad_BB2__A3__nn)*B_dot_U0 &
                    +(B_dot_grad__BdotU0_A3__nn) &
                    )*(B_dot_grad_BRVstar__p) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BRVstar__A3__nn)) &
                    )
                  Qjac_UR_visco_par_A3_knn = ( &
                    +1.d0/BB2*( &
                    +((grad_B_A3_dot_U__phi__nn) &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_A3__phi__nn))*(grad_BRVstar_phi_k)/R**2 &
                    ) &
                    )
                  Qjac_UR_visco_parpar_A3_knn = ( &
                    +1.d0/BB2**2*(( &
                    -1.d0/BB2*(B_dot_grad_BB2__A3__nn)*B_dot_U0 &
                    +(B_dot_grad__BdotU0_A3__nn) &
                    )*(B_dot_grad_BRVstar__k) &
                    ) &
                    )

                  Qjac_UR_visco_par_A3_for_Uphi_p = &
                    - v * 1.d0/BB2**2*BB2_A3*( &
                     (grad_BdotU0_R-1.d0/BB2*B_dot_U0*BB2_R) * BR0_R &
                    +(grad_BdotU0_Z-1.d0/BB2*B_dot_U0*BB2_Z) * BR0_Z &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi) * BR0_p/R**2) &
                    + v * 1.d0/BB2*( &
                     ((grad_B_A3_dot_U__R__p)+1.d0/BB2**2*BB2_A3*B_dot_U0*BB2_R &
                    -1.d0/BB2*(B_A3_dot_U__p)*BB2_R-1.d0/BB2*B_dot_U0*(grad_BB2_A3__R__p))*BR0_R &
                    +(grad_BdotU0_R-1.d0/BB2*B_dot_U0*BB2_R)*grad_BR_A3__R &
                    +((grad_B_A3_dot_U__Z__p)+1.d0/BB2**2*BB2_A3*B_dot_U0*BB2_Z &
                    -1.d0/BB2*(B_A3_dot_U__p)*BB2_Z-1.d0/BB2*B_dot_U0*(grad_BB2_A3__Z__p))*BR0_Z &
                    +(grad_BdotU0_Z-1.d0/BB2*B_dot_U0*BB2_Z)*grad_BR_A3__Z &
                    +((grad_B_A3_dot_U__phi__p)+1.d0/BB2**2*BB2_A3*B_dot_U0*BB2_phi &
                    -1.d0/BB2*(B_A3_dot_U__p)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_A3__phi__p))*BR0_p/R**2 &
                    ) 

                  Qjac_UR_visco_par_A3_for_Uphi_n = &
                    + v * 1.d0/BB2*( &
                     ((grad_B_A3_dot_U__R__n) &
                    -1.d0/BB2*(B_A3_dot_U__n)*BB2_R-1.d0/BB2*B_dot_U0*(grad_BB2_A3__R__n))*BR0_R &
                    +((grad_B_A3_dot_U__Z__n) &
                    -1.d0/BB2*(B_A3_dot_U__n)*BB2_Z-1.d0/BB2*B_dot_U0*(grad_BB2_A3__Z__n))*BR0_Z &
                    +((grad_B_A3_dot_U__phi__n) &
                    -1.d0/BB2*(B_A3_dot_U__n)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_A3__phi__n))*BR0_p/R**2 &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*grad_BR_A3__phi__n/R**2) 

                  Qjac_UR_visco_par_A3_for_Uphi_pnn = &
                    + v * 1.d0/BB2*( &
                     ( &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_A3__R__nn))*BR0_R &
                    +( &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_A3__Z__nn))*BR0_Z &
                    +((grad_B_A3_dot_U__phi__nn) &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_A3__phi__nn))*BR0_p/R**2 &
                    ) 

                  Qjac_UR_visco_parpar_A3_for_Uphi_p = &
                    - v * 2.d0/BB2**3*BB2_A3*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*B_dot_grad_BR &
                    + v * 1.d0/BB2**2*(( &
                     1.d0/BB2**2*BB2_A3*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__A3__p)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_A3_dot_U__p)+(B_A3_dot_grad_BdotU0__p)+(B_dot_grad__BdotU0_A3__p) &
                    )*B_dot_grad_BR &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BR__A3__p)) 

                  Qjac_UR_visco_parpar_A3_for_Uphi_n = &
                    + v * 1.d0/BB2**2*(( &
                    -1.d0/BB2*(B_dot_grad_BB2__A3__n)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_A3_dot_U__n)+(B_A3_dot_grad_BdotU0__n)+(B_dot_grad__BdotU0_A3__n) &
                    )*B_dot_grad_BR &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BR__A3__n)) 

                  Qjac_UR_visco_parpar_A3_for_Uphi_pnn = &
                    + v * 1.d0/BB2**2*(( &
                    -1.d0/BB2*(B_dot_grad_BB2__A3__nn)*B_dot_U0 &
                    +(B_dot_grad__BdotU0_A3__nn) &
                    )*B_dot_grad_BR &
                    ) 
                  
                  ! --- var_UR,var_UR
                  Qjac_UR_visco_par_UR_p = ( &
                    1.d0/BB2*( &
                    (BR0_R*UR+BR0*UR_R-1.d0/BB2*BR0*UR*BB2_R)*(v*BR0_R+BR0*v_R) &
                    +(BR0_Z*UR+BR0*UR_Z-1.d0/BB2*BR0*UR*BB2_Z)*(v*BR0_Z+BR0*v_Z) &
                    +(BR0_p*UR-1.d0/BB2*BR0*UR*BB2_phi)*(v*BR0_p)/R**2) &
                    )
                  Qjac_UR_visco_parpar_UR_p = ( &
                    1.d0/BB2**2*(-BR0/BB2*B_dot_grad_BB2*UR+B_dot_grad_BR*UR+BR0*(B_dot_grad__UR__p)) &
                    *(B_dot_grad_BRVstar__p) &
                    )
                  Qjac_UR_visco_par_UR_n = ( &
                    1.d0/BB2*( &
                    +(BR0*UR_p__n)*(v*BR0_p)/R**2) &
                    )
                  Qjac_UR_visco_parpar_UR_n = ( &
                    1.d0/BB2**2*(BR0*(B_dot_grad__UR__n)) &
                    *(B_dot_grad_BRVstar__p) &
                    )
                  Qjac_UR_visco_par_UR_k = ( &
                    1.d0/BB2*( &
                    +(BR0_p*UR-1.d0/BB2*BR0*UR*BB2_phi)*(BR0*v_p__k)/R**2) &
                    )
                  Qjac_UR_visco_parpar_UR_k = ( &
                    1.d0/BB2**2*(-BR0/BB2*B_dot_grad_BB2*UR+B_dot_grad_BR*UR+BR0*(B_dot_grad__UR__p)) &
                    *(B_dot_grad_BRVstar__k) &
                    )
                  Qjac_UR_visco_par_UR_kn = ( &
                    1.d0/BB2*( &
                    +(BR0*UR_p__n)*(BR0*v_p__k)/R**2) &
                    )
                  Qjac_UR_visco_parpar_UR_kn = ( &
                    1.d0/BB2**2*(BR0*(B_dot_grad__UR__n)) &
                    *(B_dot_grad_BRVstar__k) &
                    )
                  Qjac_UR_visco_par_UR_pnn = 0.d0
                  Qjac_UR_visco_parpar_UR_pnn = 0.d0
                  Qjac_UR_visco_par_UR_knn = 0.d0
                  Qjac_UR_visco_parpar_UR_knn = 0.d0

                  Qjac_UR_visco_par_UR_for_Uphi_p = &
                    v * 1.d0/BB2*( &
                     (BR0_R*UR+BR0*UR_R-1.d0/BB2*BR0*UR*BB2_R)*BR0_R &
                    +(BR0_Z*UR+BR0*UR_Z-1.d0/BB2*BR0*UR*BB2_Z)*BR0_Z &
                    +(BR0_p*UR-1.d0/BB2*BR0*UR*BB2_phi)*BR0_p/R**2)

                  Qjac_UR_visco_par_UR_for_Uphi_n = &
                    v * 1.d0/BB2*( &
                    +(BR0*UR_p__n)*BR0_p/R**2)

                  Qjac_UR_visco_par_UR_for_Uphi_pnn = 0.d0

                  Qjac_UR_visco_parpar_UR_for_Uphi_p = &
                    v * 1.d0/BB2**2*(-BR0/BB2*B_dot_grad_BB2*UR+B_dot_grad_BR*UR+BR0*(B_dot_grad__UR__p)) &
                    *B_dot_grad_BR 

                  Qjac_UR_visco_parpar_UR_for_Uphi_n = &
                    v * 1.d0/BB2**2*(BR0*(B_dot_grad__UR__n)) &
                    *B_dot_grad_BR 

                  Qjac_UR_visco_parpar_UR_for_Uphi_pnn = 0.d0

                  ! --- var_UR,var_UZ
                  Qjac_UR_visco_par_UZ_p = ( &
                    1.d0/BB2*( &
                    (BZ0_R*UZ+BZ0*UZ_R-1.d0/BB2*BZ0*UZ*BB2_R)*(v*BR0_R+BR0*v_R) &
                    +(BZ0_Z*UZ+BZ0*UZ_Z-1.d0/BB2*BZ0*UZ*BB2_Z)*(v*BR0_Z+BR0*v_Z) &
                    +(BZ0_p*UZ-1.d0/BB2*BZ0*UZ*BB2_phi)*(v*BR0_p)/R**2) &
                    )
                  Qjac_UR_visco_parpar_UZ_p = ( &
                    1.d0/BB2**2*(-BZ0/BB2*B_dot_grad_BB2*UZ+B_dot_grad_BZ*UZ+BZ0*(B_dot_grad__UZ__p)) &
                    *(B_dot_grad_BRVstar__p) &
                    )
                  Qjac_UR_visco_par_UZ_n = ( &
                    1.d0/BB2*( &
                    +(BZ0*UZ_p__n)*(v*BR0_p)/R**2) &
                    )
                  Qjac_UR_visco_parpar_UZ_n = ( &
                    1.d0/BB2**2*(BZ0*(B_dot_grad__UZ__n)) &
                    *(B_dot_grad_BRVstar__p) &
                    )
                  Qjac_UR_visco_par_UZ_k = ( &
                    1.d0/BB2*( &
                    +(BZ0_p*UZ-1.d0/BB2*BZ0*UZ*BB2_phi)*(BR0*v_p__k)/R**2) &
                    )
                  Qjac_UR_visco_parpar_UZ_k = ( &
                    1.d0/BB2**2*(-BZ0/BB2*B_dot_grad_BB2*UZ+B_dot_grad_BZ*UZ+BZ0*(B_dot_grad__UZ__p)) &
                    *(B_dot_grad_BRVstar__k) &
                    )
                  Qjac_UR_visco_par_UZ_kn = ( &
                    1.d0/BB2*( &
                    +(BZ0*UZ_p__n)*(BR0*v_p__k)/R**2) &
                    )
                  Qjac_UR_visco_parpar_UZ_kn = ( &
                    1.d0/BB2**2*(BZ0*(B_dot_grad__UZ__n)) &
                    *(B_dot_grad_BRVstar__k) &
                    )
                  Qjac_UR_visco_par_UZ_pnn = 0.d0
                  Qjac_UR_visco_parpar_UZ_pnn = 0.d0
                  Qjac_UR_visco_par_UZ_knn = 0.d0
                  Qjac_UR_visco_parpar_UZ_knn = 0.d0

                  Qjac_UR_visco_par_UZ_for_Uphi_p = &
                    v * 1.d0/BB2*( &
                    (BZ0_R*UZ+BZ0*UZ_R-1.d0/BB2*BZ0*UZ*BB2_R)*BR0_R &
                    +(BZ0_Z*UZ+BZ0*UZ_Z-1.d0/BB2*BZ0*UZ*BB2_Z)*BR0_Z &
                    +(BZ0_p*UZ-1.d0/BB2*BZ0*UZ*BB2_phi)*BR0_p/R**2) 

                  Qjac_UR_visco_par_UZ_for_Uphi_n = &
                    v * 1.d0/BB2*( &
                    +(BZ0*UZ_p__n)*BR0_p/R**2) 

                  Qjac_UR_visco_par_UZ_for_Uphi_pnn = 0.d0

                  Qjac_UR_visco_parpar_UZ_for_Uphi_p = &
                    v * 1.d0/BB2**2*(-BZ0/BB2*B_dot_grad_BB2*UZ+B_dot_grad_BZ*UZ+BZ0*(B_dot_grad__UZ__p)) &
                    *B_dot_grad_BR 

                  Qjac_UR_visco_parpar_UZ_for_Uphi_n = &
                    v * 1.d0/BB2**2*(BZ0*(B_dot_grad__UZ__n)) &
                    *B_dot_grad_BR 

                  Qjac_UR_visco_parpar_UZ_for_Uphi_pnn = 0.d0

                  ! --- var_UR,var_Up
                  Qjac_UR_visco_par_Up_p = ( &
                    1.d0/BB2*( &
                    (Bp0_R*Up+Bp0*Up_R-1.d0/BB2*Bp0*Up*BB2_R)*(v*BR0_R+BR0*v_R) &
                    +(Bp0_Z*Up+Bp0*Up_Z-1.d0/BB2*Bp0*Up*BB2_Z)*(v*BR0_Z+BR0*v_Z) &
                    +(Bp0_p*Up-1.d0/BB2*Bp0*Up*BB2_phi)*(v*BR0_p)/R**2) &
                    )
                  Qjac_UR_visco_parpar_Up_p = ( &
                    1.d0/BB2**2*(-Bp0/BB2*B_dot_grad_BB2*Up+B_dot_grad_Bp*Up+Bp0*(B_dot_grad__Up__p)) &
                    *(B_dot_grad_BRVstar__p) &
                    )
                  Qjac_UR_visco_par_Up_n = ( &
                    1.d0/BB2*( &
                    +(Bp0*Up_p__n)*(v*BR0_p)/R**2) &
                    )
                  Qjac_UR_visco_parpar_Up_n = ( &
                    1.d0/BB2**2*(Bp0*(B_dot_grad__Up__n)) &
                    *(B_dot_grad_BRVstar__p) &
                    )
                  Qjac_UR_visco_par_Up_k = ( &
                    1.d0/BB2*( &
                    +(-1.d0/BB2*Bp0*Up*BB2_phi)*(BR0*v_p__k)/R**2) &
                    )
                  Qjac_UR_visco_parpar_Up_k = ( &
                    1.d0/BB2**2*(-Bp0/BB2*B_dot_grad_BB2*Up+B_dot_grad_Bp*Up+Bp0*(B_dot_grad__Up__p)) &
                    *(B_dot_grad_BRVstar__k) &
                    )
                  Qjac_UR_visco_par_Up_kn = ( &
                    1.d0/BB2*( &
                    +(Bp0*Up_p__n)*(BR0*v_p__k)/R**2) &
                    )
                  Qjac_UR_visco_parpar_Up_kn = ( &
                    1.d0/BB2**2*(Bp0*(B_dot_grad__Up__n)) &
                    *(B_dot_grad_BRVstar__k) &
                    )
                  Qjac_UR_visco_par_Up_pnn = 0.d0
                  Qjac_UR_visco_parpar_Up_pnn = 0.d0
                  Qjac_UR_visco_par_Up_knn = 0.d0
                  Qjac_UR_visco_parpar_Up_knn = 0.d0

                  Qjac_UR_visco_par_Up_for_Uphi_p = &
                    v * 1.d0/BB2*( &
                     (Bp0_R*Up+Bp0*Up_R-1.d0/BB2*Bp0*Up*BB2_R)*BR0_R &
                    +(Bp0_Z*Up+Bp0*Up_Z-1.d0/BB2*Bp0*Up*BB2_Z)*BR0_Z &
                    +(Bp0_p*Up-1.d0/BB2*Bp0*Up*BB2_phi)*BR0_p/R**2)

                  Qjac_UR_visco_par_Up_for_Uphi_n = &
                    v * 1.d0/BB2*( &
                    +(Bp0*Up_p__n)*BR0_p/R**2)

                  Qjac_UR_visco_par_Up_for_Uphi_pnn = 0.d0
                   
                  Qjac_UR_visco_parpar_Up_for_Uphi_p = &
                    v * 1.d0/BB2**2*(-Bp0/BB2*B_dot_grad_BB2*Up+B_dot_grad_Bp*Up+Bp0*(B_dot_grad__Up__p)) &
                    *B_dot_grad_BR

                  Qjac_UR_visco_parpar_Up_for_Uphi_n = &
                    v * 1.d0/BB2**2*(Bp0*(B_dot_grad__Up__n)) &
                    *B_dot_grad_BR

                  Qjac_UR_visco_parpar_Up_for_Uphi_pnn = 0.d0

                  ! --- var_UZ,var_AR
                  Qjac_UZ_visco_par_AR_p = ( &
                    -1.d0/BB2**2*(BB2_AR__p)*( &
                      (grad_BdotU0_R - 1.d0/BB2*B_dot_U0*BB2_R)*grad_BZVstar_R  &
                    + (grad_BdotU0_Z - 1.d0/BB2*B_dot_U0*BB2_Z)*grad_BZVstar_Z  &
                    + (grad_BdotU0_phi - 1.d0/BB2*B_dot_U0*BB2_phi)*(grad_BZVstar_phi)/R**2)  &
                    + 1.d0/BB2*( &
                      ((grad_B_AR_dot_U__R__p)+1.d0/BB2**2*(BB2_AR__p)*B_dot_U0*BB2_R &
                    - 1.d0/BB2*(B_AR_dot_U__p)*BB2_R-1.d0/BB2*B_dot_U0*(grad_BB2_AR__R__p))*grad_BZVstar_R  &
                    + ((grad_B_AR_dot_U__Z__p)+1.d0/BB2**2*(BB2_AR__p)*B_dot_U0*BB2_Z &
                    - 1.d0/BB2*(B_AR_dot_U__p)*BB2_Z-1.d0/BB2*B_dot_U0*(grad_BB2_AR__Z__p))*grad_BZVstar_Z  &
                    + ((grad_B_AR_dot_U__phi__p)+1.d0/BB2**2*(BB2_AR__p)*B_dot_U0*BB2_phi  &
                    - 1.d0/BB2*(B_AR_dot_U__p)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_AR__phi__p))*(grad_BZVstar_phi)/R**2  &
                    ) &
                    )
                  Qjac_UZ_visco_parpar_AR_p = ( &
                    -2.d0/BB2**3*(BB2_AR__p)*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BZVstar__p) &
                    +1.d0/BB2**2*(( &
                    1.d0/BB2**2*(BB2_AR__p)*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__AR__p)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_AR_dot_U__p)+(B_AR_dot_grad_BdotU0__p)+(B_dot_grad__BdotU0_AR__p) &
                    )*(B_dot_grad_BZVstar__p) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BZVstar__AR__p)) &
                    )
                  Qjac_UZ_visco_par_AR_n = ( &
                    -1.d0/BB2**2*(BB2_AR__n)*( &
                      (grad_BdotU0_R - 1.d0/BB2*B_dot_U0*BB2_R)*grad_BZVstar_R  &
                    + (grad_BdotU0_Z - 1.d0/BB2*B_dot_U0*BB2_Z)*grad_BZVstar_Z  &
                    + (grad_BdotU0_phi - 1.d0/BB2*B_dot_U0*BB2_phi)*(grad_BZVstar_phi)/R**2)  &
                    + 1.d0/BB2*( &
                      ((grad_B_AR_dot_U__R__n)+1.d0/BB2**2*(BB2_AR__n)*B_dot_U0*BB2_R &
                    - 1.d0/BB2*(B_AR_dot_U__n)*BB2_R-1.d0/BB2*B_dot_U0*(grad_BB2_AR__R__n))*grad_BZVstar_R  &
                    + (grad_BdotU0_R-1.d0/BB2*B_dot_U0*BB2_R)*(v*grad_BZ_AR__R__n+BZ0_AR__n*v_R) &
                    + ((grad_B_AR_dot_U__Z__n)+1.d0/BB2**2*(BB2_AR__n)*B_dot_U0*BB2_Z &
                    - 1.d0/BB2*(B_AR_dot_U__n)*BB2_Z-1.d0/BB2*B_dot_U0*(grad_BB2_AR__Z__n))*grad_BZVstar_Z  &
                    + (grad_BdotU0_Z-1.d0/BB2*B_dot_U0*BB2_Z)*(v*grad_BZ_AR__Z__n+BZ0_AR__n*v_Z) &
                    + ((grad_B_AR_dot_U__phi__n)+1.d0/BB2**2*(BB2_AR__n)*B_dot_U0*BB2_phi  &
                    - 1.d0/BB2*(B_AR_dot_U__n)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_AR__phi__n))*(grad_BZVstar_phi)/R**2  &
                    ) &
                    )
                  Qjac_UZ_visco_parpar_AR_n = ( &
                    -2.d0/BB2**3*(BB2_AR__n)*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BZVstar__p) &
                    +1.d0/BB2**2*(( &
                    1.d0/BB2**2*(BB2_AR__n)*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__AR__n)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_AR_dot_U__n)+(B_AR_dot_grad_BdotU0__n)+(B_dot_grad__BdotU0_AR__n) &
                    )*(B_dot_grad_BZVstar__p) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BZVstar__AR__n)) &
                    )
                  Qjac_UZ_visco_par_AR_k = ( &
                    -1.d0/BB2**2*(BB2_AR__p)*( &
                    + (grad_BdotU0_phi - 1.d0/BB2*B_dot_U0*BB2_phi)*(grad_BZVstar_phi_k)/R**2)  &
                    + 1.d0/BB2*( &
                    + ((grad_B_AR_dot_U__phi__p)+1.d0/BB2**2*(BB2_AR__p)*B_dot_U0*BB2_phi  &
                    - 1.d0/BB2*(B_AR_dot_U__p)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_AR__phi__p))*(grad_BZVstar_phi_k)/R**2  &
                    ) &
                    )
                  Qjac_UZ_visco_parpar_AR_k = ( &
                    -2.d0/BB2**3*(BB2_AR__p)*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BZVstar__k) &
                    +1.d0/BB2**2*(( &
                    1.d0/BB2**2*(BB2_AR__p)*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__AR__p)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_AR_dot_U__p)+(B_AR_dot_grad_BdotU0__p)+(B_dot_grad__BdotU0_AR__p) &
                    )*(B_dot_grad_BZVstar__k) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BZVstar__AR__k)) &
                    )
                  Qjac_UZ_visco_par_AR_kn = ( &
                    -1.d0/BB2**2*(BB2_AR__n)*( &
                    + (grad_BdotU0_phi - 1.d0/BB2*B_dot_U0*BB2_phi)*(grad_BZVstar_phi_k)/R**2)  &
                    + 1.d0/BB2*( &
                    + ((grad_B_AR_dot_U__phi__n)+1.d0/BB2**2*(BB2_AR__n)*B_dot_U0*BB2_phi  &
                    - 1.d0/BB2*(B_AR_dot_U__n)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_AR__phi__n))*(grad_BZVstar_phi_k)/R**2  &
                    + (grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*(BZ0_AR__n*v_p__k)/R**2 ) &
                    )
                  Qjac_UZ_visco_parpar_AR_kn = ( &
                    -2.d0/BB2**3*(BB2_AR__n)*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BZVstar__k) &
                    +1.d0/BB2**2*(( &
                    1.d0/BB2**2*(BB2_AR__n)*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__AR__n)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_AR_dot_U__n)+(B_AR_dot_grad_BdotU0__n)+(B_dot_grad__BdotU0_AR__n) &
                    )*(B_dot_grad_BZVstar__k) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BZVstar__AR__kn)) &
                    )
                  Qjac_UZ_visco_par_AR_pnn = ( &
                    + 1.d0/BB2*( &
                      ( &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AR__R__nn))*grad_BZVstar_R  &
                    + ( &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AR__Z__nn))*grad_BZVstar_Z  &
                    + ((grad_B_AR_dot_U__phi__nn)  &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AR__phi__nn))*(grad_BZVstar_phi)/R**2  &
                    + (grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*(v*grad_BZ_AR__phi__nn)/R**2 ) &
                    )
                  Qjac_UZ_visco_parpar_AR_pnn = ( &
                    +1.d0/BB2**2*(( &
                    -1.d0/BB2*(B_dot_grad_BB2__AR__nn)*B_dot_U0 &
                    +(B_dot_grad__BdotU0_AR__nn) &
                    )*(B_dot_grad_BZVstar__p) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BZVstar__AR__nn)) &
                    )
                  Qjac_UZ_visco_par_AR_knn = ( &
                    + 1.d0/BB2*( &
                    + ((grad_B_AR_dot_U__phi__nn)  &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AR__phi__nn))*(grad_BZVstar_phi_k)/R**2  &
                    ) &
                    )
                  Qjac_UZ_visco_parpar_AR_knn = ( &
                    +1.d0/BB2**2*(( &
                    -1.d0/BB2*(B_dot_grad_BB2__AR__nn)*B_dot_U0 &
                    +(B_dot_grad__BdotU0_AR__nn) &
                    )*(B_dot_grad_BZVstar__k) &
                    ) &
                    )

                  Qjac_UZ_visco_par_AR_for_Uphi_p = &
                    - v * 1.d0/BB2**2*(BB2_AR__p)*( &
                      (grad_BdotU0_R - 1.d0/BB2*B_dot_U0*BB2_R)*BZ0_R  &
                    + (grad_BdotU0_Z - 1.d0/BB2*B_dot_U0*BB2_Z)*BZ0_Z  &
                    + (grad_BdotU0_phi - 1.d0/BB2*B_dot_U0*BB2_phi)*BZ0_p/R**2)  &
                    + v * 1.d0/BB2*( &
                      ((grad_B_AR_dot_U__R__p)+1.d0/BB2**2*(BB2_AR__p)*B_dot_U0*BB2_R &
                    - 1.d0/BB2*(B_AR_dot_U__p)*BB2_R-1.d0/BB2*B_dot_U0*(grad_BB2_AR__R__p))*BZ0_R  &
                    + ((grad_B_AR_dot_U__Z__p)+1.d0/BB2**2*(BB2_AR__p)*B_dot_U0*BB2_Z &
                    - 1.d0/BB2*(B_AR_dot_U__p)*BB2_Z-1.d0/BB2*B_dot_U0*(grad_BB2_AR__Z__p))*BZ0_Z  &
                    + ((grad_B_AR_dot_U__phi__p)+1.d0/BB2**2*(BB2_AR__p)*B_dot_U0*BB2_phi  &
                    - 1.d0/BB2*(B_AR_dot_U__p)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_AR__phi__p))*BZ0_p/R**2  &
                     ) 

                  Qjac_UZ_visco_par_AR_for_Uphi_n = &
                    - v * 1.d0/BB2**2*(BB2_AR__n)*( &
                      (grad_BdotU0_R - 1.d0/BB2*B_dot_U0*BB2_R)*BZ0_R  &
                    + (grad_BdotU0_Z - 1.d0/BB2*B_dot_U0*BB2_Z)*BZ0_Z  &
                    + (grad_BdotU0_phi - 1.d0/BB2*B_dot_U0*BB2_phi)*BZ0_p/R**2)  &
                    + v * 1.d0/BB2*( &
                      ((grad_B_AR_dot_U__R__n)+1.d0/BB2**2*(BB2_AR__n)*B_dot_U0*BB2_R &
                    - 1.d0/BB2*(B_AR_dot_U__n)*BB2_R-1.d0/BB2*B_dot_U0*(grad_BB2_AR__R__n))*BZ0_R  &
                    + (grad_BdotU0_R-1.d0/BB2*B_dot_U0*BB2_R)*grad_BZ_AR__R__n &
                    + ((grad_B_AR_dot_U__Z__n)+1.d0/BB2**2*(BB2_AR__n)*B_dot_U0*BB2_Z &
                    - 1.d0/BB2*(B_AR_dot_U__n)*BB2_Z-1.d0/BB2*B_dot_U0*(grad_BB2_AR__Z__n))*BZ0_Z  &
                    + (grad_BdotU0_Z-1.d0/BB2*B_dot_U0*BB2_Z)*grad_BZ_AR__Z__n &
                    + ((grad_B_AR_dot_U__phi__n)+1.d0/BB2**2*(BB2_AR__n)*B_dot_U0*BB2_phi  &
                    - 1.d0/BB2*(B_AR_dot_U__n)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_AR__phi__n))*BZ0_p/R**2  &
                    ) 

                  Qjac_UZ_visco_par_AR_for_Uphi_pnn = &
                    + v * 1.d0/BB2*( &
                      ( &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AR__R__nn))*BZ0_R  &
                    + ( &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AR__Z__nn))*BZ0_Z  &
                    + ((grad_B_AR_dot_U__phi__nn)  &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AR__phi__nn))*BZ0_p/R**2  &
                    + (grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*grad_BZ_AR__phi__nn/R**2 ) 

                  Qjac_UZ_visco_parpar_AR_for_Uphi_p = &
                    - v * 2.d0/BB2**3*(BB2_AR__p)*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*B_dot_grad_BZ &
                    + v * 1.d0/BB2**2*(( &
                    1.d0/BB2**2*(BB2_AR__p)*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__AR__p)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_AR_dot_U__p)+(B_AR_dot_grad_BdotU0__p)+(B_dot_grad__BdotU0_AR__p) &
                    )*B_dot_grad_BZ &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BZ__AR__p)) 

                  Qjac_UZ_visco_parpar_AR_for_Uphi_n = &
                    - v * 2.d0/BB2**3*(BB2_AR__n)*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*B_dot_grad_BZ &
                    + v * 1.d0/BB2**2*(( &
                    1.d0/BB2**2*(BB2_AR__n)*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__AR__n)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_AR_dot_U__n)+(B_AR_dot_grad_BdotU0__n)+(B_dot_grad__BdotU0_AR__n) &
                    )*B_dot_grad_BZ &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BZ__AR__n)) 

                  Qjac_UZ_visco_parpar_AR_for_Uphi_pnn = &
                    + v * 1.d0/BB2**2*(( &
                    -1.d0/BB2*(B_dot_grad_BB2__AR__nn)*B_dot_U0 &
                    +(B_dot_grad__BdotU0_AR__nn) &
                    )*B_dot_grad_BZ &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BZ__AR__nn)) 

                  ! --- var_UZ,var_AZ
                  Qjac_UZ_visco_par_AZ_p = ( &
                    -1.d0/BB2**2*(BB2_AZ__p)*( &
                    (grad_BdotU0_R-1.d0/BB2*B_dot_U0*BB2_R)*grad_BZVstar_R &
                    +(grad_BdotU0_Z-1.d0/BB2*B_dot_U0*BB2_Z)*grad_BZVstar_Z &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*(grad_BZVstar_phi)/R**2) &
                    +1.d0/BB2*( &
                    ((grad_B_AZ_dot_U__R__p) + 1.d0/BB2**2*(BB2_AZ__p)*B_dot_U0*BB2_R &
                    -1.d0/BB2*(B_AZ_dot_U__p)*BB2_R-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__R__p))*grad_BZVstar_R &
                    +(grad_BdotU0_R-1.d0/BB2*B_dot_U0*BB2_R)*(v*grad_BZ_AZ__R+BZ0_AZ*v_R) &
                    +((grad_B_AZ_dot_U__Z__p)+1.d0/BB2**2*(BB2_AZ__p)*B_dot_U0*BB2_Z &
                    -1.d0/BB2*(B_AZ_dot_U__p)*BB2_Z-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__Z__p))*grad_BZVstar_Z &
                    +(grad_BdotU0_Z-1.d0/BB2*B_dot_U0*BB2_Z)*(v*grad_BZ_AZ__Z+BZ0_AZ*v_Z) &
                    +((grad_B_AZ_dot_U__phi__p)+1.d0/BB2**2*(BB2_AZ__p)*B_dot_U0*BB2_phi &
                    -1.d0/BB2*(B_AZ_dot_U__p)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__phi__p))*(grad_BZVstar_phi)/R**2 &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*(v*grad_BZ_AZ__phi)/R**2) &
                    )
                  Qjac_UZ_visco_parpar_AZ_p = ( &
                    -2.d0/BB2**3*(BB2_AZ__p)*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BZVstar__p) &
                    +1.d0/BB2**2*(( &
                    1.d0/BB2**2*(BB2_AZ__p)*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__AZ__p)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_AZ_dot_U__p)+(B_AZ_dot_grad_BdotU0__p)+(B_dot_grad__BdotU0_AZ__p) &
                    )*(B_dot_grad_BZVstar__p) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BZVstar__AZ__p)) &
                    )
                  Qjac_UZ_visco_par_AZ_n = ( &
                    -1.d0/BB2**2*(BB2_AZ__n)*( &
                    (grad_BdotU0_R-1.d0/BB2*B_dot_U0*BB2_R)*grad_BZVstar_R &
                    +(grad_BdotU0_Z-1.d0/BB2*B_dot_U0*BB2_Z)*grad_BZVstar_Z &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*(grad_BZVstar_phi)/R**2) &
                    +1.d0/BB2*( &
                    ((grad_B_AZ_dot_U__R__n) + 1.d0/BB2**2*(BB2_AZ__n)*B_dot_U0*BB2_R &
                    -1.d0/BB2*(B_AZ_dot_U__n)*BB2_R-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__R__n))*grad_BZVstar_R &
                    +((grad_B_AZ_dot_U__Z__n)+1.d0/BB2**2*(BB2_AZ__n)*B_dot_U0*BB2_Z &
                    -1.d0/BB2*(B_AZ_dot_U__n)*BB2_Z-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__Z__n))*grad_BZVstar_Z &
                    +((grad_B_AZ_dot_U__phi__n)+1.d0/BB2**2*(BB2_AZ__n)*B_dot_U0*BB2_phi &
                    -1.d0/BB2*(B_AZ_dot_U__n)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__phi__n))*(grad_BZVstar_phi)/R**2 &
                    ) &
                    )
                  Qjac_UZ_visco_parpar_AZ_n = ( &
                    -2.d0/BB2**3*(BB2_AZ__n)*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BZVstar__p) &
                    +1.d0/BB2**2*(( &
                    1.d0/BB2**2*(BB2_AZ__n)*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__AZ__n)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_AZ_dot_U__n)+(B_AZ_dot_grad_BdotU0__n)+(B_dot_grad__BdotU0_AZ__n) &
                    )*(B_dot_grad_BZVstar__p) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BZVstar__AZ__n)) &
                    )
                  Qjac_UZ_visco_par_AZ_k = ( &
                    -1.d0/BB2**2*(BB2_AZ__p)*( &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*(grad_BZVstar_phi_k)/R**2) &
                    +1.d0/BB2*( &
                    +((grad_B_AZ_dot_U__phi__p)+1.d0/BB2**2*(BB2_AZ__p)*B_dot_U0*BB2_phi &
                    -1.d0/BB2*(B_AZ_dot_U__p)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__phi__p))*(grad_BZVstar_phi_k)/R**2 &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*(BZ0_AZ*v_p__k)/R**2) &
                    )
                  Qjac_UZ_visco_parpar_AZ_k = ( &
                    -2.d0/BB2**3*(BB2_AZ__p)*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BZVstar__k) &
                    +1.d0/BB2**2*(( &
                    1.d0/BB2**2*(BB2_AZ__p)*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__AZ__p)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_AZ_dot_U__p)+(B_AZ_dot_grad_BdotU0__p)+(B_dot_grad__BdotU0_AZ__p) &
                    )*(B_dot_grad_BZVstar__k) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BZVstar__AZ__k)) &
                    )
                  Qjac_UZ_visco_par_AZ_kn = ( &
                    -1.d0/BB2**2*(BB2_AZ__n)*( &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*(grad_BZVstar_phi_k)/R**2) &
                    +1.d0/BB2*( &
                    +((grad_B_AZ_dot_U__phi__n)+1.d0/BB2**2*(BB2_AZ__n)*B_dot_U0*BB2_phi &
                    -1.d0/BB2*(B_AZ_dot_U__n)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__phi__n))*(grad_BZVstar_phi_k)/R**2 &
                    ) &
                    )
                  Qjac_UZ_visco_parpar_AZ_kn = ( &
                    -2.d0/BB2**3*(BB2_AZ__n)*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BZVstar__k) &
                    +1.d0/BB2**2*(( &
                    1.d0/BB2**2*(BB2_AZ__n)*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__AZ__n)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_AZ_dot_U__n)+(B_AZ_dot_grad_BdotU0__n)+(B_dot_grad__BdotU0_AZ__n) &
                    )*(B_dot_grad_BZVstar__k) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BZVstar__AZ__kn)) &
                    )
                  Qjac_UZ_visco_par_AZ_pnn = ( &
                    +1.d0/BB2*( &
                    ( &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AZ__R__nn))*grad_BZVstar_R &
                    +( &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AZ__Z__nn))*grad_BZVstar_Z &
                    +((grad_B_AZ_dot_U__phi__nn) &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AZ__phi__nn))*(grad_BZVstar_phi)/R**2 &
                    ) &
                    )
                  Qjac_UZ_visco_parpar_AZ_pnn = ( &
                    +1.d0/BB2**2*(( &
                    -1.d0/BB2*(B_dot_grad_BB2__AZ__nn)*B_dot_U0 &
                    +(B_dot_grad__BdotU0_AZ__nn) &
                    )*(B_dot_grad_BZVstar__p) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BZVstar__AZ__nn)) &
                    )
                  Qjac_UZ_visco_par_AZ_knn = ( &
                    +1.d0/BB2*( &
                    +((grad_B_AZ_dot_U__phi__nn) &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AZ__phi__nn))*(grad_BZVstar_phi_k)/R**2 &
                    ) &
                    )
                  Qjac_UZ_visco_parpar_AZ_knn = ( &
                    +1.d0/BB2**2*(( &
                    -1.d0/BB2*(B_dot_grad_BB2__AZ__nn)*B_dot_U0 &
                    +(B_dot_grad__BdotU0_AZ__nn) &
                    )*(B_dot_grad_BZVstar__k) &
                    ) &
                    )

                  Qjac_UZ_visco_par_AZ_for_Uphi_p = &
                    - v * 1.d0/BB2**2*(BB2_AZ__p)*( &
                    (grad_BdotU0_R-1.d0/BB2*B_dot_U0*BB2_R)*BZ0_R &
                    +(grad_BdotU0_Z-1.d0/BB2*B_dot_U0*BB2_Z)*BZ0_Z &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*BZ0_p/R**2) &
                    + v * 1.d0/BB2*( &
                    ((grad_B_AZ_dot_U__R__p) + 1.d0/BB2**2*(BB2_AZ__p)*B_dot_U0*BB2_R &
                    -1.d0/BB2*(B_AZ_dot_U__p)*BB2_R-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__R__p))*BZ0_R &
                    +(grad_BdotU0_R-1.d0/BB2*B_dot_U0*BB2_R)*grad_BZ_AZ__R &
                    +((grad_B_AZ_dot_U__Z__p)+1.d0/BB2**2*(BB2_AZ__p)*B_dot_U0*BB2_Z &
                    -1.d0/BB2*(B_AZ_dot_U__p)*BB2_Z-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__Z__p))*BZ0_Z &
                    +(grad_BdotU0_Z-1.d0/BB2*B_dot_U0*BB2_Z)*grad_BZ_AZ__Z &
                    +((grad_B_AZ_dot_U__phi__p)+1.d0/BB2**2*(BB2_AZ__p)*B_dot_U0*BB2_phi &
                    -1.d0/BB2*(B_AZ_dot_U__p)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__phi__p))*BZ0_p/R**2 &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*grad_BZ_AZ__phi/R**2) 

                  Qjac_UZ_visco_par_AZ_for_Uphi_n = &
                    - v * 1.d0/BB2**2*(BB2_AZ__n)*( &
                    (grad_BdotU0_R-1.d0/BB2*B_dot_U0*BB2_R)*BZ0_R &
                    +(grad_BdotU0_Z-1.d0/BB2*B_dot_U0*BB2_Z)*BZ0_Z &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*BZ0_p/R**2) &
                    + v * 1.d0/BB2*( &
                    ((grad_B_AZ_dot_U__R__n) + 1.d0/BB2**2*(BB2_AZ__n)*B_dot_U0*BB2_R &
                    -1.d0/BB2*(B_AZ_dot_U__n)*BB2_R-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__R__n))*BZ0_R &
                    +((grad_B_AZ_dot_U__Z__n)+1.d0/BB2**2*(BB2_AZ__n)*B_dot_U0*BB2_Z &
                    -1.d0/BB2*(B_AZ_dot_U__n)*BB2_Z-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__Z__n))*BZ0_Z &
                    +((grad_B_AZ_dot_U__phi__n)+1.d0/BB2**2*(BB2_AZ__n)*B_dot_U0*BB2_phi &
                    -1.d0/BB2*(B_AZ_dot_U__n)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__phi__n))*BZ0_p/R**2 &
                    ) 

                  Qjac_UZ_visco_par_AZ_for_Uphi_pnn = &
                    + v * 1.d0/BB2*( &
                    ( &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AZ__R__nn))*BZ0_R &
                    +( &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AZ__Z__nn))*BZ0_Z &
                    +((grad_B_AZ_dot_U__phi__nn) &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AZ__phi__nn))*BZ0_p/R**2 &
                    ) 

                  Qjac_UZ_visco_parpar_AZ_for_Uphi_p = &
                    - v * 2.d0/BB2**3*(BB2_AZ__p)*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*B_dot_grad_BZ &
                    + v * 1.d0/BB2**2*(( &
                    1.d0/BB2**2*(BB2_AZ__p)*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__AZ__p)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_AZ_dot_U__p)+(B_AZ_dot_grad_BdotU0__p)+(B_dot_grad__BdotU0_AZ__p) &
                    )*B_dot_grad_BZ &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BZ__AZ__p)) 

                  Qjac_UZ_visco_parpar_AZ_for_Uphi_n = &
                    - v * 2.d0/BB2**3*(BB2_AZ__n)*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*B_dot_grad_BZ &
                    + v * 1.d0/BB2**2*(( &
                    1.d0/BB2**2*(BB2_AZ__n)*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__AZ__n)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_AZ_dot_U__n)+(B_AZ_dot_grad_BdotU0__n)+(B_dot_grad__BdotU0_AZ__n) &
                    )*B_dot_grad_BZ &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BZ__AZ__n)) 

                  Qjac_UZ_visco_parpar_AZ_for_Uphi_pnn = &
                    + v * 1.d0/BB2**2*(( &
                    -1.d0/BB2*(B_dot_grad_BB2__AZ__nn)*B_dot_U0 &
                    +(B_dot_grad__BdotU0_AZ__nn) &
                    )*B_dot_grad_BZ &
                    ) 

                  ! --- var_UZ,var_A3
                  Qjac_UZ_visco_par_A3_p = ( &
                    -1.d0/BB2**2*BB2_A3*( &
                    (grad_BdotU0_R-1.d0/BB2*B_dot_U0*BB2_R) * grad_BZVstar_R &
                    +(grad_BdotU0_Z-1.d0/BB2*B_dot_U0*BB2_Z) * grad_BZVstar_Z &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi) * (grad_BZVstar_phi)/R**2) &
                    +1.d0/BB2*( &
                    ((grad_B_A3_dot_U__R__p)+1.d0/BB2**2*BB2_A3*B_dot_U0*BB2_R &
                    -1.d0/BB2*(B_A3_dot_U__p)*BB2_R-1.d0/BB2*B_dot_U0*(grad_BB2_A3__R__p))*grad_BZVstar_R &
                    +(grad_BdotU0_R-1.d0/BB2*B_dot_U0*BB2_R)*(v*grad_BZ_A3__R+BZ0_A3*v_R) &
                    +((grad_B_A3_dot_U__Z__p)+1.d0/BB2**2*BB2_A3*B_dot_U0*BB2_Z &
                    -1.d0/BB2*(B_A3_dot_U__p)*BB2_Z-1.d0/BB2*B_dot_U0*(grad_BB2_A3__Z__p))*grad_BZVstar_Z &
                    +(grad_BdotU0_Z-1.d0/BB2*B_dot_U0*BB2_Z)*(v*grad_BZ_A3__Z+BZ0_A3*v_Z) &
                    +((grad_B_A3_dot_U__phi__p)+1.d0/BB2**2*BB2_A3*B_dot_U0*BB2_phi &
                    -1.d0/BB2*(B_A3_dot_U__p)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_A3__phi__p))*(grad_BZVstar_phi)/R**2 &
                    ) &
                    )
                  Qjac_UZ_visco_parpar_A3_p = ( &
                    -2.d0/BB2**3*BB2_A3*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BZVstar__p) &
                    +1.d0/BB2**2*(( &
                    1.d0/BB2**2*BB2_A3*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__A3__p)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_A3_dot_U__p)+(B_A3_dot_grad_BdotU0__p)+(B_dot_grad__BdotU0_A3__p) &
                    )*(B_dot_grad_BZVstar__p) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BZVstar__A3__p)) &
                    )
                  Qjac_UZ_visco_par_A3_n = ( &
                    +1.d0/BB2*( &
                    ((grad_B_A3_dot_U__R__n) &
                    -1.d0/BB2*(B_A3_dot_U__n)*BB2_R-1.d0/BB2*B_dot_U0*(grad_BB2_A3__R__n))*grad_BZVstar_R &
                    +((grad_B_A3_dot_U__Z__n) &
                    -1.d0/BB2*(B_A3_dot_U__n)*BB2_Z-1.d0/BB2*B_dot_U0*(grad_BB2_A3__Z__n))*grad_BZVstar_Z &
                    +((grad_B_A3_dot_U__phi__n) &
                    -1.d0/BB2*(B_A3_dot_U__n)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_A3__phi__n))*(grad_BZVstar_phi)/R**2 &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*(v*grad_BZ_A3__phi__n)/R**2) &
                    )
                  Qjac_UZ_visco_parpar_A3_n = ( &
                    +1.d0/BB2**2*(( &
                    -1.d0/BB2*(B_dot_grad_BB2__A3__n)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_A3_dot_U__n)+(B_A3_dot_grad_BdotU0__n)+(B_dot_grad__BdotU0_A3__n) &
                    )*(B_dot_grad_BZVstar__p) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BZVstar__A3__n)) &
                    )
                  Qjac_UZ_visco_par_A3_k = ( &
                    -1.d0/BB2**2*BB2_A3*( &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi) * (grad_BZVstar_phi_k)/R**2) &
                    +1.d0/BB2*( &
                    +((grad_B_A3_dot_U__phi__p)+1.d0/BB2**2*BB2_A3*B_dot_U0*BB2_phi &
                    -1.d0/BB2*(B_A3_dot_U__p)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_A3__phi__p))*(grad_BZVstar_phi_k)/R**2 &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*(BZ0_A3*v_p__k)/R**2) &
                    )
                  Qjac_UZ_visco_parpar_A3_k = ( &
                    -2.d0/BB2**3*BB2_A3*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BZVstar__k) &
                    +1.d0/BB2**2*(( &
                    1.d0/BB2**2*BB2_A3*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__A3__p)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_A3_dot_U__p)+(B_A3_dot_grad_BdotU0__p)+(B_dot_grad__BdotU0_A3__p) &
                    )*(B_dot_grad_BZVstar__k) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BZVstar__A3__k)) &
                    )
                  Qjac_UZ_visco_par_A3_kn = ( &
                    +1.d0/BB2*( &
                    +((grad_B_A3_dot_U__phi__n) &
                    -1.d0/BB2*(B_A3_dot_U__n)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_A3__phi__n))*(grad_BZVstar_phi_k)/R**2 &
                    ) &
                    )
                  Qjac_UZ_visco_parpar_A3_kn = ( &
                    +1.d0/BB2**2*(( &
                    -1.d0/BB2*(B_dot_grad_BB2__A3__n)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_A3_dot_U__n)+(B_A3_dot_grad_BdotU0__n)+(B_dot_grad__BdotU0_A3__n) &
                    )*(B_dot_grad_BZVstar__k) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BZVstar__A3__kn)) &
                    )
                  Qjac_UZ_visco_par_A3_pnn = ( &
                    +1.d0/BB2*( &
                    ( &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_A3__R__nn))*grad_BZVstar_R &
                    +( &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_A3__Z__nn))*grad_BZVstar_Z &
                    +((grad_B_A3_dot_U__phi__nn) &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_A3__phi__nn))*(grad_BZVstar_phi)/R**2 &
                    ) &
                    )
                  Qjac_UZ_visco_parpar_A3_pnn = ( &
                    +1.d0/BB2**2*(( &
                    -1.d0/BB2*(B_dot_grad_BB2__A3__nn)*B_dot_U0 &
                    +(B_dot_grad__BdotU0_A3__nn) &
                    )*(B_dot_grad_BZVstar__p) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BZVstar__A3__nn)) &
                    )
                  Qjac_UZ_visco_par_A3_knn = ( &
                    +1.d0/BB2*( &
                    +((grad_B_A3_dot_U__phi__nn) &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_A3__phi__nn))*(grad_BZVstar_phi_k)/R**2 &
                    ) &
                    )
                  Qjac_UZ_visco_parpar_A3_knn = ( &
                    +1.d0/BB2**2*(( &
                    -1.d0/BB2*(B_dot_grad_BB2__A3__nn)*B_dot_U0 &
                    +(B_dot_grad__BdotU0_A3__nn) &
                    )*(B_dot_grad_BZVstar__k) &
                    ) &
                    )

                  Qjac_UZ_visco_par_A3_for_Uphi_p = &
                    - v * 1.d0/BB2**2*BB2_A3*( &
                     (grad_BdotU0_R-1.d0/BB2*B_dot_U0*BB2_R) * BZ0_R &
                    +(grad_BdotU0_Z-1.d0/BB2*B_dot_U0*BB2_Z) * BZ0_Z &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi) * BZ0_p/R**2) &
                    + v * 1.d0/BB2*( &
                    ((grad_B_A3_dot_U__R__p)+1.d0/BB2**2*BB2_A3*B_dot_U0*BB2_R &
                    -1.d0/BB2*(B_A3_dot_U__p)*BB2_R-1.d0/BB2*B_dot_U0*(grad_BB2_A3__R__p))*BZ0_R &
                    +(grad_BdotU0_R-1.d0/BB2*B_dot_U0*BB2_R)*grad_BZ_A3__R &
                    +((grad_B_A3_dot_U__Z__p)+1.d0/BB2**2*BB2_A3*B_dot_U0*BB2_Z &
                    -1.d0/BB2*(B_A3_dot_U__p)*BB2_Z-1.d0/BB2*B_dot_U0*(grad_BB2_A3__Z__p))*BZ0_Z &
                    +(grad_BdotU0_Z-1.d0/BB2*B_dot_U0*BB2_Z)*grad_BZ_A3__Z &
                    +((grad_B_A3_dot_U__phi__p)+1.d0/BB2**2*BB2_A3*B_dot_U0*BB2_phi &
                    -1.d0/BB2*(B_A3_dot_U__p)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_A3__phi__p))*BZ0_p/R**2 &
                    ) 

                  Qjac_UZ_visco_par_A3_for_Uphi_n = &
                    + v * 1.d0/BB2*( &
                    ((grad_B_A3_dot_U__R__n) &
                    -1.d0/BB2*(B_A3_dot_U__n)*BB2_R-1.d0/BB2*B_dot_U0*(grad_BB2_A3__R__n))*BZ0_R &
                    +((grad_B_A3_dot_U__Z__n) &
                    -1.d0/BB2*(B_A3_dot_U__n)*BB2_Z-1.d0/BB2*B_dot_U0*(grad_BB2_A3__Z__n))*BZ0_Z &
                    +((grad_B_A3_dot_U__phi__n) &
                    -1.d0/BB2*(B_A3_dot_U__n)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_A3__phi__n))*BZ0_p/R**2 &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*grad_BZ_A3__phi__n/R**2) 

                  Qjac_UZ_visco_par_A3_for_Uphi_pnn = &
                    + v * 1.d0/BB2*( &
                    ( &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_A3__R__nn))*BZ0_R &
                    +( &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_A3__Z__nn))*BZ0_Z &
                    +((grad_B_A3_dot_U__phi__nn) &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_A3__phi__nn))*BZ0_p/R**2 &
                    ) 

                  Qjac_UZ_visco_parpar_A3_for_Uphi_p = &
                    - v * 2.d0/BB2**3*BB2_A3*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*B_dot_grad_BZ &
                    + v * 1.d0/BB2**2*(( &
                    1.d0/BB2**2*BB2_A3*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__A3__p)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_A3_dot_U__p)+(B_A3_dot_grad_BdotU0__p)+(B_dot_grad__BdotU0_A3__p) &
                    )*B_dot_grad_BZ &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BZ__A3__p)) 

                  Qjac_UZ_visco_parpar_A3_for_Uphi_n = &
                    + v * 1.d0/BB2**2*(( &
                    -1.d0/BB2*(B_dot_grad_BB2__A3__n)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_A3_dot_U__n)+(B_A3_dot_grad_BdotU0__n)+(B_dot_grad__BdotU0_A3__n) &
                    )*B_dot_grad_BZ &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BZ__A3__n)) 

                  Qjac_UZ_visco_parpar_A3_for_Uphi_pnn = &
                    + v * 1.d0/BB2**2*(( &
                    -1.d0/BB2*(B_dot_grad_BB2__A3__nn)*B_dot_U0 &
                    +(B_dot_grad__BdotU0_A3__nn) &
                    )*B_dot_grad_BZ &
                    ) 

                  ! --- var_UZ,var_UR
                  Qjac_UZ_visco_par_UR_p = ( &
                    1.d0/BB2*( &
                    (BR0_R*UR+BR0*UR_R-1.d0/BB2*BR0*UR*BB2_R)*(v*BZ0_R+BZ0*v_R) &
                    +(BR0_Z*UR+BR0*UR_Z-1.d0/BB2*BR0*UR*BB2_Z)*(v*BZ0_Z+BZ0*v_Z) &
                    +(BR0_p*UR-1.d0/BB2*BR0*UR*BB2_phi)*(v*BZ0_p)/R**2) &
                    )
                  Qjac_UZ_visco_parpar_UR_p = ( &
                    1.d0/BB2**2*(-BR0/BB2*B_dot_grad_BB2*UR+B_dot_grad_BR*UR+BR0*(B_dot_grad__UR__p)) &
                    *(B_dot_grad_BZVstar__p) &
                    )
                  Qjac_UZ_visco_par_UR_n = ( &
                    1.d0/BB2*( &
                    +(BR0*UR_p__n)*(v*BZ0_p)/R**2) &
                    )
                  Qjac_UZ_visco_parpar_UR_n = ( &
                    1.d0/BB2**2*(BR0*(B_dot_grad__UR__n)) &
                    *(B_dot_grad_BZVstar__p) &
                    )
                  Qjac_UZ_visco_par_UR_k = ( &
                    1.d0/BB2*( &
                    +(BR0_p*UR-1.d0/BB2*BR0*UR*BB2_phi)*(BZ0*v_p__k)/R**2) &
                    )
                  Qjac_UZ_visco_parpar_UR_k = ( &
                    1.d0/BB2**2*(-BR0/BB2*B_dot_grad_BB2*UR+B_dot_grad_BR*UR+BR0*(B_dot_grad__UR__p)) &
                    *(B_dot_grad_BZVstar__k) &
                    )
                  Qjac_UZ_visco_par_UR_kn = ( &
                    1.d0/BB2*( &
                    +(BR0*UR_p__n)*(BZ0*v_p__k)/R**2) &
                    )
                  Qjac_UZ_visco_parpar_UR_kn = ( &
                    1.d0/BB2**2*(BR0*(B_dot_grad__UR__n)) &
                    *(B_dot_grad_BZVstar__k) &
                    )
                  Qjac_UZ_visco_par_UR_pnn = 0.d0
                  Qjac_UZ_visco_parpar_UR_pnn = 0.d0
                  Qjac_UZ_visco_par_UR_knn = 0.d0
                  Qjac_UZ_visco_parpar_UR_knn = 0.d0

                  Qjac_UZ_visco_par_UR_for_Uphi_p = &
                    v * 1.d0/BB2*( &
                     (BR0_R*UR+BR0*UR_R-1.d0/BB2*BR0*UR*BB2_R)*BZ0_R &
                    +(BR0_Z*UR+BR0*UR_Z-1.d0/BB2*BR0*UR*BB2_Z)*BZ0_Z &
                    +(BR0_p*UR-1.d0/BB2*BR0*UR*BB2_phi)*BZ0_p/R**2) 

                  Qjac_UZ_visco_par_UR_for_Uphi_n = &
                    v * 1.d0/BB2*( &
                    +(BR0*UR_p__n)*BZ0_p/R**2) 

                  Qjac_UZ_visco_par_UR_for_Uphi_pnn = 0.d0

                  Qjac_UZ_visco_parpar_UR_for_Uphi_p = &
                     v * 1.d0/BB2**2*(-BR0/BB2*B_dot_grad_BB2*UR+B_dot_grad_BR*UR+BR0*(B_dot_grad__UR__p)) &
                    *B_dot_grad_BZ

                  Qjac_UZ_visco_parpar_UR_for_Uphi_n = &
                     v * 1.d0/BB2**2*(BR0*(B_dot_grad__UR__n)) &
                    *B_dot_grad_BZ

                  Qjac_UZ_visco_parpar_UR_for_Uphi_pnn = 0.d0
                  
                  ! --- var_UZ,var_UZ
                  Qjac_UZ_visco_par_UZ_p = ( &
                    1.d0/BB2*( &
                    (BZ0_R*UZ+BZ0*UZ_R-1.d0/BB2*BZ0*UZ*BB2_R)*(v*BZ0_R+BZ0*v_R) &
                    +(BZ0_Z*UZ+BZ0*UZ_Z-1.d0/BB2*BZ0*UZ*BB2_Z)*(v*BZ0_Z+BZ0*v_Z) &
                    +(BZ0_p*UZ-1.d0/BB2*BZ0*UZ*BB2_phi)*(v*BZ0_p)/R**2) &
                    )
                  Qjac_UZ_visco_parpar_UZ_p = ( &
                    1.d0/BB2**2*(-BZ0/BB2*B_dot_grad_BB2*UZ+B_dot_grad_BZ*UZ+BZ0*(B_dot_grad__UZ__p)) &
                    *(B_dot_grad_BZVstar__p) &
                    )
                  Qjac_UZ_visco_par_UZ_n = ( &
                    1.d0/BB2*( &
                    +(BZ0*UZ_p__n)*(v*BZ0_p)/R**2) &
                    )
                  Qjac_UZ_visco_parpar_UZ_n = ( &
                    1.d0/BB2**2*(BZ0*(B_dot_grad__UZ__n)) &
                    *(B_dot_grad_BZVstar__p) &
                    )
                  Qjac_UZ_visco_par_UZ_k = ( &
                    1.d0/BB2*( &
                    +(BZ0_p*UZ-1.d0/BB2*BZ0*UZ*BB2_phi)*(BZ0*v_p__k)/R**2) &
                    )
                  Qjac_UZ_visco_parpar_UZ_k = ( &
                    1.d0/BB2**2*(-BZ0/BB2*B_dot_grad_BB2*UZ+B_dot_grad_BZ*UZ+BZ0*(B_dot_grad__UZ__p)) &
                    *(B_dot_grad_BZVstar__k) &
                    )
                  Qjac_UZ_visco_par_UZ_kn = ( &
                    1.d0/BB2*( &
                    +(BZ0*UZ_p__n)*(BZ0*v_p__k)/R**2) &
                    )
                  Qjac_UZ_visco_parpar_UZ_kn = ( &
                    1.d0/BB2**2*(BZ0*(B_dot_grad__UZ__n)) &
                    *(B_dot_grad_BZVstar__k) &
                    )
                  Qjac_UZ_visco_par_UZ_pnn = 0.d0
                  Qjac_UZ_visco_parpar_UZ_pnn = 0.d0
                  Qjac_UZ_visco_par_UZ_knn = 0.d0
                  Qjac_UZ_visco_parpar_UZ_knn = 0.d0

                  Qjac_UZ_visco_par_UZ_for_Uphi_p = &
                    v * 1.d0/BB2*( &
                    (BZ0_R*UZ+BZ0*UZ_R-1.d0/BB2*BZ0*UZ*BB2_R)*BZ0_R &
                    +(BZ0_Z*UZ+BZ0*UZ_Z-1.d0/BB2*BZ0*UZ*BB2_Z)*BZ0_Z &
                    +(BZ0_p*UZ-1.d0/BB2*BZ0*UZ*BB2_phi)*BZ0_p/R**2) 

                  Qjac_UZ_visco_par_UZ_for_Uphi_n = &
                    v * 1.d0/BB2*( &
                    +(BZ0*UZ_p__n)*BZ0_p/R**2) 

                  Qjac_UZ_visco_par_UZ_for_Uphi_pnn = 0.d0

                  Qjac_UZ_visco_parpar_UZ_for_Uphi_p = &
                    v * 1.d0/BB2**2*(-BZ0/BB2*B_dot_grad_BB2*UZ+B_dot_grad_BZ*UZ+BZ0*(B_dot_grad__UZ__p)) &
                    *B_dot_grad_BZ

                  Qjac_UZ_visco_parpar_UZ_for_Uphi_n = &
                    v * 1.d0/BB2**2*(BZ0*(B_dot_grad__UZ__n)) &
                    *B_dot_grad_BZ

                  Qjac_UZ_visco_parpar_UZ_for_Uphi_pnn = 0.d0

                  ! --- var_UZ,var_Up
                  Qjac_UZ_visco_par_Up_p = ( &
                    1.d0/BB2*( &
                    (Bp0_R*Up+Bp0*Up_R-1.d0/BB2*Bp0*Up*BB2_R)*(v*BZ0_R+BZ0*v_R) &
                    +(Bp0_Z*Up+Bp0*Up_Z-1.d0/BB2*Bp0*Up*BB2_Z)*(v*BZ0_Z+BZ0*v_Z) &
                    +(Bp0_p*Up-1.d0/BB2*Bp0*Up*BB2_phi)*(v*BZ0_p)/R**2) &
                    )
                  Qjac_UZ_visco_parpar_Up_p = ( &
                    1.d0/BB2**2*(-Bp0/BB2*B_dot_grad_BB2*Up+B_dot_grad_Bp*Up+Bp0*(B_dot_grad__Up__p)) &
                    *(B_dot_grad_BZVstar__p) &
                    )
                  Qjac_UZ_visco_par_Up_n = ( &
                    1.d0/BB2*( &
                    +(Bp0*Up_p__n)*(v*BZ0_p)/R**2) &
                    )
                  Qjac_UZ_visco_parpar_Up_n = ( &
                    1.d0/BB2**2*(Bp0*(B_dot_grad__Up__n)) &
                    *(B_dot_grad_BZVstar__p) &
                    )
                  Qjac_UZ_visco_par_Up_k = ( &
                    1.d0/BB2*( &
                    +(Bp0_p*Up-1.d0/BB2*Bp0*Up*BB2_phi)*(BZ0*v_p__k)/R**2) &
                    )
                  Qjac_UZ_visco_parpar_Up_k = ( &
                    1.d0/BB2**2*(-Bp0/BB2*B_dot_grad_BB2*Up+B_dot_grad_Bp*Up+Bp0*(B_dot_grad__Up__p)) &
                    *(B_dot_grad_BZVstar__k) &
                    )
                  Qjac_UZ_visco_par_Up_kn = ( &
                    1.d0/BB2*( &
                    +(Bp0*Up_p__n)*(BZ0*v_p__k)/R**2) &
                    )
                  Qjac_UZ_visco_parpar_Up_kn = ( &
                    1.d0/BB2**2*(Bp0*(B_dot_grad__Up__n)) &
                    *(B_dot_grad_BZVstar__k) &
                    )
                  Qjac_UZ_visco_par_Up_pnn = 0.d0
                  Qjac_UZ_visco_parpar_Up_pnn = 0.d0
                  Qjac_UZ_visco_par_Up_knn = 0.d0
                  Qjac_UZ_visco_parpar_Up_knn = 0.d0

                  Qjac_UZ_visco_par_Up_for_Uphi_p = &
                    v * 1.d0/BB2*( &
                    (Bp0_R*Up+Bp0*Up_R-1.d0/BB2*Bp0*Up*BB2_R)*BZ0_R &
                    +(Bp0_Z*Up+Bp0*Up_Z-1.d0/BB2*Bp0*Up*BB2_Z)*BZ0_Z &
                    +(Bp0_p*Up-1.d0/BB2*Bp0*Up*BB2_phi)*BZ0_p/R**2) 

                  Qjac_UZ_visco_par_Up_for_Uphi_n = &
                    v * 1.d0/BB2*( &
                    +(Bp0*Up_p__n)*BZ0_p/R**2) 

                  Qjac_UZ_visco_par_Up_for_Uphi_pnn = 0.d0

                  Qjac_UZ_visco_parpar_Up_for_Uphi_p = &
                    v * 1.d0/BB2**2*(-Bp0/BB2*B_dot_grad_BB2*Up+B_dot_grad_Bp*Up+Bp0*(B_dot_grad__Up__p)) &
                    *B_dot_grad_BZ

                  Qjac_UZ_visco_parpar_Up_for_Uphi_n = &
                    v * 1.d0/BB2**2*(Bp0*(B_dot_grad__Up__n)) &
                    *B_dot_grad_BZ

                  Qjac_UZ_visco_parpar_Up_for_Uphi_pnn = 0.d0

                  ! --- var_Up,var_AR
                  Qjac_Up_visco_par_AR_p = ( &
                    -1.d0/BB2**2*(BB2_AR__p)*( &
                      (grad_BdotU0_R - 1.d0/BB2*B_dot_U0*BB2_R)*grad_BpVstar_R  &
                    + (grad_BdotU0_Z - 1.d0/BB2*B_dot_U0*BB2_Z)*grad_BpVstar_Z  &
                    + (grad_BdotU0_phi - 1.d0/BB2*B_dot_U0*BB2_phi)*(grad_BpVstar_phi)/R**2)  &
                    + 1.d0/BB2*( &
                      ((grad_B_AR_dot_U__R__p)+1.d0/BB2**2*(BB2_AR__p)*B_dot_U0*BB2_R &
                    - 1.d0/BB2*(B_AR_dot_U__p)*BB2_R-1.d0/BB2*B_dot_U0*(grad_BB2_AR__R__p))*grad_BpVstar_R  &
                    + (grad_BdotU0_R-1.d0/BB2*B_dot_U0*BB2_R)*(v*grad_Bp_AR__R+Bp0_AR*v_R) &
                    + ((grad_B_AR_dot_U__Z__p)+1.d0/BB2**2*(BB2_AR__p)*B_dot_U0*BB2_Z &
                    - 1.d0/BB2*(B_AR_dot_U__p)*BB2_Z-1.d0/BB2*B_dot_U0*(grad_BB2_AR__Z__p))*grad_BpVstar_Z  &
                    + (grad_BdotU0_Z-1.d0/BB2*B_dot_U0*BB2_Z)*(v*grad_Bp_AR__Z+Bp0_AR*v_Z) &
                    + ((grad_B_AR_dot_U__phi__p)+1.d0/BB2**2*(BB2_AR__p)*B_dot_U0*BB2_phi  &
                    - 1.d0/BB2*(B_AR_dot_U__p)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_AR__phi__p))*(grad_BpVstar_phi)/R**2  &
                    ) &
                    )
                  Qjac_Up_visco_parpar_AR_p = ( &
                    -2.d0/BB2**3*(BB2_AR__p)*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BpVstar__p) &
                    +1.d0/BB2**2*(( &
                    1.d0/BB2**2*(BB2_AR__p)*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__AR__p)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_AR_dot_U__p)+(B_AR_dot_grad_BdotU0__p)+(B_dot_grad__BdotU0_AR__p) &
                    )*(B_dot_grad_BpVstar__p) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BpVstar__AR__p)) &
                    )
                  Qjac_Up_visco_par_AR_n = ( &
                    -1.d0/BB2**2*(BB2_AR__n)*( &
                      (grad_BdotU0_R - 1.d0/BB2*B_dot_U0*BB2_R)*grad_BpVstar_R  &
                    + (grad_BdotU0_Z - 1.d0/BB2*B_dot_U0*BB2_Z)*grad_BpVstar_Z  &
                    + (grad_BdotU0_phi - 1.d0/BB2*B_dot_U0*BB2_phi)*(grad_BpVstar_phi)/R**2)  &
                    + 1.d0/BB2*( &
                      ((grad_B_AR_dot_U__R__n)+1.d0/BB2**2*(BB2_AR__n)*B_dot_U0*BB2_R &
                    - 1.d0/BB2*(B_AR_dot_U__n)*BB2_R-1.d0/BB2*B_dot_U0*(grad_BB2_AR__R__n))*grad_BpVstar_R  &
                    + ((grad_B_AR_dot_U__Z__n)+1.d0/BB2**2*(BB2_AR__n)*B_dot_U0*BB2_Z &
                    - 1.d0/BB2*(B_AR_dot_U__n)*BB2_Z-1.d0/BB2*B_dot_U0*(grad_BB2_AR__Z__n))*grad_BpVstar_Z  &
                    + ((grad_B_AR_dot_U__phi__n)+1.d0/BB2**2*(BB2_AR__n)*B_dot_U0*BB2_phi  &
                    - 1.d0/BB2*(B_AR_dot_U__n)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_AR__phi__n))*(grad_BpVstar_phi)/R**2  &
                    + (grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*(v*grad_Bp_AR__phi__n)/R**2 ) &
                    )
                  Qjac_Up_visco_parpar_AR_n = ( &
                    -2.d0/BB2**3*(BB2_AR__n)*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BpVstar__p) &
                    +1.d0/BB2**2*(( &
                    1.d0/BB2**2*(BB2_AR__n)*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__AR__n)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_AR_dot_U__n)+(B_AR_dot_grad_BdotU0__n)+(B_dot_grad__BdotU0_AR__n) &
                    )*(B_dot_grad_BpVstar__p) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BpVstar__AR__n)) &
                    )
                  Qjac_Up_visco_par_AR_k = ( &
                    -1.d0/BB2**2*(BB2_AR__p)*( &
                    + (grad_BdotU0_phi - 1.d0/BB2*B_dot_U0*BB2_phi)*(grad_BpVstar_phi_k)/R**2)  &
                    + 1.d0/BB2*( &
                    + ((grad_B_AR_dot_U__phi__p)+1.d0/BB2**2*(BB2_AR__p)*B_dot_U0*BB2_phi  &
                    - 1.d0/BB2*(B_AR_dot_U__p)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_AR__phi__p))*(grad_BpVstar_phi_k)/R**2  &
                    + (grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*(Bp0_AR*v_p__k)/R**2 ) &
                    )
                  Qjac_Up_visco_parpar_AR_k = ( &
                    -2.d0/BB2**3*(BB2_AR__p)*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BpVstar__k) &
                    +1.d0/BB2**2*(( &
                    1.d0/BB2**2*(BB2_AR__p)*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__AR__p)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_AR_dot_U__p)+(B_AR_dot_grad_BdotU0__p)+(B_dot_grad__BdotU0_AR__p) &
                    )*(B_dot_grad_BpVstar__k) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BpVstar__AR__k)) &
                    )
                  Qjac_Up_visco_par_AR_kn = ( &
                    -1.d0/BB2**2*(BB2_AR__n)*( &
                    + (grad_BdotU0_phi - 1.d0/BB2*B_dot_U0*BB2_phi)*(grad_BpVstar_phi_k)/R**2)  &
                    + 1.d0/BB2*( &
                    + ((grad_B_AR_dot_U__phi__n)+1.d0/BB2**2*(BB2_AR__n)*B_dot_U0*BB2_phi  &
                    - 1.d0/BB2*(B_AR_dot_U__n)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_AR__phi__n))*(grad_BpVstar_phi_k)/R**2  &
                    ) &
                    )
                  Qjac_Up_visco_parpar_AR_kn = ( &
                    -2.d0/BB2**3*(BB2_AR__n)*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BpVstar__k) &
                    +1.d0/BB2**2*(( &
                    1.d0/BB2**2*(BB2_AR__n)*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__AR__n)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_AR_dot_U__n)+(B_AR_dot_grad_BdotU0__n)+(B_dot_grad__BdotU0_AR__n) &
                    )*(B_dot_grad_BpVstar__k) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BpVstar__AR__kn)) &
                    )
                  Qjac_Up_visco_par_AR_pnn = ( &
                    + 1.d0/BB2*( &
                      ( &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AR__R__nn))*grad_BpVstar_R  &
                    + ( &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AR__Z__nn))*grad_BpVstar_Z  &
                    + ((grad_B_AR_dot_U__phi__nn)  &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AR__phi__nn))*(grad_BpVstar_phi)/R**2  &
                    ) &
                    )
                  Qjac_Up_visco_parpar_AR_pnn = ( &
                    +1.d0/BB2**2*(( &
                    -1.d0/BB2*(B_dot_grad_BB2__AR__nn)*B_dot_U0 &
                    +(B_dot_grad__BdotU0_AR__nn) &
                    )*(B_dot_grad_BpVstar__p) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BpVstar__AR__nn)) &
                    )
                  Qjac_Up_visco_par_AR_knn = ( &
                    + 1.d0/BB2*( &
                    + ((grad_B_AR_dot_U__phi__nn)  &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AR__phi__nn))*(grad_BpVstar_phi_k)/R**2  &
                    ) &
                    )
                  Qjac_Up_visco_parpar_AR_knn = ( &
                    +1.d0/BB2**2*(( &
                    -1.d0/BB2*(B_dot_grad_BB2__AR__nn)*B_dot_U0 &
                    +(B_dot_grad__BdotU0_AR__nn) &
                    )*(B_dot_grad_BpVstar__k) &
                    ) &
                    )

                  Qjac_Up_visco_par_AR_for_Uphi_p = &
                    - v * 1.d0/BB2**2*(BB2_AR__p)*( &
                      (grad_BdotU0_R - 1.d0/BB2*B_dot_U0*BB2_R)*Bp0_R  &
                    + (grad_BdotU0_Z - 1.d0/BB2*B_dot_U0*BB2_Z)*Bp0_Z  &
                    + (grad_BdotU0_phi - 1.d0/BB2*B_dot_U0*BB2_phi)*Bp0_p/R**2)  &
                    + v * 1.d0/BB2*( &
                      ((grad_B_AR_dot_U__R__p)+1.d0/BB2**2*(BB2_AR__p)*B_dot_U0*BB2_R &
                    - 1.d0/BB2*(B_AR_dot_U__p)*BB2_R-1.d0/BB2*B_dot_U0*(grad_BB2_AR__R__p))*Bp0_R  &
                    + (grad_BdotU0_R-1.d0/BB2*B_dot_U0*BB2_R)*grad_Bp_AR__R &
                    + ((grad_B_AR_dot_U__Z__p)+1.d0/BB2**2*(BB2_AR__p)*B_dot_U0*BB2_Z &
                    - 1.d0/BB2*(B_AR_dot_U__p)*BB2_Z-1.d0/BB2*B_dot_U0*(grad_BB2_AR__Z__p))*Bp0_Z  &
                    + (grad_BdotU0_Z-1.d0/BB2*B_dot_U0*BB2_Z)*grad_Bp_AR__Z &
                    + ((grad_B_AR_dot_U__phi__p)+1.d0/BB2**2*(BB2_AR__p)*B_dot_U0*BB2_phi  &
                    - 1.d0/BB2*(B_AR_dot_U__p)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_AR__phi__p))*Bp0_p/R**2  &
                     ) 

                  Qjac_Up_visco_par_AR_for_Uphi_n = &
                    - v * 1.d0/BB2**2*(BB2_AR__n)*( &
                      (grad_BdotU0_R - 1.d0/BB2*B_dot_U0*BB2_R)*Bp0_R  &
                    + (grad_BdotU0_Z - 1.d0/BB2*B_dot_U0*BB2_Z)*Bp0_Z  &
                    + (grad_BdotU0_phi - 1.d0/BB2*B_dot_U0*BB2_phi)*Bp0_p/R**2)  &
                    + v * 1.d0/BB2*( &
                      ((grad_B_AR_dot_U__R__n)+1.d0/BB2**2*(BB2_AR__n)*B_dot_U0*BB2_R &
                    - 1.d0/BB2*(B_AR_dot_U__n)*BB2_R-1.d0/BB2*B_dot_U0*(grad_BB2_AR__R__n))*Bp0_R  &
                    + ((grad_B_AR_dot_U__Z__n)+1.d0/BB2**2*(BB2_AR__n)*B_dot_U0*BB2_Z &
                    - 1.d0/BB2*(B_AR_dot_U__n)*BB2_Z-1.d0/BB2*B_dot_U0*(grad_BB2_AR__Z__n))*Bp0_Z  &
                    + ((grad_B_AR_dot_U__phi__n)+1.d0/BB2**2*(BB2_AR__n)*B_dot_U0*BB2_phi  &
                    - 1.d0/BB2*(B_AR_dot_U__n)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_AR__phi__n))*Bp0_p/R**2  &
                    + (grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*grad_Bp_AR__phi__n/R**2 ) 

                  Qjac_Up_visco_par_AR_for_Uphi_pnn = &
                    + v * 1.d0/BB2*( &
                      ( &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AR__R__nn))*Bp0_R  &
                    + ( &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AR__Z__nn))*Bp0_Z  &
                    + ((grad_B_AR_dot_U__phi__nn)  &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AR__phi__nn))*Bp0_p/R**2  &
                    )

                  Qjac_Up_visco_parpar_AR_for_Uphi_p = &
                    - v * 2.d0/BB2**3*(BB2_AR__p)*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*B_dot_grad_Bp &
                    + v * 1.d0/BB2**2*(( &
                    1.d0/BB2**2*(BB2_AR__p)*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__AR__p)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_AR_dot_U__p)+(B_AR_dot_grad_BdotU0__p)+(B_dot_grad__BdotU0_AR__p) &
                    )*B_dot_grad_Bp &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_Bp__AR__p)) 

                  Qjac_Up_visco_parpar_AR_for_Uphi_n = &
                    - v * 2.d0/BB2**3*(BB2_AR__n)*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*B_dot_grad_Bp &
                    + v * 1.d0/BB2**2*(( &
                    1.d0/BB2**2*(BB2_AR__n)*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__AR__n)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_AR_dot_U__n)+(B_AR_dot_grad_BdotU0__n)+(B_dot_grad__BdotU0_AR__n) &
                    )*B_dot_grad_Bp &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_Bp__AR__n)) 

                  Qjac_Up_visco_parpar_AR_for_Uphi_pnn = &
                    + v * 1.d0/BB2**2*(( &
                    -1.d0/BB2*(B_dot_grad_BB2__AR__nn)*B_dot_U0 &
                    +(B_dot_grad__BdotU0_AR__nn) &
                    )*B_dot_grad_Bp &
                    ) 

                  ! --- var_Up,var_AZ
                  Qjac_Up_visco_par_AZ_p = ( &
                    -1.d0/BB2**2*(BB2_AZ__p)*( &
                    (grad_BdotU0_R-1.d0/BB2*B_dot_U0*BB2_R)*grad_BpVstar_R &
                    +(grad_BdotU0_Z-1.d0/BB2*B_dot_U0*BB2_Z)*grad_BpVstar_Z &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*(grad_BpVstar_phi)/R**2) &
                    +1.d0/BB2*( &
                    ((grad_B_AZ_dot_U__R__p) + 1.d0/BB2**2*(BB2_AZ__p)*B_dot_U0*BB2_R &
                    -1.d0/BB2*(B_AZ_dot_U__p)*BB2_R-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__R__p))*grad_BpVstar_R &
                    +(grad_BdotU0_R-1.d0/BB2*B_dot_U0*BB2_R)*(v*grad_Bp_AZ__R+Bp0_AZ*v_R) &
                    +((grad_B_AZ_dot_U__Z__p)+1.d0/BB2**2*(BB2_AZ__p)*B_dot_U0*BB2_Z &
                    -1.d0/BB2*(B_AZ_dot_U__p)*BB2_Z-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__Z__p))*grad_BpVstar_Z &
                    +(grad_BdotU0_Z-1.d0/BB2*B_dot_U0*BB2_Z)*(v*grad_Bp_AZ__Z+Bp0_AZ*v_Z) &
                    +((grad_B_AZ_dot_U__phi__p)+1.d0/BB2**2*(BB2_AZ__p)*B_dot_U0*BB2_phi &
                    -1.d0/BB2*(B_AZ_dot_U__p)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__phi__p))*(grad_BpVstar_phi)/R**2 &
                    ) &
                    )
                  Qjac_Up_visco_parpar_AZ_p = ( &
                    -2.d0/BB2**3*(BB2_AZ__p)*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BpVstar__p) &
                    +1.d0/BB2**2*(( &
                    1.d0/BB2**2*(BB2_AZ__p)*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__AZ__p)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_AZ_dot_U__p)+(B_AZ_dot_grad_BdotU0__p)+(B_dot_grad__BdotU0_AZ__p) &
                    )*(B_dot_grad_BpVstar__p) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BpVstar__AZ__p)) &
                    )
                  Qjac_Up_visco_par_AZ_n = ( &
                    -1.d0/BB2**2*(BB2_AZ__n)*( &
                    (grad_BdotU0_R-1.d0/BB2*B_dot_U0*BB2_R)*grad_BpVstar_R &
                    +(grad_BdotU0_Z-1.d0/BB2*B_dot_U0*BB2_Z)*grad_BpVstar_Z &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*(grad_BpVstar_phi)/R**2) &
                    +1.d0/BB2*( &
                    ((grad_B_AZ_dot_U__R__n) + 1.d0/BB2**2*(BB2_AZ__n)*B_dot_U0*BB2_R &
                    -1.d0/BB2*(B_AZ_dot_U__n)*BB2_R-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__R__n))*grad_BpVstar_R &
                    +((grad_B_AZ_dot_U__Z__n)+1.d0/BB2**2*(BB2_AZ__n)*B_dot_U0*BB2_Z &
                    -1.d0/BB2*(B_AZ_dot_U__n)*BB2_Z-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__Z__n))*grad_BpVstar_Z &
                    +((grad_B_AZ_dot_U__phi__n)+1.d0/BB2**2*(BB2_AZ__n)*B_dot_U0*BB2_phi &
                    -1.d0/BB2*(B_AZ_dot_U__n)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__phi__n))*(grad_BpVstar_phi)/R**2 &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*(v*grad_Bp_AZ__phi__n)/R**2) &
                    )
                  Qjac_Up_visco_parpar_AZ_n = ( &
                    -2.d0/BB2**3*(BB2_AZ__n)*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BpVstar__p) &
                    +1.d0/BB2**2*(( &
                    1.d0/BB2**2*(BB2_AZ__n)*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__AZ__n)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_AZ_dot_U__n)+(B_AZ_dot_grad_BdotU0__n)+(B_dot_grad__BdotU0_AZ__n) &
                    )*(B_dot_grad_BpVstar__p) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BpVstar__AZ__n)) &
                    )
                  Qjac_Up_visco_par_AZ_k = ( &
                    -1.d0/BB2**2*(BB2_AZ__p)*( &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*(grad_BpVstar_phi_k)/R**2) &
                    +1.d0/BB2*( &
                    +((grad_B_AZ_dot_U__phi__p)+1.d0/BB2**2*(BB2_AZ__p)*B_dot_U0*BB2_phi &
                    -1.d0/BB2*(B_AZ_dot_U__p)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__phi__p))*(grad_BpVstar_phi_k)/R**2 &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*(Bp0_AZ*v_p__k)/R**2) &
                    )
                  Qjac_Up_visco_parpar_AZ_k = ( &
                    -2.d0/BB2**3*(BB2_AZ__p)*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BpVstar__k) &
                    +1.d0/BB2**2*(( &
                    1.d0/BB2**2*(BB2_AZ__p)*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__AZ__p)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_AZ_dot_U__p)+(B_AZ_dot_grad_BdotU0__p)+(B_dot_grad__BdotU0_AZ__p) &
                    )*(B_dot_grad_BpVstar__k) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BpVstar__AZ__k)) &
                    )
                  Qjac_Up_visco_par_AZ_kn = ( &
                    -1.d0/BB2**2*(BB2_AZ__n)*( &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*(grad_BpVstar_phi_k)/R**2) &
                    +1.d0/BB2*( &
                    +((grad_B_AZ_dot_U__phi__n)+1.d0/BB2**2*(BB2_AZ__n)*B_dot_U0*BB2_phi &
                    -1.d0/BB2*(B_AZ_dot_U__n)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__phi__n))*(grad_BpVstar_phi_k)/R**2 &
                    ) &
                    )
                  Qjac_Up_visco_parpar_AZ_kn = ( &
                    -2.d0/BB2**3*(BB2_AZ__n)*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BpVstar__k) &
                    +1.d0/BB2**2*(( &
                    1.d0/BB2**2*(BB2_AZ__n)*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__AZ__n)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_AZ_dot_U__n)+(B_AZ_dot_grad_BdotU0__n)+(B_dot_grad__BdotU0_AZ__n) &
                    )*(B_dot_grad_BpVstar__k) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BpVstar__AZ__kn)) &
                    )
                  Qjac_Up_visco_par_AZ_pnn = ( &
                    +1.d0/BB2*( &
                    ( &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AZ__R__nn))*grad_BpVstar_R &
                    +( &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AZ__Z__nn))*grad_BpVstar_Z &
                    +((grad_B_AZ_dot_U__phi__nn) &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AZ__phi__nn))*(grad_BpVstar_phi)/R**2 &
                    ) &
                    )
                  Qjac_Up_visco_parpar_AZ_pnn = ( &
                    +1.d0/BB2**2*(( &
                    -1.d0/BB2*(B_dot_grad_BB2__AZ__nn)*B_dot_U0 &
                    +(B_dot_grad__BdotU0_AZ__nn) &
                    )*(B_dot_grad_BpVstar__p) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BpVstar__AZ__nn)) &
                    )
                  Qjac_Up_visco_par_AZ_knn = ( &
                    +1.d0/BB2*( &
                    +((grad_B_AZ_dot_U__phi__nn) &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AZ__phi__nn))*(grad_BpVstar_phi_k)/R**2 &
                    ) &
                    )
                  Qjac_Up_visco_parpar_AZ_knn = ( &
                    +1.d0/BB2**2*(( &
                    -1.d0/BB2*(B_dot_grad_BB2__AZ__nn)*B_dot_U0 &
                    +(B_dot_grad__BdotU0_AZ__nn) &
                    )*(B_dot_grad_BpVstar__k) &
                    ) &
                    )

                  Qjac_Up_visco_par_AZ_for_Uphi_p = &
                    - v * 1.d0/BB2**2*(BB2_AZ__p)*( &
                     (grad_BdotU0_R-1.d0/BB2*B_dot_U0*BB2_R)*Bp0_R &
                    +(grad_BdotU0_Z-1.d0/BB2*B_dot_U0*BB2_Z)*Bp0_Z &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*Bp0_p/R**2) &
                    + v * 1.d0/BB2*( &
                     ((grad_B_AZ_dot_U__R__p) + 1.d0/BB2**2*(BB2_AZ__p)*B_dot_U0*BB2_R &
                    -1.d0/BB2*(B_AZ_dot_U__p)*BB2_R-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__R__p))*Bp0_R &
                    +(grad_BdotU0_R-1.d0/BB2*B_dot_U0*BB2_R)*grad_Bp_AZ__R &
                    +((grad_B_AZ_dot_U__Z__p)+1.d0/BB2**2*(BB2_AZ__p)*B_dot_U0*BB2_Z &
                    -1.d0/BB2*(B_AZ_dot_U__p)*BB2_Z-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__Z__p))*Bp0_Z &
                    +(grad_BdotU0_Z-1.d0/BB2*B_dot_U0*BB2_Z)*grad_Bp_AZ__Z &
                    +((grad_B_AZ_dot_U__phi__p)+1.d0/BB2**2*(BB2_AZ__p)*B_dot_U0*BB2_phi &
                    -1.d0/BB2*(B_AZ_dot_U__p)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__phi__p))*Bp0_p/R**2 &
                    ) 

                  Qjac_Up_visco_par_AZ_for_Uphi_n = &
                    - v * 1.d0/BB2**2*(BB2_AZ__n)*( &
                     (grad_BdotU0_R-1.d0/BB2*B_dot_U0*BB2_R)*Bp0_R &
                    +(grad_BdotU0_Z-1.d0/BB2*B_dot_U0*BB2_Z)*Bp0_Z &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*Bp0_p/R**2) &
                    + v * 1.d0/BB2*( &
                     ((grad_B_AZ_dot_U__R__n) + 1.d0/BB2**2*(BB2_AZ__n)*B_dot_U0*BB2_R &
                    -1.d0/BB2*(B_AZ_dot_U__n)*BB2_R-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__R__n))*Bp0_R &
                    +((grad_B_AZ_dot_U__Z__n)+1.d0/BB2**2*(BB2_AZ__n)*B_dot_U0*BB2_Z &
                    -1.d0/BB2*(B_AZ_dot_U__n)*BB2_Z-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__Z__n))*Bp0_Z &
                    +((grad_B_AZ_dot_U__phi__n)+1.d0/BB2**2*(BB2_AZ__n)*B_dot_U0*BB2_phi &
                    -1.d0/BB2*(B_AZ_dot_U__n)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_AZ__phi__n))*Bp0_p/R**2 &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*grad_Bp_AZ__phi__n/R**2) 

                  Qjac_Up_visco_par_AZ_for_Uphi_pnn = &
                    + v * 1.d0/BB2*( &
                     ( &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AZ__R__nn))*Bp0_R &
                    +( &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AZ__Z__nn))*Bp0_Z &
                    +((grad_B_AZ_dot_U__phi__nn) &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_AZ__phi__nn))*Bp0_p/R**2 &
                    ) 

                  Qjac_Up_visco_parpar_AZ_for_Uphi_p = &
                    - v * 2.d0/BB2**3*(BB2_AZ__p)*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*B_dot_grad_Bp &
                    + v * 1.d0/BB2**2*(( &
                     1.d0/BB2**2*(BB2_AZ__p)*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__AZ__p)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_AZ_dot_U__p)+(B_AZ_dot_grad_BdotU0__p)+(B_dot_grad__BdotU0_AZ__p) &
                    )*B_dot_grad_Bp &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_Bp__AZ__p)) 

                  Qjac_Up_visco_parpar_AZ_for_Uphi_n = &
                    - v * 2.d0/BB2**3*(BB2_AZ__n)*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*B_dot_grad_Bp &
                    + v * 1.d0/BB2**2*(( &
                     1.d0/BB2**2*(BB2_AZ__n)*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__AZ__n)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_AZ_dot_U__n)+(B_AZ_dot_grad_BdotU0__n)+(B_dot_grad__BdotU0_AZ__n) &
                    )*B_dot_grad_Bp &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_Bp__AZ__n))

                  Qjac_Up_visco_parpar_AZ_for_Uphi_pnn = &
                    + v * 1.d0/BB2**2*(( &
                    -1.d0/BB2*(B_dot_grad_BB2__AZ__nn)*B_dot_U0 &
                    +(B_dot_grad__BdotU0_AZ__nn) &
                    )*B_dot_grad_Bp &
                    ) 

                  ! --- var_Up,var_A3
                  Qjac_Up_visco_par_A3_p = ( &
                    -1.d0/BB2**2*BB2_A3*( &
                    (grad_BdotU0_R-1.d0/BB2*B_dot_U0*BB2_R) * grad_BpVstar_R &
                    +(grad_BdotU0_Z-1.d0/BB2*B_dot_U0*BB2_Z) * grad_BpVstar_Z &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi) * (grad_BpVstar_phi)/R**2) &
                    +1.d0/BB2*( &
                    ((grad_B_A3_dot_U__R__p)+1.d0/BB2**2*BB2_A3*B_dot_U0*BB2_R &
                    -1.d0/BB2*(B_A3_dot_U__p)*BB2_R-1.d0/BB2*B_dot_U0*(grad_BB2_A3__R__p))*grad_BpVstar_R &
                    +(grad_BdotU0_R-1.d0/BB2*B_dot_U0*BB2_R)*(v*grad_Bp_A3__R+Bp0_A3*v_R) &
                    +((grad_B_A3_dot_U__Z__p)+1.d0/BB2**2*BB2_A3*B_dot_U0*BB2_Z &
                    -1.d0/BB2*(B_A3_dot_U__p)*BB2_Z-1.d0/BB2*B_dot_U0*(grad_BB2_A3__Z__p))*grad_BpVstar_Z &
                    +(grad_BdotU0_Z-1.d0/BB2*B_dot_U0*BB2_Z)*(v*grad_Bp_A3__Z+Bp0_A3*v_Z) &
                    +((grad_B_A3_dot_U__phi__p)+1.d0/BB2**2*BB2_A3*B_dot_U0*BB2_phi &
                    -1.d0/BB2*(B_A3_dot_U__p)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_A3__phi__p))*(grad_BpVstar_phi)/R**2 &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*(v*grad_Bp_A3__phi)/R**2) &
                    )
                  Qjac_Up_visco_parpar_A3_p = ( &
                    -2.d0/BB2**3*BB2_A3*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BpVstar__p) &
                    +1.d0/BB2**2*(( &
                    1.d0/BB2**2*BB2_A3*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__A3__p)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_A3_dot_U__p)+(B_A3_dot_grad_BdotU0__p)+(B_dot_grad__BdotU0_A3__p) &
                    )*(B_dot_grad_BpVstar__p) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BpVstar__A3__p)) &
                    )
                  Qjac_Up_visco_par_A3_n = ( &
                    +1.d0/BB2*( &
                    ((grad_B_A3_dot_U__R__n) &
                    -1.d0/BB2*(B_A3_dot_U__n)*BB2_R-1.d0/BB2*B_dot_U0*(grad_BB2_A3__R__n))*grad_BpVstar_R &
                    +((grad_B_A3_dot_U__Z__n) &
                    -1.d0/BB2*(B_A3_dot_U__n)*BB2_Z-1.d0/BB2*B_dot_U0*(grad_BB2_A3__Z__n))*grad_BpVstar_Z &
                    +((grad_B_A3_dot_U__phi__n) &
                    -1.d0/BB2*(B_A3_dot_U__n)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_A3__phi__n))*(grad_BpVstar_phi)/R**2 &
                    ) &
                    )
                  Qjac_Up_visco_parpar_A3_n = ( &
                    +1.d0/BB2**2*(( &
                    -1.d0/BB2*(B_dot_grad_BB2__A3__n)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_A3_dot_U__n)+(B_A3_dot_grad_BdotU0__n)+(B_dot_grad__BdotU0_A3__n) &
                    )*(B_dot_grad_BpVstar__p) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BpVstar__A3__n)) &
                    )
                  Qjac_Up_visco_par_A3_k = ( &
                    -1.d0/BB2**2*BB2_A3*( &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi) * (grad_BpVstar_phi_k)/R**2) &
                    +1.d0/BB2*( &
                    +((grad_B_A3_dot_U__phi__p)+1.d0/BB2**2*BB2_A3*B_dot_U0*BB2_phi &
                    -1.d0/BB2*(B_A3_dot_U__p)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_A3__phi__p))*(grad_BpVstar_phi_k)/R**2 &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*(Bp0_A3*v_p__k)/R**2) &
                    )
                  Qjac_Up_visco_parpar_A3_k = ( &
                    -2.d0/BB2**3*BB2_A3*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BpVstar__k) &
                    +1.d0/BB2**2*(( &
                    1.d0/BB2**2*BB2_A3*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__A3__p)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_A3_dot_U__p)+(B_A3_dot_grad_BdotU0__p)+(B_dot_grad__BdotU0_A3__p) &
                    )*(B_dot_grad_BpVstar__k) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BpVstar__A3__k)) &
                    )
                  Qjac_Up_visco_par_A3_kn = ( &
                    +1.d0/BB2*( &
                    +((grad_B_A3_dot_U__phi__n) &
                    -1.d0/BB2*(B_A3_dot_U__n)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_A3__phi__n))*(grad_BpVstar_phi_k)/R**2 &
                    ) &
                    )
                  Qjac_Up_visco_parpar_A3_kn = ( &
                    +1.d0/BB2**2*(( &
                    -1.d0/BB2*(B_dot_grad_BB2__A3__n)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_A3_dot_U__n)+(B_A3_dot_grad_BdotU0__n)+(B_dot_grad__BdotU0_A3__n) &
                    )*(B_dot_grad_BpVstar__k) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BpVstar__A3__kn)) &
                    )
                  Qjac_Up_visco_par_A3_pnn = ( &
                    +1.d0/BB2*( &
                    ( &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_A3__R__nn))*grad_BpVstar_R &
                    +( &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_A3__Z__nn))*grad_BpVstar_Z &
                    +((grad_B_A3_dot_U__phi__nn) &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_A3__phi__nn))*(grad_BpVstar_phi)/R**2 &
                    ) &
                    )
                  Qjac_Up_visco_parpar_A3_pnn = ( &
                    +1.d0/BB2**2*(( &
                    -1.d0/BB2*(B_dot_grad_BB2__A3__nn)*B_dot_U0 &
                    +(B_dot_grad__BdotU0_A3__nn) &
                    )*(B_dot_grad_BpVstar__p) &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*(B_dot_grad_BpVstar__A3__nn)) &
                    )
                  Qjac_Up_visco_par_A3_knn = ( &
                    +1.d0/BB2*( &
                    +((grad_B_A3_dot_U__phi__nn) &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_A3__phi__nn))*(grad_BpVstar_phi_k)/R**2 &
                    ) &
                    )
                  Qjac_Up_visco_parpar_A3_knn = ( &
                    +1.d0/BB2**2*(( &
                    -1.d0/BB2*(B_dot_grad_BB2__A3__nn)*B_dot_U0 &
                    +(B_dot_grad__BdotU0_A3__nn) &
                    )*(B_dot_grad_BpVstar__k) &
                    ) &
                    )

                  Qjac_Up_visco_par_A3_for_Uphi_p = &
                    - v * 1.d0/BB2**2*BB2_A3*( &
                    (grad_BdotU0_R-1.d0/BB2*B_dot_U0*BB2_R) * Bp0_R &
                    +(grad_BdotU0_Z-1.d0/BB2*B_dot_U0*BB2_Z) * Bp0_Z &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi) * Bp0_p/R**2) &
                    + v * 1.d0/BB2*( &
                    ((grad_B_A3_dot_U__R__p)+1.d0/BB2**2*BB2_A3*B_dot_U0*BB2_R &
                    -1.d0/BB2*(B_A3_dot_U__p)*BB2_R-1.d0/BB2*B_dot_U0*(grad_BB2_A3__R__p))*Bp0_R &
                    +(grad_BdotU0_R-1.d0/BB2*B_dot_U0*BB2_R)*grad_Bp_A3__R &
                    +((grad_B_A3_dot_U__Z__p)+1.d0/BB2**2*BB2_A3*B_dot_U0*BB2_Z &
                    -1.d0/BB2*(B_A3_dot_U__p)*BB2_Z-1.d0/BB2*B_dot_U0*(grad_BB2_A3__Z__p))*Bp0_Z &
                    +(grad_BdotU0_Z-1.d0/BB2*B_dot_U0*BB2_Z)*grad_Bp_A3__Z &
                    +((grad_B_A3_dot_U__phi__p)+1.d0/BB2**2*BB2_A3*B_dot_U0*BB2_phi &
                    -1.d0/BB2*(B_A3_dot_U__p)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_A3__phi__p))*Bp0_p/R**2 &
                    +(grad_BdotU0_phi-1.d0/BB2*B_dot_U0*BB2_phi)*grad_Bp_A3__phi/R**2) 

                  Qjac_Up_visco_par_A3_for_Uphi_n = &
                    + v * 1.d0/BB2*( &
                    ((grad_B_A3_dot_U__R__n) &
                    -1.d0/BB2*(B_A3_dot_U__n)*BB2_R-1.d0/BB2*B_dot_U0*(grad_BB2_A3__R__n))*Bp0_R &
                    +((grad_B_A3_dot_U__Z__n) &
                    -1.d0/BB2*(B_A3_dot_U__n)*BB2_Z-1.d0/BB2*B_dot_U0*(grad_BB2_A3__Z__n))*Bp0_Z &
                    +((grad_B_A3_dot_U__phi__n) &
                    -1.d0/BB2*(B_A3_dot_U__n)*BB2_phi-1.d0/BB2*B_dot_U0*(grad_BB2_A3__phi__n))*Bp0_p/R**2 &
                    ) 

                  Qjac_Up_visco_par_A3_for_Uphi_pnn = &
                    + v * 1.d0/BB2*( &
                    ( &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_A3__R__nn))*Bp0_R &
                    +( &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_A3__Z__nn))*Bp0_Z &
                    +((grad_B_A3_dot_U__phi__nn) &
                    -1.d0/BB2*B_dot_U0*(grad_BB2_A3__phi__nn))*Bp0_p/R**2 &
                    ) 

                  Qjac_Up_visco_parpar_A3_for_Uphi_p = &
                    - v * 2.d0/BB2**3*BB2_A3*(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*B_dot_grad_Bp &
                    + v * 1.d0/BB2**2*(( &
                    1.d0/BB2**2*BB2_A3*B_dot_grad_BB2*B_dot_U0-1.d0/BB2*(B_dot_grad_BB2__A3__p)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_A3_dot_U__p)+(B_A3_dot_grad_BdotU0__p)+(B_dot_grad__BdotU0_A3__p) &
                    )*B_dot_grad_Bp &
                    +(-1.d0/BB2*B_dot_grad_BB2*B_dot_U0+B_dot_grad_BdotU0)*B_dot_grad_Bp__A3__p) 

                  Qjac_Up_visco_parpar_A3_for_Uphi_n = &
                    + v * 1.d0/BB2**2*(( &
                    -1.d0/BB2*(B_dot_grad_BB2__A3__n)*B_dot_U0 &
                    -1.d0/BB2*B_dot_grad_BB2*(B_A3_dot_U__n)+(B_A3_dot_grad_BdotU0__n)+(B_dot_grad__BdotU0_A3__n) &
                    )*B_dot_grad_Bp &
                    ) 

                  Qjac_Up_visco_parpar_A3_for_Uphi_pnn = &
                    + v * 1.d0/BB2**2*(( &
                    -1.d0/BB2*(B_dot_grad_BB2__A3__nn)*B_dot_U0 &
                    +(B_dot_grad__BdotU0_A3__nn) &
                    )*B_dot_grad_Bp &
                    ) 

                  ! --- var_Up,var_UR
                  Qjac_Up_visco_par_UR_p = ( &
                    1.d0/BB2*( &
                    (BR0_R*UR+BR0*UR_R-1.d0/BB2*BR0*UR*BB2_R)*(v*Bp0_R+Bp0*v_R) &
                    +(BR0_Z*UR+BR0*UR_Z-1.d0/BB2*BR0*UR*BB2_Z)*(v*Bp0_Z+Bp0*v_Z) &
                    +(BR0_p*UR-1.d0/BB2*BR0*UR*BB2_phi)*(v*Bp0_p)/R**2) &
                    )
                  Qjac_Up_visco_parpar_UR_p = ( &
                    1.d0/BB2**2*(-BR0/BB2*B_dot_grad_BB2*UR+B_dot_grad_BR*UR+BR0*(B_dot_grad__UR__p)) &
                    *(B_dot_grad_BpVstar__p) &
                    )
                  Qjac_Up_visco_par_UR_n = ( &
                    1.d0/BB2*( &
                    +(BR0*UR_p__n)*(v*Bp0_p)/R**2) &
                    )
                  Qjac_Up_visco_parpar_UR_n = ( &
                    1.d0/BB2**2*(BR0*(B_dot_grad__UR__n)) &
                    *(B_dot_grad_BpVstar__p) &
                    )
                  Qjac_Up_visco_par_UR_k = ( &
                    1.d0/BB2*( &
                    +(BR0_p*UR-1.d0/BB2*BR0*UR*BB2_phi)*(Bp0*v_p__k)/R**2) &
                    )
                  Qjac_Up_visco_parpar_UR_k = ( &
                    1.d0/BB2**2*(-BR0/BB2*B_dot_grad_BB2*UR+B_dot_grad_BR*UR+BR0*(B_dot_grad__UR__p)) &
                    *(B_dot_grad_BpVstar__k) &
                    )
                  Qjac_Up_visco_par_UR_kn = ( &
                    1.d0/BB2*( &
                    +(BR0*UR_p__n)*(Bp0*v_p__k)/R**2) &
                    )
                  Qjac_Up_visco_parpar_UR_kn = ( &
                    1.d0/BB2**2*(BR0*(B_dot_grad__UR__n)) &
                    *(B_dot_grad_BpVstar__k) &
                    )
                  Qjac_Up_visco_par_UR_pnn = 0.d0
                  Qjac_Up_visco_parpar_UR_pnn = 0.d0
                  Qjac_Up_visco_par_UR_knn = 0.d0
                  Qjac_Up_visco_parpar_UR_knn = 0.d0

                  Qjac_Up_visco_par_UR_for_Uphi_p = &
                    v * 1.d0/BB2*( &
                    (BR0_R*UR+BR0*UR_R-1.d0/BB2*BR0*UR*BB2_R)*Bp0_R &
                    +(BR0_Z*UR+BR0*UR_Z-1.d0/BB2*BR0*UR*BB2_Z)*Bp0_Z &
                    +(BR0_p*UR-1.d0/BB2*BR0*UR*BB2_phi)*Bp0_p/R**2) 

                  Qjac_Up_visco_par_UR_for_Uphi_n = &
                    v * 1.d0/BB2*( &
                    +(BR0*UR_p__n)*Bp0_p/R**2) 

                  Qjac_Up_visco_par_UR_for_Uphi_pnn = 0.d0

                  Qjac_Up_visco_parpar_UR_for_Uphi_p = &
                    v * 1.d0/BB2**2*(-BR0/BB2*B_dot_grad_BB2*UR+B_dot_grad_BR*UR+BR0*(B_dot_grad__UR__p)) &
                    *B_dot_grad_Bp

                  Qjac_Up_visco_parpar_UR_for_Uphi_n = &
                    v * 1.d0/BB2**2*(BR0*(B_dot_grad__UR__n)) &
                    *B_dot_grad_Bp

                  Qjac_Up_visco_parpar_UR_for_Uphi_pnn = 0.d0

                  ! --- var_Up,var_UZ
                  Qjac_Up_visco_par_UZ_p = ( &
                    1.d0/BB2*( &
                    (BZ0_R*UZ+BZ0*UZ_R-1.d0/BB2*BZ0*UZ*BB2_R)*(v*Bp0_R+Bp0*v_R) &
                    +(BZ0_Z*UZ+BZ0*UZ_Z-1.d0/BB2*BZ0*UZ*BB2_Z)*(v*Bp0_Z+Bp0*v_Z) &
                    +(BZ0_p*UZ-1.d0/BB2*BZ0*UZ*BB2_phi)*(v*Bp0_p)/R**2) &
                    )
                  Qjac_Up_visco_parpar_UZ_p = ( &
                    1.d0/BB2**2*(-BZ0/BB2*B_dot_grad_BB2*UZ+B_dot_grad_BZ*UZ+BZ0*(B_dot_grad__UZ__p)) &
                    *(B_dot_grad_BpVstar__p) &
                    )
                  Qjac_Up_visco_par_UZ_n = ( &
                    1.d0/BB2*( &
                    +(BZ0*UZ_p__n)*(v*Bp0_p)/R**2) &
                    )
                  Qjac_Up_visco_parpar_UZ_n = ( &
                    1.d0/BB2**2*(BZ0*(B_dot_grad__UZ__n)) &
                    *(B_dot_grad_BpVstar__p) &
                    )
                  Qjac_Up_visco_par_UZ_k = ( &
                    1.d0/BB2*( &
                    +(BZ0_p*UZ-1.d0/BB2*BZ0*UZ*BB2_phi)*(Bp0*v_p__k)/R**2) &
                    )
                  Qjac_Up_visco_parpar_UZ_k = ( &
                    1.d0/BB2**2*(-BZ0/BB2*B_dot_grad_BB2*UZ+B_dot_grad_BZ*UZ+BZ0*(B_dot_grad__UZ__p)) &
                    *(B_dot_grad_BpVstar__k) &
                    )
                  Qjac_Up_visco_par_UZ_kn = ( &
                    1.d0/BB2*( &
                    +(BZ0*UZ_p__n)*(Bp0*v_p__k)/R**2) &
                    )
                  Qjac_Up_visco_parpar_UZ_kn = ( &
                    1.d0/BB2**2*(BZ0*(B_dot_grad__UZ__n)) &
                    *(B_dot_grad_BpVstar__k) &
                    )
                  Qjac_Up_visco_par_UZ_pnn = 0.d0
                  Qjac_Up_visco_parpar_UZ_pnn = 0.d0
                  Qjac_Up_visco_par_UZ_knn = 0.d0
                  Qjac_Up_visco_parpar_UZ_knn = 0.d0

                  Qjac_Up_visco_par_UZ_for_Uphi_p = &
                    v * 1.d0/BB2*( &
                    (BZ0_R*UZ+BZ0*UZ_R-1.d0/BB2*BZ0*UZ*BB2_R)*Bp0_R &
                    +(BZ0_Z*UZ+BZ0*UZ_Z-1.d0/BB2*BZ0*UZ*BB2_Z)*Bp0_Z &
                    +(BZ0_p*UZ-1.d0/BB2*BZ0*UZ*BB2_phi)*Bp0_p/R**2) 

                  Qjac_Up_visco_par_UZ_for_Uphi_n = &
                    v * 1.d0/BB2*( &
                    +(BZ0*UZ_p__n)*Bp0_p/R**2) 

                  Qjac_Up_visco_par_UZ_for_Uphi_pnn = 0.d0

                  Qjac_Up_visco_parpar_UZ_for_Uphi_p = &
                    v * 1.d0/BB2**2*(-BZ0/BB2*B_dot_grad_BB2*UZ+B_dot_grad_BZ*UZ+BZ0*(B_dot_grad__UZ__p)) &
                    *B_dot_grad_Bp

                  Qjac_Up_visco_parpar_UZ_for_Uphi_n = &
                    v * 1.d0/BB2**2*(BZ0*(B_dot_grad__UZ__n)) &
                    *B_dot_grad_Bp

                  Qjac_Up_visco_parpar_UZ_for_Uphi_pnn = 0.d0

                  ! --- var_Up,var_Up
                  Qjac_Up_visco_par_Up_p = ( &
                    1.d0/BB2*( &
                    (Bp0_R*Up+Bp0*Up_R-1.d0/BB2*Bp0*Up*BB2_R)*(v*Bp0_R+Bp0*v_R) &
                    +(Bp0_Z*Up+Bp0*Up_Z-1.d0/BB2*Bp0*Up*BB2_Z)*(v*Bp0_Z+Bp0*v_Z) &
                    +(Bp0_p*Up-1.d0/BB2*Bp0*Up*BB2_phi)*(v*Bp0_p)/R**2) &
                    )
                  Qjac_Up_visco_parpar_Up_p = ( &
                    1.d0/BB2**2*(-Bp0/BB2*B_dot_grad_BB2*Up+B_dot_grad_Bp*Up+Bp0*(B_dot_grad__Up__p)) &
                    *(B_dot_grad_BpVstar__p) &
                    )
                  Qjac_Up_visco_par_Up_n = ( &
                    1.d0/BB2*( &
                    +(Bp0*Up_p__n)*(v*Bp0_p)/R**2) &
                    )
                  Qjac_Up_visco_parpar_Up_n = ( &
                    1.d0/BB2**2*(Bp0*(B_dot_grad__Up__n)) &
                    *(B_dot_grad_BpVstar__p) &
                    )
                  Qjac_Up_visco_par_Up_k = ( &
                    1.d0/BB2*( &
                    +(Bp0_p*Up-1.d0/BB2*Bp0*Up*BB2_phi)*(Bp0*v_p__k)/R**2) &
                    )
                  Qjac_Up_visco_parpar_Up_k = ( &
                    1.d0/BB2**2*(-Bp0/BB2*B_dot_grad_BB2*Up+B_dot_grad_Bp*Up+Bp0*(B_dot_grad__Up__p)) &
                    *(B_dot_grad_BpVstar__k) &
                    )
                  Qjac_Up_visco_par_Up_kn = ( &
                    1.d0/BB2*( &
                    +(Bp0*Up_p__n)*(Bp0*v_p__k)/R**2) &
                    )
                  Qjac_Up_visco_parpar_Up_kn = ( &
                    1.d0/BB2**2*(Bp0*(B_dot_grad__Up__n)) &
                    *(B_dot_grad_BpVstar__k) &
                    )
                  Qjac_Up_visco_par_Up_pnn = 0.d0
                  Qjac_Up_visco_parpar_Up_pnn = 0.d0
                  Qjac_Up_visco_par_Up_knn = 0.d0
                  Qjac_Up_visco_parpar_Up_knn = 0.d0

                  Qjac_Up_visco_par_Up_for_Uphi_p = &
                    v * 1.d0/BB2*( &
                    (Bp0_R*Up+Bp0*Up_R-1.d0/BB2*Bp0*Up*BB2_R)*Bp0_R &
                    +(Bp0_Z*Up+Bp0*Up_Z-1.d0/BB2*Bp0*Up*BB2_Z)*Bp0_Z &
                    +(Bp0_p*Up-1.d0/BB2*Bp0*Up*BB2_phi)*Bp0_p/R**2) 

                  Qjac_Up_visco_par_Up_for_Uphi_n = &
                    v * 1.d0/BB2*( &
                    +(Bp0*Up_p__n)*Bp0_p/R**2) 

                  Qjac_Up_visco_par_Up_for_Uphi_pnn = 0.d0

                  Qjac_Up_visco_parpar_Up_for_Uphi_p = &
                    v * 1.d0/BB2**2*(-Bp0/BB2*B_dot_grad_BB2*Up+B_dot_grad_Bp*Up+Bp0*(B_dot_grad__Up__p)) &
                    *B_dot_grad_Bp
                  
                  Qjac_Up_visco_parpar_Up_for_Uphi_n = &
                    v * 1.d0/BB2**2*(Bp0*(B_dot_grad__Up__n)) &
                    *B_dot_grad_Bp   
                    
                  Qjac_Up_visco_parpar_Up_for_Uphi_pnn = 0.d0
#endif

                  ! --- Viscous terms
                  Qvisc_UR_UR__p = 0.d0 ; Qvisc_UR_UR__n = 0.d0 ; Qvisc_UR_UR__k = 0.d0 ; Qvisc_UR_UR__kn = 0.d0
                  Qvisc_UR_UZ__p = 0.d0 ; Qvisc_UR_UZ__n = 0.d0 ; Qvisc_UR_UZ__k = 0.d0 ; Qvisc_UR_UZ__kn = 0.d0
                  Qvisc_UR_Up__p = 0.d0 ; Qvisc_UR_Up__n = 0.d0 ; Qvisc_UR_Up__k = 0.d0 ; Qvisc_UR_Up__kn = 0.d0
                  Qvisc_UZ_UR__p = 0.d0 ; Qvisc_UZ_UR__n = 0.d0 ; Qvisc_UZ_UR__k = 0.d0 ; Qvisc_UZ_UR__kn = 0.d0
                  Qvisc_UZ_UZ__p = 0.d0 ; Qvisc_UZ_UZ__n = 0.d0 ; Qvisc_UZ_UZ__k = 0.d0 ; Qvisc_UZ_UZ__kn = 0.d0
                  Qvisc_UZ_Up__p = 0.d0 ; Qvisc_UZ_Up__n = 0.d0 ; Qvisc_UZ_Up__k = 0.d0 ; Qvisc_UZ_Up__kn = 0.d0
                  Qvisc_Up_UR__p = 0.d0 ; Qvisc_Up_UR__n = 0.d0 ; Qvisc_Up_UR__k = 0.d0 ; Qvisc_Up_UR__kn = 0.d0
                  Qvisc_Up_UZ__p = 0.d0 ; Qvisc_Up_UZ__n = 0.d0 ; Qvisc_Up_UZ__k = 0.d0 ; Qvisc_Up_UZ__kn = 0.d0
                  Qvisc_Up_Up__p = 0.d0 ; Qvisc_Up_Up__n = 0.d0 ; Qvisc_Up_Up__k = 0.d0 ; Qvisc_Up_Up__kn = 0.d0
                  Qvisc_T_UR__p  = 0.d0 ; Qvisc_T_UR__n  = 0.d0
                  Qvisc_T_UZ__p  = 0.d0 ; Qvisc_T_UZ__n  = 0.d0
                  Qvisc_T_Up__p  = 0.d0 ; Qvisc_T_Up__n  = 0.d0
                  Qvisc_T_T__p   = 0.d0 ; Qvisc_T_T__n   = 0.d0

                  Qgrad_visc_dot_grad_UR__UR__p = 0.d0; Qgrad_visc_dot_grad_UR__UR__n = 0.d0
                  Qgrad_visc_dot_grad_UZ__UZ__p = 0.d0; Qgrad_visc_dot_grad_UZ__UZ__n = 0.d0
                  Qgrad_visc_dot_grad_Up__Up__p = 0.d0; Qgrad_visc_dot_grad_Up__Up__n = 0.d0
                  Qgrad_visc_dot_grad_UR__Te__p = 0.d0; Qgrad_visc_dot_grad_UR__Te__n = 0.d0
                  Qgrad_visc_dot_grad_UZ__Te__p = 0.d0; Qgrad_visc_dot_grad_UZ__Te__n = 0.d0
                  Qgrad_visc_dot_grad_Up__Te__p = 0.d0; Qgrad_visc_dot_grad_Up__Te__n = 0.d0
                  
                  Qgrad_visc_dot_grad_UR__UZ__p = 0.d0; Qgrad_visc_dot_grad_UR__Up__p = 0.d0; Qgrad_visc_dot_grad_UR__Up__n = 0.d0;
                  Qgrad_visc_dot_grad_UZ__UR__p = 0.d0; Qgrad_visc_dot_grad_UZ__Up__p = 0.d0; Qgrad_visc_dot_grad_UZ__Up__n = 0.d0;
                  Qgrad_visc_dot_grad_Up__UR__p = 0.d0; Qgrad_visc_dot_grad_Up__UR__n = 0.d0; Qgrad_visc_dot_grad_Up__UZ__p = 0.d0; Qgrad_visc_dot_grad_Up__UZ__n = 0.d0
                  
                  QdivU_BdotGradVisc__UR__p = 0.d0; QdivU_BdotGradVisc__UZ__p = 0.d0; QdivU_BdotGradVisc__Up__n = 0.d0
                  QdivU_BdotGradVisc__AR__p = 0.d0
                  QdivU_BdotGradVisc__AR__n = 0.d0
                  QdivU_BdotGradVisc__AZ__p = 0.d0
                  QdivU_BdotGradVisc__AZ__n = 0.d0
                  QdivU_BdotGradVisc__A3__p = 0.d0
                  QdivU_BdotGradVisc__Te__p = 0.d0
                  QdivU_BdotGradVisc__Te__n = 0.d0
                  
                  curl_U0_dot_gradmu_x_B0__UR__p = 0.d0; curl_U0_dot_gradmu_x_B0__UR__n = 0.d0; curl_U0_dot_gradmu_x_B0__UZ__p = 0.d0; curl_U0_dot_gradmu_x_B0__UZ__n = 0.d0; curl_U0_dot_gradmu_x_B0__Up__p = 0.d0
                  curl_U0_dot_gradmu_x_B0__AR__p = 0.d0; curl_U0_dot_gradmu_x_B0__AR__n = 0.d0; curl_U0_dot_gradmu_x_B0__AZ__p = 0.d0; curl_U0_dot_gradmu_x_B0__AZ__n = 0.d0; curl_U0_dot_gradmu_x_B0__A3__p = 0.d0
                  curl_U0_dot_gradmu_x_B0__Te__p = 0.d0; curl_U0_dot_gradmu_x_B0__Te__n = 0.d0;
                  
                  curl_U0_dot_curl_B0__UR__p = 0.d0; curl_U0_dot_curl_B0__UR__n = 0.d0
                  curl_U0_dot_curl_B0__UZ__p = 0.d0; curl_U0_dot_curl_B0__UZ__n = 0.d0
                  curl_U0_dot_curl_B0__Up__p = 0.d0;
                  
                  curl_U0_dot_curl_B0__AR__p = 0.d0; curl_U0_dot_curl_B0__AR__n = 0.d0; curl_U0_dot_curl_B0__AR__nn = 0.d0
                  curl_U0_dot_curl_B0__AZ__p = 0.d0; curl_U0_dot_curl_B0__AZ__n = 0.d0; curl_U0_dot_curl_B0__AZ__nn = 0.d0
                  curl_U0_dot_curl_B0__A3__p = 0.d0; curl_U0_dot_curl_B0__A3__n = 0.d0
                  
                  Qvisc_grad_BR_dot_grad_UR__UR__p = 0.d0; Qvisc_grad_BR_dot_grad_UR__UR__n = 0.d0
                  Qvisc_grad_BZ_dot_grad_UZ__UZ__p = 0.d0; Qvisc_grad_BZ_dot_grad_UZ__UZ__n = 0.d0
                  Qvisc_grad_Bp_dot_grad_Up__Up__p = 0.d0; Qvisc_grad_Bp_dot_grad_Up__Up__n = 0.d0
                  Qvisc_grad_BR_dot_grad_UR__AR__p = 0.d0;
                  Qvisc_grad_BR_dot_grad_UR__AZ__n = 0.d0; Qvisc_grad_BR_dot_grad_UR__AZ__nn = 0.d0
                  Qvisc_grad_BR_dot_grad_UR__A3__p = 0.d0; Qvisc_grad_BR_dot_grad_UR__A3__n = 0.d0
                  Qvisc_grad_BZ_dot_grad_UZ__AR__p = 0.d0; Qvisc_grad_BZ_dot_grad_UZ__AR__nn = 0.d0
                  Qvisc_grad_BZ_dot_grad_UZ__AZ__p = 0.d0;
                  Qvisc_grad_BZ_dot_grad_UZ__A3__p = 0.d0; Qvisc_grad_BZ_dot_grad_UZ__A3__n = 0.d0
                  Qvisc_grad_Bp_dot_grad_Up__AR__p = 0.d0; Qvisc_grad_Bp_dot_grad_Up__AR__n = 0.d0
                  Qvisc_grad_Bp_dot_grad_Up__AZ__p = 0.d0; Qvisc_grad_Bp_dot_grad_Up__AZ__n = 0.d0
                  Qvisc_grad_Bp_dot_grad_Up__A3__p = 0.d0;

                  Qvisc_gradmu_dot_grad_vecV_UR__UR__p = 0.d0
                  Qvisc_gradmu_dot_grad_vecV_UR__UR__n = 0.d0
                  Qvisc_gradmu_dot_grad_vecV_UR__Up__p = 0.d0
                  Qvisc_gradmu_dot_grad_vecV_UR__Te__p = 0.d0
                  Qvisc_gradmu_dot_grad_vecV_UR__Te__n = 0.d0

                  Qvisc_gradmu_dot_grad_vecV_UZ__UZ__p = 0.d0
                  Qvisc_gradmu_dot_grad_vecV_UZ__UZ__n = 0.d0
                  Qvisc_gradmu_dot_grad_vecV_UZ__Te__p = 0.d0
                  Qvisc_gradmu_dot_grad_vecV_UZ__Te__n = 0.d0

                  Qvisc_gradmu_dot_grad_vecV_Up__UR__p = 0.d0
                  Qvisc_gradmu_dot_grad_vecV_Up__Up__p = 0.d0
                  Qvisc_gradmu_dot_grad_vecV_Up__Up__n = 0.d0
                  Qvisc_gradmu_dot_grad_vecV_Up__Te__p = 0.d0
                  Qvisc_gradmu_dot_grad_vecV_Up__Te__n = 0.d0

                  Qgrad_visc_dot_gradU_T_UR__UR__p = 0.d0
                  Qgrad_visc_dot_gradU_T_UR__UZ__p = 0.d0
                  Qgrad_visc_dot_gradU_T_UR__Up__p = 0.d0
                  Qgrad_visc_dot_gradU_T_UR__Te__p = 0.d0
                  Qgrad_visc_dot_gradU_T_UR__Te__n = 0.d0

                  Qgrad_visc_dot_gradU_T_UZ__UR__p = 0.d0
                  Qgrad_visc_dot_gradU_T_UZ__UZ__p = 0.d0
                  Qgrad_visc_dot_gradU_T_UZ__Up__p = 0.d0
                  Qgrad_visc_dot_gradU_T_UZ__Te__p = 0.d0
                  Qgrad_visc_dot_gradU_T_UZ__Te__n = 0.d0

                  Qgrad_visc_dot_gradU_T_Up__UR__p = 0.d0
                  Qgrad_visc_dot_gradU_T_Up__UR__n = 0.d0
                  Qgrad_visc_dot_gradU_T_Up__UZ__n = 0.d0
                  Qgrad_visc_dot_gradU_T_Up__Up__p = 0.d0
                  Qgrad_visc_dot_gradU_T_Up__Up__n = 0.d0
                  Qgrad_visc_dot_gradU_T_Up__Te__p = 0.d0
                  Qgrad_visc_dot_gradU_T_Up__Te__n = 0.d0

                  Q_minus_2_over_3_divU_grad_visc_UR__UR__p = 0.d0
                  Q_minus_2_over_3_divU_grad_visc_UR__UZ__p = 0.d0
                  Q_minus_2_over_3_divU_grad_visc_UR__Up__n = 0.d0
                  Q_minus_2_over_3_divU_grad_visc_UR__Te__p = 0.d0

                  Q_minus_2_over_3_divU_grad_visc_UZ__UR__p = 0.d0
                  Q_minus_2_over_3_divU_grad_visc_UZ__UZ__p = 0.d0
                  Q_minus_2_over_3_divU_grad_visc_UZ__Up__n = 0.d0
                  Q_minus_2_over_3_divU_grad_visc_UZ__Te__p = 0.d0

                  Q_minus_2_over_3_divU_grad_visc_Up__UR__p = 0.d0
                  Q_minus_2_over_3_divU_grad_visc_Up__UZ__p = 0.d0
                  Q_minus_2_over_3_divU_grad_visc_Up__Up__n = 0.d0
                  Q_minus_2_over_3_divU_grad_visc_Up__Te__p = 0.d0
                  Q_minus_2_over_3_divU_grad_visc_Up__Te__n = 0.d0

                  Q_1_over_3_visc_grad_divU_UR__UR__p = 0.d0
                  Q_1_over_3_visc_grad_divU_UR__UZ__p = 0.d0
                  Q_1_over_3_visc_grad_divU_UR__Up__n = 0.d0
                  Q_1_over_3_visc_grad_divU_UR__Te__p = 0.d0

                  Q_1_over_3_visc_grad_divU_UZ__UR__p = 0.d0
                  Q_1_over_3_visc_grad_divU_UZ__UZ__p = 0.d0
                  Q_1_over_3_visc_grad_divU_UZ__Up__n = 0.d0
                  Q_1_over_3_visc_grad_divU_UZ__Te__p = 0.d0

                  Q_1_over_3_visc_grad_divU_Up__UR__p = 0.d0
                  Q_1_over_3_visc_grad_divU_Up__UZ__p = 0.d0
                  Q_1_over_3_visc_grad_divU_Up__Up__n = 0.d0
                  Q_1_over_3_visc_grad_divU_Up__Te__p = 0.d0
                  Q_1_over_3_visc_grad_divU_Up__Te__n = 0.d0
                  Q_1_over_3_visc_grad_divU_Up__UR__k = 0.d0
                  Q_1_over_3_visc_grad_divU_Up__UZ__k = 0.d0
                  Q_1_over_3_visc_grad_divU_Up__Up__kn = 0.d0
                  Q_1_over_3_visc_grad_divU_Up__Te__k = 0.d0

                  Q_1_over_3_B_dot_visc_grad_divU__AR__p = 0.d0
                  Q_1_over_3_B_dot_visc_grad_divU__AR__n = 0.d0
                  Q_1_over_3_B_dot_visc_grad_divU__AZ__p = 0.d0
                  Q_1_over_3_B_dot_visc_grad_divU__AZ__n = 0.d0
                  Q_1_over_3_B_dot_visc_grad_divU__A3__p = 0.d0
                  Q_1_over_3_B_dot_visc_grad_divU__UR__p = 0.d0
                  Q_1_over_3_B_dot_visc_grad_divU__UZ__p = 0.d0
                  Q_1_over_3_B_dot_visc_grad_divU__Up__n = 0.d0
                  Q_1_over_3_B_dot_visc_grad_divU__Te__p = 0.d0
                  Q_1_over_3_B_dot_visc_grad_divU__Te__n = 0.d0

                  Q_1_over_3_B_dot_visc_grad_divU__AR__k = 0.d0
                  Q_1_over_3_B_dot_visc_grad_divU__AZ__k = 0.d0
                  Q_1_over_3_B_dot_visc_grad_divU__A3__k = 0.d0
                  Q_1_over_3_B_dot_visc_grad_divU__UR__k = 0.d0
                  Q_1_over_3_B_dot_visc_grad_divU__UZ__k = 0.d0
                  Q_1_over_3_B_dot_visc_grad_divU__Up__kn = 0.d0
                  Q_1_over_3_B_dot_visc_grad_divU__Te__k = 0.d0

#ifdef VISCOSITY_TENSOR_FORM                  
                  ! --- The correct viscosity (there is only one...)
                  Qvisc_UR_UR__p  = - (v_R + v/R) * divU_UR + v_Z * (- UR_Z)
                  Qvisc_UR_UR__k  = 0.d0
                  Qvisc_UR_UR__n  = 0.d0
                  Qvisc_UR_UR__kn = - v_p__k/R**2 * (UR_p__n)

                  Qvisc_UR_UZ__p  = - (v_R + v/R) * divU_UZ + v_Z * (UZ_R)
                  Qvisc_UR_UZ__k  = 0.d0
                  Qvisc_UR_UZ__n  = 0.d0
                  Qvisc_UR_UZ__kn = 0.d0

                  Qvisc_UR_Up__p  = 0.d0
                  Qvisc_UR_Up__k  = - v_p__k/R**2 * (- Up - R*Up_R)
                  Qvisc_UR_Up__n  = - (v_R + v/R) * divU_Up__n
                  Qvisc_UR_Up__kn = 0.d0

                  Qvisc_UZ_UR__p  = - v_Z * divU_UR - v_R * (- UR_Z)
                  Qvisc_UZ_UR__k  = 0.d0
                  Qvisc_UZ_UR__n  = 0.d0
                  Qvisc_UZ_UR__kn = 0.d0

                  Qvisc_UZ_UZ__p  = - v_Z * divU_UZ - v_R * (UZ_R)
                  Qvisc_UZ_UZ__k  = 0.d0
                  Qvisc_UZ_UZ__n  = 0.d0
                  Qvisc_UZ_UZ__kn = + v_p__k/R**2 * (- UZ_p__n)

                  Qvisc_UZ_Up__p  = 0.d0
                  Qvisc_UZ_Up__k  = + v_p__k/R**2 * (R*Up_Z)
                  Qvisc_UZ_Up__n  = - v_Z * divU_Up__n
                  Qvisc_UZ_Up__kn = 0.d0

                  Qvisc_Up_UR__p  = 0.d0
                  Qvisc_Up_UR__k  = - v_p__k/R * divU_UR
                  Qvisc_Up_UR__n  = + (v_R + v/R)/R * (UR_p__n)
                  Qvisc_Up_UR__kn = 0.d0

                  Qvisc_Up_UZ__p  = 0.d0
                  Qvisc_Up_UZ__k  = - v_p__k/R * divU_UZ
                  Qvisc_Up_UZ__n  = - v_Z/R * (- UZ_p__n)
                  Qvisc_Up_UZ__kn = 0.d0

                  Qvisc_Up_Up__p  = - v_Z/R * (R*Up_Z) + (v_R + v/R)/R * (- Up - R*Up_R)
                  Qvisc_Up_Up__k  = 0.d0
                  Qvisc_Up_Up__n  = 0.d0
                  Qvisc_Up_Up__kn = - v_p__k/R * divU_Up__n

                  Qgrad_visc_dot_grad_UR__UR__p = - (visco_R + visco_T / R) * divU_UR + visco_Z * (-UR_Z)
                  Qgrad_visc_dot_grad_UR__UR__n = - 1.d0 / R**2 * visco_p * (UR_p__n)
                  Qgrad_visc_dot_grad_UR__UZ__p = - (visco_R + visco_T / R) * divU_UZ + visco_Z * (UZ_R)
                  Qgrad_visc_dot_grad_UR__Up__p = - 1.d0 / R**2 * visco_p * (- Up - R * Up_R)
                  Qgrad_visc_dot_grad_UR__Up__n = - (visco_R + visco_T / R) * divU_Up__n

                  Qgrad_visc_dot_grad_UR__Te__p = - (visco_R_T + visco_T_T / R) * divU + visco_Z_T * (UZ0_R-UR0_Z) - 1.d0 / R**2 * visco_p_T__p * (UR0_p - Up0 - R * Up0_R)
                  Qgrad_visc_dot_grad_UR__Te__n = - 1.d0 / R**2 * visco_p_T__n * (UR0_p - Up0 - R * Up0_R)

                  Qgrad_visc_dot_grad_UZ__UR__p = - visco_Z * divU_UR - visco_R * (- UR_Z)
                  Qgrad_visc_dot_grad_UZ__UZ__p = - visco_Z * divU_UZ - visco_R * (UZ_R)
                  Qgrad_visc_dot_grad_UZ__UZ__n = + 1.d0 / R**2 * visco_p * (-UZ_p__n)
                  Qgrad_visc_dot_grad_UZ__Up__p = + 1.d0 / R**2 * visco_p * (R*Up_Z)
                  Qgrad_visc_dot_grad_UZ__Up__n = - visco_Z * divU_Up__n

                  Qgrad_visc_dot_grad_UZ__Te__p = - visco_Z_T * divU + 1.d0 / R**2 * visco_p_T__p * (R*Up0_Z-UZ0_p) - visco_R_T * (UZ0_R - UR0_Z)
                  Qgrad_visc_dot_grad_UZ__Te__n = + 1.d0 / R**2 * visco_p_T__n * (R*Up0_Z-UZ0_p)

                  Qgrad_visc_dot_grad_Up__UR__p = - 1.d0 / R * visco_p * divU_UR
                  Qgrad_visc_dot_grad_Up__UR__n = + 1.d0 / R * (visco_R + visco_T / R) * (UR_p__n)
                  Qgrad_visc_dot_grad_Up__UZ__p = - 1.d0 / R * visco_p * divU_UZ
                  Qgrad_visc_dot_grad_Up__UZ__n = - 1.d0 / R * visco_Z * (- UZ_p__n)
                  Qgrad_visc_dot_grad_Up__Up__p = + 1.d0 / R * (visco_R + visco_T / R) * (- Up - R*Up_R) - 1.d0 / R * visco_Z * (R * Up_Z)
                  Qgrad_visc_dot_grad_Up__Up__n = - 1.d0 / R * visco_p * divU_Up__n

                  Qgrad_visc_dot_grad_Up__Te__p = - 1.d0 / R * visco_p_T__p * divU + 1.d0 / R * (visco_R_T + visco_T_T / R) * (UR0_p - Up0 - R*Up0_R) - 1.d0 / R * visco_Z_T * (R * Up0_Z - UZ0_p)
                  Qgrad_visc_dot_grad_Up__Te__n = - 1.d0 / R * visco_p_T__n * divU

                  QdivU_BdotGradVisc__UR__p = divU_UR * (BR0 * visco_R + BZ0 * visco_Z + Bp0 * visco_p / R)
                  QdivU_BdotGradVisc__UZ__p = divU_UZ * (BR0 * visco_R + BZ0 * visco_Z + Bp0 * visco_p / R)
                  QdivU_BdotGradVisc__Up__n = divU_Up__n * (BR0 * visco_R + BZ0 * visco_Z + Bp0 * visco_p / R)

                  QdivU_BdotGradVisc__AR__p = divU * (BR0_AR * visco_R + Bp0_AR * visco_p / R)
                  QdivU_BdotGradVisc__AR__n = divU * (+ BZ0_AR__n * visco_Z)
                  QdivU_BdotGradVisc__AZ__p = divU * (+ BZ0_AZ * visco_Z + Bp0_AZ * visco_p / R)
                  QdivU_BdotGradVisc__AZ__n = divU * (BR0_AZ__n * visco_R)
                  QdivU_BdotGradVisc__A3__p = divU * (BR0_A3 * visco_R + BZ0_A3 * visco_Z + Bp0_A3 * visco_p / R)

                  QdivU_BdotGradVisc__Te__p = divU * (BR0 * visco_R_T + BZ0 * visco_Z_T + Bp0 * visco_p_T__p / R)
                  QdivU_BdotGradVisc__Te__n = divU * (+ Bp0 * visco_p_T__n / R)

                  curl_U0_dot_gradmu_x_B0__UR__p = &
                                    + (- UR_Z) * (visco_R * BZ0 - visco_Z * BR0)
                  curl_U0_dot_gradmu_x_B0__UR__n = &
                                    + 1.d0 / R**2 * (UR_p__n) * (visco_p * BR0 - R * visco_R * Bp0)
                  curl_U0_dot_gradmu_x_B0__UZ__p = &
                                    + (UZ_R) * (visco_R * BZ0 - visco_Z * BR0)
                  curl_U0_dot_gradmu_x_B0__UZ__n = (- UZ_p__n / R) * (visco_Z * Bp0 - 1.d0 / R * visco_p * BZ0) 
                  curl_U0_dot_gradmu_x_B0__Up__p = (Up_Z) * (visco_Z * Bp0 - 1.d0 / R * visco_p * BZ0) &
                                    + 1.d0 / R**2 * (- Up - R*Up_R) * (visco_p * BR0 - R * visco_R * Bp0) 

                  curl_U0_dot_gradmu_x_B0__AR__p = (Up0_Z - UZ0_p / R) * (visco_Z * Bp0_AR) &
                                    + 1.d0 / R**2 * (UR0_p - Up0 - R*Up0_R) * (visco_p * BR0_AR - R * visco_R * Bp0_AR) &
                                    + (UZ0_R - UR0_Z) * ( - visco_Z * BR0_AR)
                  curl_U0_dot_gradmu_x_B0__AR__n = (Up0_Z - UZ0_p / R) * (- 1.d0 / R * visco_p * BZ0_AR__n) &
                                    + (UZ0_R - UR0_Z) * (visco_R * BZ0_AR__n)
                  curl_U0_dot_gradmu_x_B0__AZ__p = (Up0_Z - UZ0_p / R) * (visco_Z * Bp0_AZ - 1.d0 / R * visco_p * BZ0_AZ) &
                                    + 1.d0 / R**2 * (UR0_p - Up0 - R*Up0_R) * (- R * visco_R * Bp0_AZ) &
                                    + (UZ0_R - UR0_Z) * (visco_R * BZ0_AZ)
                  curl_U0_dot_gradmu_x_B0__AZ__n = &
                                    + 1.d0 / R**2 * (UR0_p - Up0 - R*Up0_R) * (visco_p * BR0_AZ__n) &
                                    + (UZ0_R - UR0_Z) * (- visco_Z * BR0_AZ__n)
                  curl_U0_dot_gradmu_x_B0__A3__p = (Up0_Z - UZ0_p / R) * (visco_Z * Bp0_A3 - 1.d0 / R * visco_p * BZ0_A3) &
                                    + 1.d0 / R**2 * (UR0_p - Up0 - R*Up0_R) * (visco_p * BR0_A3 - R * visco_R * Bp0_A3) &
                                    + (UZ0_R - UR0_Z) * (visco_R * BZ0_A3 - visco_Z * BR0_A3)

                  curl_U0_dot_gradmu_x_B0__Te__p = (Up0_Z - UZ0_p / R) * (visco_Z_T * Bp0 - 1.d0 / R * visco_p_T__p * BZ0) &
                                    + 1.d0 / R**2 * (UR0_p - Up0 - R*Up0_R) * (visco_p_T__p * BR0 - R * visco_R_T * Bp0) &
                                    + (UZ0_R - UR0_Z) * (visco_R_T * BZ0 - visco_Z_T * BR0)
                  curl_U0_dot_gradmu_x_B0__Te__n = (Up0_Z - UZ0_p / R) * (- 1.d0 / R * visco_p_T__n * BZ0) &
                                    + 1.d0 / R**2 * (UR0_p - Up0 - R*Up0_R) * (visco_p_T__n * BR0)

                  curl_U0_dot_curl_B0__UR__p = &
                                + (- UR_Z) * (BZ0_R - BR0_Z)
                  curl_U0_dot_curl_B0__UR__n = &
                                + 1.d0 / R**2 * (UR_p__n) * (BR0_p - Bp0 - R * Bp0_R) 

                  curl_U0_dot_curl_B0__UZ__p =  &
                                + (UZ_R) * (BZ0_R - BR0_Z)
                  curl_U0_dot_curl_B0__UZ__n = (- UZ_p__n / R) * (Bp0_Z - 1.d0 / R * BZ0_p) 

                  curl_U0_dot_curl_B0__Up__p = (Up_Z) * (Bp0_Z - 1.d0 / R * BZ0_p) &
                                + 1.d0 / R**2 * (- Up - R*Up_R) * (BR0_p - Bp0 - R * Bp0_R) 

                  curl_U0_dot_curl_B0__AR__p = (Up0_Z - UZ0_p / R) * (grad_Bp_AR__Z) &
                                + 1.d0 / R**2 * (UR0_p - Up0 - R*Up0_R) * (grad_BR_AR__phi - Bp0_AR - R * grad_Bp_AR__R) &
                                + (UZ0_R - UR0_Z) * (- grad_BR_AR__Z)
                  curl_U0_dot_curl_B0__AR__n = &
                                + (UZ0_R - UR0_Z) * (grad_BZ_AR__R__n)
                  curl_U0_dot_curl_B0__AR__nn = (Up0_Z - UZ0_p / R) * (- 1.d0 / R * grad_BZ_AR__phi__nn) 

                  curl_U0_dot_curl_B0__AZ__p = (Up0_Z - UZ0_p / R) * (grad_Bp_AZ__Z - 1.d0 / R * grad_BZ_AZ__phi) &
                                + 1.d0 / R**2 * (UR0_p - Up0 - R*Up0_R) * (- Bp0_AZ - R * grad_Bp_AZ__R) &
                                + (UZ0_R - UR0_Z) * (grad_BZ_AZ__R)
                  curl_U0_dot_curl_B0__AZ__n = &
                                + (UZ0_R - UR0_Z) * (- grad_BR_AZ__Z__n)
                  curl_U0_dot_curl_B0__AZ__nn = &
                                + 1.d0 / R**2 * (UR0_p - Up0 - R*Up0_R) * (grad_BR_AZ__phi__nn)

                  curl_U0_dot_curl_B0__A3__p = (Up0_Z - UZ0_p / R) * (grad_Bp_A3__Z) &
                                + 1.d0 / R**2 * (UR0_p - Up0 - R*Up0_R) * (- Bp0_A3 - R * grad_Bp_A3__R) &
                                + (UZ0_R - UR0_Z) * (grad_BZ_A3__R - grad_BR_A3__Z)
                  curl_U0_dot_curl_B0__A3__n = (Up0_Z - UZ0_p / R) * (- 1.d0 / R * grad_BZ_A3__phi__n) &
                                + 1.d0 / R**2 * (UR0_p - Up0 - R*Up0_R) * (grad_BR_A3__phi__n) 
#else
                  ! --- full MHD viscosity
                  Qvisc_UR_UR__p  = - (v_R * UR_R + v_Z * UR_Z) - v * ( UR / R**2 ) 
                  Qvisc_UR_Up__n  = - v * (+ 2.d0 * Up_p__n / R**2 ) 
                  Qvisc_UR_UR__kn = - (v_p__k * UR_p__n / R**2)

                  Qvisc_UZ_UZ__p  = - (v_R * UZ_R + v_Z * UZ_Z)
                  Qvisc_UZ_UZ__kn = - (v_p__k * UZ_p__n / R**2)

                  Qvisc_Up_Up__p  = - (v_R * Up_R + v_Z * Up_Z) - v * ( Up / R**2 )
                  Qvisc_Up_UR__n  = - v * ( - 2.d0 * UR_p__n / R**2 )
                  Qvisc_Up_Up__kn = - (v_p__k * Up_p__n / R**2)

                  Qgrad_visc_dot_grad_UR__UR__p = visco_R * UR_R + visco_Z * UR_Z
                  Qgrad_visc_dot_grad_UR__UR__n = visco_p * UR_p__n / R**2

                  Qgrad_visc_dot_grad_UZ__UZ__p = visco_R * UZ_R + visco_Z * UZ_Z
                  Qgrad_visc_dot_grad_UZ__UZ__n = visco_p * UZ_p__n / R**2

                  Qgrad_visc_dot_grad_Up__Up__p = visco_R * Up_R + visco_Z * Up_Z
                  Qgrad_visc_dot_grad_Up__Up__n = visco_p * Up_p__n / R**2

                  Qgrad_visc_dot_grad_UR__Te__p = visco_R_T * UR0_R + visco_Z_T * UR0_Z + visco_p_T__p * UR0_p / R**2
                  Qgrad_visc_dot_grad_UR__Te__n = visco_p_T__n * UR0_p / R**2

                  Qgrad_visc_dot_grad_UZ__Te__p = visco_R_T * UZ0_R + visco_Z_T * UZ0_Z + visco_p_T__p * UZ0_p / R**2
                  Qgrad_visc_dot_grad_UZ__Te__n = visco_p_T__n * UZ0_p / R**2

                  Qgrad_visc_dot_grad_Up__Te__p = visco_R_T * Up0_R + visco_Z_T * Up0_Z + visco_p_T__p * Up0_p / R**2
                  Qgrad_visc_dot_grad_Up__Te__n = visco_p_T__n * Up0_p / R**2

                  Qvisc_grad_BR_dot_grad_UR__UR__p = BR0_R * UR_R + BR0_Z * UR_Z
                  Qvisc_grad_BR_dot_grad_UR__UR__n = BR0_p * UR_p__n / R**2
                  Qvisc_grad_BZ_dot_grad_UZ__UZ__p = BZ0_R * UZ_R + BZ0_Z * UZ_Z
                  Qvisc_grad_BZ_dot_grad_UZ__UZ__n = BZ0_p * UZ_p__n / R**2
                  Qvisc_grad_Bp_dot_grad_Up__Up__p = Bp0_R * Up_R + Bp0_Z * Up_Z
                  Qvisc_grad_Bp_dot_grad_Up__Up__n = Bp0_p * Up_p__n / R**2

                  Qvisc_grad_BR_dot_grad_UR__AR__p = grad_BR_AR__R * UR0_R + grad_BR_AR__Z * UR0_Z + grad_BR_AR__phi * UR0_p / R**2
                  Qvisc_grad_BR_dot_grad_UR__AZ__n = grad_BR_AZ__R__n * UR0_R + grad_BR_AZ__Z__n * UR0_Z
                  Qvisc_grad_BR_dot_grad_UR__AZ__nn = grad_BR_AZ__phi__nn * UR0_p / R**2
                  Qvisc_grad_BR_dot_grad_UR__A3__p = grad_BR_A3__R * UR0_R + grad_BR_A3__Z * UR0_Z 
                  Qvisc_grad_BR_dot_grad_UR__A3__n = grad_BR_A3__phi__n * UR0_p / R**2

                  Qvisc_grad_BZ_dot_grad_UZ__AR__p = grad_BZ_AR__R__n * UZ0_R + grad_BZ_AR__Z__n * UZ0_Z 
                  Qvisc_grad_BZ_dot_grad_UZ__AR__nn = grad_BZ_AR__phi__nn * UZ0_p / R**2
                  Qvisc_grad_BZ_dot_grad_UZ__AZ__p = grad_BZ_AZ__R * UZ0_R + grad_BZ_AZ__Z * UZ0_Z + grad_BZ_AZ__phi * UZ0_p / R**2
                  Qvisc_grad_BZ_dot_grad_UZ__A3__p = grad_BZ_A3__R * UZ0_R + grad_BZ_A3__Z * UZ0_Z 
                  Qvisc_grad_BZ_dot_grad_UZ__A3__n = grad_BZ_A3__phi__n * UZ0_p / R**2
                  
                  Qvisc_grad_Bp_dot_grad_Up__AR__p = grad_Bp_AR__R * Up0_R + grad_Bp_AR__Z * Up0_Z 
                  Qvisc_grad_Bp_dot_grad_Up__AR__n = grad_Bp_AR__phi__n * Up0_p / R**2
                  Qvisc_grad_Bp_dot_grad_Up__AZ__p = grad_Bp_AZ__R * Up0_R + grad_Bp_AZ__Z * Up0_Z 
                  Qvisc_grad_Bp_dot_grad_Up__AZ__n = grad_Bp_AZ__phi__n * Up0_p / R**2
                  Qvisc_grad_Bp_dot_grad_Up__A3__p = grad_Bp_A3__R * Up0_R + grad_Bp_A3__Z * Up0_Z + grad_Bp_A3__phi * Up0_p / R**2

#ifdef VISCOSITY_INSIDE_LAPLACIAN
                  Qvisc_gradmu_dot_grad_vecV_UR__UR__p = visco_R * UR_R + visco_Z * UR_Z
                  Qvisc_gradmu_dot_grad_vecV_UR__UR__n = visco_p * UR_p__n / R**2

                  Qvisc_gradmu_dot_grad_vecV_UR__Up__p = - visco_p * Up / R**2

                  Qvisc_gradmu_dot_grad_vecV_UR__Te__p = visco_R_T * UR0_R + visco_Z_T * UR0_Z + visco_p_T__p * UR0_p / R**2 - visco_p_T__p * Up0 / R**2
                  Qvisc_gradmu_dot_grad_vecV_UR__Te__n =                                         visco_p_T__n * UR0_p / R**2 - visco_p_T__n * Up0 / R**2


                  Qvisc_gradmu_dot_grad_vecV_UZ__UZ__p = visco_R * UZ_R + visco_Z * UZ_Z
                  Qvisc_gradmu_dot_grad_vecV_UZ__UZ__n = visco_p * UZ_p__n / R**2

                  Qvisc_gradmu_dot_grad_vecV_UZ__Te__p = visco_R_T * UZ0_R + visco_Z_T * UZ0_Z + visco_p_T__p * UZ0_p / R**2
                  Qvisc_gradmu_dot_grad_vecV_UZ__Te__n =                                         visco_p_T__n * UZ0_p / R**2

                  Qvisc_gradmu_dot_grad_vecV_Up__UR__p = visco_p * UR / R**2

                  Qvisc_gradmu_dot_grad_vecV_Up__Up__p = visco_R * Up_R + visco_Z * Up_Z
                  Qvisc_gradmu_dot_grad_vecV_Up__Up__n = visco_p * Up_p__n / R**2
                  Qvisc_gradmu_dot_grad_vecV_Up__Te__p = visco_R_T * Up0_R + visco_Z_T * Up0_Z + visco_p_T__p * Up0_p / R**2 + visco_p_T__p * UR0 / R**2
                  Qvisc_gradmu_dot_grad_vecV_Up__Te__n =                                         visco_p_T__n * Up0_p / R**2 + visco_p_T__n * UR0 / R**2

#endif
#ifdef VISCOSITY_COMPLETE_TERMS
                  Qgrad_visc_dot_gradU_T_UR__UR__p = visco_R * UR_R
                  Qgrad_visc_dot_gradU_T_UR__UZ__p = visco_Z * UZ_R
                  Qgrad_visc_dot_gradU_T_UR__Up__p = visco_p / R * Up_R
                  Qgrad_visc_dot_gradU_T_UR__Te__p = visco_R_T * UR0_R + visco_Z_T * UZ0_R + visco_p_T__p / R * Up0_R
                  Qgrad_visc_dot_gradU_T_UR__Te__n = visco_p_T__n / R * Up0_R

                  Qgrad_visc_dot_gradU_T_UZ__UR__p = visco_R * UR_Z
                  Qgrad_visc_dot_gradU_T_UZ__UZ__p = visco_Z * UZ_Z
                  Qgrad_visc_dot_gradU_T_UZ__Up__p = visco_p / R * Up_Z
                  Qgrad_visc_dot_gradU_T_UZ__Te__p = visco_R_T * UR0_Z + visco_Z_T * UZ0_Z + visco_p_T__p / R * Up0_Z
                  Qgrad_visc_dot_gradU_T_UZ__Te__n = visco_p_T__n / R * Up0_Z

                  Qgrad_visc_dot_gradU_T_Up__UR__p = visco_p / R**2 * UR
                  Qgrad_visc_dot_gradU_T_Up__UR__n = visco_R / R * UR_p__n
                  Qgrad_visc_dot_gradU_T_Up__UZ__n = visco_Z / R * UZ_p__n
                  Qgrad_visc_dot_gradU_T_Up__Up__p = visco_R / R * (- Up)
                  Qgrad_visc_dot_gradU_T_Up__Up__n = visco_p / R**2 * Up_p__n
                  Qgrad_visc_dot_gradU_T_Up__Te__p = visco_R_T / R * (UR0_p - Up0) + visco_Z_T / R * UZ0_p + visco_p_T__p / R**2 * (Up0_p + UR0)
                  Qgrad_visc_dot_gradU_T_Up__Te__n = visco_p_T__n / R**2 * (Up0_p + UR0)

                  Q_minus_2_over_3_divU_grad_visc_UR__UR__p = - 2.d0/3.d0 * visco_R * divU_UR
                  Q_minus_2_over_3_divU_grad_visc_UR__UZ__p = - 2.d0/3.d0 * visco_R * divU_UZ
                  Q_minus_2_over_3_divU_grad_visc_UR__Up__n = - 2.d0/3.d0 * visco_R * divU_Up__n
                  Q_minus_2_over_3_divU_grad_visc_UR__Te__p = - 2.d0/3.d0 * visco_R_T * divU

                  Q_minus_2_over_3_divU_grad_visc_UZ__UR__p = - 2.d0/3.d0 * visco_Z * divU_UR
                  Q_minus_2_over_3_divU_grad_visc_UZ__UZ__p = - 2.d0/3.d0 * visco_Z * divU_UZ
                  Q_minus_2_over_3_divU_grad_visc_UZ__Up__n = - 2.d0/3.d0 * visco_Z * divU_Up__n
                  Q_minus_2_over_3_divU_grad_visc_UZ__Te__p = - 2.d0/3.d0 * visco_Z_T * divU

                  Q_minus_2_over_3_divU_grad_visc_Up__UR__p = - 2.d0/3.d0 * visco_p * divU_UR / R
                  Q_minus_2_over_3_divU_grad_visc_Up__UZ__p = - 2.d0/3.d0 * visco_p * divU_UZ / R
                  Q_minus_2_over_3_divU_grad_visc_Up__Up__n = - 2.d0/3.d0 * visco_p * divU_Up__n / R
                  Q_minus_2_over_3_divU_grad_visc_Up__Te__p = - 2.d0/3.d0 * visco_p_T__p * divU / R
                  Q_minus_2_over_3_divU_grad_visc_Up__Te__n = - 2.d0/3.d0 * visco_p_T__n * divU / R

                  Q_1_over_3_visc_grad_divU_UR__UR__p = -1.d0/3.d0 * divU_UR * (visco * (v_R + v/R) + v * (visco_R + visco/R))
                  Q_1_over_3_visc_grad_divU_UR__UZ__p = -1.d0/3.d0 * divU_UZ * (visco * (v_R + v/R) + v * (visco_R + visco/R))
                  Q_1_over_3_visc_grad_divU_UR__Up__n = -1.d0/3.d0 * divU_Up__n * (visco * (v_R + v/R) + v * (visco_R + visco/R))
                  Q_1_over_3_visc_grad_divU_UR__Te__p = -1.d0/3.d0 * divU * (visco_T * (v_R + v/R) + v * (visco_R_T + visco_T/R))

                  Q_1_over_3_visc_grad_divU_UZ__UR__p = -1.d0/3.d0 * divU_UR * (visco * v_Z         + v * visco_Z)
                  Q_1_over_3_visc_grad_divU_UZ__UZ__p = -1.d0/3.d0 * divU_UZ * (visco * v_Z         + v * visco_Z)
                  Q_1_over_3_visc_grad_divU_UZ__Up__n = -1.d0/3.d0 * divU_Up__n * (visco * v_Z         + v * visco_Z)
                  Q_1_over_3_visc_grad_divU_UZ__Te__p = -1.d0/3.d0 * divU * (visco_T * v_Z         + v * visco_Z_T)

                  Q_1_over_3_visc_grad_divU_Up__UR__p = -1.d0/3.d0 * divU_UR * (                    + v * visco_p) / R
                  Q_1_over_3_visc_grad_divU_Up__UZ__p = -1.d0/3.d0 * divU_UZ * (                    + v * visco_p) / R
                  Q_1_over_3_visc_grad_divU_Up__Up__n = -1.d0/3.d0 * divU_Up__n * (                    + v * visco_p) / R
                  Q_1_over_3_visc_grad_divU_Up__Te__p = -1.d0/3.d0 * divU * (                    + v * visco_p_T__p) / R
                  Q_1_over_3_visc_grad_divU_Up__Te__n = -1.d0/3.d0 * divU * (                    + v * visco_p_T__n) / R
                  Q_1_over_3_visc_grad_divU_Up__UR__k = -1.d0/3.d0 * divU_UR * (visco * v_p__k                   ) / R
                  Q_1_over_3_visc_grad_divU_Up__UZ__k = -1.d0/3.d0 * divU_UZ * (visco * v_p__k                   ) / R
                  Q_1_over_3_visc_grad_divU_Up__Up__kn = -1.d0/3.d0 * divU_Up__n * (visco * v_p__k                   ) / R
                  Q_1_over_3_visc_grad_divU_Up__Te__k = -1.d0/3.d0 * divU * (visco_T * v_p__k                   ) / R

                  Q_1_over_3_B_dot_visc_grad_divU__AR__p = -1.d0/3.d0 * divU * (visco * (BR0_AR * v_R) + v * (BR0_AR * visco_R + Bp0_AR * visco_p / R))
                  Q_1_over_3_B_dot_visc_grad_divU__AR__n = -1.d0/3.d0 * divU * (visco * (BZ0_AR__n * v_Z) + v * (BZ0_AR__n * visco_Z))
                  Q_1_over_3_B_dot_visc_grad_divU__AZ__p = -1.d0/3.d0 * divU * (visco * (BZ0_AZ * v_Z) + v * (BZ0_AZ * visco_Z + Bp0_AZ * visco_p / R))
                  Q_1_over_3_B_dot_visc_grad_divU__AZ__n = -1.d0/3.d0 * divU * (visco * (BR0_AZ__n * v_R) + v * (BR0_AZ__n * visco_R))
                  Q_1_over_3_B_dot_visc_grad_divU__A3__p = -1.d0/3.d0 * divU * (visco * (BR0_A3 * v_R + BZ0_A3 * v_Z) + v * (BR0_A3 * visco_R + BZ0_A3 * visco_Z + Bp0_A3 * visco_p / R))
                  Q_1_over_3_B_dot_visc_grad_divU__UR__p = -1.d0/3.d0 * divU_UR * (visco * (BR0 * v_R + BZ0 * v_Z) + v * (BR0 * visco_R + BZ0 * visco_Z + Bp0 * visco_p / R))
                  Q_1_over_3_B_dot_visc_grad_divU__UZ__p = -1.d0/3.d0 * divU_UZ * (visco * (BR0 * v_R + BZ0 * v_Z) + v * (BR0 * visco_R + BZ0 * visco_Z + Bp0 * visco_p / R))
                  Q_1_over_3_B_dot_visc_grad_divU__Up__n = -1.d0/3.d0 * divU_Up__n * (visco * (BR0 * v_R + BZ0 * v_Z) + v * (BR0 * visco_R + BZ0 * visco_Z + Bp0 * visco_p / R))
                  Q_1_over_3_B_dot_visc_grad_divU__Te__p = -1.d0/3.d0 * divU * (visco_T * (BR0 * v_R + BZ0 * v_Z) + v * (BR0 * visco_R_T + BZ0 * visco_Z_T + Bp0 * visco_p_T__p / R))
                  Q_1_over_3_B_dot_visc_grad_divU__Te__n = -1.d0/3.d0 * divU * (v * (Bp0 * visco_p_T__n / R))

                  Q_1_over_3_B_dot_visc_grad_divU__AR__k = -1.d0/3.d0 * divU * (visco * (Bp0_AR * v_p__k / R))
                  Q_1_over_3_B_dot_visc_grad_divU__AZ__k = -1.d0/3.d0 * divU * (visco * (Bp0_AZ * v_p__k / R))
                  Q_1_over_3_B_dot_visc_grad_divU__A3__k = -1.d0/3.d0 * divU * (visco * (Bp0_A3 * v_p__k / R))
                  Q_1_over_3_B_dot_visc_grad_divU__UR__k = -1.d0/3.d0 * divU_UR * (visco * (Bp0 * v_p__k / R))
                  Q_1_over_3_B_dot_visc_grad_divU__UZ__k = -1.d0/3.d0 * divU_UZ * (visco * (Bp0 * v_p__k / R))
                  Q_1_over_3_B_dot_visc_grad_divU__Up__kn = -1.d0/3.d0 * divU_Up__n * (visco * (Bp0 * v_p__k / R))
                  Q_1_over_3_B_dot_visc_grad_divU__Te__k = -1.d0/3.d0 * divU * (visco_T * (Bp0 * v_p__k / R))
#endif

#endif

                  lap_bf = bf_R / R + bf_RR + bf_ZZ

                  amat      = 0.d0
                  Pjac      = 0.d0 ! time derivative part
                  Qjac_p    = 0.d0 ! rest of the LHS (poloidal part)
                  Qjac_k    = 0.d0 ! rest of the LHS (toroidal part with phi-derivatives of the test-function)
                  Qjac_n    = 0.d0 ! rest of the LHS (toroidal part with phi-derivatives of the basis-functions)
                  Qjac_kn   = 0.d0 ! rest of the LHS (toroidal part with phi-derivatives of the basis-functions and the test-functions)
                  Qjac_pnn  = 0.d0 ! rest of the LHS (toroidal part with second-order phi-derivatives of the basis-functions)
#ifdef ENABLE_VISCO_PAR_AND_PARPAR
                  Qjac_knn  = 0.d0 ! rest of the LHS (toroidal part with second-order phi-derivatives of the basis-functions and phi-derivatives of the test-function)
#endif

                  !###################################################################################################
                  !#  equation 1   (R component induction equation)                                                  #
                  !###################################################################################################
                  Pjac   (var_AR,var_AR) =   v * AR 

                  Qjac_p (var_AR,var_AR) = + v * (UZ0   * Bp0_AR ) &
                                           + tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2    * BR0 * BgradPe_AR__p &
                                           - tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2**2 * BR0 * BgradPe * BB2_AR__p &
                                           + eta_ARAZ * v * (eta_Z * Bp0_AR ) &
                                           - eta_ARAZ * eta_T * ( - v_Z * Bp0_AR ) &
                                           - eta_ARAZ_const   * ( - v_Z * Bp0_AR ) &
                                           - eta_ARAZ * eta_num * lap_Vstar * lap_bf
                  Qjac_n (var_AR,var_AR) = + v * (- Up0   * BZ0_AR__n     ) &
                                           + tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2    * BR0 * BgradPe_AR__n &
                                           - tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2**2 * BR0 * BgradPe * BB2_AR__n &
                                           + eta_ARAZ * v * (- eta_p * BZ0_AR__n / R )
                  Qjac_kn(var_AR,var_AR) = - eta_ARAZ * eta_T * ( + v_p__k * BZ0_AR__n / R)

                  Qjac_p (var_AR,var_AZ) = + v * (UZ0 * Bp0_AZ - Up0 * BZ0_AZ)          &
                                           + tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2    * BR0 * BgradPe_AZ__p &
                                           - tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2**2 * BR0 * BgradPe * BB2_AZ__p &
                                           + eta_ARAZ * v * (eta_Z * Bp0_AZ - eta_p * BZ0_AZ / R ) &
                                           - eta_ARAZ * eta_T * ( - v_Z * Bp0_AZ ) &
                                           - eta_ARAZ_const   * ( - v_Z * Bp0_AZ )
                  Qjac_n (var_AR,var_AZ) = + tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2    * BR0_AZ__n * BgradPe &
                                           + tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2    * BR0       * BgradPe_AZ__n &
                                           - tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2**2 * BR0       * BgradPe * BB2_AZ__n
                  Qjac_k (var_AR,var_AZ) = - eta_ARAZ * eta_T * ( + v_p__k * BZ0_AZ / R)

                  Qjac_p (var_AR,var_A3) = + v * (UZ0 * Bp0_A3 - Up0 * BZ0_A3)          &
                                           + tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2    * BR0_A3 * BgradPe &
                                           + tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2    * BR0    * BgradPe_A3 &
                                           - tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2**2 * BR0    * BgradPe * BB2_A3 &
                                           + eta_ARAZ * v * (eta_Z * Bp0_A3 - eta_p * BZ0_A3 / R ) &
                                           - eta_ARAZ * eta_T * ( - v_Z * Bp0_A3 ) &
                                           - eta_ARAZ_const   * ( - v_Z * Bp0_A3 )
                  Qjac_k (var_AR,var_A3) = - eta_ARAZ * eta_T * ( + v_p__k * BZ0_A3 / R)

                  Qjac_p (var_AR,var_UZ) = + v * (  UZ * Bp0)
                  Qjac_p (var_AR,var_Up) = + v * (- Up * BZ0)

                  Qjac_p (var_AR,var_rho)= + tauIC_ARAZ * v * tau_IC*F0/rho0_corr   /BB2 * BR0 * BgradPe_rho__p &
                                           - tauIC_ARAZ * v * tau_IC*F0/rho0_corr**2/BB2 * BR0 * BgradPe * rho
                  Qjac_n (var_AR,var_rho)= - tauIC_ARAZ * v * tau_IC*F0/rho0_corr   /BB2 * BR0 * BgradPe_rho__n

                  Qjac_p (var_AR,var_Te )= + tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2 * BR0 * BgradPe_Te__p &
                                           + eta_ARAZ * v * (eta_Z_T * Bp00 - eta_p_T__p * BZ0 / R ) &
                                           - eta_ARAZ * eta_T_T * ( - v_Z * Bp00 )                   &
                                           + eta_ARAZ * eta_T_T * v * current_source_JR(ms,mt)
                  Qjac_n (var_AR,var_Te )= + tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2 * BR0 * BgradPe_Te__n &
                                           + eta_ARAZ * v * (              - eta_p_T__n * BZ0 / R )
                  Qjac_k (var_AR,var_Te )= - eta_ARAZ * eta_T_T * ( + v_p__k * BZ0 / R)

                  !###################################################################################################
                  !#  equation 2   (Z component induction equation)                                                  #
                  !###################################################################################################
                  Pjac   (var_AZ,var_AZ) =   v * AZ 

                  Qjac_p (var_AZ,var_AR) = + v * (Up0 * BR0_AR - UR0 * Bp0_AR)         &
                                           + tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2    * BZ0 * BgradPe_AR__p &
                                           - tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2**2 * BZ0 * BgradPe * BB2_AR__p &
                                           + eta_ARAZ * v * (eta_p / R * BR0_AR - eta_R * Bp0_AR) &
                                           - eta_ARAZ * eta_T * ( + v_R * Bp0_AR ) &
                                           - eta_ARAZ_const   * ( + v_R * Bp0_AR )
                  Qjac_n (var_AZ,var_AR) = + tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2    * BZ0_AR__n * BgradPe &
                                           + tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2    * BZ0       * BgradPe_AR__n &
                                           - tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2**2 * BZ0       * BgradPe * BB2_AR__n
                  Qjac_k (var_AZ,var_AR) = - eta_ARAZ * eta_T * ( - v_p__k * BR0_AR / R)

                  Qjac_p (var_AZ,var_AZ) = + v * (- UR0   * Bp0_AZ) &
                                           + tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2    * BZ0 * BgradPe_AZ__p &
                                           - tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2**2 * BZ0 * BgradPe * BB2_AZ__p &
                                           + eta_ARAZ * v * (- eta_R * Bp0_AZ) &
                                           - eta_ARAZ * eta_T * ( + v_R * Bp0_AZ ) &
                                           - eta_ARAZ_const   * ( + v_R * Bp0_AZ ) &
                                           - eta_ARAZ * eta_num * lap_Vstar * lap_bf
                  Qjac_n (var_AZ,var_AZ) = + v * (Up0 * BR0_AZ__n)       &
                                           + tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2    * BZ0 * BgradPe_AZ__n &
                                           - tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2**2 * BZ0 * BgradPe * BB2_AZ__n &
                                           + eta_ARAZ * v * (eta_p / R * BR0_AZ__n)
                  Qjac_kn(var_AZ,var_AZ) = - eta_ARAZ * eta_T * ( - v_p__k * BR0_AZ__n / R)

                  Qjac_p (var_AZ,var_A3) = + v * (Up0 * BR0_A3 - UR0 * Bp0_A3)         &
                                           + tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2    * BZ0_A3 * BgradPe &
                                           + tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2    * BZ0    * BgradPe_A3 &
                                           - tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2**2 * BZ0    * BgradPe * BB2_A3 &
                                           + eta_ARAZ * v * (eta_p / R * BR0_A3 - eta_R * Bp0_A3) &
                                           - eta_ARAZ * eta_T * ( + v_R * Bp0_A3 ) &
                                           - eta_ARAZ_const   * ( + v_R * Bp0_A3 )
                  Qjac_k (var_AZ,var_A3) = - eta_ARAZ * eta_T * ( - v_p__k * BR0_A3 / R)

                  Qjac_p (var_AZ,var_UR) = + v * (- UR * Bp0)
                  Qjac_p (var_AZ,var_Up) = + v * (  Up * BR0)

                  Qjac_p (var_AZ,var_rho)= + tauIC_ARAZ * v * tau_IC*F0/rho0_corr   /BB2 * BZ0 * BgradPe_rho__p &
                                           - tauIC_ARAZ * v * tau_IC*F0/rho0_corr**2/BB2 * BZ0 * BgradPe * rho
                  Qjac_n (var_AZ,var_rho)= + tauIC_ARAZ * v * tau_IC*F0/rho0_corr   /BB2 * BZ0 * BgradPe_rho__n

                  Qjac_p (var_AZ,var_Te )= + tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2 * BZ0 * BgradPe_Te__p &
                                           + eta_ARAZ * v * (eta_p_T__p / R * BR0 - eta_R_T * Bp00) &
                                           - eta_ARAZ * eta_T_T * ( + v_R * Bp00 )                  &
                                           + eta_ARAZ * eta_T_T * v * current_source_JZ(ms,mt)
                  Qjac_n (var_AZ,var_Te )= + tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2 * BZ0 * BgradPe_Te__n &
                                           + eta_ARAZ * v * (eta_p_T__n / R * BR0 )
                  Qjac_k (var_AZ,var_Te )= - eta_ARAZ * eta_T_T * ( - v_p__k * BR0 / R)

                  !###################################################################################################
                  !#  equation 3   (Phi component induction equation)                                                #
                  !###################################################################################################
                  Pjac   (var_A3,var_A3) =   v * A3

                  Qjac_p (var_A3,var_AR) = + R * v * (- UZ0 * BR0_AR)                &
                                           + tauIC_AP * v * tau_IC*F0/rho0_corr/BB2    * R*Bp0_AR * BgradPe &
                                           + tauIC_AP * v * tau_IC*F0/rho0_corr/BB2    * R*Bp0    * BgradPe_AR__p &
                                           - tauIC_AP * v * tau_IC*F0/rho0_corr/BB2**2 * R*Bp0    * BgradPe * BB2_AR__p &
                                           - eta_T * v_Z * R                * BR0_AR &
                                           + R * v * (- eta_Z * BR0_AR)
                  Qjac_n (var_A3,var_AR) = + R * v * (UR0 * BZ0_AR__n)                  &
                                           + tauIC_AP * v * tau_IC*F0/rho0_corr/BB2    * R*Bp0    * BgradPe_AR__n &
                                           - tauIC_AP * v * tau_IC*F0/rho0_corr/BB2**2 * R*Bp0    * BgradPe * BB2_AR__n &
                                           + eta_T * ( 2.d0 * v + R * v_R ) * BZ0_AR__n &
                                           + R * v * (eta_R * BZ0_AR__n)

                  Qjac_p (var_A3,var_AZ) = + R * v * (UR0 * BZ0_AZ)                  &
                                           + tauIC_AP * v * tau_IC*F0/rho0_corr/BB2    * R*Bp0_AZ * BgradPe &
                                           + tauIC_AP * v * tau_IC*F0/rho0_corr/BB2    * R*Bp0    * BgradPe_AZ__p &
                                           - tauIC_AP * v * tau_IC*F0/rho0_corr/BB2**2 * R*Bp0    * BgradPe * BB2_AZ__p &
                                           + eta_T * ( 2.d0 * v + R * v_R ) * BZ0_AZ &
                                            + R * v * (eta_R * BZ0_AZ)
                  Qjac_n (var_A3,var_AZ) = + R * v * (- UZ0 * BR0_AZ__n)                &
                                           + tauIC_AP * v * tau_IC*F0/rho0_corr/BB2    * R*Bp0    * BgradPe_AZ__n &
                                           - tauIC_AP * v * tau_IC*F0/rho0_corr/BB2**2 * R*Bp0    * BgradPe * BB2_AZ__n &
                                           - eta_T * v_Z * R                * BR0_AZ__n &
                                           + R * v * (- eta_Z * BR0_AZ__n)

                  Qjac_p (var_A3,var_A3) = + R * v * (UR0 * BZ0_A3 - UZ0 * BR0_A3)   &
                                           + tauIC_AP * v * tau_IC*F0/rho0_corr/BB2    * R*Bp0    * BgradPe_A3 &
                                           - tauIC_AP * v * tau_IC*F0/rho0_corr/BB2**2 * R*Bp0    * BgradPe * BB2_A3 &
                                           - eta_T * v_Z * R                * BR0_A3 &
                                           + eta_T * ( 2.d0 * v + R * v_R ) * BZ0_A3 &
                                           + R * v * (eta_R * BZ0_A3 - eta_Z * BR0_A3) &
                                           - eta_num * lap_Vstar * lap_bf

                  Qjac_p (var_A3,var_UR) = + R * v * (  UR * BZ0)
                  Qjac_p (var_A3,var_UZ) = + R * v * (- UZ * BR0)

                  Qjac_p (var_A3,var_rho)= + tauIC_AP * v * tau_IC*F0/rho0_corr   /BB2 * R*Bp0 * BgradPe_rho__p &
                                           - tauIC_AP * v * tau_IC*F0/rho0_corr**2/BB2 * R*Bp0 * BgradPe * rho
                  Qjac_n (var_A3,var_rho)= + tauIC_AP * v * tau_IC*F0/rho0_corr/BB2 * R*Bp0 * BgradPe_rho__n

                  Qjac_p (var_A3,var_Te )= + tauIC_AP * v * tau_IC*F0/rho0_corr/BB2 * R*Bp0 * BgradPe_Te__p &
                                           - eta_T_T * v_Z * R                * BR0 &
                                           + eta_T_T * ( 2.d0 * v + R * v_R ) * BZ0 &
                                           + eta_T_T * v * current_source_Jp(ms,mt) &
                                           + R * v * (eta_R_T * BZ0 - eta_Z_T * BR0)
                  Qjac_n (var_A3,var_Te )= + tauIC_AP * v * tau_IC*F0/rho0_corr/BB2 * R*Bp0 * BgradPe_Te__n

                  !###################################################################################################
                  !#  equation 4   (R component momentum equation)                                                   #
                  !###################################################################################################
                  Pjac(var_UR,var_UR )    =   v * rho0_corr * UR

                  Qjac_p (var_UR,var_AR ) =  Qconv_UR_AR__p + JxB_UR_AR__p  - v * PneoR_AR__p
                  Qjac_n (var_UR,var_AR ) =  Qconv_UR_AR__n + JxB_UR_AR__n  - v * PneoR_AR__n
                  Qjac_k (var_UR,var_AR ) =                 + JxB_UR_AR__k
                  Qjac_kn(var_UR,var_AR ) =                 + JxB_UR_AR__kn
#ifdef ENABLE_VISCO_PAR_AND_PARPAR
                  Qjac_p(var_UR, var_AR) = Qjac_p(var_UR, var_AR) + &
                    -visco_par_prof * Qjac_UR_visco_par_AR_p &
                    -visco_pp_diff_prof * Qjac_UR_visco_parpar_AR_p

                  Qjac_n(var_UR, var_AR) = Qjac_n(var_UR, var_AR) + &
                    -visco_par_prof * Qjac_UR_visco_par_AR_n &
                    -visco_pp_diff_prof * Qjac_UR_visco_parpar_AR_n

                  Qjac_k(var_UR, var_AR) = Qjac_k(var_UR, var_AR) + &
                    -visco_par_prof * Qjac_UR_visco_par_AR_k &
                    -visco_pp_diff_prof * Qjac_UR_visco_parpar_AR_k

                  Qjac_kn(var_UR, var_AR) = Qjac_kn(var_UR, var_AR) + &
                    -visco_par_prof * Qjac_UR_visco_par_AR_kn &
                    -visco_pp_diff_prof * Qjac_UR_visco_parpar_AR_kn

                  Qjac_pnn(var_UR, var_AR) = Qjac_pnn(var_UR, var_AR) + &
                    -visco_par_prof * Qjac_UR_visco_par_AR_pnn &
                    -visco_pp_diff_prof * Qjac_UR_visco_parpar_AR_pnn

                  Qjac_knn(var_UR, var_AR) = Qjac_knn(var_UR, var_AR) + &
                    -visco_par_prof * Qjac_UR_visco_par_AR_knn &
                    -visco_pp_diff_prof * Qjac_UR_visco_parpar_AR_knn
#endif

                  Qjac_p (var_UR,var_AZ ) =  Qconv_UR_AZ__p + JxB_UR_AZ__p  - v * PneoR_AZ__p
                  Qjac_n (var_UR,var_AZ ) =  Qconv_UR_AZ__n + JxB_UR_AZ__n  - v * PneoR_AZ__n
                  Qjac_k (var_UR,var_AZ ) =                 + JxB_UR_AZ__k
                  Qjac_kn(var_UR,var_AZ ) =                 + JxB_UR_AZ__kn
#ifdef ENABLE_VISCO_PAR_AND_PARPAR
                  Qjac_p(var_UR, var_AZ) = Qjac_p(var_UR, var_AZ) + &
                    -visco_par_prof * Qjac_UR_visco_par_AZ_p &
                    -visco_pp_diff_prof * Qjac_UR_visco_parpar_AZ_p

                  Qjac_n(var_UR, var_AZ) = Qjac_n(var_UR, var_AZ) + &
                    -visco_par_prof * Qjac_UR_visco_par_AZ_n &
                    -visco_pp_diff_prof * Qjac_UR_visco_parpar_AZ_n

                  Qjac_k(var_UR, var_AZ) = Qjac_k(var_UR, var_AZ) + &
                    -visco_par_prof * Qjac_UR_visco_par_AZ_k &
                    -visco_pp_diff_prof * Qjac_UR_visco_parpar_AZ_k

                  Qjac_kn(var_UR, var_AZ) = Qjac_kn(var_UR, var_AZ) + &
                    -visco_par_prof * Qjac_UR_visco_par_AZ_kn &
                    -visco_pp_diff_prof * Qjac_UR_visco_parpar_AZ_kn

                  Qjac_pnn(var_UR, var_AZ) = Qjac_pnn(var_UR, var_AZ) + &
                    -visco_par_prof * Qjac_UR_visco_par_AZ_pnn &
                    -visco_pp_diff_prof * Qjac_UR_visco_parpar_AZ_pnn

                  Qjac_knn(var_UR, var_AZ) = Qjac_knn(var_UR, var_AZ) + &
                    -visco_par_prof * Qjac_UR_visco_par_AZ_knn &
                    -visco_pp_diff_prof * Qjac_UR_visco_parpar_AZ_knn
#endif

                  Qjac_p (var_UR,var_A3 ) =  Qconv_UR_A3__p + JxB_UR_A3__p  - v * PneoR_A3__p
                  Qjac_n (var_UR,var_A3 ) =  Qconv_UR_A3__n + JxB_UR_A3__n  - v * PneoR_A3__n
                  Qjac_k (var_UR,var_A3 ) =                 + JxB_UR_A3__k
                  Qjac_kn(var_UR,var_A3 ) =                 + JxB_UR_A3__kn
#ifdef ENABLE_VISCO_PAR_AND_PARPAR
                  Qjac_p(var_UR, var_A3) = Qjac_p(var_UR, var_A3) + &
                    -visco_par_prof * Qjac_UR_visco_par_A3_p &
                    -visco_pp_diff_prof * Qjac_UR_visco_parpar_A3_p

                  Qjac_n(var_UR, var_A3) = Qjac_n(var_UR, var_A3) + &
                    -visco_par_prof * Qjac_UR_visco_par_A3_n &
                    -visco_pp_diff_prof * Qjac_UR_visco_parpar_A3_n

                  Qjac_k(var_UR, var_A3) = Qjac_k(var_UR, var_A3) + &
                    -visco_par_prof * Qjac_UR_visco_par_A3_k &
                    -visco_pp_diff_prof * Qjac_UR_visco_parpar_A3_k

                  Qjac_kn(var_UR, var_A3) = Qjac_kn(var_UR, var_A3) + &
                    -visco_par_prof * Qjac_UR_visco_par_A3_kn &
                    -visco_pp_diff_prof * Qjac_UR_visco_parpar_A3_kn

                  Qjac_pnn(var_UR, var_A3) = Qjac_pnn(var_UR, var_A3) + &
                    -visco_par_prof * Qjac_UR_visco_par_A3_pnn &
                    -visco_pp_diff_prof * Qjac_UR_visco_parpar_A3_pnn

                  Qjac_knn(var_UR, var_A3) = Qjac_knn(var_UR, var_A3) + &
                    -visco_par_prof * Qjac_UR_visco_par_A3_knn &
                    -visco_pp_diff_prof * Qjac_UR_visco_parpar_A3_knn
#endif

                  Qjac_p (var_UR,var_UR ) =  Qconv_UR_UR__p + visco_T * Qvisc_UR_UR__p  - v * PneoR_UR - v * particle_source(ms,mt) * UR &
                                            - visco_num * lap_Vstar * lap_bf
                  Qjac_n (var_UR,var_UR ) =  Qconv_UR_UR__n + visco_T * Qvisc_UR_UR__n
                  Qjac_k (var_UR,var_UR ) =                 + visco_T * Qvisc_UR_UR__k
                  Qjac_kn(var_UR,var_UR ) =                 + visco_T * Qvisc_UR_UR__kn

#ifdef VISCOSITY_TENSOR_FORM
                  Qjac_p (var_UR,var_UR ) = Qjac_p (var_UR,var_UR ) + v * Qgrad_visc_dot_grad_UR__UR__p
                  Qjac_n (var_UR,var_UR ) = Qjac_n (var_UR,var_UR ) + v * Qgrad_visc_dot_grad_UR__UR__n
#else
                  ! --- the missed term in full MHD viscosity
                  Qjac_p (var_UR,var_UR ) = Qjac_p (var_UR,var_UR ) - v * Qgrad_visc_dot_grad_UR__UR__p
                  Qjac_n (var_UR,var_UR ) = Qjac_n (var_UR,var_UR ) - v * Qgrad_visc_dot_grad_UR__UR__n
#ifdef VISCOSITY_INSIDE_LAPLACIAN
                  Qjac_p (var_UR,var_UR ) = Qjac_p (var_UR,var_UR ) + v * Qvisc_gradmu_dot_grad_vecV_UR__UR__p
                  Qjac_n (var_UR,var_UR ) = Qjac_n (var_UR,var_UR ) + v * Qvisc_gradmu_dot_grad_vecV_UR__UR__n
#endif
#ifdef VISCOSITY_COMPLETE_TERMS
                  Qjac_p (var_UR,var_UR ) = Qjac_p (var_UR,var_UR ) + v * Qgrad_visc_dot_gradU_T_UR__UR__p          &
                                          + v * Q_minus_2_over_3_divU_grad_visc_UR__UR__p &
                                          +     Q_1_over_3_visc_grad_divU_UR__UR__p
#endif
#endif

#ifdef ENABLE_VISCO_PAR_AND_PARPAR
                  Qjac_p(var_UR, var_UR) = Qjac_p(var_UR, var_UR) + &
                    -visco_par_prof * Qjac_UR_visco_par_UR_p &
                    -visco_pp_diff_prof * Qjac_UR_visco_parpar_UR_p

                  Qjac_n(var_UR, var_UR) = Qjac_n(var_UR, var_UR) + &
                    -visco_par_prof * Qjac_UR_visco_par_UR_n &
                    -visco_pp_diff_prof * Qjac_UR_visco_parpar_UR_n

                  Qjac_k(var_UR, var_UR) = Qjac_k(var_UR, var_UR) + &
                    -visco_par_prof * Qjac_UR_visco_par_UR_k &
                    -visco_pp_diff_prof * Qjac_UR_visco_parpar_UR_k

                  Qjac_kn(var_UR, var_UR) = Qjac_kn(var_UR, var_UR) + &
                    -visco_par_prof * Qjac_UR_visco_par_UR_kn &
                    -visco_pp_diff_prof * Qjac_UR_visco_parpar_UR_kn

                  Qjac_pnn(var_UR, var_UR) = Qjac_pnn(var_UR, var_UR) + &
                    -visco_par_prof * Qjac_UR_visco_par_UR_pnn &
                    -visco_pp_diff_prof * Qjac_UR_visco_parpar_UR_pnn

                  Qjac_knn(var_UR, var_UR) = Qjac_knn(var_UR, var_UR) + &
                    -visco_par_prof * Qjac_UR_visco_par_UR_knn &
                    -visco_pp_diff_prof * Qjac_UR_visco_parpar_UR_knn
#endif

                  Qjac_p (var_UR,var_UZ ) =  Qconv_UR_UZ__p + visco_T * Qvisc_UR_UZ__p  - v * PneoR_UZ
                  Qjac_n (var_UR,var_UZ ) =  Qconv_UR_UZ__n + visco_T * Qvisc_UR_UZ__n
                  Qjac_k (var_UR,var_UZ ) =                 + visco_T * Qvisc_UR_UZ__k
                  Qjac_kn(var_UR,var_UZ ) =                 + visco_T * Qvisc_UR_UZ__kn
#ifdef VISCOSITY_TENSOR_FORM
                  Qjac_p (var_UR,var_UZ ) = Qjac_p (var_UR,var_UZ ) + v * Qgrad_visc_dot_grad_UR__UZ__p
#else
#ifdef VISCOSITY_COMPLETE_TERMS
                  Qjac_p (var_UR,var_UZ ) = Qjac_p (var_UR,var_UZ ) + v * Qgrad_visc_dot_gradU_T_UR__UZ__p          &
                                          + v * Q_minus_2_over_3_divU_grad_visc_UR__UZ__p &
                                          +     Q_1_over_3_visc_grad_divU_UR__UZ__p
#endif
#endif

#ifdef ENABLE_VISCO_PAR_AND_PARPAR
                  Qjac_p(var_UR, var_UZ) = Qjac_p(var_UR, var_UZ) + &
                    -visco_par_prof * Qjac_UR_visco_par_UZ_p &
                    -visco_pp_diff_prof * Qjac_UR_visco_parpar_UZ_p

                  Qjac_n(var_UR, var_UZ) = Qjac_n(var_UR, var_UZ) + &
                    -visco_par_prof * Qjac_UR_visco_par_UZ_n &
                    -visco_pp_diff_prof * Qjac_UR_visco_parpar_UZ_n

                  Qjac_k(var_UR, var_UZ) = Qjac_k(var_UR, var_UZ) + &
                    -visco_par_prof * Qjac_UR_visco_par_UZ_k &
                    -visco_pp_diff_prof * Qjac_UR_visco_parpar_UZ_k

                  Qjac_kn(var_UR, var_UZ) = Qjac_kn(var_UR, var_UZ) + &
                    -visco_par_prof * Qjac_UR_visco_par_UZ_kn &
                    -visco_pp_diff_prof * Qjac_UR_visco_parpar_UZ_kn

                  Qjac_pnn(var_UR, var_UZ) = Qjac_pnn(var_UR, var_UZ) + &
                    -visco_par_prof * Qjac_UR_visco_par_UZ_pnn &
                    -visco_pp_diff_prof * Qjac_UR_visco_parpar_UZ_pnn

                  Qjac_knn(var_UR, var_UZ) = Qjac_knn(var_UR, var_UZ) + &
                    -visco_par_prof * Qjac_UR_visco_par_UZ_knn &
                    -visco_pp_diff_prof * Qjac_UR_visco_parpar_UZ_knn
#endif

                  Qjac_p (var_UR,var_Up ) =  Qconv_UR_Up__p + visco_T * Qvisc_UR_Up__p
                  Qjac_n (var_UR,var_Up ) =  Qconv_UR_Up__n + visco_T * Qvisc_UR_Up__n
                  Qjac_k (var_UR,var_Up ) =                 + visco_T * Qvisc_UR_Up__k
                  Qjac_kn(var_UR,var_Up ) =                 + visco_T * Qvisc_UR_Up__kn
#ifdef VISCOSITY_TENSOR_FORM
                  Qjac_p (var_UR,var_Up ) = Qjac_p (var_UR,var_Up ) + v * Qgrad_visc_dot_grad_UR__Up__p
                  Qjac_n (var_UR,var_Up ) = Qjac_n (var_UR,var_Up ) + v * Qgrad_visc_dot_grad_UR__Up__n
#else
#ifdef VISCOSITY_INSIDE_LAPLACIAN
                  Qjac_p (var_UR,var_Up ) = Qjac_p (var_UR,var_Up ) + v * Qvisc_gradmu_dot_grad_vecV_UR__Up__p
#endif
#ifdef VISCOSITY_COMPLETE_TERMS
                  Qjac_p (var_UR,var_Up ) = Qjac_p (var_UR,var_Up ) + v * Qgrad_visc_dot_gradU_T_UR__Up__p
                  Qjac_n (var_UR,var_Up ) = Qjac_n (var_UR,var_Up )                       &
                                          + v * Q_minus_2_over_3_divU_grad_visc_UR__Up__n &
                                          +     Q_1_over_3_visc_grad_divU_UR__Up__n
#endif
#endif
#ifdef ENABLE_VISCO_PAR_AND_PARPAR
                  Qjac_p(var_UR, var_Up) = Qjac_p(var_UR, var_Up) + &
                    -visco_par_prof * Qjac_UR_visco_par_Up_p &
                    -visco_pp_diff_prof * Qjac_UR_visco_parpar_Up_p

                  Qjac_n(var_UR, var_Up) = Qjac_n(var_UR, var_Up) + &
                    -visco_par_prof * Qjac_UR_visco_par_Up_n &
                    -visco_pp_diff_prof * Qjac_UR_visco_parpar_Up_n

                  Qjac_k(var_UR, var_Up) = Qjac_k(var_UR, var_Up) + &
                    -visco_par_prof * Qjac_UR_visco_par_Up_k &
                    -visco_pp_diff_prof * Qjac_UR_visco_parpar_Up_k

                  Qjac_kn(var_UR, var_Up) = Qjac_kn(var_UR, var_Up) + &
                    -visco_par_prof * Qjac_UR_visco_par_Up_kn &
                    -visco_pp_diff_prof * Qjac_UR_visco_parpar_Up_kn

                  Qjac_pnn(var_UR, var_Up) = Qjac_pnn(var_UR, var_Up) + &
                    -visco_par_prof * Qjac_UR_visco_par_Up_pnn &
                    -visco_pp_diff_prof * Qjac_UR_visco_parpar_Up_pnn

                  Qjac_knn(var_UR, var_Up) = Qjac_knn(var_UR, var_Up) + &
                    -visco_par_prof * Qjac_UR_visco_par_Up_knn &
                    -visco_pp_diff_prof * Qjac_UR_visco_parpar_Up_knn
#endif

                  Qjac_p (var_UR,var_rho) =  Qconv_UR_rho__p + ( v_R + v / R ) * (rho*(Ti0+Te0)) - v * PneoR_rho__p
                  Qjac_n (var_UR,var_rho) =  Qconv_UR_rho__n                                     - v * PneoR_rho__n

                  Qjac_p (var_UR,var_Ti ) =  Qconv_UR_Ti__p  + ( v_R + v / R ) * (rho0*Ti      ) - v * PneoR_Ti__p
                  Qjac_n (var_UR,var_Ti ) =  Qconv_UR_Ti__n                                      - v * PneoR_Ti__n

                  Qjac_p (var_UR,var_Te ) =                  + ( v_R + v / R ) * (rho0*Te      ) + visco_T_T * Qvisc_UR__p
                  Qjac_k (var_UR,var_Te ) =                                                      + visco_T_T * Qvisc_UR__k

#ifdef VISCOSITY_TENSOR_FORM
                  Qjac_p (var_UR,var_Te ) = Qjac_p (var_UR,var_Te ) + v * Qgrad_visc_dot_grad_UR__Te__p
                  Qjac_n (var_UR,var_Te ) = Qjac_n (var_UR,var_Te ) + v * Qgrad_visc_dot_grad_UR__Te__n
#else
                  ! --- the missed term in full MHD viscosity
                  Qjac_p (var_UR,var_Te ) = Qjac_p (var_UR,var_Te ) - v * Qgrad_visc_dot_grad_UR__Te__p
                  Qjac_n (var_UR,var_Te ) = Qjac_n (var_UR,var_Te ) - v * Qgrad_visc_dot_grad_UR__Te__n
#ifdef VISCOSITY_INSIDE_LAPLACIAN
                  Qjac_p (var_UR,var_Te ) = Qjac_p (var_UR,var_Te ) + v * Qvisc_gradmu_dot_grad_vecV_UR__Te__p
                  Qjac_n (var_UR,var_Te ) = Qjac_n (var_UR,var_Te ) + v * Qvisc_gradmu_dot_grad_vecV_UR__Te__n
#endif
#ifdef VISCOSITY_COMPLETE_TERMS
                  Qjac_p (var_UR,var_Te ) = Qjac_p (var_UR,var_Te ) + v * Qgrad_visc_dot_gradU_T_UR__Te__p          &
                                + v * Q_minus_2_over_3_divU_grad_visc_UR__Te__p &
                                +     Q_1_over_3_visc_grad_divU_UR__Te__p
                  Qjac_n (var_UR,var_Te ) = Qjac_n (var_UR,var_Te ) + v * Qgrad_visc_dot_gradU_T_UR__Te__n
#endif
#endif

                  !###################################################################################################
                  !#  equation 5   (Z component momentum equation)                                                   #
                  !###################################################################################################
                  Pjac(var_UZ,var_UZ)      =   v * rho0_corr * UZ

                  Qjac_p (var_UZ,var_AR )  =  Qconv_UZ_AR__p + JxB_UZ_AR__p  - v * PneoZ_AR__p
                  Qjac_n (var_UZ,var_AR )  =  Qconv_UZ_AR__n + JxB_UZ_AR__n  - v * PneoZ_AR__n
                  Qjac_k (var_UZ,var_AR )  =                 + JxB_UZ_AR__k
                  Qjac_kn(var_UZ,var_AR )  =                 + JxB_UZ_AR__kn
#ifdef ENABLE_VISCO_PAR_AND_PARPAR
                  Qjac_p(var_UZ, var_AR) = Qjac_p(var_UZ, var_AR) + &
                    -visco_par_prof * Qjac_UZ_visco_par_AR_p &
                    -visco_pp_diff_prof * Qjac_UZ_visco_parpar_AR_p

                  Qjac_n(var_UZ, var_AR) = Qjac_n(var_UZ, var_AR) + &
                    -visco_par_prof * Qjac_UZ_visco_par_AR_n &
                    -visco_pp_diff_prof * Qjac_UZ_visco_parpar_AR_n

                  Qjac_k(var_UZ, var_AR) = Qjac_k(var_UZ, var_AR) + &
                    -visco_par_prof * Qjac_UZ_visco_par_AR_k &
                    -visco_pp_diff_prof * Qjac_UZ_visco_parpar_AR_k

                  Qjac_kn(var_UZ, var_AR) = Qjac_kn(var_UZ, var_AR) + &
                    -visco_par_prof * Qjac_UZ_visco_par_AR_kn &
                    -visco_pp_diff_prof * Qjac_UZ_visco_parpar_AR_kn

                  Qjac_pnn(var_UZ, var_AR) = Qjac_pnn(var_UZ, var_AR) + &
                    -visco_par_prof * Qjac_UZ_visco_par_AR_pnn &
                    -visco_pp_diff_prof * Qjac_UZ_visco_parpar_AR_pnn

                  Qjac_knn(var_UZ, var_AR) = Qjac_knn(var_UZ, var_AR) + &
                    -visco_par_prof * Qjac_UZ_visco_par_AR_knn &
                    -visco_pp_diff_prof * Qjac_UZ_visco_parpar_AR_knn
#endif

                  Qjac_p (var_UZ,var_AZ )  =  Qconv_UZ_AZ__p + JxB_UZ_AZ__p  - v * PneoZ_AZ__p
                  Qjac_n (var_UZ,var_AZ )  =  Qconv_UZ_AZ__n + JxB_UZ_AZ__n  - v * PneoZ_AZ__n
                  Qjac_k (var_UZ,var_AZ )  =                 + JxB_UZ_AZ__k
                  Qjac_kn(var_UZ,var_AZ )  =                 + JxB_UZ_AZ__kn
#ifdef ENABLE_VISCO_PAR_AND_PARPAR
                  Qjac_p(var_UZ, var_AZ) = Qjac_p(var_UZ, var_AZ) + &
                    -visco_par_prof * Qjac_UZ_visco_par_AZ_p &
                    -visco_pp_diff_prof * Qjac_UZ_visco_parpar_AZ_p

                  Qjac_n(var_UZ, var_AZ) = Qjac_n(var_UZ, var_AZ) + &
                    -visco_par_prof * Qjac_UZ_visco_par_AZ_n &
                    -visco_pp_diff_prof * Qjac_UZ_visco_parpar_AZ_n

                  Qjac_k(var_UZ, var_AZ) = Qjac_k(var_UZ, var_AZ) + &
                    -visco_par_prof * Qjac_UZ_visco_par_AZ_k &
                    -visco_pp_diff_prof * Qjac_UZ_visco_parpar_AZ_k

                  Qjac_kn(var_UZ, var_AZ) = Qjac_kn(var_UZ, var_AZ) + &
                    -visco_par_prof * Qjac_UZ_visco_par_AZ_kn &
                    -visco_pp_diff_prof * Qjac_UZ_visco_parpar_AZ_kn

                  Qjac_pnn(var_UZ, var_AZ) = Qjac_pnn(var_UZ, var_AZ) + &
                    -visco_par_prof * Qjac_UZ_visco_par_AZ_pnn &
                    -visco_pp_diff_prof * Qjac_UZ_visco_parpar_AZ_pnn

                  Qjac_knn(var_UZ, var_AZ) = Qjac_knn(var_UZ, var_AZ) + &
                    -visco_par_prof * Qjac_UZ_visco_par_AZ_knn &
                    -visco_pp_diff_prof * Qjac_UZ_visco_parpar_AZ_knn
#endif

                  Qjac_p (var_UZ,var_A3 )  =  Qconv_UZ_A3__p + JxB_UZ_A3__p  - v * PneoZ_A3__p
                  Qjac_n (var_UZ,var_A3 )  =  Qconv_UZ_A3__n + JxB_UZ_A3__n  - v * PneoZ_A3__n
                  Qjac_k (var_UZ,var_A3 )  =                 + JxB_UZ_A3__k
                  Qjac_kn(var_UZ,var_A3 )  =                 + JxB_UZ_A3__kn
#ifdef ENABLE_VISCO_PAR_AND_PARPAR
                  Qjac_p(var_UZ, var_A3) = Qjac_p(var_UZ, var_A3) + &
                    -visco_par_prof * Qjac_UZ_visco_par_A3_p &
                    -visco_pp_diff_prof * Qjac_UZ_visco_parpar_A3_p

                  Qjac_n(var_UZ, var_A3) = Qjac_n(var_UZ, var_A3) + &
                    -visco_par_prof * Qjac_UZ_visco_par_A3_n &
                    -visco_pp_diff_prof * Qjac_UZ_visco_parpar_A3_n

                  Qjac_k(var_UZ, var_A3) = Qjac_k(var_UZ, var_A3) + &
                    -visco_par_prof * Qjac_UZ_visco_par_A3_k &
                    -visco_pp_diff_prof * Qjac_UZ_visco_parpar_A3_k

                  Qjac_kn(var_UZ, var_A3) = Qjac_kn(var_UZ, var_A3) + &
                    -visco_par_prof * Qjac_UZ_visco_par_A3_kn &
                    -visco_pp_diff_prof * Qjac_UZ_visco_parpar_A3_kn

                  Qjac_pnn(var_UZ, var_A3) = Qjac_pnn(var_UZ, var_A3) + &
                    -visco_par_prof * Qjac_UZ_visco_par_A3_pnn &
                    -visco_pp_diff_prof * Qjac_UZ_visco_parpar_A3_pnn

                  Qjac_knn(var_UZ, var_A3) = Qjac_knn(var_UZ, var_A3) + &
                    -visco_par_prof * Qjac_UZ_visco_par_A3_knn &
                    -visco_pp_diff_prof * Qjac_UZ_visco_parpar_A3_knn
#endif    

                  Qjac_p (var_UZ,var_UR )  =  Qconv_UZ_UR__p + visco_T * Qvisc_UZ_UR__p - v * PneoZ_UR
                  Qjac_n (var_UZ,var_UR )  =  Qconv_UZ_UR__n + visco_T * Qvisc_UZ_UR__n
                  Qjac_k (var_UZ,var_UR )  =                 + visco_T * Qvisc_UZ_UR__k
                  Qjac_kn(var_UZ,var_UR )  =                 + visco_T * Qvisc_UZ_UR__kn
#ifdef VISCOSITY_TENSOR_FORM
                  Qjac_p (var_UZ,var_UR ) = Qjac_p (var_UZ,var_UR ) + v * Qgrad_visc_dot_grad_UZ__UR__p
#else
#ifdef VISCOSITY_COMPLETE_TERMS
                  Qjac_p (var_UZ,var_UR ) = Qjac_p (var_UZ,var_UR ) + v * Qgrad_visc_dot_gradU_T_UZ__UR__p          &
                                + v * Q_minus_2_over_3_divU_grad_visc_UZ__UR__p &
                                +     Q_1_over_3_visc_grad_divU_UZ__UR__p
#endif
#endif
#ifdef ENABLE_VISCO_PAR_AND_PARPAR
                  Qjac_p(var_UZ, var_UR) = Qjac_p(var_UZ, var_UR) + &
                    -visco_par_prof * Qjac_UZ_visco_par_UR_p &
                    -visco_pp_diff_prof * Qjac_UZ_visco_parpar_UR_p

                  Qjac_n(var_UZ, var_UR) = Qjac_n(var_UZ, var_UR) + &
                    -visco_par_prof * Qjac_UZ_visco_par_UR_n &
                    -visco_pp_diff_prof * Qjac_UZ_visco_parpar_UR_n

                  Qjac_k(var_UZ, var_UR) = Qjac_k(var_UZ, var_UR) + &
                    -visco_par_prof * Qjac_UZ_visco_par_UR_k &
                    -visco_pp_diff_prof * Qjac_UZ_visco_parpar_UR_k

                  Qjac_kn(var_UZ, var_UR) = Qjac_kn(var_UZ, var_UR) + &
                    -visco_par_prof * Qjac_UZ_visco_par_UR_kn &
                    -visco_pp_diff_prof * Qjac_UZ_visco_parpar_UR_kn

                  Qjac_pnn(var_UZ, var_UR) = Qjac_pnn(var_UZ, var_UR) + &
                    -visco_par_prof * Qjac_UZ_visco_par_UR_pnn &
                    -visco_pp_diff_prof * Qjac_UZ_visco_parpar_UR_pnn

                  Qjac_knn(var_UZ, var_UR) = Qjac_knn(var_UZ, var_UR) + &
                    -visco_par_prof * Qjac_UZ_visco_par_UR_knn &
                    -visco_pp_diff_prof * Qjac_UZ_visco_parpar_UR_knn
#endif

                  Qjac_p (var_UZ,var_UZ )  =  Qconv_UZ_UZ__p + visco_T * Qvisc_UZ_UZ__p - v * PneoZ_UZ - v * particle_source(ms,mt) * UZ &
                                            - visco_num * lap_Vstar * lap_bf
                  Qjac_n (var_UZ,var_UZ )  =  Qconv_UZ_UZ__n + visco_T * Qvisc_UZ_UZ__n
                  Qjac_k (var_UZ,var_UZ )  =                 + visco_T * Qvisc_UZ_UZ__k
                  Qjac_kn(var_UZ,var_UZ )  =                 + visco_T * Qvisc_UZ_UZ__kn

#ifdef VISCOSITY_TENSOR_FORM
                  Qjac_p (var_UZ,var_UZ ) = Qjac_p (var_UZ,var_UZ ) + v * Qgrad_visc_dot_grad_UZ__UZ__p
                  Qjac_n (var_UZ,var_UZ ) = Qjac_n (var_UZ,var_UZ ) + v * Qgrad_visc_dot_grad_UZ__UZ__n
#else
                  ! --- the missed term in full MHD viscosity
                  Qjac_p (var_UZ,var_UZ ) = Qjac_p (var_UZ,var_UZ ) - v * Qgrad_visc_dot_grad_UZ__UZ__p
                  Qjac_n (var_UZ,var_UZ ) = Qjac_n (var_UZ,var_UZ ) - v * Qgrad_visc_dot_grad_UZ__UZ__n
#ifdef VISCOSITY_INSIDE_LAPLACIAN
                  Qjac_p (var_UZ,var_UZ ) = Qjac_p (var_UZ,var_UZ ) + v * Qvisc_gradmu_dot_grad_vecV_UZ__UZ__p
                  Qjac_n (var_UZ,var_UZ ) = Qjac_n (var_UZ,var_UZ ) + v * Qvisc_gradmu_dot_grad_vecV_UZ__UZ__n
#endif
#ifdef VISCOSITY_COMPLETE_TERMS
                  Qjac_p (var_UZ,var_UZ ) = Qjac_p (var_UZ,var_UZ ) + v * Qgrad_visc_dot_gradU_T_UZ__UZ__p          &
                                + v * Q_minus_2_over_3_divU_grad_visc_UZ__UZ__p &
                                +     Q_1_over_3_visc_grad_divU_UZ__UZ__p
#endif
#endif                  

#ifdef ENABLE_VISCO_PAR_AND_PARPAR
                  Qjac_p(var_UZ, var_UZ) = Qjac_p(var_UZ, var_UZ) + &
                    -visco_par_prof * Qjac_UZ_visco_par_UZ_p &
                    -visco_pp_diff_prof * Qjac_UZ_visco_parpar_UZ_p

                  Qjac_n(var_UZ, var_UZ) = Qjac_n(var_UZ, var_UZ) + &
                    -visco_par_prof * Qjac_UZ_visco_par_UZ_n &
                    -visco_pp_diff_prof * Qjac_UZ_visco_parpar_UZ_n

                  Qjac_k(var_UZ, var_UZ) = Qjac_k(var_UZ, var_UZ) + &
                    -visco_par_prof * Qjac_UZ_visco_par_UZ_k &
                    -visco_pp_diff_prof * Qjac_UZ_visco_parpar_UZ_k

                  Qjac_kn(var_UZ, var_UZ) = Qjac_kn(var_UZ, var_UZ) + &
                    -visco_par_prof * Qjac_UZ_visco_par_UZ_kn &
                    -visco_pp_diff_prof * Qjac_UZ_visco_parpar_UZ_kn

                  Qjac_pnn(var_UZ, var_UZ) = Qjac_pnn(var_UZ, var_UZ) + &
                    -visco_par_prof * Qjac_UZ_visco_par_UZ_pnn &
                    -visco_pp_diff_prof * Qjac_UZ_visco_parpar_UZ_pnn

                  Qjac_knn(var_UZ, var_UZ) = Qjac_knn(var_UZ, var_UZ) + &
                    -visco_par_prof * Qjac_UZ_visco_par_UZ_knn &
                    -visco_pp_diff_prof * Qjac_UZ_visco_parpar_UZ_knn
#endif

                  Qjac_p (var_UZ,var_Up )  =  Qconv_UZ_Up__p + visco_T * Qvisc_UZ_Up__p
                  Qjac_n (var_UZ,var_Up )  =  Qconv_UZ_Up__n + visco_T * Qvisc_UZ_Up__n
                  Qjac_k (var_UZ,var_Up )  =                 + visco_T * Qvisc_UZ_Up__k
                  Qjac_kn(var_UZ,var_Up )  =                 + visco_T * Qvisc_UZ_Up__kn
#ifdef VISCOSITY_TENSOR_FORM
                  Qjac_p (var_UZ,var_Up ) = Qjac_p (var_UZ,var_Up ) + v * Qgrad_visc_dot_grad_UZ__Up__p
                  Qjac_n (var_UZ,var_Up ) = Qjac_n (var_UZ,var_Up ) + v * Qgrad_visc_dot_grad_UZ__Up__n
#else
#ifdef VISCOSITY_COMPLETE_TERMS
                  Qjac_p (var_UZ,var_Up ) = Qjac_p (var_UZ,var_Up ) + v * Qgrad_visc_dot_gradU_T_UZ__Up__p
                  Qjac_n (var_UZ,var_Up ) = Qjac_n (var_UZ,var_Up )             &
                                + v * Q_minus_2_over_3_divU_grad_visc_UZ__Up__n &
                                +     Q_1_over_3_visc_grad_divU_UZ__Up__n
#endif
#endif
#ifdef ENABLE_VISCO_PAR_AND_PARPAR
                  Qjac_p(var_UZ, var_Up) = Qjac_p(var_UZ, var_Up) + &
                    -visco_par_prof * Qjac_UZ_visco_par_Up_p &
                    -visco_pp_diff_prof * Qjac_UZ_visco_parpar_Up_p

                  Qjac_n(var_UZ, var_Up) = Qjac_n(var_UZ, var_Up) + &
                    -visco_par_prof * Qjac_UZ_visco_par_Up_n &
                    -visco_pp_diff_prof * Qjac_UZ_visco_parpar_Up_n

                  Qjac_k(var_UZ, var_Up) = Qjac_k(var_UZ, var_Up) + &
                    -visco_par_prof * Qjac_UZ_visco_par_Up_k &
                    -visco_pp_diff_prof * Qjac_UZ_visco_parpar_Up_k

                  Qjac_kn(var_UZ, var_Up) = Qjac_kn(var_UZ, var_Up) + &
                    -visco_par_prof * Qjac_UZ_visco_par_Up_kn &
                    -visco_pp_diff_prof * Qjac_UZ_visco_parpar_Up_kn

                  Qjac_pnn(var_UZ, var_Up) = Qjac_pnn(var_UZ, var_Up) + &
                    -visco_par_prof * Qjac_UZ_visco_par_Up_pnn &
                    -visco_pp_diff_prof * Qjac_UZ_visco_parpar_Up_pnn

                  Qjac_knn(var_UZ, var_Up) = Qjac_knn(var_UZ, var_Up) + &
                    -visco_par_prof * Qjac_UZ_visco_par_Up_knn &
                    -visco_pp_diff_prof * Qjac_UZ_visco_parpar_Up_knn
#endif

                  Qjac_p (var_UZ,var_rho)  =  Qconv_UZ_rho__p + v_Z * (rho*(Ti0+Te0)) - v * PneoZ_rho__p
                  Qjac_n (var_UZ,var_rho)  =  Qconv_UZ_rho__n                         - v * PneoZ_rho__n

                  Qjac_p (var_UZ,var_Ti )  =  Qconv_UZ_Ti__p  + v_Z * (rho0*Ti      ) - v * PneoZ_Ti__p
                  Qjac_n (var_UZ,var_Ti )  =  Qconv_UZ_Ti__n                          - v * PneoZ_Ti__n

                  Qjac_p (var_UZ,var_Te )  =                  + v_Z * (rho0*Te      ) + visco_T_T * Qvisc_UZ__p
                  Qjac_k (var_UZ,var_Te )  =                                          + visco_T_T * Qvisc_UZ__k

#ifdef VISCOSITY_TENSOR_FORM
                  Qjac_p (var_UZ,var_Te ) = Qjac_p (var_UZ,var_Te ) + v * Qgrad_visc_dot_grad_UZ__Te__p
                  Qjac_n (var_UZ,var_Te ) = Qjac_n (var_UZ,var_Te ) + v * Qgrad_visc_dot_grad_UZ__Te__n
#else
                  ! --- the missed term in full MHD viscosity
                  Qjac_p (var_UZ,var_Te ) = Qjac_p (var_UZ,var_Te ) - v * Qgrad_visc_dot_grad_UZ__Te__p
                  Qjac_n (var_UZ,var_Te ) = Qjac_n (var_UZ,var_Te ) - v * Qgrad_visc_dot_grad_UZ__Te__n
#ifdef VISCOSITY_INSIDE_LAPLACIAN
                  Qjac_p (var_UZ,var_Te ) = Qjac_p (var_UZ,var_Te ) + v * Qvisc_gradmu_dot_grad_vecV_UZ__Te__p
                  Qjac_n (var_UZ,var_Te ) = Qjac_n (var_UZ,var_Te ) + v * Qvisc_gradmu_dot_grad_vecV_UZ__Te__n
#endif
#ifdef VISCOSITY_COMPLETE_TERMS
                  Qjac_p (var_UZ,var_Te ) = Qjac_p (var_UZ,var_Te ) + v * Qgrad_visc_dot_gradU_T_UZ__Te__p          &
                                + v * Q_minus_2_over_3_divU_grad_visc_UZ__Te__p &
                                +     Q_1_over_3_visc_grad_divU_UZ__Te__p
                  Qjac_n (var_UZ,var_Te ) = Qjac_n (var_UZ,var_Te ) + v * Qgrad_visc_dot_gradU_T_UZ__Te__n
#endif
#endif                  

                  !###################################################################################################
                  !#  equation 6   (Phi component momentum equation)                                                 #
                  !###################################################################################################
#ifdef MOMENTUM_UPHI_IN_PARALLEL
                  Pjac(var_Up,var_UR)    = v * BR0 * rho0_corr * UR
                  Pjac(var_Up,var_UZ)    = v * BZ0 * rho0_corr * UZ
                  Pjac(var_Up,var_Up)    = v * Bp0 * rho0_corr * Up 

                  Qjac_p (var_Up,var_AR)  = + BR0_AR * Qconv_UR              &
                                            + Bp0_AR * Qconv_Up              &
                                            + BR0    * Qconv_UR_AR__p        &
                                            + BZ0    * Qconv_UZ_AR__p        &
                                            + Bp0    * Qconv_Up_AR__p        &
                                            + BR0_AR * visco_T * Qvisc_UR__p &
                                            + Bp0_AR * visco_T * Qvisc_Up__p &
                                            + p0 * BgradVstar_AR__p          &
                                            - v * (BR0*PneoR_AR__p + BZ0*PneoZ_AR__p) &
                                            - v * particle_source(ms,mt) * (Bp0_AR*Up0) &
                                            - visco_num * lap_Vstar * (BR0_AR * lap_UR               + Bp0_AR * lap_Up)
                  Qjac_n (var_Up,var_AR)  = + BZ0_AR__n * Qconv_UZ              &
                                            + BR0       * Qconv_UR_AR__n        &
                                            + BZ0       * Qconv_UZ_AR__n        &
                                            + Bp0       * Qconv_Up_AR__n        &
                                            + BZ0_AR__n * visco_T * Qvisc_UZ__p &
                                            + p0 * BgradVstar_AR__n             &
                                            - v * (BR0*PneoR_AR__n + BZ0*PneoZ_AR__n) &
                                            - v * (                + BZ0_AR__n*PneoZ) &
                                            - v * particle_source(ms,mt) * (BZ0_AR__n*UZ0) &
                                            - visco_num * lap_Vstar * (               + BZ0_AR__n * lap_UZ)
                  Qjac_k (var_Up,var_AR)  = + BR0_AR * visco_T * Qvisc_UR__k &
                                            + Bp0_AR * visco_T * Qvisc_Up__k &
                                            + p0 * BgradVstar_AR__k
                  Qjac_kn(var_Up,var_AR)  = + BZ0_AR__n * visco_T * Qvisc_UZ__k

#ifdef VISCOSITY_TENSOR_FORM
                  Qjac_p (var_Up,var_AR) = Qjac_p (var_Up,var_AR) + &
                                         - v * QdivU_BdotGradVisc__AR__p &
                                         - v * curl_U0_dot_gradmu_x_B0__AR__p &
                                         - v * visco_T * curl_U0_dot_curl_B0__AR__p
                  Qjac_n (var_Up,var_AR) = Qjac_n (var_Up,var_AR) + &
                                         - v * QdivU_BdotGradVisc__AR__n &
                                         - v * curl_U0_dot_gradmu_x_B0__AR__n &
                                         - v * visco_T * curl_U0_dot_curl_B0__AR__n
                  Qjac_pnn (var_Up,var_AR) = Qjac_pnn (var_Up,var_AR) + &
                                         - v * visco_T * curl_U0_dot_curl_B0__AR__nn             
#else
                  ! --- the missed term in full MHD viscosity, the B\cdot operation should be considerd when perform integration-by-parts
                  Qjac_p (var_Up,var_AR) = Qjac_p (var_Up,var_AR) + &
                                         - BR0_AR * v * Qgrad_visc_dot_grad_UR_p &
                                         - Bp0_AR * v * Qgrad_visc_dot_grad_Up_p &
                                         - visco_T * v * Qvisc_grad_BR_dot_grad_UR__AR__p   &
                                         - visco_T * v * Qvisc_grad_BZ_dot_grad_UZ__AR__p   &
                                         - visco_T * v * Qvisc_grad_Bp_dot_grad_Up__AR__p 

                  Qjac_n (var_Up,var_AR) = Qjac_n (var_Up,var_AR) + &
                                         - BZ0_AR__n * v * Qgrad_visc_dot_grad_UZ_p &
                                         - visco_T * v * Qvisc_grad_Bp_dot_grad_Up__AR__n 

                  Qjac_pnn (var_Up,var_AR) = Qjac_pnn (var_Up,var_AR) + &
                                         - visco_T * v * Qvisc_grad_BZ_dot_grad_UZ__AR__nn
#ifdef VISCOSITY_INSIDE_LAPLACIAN
                  Qjac_p (var_Up,var_AR) = Qjac_p (var_Up,var_AR) + &
                                         + BR0_AR * v * Qvisc_gradmu_dot_grad_vecV_UR_p &
                                         + Bp0_AR * v * Qvisc_gradmu_dot_grad_vecV_Up_p
                  Qjac_n (var_Up,var_AR) = Qjac_n (var_Up,var_AR) + &
                                         + BZ0_AR__n * v * Qvisc_gradmu_dot_grad_vecV_UZ_p
#endif
#ifdef VISCOSITY_COMPLETE_TERMS
                  Qjac_p (var_Up,var_AR) = Qjac_p (var_Up,var_AR) + &
                                + BR0_AR * v * Qgrad_visc_dot_gradU_T_UR__p &
                                + Bp0_AR * v * Qgrad_visc_dot_gradU_T_Up__p &
                                + BR0_AR * v * Q_minus_2_over_3_divU_grad_visc_UR__p &
                                + Bp0_AR * v * Q_minus_2_over_3_divU_grad_visc_Up__p &
                                + Q_1_over_3_B_dot_visc_grad_divU__AR__p
                  Qjac_n (var_Up,var_AR) = Qjac_n (var_Up,var_AR) + &
                                + BZ0_AR__n * v * Qgrad_visc_dot_gradU_T_UZ__p &
                                + BZ0_AR__n * v * Q_minus_2_over_3_divU_grad_visc_UZ__p &
                                + Q_1_over_3_B_dot_visc_grad_divU__AR__n
                  Qjac_k (var_Up,var_AR) = Qjac_k (var_Up,var_AR) + &
                                + Q_1_over_3_B_dot_visc_grad_divU__AR__k
#endif
#endif

#ifdef ENABLE_VISCO_PAR_AND_PARPAR
                  Qjac_p(var_Up, var_AR) = Qjac_p(var_Up, var_AR) + &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_AR_p + BZ0 * Qjac_UZ_visco_par_AR_p + Bp0 * Qjac_Up_visco_par_AR_p) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_AR_p + BZ0 * Qjac_UZ_visco_parpar_AR_p + Bp0 * Qjac_Up_visco_parpar_AR_p) &
                    -visco_par_prof * (BR0_AR * Qvec_UR_visco_par__p + Bp0_AR * Qvec_Up_visco_par__p) &
                    -visco_pp_diff_prof * (BR0_AR * Qvec_UR_visco_parpar__p + Bp0_AR * Qvec_Up_visco_parpar__p) &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_AR_for_Uphi_p + BZ0 * Qjac_UZ_visco_par_AR_for_Uphi_p + Bp0 * Qjac_Up_visco_par_AR_for_Uphi_p) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_AR_for_Uphi_p + BZ0 * Qjac_UZ_visco_parpar_AR_for_Uphi_p + Bp0 * Qjac_Up_visco_parpar_AR_for_Uphi_p) &
                    -visco_par_prof * (BR0_AR * Qvec_UR_visco_par_for_Uphi__p + Bp0_AR * Qvec_Up_visco_par_for_Uphi__p) &
                    -visco_pp_diff_prof * (BR0_AR * Qvec_UR_visco_parpar_for_Uphi__p + Bp0_AR * Qvec_Up_visco_parpar_for_Uphi__p)


                  Qjac_n(var_Up, var_AR) = Qjac_n(var_Up, var_AR) + &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_AR_n + BZ0 * Qjac_UZ_visco_par_AR_n + Bp0 * Qjac_Up_visco_par_AR_n) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_AR_n + BZ0 * Qjac_UZ_visco_parpar_AR_n + Bp0 * Qjac_Up_visco_parpar_AR_n) &
                    -visco_par_prof * (BZ0_AR__n * Qvec_UZ_visco_par__p) &
                    -visco_pp_diff_prof * (BZ0_AR__n * Qvec_UZ_visco_parpar__p) &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_AR_for_Uphi_n + BZ0 * Qjac_UZ_visco_par_AR_for_Uphi_n + Bp0 * Qjac_Up_visco_par_AR_for_Uphi_n) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_AR_for_Uphi_n + BZ0 * Qjac_UZ_visco_parpar_AR_for_Uphi_n + Bp0 * Qjac_Up_visco_parpar_AR_for_Uphi_n) &
                    -visco_par_prof * (BZ0_AR__n * Qvec_UZ_visco_par_for_Uphi__p) &
                    -visco_pp_diff_prof * (BZ0_AR__n * Qvec_UZ_visco_parpar_for_Uphi__p)

                  Qjac_k(var_Up, var_AR) = Qjac_k(var_Up, var_AR) + &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_AR_k + BZ0 * Qjac_UZ_visco_par_AR_k + Bp0 * Qjac_Up_visco_par_AR_k) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_AR_k + BZ0 * Qjac_UZ_visco_parpar_AR_k + Bp0 * Qjac_Up_visco_parpar_AR_k) &
                    -visco_par_prof * (BR0_AR * Qvec_UR_visco_par__k + Bp0_AR * Qvec_Up_visco_par__k) &
                    -visco_pp_diff_prof * (BR0_AR * Qvec_UR_visco_parpar__k + Bp0_AR * Qvec_Up_visco_parpar__k)                    

                  Qjac_kn(var_Up, var_AR) = Qjac_kn(var_Up, var_AR) + &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_AR_kn + BZ0 * Qjac_UZ_visco_par_AR_kn + Bp0 * Qjac_Up_visco_par_AR_kn) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_AR_kn + BZ0 * Qjac_UZ_visco_parpar_AR_kn + Bp0 * Qjac_Up_visco_parpar_AR_kn) &
                    -visco_par_prof * (BZ0_AR__n * Qvec_UZ_visco_par__k) &
                    -visco_pp_diff_prof * (BZ0_AR__n * Qvec_UZ_visco_parpar__k)

                  Qjac_pnn(var_Up, var_AR) = Qjac_pnn(var_Up, var_AR) + &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_AR_pnn + BZ0 * Qjac_UZ_visco_par_AR_pnn + Bp0 * Qjac_Up_visco_par_AR_pnn) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_AR_pnn + BZ0 * Qjac_UZ_visco_parpar_AR_pnn + Bp0 * Qjac_Up_visco_parpar_AR_pnn) &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_AR_for_Uphi_pnn + BZ0 * Qjac_UZ_visco_par_AR_for_Uphi_pnn + Bp0 * Qjac_Up_visco_par_AR_for_Uphi_pnn) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_AR_for_Uphi_pnn + BZ0 * Qjac_UZ_visco_parpar_AR_for_Uphi_pnn + Bp0 * Qjac_Up_visco_parpar_AR_for_Uphi_pnn)

                  Qjac_knn(var_Up, var_AR) = Qjac_knn(var_Up, var_AR) + &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_AR_knn + BZ0 * Qjac_UZ_visco_par_AR_knn + Bp0 * Qjac_Up_visco_par_AR_knn) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_AR_knn + BZ0 * Qjac_UZ_visco_parpar_AR_knn + Bp0 * Qjac_Up_visco_parpar_AR_knn)
#endif

                  Qjac_p (var_Up,var_AZ)  = + BZ0_AZ * Qconv_UZ              &
                                            + Bp0_AZ * Qconv_Up              &
                                            + BR0    * Qconv_UR_AZ__p        &
                                            + BZ0    * Qconv_UZ_AZ__p        &
                                            + Bp0    * Qconv_Up_AZ__p        &
                                            + BZ0_AZ * visco_T * Qvisc_UZ__p &
                                            + Bp0_AZ * visco_T * Qvisc_Up__p &
                                            + p0 * BgradVstar_AZ__p          &
                                            - v * (BR0*PneoR_AZ__p + BZ0*PneoZ_AZ__p) &
                                            - v * particle_source(ms,mt) * (Bp0_AZ*Up0) &
                                            - visco_num * lap_Vstar * (              + BZ0_AZ * lap_UZ + Bp0_AZ * lap_Up)
                  Qjac_n (var_Up,var_AZ)  = + BR0_AZ__n * Qconv_UR              &
                                            + BR0       * Qconv_UR_AZ__n        &
                                            + BZ0       * Qconv_UZ_AZ__n        &
                                            + Bp0       * Qconv_Up_AZ__n        &
                                            + BR0_AZ__n * visco_T * Qvisc_UR__p &
                                            + p0 * BgradVstar_AZ__n             &
                                            - v * (BR0*PneoR_AZ__n + BZ0*PneoZ_AZ__n) &
                                            - v * (BR0_AZ__n*PneoR                  ) &
                                            - v * particle_source(ms,mt) * (BR0_AZ__n*UR0) &
                                            - visco_num * lap_Vstar * (BR0_AZ__n * lap_UR )
                  Qjac_k (var_Up,var_AZ)  = + BZ0_AZ * visco_T * Qvisc_UZ__k &
                                            + Bp0_AZ * visco_T * Qvisc_Up__k &
                                            + p0 * BgradVstar_AZ__k
                  Qjac_kn(var_Up,var_AZ)  = + BR0_AZ__n * visco_T * Qvisc_UR__k

#ifdef VISCOSITY_TENSOR_FORM
                  Qjac_p (var_Up,var_AZ) = Qjac_p (var_Up,var_AZ) + &
                                         - v * QdivU_BdotGradVisc__AZ__p &
                                         - v * curl_U0_dot_gradmu_x_B0__AZ__p &
                                         - v * visco_T * curl_U0_dot_curl_B0__AZ__p
                  Qjac_n (var_Up,var_AZ) = Qjac_n (var_Up,var_AZ) + &
                                         - v * QdivU_BdotGradVisc__AZ__n &
                                         - v * curl_U0_dot_gradmu_x_B0__AZ__n &
                                         - v * visco_T * curl_U0_dot_curl_B0__AZ__n
                  Qjac_pnn (var_Up,var_AZ) = Qjac_pnn (var_Up,var_AZ) + &
                                         - v * visco_T * curl_U0_dot_curl_B0__AZ__nn
#else
                  ! --- the missed term in full MHD viscosity, the B\cdot operation should be considerd when perform integration-by-parts
                  Qjac_p (var_Up,var_AZ) = Qjac_p (var_Up,var_AZ) + &
                                         - BZ0_AZ * v * Qgrad_visc_dot_grad_UZ_p &
                                         - Bp0_AZ * v * Qgrad_visc_dot_grad_Up_p &
                                         - visco_T * v * Qvisc_grad_BZ_dot_grad_UZ__AZ__p   &
                                         - visco_T * v * Qvisc_grad_Bp_dot_grad_Up__AZ__p 

                  Qjac_n (var_Up,var_AZ) = Qjac_n (var_Up,var_AZ) + &
                                         - BR0_AZ__n * v * Qgrad_visc_dot_grad_UR_p &
                                         - visco_T * v * Qvisc_grad_BR_dot_grad_UR__AZ__n   &
                                         - visco_T * v * Qvisc_grad_Bp_dot_grad_Up__AZ__n

                  Qjac_pnn (var_Up,var_AZ) = Qjac_pnn (var_Up,var_AZ) + &
                                         - visco_T * v * Qvisc_grad_BR_dot_grad_UR__AZ__nn 
#ifdef VISCOSITY_INSIDE_LAPLACIAN
                  Qjac_p (var_Up,var_AZ) = Qjac_p (var_Up,var_AZ) + &
                                         + BZ0_AZ * v * Qvisc_gradmu_dot_grad_vecV_UZ_p &
                                         + Bp0_AZ * v * Qvisc_gradmu_dot_grad_vecV_Up_p
                  Qjac_n (var_Up,var_AZ) = Qjac_n (var_Up,var_AZ) + &
                                         + BR0_AZ__n * v * Qvisc_gradmu_dot_grad_vecV_UR_p 
#endif
#ifdef VISCOSITY_COMPLETE_TERMS
                  Qjac_p (var_Up,var_AZ) = Qjac_p (var_Up,var_AZ) + &
                                + BZ0_AZ * v * Qgrad_visc_dot_gradU_T_UZ__p &
                                + Bp0_AZ * v * Qgrad_visc_dot_gradU_T_Up__p &
                                + BZ0_AZ * v * Q_minus_2_over_3_divU_grad_visc_UZ__p &
                                + Bp0_AZ * v * Q_minus_2_over_3_divU_grad_visc_Up__p &
                                + Q_1_over_3_B_dot_visc_grad_divU__AZ__p
                  Qjac_n (var_Up,var_AZ) = Qjac_n (var_Up,var_AZ) + &
                                + BR0_AZ__n * v * Qgrad_visc_dot_gradU_T_UR__p &
                                + BR0_AZ__n * v * Q_minus_2_over_3_divU_grad_visc_UR__p &
                                + Q_1_over_3_B_dot_visc_grad_divU__AZ__n
                  Qjac_k (var_Up,var_AZ) = Qjac_k (var_Up,var_AZ) + &
                                + Q_1_over_3_B_dot_visc_grad_divU__AZ__k
#endif
#endif

#ifdef ENABLE_VISCO_PAR_AND_PARPAR
                  Qjac_p(var_Up, var_AZ) = Qjac_p(var_Up, var_AZ) + &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_AZ_p + BZ0 * Qjac_UZ_visco_par_AZ_p + Bp0 * Qjac_Up_visco_par_AZ_p) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_AZ_p + BZ0 * Qjac_UZ_visco_parpar_AZ_p + Bp0 * Qjac_Up_visco_parpar_AZ_p) &
                    -visco_par_prof * (BZ0_AZ * Qvec_UZ_visco_par__p + Bp0_AZ * Qvec_Up_visco_par__p) &
                    -visco_pp_diff_prof * (BZ0_AZ * Qvec_UZ_visco_parpar__p + Bp0_AZ * Qvec_Up_visco_parpar__p) &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_AZ_for_Uphi_p + BZ0 * Qjac_UZ_visco_par_AZ_for_Uphi_p + Bp0 * Qjac_Up_visco_par_AZ_for_Uphi_p) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_AZ_for_Uphi_p + BZ0 * Qjac_UZ_visco_parpar_AZ_for_Uphi_p + Bp0 * Qjac_Up_visco_parpar_AZ_for_Uphi_p) &
                    -visco_par_prof * (BZ0_AZ * Qvec_UZ_visco_par_for_Uphi__p + Bp0_AZ * Qvec_Up_visco_par_for_Uphi__p) &
                    -visco_pp_diff_prof * (BZ0_AZ * Qvec_UZ_visco_parpar_for_Uphi__p + Bp0_AZ * Qvec_Up_visco_parpar_for_Uphi__p)

                  Qjac_n(var_Up, var_AZ) = Qjac_n(var_Up, var_AZ) + &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_AZ_n + BZ0 * Qjac_UZ_visco_par_AZ_n + Bp0 * Qjac_Up_visco_par_AZ_n) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_AZ_n + BZ0 * Qjac_UZ_visco_parpar_AZ_n + Bp0 * Qjac_Up_visco_parpar_AZ_n) &
                    -visco_par_prof * (BR0_AZ__n * Qvec_UR_visco_par__p) &
                    -visco_pp_diff_prof * (BR0_AZ__n * Qvec_UR_visco_parpar__p) &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_AZ_for_Uphi_n + BZ0 * Qjac_UZ_visco_par_AZ_for_Uphi_n + Bp0 * Qjac_Up_visco_par_AZ_for_Uphi_n) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_AZ_for_Uphi_n + BZ0 * Qjac_UZ_visco_parpar_AZ_for_Uphi_n + Bp0 * Qjac_Up_visco_parpar_AZ_for_Uphi_n) &
                    -visco_par_prof * (BR0_AZ__n * Qvec_UR_visco_par_for_Uphi__p) &
                    -visco_pp_diff_prof * (BR0_AZ__n * Qvec_UR_visco_parpar_for_Uphi__p)

                  Qjac_k(var_Up, var_AZ) = Qjac_k(var_Up, var_AZ) + &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_AZ_k + BZ0 * Qjac_UZ_visco_par_AZ_k + Bp0 * Qjac_Up_visco_par_AZ_k) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_AZ_k + BZ0 * Qjac_UZ_visco_parpar_AZ_k + Bp0 * Qjac_Up_visco_parpar_AZ_k) &
                    -visco_par_prof * (BZ0_AZ * Qvec_UZ_visco_par__k + Bp0_AZ * Qvec_Up_visco_par__k) &
                    -visco_pp_diff_prof * (BZ0_AZ * Qvec_UZ_visco_parpar__k + Bp0_AZ * Qvec_Up_visco_parpar__k)

                  Qjac_kn(var_Up, var_AZ) = Qjac_kn(var_Up, var_AZ) + &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_AZ_kn + BZ0 * Qjac_UZ_visco_par_AZ_kn + Bp0 * Qjac_Up_visco_par_AZ_kn) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_AZ_kn + BZ0 * Qjac_UZ_visco_parpar_AZ_kn + Bp0 * Qjac_Up_visco_parpar_AZ_kn) &
                    -visco_par_prof * (BR0_AZ__n * Qvec_UR_visco_par__k) &
                    -visco_pp_diff_prof * (BR0_AZ__n * Qvec_UR_visco_parpar__k)

                  Qjac_pnn(var_Up, var_AZ) = Qjac_pnn(var_Up, var_AZ) + &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_AZ_pnn + BZ0 * Qjac_UZ_visco_par_AZ_pnn + Bp0 * Qjac_Up_visco_par_AZ_pnn) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_AZ_pnn + BZ0 * Qjac_UZ_visco_parpar_AZ_pnn + Bp0 * Qjac_Up_visco_parpar_AZ_pnn) &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_AZ_for_Uphi_pnn + BZ0 * Qjac_UZ_visco_par_AZ_for_Uphi_pnn + Bp0 * Qjac_Up_visco_par_AZ_for_Uphi_pnn) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_AZ_for_Uphi_pnn + BZ0 * Qjac_UZ_visco_parpar_AZ_for_Uphi_pnn + Bp0 * Qjac_Up_visco_parpar_AZ_for_Uphi_pnn)

                  Qjac_knn(var_Up, var_AZ) = Qjac_knn(var_Up, var_AZ) + &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_AZ_knn + BZ0 * Qjac_UZ_visco_par_AZ_knn + Bp0 * Qjac_Up_visco_par_AZ_knn) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_AZ_knn + BZ0 * Qjac_UZ_visco_parpar_AZ_knn + Bp0 * Qjac_Up_visco_parpar_AZ_knn)
#endif

                  Qjac_p (var_Up,var_A3)  = + BR0_A3 * Qconv_UR              &
                                            + BZ0_A3 * Qconv_UZ              &
                                            + Bp0_A3 * Qconv_Up              &
                                            + BR0    * Qconv_UR_A3__p        &
                                            + BZ0    * Qconv_UZ_A3__p        &
                                            + Bp0    * Qconv_Up_A3__p        &
                                            + BR0_A3 * visco_T * Qvisc_UR__p &
                                            + BZ0_A3 * visco_T * Qvisc_UZ__p &
                                            + Bp0_A3 * visco_T * Qvisc_Up__p &
                                            + p0 * BgradVstar_A3__p          &
                                            - v * (BR0*PneoR_A3__p + BZ0*PneoZ_A3__p) &
                                            - v * (BR0_A3*PneoR    + BZ0_A3*PneoZ   ) &
                                            - v * particle_source(ms,mt) * (BR0_A3*UR0 + BZ0_A3*UZ0) &
                                            - visco_num * lap_Vstar * (BR0_A3 * lap_UR + BZ0_A3 * lap_UZ + Bp0_A3 * lap_Up)
                  Qjac_n (var_Up,var_A3)  = + BR0    * Qconv_UR_A3__n        &
                                            + BZ0    * Qconv_UZ_A3__n        &
                                            + Bp0    * Qconv_Up_A3__n        &
                                            - v * (BR0*PneoR_A3__n + BZ0*PneoZ_A3__n)
                  Qjac_k (var_Up,var_A3)  = + BR0_A3 * visco_T * Qvisc_UR__k &
                                            + BZ0_A3 * visco_T * Qvisc_UZ__k &
                                            + Bp0_A3 * visco_T * Qvisc_Up__k &
                                            + p0 * BgradVstar_A3__k

#ifdef VISCOSITY_TENSOR_FORM
                  Qjac_p (var_Up,var_A3) = Qjac_p (var_Up,var_A3) + &
                                         - v * QdivU_BdotGradVisc__A3__p &
                                         - v * curl_U0_dot_gradmu_x_B0__A3__p &
                                         - v * visco_T * curl_U0_dot_curl_B0__A3__p
                  Qjac_n (var_Up,var_A3) = Qjac_n (var_Up,var_A3) + &
                                         - v * visco_T * curl_U0_dot_curl_B0__A3__n
#else
                  ! --- the missed term in full MHD viscosity, the B\cdot operation should be considerd when perform integration-by-parts
                  Qjac_p (var_Up,var_A3) = Qjac_p (var_Up,var_A3) + &
                                         - BR0_A3 * v * Qgrad_visc_dot_grad_UR_p &
                                         - BZ0_A3 * v * Qgrad_visc_dot_grad_UZ_p &
                                         - Bp0_A3 * v * Qgrad_visc_dot_grad_Up_p &
                                         - visco_T * v * Qvisc_grad_BR_dot_grad_UR__A3__p   &
                                         - visco_T * v * Qvisc_grad_BZ_dot_grad_UZ__A3__p   &
                                         - visco_T * v * Qvisc_grad_Bp_dot_grad_Up__A3__p 

                  Qjac_n (var_Up,var_A3) = Qjac_n (var_Up,var_A3) + &
                                         - visco_T * v * Qvisc_grad_BR_dot_grad_UR__A3__n   &
                                         - visco_T * v * Qvisc_grad_BZ_dot_grad_UZ__A3__n
#ifdef VISCOSITY_INSIDE_LAPLACIAN
                  Qjac_p (var_Up,var_A3) = Qjac_p (var_Up,var_A3) + &
                                         + BR0_A3 * v * Qvisc_gradmu_dot_grad_vecV_UR_p &
                                         + BZ0_A3 * v * Qvisc_gradmu_dot_grad_vecV_UZ_p &
                                         + Bp0_A3 * v * Qvisc_gradmu_dot_grad_vecV_Up_p
#endif
#ifdef VISCOSITY_COMPLETE_TERMS
                  Qjac_p (var_Up,var_A3) = Qjac_p (var_Up,var_A3) + &
                                + BR0_A3 * v * Qgrad_visc_dot_gradU_T_UR__p &
                                + BZ0_A3 * v * Qgrad_visc_dot_gradU_T_UZ__p &
                                + Bp0_A3 * v * Qgrad_visc_dot_gradU_T_Up__p &
                                + BR0_A3 * v * Q_minus_2_over_3_divU_grad_visc_UR__p &
                                + BZ0_A3 * v * Q_minus_2_over_3_divU_grad_visc_UZ__p &
                                + Bp0_A3 * v * Q_minus_2_over_3_divU_grad_visc_Up__p &
                                + Q_1_over_3_B_dot_visc_grad_divU__A3__p
                  Qjac_k (var_Up,var_A3) = Qjac_k (var_Up,var_A3) + &
                                + Q_1_over_3_B_dot_visc_grad_divU__A3__k
#endif
#endif

#ifdef ENABLE_VISCO_PAR_AND_PARPAR
                  Qjac_p(var_Up, var_A3) = Qjac_p(var_Up, var_A3) + &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_A3_p + BZ0 * Qjac_UZ_visco_par_A3_p + Bp0 * Qjac_Up_visco_par_A3_p) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_A3_p + BZ0 * Qjac_UZ_visco_parpar_A3_p + Bp0 * Qjac_Up_visco_parpar_A3_p) &
                    -visco_par_prof * (BR0_A3 * Qvec_UR_visco_par__p + BZ0_A3 * Qvec_UZ_visco_par__p + Bp0_A3 * Qvec_Up_visco_par__p) &
                    -visco_pp_diff_prof * (BR0_A3 * Qvec_UR_visco_parpar__p + BZ0_A3 * Qvec_UZ_visco_parpar__p + Bp0_A3 * Qvec_Up_visco_parpar__p) &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_A3_for_Uphi_p + BZ0 * Qjac_UZ_visco_par_A3_for_Uphi_p + Bp0 * Qjac_Up_visco_par_A3_for_Uphi_p) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_A3_for_Uphi_p + BZ0 * Qjac_UZ_visco_parpar_A3_for_Uphi_p + Bp0 * Qjac_Up_visco_parpar_A3_for_Uphi_p) &
                    -visco_par_prof * (BR0_A3 * Qvec_UR_visco_par_for_Uphi__p + BZ0_A3 * Qvec_UZ_visco_par_for_Uphi__p + Bp0_A3 * Qvec_Up_visco_par_for_Uphi__p) &
                    -visco_pp_diff_prof * (BR0_A3 * Qvec_UR_visco_parpar_for_Uphi__p + BZ0_A3 * Qvec_UZ_visco_parpar_for_Uphi__p + Bp0_A3 * Qvec_Up_visco_parpar_for_Uphi__p)

                  Qjac_n(var_Up, var_A3) = Qjac_n(var_Up, var_A3) + &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_A3_n + BZ0 * Qjac_UZ_visco_par_A3_n + Bp0 * Qjac_Up_visco_par_A3_n) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_A3_n + BZ0 * Qjac_UZ_visco_parpar_A3_n + Bp0 * Qjac_Up_visco_parpar_A3_n) &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_A3_for_Uphi_n + BZ0 * Qjac_UZ_visco_par_A3_for_Uphi_n + Bp0 * Qjac_Up_visco_par_A3_for_Uphi_n) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_A3_for_Uphi_n + BZ0 * Qjac_UZ_visco_parpar_A3_for_Uphi_n + Bp0 * Qjac_Up_visco_parpar_A3_for_Uphi_n) 

                  Qjac_k(var_Up, var_A3) = Qjac_k(var_Up, var_A3) + &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_A3_k + BZ0 * Qjac_UZ_visco_par_A3_k + Bp0 * Qjac_Up_visco_par_A3_k) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_A3_k + BZ0 * Qjac_UZ_visco_parpar_A3_k + Bp0 * Qjac_Up_visco_parpar_A3_k) &
                    -visco_par_prof * (BR0_A3 * Qvec_UR_visco_par__k + BZ0_A3 * Qvec_UZ_visco_par__k + Bp0_A3 * Qvec_Up_visco_par__k) &
                    -visco_pp_diff_prof * (BR0_A3 * Qvec_UR_visco_parpar__k + BZ0_A3 * Qvec_UZ_visco_parpar__k + Bp0_A3 * Qvec_Up_visco_parpar__k)

                  Qjac_kn(var_Up, var_A3) = Qjac_kn(var_Up, var_A3) + &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_A3_kn + BZ0 * Qjac_UZ_visco_par_A3_kn + Bp0 * Qjac_Up_visco_par_A3_kn) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_A3_kn + BZ0 * Qjac_UZ_visco_parpar_A3_kn + Bp0 * Qjac_Up_visco_parpar_A3_kn)

                  Qjac_pnn(var_Up, var_A3) = Qjac_pnn(var_Up, var_A3) + &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_A3_pnn + BZ0 * Qjac_UZ_visco_par_A3_pnn + Bp0 * Qjac_Up_visco_par_A3_pnn) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_A3_pnn + BZ0 * Qjac_UZ_visco_parpar_A3_pnn + Bp0 * Qjac_Up_visco_parpar_A3_pnn) &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_A3_for_Uphi_pnn + BZ0 * Qjac_UZ_visco_par_A3_for_Uphi_pnn + Bp0 * Qjac_Up_visco_par_A3_for_Uphi_pnn) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_A3_for_Uphi_pnn + BZ0 * Qjac_UZ_visco_parpar_A3_for_Uphi_pnn + Bp0 * Qjac_Up_visco_parpar_A3_for_Uphi_pnn)

                  Qjac_knn(var_Up, var_A3) = Qjac_knn(var_Up, var_A3) + &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_A3_knn + BZ0 * Qjac_UZ_visco_par_A3_knn + Bp0 * Qjac_Up_visco_par_A3_knn) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_A3_knn + BZ0 * Qjac_UZ_visco_parpar_A3_knn + Bp0 * Qjac_Up_visco_parpar_A3_knn)
#endif    

                  Qjac_p (var_Up,var_UR)  = + BR0 * Qconv_UR_UR__p           &
                                            + BZ0 * Qconv_UZ_UR__p           &
                                            + Bp0 * Qconv_Up_UR__p           &
                                            + BR0 * visco_T * Qvisc_UR_UR__p &
                                            + BZ0 * visco_T * Qvisc_UZ_UR__p &
                                            + Bp0 * visco_T * Qvisc_Up_UR__p &
                                            - v * (BR0*PneoR_UR + BZ0*PneoZ_UR) &
                                            - v * particle_source(ms,mt) * (BR0*UR) &
                                            - visco_num * lap_Vstar * BR0 * lap_bf
                  Qjac_n (var_Up,var_UR)  = + BR0 * Qconv_UR_UR__n           &
                                            + BZ0 * Qconv_UZ_UR__n           &
                                            + Bp0 * Qconv_Up_UR__n           &
                                            + BR0 * visco_T * Qvisc_UR_UR__n &
                                            + BZ0 * visco_T * Qvisc_UZ_UR__n &
                                            + Bp0 * visco_T * Qvisc_Up_UR__n
                  Qjac_k (var_Up,var_UR)  = + BR0 * visco_T * Qvisc_UR_UR__k &
                                            + BZ0 * visco_T * Qvisc_UZ_UR__k &
                                            + Bp0 * visco_T * Qvisc_Up_UR__k
                  Qjac_kn(var_Up,var_UR)  = + BR0 * visco_T * Qvisc_UR_UR__kn &
                                            + BZ0 * visco_T * Qvisc_UZ_UR__kn &
                                            + Bp0 * visco_T * Qvisc_Up_UR__kn

#ifdef VISCOSITY_TENSOR_FORM
                  Qjac_p (var_Up,var_UR) = Qjac_p (var_Up,var_UR) + &
                                         - v * QdivU_BdotGradVisc__UR__p      &
                                         - v * curl_U0_dot_gradmu_x_B0__UR__p &
                                         - v * visco_T * curl_U0_dot_curl_B0__UR__p
                  Qjac_n (var_Up,var_UR) = Qjac_n (var_Up,var_UR) + &
                                         - v * curl_U0_dot_gradmu_x_B0__UR__n &
                                         - v * visco_T * curl_U0_dot_curl_B0__UR__n
#else
                  ! --- the missed term in full MHD viscosity, the B\cdot operation should be considerd when perform integration-by-parts
                  Qjac_p (var_Up,var_UR) = Qjac_p (var_Up,var_UR) + &
                                         - BR0 * v * Qgrad_visc_dot_grad_UR__UR__p &
                                         - visco_T * v * Qvisc_grad_BR_dot_grad_UR__UR__p

                  Qjac_n (var_Up,var_UR) = Qjac_n (var_Up,var_UR) + &
                                         - BR0 * v * Qgrad_visc_dot_grad_UR__UR__n &
                                         - visco_T * v * Qvisc_grad_BR_dot_grad_UR__UR__n
#ifdef VISCOSITY_INSIDE_LAPLACIAN
                  Qjac_p (var_Up,var_UR) = Qjac_p (var_Up,var_UR) + &
                                         + BR0 * v * Qvisc_gradmu_dot_grad_vecV_UR__UR__p &
                                         + Bp0 * v * Qvisc_gradmu_dot_grad_vecV_Up__UR__p
                  Qjac_n (var_Up,var_UR) = Qjac_n (var_Up,var_UR) + &
                                         + BR0 * v * Qvisc_gradmu_dot_grad_vecV_UR__UR__n
#endif
#ifdef VISCOSITY_COMPLETE_TERMS
                  Qjac_p (var_Up,var_UR) = Qjac_p (var_Up,var_UR) + &
                                + BR0 * v * Qgrad_visc_dot_gradU_T_UR__UR__p &
                                + BZ0 * v * Qgrad_visc_dot_gradU_T_UZ__UR__p &
                                + Bp0 * v * Qgrad_visc_dot_gradU_T_Up__UR__p &
                                + BR0 * v * Q_minus_2_over_3_divU_grad_visc_UR__UR__p &
                                + BZ0 * v * Q_minus_2_over_3_divU_grad_visc_UZ__UR__p &
                                + Bp0 * v * Q_minus_2_over_3_divU_grad_visc_Up__UR__p &
                                + Q_1_over_3_B_dot_visc_grad_divU__UR__p
                  Qjac_n (var_Up,var_UR) = Qjac_n (var_Up,var_UR) + &
                                + Bp0 * v * Qgrad_visc_dot_gradU_T_Up__UR__n
                  Qjac_k (var_Up,var_UR) = Qjac_k (var_Up,var_UR) + &
                                + Q_1_over_3_B_dot_visc_grad_divU__UR__k
#endif
#endif

#ifdef ENABLE_VISCO_PAR_AND_PARPAR
                  Qjac_p(var_Up, var_UR) = Qjac_p(var_Up, var_UR) + &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_UR_p + BZ0 * Qjac_UZ_visco_par_UR_p + Bp0 * Qjac_Up_visco_par_UR_p) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_UR_p + BZ0 * Qjac_UZ_visco_parpar_UR_p + Bp0 * Qjac_Up_visco_parpar_UR_p) &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_UR_for_Uphi_p + BZ0 * Qjac_UZ_visco_par_UR_for_Uphi_p + Bp0 * Qjac_Up_visco_par_UR_for_Uphi_p) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_UR_for_Uphi_p + BZ0 * Qjac_UZ_visco_parpar_UR_for_Uphi_p + Bp0 * Qjac_Up_visco_parpar_UR_for_Uphi_p) 

                  Qjac_n(var_Up, var_UR) = Qjac_n(var_Up, var_UR) + &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_UR_n + BZ0 * Qjac_UZ_visco_par_UR_n + Bp0 * Qjac_Up_visco_par_UR_n) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_UR_n + BZ0 * Qjac_UZ_visco_parpar_UR_n + Bp0 * Qjac_Up_visco_parpar_UR_n) &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_UR_for_Uphi_n + BZ0 * Qjac_UZ_visco_par_UR_for_Uphi_n + Bp0 * Qjac_Up_visco_par_UR_for_Uphi_n) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_UR_for_Uphi_n + BZ0 * Qjac_UZ_visco_parpar_UR_for_Uphi_n + Bp0 * Qjac_Up_visco_parpar_UR_for_Uphi_n) 

                  Qjac_k(var_Up, var_UR) = Qjac_k(var_Up, var_UR) + &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_UR_k + BZ0 * Qjac_UZ_visco_par_UR_k + Bp0 * Qjac_Up_visco_par_UR_k) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_UR_k + BZ0 * Qjac_UZ_visco_parpar_UR_k + Bp0 * Qjac_Up_visco_parpar_UR_k)

                  Qjac_kn(var_Up, var_UR) = Qjac_kn(var_Up, var_UR) + &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_UR_kn + BZ0 * Qjac_UZ_visco_par_UR_kn + Bp0 * Qjac_Up_visco_par_UR_kn) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_UR_kn + BZ0 * Qjac_UZ_visco_parpar_UR_kn + Bp0 * Qjac_Up_visco_parpar_UR_kn)

                  Qjac_pnn(var_Up, var_UR) = Qjac_pnn(var_Up, var_UR) + &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_UR_pnn + BZ0 * Qjac_UZ_visco_par_UR_pnn + Bp0 * Qjac_Up_visco_par_UR_pnn) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_UR_pnn + BZ0 * Qjac_UZ_visco_parpar_UR_pnn + Bp0 * Qjac_Up_visco_parpar_UR_pnn) &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_UR_for_Uphi_pnn + BZ0 * Qjac_UZ_visco_par_UR_for_Uphi_pnn + Bp0 * Qjac_Up_visco_par_UR_for_Uphi_pnn) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_UR_for_Uphi_pnn + BZ0 * Qjac_UZ_visco_parpar_UR_for_Uphi_pnn + Bp0 * Qjac_Up_visco_parpar_UR_for_Uphi_pnn) 

                  Qjac_knn(var_Up, var_UR) = Qjac_knn(var_Up, var_UR) + &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_UR_knn + BZ0 * Qjac_UZ_visco_par_UR_knn + Bp0 * Qjac_Up_visco_par_UR_knn) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_UR_knn + BZ0 * Qjac_UZ_visco_parpar_UR_knn + Bp0 * Qjac_Up_visco_parpar_UR_knn)
#endif

                  Qjac_p (var_Up,var_UZ)  = + BR0 * Qconv_UR_UZ__p           &
                                            + BZ0 * Qconv_UZ_UZ__p           &
                                            + Bp0 * Qconv_Up_UZ__p           &
                                            + BR0 * visco_T * Qvisc_UR_UZ__p &
                                            + BZ0 * visco_T * Qvisc_UZ_UZ__p &
                                            + Bp0 * visco_T * Qvisc_Up_UZ__p &
                                            - v * (BR0*PneoR_UZ + BZ0*PneoZ_UZ) &
                                            - v * particle_source(ms,mt) * (BZ0*UZ) &
                                            - visco_num * lap_Vstar * BZ0 * lap_bf
                  Qjac_n (var_Up,var_UZ)  = + BR0 * Qconv_UR_UZ__n           &
                                            + BZ0 * Qconv_UZ_UZ__n           &
                                            + Bp0 * Qconv_Up_UZ__n           &
                                            + BR0 * visco_T * Qvisc_UR_UZ__n &
                                            + BZ0 * visco_T * Qvisc_UZ_UZ__n &
                                            + Bp0 * visco_T * Qvisc_Up_UZ__n
                  Qjac_k (var_Up,var_UZ)  = + BR0 * visco_T * Qvisc_UR_UZ__k &
                                            + BZ0 * visco_T * Qvisc_UZ_UZ__k &
                                            + Bp0 * visco_T * Qvisc_Up_UZ__k
                  Qjac_kn(var_Up,var_UZ)  = + BR0 * visco_T * Qvisc_UR_UZ__kn &
                                            + BZ0 * visco_T * Qvisc_UZ_UZ__kn &
                                            + Bp0 * visco_T * Qvisc_Up_UZ__kn
#ifdef VISCOSITY_TENSOR_FORM
                  Qjac_p (var_Up,var_UZ) = Qjac_p (var_Up,var_UZ) + &
                                         - v * QdivU_BdotGradVisc__UZ__p      &
                                         - v * curl_U0_dot_gradmu_x_B0__UZ__p &
                                         - v * visco_T * curl_U0_dot_curl_B0__UZ__p
                  Qjac_n (var_Up,var_UZ) = Qjac_n (var_Up,var_UZ) + &
                                         - v * curl_U0_dot_gradmu_x_B0__UZ__n &
                                         - v * visco_T * curl_U0_dot_curl_B0__UZ__n
#else
                  ! --- the missed term in full MHD viscosity, the B\cdot operation should be considerd when perform integration-by-parts
                  Qjac_p (var_Up,var_UZ) = Qjac_p (var_Up,var_UZ) + &
                                         - BZ0 * v * Qgrad_visc_dot_grad_UZ__UZ__p &
                                         - visco_T * v * Qvisc_grad_BZ_dot_grad_UZ__UZ__p

                  Qjac_n (var_Up,var_UZ) = Qjac_n (var_Up,var_UZ) + &
                                         - BZ0 * v * Qgrad_visc_dot_grad_UZ__UZ__n &
                                         - visco_T * v * Qvisc_grad_BZ_dot_grad_UZ__UZ__n
#ifdef VISCOSITY_INSIDE_LAPLACIAN
                  Qjac_p (var_Up,var_UZ) = Qjac_p (var_Up,var_UZ) + &
                                         + BZ0 * v * Qvisc_gradmu_dot_grad_vecV_UZ__UZ__p
                  Qjac_n (var_Up,var_UZ) = Qjac_n (var_Up,var_UZ) + &
                                         + BZ0 * v * Qvisc_gradmu_dot_grad_vecV_UZ__UZ__n
#endif
#ifdef VISCOSITY_COMPLETE_TERMS
                  Qjac_p (var_Up,var_UZ) = Qjac_p (var_Up,var_UZ) + &
                                + BR0 * v * Qgrad_visc_dot_gradU_T_UR__UZ__p &
                                + BZ0 * v * Qgrad_visc_dot_gradU_T_UZ__UZ__p &
                                + BR0 * v * Q_minus_2_over_3_divU_grad_visc_UR__UZ__p &
                                + BZ0 * v * Q_minus_2_over_3_divU_grad_visc_UZ__UZ__p &
                                + Bp0 * v * Q_minus_2_over_3_divU_grad_visc_Up__UZ__p &
                                + Q_1_over_3_B_dot_visc_grad_divU__UZ__p
                  Qjac_n (var_Up,var_UZ) = Qjac_n (var_Up,var_UZ) + &
                                + Bp0 * v * Qgrad_visc_dot_gradU_T_Up__UZ__n
                  Qjac_k (var_Up,var_UZ) = Qjac_k (var_Up,var_UZ) + &
                                + Q_1_over_3_B_dot_visc_grad_divU__UZ__k
#endif
#endif
#ifdef ENABLE_VISCO_PAR_AND_PARPAR
                  Qjac_p(var_Up, var_UZ) = Qjac_p(var_Up, var_UZ) + &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_UZ_p + BZ0 * Qjac_UZ_visco_par_UZ_p + Bp0 * Qjac_Up_visco_par_UZ_p) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_UZ_p + BZ0 * Qjac_UZ_visco_parpar_UZ_p + Bp0 * Qjac_Up_visco_parpar_UZ_p) &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_UZ_for_Uphi_p + BZ0 * Qjac_UZ_visco_par_UZ_for_Uphi_p + Bp0 * Qjac_Up_visco_par_UZ_for_Uphi_p) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_UZ_for_Uphi_p + BZ0 * Qjac_UZ_visco_parpar_UZ_for_Uphi_p + Bp0 * Qjac_Up_visco_parpar_UZ_for_Uphi_p) 

                  Qjac_n(var_Up, var_UZ) = Qjac_n(var_Up, var_UZ) + &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_UZ_n + BZ0 * Qjac_UZ_visco_par_UZ_n + Bp0 * Qjac_Up_visco_par_UZ_n) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_UZ_n + BZ0 * Qjac_UZ_visco_parpar_UZ_n + Bp0 * Qjac_Up_visco_parpar_UZ_n) &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_UZ_for_Uphi_n + BZ0 * Qjac_UZ_visco_par_UZ_for_Uphi_n + Bp0 * Qjac_Up_visco_par_UZ_for_Uphi_n) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_UZ_for_Uphi_n + BZ0 * Qjac_UZ_visco_parpar_UZ_for_Uphi_n + Bp0 * Qjac_Up_visco_parpar_UZ_for_Uphi_n) 

                  Qjac_k(var_Up, var_UZ) = Qjac_k(var_Up, var_UZ) + &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_UZ_k + BZ0 * Qjac_UZ_visco_par_UZ_k + Bp0 * Qjac_Up_visco_par_UZ_k) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_UZ_k + BZ0 * Qjac_UZ_visco_parpar_UZ_k + Bp0 * Qjac_Up_visco_parpar_UZ_k)

                  Qjac_kn(var_Up, var_UZ) = Qjac_kn(var_Up, var_UZ) + &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_UZ_kn + BZ0 * Qjac_UZ_visco_par_UZ_kn + Bp0 * Qjac_Up_visco_par_UZ_kn) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_UZ_kn + BZ0 * Qjac_UZ_visco_parpar_UZ_kn + Bp0 * Qjac_Up_visco_parpar_UZ_kn)

                  Qjac_pnn(var_Up, var_UZ) = Qjac_pnn(var_Up, var_UZ) + &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_UZ_pnn + BZ0 * Qjac_UZ_visco_par_UZ_pnn + Bp0 * Qjac_Up_visco_par_UZ_pnn) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_UZ_pnn + BZ0 * Qjac_UZ_visco_parpar_UZ_pnn + Bp0 * Qjac_Up_visco_parpar_UZ_pnn) &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_UZ_for_Uphi_pnn + BZ0 * Qjac_UZ_visco_par_UZ_for_Uphi_pnn + Bp0 * Qjac_Up_visco_par_UZ_for_Uphi_pnn) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_UZ_for_Uphi_pnn + BZ0 * Qjac_UZ_visco_parpar_UZ_for_Uphi_pnn + Bp0 * Qjac_Up_visco_parpar_UZ_for_Uphi_pnn) 

                  Qjac_knn(var_Up, var_UZ) = Qjac_knn(var_Up, var_UZ) + &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_UZ_knn + BZ0 * Qjac_UZ_visco_par_UZ_knn + Bp0 * Qjac_Up_visco_par_UZ_knn) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_UZ_knn + BZ0 * Qjac_UZ_visco_parpar_UZ_knn + Bp0 * Qjac_Up_visco_parpar_UZ_knn)
#endif

                  Qjac_p (var_Up,var_Up)  = + BR0 * Qconv_UR_Up__p           &
                                            + BZ0 * Qconv_UZ_Up__p           &
                                            + Bp0 * Qconv_Up_Up__p           &
                                            + BR0 * visco_T * Qvisc_UR_Up__p &
                                            + BZ0 * visco_T * Qvisc_UZ_Up__p &
                                            + Bp0 * visco_T * Qvisc_Up_Up__p &
                                            - v * particle_source(ms,mt) * (Bp0*Up) &
                                            - visco_num * lap_Vstar * Bp0 * lap_bf
                  Qjac_n (var_Up,var_Up)  = + BR0 * Qconv_UR_Up__n              &
                                            + BZ0 * Qconv_UZ_Up__n              &
                                            + Bp0 * Qconv_Up_Up__n              &
                                            + BR0 * visco_T * Qvisc_UR_Up__n &
                                            + BZ0 * visco_T * Qvisc_UZ_Up__n &
                                            + Bp0 * visco_T * Qvisc_Up_Up__n
                  Qjac_k (var_Up,var_Up)  = + BR0 * visco_T * Qvisc_UR_Up__k &
                                            + BZ0 * visco_T * Qvisc_UZ_Up__k &
                                            + Bp0 * visco_T * Qvisc_Up_Up__k
                  Qjac_kn(var_Up,var_Up)  = + BR0 * visco_T * Qvisc_UR_Up__kn &
                                            + BZ0 * visco_T * Qvisc_UZ_Up__kn &
                                            + Bp0 * visco_T * Qvisc_Up_Up__kn
#ifdef VISCOSITY_TENSOR_FORM
                  Qjac_p (var_Up,var_Up) = Qjac_p (var_Up,var_Up) + &
                                         - v * curl_U0_dot_gradmu_x_B0__Up__p &
                                         - v * visco_T * curl_U0_dot_curl_B0__Up__p
                  Qjac_n (var_Up,var_Up) = Qjac_n (var_Up,var_Up) + &
                                         - v * QdivU_BdotGradVisc__Up__n
#else
                  ! --- the missed term in full MHD viscosity, the B\cdot operation should be considerd when perform integration-by-parts
                  Qjac_p (var_Up,var_Up) = Qjac_p (var_Up,var_Up) + &
                                         - Bp0 * v * Qgrad_visc_dot_grad_Up__Up__p &
                                         - visco_T * v * Qvisc_grad_Bp_dot_grad_Up__Up__p

                  Qjac_n (var_Up,var_Up) = Qjac_n (var_Up,var_Up) + &
                                         - Bp0 * v * Qgrad_visc_dot_grad_Up__Up__n &
                                         - visco_T * v * Qvisc_grad_Bp_dot_grad_Up__Up__n
#ifdef VISCOSITY_INSIDE_LAPLACIAN
                  Qjac_p (var_Up,var_Up) = Qjac_p (var_Up,var_Up) + &
                                         + BR0 * v * Qvisc_gradmu_dot_grad_vecV_UR__Up__p &
                                         + Bp0 * v * Qvisc_gradmu_dot_grad_vecV_Up__Up__p
                  Qjac_n (var_Up,var_Up) = Qjac_n (var_Up,var_Up) + &
                                         + Bp0 * v * Qvisc_gradmu_dot_grad_vecV_Up__Up__n
#endif
#ifdef VISCOSITY_COMPLETE_TERMS
                  Qjac_p (var_Up,var_Up) = Qjac_p (var_Up,var_Up) + &
                                + BR0 * v * Qgrad_visc_dot_gradU_T_UR__Up__p &
                                + BZ0 * v * Qgrad_visc_dot_gradU_T_UZ__Up__p &
                                + Bp0 * v * Qgrad_visc_dot_gradU_T_Up__Up__p
                  Qjac_n (var_Up,var_Up) = Qjac_n (var_Up,var_Up) + &
                                + Bp0 * v * Qgrad_visc_dot_gradU_T_Up__Up__n &
                                + BR0 * v * Q_minus_2_over_3_divU_grad_visc_UR__Up__n &
                                + BZ0 * v * Q_minus_2_over_3_divU_grad_visc_UZ__Up__n &
                                + Bp0 * v * Q_minus_2_over_3_divU_grad_visc_Up__Up__n &
                                + Q_1_over_3_B_dot_visc_grad_divU__Up__n
                  Qjac_kn(var_Up,var_Up) = Qjac_kn(var_Up,var_Up) + &
                                + Q_1_over_3_B_dot_visc_grad_divU__Up__kn
#endif
#endif
#ifdef ENABLE_VISCO_PAR_AND_PARPAR
                  Qjac_p(var_Up, var_Up) = Qjac_p(var_Up, var_Up) + &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_Up_p + BZ0 * Qjac_UZ_visco_par_Up_p + Bp0 * Qjac_Up_visco_par_Up_p) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_Up_p + BZ0 * Qjac_UZ_visco_parpar_Up_p + Bp0 * Qjac_Up_visco_parpar_Up_p) &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_Up_for_Uphi_p + BZ0 * Qjac_UZ_visco_par_Up_for_Uphi_p + Bp0 * Qjac_Up_visco_par_Up_for_Uphi_p) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_Up_for_Uphi_p + BZ0 * Qjac_UZ_visco_parpar_Up_for_Uphi_p + Bp0 * Qjac_Up_visco_parpar_Up_for_Uphi_p) 

                  Qjac_n(var_Up, var_Up) = Qjac_n(var_Up, var_Up) + &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_Up_n + BZ0 * Qjac_UZ_visco_par_Up_n + Bp0 * Qjac_Up_visco_par_Up_n) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_Up_n + BZ0 * Qjac_UZ_visco_parpar_Up_n + Bp0 * Qjac_Up_visco_parpar_Up_n) &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_Up_for_Uphi_n + BZ0 * Qjac_UZ_visco_par_Up_for_Uphi_n + Bp0 * Qjac_Up_visco_par_Up_for_Uphi_n) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_Up_for_Uphi_n + BZ0 * Qjac_UZ_visco_parpar_Up_for_Uphi_n + Bp0 * Qjac_Up_visco_parpar_Up_for_Uphi_n) 

                  Qjac_k(var_Up, var_Up) = Qjac_k(var_Up, var_Up) + &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_Up_k + BZ0 * Qjac_UZ_visco_par_Up_k + Bp0 * Qjac_Up_visco_par_Up_k) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_Up_k + BZ0 * Qjac_UZ_visco_parpar_Up_k + Bp0 * Qjac_Up_visco_parpar_Up_k)

                  Qjac_kn(var_Up, var_Up) = Qjac_kn(var_Up, var_Up) + &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_Up_kn + BZ0 * Qjac_UZ_visco_par_Up_kn + Bp0 * Qjac_Up_visco_par_Up_kn) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_Up_kn + BZ0 * Qjac_UZ_visco_parpar_Up_kn + Bp0 * Qjac_Up_visco_parpar_Up_kn)

                  Qjac_pnn(var_Up, var_Up) = Qjac_pnn(var_Up, var_Up) + &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_Up_pnn + BZ0 * Qjac_UZ_visco_par_Up_pnn + Bp0 * Qjac_Up_visco_par_Up_pnn) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_Up_pnn + BZ0 * Qjac_UZ_visco_parpar_Up_pnn + Bp0 * Qjac_Up_visco_parpar_Up_pnn) &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_Up_for_Uphi_pnn + BZ0 * Qjac_UZ_visco_par_Up_for_Uphi_pnn + Bp0 * Qjac_Up_visco_par_Up_for_Uphi_pnn) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_Up_for_Uphi_pnn + BZ0 * Qjac_UZ_visco_parpar_Up_for_Uphi_pnn + Bp0 * Qjac_Up_visco_parpar_Up_for_Uphi_pnn) 

                  Qjac_knn(var_Up, var_Up) = Qjac_knn(var_Up, var_Up) + &
                    -visco_par_prof * (BR0 * Qjac_UR_visco_par_Up_knn + BZ0 * Qjac_UZ_visco_par_Up_knn + Bp0 * Qjac_Up_visco_par_Up_knn) &
                    -visco_pp_diff_prof * (BR0 * Qjac_UR_visco_parpar_Up_knn + BZ0 * Qjac_UZ_visco_parpar_Up_knn + Bp0 * Qjac_Up_visco_parpar_Up_knn)
#endif

                  Qjac_p (var_Up,var_rho) = + BR0 * Qconv_UR_rho__p              &
                                            + BZ0 * Qconv_UZ_rho__p              &
                                            + Bp0 * Qconv_Up_rho__p              &
                                            + (rho*(Ti0+Te0)) * BgradVstar__p    &
                                            - v * (BR0*PneoR_rho__p + BZ0*PneoZ_rho__p)
                  Qjac_n (var_Up,var_rho) = + BR0 * Qconv_UR_rho__n              &
                                            + BZ0 * Qconv_UZ_rho__n              &
                                            + Bp0 * Qconv_Up_rho__n              &
                                            - v * (BR0*PneoR_rho__n + BZ0*PneoZ_rho__n)
                  Qjac_k (var_Up,var_rho) = + (rho*(Ti0+Te0)) * BgradVstar__k

                  Qjac_p (var_Up,var_Ti)  = + BR0 * Qconv_UR_Ti__p              &
                                            + BZ0 * Qconv_UZ_Ti__p              &
                                            + Bp0 * Qconv_Up_Ti__p              &
                                            + (rho0*Ti) * BgradVstar__p         &
                                            - v * (BR0*PneoR_Ti__p + BZ0*PneoZ_Ti__p)
                  Qjac_n (var_Up,var_Ti)  = + BR0 * Qconv_UR_Ti__n              &
                                            + BZ0 * Qconv_UZ_Ti__n              &
                                            + Bp0 * Qconv_Up_Ti__n              &
                                            - v * (BR0*PneoR_Ti__n + BZ0*PneoZ_Ti__n)
                  Qjac_k (var_Up,var_Ti)  = + (rho0*Ti) * BgradVstar__k

                  Qjac_p (var_Up,var_Te)  = + BR0 * visco_T_T * Qvisc_UR__p &
                                            + BZ0 * visco_T_T * Qvisc_UZ__p &
                                            + Bp0 * visco_T_T * Qvisc_Up__p &
                                            + (rho0*Te) * BgradVstar__p
                  Qjac_k (var_Up,var_Te)  = + BR0 * visco_T_T * Qvisc_UR__k &
                                            + BZ0 * visco_T_T * Qvisc_UZ__k &
                                            + Bp0 * visco_T_T * Qvisc_Up__k &
                                            + (rho0*Te) * BgradVstar__k
#ifdef VISCOSITY_TENSOR_FORM
                  Qjac_p (var_Up,var_Te) = Qjac_p (var_Up,var_Te) + &
                                         - v * QdivU_BdotGradVisc__Te__p      &
                                         - v * curl_U0_dot_gradmu_x_B0__Te__p &
                                         - v * visco_T_T * curl_U0_dot_curl_B0
                  Qjac_n (var_Up,var_Te) = Qjac_n (var_Up,var_Te) + &
                                         - v * QdivU_BdotGradVisc__Te__n      &
                                         - v * curl_U0_dot_gradmu_x_B0__Te__n 
#else
                  ! --- the missed term in full MHD viscosity, the B\cdot operation should be considerd when perform integration-by-parts
                  Qjac_p (var_Up,var_Te) = Qjac_p (var_Up,var_Te) + &
                                         - BR0 * v * Qgrad_visc_dot_grad_UR__Te__p &
                                         - BZ0 * v * Qgrad_visc_dot_grad_UZ__Te__p &
                                         - Bp0 * v * Qgrad_visc_dot_grad_Up__Te__p &
                                         - visco_T_T * v * Qvisc_grad_BR_dot_grad_UR_p   &
                                         - visco_T_T * v * Qvisc_grad_BZ_dot_grad_UZ_p   &
                                         - visco_T_T * v * Qvisc_grad_Bp_dot_grad_Up_p

                  Qjac_n (var_Up,var_Te) = Qjac_n (var_Up,var_Te) + &
                                         - BR0 * v * Qgrad_visc_dot_grad_UR__Te__n &
                                         - BZ0 * v * Qgrad_visc_dot_grad_UZ__Te__n &
                                         - Bp0 * v * Qgrad_visc_dot_grad_Up__Te__n

#ifdef VISCOSITY_INSIDE_LAPLACIAN
                  Qjac_p (var_Up,var_Te) = Qjac_p (var_Up,var_Te) + &
                                         + BR0 * v * Qvisc_gradmu_dot_grad_vecV_UR__Te__p &
                                         + BZ0 * v * Qvisc_gradmu_dot_grad_vecV_UZ__Te__p &
                                         + Bp0 * v * Qvisc_gradmu_dot_grad_vecV_Up__Te__p
                  Qjac_n (var_Up,var_Te) = Qjac_n (var_Up,var_Te) + &
                                         + BR0 * v * Qvisc_gradmu_dot_grad_vecV_UR__Te__n &
                                         + BZ0 * v * Qvisc_gradmu_dot_grad_vecV_UZ__Te__n &
                                         + Bp0 * v * Qvisc_gradmu_dot_grad_vecV_Up__Te__n
#endif
#ifdef VISCOSITY_COMPLETE_TERMS
                  Qjac_p (var_Up,var_Te) = Qjac_p (var_Up,var_Te) + &
                                + BR0 * v * Qgrad_visc_dot_gradU_T_UR__Te__p &
                                + BZ0 * v * Qgrad_visc_dot_gradU_T_UZ__Te__p &
                                + Bp0 * v * Qgrad_visc_dot_gradU_T_Up__Te__p &
                                + BR0 * v * Q_minus_2_over_3_divU_grad_visc_UR__Te__p &
                                + BZ0 * v * Q_minus_2_over_3_divU_grad_visc_UZ__Te__p &
                                + Bp0 * v * Q_minus_2_over_3_divU_grad_visc_Up__Te__p &
                                + Q_1_over_3_B_dot_visc_grad_divU__Te__p
                  Qjac_n (var_Up,var_Te) = Qjac_n (var_Up,var_Te) + &
                                + BR0 * v * Qgrad_visc_dot_gradU_T_UR__Te__n &
                                + BZ0 * v * Qgrad_visc_dot_gradU_T_UZ__Te__n &
                                + Bp0 * v * Qgrad_visc_dot_gradU_T_Up__Te__n &
                                + Bp0 * v * Q_minus_2_over_3_divU_grad_visc_Up__Te__n &
                                + Q_1_over_3_B_dot_visc_grad_divU__Te__n
                  Qjac_k (var_Up,var_Te) = Qjac_k (var_Up,var_Te) + &
                                + Q_1_over_3_B_dot_visc_grad_divU__Te__k
#endif
#endif
#else

                  Pjac(var_Up,var_Up)    = v * rho0_corr * Up 

                  Qjac_p (var_Up,var_AR)  = Qconv_Up_AR__p + JxB_Up_AR__p
                  Qjac_n (var_Up,var_AR)  = Qconv_Up_AR__n + JxB_Up_AR__n
                  Qjac_k (var_Up,var_AR)  =                  JxB_Up_AR__k                                            
                  Qjac_kn(var_Up,var_AR)  =                  JxB_Up_AR__kn

#ifdef ENABLE_VISCO_PAR_AND_PARPAR
                  Qjac_p(var_Up, var_AR) = Qjac_p(var_Up, var_AR) + &
                    -visco_par_prof * (Qjac_Up_visco_par_AR_p) &
                    -visco_pp_diff_prof * (Qjac_Up_visco_parpar_AR_p)

                  Qjac_n(var_Up, var_AR) = Qjac_n(var_Up, var_AR) + &
                    -visco_par_prof * (Qjac_Up_visco_par_AR_n) &
                    -visco_pp_diff_prof * (Qjac_Up_visco_parpar_AR_n)

                  Qjac_k(var_Up, var_AR) = Qjac_k(var_Up, var_AR) + &
                    -visco_par_prof * (Qjac_Up_visco_par_AR_k) &
                    -visco_pp_diff_prof * (Qjac_Up_visco_parpar_AR_k)                 

                  Qjac_kn(var_Up, var_AR) = Qjac_kn(var_Up, var_AR) + &
                    -visco_par_prof * (Qjac_Up_visco_par_AR_kn) &
                    -visco_pp_diff_prof * (Qjac_Up_visco_parpar_AR_kn)

                  Qjac_pnn(var_Up, var_AR) = Qjac_pnn(var_Up, var_AR) + &
                    -visco_par_prof * (Qjac_Up_visco_par_AR_pnn) &
                    -visco_pp_diff_prof * (Qjac_Up_visco_parpar_AR_pnn)

                  Qjac_knn(var_Up, var_AR) = Qjac_knn(var_Up, var_AR) + &
                    -visco_par_prof * (Qjac_Up_visco_par_AR_knn) &
                    -visco_pp_diff_prof * (Qjac_Up_visco_parpar_AR_knn)
#endif

                  Qjac_p (var_Up,var_AZ)  = Qconv_Up_AZ__p + JxB_Up_AZ__p
                  Qjac_n (var_Up,var_AZ)  = Qconv_Up_AZ__n + JxB_Up_AZ__n
                  Qjac_k (var_Up,var_AZ)  =                  JxB_Up_AZ__k
                  Qjac_kn(var_Up,var_AZ)  =                  JxB_Up_AZ__kn
#ifdef ENABLE_VISCO_PAR_AND_PARPAR
                  Qjac_p(var_Up, var_AZ) = Qjac_p(var_Up, var_AZ) + &
                    -visco_par_prof * (Qjac_Up_visco_par_AZ_p) &
                    -visco_pp_diff_prof * (Qjac_Up_visco_parpar_AZ_p)

                  Qjac_n(var_Up, var_AZ) = Qjac_n(var_Up, var_AZ) + &
                    -visco_par_prof * (Qjac_Up_visco_par_AZ_n) &
                    -visco_pp_diff_prof * (Qjac_Up_visco_parpar_AZ_n) 

                  Qjac_k(var_Up, var_AZ) = Qjac_k(var_Up, var_AZ) + &
                    -visco_par_prof * (Qjac_Up_visco_par_AZ_k) &
                    -visco_pp_diff_prof * (Qjac_Up_visco_parpar_AZ_k)

                  Qjac_kn(var_Up, var_AZ) = Qjac_kn(var_Up, var_AZ) + &
                    -visco_par_prof * (Qjac_Up_visco_par_AZ_kn) &
                    -visco_pp_diff_prof * (Qjac_Up_visco_parpar_AZ_kn)

                  Qjac_pnn(var_Up, var_AZ) = Qjac_pnn(var_Up, var_AZ) + &
                    -visco_par_prof * (Qjac_Up_visco_par_AZ_pnn) &
                    -visco_pp_diff_prof * (Qjac_Up_visco_parpar_AZ_pnn) 

                  Qjac_knn(var_Up, var_AZ) = Qjac_knn(var_Up, var_AZ) + &
                    -visco_par_prof * (Qjac_Up_visco_par_AZ_knn) &
                    -visco_pp_diff_prof * (Qjac_Up_visco_parpar_AZ_knn)
#endif

                  Qjac_p (var_Up,var_A3)  = Qconv_Up_A3__p + JxB_Up_A3__p
                  Qjac_n (var_Up,var_A3)  = Qconv_Up_A3__n + JxB_Up_A3__n
                  Qjac_k (var_Up,var_A3)  =                  JxB_Up_A3__k
                  Qjac_kn(var_Up,var_A3)  =                  JxB_Up_A3__kn

#ifdef ENABLE_VISCO_PAR_AND_PARPAR
                  Qjac_p(var_Up, var_A3) = Qjac_p(var_Up, var_A3) + &
                    -visco_par_prof * (Qjac_Up_visco_par_A3_p) &
                    -visco_pp_diff_prof * (Qjac_Up_visco_parpar_A3_p) 

                  Qjac_n(var_Up, var_A3) = Qjac_n(var_Up, var_A3) + &
                    -visco_par_prof * (Qjac_Up_visco_par_A3_n) &
                    -visco_pp_diff_prof * (Qjac_Up_visco_parpar_A3_n)

                  Qjac_k(var_Up, var_A3) = Qjac_k(var_Up, var_A3) + &
                    -visco_par_prof * (Qjac_Up_visco_par_A3_k) &
                    -visco_pp_diff_prof * (Qjac_Up_visco_parpar_A3_k)

                  Qjac_kn(var_Up, var_A3) = Qjac_kn(var_Up, var_A3) + &
                    -visco_par_prof * (Qjac_Up_visco_par_A3_kn) &
                    -visco_pp_diff_prof * (Qjac_Up_visco_parpar_A3_kn)

                  Qjac_pnn(var_Up, var_A3) = Qjac_pnn(var_Up, var_A3) + &
                    -visco_par_prof * (Qjac_Up_visco_par_A3_pnn) &
                    -visco_pp_diff_prof * (Qjac_Up_visco_parpar_A3_pnn)

                  Qjac_knn(var_Up, var_A3) = Qjac_knn(var_Up, var_A3) + &
                    -visco_par_prof * (Qjac_Up_visco_par_A3_knn) &
                    -visco_pp_diff_prof * (Qjac_Up_visco_parpar_A3_knn)
#endif    

                  Qjac_p (var_Up,var_UR)  = Qconv_Up_UR__p + visco_T * Qvisc_Up_UR__p                                             
                  Qjac_n (var_Up,var_UR)  = Qconv_Up_UR__n + visco_T * Qvisc_Up_UR__n
                  Qjac_k (var_Up,var_UR)  =                  visco_T * Qvisc_Up_UR__k
                  Qjac_kn(var_Up,var_UR)  =                  visco_T * Qvisc_Up_UR__kn
#ifdef VISCOSITY_TENSOR_FORM
                  Qjac_p (var_Up,var_UR ) = Qjac_p (var_Up,var_UR ) + v * Qgrad_visc_dot_grad_Up__UR__p
                  Qjac_n (var_Up,var_UR ) = Qjac_n (var_Up,var_UR ) + v * Qgrad_visc_dot_grad_Up__UR__n
#else
#ifdef VISCOSITY_INSIDE_LAPLACIAN
                  Qjac_p (var_Up,var_UR ) = Qjac_p (var_Up,var_UR ) + v * Qvisc_gradmu_dot_grad_vecV_Up__UR__p
#endif
#ifdef VISCOSITY_COMPLETE_TERMS
                  Qjac_p (var_Up,var_UR ) = Qjac_p (var_Up,var_UR ) + v * Qgrad_visc_dot_gradU_T_Up__UR__p          &
                                                                    + v * Q_minus_2_over_3_divU_grad_visc_Up__UR__p &
                                                                    +     Q_1_over_3_visc_grad_divU_Up__UR__p
                  Qjac_n (var_Up,var_UR ) = Qjac_n (var_Up,var_UR ) + v * Qgrad_visc_dot_gradU_T_Up__UR__n
                  Qjac_k (var_Up,var_UR ) = Qjac_k (var_Up,var_UR )                                                 &
                                                                    +     Q_1_over_3_visc_grad_divU_Up__UR__k
#endif
#endif

#ifdef ENABLE_VISCO_PAR_AND_PARPAR
                  Qjac_p(var_Up, var_UR) = Qjac_p(var_Up, var_UR) + &
                    -visco_par_prof * (Qjac_Up_visco_par_UR_p) &
                    -visco_pp_diff_prof * (Qjac_Up_visco_parpar_UR_p) 

                  Qjac_n(var_Up, var_UR) = Qjac_n(var_Up, var_UR) + &
                    -visco_par_prof * (Qjac_Up_visco_par_UR_n) &
                    -visco_pp_diff_prof * (Qjac_Up_visco_parpar_UR_n) 

                  Qjac_k(var_Up, var_UR) = Qjac_k(var_Up, var_UR) + &
                    -visco_par_prof * (Qjac_Up_visco_par_UR_k) &
                    -visco_pp_diff_prof * (Qjac_Up_visco_parpar_UR_k)

                  Qjac_kn(var_Up, var_UR) = Qjac_kn(var_Up, var_UR) + &
                    -visco_par_prof * (Qjac_Up_visco_par_UR_kn) &
                    -visco_pp_diff_prof * (Qjac_Up_visco_parpar_UR_kn)

                  Qjac_pnn(var_Up, var_UR) = Qjac_pnn(var_Up, var_UR) + &
                    -visco_par_prof * (Qjac_Up_visco_par_UR_pnn) &
                    -visco_pp_diff_prof * (Qjac_Up_visco_parpar_UR_pnn) 

                  Qjac_knn(var_Up, var_UR) = Qjac_knn(var_Up, var_UR) + &
                    -visco_par_prof * (Qjac_Up_visco_par_UR_knn) &
                    -visco_pp_diff_prof * (Qjac_Up_visco_parpar_UR_knn)
#endif

                  Qjac_p (var_Up,var_UZ)  = Qconv_Up_UZ__p + visco_T * Qvisc_Up_UZ__p 
                  Qjac_n (var_Up,var_UZ)  = Qconv_Up_UZ__n + visco_T * Qvisc_Up_UZ__n
                  Qjac_k (var_Up,var_UZ)  =                  visco_T * Qvisc_Up_UZ__k
                  Qjac_kn(var_Up,var_UZ)  =                  visco_T * Qvisc_Up_UZ__kn
#ifdef VISCOSITY_TENSOR_FORM
                  Qjac_p (var_Up,var_UZ ) = Qjac_p (var_Up,var_UZ ) + v * Qgrad_visc_dot_grad_Up__UZ__p
                  Qjac_n (var_Up,var_UZ ) = Qjac_n (var_Up,var_UZ ) + v * Qgrad_visc_dot_grad_Up__UZ__n
#else
#ifdef VISCOSITY_COMPLETE_TERMS
                  Qjac_p (var_Up,var_UZ ) = Qjac_p (var_Up,var_UZ )                                                 &
                                                                    + v * Q_minus_2_over_3_divU_grad_visc_Up__UZ__p &
                                                                    +     Q_1_over_3_visc_grad_divU_Up__UZ__p
                  Qjac_n (var_Up,var_UZ ) = Qjac_n (var_Up,var_UZ ) + v * Qgrad_visc_dot_gradU_T_Up__UZ__n
                  Qjac_k (var_Up,var_UZ ) = Qjac_k (var_Up,var_UZ )                                                 &
                                                                    +     Q_1_over_3_visc_grad_divU_Up__UZ__k
#endif
#endif
#ifdef ENABLE_VISCO_PAR_AND_PARPAR
                  Qjac_p(var_Up, var_UZ) = Qjac_p(var_Up, var_UZ) + &
                    -visco_par_prof * (Qjac_Up_visco_par_UZ_p) &
                    -visco_pp_diff_prof * (Qjac_Up_visco_parpar_UZ_p) 

                  Qjac_n(var_Up, var_UZ) = Qjac_n(var_Up, var_UZ) + &
                    -visco_par_prof * (Qjac_Up_visco_par_UZ_n) &
                    -visco_pp_diff_prof * (Qjac_Up_visco_parpar_UZ_n) 

                  Qjac_k(var_Up, var_UZ) = Qjac_k(var_Up, var_UZ) + &
                    -visco_par_prof * (Qjac_Up_visco_par_UZ_k) &
                    -visco_pp_diff_prof * (Qjac_Up_visco_parpar_UZ_k)

                  Qjac_kn(var_Up, var_UZ) = Qjac_kn(var_Up, var_UZ) + &
                    -visco_par_prof * (Qjac_Up_visco_par_UZ_kn) &
                    -visco_pp_diff_prof * (Qjac_Up_visco_parpar_UZ_kn)

                  Qjac_pnn(var_Up, var_UZ) = Qjac_pnn(var_Up, var_UZ) + &
                    -visco_par_prof * (Qjac_Up_visco_par_UZ_pnn) &
                    -visco_pp_diff_prof * (Qjac_Up_visco_parpar_UZ_pnn) 

                  Qjac_knn(var_Up, var_UZ) = Qjac_knn(var_Up, var_UZ) + &
                    -visco_par_prof * (Qjac_Up_visco_par_UZ_knn) &
                    -visco_pp_diff_prof * (Qjac_Up_visco_parpar_UZ_knn)
#endif

                  Qjac_p (var_Up,var_Up)  = Qconv_Up_Up__p + visco_T * Qvisc_Up_Up__p &
                                            - v * particle_source(ms,mt) * Up &
                                            - visco_num * lap_Vstar * lap_bf
                  Qjac_n (var_Up,var_Up)  = Qconv_Up_Up__n + visco_T * Qvisc_Up_Up__n
                  Qjac_k (var_Up,var_Up)  =                  visco_T * Qvisc_Up_Up__k
                  Qjac_kn(var_Up,var_Up)  =                  visco_T * Qvisc_Up_Up__kn

#ifdef VISCOSITY_TENSOR_FORM
                  Qjac_p (var_Up,var_Up ) = Qjac_p (var_Up,var_Up ) + v * Qgrad_visc_dot_grad_Up__Up__p
                  Qjac_n (var_Up,var_Up ) = Qjac_n (var_Up,var_Up ) + v * Qgrad_visc_dot_grad_Up__Up__n
#else
                  Qjac_p (var_Up,var_Up) = Qjac_p (var_Up,var_Up) + &
                                           - v * Qgrad_visc_dot_grad_Up__Up__p
                  Qjac_n (var_Up,var_Up) = Qjac_n (var_Up,var_Up) + &
                                           - v * Qgrad_visc_dot_grad_Up__Up__n
#ifdef VISCOSITY_INSIDE_LAPLACIAN
                  Qjac_p (var_Up,var_Up ) = Qjac_p (var_Up,var_Up ) + v * Qvisc_gradmu_dot_grad_vecV_Up__Up__p
                  Qjac_n (var_Up,var_Up ) = Qjac_n (var_Up,var_Up ) + v * Qvisc_gradmu_dot_grad_vecV_Up__Up__n
#endif
#ifdef VISCOSITY_COMPLETE_TERMS
                  Qjac_p (var_Up,var_Up ) = Qjac_p (var_Up,var_Up ) + v * Qgrad_visc_dot_gradU_T_Up__Up__p
                  Qjac_n (var_Up,var_Up ) = Qjac_n (var_Up,var_Up ) + v * Qgrad_visc_dot_gradU_T_Up__Up__n          &
                                                                    + v * Q_minus_2_over_3_divU_grad_visc_Up__Up__n &
                                                                    +     Q_1_over_3_visc_grad_divU_Up__Up__n
                  Qjac_kn(var_Up,var_Up ) = Qjac_kn(var_Up,var_Up )                                                 &
                                                                    +     Q_1_over_3_visc_grad_divU_Up__Up__kn

#endif
#endif

#ifdef ENABLE_VISCO_PAR_AND_PARPAR
                  Qjac_p(var_Up, var_Up) = Qjac_p(var_Up, var_Up) + &
                    -visco_par_prof * (Qjac_Up_visco_par_Up_p) &
                    -visco_pp_diff_prof * (Qjac_Up_visco_parpar_Up_p) 

                  Qjac_n(var_Up, var_Up) = Qjac_n(var_Up, var_Up) + &
                    -visco_par_prof * (Qjac_Up_visco_par_Up_n) &
                    -visco_pp_diff_prof * (Qjac_Up_visco_parpar_Up_n) 

                  Qjac_k(var_Up, var_Up) = Qjac_k(var_Up, var_Up) + &
                    -visco_par_prof * (Qjac_Up_visco_par_Up_k) &
                    -visco_pp_diff_prof * (Qjac_Up_visco_parpar_Up_k)

                  Qjac_kn(var_Up, var_Up) = Qjac_kn(var_Up, var_Up) + &
                    -visco_par_prof * (Qjac_Up_visco_par_Up_kn) &
                    -visco_pp_diff_prof * (Qjac_Up_visco_parpar_Up_kn)

                  Qjac_pnn(var_Up, var_Up) = Qjac_pnn(var_Up, var_Up) + &
                    -visco_par_prof * (Qjac_Up_visco_par_Up_pnn) &
                    -visco_pp_diff_prof * (Qjac_Up_visco_parpar_Up_pnn) 

                  Qjac_knn(var_Up, var_Up) = Qjac_knn(var_Up, var_Up) + &
                    -visco_par_prof * (Qjac_Up_visco_par_Up_knn) &
                    -visco_pp_diff_prof * (Qjac_Up_visco_parpar_Up_knn)
#endif

                  Qjac_p (var_Up,var_rho) = Qconv_Up_rho__p  
                  Qjac_n (var_Up,var_rho) = Qconv_Up_rho__n
                  Qjac_k (var_Up,var_rho) =                 + (rho*(Ti0+Te0)) * v_p__k / R

                  Qjac_p (var_Up,var_Ti)  = Qconv_Up_Ti__p
                  Qjac_n (var_Up,var_Ti)  = Qconv_Up_Ti__n
                  Qjac_k (var_Up,var_Ti)  =                 + (rho0*Ti) * v_p__k / R

                  Qjac_p (var_Up,var_Te)  = visco_T_T * Qvisc_Up__p
                  Qjac_k (var_Up,var_Te)  = visco_T_T * Qvisc_Up__k + (rho0*Te) * v_p__k / R

#ifdef VISCOSITY_TENSOR_FORM
                  Qjac_p (var_Up,var_Te) = Qjac_p (var_Up,var_Te) + &
                                           + v * Qgrad_visc_dot_grad_Up__Te__p
                  Qjac_n (var_Up,var_Te) = Qjac_n (var_Up,var_Te) + &
                                           + v * Qgrad_visc_dot_grad_Up__Te__n
#else
                  Qjac_p (var_Up,var_Te) = Qjac_p (var_Up,var_Te) + &
                                           - v * Qgrad_visc_dot_grad_Up__Te__p
                  Qjac_n (var_Up,var_Te) = Qjac_n (var_Up,var_Te) + &
                                           - v * Qgrad_visc_dot_grad_Up__Te__n
#ifdef VISCOSITY_INSIDE_LAPLACIAN
                  Qjac_p (var_Up,var_Te ) = Qjac_p (var_Up,var_Te ) + v * Qvisc_gradmu_dot_grad_vecV_Up__Te__p
                  Qjac_n (var_Up,var_Te ) = Qjac_n (var_Up,var_Te ) + v * Qvisc_gradmu_dot_grad_vecV_Up__Te__n
#endif
#ifdef VISCOSITY_COMPLETE_TERMS
                  Qjac_p (var_Up,var_Te ) = Qjac_p (var_Up,var_Te ) + v * Qgrad_visc_dot_gradU_T_Up__Te__p          &
                                                                    + v * Q_minus_2_over_3_divU_grad_visc_Up__Te__p &
                                                                    +     Q_1_over_3_visc_grad_divU_Up__Te__p
                  Qjac_n (var_Up,var_Te ) = Qjac_n (var_Up,var_Te ) + v * Qgrad_visc_dot_gradU_T_Up__Te__n          &
                                                                    + v * Q_minus_2_over_3_divU_grad_visc_Up__Te__n &
                                                                    +     Q_1_over_3_visc_grad_divU_Up__Te__n
                  Qjac_k (var_Up,var_Te ) = Qjac_k (var_Up,var_Te )                                                 &
                                                                    +     Q_1_over_3_visc_grad_divU_Up__Te__k
#endif
#endif
#endif

                  !###################################################################################################
                  !#  equation 7   (Density equation)                                                                #
                  !###################################################################################################
                  Pjac   (var_rho,var_rho) =   v * rho

                  Qjac_p (var_rho,var_AR)  = + tauIC_rho * rho0 * VdiaGradVstar_AR__p &
                                             - (D_par_prof-D_prof) * BgradVstar_AR__p * BgradRho       / BB2 &
                                             - (D_par_prof-D_prof) * BgradVstar__p    * BgradRho_AR__p / BB2 &
                                             + (D_par_prof-D_prof) * BgradVstar__p    * BgradRho       / BB2**2 * BB2_AR__p
                  Qjac_n (var_rho,var_AR)  = + tauIC_rho * rho0 * VdiaGradVstar_AR__n &
                                             - (D_par_prof-D_prof) * BgradVstar_AR__n * BgradRho       / BB2 &
                                             - (D_par_prof-D_prof) * BgradVstar__p    * BgradRho_AR__n / BB2 &
                                             + (D_par_prof-D_prof) * BgradVstar__p    * BgradRho       / BB2**2 * BB2_AR__n
                  Qjac_k (var_rho,var_AR)  = + tauIC_rho * rho0 * VdiaGradVstar_AR__k &
                                             - (D_par_prof-D_prof) * BgradVstar_AR__k * BgradRho       / BB2 &
                                             - (D_par_prof-D_prof) * BgradVstar__k    * BgradRho_AR__p / BB2 &
                                             + (D_par_prof-D_prof) * BgradVstar__k    * BgradRho       / BB2**2 * BB2_AR__p
                  Qjac_kn(var_rho,var_AR)  = + tauIC_rho * rho0 * VdiaGradVstar_AR__kn &
                                             - (D_par_prof-D_prof) * BgradVstar__k    * BgradRho_AR__n / BB2 &
                                             + (D_par_prof-D_prof) * BgradVstar__k    * BgradRho       / BB2**2 * BB2_AR__n

                  Qjac_p (var_rho,var_AZ)  = + tauIC_rho * rho0 * VdiaGradVstar_AZ__p &
                                             - (D_par_prof-D_prof) * BgradVstar_AZ__p * BgradRho       / BB2 &
                                             - (D_par_prof-D_prof) * BgradVstar__p    * BgradRho_AZ__p / BB2 &
                                             + (D_par_prof-D_prof) * BgradVstar__p    * BgradRho       / BB2**2 * BB2_AZ__p
                  Qjac_n (var_rho,var_AZ)  = + tauIC_rho * rho0 * VdiaGradVstar_AZ__n &
                                             - (D_par_prof-D_prof) * BgradVstar_AZ__n * BgradRho       / BB2 &
                                             - (D_par_prof-D_prof) * BgradVstar__p    * BgradRho_AZ__n / BB2 &
                                             + (D_par_prof-D_prof) * BgradVstar__p    * BgradRho       / BB2**2 * BB2_AZ__n
                  Qjac_k (var_rho,var_AZ)  = + tauIC_rho * rho0 * VdiaGradVstar_AZ__k &
                                             - (D_par_prof-D_prof) * BgradVstar_AZ__k * BgradRho       / BB2 &
                                             - (D_par_prof-D_prof) * BgradVstar__k    * BgradRho_AZ__p / BB2 &
                                             + (D_par_prof-D_prof) * BgradVstar__k    * BgradRho       / BB2**2 * BB2_AZ__p
                  Qjac_kn(var_rho,var_AZ)  = + tauIC_rho * rho0 * VdiaGradVstar_AZ__kn &
                                             - (D_par_prof-D_prof) * BgradVstar__k    * BgradRho_AZ__n / BB2 &
                                             + (D_par_prof-D_prof) * BgradVstar__k    * BgradRho       / BB2**2 * BB2_AZ__n

                  Qjac_p (var_rho,var_A3)  = + tauIC_rho * rho0 * VdiaGradVstar_A3__p &
                                             - (D_par_prof-D_prof) * BgradVstar_A3__p * BgradRho    / BB2 &
                                             - (D_par_prof-D_prof) * BgradVstar__p    * BgradRho_A3 / BB2 &
                                             + (D_par_prof-D_prof) * BgradVstar__p    * BgradRho    / BB2**2 * BB2_A3
                  Qjac_n (var_rho,var_A3)  = + tauIC_rho * rho0 * VdiaGradVstar_A3__n
                  Qjac_k (var_rho,var_A3)  = + tauIC_rho * rho0 * VdiaGradVstar_A3__k &
                                             - (D_par_prof-D_prof) * BgradVstar_A3__k * BgradRho    / BB2 &
                                             - (D_par_prof-D_prof) * BgradVstar__k    * BgradRho_A3 / BB2 &
                                             + (D_par_prof-D_prof) * BgradVstar__k    * BgradRho    / BB2**2 * BB2_A3
                  Qjac_kn(var_rho,var_A3)  = + tauIC_rho * rho0 * VdiaGradVstar_A3__kn

                  Qjac_p (var_rho,var_UR)  = - v * ( rho0 * divU_UR + UgradRho_UR )

                  Qjac_p (var_rho,var_UZ)  = - v * ( rho0 * divU_UZ + UgradRho_UZ )

                  Qjac_p (var_rho,var_Up)  = - v * (                + UgradRho_Up )
                  Qjac_n (var_rho,var_Up)  = - v * ( rho0 * divU_Up__n )

                  Qjac_p (var_rho,var_rho) = - v * ( rho * divU + UgradRho_rho__p )                   &
                                             + tauIC_rho * rho  * VdiaGradVstar__p &
                                             + tauIC_rho * rho0 * VdiaGradVstar_rho__p &
                                             - D_prof * gradRho_gradVstar_rho__p                      &
                                             - (D_par_prof-D_prof) * BgradVstar__p * BgradRho_rho__p / BB2 &
                                             - D_perp_num * lap_Vstar * lap_bf
                  Qjac_n (var_rho,var_rho) = - v * (UgradRho_rho__n )                                 &
                                             + tauIC_rho * rho0 * VdiaGradVstar_rho__n &
                                             - (D_par_prof-D_prof) * BgradVstar__p * BgradRho_rho__n / BB2
                  Qjac_k (var_rho,var_rho) = + tauIC_rho * rho0 * VdiaGradVstar_rho__k &
                                             - (D_par_prof-D_prof) * BgradVstar__k * BgradRho_rho__p / BB2
                  Qjac_kn(var_rho,var_rho) = + tauIC_rho * rho0 * VdiaGradVstar_rho__kn &
                                             - D_prof * gradRho_gradVstar_rho__kn                     &
                                             - (D_par_prof-D_prof) * BgradVstar__k * BgradRho_rho__n / BB2

                  Qjac_p (var_rho,var_Ti ) = + tauIC_rho * rho0 * VdiaGradVstar_Ti__p
                  Qjac_n (var_rho,var_Ti ) = + tauIC_rho * rho0 * VdiaGradVstar_Ti__n
                  Qjac_k (var_rho,var_Ti ) = + tauIC_rho * rho0 * VdiaGradVstar_Ti__k
                  Qjac_kn(var_rho,var_Ti ) = + tauIC_rho * rho0 * VdiaGradVstar_Ti__kn

                  !###################################################################################################
                  !#  equation 8   (Ion Temperature  equation)                                                       #
                  !###################################################################################################
                  Pjac   (var_Ti,var_Ti)  =   v * rho0_corr * Ti
                  Pjac   (var_Ti,var_rho) =   v * rho       * Ti0_corr

                  Qjac_p (var_Ti,var_AR)  = - (ZKi_par_T-ZKi_prof) * BgradVstar_AR__p * BgradTi       / BB2                &
                                            - (ZKi_par_T-ZKi_prof) * BgradVstar__p    * BgradTi_AR__p / BB2                &
                                            + (ZKi_par_T-ZKi_prof) * BgradVstar__p    * BgradTi       / BB2**2 * BB2_AR__p
                  Qjac_n (var_Ti,var_AR)  = - (ZKi_par_T-ZKi_prof) * BgradVstar_AR__n * BgradTi       / BB2                &
                                            - (ZKi_par_T-ZKi_prof) * BgradVstar__p    * BgradTi_AR__n / BB2                &
                                            + (ZKi_par_T-ZKi_prof) * BgradVstar__p    * BgradTi       / BB2**2 * BB2_AR__n
                  Qjac_k (var_Ti,var_AR)  = - (ZKi_par_T-ZKi_prof) * BgradVstar_AR__k * BgradTi       / BB2                &
                                            - (ZKi_par_T-ZKi_prof) * BgradVstar__k    * BgradTi_AR__p / BB2                &
                                            + (ZKi_par_T-ZKi_prof) * BgradVstar__k    * BgradTi       / BB2**2 * BB2_AR__p
                  Qjac_kn(var_Ti,var_AR)  = - (ZKi_par_T-ZKi_prof) * BgradVstar__k    * BgradTi_AR__n / BB2                &
                                            + (ZKi_par_T-ZKi_prof) * BgradVstar__k    * BgradTi       / BB2**2 * BB2_AR__n

#ifdef ENABLE_DIA_DRIFT_TEMPERATURE
                  Qjac_p(var_Ti,var_AR) = Qjac_p(var_Ti,var_AR)                                                            &
                                          + tauIC_tem * (gamma-1.d0) * v                                                   &
                                          * (Ti0 * (VdiaR0_AR__p * rho0_R + VdiaZ0_AR__p * rho0_Z + VdiaP0_AR__p * rho0_p / R) &
                                          + rho0 * (VdiaR0_AR__p * Ti0_R + VdiaZ0_AR__p * Ti0_Z + VdiaP0_AR__p * Ti0_p / R))   &
                                          + tauIC_tem * gamma * rho0 * Ti0 * VdiaGradVstar_AR__p
                  Qjac_n(var_Ti,var_AR) = Qjac_n(var_Ti,var_AR)                                                            &
                                          + tauIC_tem * (gamma-1.d0) * v                                                   &
                                          * (Ti0 * (VdiaR0_AR__n * rho0_R + VdiaZ0_AR__n * rho0_Z + VdiaP0_AR__n * rho0_p / R) &
                                          + rho0 * (VdiaR0_AR__n * Ti0_R + VdiaZ0_AR__n * Ti0_Z + VdiaP0_AR__n * Ti0_p / R))   &
                                          + tauIC_tem * gamma * rho0 * Ti0 * VdiaGradVstar_AR__n
                  Qjac_k(var_Ti,var_AR) = Qjac_k(var_Ti,var_AR)                                                            &
                                          + tauIC_tem * gamma * rho0 * Ti0 * VdiaGradVstar_AR__k
                  Qjac_kn(var_Ti,var_AR) = Qjac_kn(var_Ti,var_AR)                                                          &
                                          + tauIC_tem * gamma * rho0 * Ti0 * VdiaGradVstar_AR__kn
#endif

                  Qjac_p (var_Ti,var_AZ)  = - (ZKi_par_T-ZKi_prof) * BgradVstar_AZ__p * BgradTi       / BB2                &
                                            - (ZKi_par_T-ZKi_prof) * BgradVstar__p    * BgradTi_AZ__p / BB2                &
                                            + (ZKi_par_T-ZKi_prof) * BgradVstar__p    * BgradTi       / BB2**2 * BB2_AZ__p
                  Qjac_n (var_Ti,var_AZ)  = - (ZKi_par_T-ZKi_prof) * BgradVstar_AZ__n * BgradTi       / BB2                &
                                            - (ZKi_par_T-ZKi_prof) * BgradVstar__p    * BgradTi_AZ__n / BB2                &
                                            + (ZKi_par_T-ZKi_prof) * BgradVstar__p    * BgradTi       / BB2**2 * BB2_AZ__n
                  Qjac_k (var_Ti,var_AZ)  = - (ZKi_par_T-ZKi_prof) * BgradVstar_AZ__k * BgradTi       / BB2                &
                                            - (ZKi_par_T-ZKi_prof) * BgradVstar__k    * BgradTi_AZ__p / BB2                &
                                            + (ZKi_par_T-ZKi_prof) * BgradVstar__k    * BgradTi       / BB2**2 * BB2_AZ__p
                  Qjac_kn(var_Ti,var_AZ)  = - (ZKi_par_T-ZKi_prof) * BgradVstar__k    * BgradTi_AZ__n / BB2                &
                                            + (ZKi_par_T-ZKi_prof) * BgradVstar__k    * BgradTi       / BB2**2 * BB2_AZ__n

#ifdef ENABLE_DIA_DRIFT_TEMPERATURE
                  Qjac_p(var_Ti,var_AZ) = Qjac_p(var_Ti,var_AZ)                                                            &
                                          + tauIC_tem * (gamma-1.d0) * v                                                   &
                                          * (Ti0 * (VdiaR0_AZ__p * rho0_R + VdiaZ0_AZ__p * rho0_Z + VdiaP0_AZ__p * rho0_p / R) &
                                          + rho0 * (VdiaR0_AZ__p * Ti0_R + VdiaZ0_AZ__p * Ti0_Z + VdiaP0_AZ__p * Ti0_p / R))   &
                                          + tauIC_tem * gamma * rho0 * Ti0 * VdiaGradVstar_AZ__p
                  Qjac_n(var_Ti,var_AZ) = Qjac_n(var_Ti,var_AZ)                                                            &
                                          + tauIC_tem * (gamma-1.d0) * v                                                   &
                                          * (Ti0 * (VdiaR0_AZ__n * rho0_R + VdiaZ0_AZ__n * rho0_Z + VdiaP0_AZ__n * rho0_p / R) &
                                          + rho0 * (VdiaR0_AZ__n * Ti0_R + VdiaZ0_AZ__n * Ti0_Z + VdiaP0_AZ__n * Ti0_p / R))   &
                                          + tauIC_tem * gamma * rho0 * Ti0 * VdiaGradVstar_AZ__n
                  Qjac_k(var_Ti,var_AZ) = Qjac_k(var_Ti,var_AZ)                                                            &
                                          + tauIC_tem * gamma * rho0 * Ti0 * VdiaGradVstar_AZ__k
                  Qjac_kn(var_Ti,var_AZ) = Qjac_kn(var_Ti,var_AZ)                                                          &
                                          + tauIC_tem * gamma * rho0 * Ti0 * VdiaGradVstar_AZ__kn
#endif

                  Qjac_p (var_Ti,var_A3)  = - (ZKi_par_T-ZKi_prof) * BgradVstar_A3__p * BgradTi    / BB2             &
                                            - (ZKi_par_T-ZKi_prof) * BgradVstar__p    * BgradTi_A3 / BB2             &
                                            + (ZKi_par_T-ZKi_prof) * BgradVstar__p    * BgradTi    / BB2**2 * BB2_A3
                  Qjac_k (var_Ti,var_A3)  = - (ZKi_par_T-ZKi_prof) * BgradVstar_A3__k * BgradTi    / BB2             &
                                            - (ZKi_par_T-ZKi_prof) * BgradVstar__k    * BgradTi_A3 / BB2             &
                                            + (ZKi_par_T-ZKi_prof) * BgradVstar__k    * BgradTi    / BB2**2 * BB2_A3
#ifdef ENABLE_DIA_DRIFT_TEMPERATURE
                  Qjac_p(var_Ti,var_A3) = Qjac_p(var_Ti,var_A3)                                                            &
                                          + tauIC_tem * (gamma-1.d0) * v                                                   &
                                          * (Ti0 * (VdiaR0_A3__p * rho0_R + VdiaZ0_A3__p * rho0_Z + VdiaP0_A3__p * rho0_p / R) &
                                          + rho0 * (VdiaR0_A3__p * Ti0_R + VdiaZ0_A3__p * Ti0_Z + VdiaP0_A3__p * Ti0_p / R))   &
                                          + tauIC_tem * gamma * rho0 * Ti0 * VdiaGradVstar_A3__p
                  Qjac_n(var_Ti,var_A3) = Qjac_n(var_Ti,var_A3)                                                            &
                                          + tauIC_tem * (gamma-1.d0) * v                                                   &
                                          * (Ti0 * (VdiaR0_A3__n * rho0_R + VdiaZ0_A3__n * rho0_Z + VdiaP0_A3__n * rho0_p / R) &
                                          + rho0 * (VdiaR0_A3__n * Ti0_R + VdiaZ0_A3__n * Ti0_Z + VdiaP0_A3__n * Ti0_p / R))   &
                                          + tauIC_tem * gamma * rho0 * Ti0 * VdiaGradVstar_A3__n
                  Qjac_k(var_Ti,var_A3) = Qjac_k(var_Ti,var_A3)                                                            &
                                          + tauIC_tem * gamma * rho0 * Ti0 * VdiaGradVstar_A3__k
                  Qjac_kn(var_Ti,var_A3) = Qjac_kn(var_Ti,var_A3)                                                          &
                                          + tauIC_tem * gamma * rho0 * Ti0 * VdiaGradVstar_A3__kn
#endif

                  Qjac_p (var_Ti,var_UR)  = + v * ( - rho0 * UgradTi_UR - Ti0 * UgradRho_UR  -  gamma * pi0 * divU_UR) &
                                            + v * (gamma-1.d0) * Qvisc_T_UR__p
                  Qjac_n (var_Ti,var_UR)  = + v * (gamma-1.d0) * Qvisc_T_UR__n

                  Qjac_p (var_Ti,var_UZ)  = + v * ( - rho0 * UgradTi_UZ - Ti0 * UgradRho_UZ  -  gamma * pi0 * divU_UZ) &
                                            + v * (gamma-1.d0) * Qvisc_T_UZ__p
                  Qjac_n (var_Ti,var_UZ)  = + v * (gamma-1.d0) * Qvisc_T_UZ__n

                  Qjac_p (var_Ti,var_Up)  = + v * ( - rho0 * UgradTi_Up - Ti0 * UgradRho_Up                             ) &
                                            + v * (gamma-1.d0) * Qvisc_T_Up__p
                  Qjac_n (var_Ti,var_Up)  = + v * (                                          -  gamma * pi0 * divU_Up__n) &
                                            + v * (gamma-1.d0) * Qvisc_T_Up__n

                  Qjac_p (var_Ti,var_rho) = + v * ( - rho * UgradTi - Ti0 * UgradRho_rho__p  -  gamma * (rho*Ti0) * divU) &
                                            + v * ddTi_e_drho
                  Qjac_n (var_Ti,var_rho) = + v * (                 - Ti0 * UgradRho_rho__n                             )

#ifdef ENABLE_DIA_DRIFT_TEMPERATURE
                  Qjac_p(var_Ti,var_rho) = Qjac_p(var_Ti,var_rho)                                                            &
                                          + tauIC_tem * (gamma-1.d0) * v                                                     &
                                          * (Ti0 * (VdiaR0_rho__p * rho0_R + VdiaZ0_rho__p * rho0_Z + VdiaP0_rho__p * rho0_p / R) &
                                          + rho0 * (VdiaR0_rho__p * Ti0_R + VdiaZ0_rho__p * Ti0_Z + VdiaP0_rho__p * Ti0_p / R))   &
                                          + tauIC_tem * (gamma-1.d0) * v                                                     &
                                          * (Ti0 * (VdiaR0 * rho_R + VdiaZ0 * rho_Z)                                         &
                                          + rho * (VdiaR0 * Ti0_R + VdiaZ0 * Ti0_Z + VdiaP0 * Ti0_p / R))                    &
                                          + tauIC_tem * gamma * rho0 * Ti0 * VdiaGradVstar_rho__p                            &
                                          + tauIC_tem * gamma * rho * Ti0 * VdiaGradVstar__p
                  Qjac_n(var_Ti,var_rho) = Qjac_n(var_Ti,var_rho)                                                            &
                                          + tauIC_tem * (gamma-1.d0) * v                                                     &
                                          * (Ti0 * (VdiaR0_rho__n * rho0_R + VdiaZ0_rho__n * rho0_Z + VdiaP0_rho__n * rho0_p / R) &
                                          + rho0 * (VdiaR0_rho__n * Ti0_R + VdiaZ0_rho__n * Ti0_Z + VdiaP0_rho__n * Ti0_p / R))   &
                                          + tauIC_tem * (gamma-1.d0) * v                                                     & 
                                          * (Ti0 * (VdiaP0 * rho_p__n / R))                                                  &                                          
                                          + tauIC_tem * gamma * rho0 * Ti0 * VdiaGradVstar_rho__n                          

                  Qjac_k(var_Ti,var_rho) = Qjac_k(var_Ti,var_rho)                                                            &
                                          + tauIC_tem * gamma * rho0 * Ti0 * VdiaGradVstar_rho__k                            &
                                          + tauIC_tem * gamma * rho * Ti0 * VdiaGradVstar__k
                  Qjac_kn(var_Ti,var_rho) = Qjac_kn(var_Ti,var_rho)                                                          &
                                          + tauIC_tem * gamma * rho0 * Ti0 * VdiaGradVstar_rho__kn
#endif

                  Qjac_p (var_Ti,var_Ti)  = + v * ( - rho0 * UgradTi_Ti__p - Ti * UgradRho - gamma * (rho0*Ti) * divU ) &
                                            + v * (gamma-1.d0) * Qvisc_T_T__p                                           &
                                            - ZKi_prof * gradTi_gradVstar_Ti__p                                         &
                                            - (ZKi_par_T-ZKi_prof) * BgradVstar__p * BgradTi_Ti__p / BB2                &
                                            - (dZKi_par_dT*Ti    ) * BgradVstar__p * BgradTi       / BB2                &
                                            + v * ddTi_e_dTi                                                            &
                                            - ZK_i_perp_num * lap_Vstar * lap_bf
                  Qjac_n (var_Ti,var_Ti)  = + v * ( - rho0 * UgradTi_Ti__n                                            ) &
                                            + v * (gamma-1.d0) * Qvisc_T_T__n                                           &
                                            - (ZKi_par_T-ZKi_prof) * BgradVstar__p * BgradTi_Ti__n / BB2
                  Qjac_k (var_Ti,var_Ti)  = - (ZKi_par_T-ZKi_prof) * BgradVstar__k * BgradTi_Ti__p / BB2                &
                                            - (dZKi_par_dT*Ti    ) * BgradVstar__k * BgradTi       / BB2
                  Qjac_kn(var_Ti,var_Ti)  = - ZKi_prof * gradTi_gradVstar_Ti__kn                                        &
                                            - (ZKi_par_T-ZKi_prof) * BgradVstar__k * BgradTi_Ti__n / BB2

#ifdef ENABLE_DIA_DRIFT_TEMPERATURE
                  Qjac_p(var_Ti,var_Ti) = Qjac_p(var_Ti,var_Ti)                                                              &
                                          + tauIC_tem * (gamma-1.d0) * v                                                     &
                                          * (Ti0 * (VdiaR0_Ti__p * rho0_R + VdiaZ0_Ti__p * rho0_Z + VdiaP0_Ti__p * rho0_p / R) &
                                          + rho0 * (VdiaR0_Ti__p * Ti0_R + VdiaZ0_Ti__p * Ti0_Z + VdiaP0_Ti__p * Ti0_p / R))   &
                                          + tauIC_tem * (gamma-1.d0) * v                                                     &
                                          * (Ti * (VdiaR0 * rho0_R + VdiaZ0 * rho0_Z + VdiaP0 * rho0_p / R)                  &
                                          + rho0 * (VdiaR0 * Ti_R + VdiaZ0 * Ti_Z))                                          &
                                          + tauIC_tem * gamma * rho0 * Ti0 * VdiaGradVstar_Ti__p                             &
                                          + tauIC_tem * gamma * rho0 * Ti * VdiaGradVstar__p
                  Qjac_n(var_Ti,var_Ti) = Qjac_n(var_Ti,var_Ti)                                                              &
                                          + tauIC_tem * (gamma-1.d0) * v                                                     &
                                          * (Ti0 * (VdiaR0_Ti__n * rho0_R + VdiaZ0_Ti__n * rho0_Z + VdiaP0_Ti__n * rho0_p / R) &
                                          + rho0 * (VdiaR0_Ti__n * Ti0_R + VdiaZ0_Ti__n * Ti0_Z + VdiaP0_Ti__n * Ti0_p / R))   &
                                          + tauIC_tem * (gamma-1.d0) * v                                                     &
                                          * (rho0 * (VdiaP0 * Ti_p__n / R))                                                  &
                                          + tauIC_tem * gamma * rho0 * Ti0 * VdiaGradVstar_Ti__n                                            

                  Qjac_k(var_Ti,var_Ti) = Qjac_k(var_Ti,var_Ti)                                                              &
                                          + tauIC_tem * gamma * rho0 * Ti0 * VdiaGradVstar_Ti__k                             &
                                          + tauIC_tem * gamma * rho0 * Ti * VdiaGradVstar__k
                  Qjac_kn(var_Ti,var_Ti) = Qjac_kn(var_Ti,var_Ti)                                                            &
                                          + tauIC_tem * gamma * rho0 * Ti0 * VdiaGradVstar_Ti__kn
#endif

                  Qjac_p (var_Ti,var_Te)  = + v * ddTi_e_dTe

                  !###################################################################################################
                  !#  equation 9   (Electron Temperature  equation)                                                  #
                  !###################################################################################################
                  Pjac   (var_Te,var_Te)  =   v * rho0_corr * Te
                  Pjac   (var_Te,var_rho) =   v * rho       * Te0_corr

                  Qjac_p (var_Te,var_AR)  = - (ZKe_par_T-ZKe_prof) * BgradVstar_AR__p * BgradTe       / BB2                &
                                            - (ZKe_par_T-ZKe_prof) * BgradVstar__p    * BgradTe_AR__p / BB2                &
                                            + (ZKe_par_T-ZKe_prof) * BgradVstar__p    * BgradTe       / BB2**2 * BB2_AR__p &
                                            + v * (gamma-1.0d0) * eta_T_ohm * JJ2_AR__p
                  Qjac_n (var_Te,var_AR)  = - (ZKe_par_T-ZKe_prof) * BgradVstar_AR__n * BgradTe       / BB2                &
                                            - (ZKe_par_T-ZKe_prof) * BgradVstar__p    * BgradTe_AR__n / BB2                &
                                            + (ZKe_par_T-ZKe_prof) * BgradVstar__p    * BgradTe       / BB2**2 * BB2_AR__n &
                                            + v * (gamma-1.0d0) * eta_T_ohm * JJ2_AR__n
                  Qjac_k (var_Te,var_AR)  = - (ZKe_par_T-ZKe_prof) * BgradVstar_AR__k * BgradTe       / BB2                &
                                            - (ZKe_par_T-ZKe_prof) * BgradVstar__k    * BgradTe_AR__p / BB2                &
                                            + (ZKe_par_T-ZKe_prof) * BgradVstar__k    * BgradTe       / BB2**2 * BB2_AR__p
                  Qjac_kn(var_Te,var_AR)  = - (ZKe_par_T-ZKe_prof) * BgradVstar__k    * BgradTe_AR__n / BB2                &
                                            + (ZKe_par_T-ZKe_prof) * BgradVstar__k    * BgradTe       / BB2**2 * BB2_AR__n
                  Qjac_pnn(var_Te,var_AR) = v * (gamma-1.0d0) * eta_T_ohm * JJ2_AR__nn

#ifdef ENABLE_DIA_DRIFT_TEMPERATURE
                  Qjac_p(var_Te,var_AR) = Qjac_p(var_Te,var_AR)                                                            &
                                          + tauIC_tem * (gamma-1.d0) * v                                                   &
                                          * (Te0 * (VdiaR0_AR__p * rho0_R + VdiaZ0_AR__p * rho0_Z + VdiaP0_AR__p * rho0_p / R) &
                                          + rho0 * (VdiaR0_AR__p * Te0_R + VdiaZ0_AR__p * Te0_Z + VdiaP0_AR__p * Te0_p / R))   &
                                          + tauIC_tem * gamma * rho0 * Te0 * VdiaGradVstar_AR__p
                  Qjac_n(var_Te,var_AR) = Qjac_n(var_Te,var_AR)                                                            &
                                          + tauIC_tem * (gamma-1.d0) * v                                                   &
                                          * (Te0 * (VdiaR0_AR__n * rho0_R + VdiaZ0_AR__n * rho0_Z + VdiaP0_AR__n * rho0_p / R) &
                                          + rho0 * (VdiaR0_AR__n * Te0_R + VdiaZ0_AR__n * Te0_Z + VdiaP0_AR__n * Te0_p / R))   &
                                          + tauIC_tem * gamma * rho0 * Te0 * VdiaGradVstar_AR__n
                  Qjac_k(var_Te,var_AR) = Qjac_k(var_Te,var_AR)                                                            &
                                          + tauIC_tem * gamma * rho0 * Te0 * VdiaGradVstar_AR__k
                  Qjac_kn(var_Te,var_AR) = Qjac_kn(var_Te,var_AR)                                                          &
                                          + tauIC_tem * gamma * rho0 * Te0 * VdiaGradVstar_AR__kn
#endif

                  Qjac_p (var_Te,var_AZ)  = - (ZKe_par_T-ZKe_prof) * BgradVstar_AZ__p * BgradTe       / BB2                &
                                            - (ZKe_par_T-ZKe_prof) * BgradVstar__p    * BgradTe_AZ__p / BB2                &
                                            + (ZKe_par_T-ZKe_prof) * BgradVstar__p    * BgradTe       / BB2**2 * BB2_AZ__p &
                                            + v * (gamma-1.0d0) * eta_T_ohm * JJ2_AZ__p
                  Qjac_n (var_Te,var_AZ)  = - (ZKe_par_T-ZKe_prof) * BgradVstar_AZ__n * BgradTe       / BB2                &
                                            - (ZKe_par_T-ZKe_prof) * BgradVstar__p    * BgradTe_AZ__n / BB2                &
                                            + (ZKe_par_T-ZKe_prof) * BgradVstar__p    * BgradTe       / BB2**2 * BB2_AZ__n &
                                            + v * (gamma-1.0d0) * eta_T_ohm * JJ2_AZ__n
                  Qjac_k (var_Te,var_AZ)  = - (ZKe_par_T-ZKe_prof) * BgradVstar_AZ__k * BgradTe       / BB2                &
                                            - (ZKe_par_T-ZKe_prof) * BgradVstar__k    * BgradTe_AZ__p / BB2                &
                                            + (ZKe_par_T-ZKe_prof) * BgradVstar__k    * BgradTe       / BB2**2 * BB2_AZ__p
                  Qjac_kn(var_Te,var_AZ)  = - (ZKe_par_T-ZKe_prof) * BgradVstar__k    * BgradTe_AZ__n / BB2                &
                                            + (ZKe_par_T-ZKe_prof) * BgradVstar__k    * BgradTe       / BB2**2 * BB2_AZ__n
                  Qjac_pnn(var_Te,var_AZ) = v * (gamma-1.0d0) * eta_T_ohm * JJ2_AZ__nn

#ifdef ENABLE_DIA_DRIFT_TEMPERATURE
                  Qjac_p(var_Te,var_AZ) = Qjac_p(var_Te,var_AZ)                                                            &
                                          + tauIC_tem * (gamma-1.d0) * v                                                   &
                                          * (Te0 * (VdiaR0_AZ__p * rho0_R + VdiaZ0_AZ__p * rho0_Z + VdiaP0_AZ__p * rho0_p / R) &
                                          + rho0 * (VdiaR0_AZ__p * Te0_R + VdiaZ0_AZ__p * Te0_Z + VdiaP0_AZ__p * Te0_p / R))   &
                                          + tauIC_tem * gamma * rho0 * Te0 * VdiaGradVstar_AZ__p
                  Qjac_n(var_Te,var_AZ) = Qjac_n(var_Te,var_AZ)                                                            &
                                          + tauIC_tem * (gamma-1.d0) * v                                                   &
                                          * (Te0 * (VdiaR0_AZ__n * rho0_R + VdiaZ0_AZ__n * rho0_Z + VdiaP0_AZ__n * rho0_p / R) &
                                          + rho0 * (VdiaR0_AZ__n * Te0_R + VdiaZ0_AZ__n * Te0_Z + VdiaP0_AZ__n * Te0_p / R))   &
                                          + tauIC_tem * gamma * rho0 * Te0 * VdiaGradVstar_AZ__n
                  Qjac_k(var_Te,var_AZ) = Qjac_k(var_Te,var_AZ)                                                            &
                                          + tauIC_tem * gamma * rho0 * Te0 * VdiaGradVstar_AZ__k
                  Qjac_kn(var_Te,var_AZ) = Qjac_kn(var_Te,var_AZ)                                                          &
                                          + tauIC_tem * gamma * rho0 * Te0 * VdiaGradVstar_AZ__kn
#endif

                  Qjac_p (var_Te,var_A3)  = - (ZKe_par_T-ZKe_prof) * BgradVstar_A3__p * BgradTe    / BB2             &
                                            - (ZKe_par_T-ZKe_prof) * BgradVstar__p    * BgradTe_A3 / BB2             &
                                            + (ZKe_par_T-ZKe_prof) * BgradVstar__p    * BgradTe    / BB2**2 * BB2_A3 &
                                            + v * (gamma-1.0d0) * eta_T_ohm * JJ2_A3__p
                  Qjac_k (var_Te,var_A3)  = - (ZKe_par_T-ZKe_prof) * BgradVstar_A3__k * BgradTe    / BB2             &
                                            - (ZKe_par_T-ZKe_prof) * BgradVstar__k    * BgradTe_A3 / BB2             &
                                            + (ZKe_par_T-ZKe_prof) * BgradVstar__k    * BgradTe    / BB2**2 * BB2_A3
                  Qjac_n (var_Te,var_A3)  = v * (gamma-1.0d0) * eta_T_ohm * JJ2_A3__n
                  Qjac_pnn(var_Te,var_A3) = v * (gamma-1.0d0) * eta_T_ohm * JJ2_A3__nn

#ifdef ENABLE_DIA_DRIFT_TEMPERATURE
                  Qjac_p(var_Te,var_A3) = Qjac_p(var_Te,var_A3)                                                            &
                                          + tauIC_tem * (gamma-1.d0) * v                                                   &
                                          * (Te0 * (VdiaR0_A3__p * rho0_R + VdiaZ0_A3__p * rho0_Z + VdiaP0_A3__p * rho0_p / R) &
                                          + rho0 * (VdiaR0_A3__p * Te0_R + VdiaZ0_A3__p * Te0_Z + VdiaP0_A3__p * Te0_p / R))   &
                                          + tauIC_tem * gamma * rho0 * Te0 * VdiaGradVstar_A3__p
                  Qjac_n(var_Te,var_A3) = Qjac_n(var_Te,var_A3)                                                            &
                                          + tauIC_tem * (gamma-1.d0) * v                                                   &
                                          * (Te0 * (VdiaR0_A3__n * rho0_R + VdiaZ0_A3__n * rho0_Z + VdiaP0_A3__n * rho0_p / R) &
                                          + rho0 * (VdiaR0_A3__n * Te0_R + VdiaZ0_A3__n * Te0_Z + VdiaP0_A3__n * Te0_p / R))   &
                                          + tauIC_tem * gamma * rho0 * Te0 * VdiaGradVstar_A3__n
                  Qjac_k(var_Te,var_A3) = Qjac_k(var_Te,var_A3)                                                            &
                                          + tauIC_tem * gamma * rho0 * Te0 * VdiaGradVstar_A3__k
                  Qjac_kn(var_Te,var_A3) = Qjac_kn(var_Te,var_A3)                                                          &
                                          + tauIC_tem * gamma * rho0 * Te0 * VdiaGradVstar_A3__kn
#endif

                  Qjac_p (var_Te,var_UR)  = + v * ( - rho0 * UgradTe_UR - Te0 * UgradRho_UR - gamma * pe0 * divU_UR  ) &
                                            + v * (gamma-1.d0) * Qvisc_T_UR__p
                  Qjac_n (var_Te,var_UR)  = + v * (gamma-1.d0) * Qvisc_T_UR__n

                  Qjac_p (var_Te,var_UZ)  = + v * ( - rho0 * UgradTe_UZ - Te0 * UgradRho_UZ - gamma * pe0 * divU_UZ  ) &
                                            + v * (gamma-1.d0) * Qvisc_T_UZ__p
                  Qjac_n (var_Te,var_UZ)  = + v * (gamma-1.d0) * Qvisc_T_UZ__n

                  Qjac_p (var_Te,var_Up)  = + v * ( - rho0 * UgradTe_Up - Te0 * UgradRho_Up                             ) &
                                            + v * (gamma-1.d0) * Qvisc_T_Up__p
                  Qjac_n (var_Te,var_Up)  = + v * (                                          -  gamma * pe0 * divU_Up__n) &
                                            + v * (gamma-1.d0) * Qvisc_T_Up__n

                  Qjac_p (var_Te,var_rho) = + v * ( - rho * UgradTe - Te0 * UgradRho_rho__p  - gamma * (rho*Te0) * divU ) &
                                            + v * ddTe_i_drho
                  Qjac_n (var_Te,var_rho) = + v * (                 - Te0 * UgradRho_rho__n                             )

#ifdef ENABLE_DIA_DRIFT_TEMPERATURE
                  Qjac_p(var_Te,var_rho) = Qjac_p(var_Te,var_rho)                                                            &
                                          + tauIC_tem * (gamma-1.d0) * v                                                     &
                                          * (Te0 * (VdiaR0_rho__p * rho0_R + VdiaZ0_rho__p * rho0_Z + VdiaP0_rho__p * rho0_p / R) &
                                          + rho0 * (VdiaR0_rho__p * Te0_R + VdiaZ0_rho__p * Te0_Z + VdiaP0_rho__p * Te0_p / R))   &
                                          + tauIC_tem * (gamma-1.d0) * v                                                     &
                                          * (Te0 * (VdiaR0 * rho_R + VdiaZ0 * rho_Z)                                         &
                                          + rho * (VdiaR0 * Te0_R + VdiaZ0 * Te0_Z + VdiaP0 * Te0_p / R))                    &
                                          + tauIC_tem * gamma * rho0 * Te0 * VdiaGradVstar_rho__p                            &
                                          + tauIC_tem * gamma * rho * Te0 * VdiaGradVstar__p
                  Qjac_n(var_Te,var_rho) = Qjac_n(var_Te,var_rho)                                                            &
                                          + tauIC_tem * (gamma-1.d0) * v                                                     &
                                          * (Te0 * (VdiaR0_rho__n * rho0_R + VdiaZ0_rho__n * rho0_Z + VdiaP0_rho__n * rho0_p / R) &
                                          + rho0 * (VdiaR0_rho__n * Te0_R + VdiaZ0_rho__n * Te0_Z + VdiaP0_rho__n * Te0_p / R))   &
                                          + tauIC_tem * (gamma-1.d0) * v                                                     & 
                                          * (Te0 * (VdiaP0 * rho_p__n / R))                                                  &                                          
                                          + tauIC_tem * gamma * rho0 * Te0 * VdiaGradVstar_rho__n                          

                  Qjac_k(var_Te,var_rho) = Qjac_k(var_Te,var_rho)                                                            &
                                          + tauIC_tem * gamma * rho0 * Te0 * VdiaGradVstar_rho__k                            &
                                          + tauIC_tem * gamma * rho * Te0 * VdiaGradVstar__k
                  Qjac_kn(var_Te,var_rho) = Qjac_kn(var_Te,var_rho)                                                          &
                                          + tauIC_tem * gamma * rho0 * Te0 * VdiaGradVstar_rho__kn
#endif

                  Qjac_p (var_Te,var_Te)  = + v * ( - rho0 * UgradTe_Te__p - Te * UgradRho - gamma * (rho0*Te) * divU ) &
                                            + v * (gamma-1.d0) * Qvisc_T_T__p                                           &
                                            - ZKe_prof * gradTe_gradVstar_Te__p                                         &
                                            - (ZKe_par_T-ZKe_prof) * BgradVstar__p * BgradTe_Te__p / BB2                &
                                            - (dZKe_par_dT*Te    ) * BgradVstar__p * BgradTe       / BB2                &
                                            + v * ddTe_i_dTe                                                            &
                                            + v * (gamma-1.0d0) * deta_dT_ohm * Te * JJ2                                &
                                            - ZK_e_perp_num * lap_Vstar * lap_bf
                  Qjac_n (var_Te,var_Te)  = + v * ( - rho0 * UgradTe_Te__n                                            ) &
                                            + v * (gamma-1.d0) * Qvisc_T_T__n                                           &
                                            - (ZKe_par_T-ZKe_prof) * BgradVstar__p * BgradTe_Te__n / BB2
                  Qjac_k (var_Te,var_Te)  = - (ZKe_par_T-ZKe_prof) * BgradVstar__k * BgradTe_Te__p / BB2                &
                                            - (dZKe_par_dT*Te    ) * BgradVstar__k * BgradTe       / BB2
                  Qjac_kn(var_Te,var_Te)  = - ZKe_prof * gradTe_gradVstar_Te__kn                                        &
                                           - (ZKe_par_T-ZKe_prof) * BgradVstar__k * BgradTe_Te__n / BB2
#ifdef ENABLE_DIA_DRIFT_TEMPERATURE
                  Qjac_p(var_Te,var_Te) = Qjac_p(var_Te,var_Te)                                                            &
                                          + tauIC_tem * (gamma-1.d0) * v                                                   &
                                          * (Te * (VdiaR0 * rho0_R + VdiaZ0 * rho0_Z + VdiaP0 * rho0_p / R)                &
                                          + rho0 * (VdiaR0 * Te_R + VdiaZ0 * Te_Z))                                        &
                                          + tauIC_tem * gamma * rho0 * Te * VdiaGradVstar__p
                  Qjac_n(var_Te,var_Te) = Qjac_n(var_Te,var_Te)                                                            &
                                          + tauIC_tem * (gamma-1.d0) * v                                                   &
                                          * (rho0 * (VdiaP0 * Te_p__n / R))

                  Qjac_k(var_Te,var_Te) = Qjac_k(var_Te,var_Te)                                                            &
                                          + tauIC_tem * gamma * rho0 * Te * VdiaGradVstar__k
#endif

                  Qjac_p (var_Te,var_Ti)  = + v * ddTe_i_dTi
#ifdef ENABLE_DIA_DRIFT_TEMPERATURE
                  Qjac_p(var_Te,var_Ti) = Qjac_p(var_Te,var_Ti)                                                              &
                                          + tauIC_tem * (gamma-1.d0) * v                                                     &
                                          * (Te0 * (VdiaR0_Ti__p * rho0_R + VdiaZ0_Ti__p * rho0_Z + VdiaP0_Ti__p * rho0_p / R) &
                                          + rho0 * (VdiaR0_Ti__p * Te0_R + VdiaZ0_Ti__p * Te0_Z + VdiaP0_Ti__p * Te0_p / R))   &
                                          + tauIC_tem * gamma * rho0 * Te0 * VdiaGradVstar_Ti__p
                  Qjac_n(var_Te,var_Ti) = Qjac_n(var_Te,var_Ti)                                                              &
                                          + tauIC_tem * (gamma-1.d0) * v                                                     &
                                          * (Te0 * (VdiaR0_Ti__n * rho0_R + VdiaZ0_Ti__n * rho0_Z + VdiaP0_Ti__n * rho0_p / R) &
                                          + rho0 * (VdiaR0_Ti__n * Te0_R + VdiaZ0_Ti__n * Te0_Z + VdiaP0_Ti__n * Te0_p / R))   &
                                          + tauIC_tem * gamma * rho0 * Te0 * VdiaGradVstar_Ti__n

                  Qjac_k(var_Te,var_Ti) = Qjac_k(var_Te,var_Ti)                                                              &
                                          + tauIC_tem * gamma * rho0 * Te0 * VdiaGradVstar_Ti__k
                  Qjac_kn(var_Te,var_Ti) = Qjac_kn(var_Te,var_Ti)                                                            &
                                          + tauIC_tem * gamma * rho0 * Te0 * VdiaGradVstar_Ti__kn
#endif

                  if (use_vms) call add_vms_to_elm()

                  if (use_fft) then
                    index_kl =       n_var*n_degrees*(k-1) +       n_var*(l-1) + 1
                    do kvar= 1,n_var
                      do ivar= 1,n_var
                        ij = ivar
                        kl = index_kl + (kvar-1)

                        amat(ivar,kvar)  = (1.d0+zeta)*Pjac(ivar,kvar) - tstep * theta * Qjac_p(ivar,kvar)
                        ELM_p(mp,kl,ij)  = ELM_p(mp,kl,ij)  +  wst * amat(ivar,kvar) * R * xjac

                        amat(ivar,kvar)  = - tstep * theta * Qjac_k(ivar,kvar)
                        ELM_k(mp,kl,ij)  = ELM_k(mp,kl,ij)  +  wst * amat(ivar,kvar) * R * xjac

                        amat(ivar,kvar)  = - tstep * theta * Qjac_n(ivar,kvar)
                        ELM_n(mp,kl,ij)  = ELM_n(mp,kl,ij)  +  wst * amat(ivar,kvar) * R * xjac

                        amat(ivar,kvar)  = - tstep * theta * Qjac_kn(ivar,kvar)
                        ELM_kn(mp,kl,ij) = ELM_kn(mp,kl,ij) +  wst * amat(ivar,kvar) * R * xjac

                        amat(ivar,kvar)  = - tstep * theta * Qjac_pnn(ivar,kvar)
                        ELM_pnn(mp,kl,ij)= ELM_pnn(mp,kl,ij) +  wst * amat(ivar,kvar) * R * xjac
#ifdef ENABLE_VISCO_PAR_AND_PARPAR
                        amat(ivar,kvar)  = - tstep * theta * Qjac_knn(ivar,kvar)
                        ELM_knn(mp,kl,ij)= ELM_knn(mp,kl,ij) +  wst * amat(ivar,kvar) * R * xjac
#endif
                      enddo
                    enddo
                  else
                    index_kl = n_tor_local*n_var*n_degrees*(k-1) + n_tor_local*n_var*(l-1) + in - n_tor_start +1
                    do kvar= 1,n_var
                      do ivar= 1,n_var
                        ij = index_ij + (ivar-1) * n_tor_local
                        kl = index_kl + (kvar-1) * n_tor_local

                        amat(ivar,kvar) = (1.d0+zeta)*Pjac(ivar,kvar) - tstep * theta * Qjac_p(ivar,kvar)
                        ELM(ij,kl)      = ELM(ij,kl) + wst * amat(ivar,kvar) * R * xjac

                        amat(ivar,kvar) = - tstep * theta * Qjac_k(ivar,kvar)
                        ELM(ij,kl)      = ELM(ij,kl) + wst * amat(ivar,kvar) * R * xjac

                        amat(ivar,kvar) = - tstep * theta * Qjac_n(ivar,kvar)
                        ELM(ij,kl)      = ELM(ij,kl) + wst * amat(ivar,kvar) * R * xjac

                        amat(ivar,kvar) = - tstep * theta * Qjac_kn(ivar,kvar)
                        ELM(ij,kl)      = ELM(ij,kl) + wst * amat(ivar,kvar) * R * xjac

                        amat(ivar,kvar) = - tstep * theta * Qjac_pnn(ivar,kvar)
                        ELM(ij,kl)      = ELM(ij,kl) + wst * amat(ivar,kvar) * R * xjac
#ifdef ENABLE_VISCO_PAR_AND_PARPAR
                        amat(ivar,kvar) = - tstep * theta * Qjac_knn(ivar,kvar)
                        ELM(ij,kl)      = ELM(ij,kl) + wst * amat(ivar,kvar) * R * xjac
#endif
                      enddo
                    enddo
                  endif

                enddo ! in loop (n_tor, or not...)

              enddo ! l loop n_degrees
            enddo ! k loop (n_vertex)

          enddo ! im loop (n_tor, or not...)

        enddo ! mp loop (n_plane)

      enddo ! ms loop
    enddo ! mt loop


    if (use_fft) then
      do i_v = 1, n_var
        do j_loc=1, n_vertex_max*n_var*n_degrees
          i_loc = n_var*n_degrees*(i-1) + n_var * (j-1) + i_v 
          in_fft =  ELM_p(1:n_plane,j_loc,i_v)
#ifdef USE_FFTW
          call dfftw_execute_dft_r2c(fftw_plan, in_fft, out_fft)
#else
          call my_fft(in_fft, out_fft, n_plane)
#endif
          do m=1,(n_tor+1)/2

            index_m = n_tor*(j_loc-1) + max(2*(m-1),1)

            do k=1,(n_tor+1)/2

              index_k = n_tor*(i_loc-1) + max(2*(k-1),1)

              l = (k-1) + (m-1)

              if (l .le. n_plane/2) then
                ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + real(out_fft(l+1))
                ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   - imag(out_fft(l+1))
                ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) - imag(out_fft(l+1))
                ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) - real(out_fft(l+1))
              endif

              l = abs((k-1) - (m-1)); sgn_l = sign(1, (k-m))

              if (l .le. n_plane/2) then
                ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   +         real(out_fft(l+1))
                ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   - sgn_l * imag(out_fft(l+1))
                ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + sgn_l * imag(out_fft(l+1))
                ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) +         real(out_fft(l+1))
              endif

            enddo

          enddo

          if (maxval(abs(ELM_n(1:n_plane,j_loc, i_v))) .ne. 0.d0) then

            in_fft =  ELM_n(1:n_plane,j_loc, i_v)
#ifdef USE_FFTW
            call dfftw_execute_dft_r2c(fftw_plan, in_fft, out_fft)
#else
            call my_fft(in_fft, out_fft, n_plane)
#endif

            do m=1,(n_tor+1)/2
              im = max(2*(m-1),1)
              index_m = n_tor*(j_loc-1) + max(2*(m-1),1)

              do k=1,(n_tor+1)/2

                index_k = n_tor*(i_loc-1) + max(2*(k-1),1)

                l = (k-1) + (m-1)

                if (l .le. n_plane/2) then
                  ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + imag(out_fft(l+1)) * float(mode(im))
                  ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + real(out_fft(l+1)) * float(mode(im))
                  ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + real(out_fft(l+1)) * float(mode(im))
                  ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) - imag(out_fft(l+1)) * float(mode(im))
                endif

                l = abs((k-1) - (m-1)); sgn_l = sign(1, (k-m))

                if (l .le. n_plane/2) then
                  ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   - sgn_l * imag(out_fft(l+1)) * float(mode(im))
                  ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   -         real(out_fft(l+1)) * float(mode(im))
                  ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) +         real(out_fft(l+1)) * float(mode(im))
                  ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) - sgn_l * imag(out_fft(l+1)) * float(mode(im))
                endif

              enddo

            enddo

          endif

          if (maxval(abs(ELM_k(1:n_plane,j_loc, i_v))) .ne. 0.d0) then

            in_fft =  ELM_k(1:n_plane,j_loc, i_v)

#ifdef USE_FFTW
            call dfftw_execute_dft_r2c(fftw_plan, in_fft, out_fft)
#else
            call my_fft(in_fft, out_fft, n_plane)
#endif

            do m=1,(n_tor+1)/2

              index_m = n_tor*(j_loc-1) + max(2*(m-1),1)

              do k=1,(n_tor+1)/2

                ik      = max(2*(k-1),1)
                index_k = n_tor*(i_loc-1) + max(2*(k-1),1)

                l = (k-1) + (m-1)

                if (l .le. n_plane/2) then
                  ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + imag(out_fft(l+1)) * float(mode(ik))
                  ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + real(out_fft(l+1)) * float(mode(ik))
                  ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + real(out_fft(l+1)) * float(mode(ik))
                  ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) - imag(out_fft(l+1)) * float(mode(ik))
                endif

                l = abs((k-1) - (m-1)); sgn_l = sign(1, (k-m))

                if (l .le. n_plane/2) then
                  ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + sgn_l * imag(out_fft(l+1)) * float(mode(ik))
                  ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   +         real(out_fft(l+1)) * float(mode(ik))
                  ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) -         real(out_fft(l+1)) * float(mode(ik))
                  ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + sgn_l * imag(out_fft(l+1)) * float(mode(ik))
                endif

              enddo

            enddo

          endif

          if (maxval(abs(ELM_kn(1:n_plane,j_loc, i_v))) .ne. 0.d0) then

            in_fft =  ELM_kn(1:n_plane,j_loc, i_v)

#ifdef USE_FFTW
            call dfftw_execute_dft_r2c(fftw_plan, in_fft, out_fft)
#else
            call my_fft(in_fft, out_fft, n_plane)
#endif

            do m=1,(n_tor+1)/2

              im      = max(2*(m-1),1)
              index_m = n_tor*(j_loc-1) + max(2*(m-1),1)

              do k=1,(n_tor+1)/2

                ik      = max(2*(k-1),1)
                index_k = n_tor*(i_loc-1) + max(2*(k-1),1)

                l = (k-1) + (m-1)

                if (l .le. n_plane/2) then
                  ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   - real(out_fft(l+1)) * float(mode(im)) * float(mode(ik))
                  ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + imag(out_fft(l+1)) * float(mode(im)) * float(mode(ik))
                  ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + imag(out_fft(l+1)) * float(mode(im)) * float(mode(ik))
                  ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + real(out_fft(l+1)) * float(mode(im)) * float(mode(ik))
                endif

                l = abs((k-1) - (m-1)); sgn_l = sign(1, (k-m))

                if (l .le. n_plane/2) then
                  ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   +         real(out_fft(l+1)) * float(mode(im)) * float(mode(ik))
                  ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   - sgn_l * imag(out_fft(l+1)) * float(mode(im)) * float(mode(ik))
                  ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + sgn_l * imag(out_fft(l+1)) * float(mode(im)) * float(mode(ik))
                  ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) +         real(out_fft(l+1)) * float(mode(im)) * float(mode(ik))
                endif

              enddo

            enddo

          endif

          if (maxval(abs(ELM_pnn(1:n_plane,j_loc, i_v))) .ne. 0.d0) then

            in_fft =  ELM_pnn(1:n_plane,j_loc, i_v)

#ifdef USE_FFTW
            call dfftw_execute_dft_r2c(fftw_plan, in_fft, out_fft)
#else
            call my_fft(in_fft, out_fft, n_plane)
#endif

            do m=1,(n_tor+1)/2

              im      = max(2*(m-1),1)
              index_m = n_tor*(j_loc-1) + max(2*(m-1),1)

              do k=1,(n_tor+1)/2

                ik      = max(2*(k-1),1)
                index_k = n_tor*(i_loc-1) + max(2*(k-1),1)

                l = (k-1) + (m-1)
                if (l .le. n_plane/2) then
                  ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   - real(out_fft(l+1)) * float(mode(im)**2)
                  ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + imag(out_fft(l+1)) * float(mode(im)**2)
                  ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + imag(out_fft(l+1)) * float(mode(im)**2)
                  ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + real(out_fft(l+1)) * float(mode(im)**2)
                endif

                l = abs((k-1) - (m-1)); sgn_l = sign(1, (k-m))

                if (l .le. n_plane/2) then
                  ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   -         real(out_fft(l+1)) * float(mode(im)**2)
                  ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + sgn_l * imag(out_fft(l+1)) * float(mode(im)**2)
                  ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) - sgn_l * imag(out_fft(l+1)) * float(mode(im)**2)
                  ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) -         real(out_fft(l+1)) * float(mode(im)**2)
                endif

              enddo

            enddo

          endif

#ifdef ENABLE_VISCO_PAR_AND_PARPAR
          if (maxval(abs(ELM_knn(1:n_plane,j_loc, i_v))) .ne. 0.d0) then

            in_fft =  ELM_knn(1:n_plane,j_loc, i_v)

#ifdef USE_FFTW
            call dfftw_execute_dft_r2c(fftw_plan, in_fft, out_fft)
#else
            call my_fft(in_fft, out_fft, n_plane)
#endif

            do m=1,(n_tor+1)/2

              im      = max(2*(m-1),1)
              index_m = n_tor*(j_loc-1) + max(2*(m-1),1)

              do k=1,(n_tor+1)/2

                ik      = max(2*(k-1),1)
                index_k = n_tor*(i_loc-1) + max(2*(k-1),1)

                l = (k-1) + (m-1)
                
                if (l .le. n_plane/2) then
                  ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   - imag(out_fft(l+1)) * float(mode(im)**2) * float(mode(ik))
                  ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   - real(out_fft(l+1)) * float(mode(im)**2) * float(mode(ik))
                  ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) - real(out_fft(l+1)) * float(mode(im)**2) * float(mode(ik))
                  ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + imag(out_fft(l+1)) * float(mode(im)**2) * float(mode(ik))
                endif

                l = abs((k-1) - (m-1)); sgn_l = sign(1, (k-m))

                if (l .le. n_plane/2) then
                  ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   - sgn_l * imag(out_fft(l+1)) * float(mode(im)**2) * float(mode(ik))
                  ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   -         real(out_fft(l+1)) * float(mode(im)**2) * float(mode(ik))
                  ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) +         real(out_fft(l+1)) * float(mode(im)**2) * float(mode(ik))
                  ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) - sgn_l * imag(out_fft(l+1)) * float(mode(im)**2) * float(mode(ik))
                endif

              enddo

            enddo

          endif
#endif

        enddo

      enddo

    endif ! apply fft (or not)

  enddo ! j loop n_degrees
enddo ! i loop (n_vertex)

if (n_tor .le. n_tor_fft_thresh) return

ELM = 0.5d0 * ELM

do j=1, n_vertex_max*n_var*n_degrees

  in_fft = RHS_p(1:n_plane,j)
#ifdef USE_FFTW
  call dfftw_execute_dft_r2c(fftw_plan, in_fft, out_fft)
#else
  call my_fft(in_fft, out_fft, n_plane)
#endif
    
  index = n_tor*(j-1) + 1
  RHS(index) = real(out_fft(1))

  do k=2,(n_tor+1)/2
    index = n_tor*(j-1) + 2*(k-1)
    RHS(index)   =   real(out_fft(k))
    RHS(index+1) = - imag(out_fft(k))
  enddo

enddo

do j=1, n_vertex_max*n_var*n_degrees

  in_fft = RHS_k(1:n_plane,j)
#ifdef USE_FFTW
  call dfftw_execute_dft_r2c(fftw_plan, in_fft, out_fft)
#else
  call my_fft(in_fft, out_fft, n_plane)
#endif
  
  index = n_tor*(j-1) + 1
  ik    = 1
  RHS(index) = RHS(index) + imag(out_fft(1)) * float(mode(ik))

  do k=2,(n_tor+1)/2
    ik    = max(2*(k-1),1)
    index = n_tor*(j-1) + 2*(k-1)
    RHS(index)   = RHS(index)   + imag(out_fft(k)) * float(mode(ik))
    RHS(index+1) = RHS(index+1) + real(out_fft(k)) * float(mode(ik))
  enddo

enddo

return


CONTAINS

! subroutine that calculates shock-capturing stabilization related terms
subroutine calculate_sc_quantities()

  d_p = 0.d0      
  ! approximate residuals
  R_rho    = rho0 * divU + UgradRho
  R_Ti     = 0.d0
  R_Te     = 0.d0
  
  ! total pressure including neutrals and impurities
  Ptot_corr= rho0_corr * (Ti0_corr + Te0_corr)
  Ptot_p   = Pi0_p + Pe0_p
  Ptot_R   = Pi0_R + Pe0_R
  Ptot_Z   = Pi0_Z + Pe0_Z

  rhoi_eff = rho0_corr 
  rhoe_eff = rho0_corr
  R_Ti = UgradTi + (gamma-1.d0) * pi0 / rhoi_eff * divU
  R_Te = UgradTe + (gamma-1.d0) * pe0 / rhoe_eff * divU

  d_p  = (Ti0 + Te0) * R_rho + rhoi_eff * R_Ti + rhoe_eff * R_Te
  
  ! Shock-detector term based on the total pressure gradient
  f_p = dsqrt( Ptot_R*Ptot_R + Ptot_Z*Ptot_Z + Ptot_p*Ptot_p/ (R*R) ) / Ptot_corr * h_e
  ! Estimation of the numerical stabilization coefficient
  tau_sc = h_e * h_e * abs(d_p) / Ptot_corr * f_p
  
  ! Use of source terms to increase the stabilization coefficients
  if(add_sources_in_sc)then
    src_rho = particle_source(ms,mt)
    src_pi  = heat_source_i(ms,mt)
    src_pe  = heat_source_e(ms,mt)
    s_p = (Ti0 + Te0) * src_rho + src_pi + src_pe
  
    tau_sc = h_e * h_e * (abs(d_p) + abs(s_p)) / Ptot_corr * f_p
  endif
  
  ! Updates in the physical diffusivities to locally add numerical
  ! stabilization.
  D_par_prof = D_par_prof + D_par_sc_num     * tau_sc ! every step D_par_prof     should be reset
  visco_T    = visco_T    + visco_sc_num     * tau_sc
  D_prof     = D_prof     + D_perp_sc_num    * tau_sc
  ZKi_par_T  = ZKi_par_T  + ZK_i_par_sc_num  * tau_sc
  ZKi_prof   = ZKi_prof   + ZK_i_perp_sc_num * tau_sc
  ZKe_par_T  = ZKe_par_T  + ZK_e_par_sc_num  * tau_sc
  ZKe_prof   = ZKe_prof   + ZK_e_perp_sc_num * tau_sc
  visco_par_prof = visco_par_prof + visco_par_sc_num * tau_sc
  ! the defination of visco_par_par_sc_num or visco_pp_diff_sc_num is strange, here used just for numerical reasons.
  visco_pp_diff_prof = visco_pp_diff_prof + visco_par_par_sc_num * tau_sc 

end subroutine calculate_sc_quantities

! subroutines to calculates VMS stabilization terms
subroutine add_vms_to_rhs()

  vms_AR__p = 0.d0 ; vms_AR__k = 0.d0
  vms_AZ__p = 0.d0 ; vms_AZ__k = 0.d0
  vms_A3__p = 0.d0 ; vms_A3__k = 0.d0
  vms_UR__p = 0.d0 ; vms_UR__k = 0.d0
  vms_UZ__p = 0.d0 ; vms_UZ__k = 0.d0
  vms_Up__p = 0.d0 ; vms_Up__k = 0.d0
  vms_rho__p= 0.d0 ; vms_rho__k= 0.d0
  vms_Ti__p  = 0.d0; vms_Ti__k  = 0.d0
  vms_Te__p  = 0.d0; vms_Te__k  = 0.d0
   
  res = 0.d0
  
  ! viscosity terms in strong form
  vsR = UR0 / R + UR0_RR + UR0_ZZ + UR0_pp / R**2 - 2.d0 * Up0_p / R - UR0 / R**2
  vsZ = UZ0 / R + UZ0_RR + UZ0_ZZ + UZ0_pp / R**2 
  vsp = Up0 / R + Up0_RR + Up0_ZZ + Up0_pp / R**2 + 2.d0 * UR0_p / R - Up0 / R**2
  
  ! residual in A
  res(var_AR) = -     (UZ0 * Bp0 - Up0 * BZ0) - tauIC_ARAZ * tau_IC * F0/rho0_corr/BB2 * BR0  * BgradPe &
                +     eta_T * (JR0 - current_source_JR(ms,mt))
  res(var_AZ) = -     (Up0 * BR0 - UR0 * Bp0) - tauIC_ARAZ * tau_IC * F0/rho0_corr/BB2 * BZ0  * BgradPe &
                +     eta_T * (JZ0 - current_source_JZ(ms,mt))
  res(var_A3) = - R * (UR0 * BZ0 - UZ0 * BR0) -              tau_IC * F0/rho0_corr/BB2 * R*Bp0* BgradPe &
                + R * eta_T * (Jp0 - current_source_Jp(ms,mt))
  
  ! residual in momentum
  res(var_UR) = rho0 * UgradUR - rho0 * Up0 * Up0/R + rho0 * VdiaGradUR - rho0 * VdiaP0*Up0 / R &
              + p0_R   - visco_T * vsR - PneoR + particle_source(ms,mt) * UR0
  res(var_UZ) = rho0 * UgradUZ                      + rho0 * VdiaGradUZ                  & 
              + p0_Z   - visco_T * vsZ - PneoZ + particle_source(ms,mt) * UZ0
  res(var_Up) = rho0 * UgradUp + rho0 * UR0 * Up0/R + rho0 * VdiaGradUp + rho0 * VdiaP0*UR0 / R &
              + p0_p/R - visco_T * vsp         + particle_source(ms,mt) * Up0
  
  ! residual in density equation
  res(var_rho)=   UgradRho + rho0 * divU - D_prof * (rho0_R / R + rho0_RR + rho0_ZZ + rho0_pp/R**2)
  
  ! residual in temperature equations
  res(var_Ti) = rho0 * UgradTi + Ti0 * UgradRho + gamma * pi0 * divU &
              - heat_source_i(ms,mt) - (gamma-1.d0) * Qvisc_T        &
              - ZKi_prof * (Ti0_R / R + Ti0_RR + Ti0_ZZ + Ti0_pp/R**2) - dTi_e
  res(var_Te) = rho0 * UgradTe + Te0 * UgradRho + gamma * pe0 * divU  &
              - heat_source_e(ms,mt) - (gamma-1.d0) * Qvisc_T       &
              - ZKe_prof * (Te0_R / R + Te0_RR + Te0_ZZ + Te0_pp/R**2) - dTe_i &
              - (gamma-1.0d0) * eta_T_ohm * JJ2  
    
  ! stabilization operators based on first derivatives only
  ! equation A_R
  vms_AR__p(var_AR) =   UZ0 * v_Z
  vms_AR__k(var_AR) =   Up0 * v_p__k / R
  
  vms_AR__p(var_AZ) = - UZ0 * v_R
  
  vms_AR__p(var_A3) = - Up0 * v_R / R
  
  ! equation A_Z
  vms_AZ__p(var_AR) = - UR0 * v_Z
  
  vms_AZ__p(var_AZ) =   UR0 * v_R
  vms_AZ__k(var_AZ) =   Up0 * v_p__k / R
  
  vms_AZ__p(var_A3) = - Up0 * v_Z / R
  
  ! equation A_3
  vms_A3__k(var_AR) = - UR0 * v_p__k
  
  vms_A3__k(var_AZ) = - UZ0 * v_p__k
  
  vms_A3__p(var_A3) = + UZ0 * v_Z + UR0 * v_R + UR0 * v / R
    
  ! equation UR
  vms_UR__p(var_UR) =  rho0 * UgradVstar__p
  vms_UR__k(var_UR) =  rho0 * UgradVstar__k
    
  vms_UR__p(var_Up) =  rho0 * (Up0 * v / R)
   
  vms_UR__p(var_rho)=  (Ti0+Te0) * (v/R + v_R) 
  vms_UR__p(var_Ti) =  rho0 * (v/R + v_R)
  vms_UR__p(var_Te) =  rho0 * (v/R + v_R)
  
  ! equation UZ
  vms_UZ__p(var_UZ) =  rho0 * UgradVstar__p
  vms_UZ__k(var_UZ) =  rho0 * UgradVstar__k
    
  vms_UZ__p(var_rho)=  (Ti0+Te0) * v_Z
  vms_UZ__p(var_Ti) =  rho0 * v_Z
  vms_UZ__p(var_Te) =  rho0 * v_Z
  
  ! equation Up
  vms_Up__p(var_UR) =  - rho0 * Up0 * v / R 
    
  vms_Up__p(var_Up) =  rho0 * UgradVstar__p
  vms_Up__k(var_Up) =  rho0 * UgradVstar__k
    
  vms_Up__k(var_rho)=  (Ti0+Te0) * v_p__k / R 
  vms_Up__k(var_Ti) =  rho0 * v_p__k / R
  vms_Up__k(var_Te) =  rho0 * v_p__k / R
  
  !! equation rho
  vms_rho__p(var_UR) =  rho0 * v_R
  
  vms_rho__p(var_UZ) =  rho0 * v_Z
    
  vms_rho__k(var_Up) = rho0 * v_p__k / R
    
  vms_rho__p(var_rho)=  UgradVstar__p
  vms_rho__k(var_rho)=  UgradVstar__k
    
  !! equation Ti
  vms_Ti__p(var_UR) = gamma * pi0 * v_R
      
  vms_Ti__p(var_UZ) = gamma * pi0 * v_Z
  
  vms_Ti__k(var_Up) = gamma * pi0 * v_p__k / R
  
  vms_Ti__p(var_rho)= Ti0 * UgradVstar__p
  vms_Ti__k(var_rho)= Ti0 * UgradVstar__k
      
  vms_Ti__p(var_Ti) = rho0 * UgradVstar__p
  vms_Ti__k(var_Ti) = rho0 * UgradVstar__k
  
  !! equation Te
  vms_Te__p(var_UR) = gamma * pe0 * v_R
    
  vms_Te__p(var_UZ) = gamma * pe0 * v_Z
  
  vms_Te__k(var_Up) = gamma * pe0 * v_p__k / R
  
  vms_Te__p(var_rho)=  Te0 * UgradVstar__p
  vms_Te__k(var_rho)=  Te0 * UgradVstar__k
    
  vms_Te__p(var_Te) =  rho0 * UgradVstar__p
  vms_Te__k(var_Te) =  rho0 * UgradVstar__k
  
#ifdef MOMENTUM_UPHI_IN_PARALLEL
  ! for B-projection
  vms_Up__p = BR0*vms_UR__p + BZ0*vms_UZ__p + Bp0*vms_Up__p
  vms_Up__k = BR0*vms_UR__k + BZ0*vms_UZ__k + Bp0*vms_Up__k
#endif
  
  Qvec_p(var_AR) = Qvec_p(var_AR) - vms_coeff_AR * tscale * dot_product(res, vms_AR__p)
  Qvec_k(var_AR) = Qvec_k(var_AR) - vms_coeff_AR * tscale * dot_product(res, vms_AR__k)
  
  Qvec_p(var_AZ) = Qvec_p(var_AZ) - vms_coeff_AZ * tscale * dot_product(res, vms_AZ__p)
  Qvec_k(var_AZ) = Qvec_k(var_AZ) - vms_coeff_AZ * tscale * dot_product(res, vms_AZ__k)
  
  Qvec_p(var_A3) = Qvec_p(var_A3) - vms_coeff_A3 * tscale * dot_product(res, vms_A3__p)
  Qvec_k(var_A3) = Qvec_k(var_A3) - vms_coeff_A3 * tscale * dot_product(res, vms_A3__k)
  
  Qvec_p(var_UR) = Qvec_p(var_UR) - vms_coeff_UR * tscale * dot_product(res, vms_UR__p)
  Qvec_k(var_UR) = Qvec_k(var_UR) - vms_coeff_UR * tscale * dot_product(res, vms_UR__k)
  
  Qvec_p(var_UZ) = Qvec_p(var_UZ) - vms_coeff_UZ * tscale * dot_product(res, vms_UZ__p)
  Qvec_k(var_UZ) = Qvec_k(var_UZ) - vms_coeff_UZ * tscale * dot_product(res, vms_UZ__k)
  
  Qvec_p(var_Up) = Qvec_p(var_Up) - vms_coeff_Up * tscale * dot_product(res, vms_Up__p)
  Qvec_k(var_Up) = Qvec_k(var_Up) - vms_coeff_Up * tscale * dot_product(res, vms_Up__k)
  
  Qvec_p(var_rho)= Qvec_p(var_rho)- vms_coeff_rho* tscale * dot_product(res, vms_rho__p)
  Qvec_k(var_rho)= Qvec_k(var_rho)- vms_coeff_rho* tscale * dot_product(res, vms_rho__k)
  
  Qvec_p(var_Ti) = Qvec_p(var_Ti) - vms_coeff_Ti  * tscale * dot_product(res, vms_Ti__p)
  Qvec_k(var_Ti) = Qvec_k(var_Ti) - vms_coeff_Ti  * tscale * dot_product(res, vms_Ti__k)
  
  Qvec_p(var_Te) = Qvec_p(var_Te) - vms_coeff_Te  * tscale * dot_product(res, vms_Te__p)
  Qvec_k(var_Te) = Qvec_k(var_Te) - vms_coeff_Te  * tscale * dot_product(res, vms_Te__k)

  
end subroutine add_vms_to_rhs

subroutine add_vms_to_elm()

  res_jac__p = 0.d0 ; res_jac__n = 0.d0 ; res_jac__nn = 0.d0
  
  vsR_UR__p = UR / R + UR_RR + UR_ZZ - UR / R**2 + UR / R**2
  vsR_UR__n = 0.d0
  vsR_UR__nn= UR_pp / R**2
  
  vsR_UZ__p = 0.d0
  vsR_UZ__n = 0.d0
  vsR_UZ__nn= 0.d0
  
  vsR_Up__p = 0.d0
  vsR_Up__n = - 2.d0 * Up_p__n / R
  vsR_Up__nn= 0.d0
  
  vsZ_UR__p = 0.d0
  vsZ_UR__n = 0.d0
  vsZ_UR__nn= 0.d0
  
  vsZ_UZ__p = UZ / R + UZ_RR + UZ_ZZ
  vsZ_UZ__n = 0.d0
  vsZ_UZ__nn= UZ_pp / R**2
  
  vsZ_Up__p = 0.d0
  vsZ_Up__n = 0.d0
  vsZ_Up__nn= 0.d0
  
  vsp_UR__p = 0.d0
  vsp_UR__n = 2.d0 * UR_p__n / R
  vsp_UR__nn= 0.d0
  
  vsp_UZ__p = 0.d0
  vsp_UZ__n = 0.d0
  vsp_UZ__nn= 0.d0
  
  vsp_Up__p = Up / R + Up_RR + Up_ZZ - Up / R**2
  vsp_Up__n = 0.d0
  vsp_Up__nn= Up_pp / R**2
  
  ! AR equation:
  res_jac__p(var_AR, var_AR) = - UZ0 * Bp0_AR &
                               - tauIC_ARAZ * tau_IC*F0/rho0_corr/BB2    * BR0 * BgradPe_AR__p &
                               + tauIC_ARAZ * tau_IC*F0/rho0_corr/BB2**2 * BR0 * BgradPe * BB2_AR__p &
                               + eta_T * JR0_AR__p
  res_jac__n(var_AR, var_AR) =   Up0 * BZ0_AR__n &
                               - tauIC_ARAZ * tau_IC*F0/rho0_corr/BB2    * BR0 * BgradPe_AR__n &
                               + tauIC_ARAZ * tau_IC*F0/rho0_corr/BB2**2 * BR0 * BgradPe * BB2_AR__n &
                               + eta_T * JR0_AR__n
  res_jac__nn(var_AR, var_AR) = + eta_T * JR0_AR__nn
  
  res_jac__p(var_AR, var_AZ) = - UZ0 * Bp0_AZ + Up0 * BZ0_AZ &
                               - tauIC_ARAZ * tau_IC*F0/rho0_corr/BB2    * BR0 * BgradPe_AZ__p &
                               + tauIC_ARAZ * tau_IC*F0/rho0_corr/BB2**2 * BR0 * BgradPe * BB2_AZ__p &
                               + eta_T * JR0_AZ__p
  res_jac__n(var_AR, var_AZ) = + tauIC_ARAZ * tau_IC*F0/rho0_corr/BB2    * BR0_AZ__n * BgradPe &
                               - tauIC_ARAZ * tau_IC*F0/rho0_corr/BB2    * BR0       * BgradPe_AZ__n &
                               + tauIC_ARAZ * tau_IC*F0/rho0_corr/BB2**2 * BR0       * BgradPe * BB2_AZ__n &
                               + eta_T * JR0_AZ__n
  res_jac__nn(var_AR, var_AZ)= + eta_T * JR0_AZ__nn
  
  res_jac__p(var_AR, var_A3) = - UZ0 * Bp0_A3 + Up0 * BZ0_A3 &
                               - tauIC_ARAZ* tau_IC*F0/rho0_corr/BB2    * BR0_A3 * BgradPe &
                               - tauIC_ARAZ* tau_IC*F0/rho0_corr/BB2    * BR0    * BgradPe_A3 &
                               + tauIC_ARAZ* tau_IC*F0/rho0_corr/BB2**2 * BR0    * BgradPe * BB2_A3 &
                               + eta_T * JR0_A3__p
  res_jac__n(var_AR, var_A3) = + eta_T * JR0_A3__n
  res_jac__nn(var_AR, var_A3)= + eta_T * JR0_A3__nn
  
  res_jac__p(var_AR, var_UZ) = - UZ * Bp0 
  res_jac__p(var_AR, var_Up) =   Up * BZ0 
  
  res_jac__p(var_AR, var_rho)= - tauIC_ARAZ * tau_IC*F0/rho0_corr   /BB2 * BR0 * BgradPe_rho__p &
                               + tauIC_ARAZ * tau_IC*F0/rho0_corr**2/BB2 * BR0 * BgradPe * rho
  
  res_jac__n(var_AR, var_rho)= + tauIC_ARAZ * tau_IC*F0/rho0_corr   /BB2 * BR0 * BgradPe_rho__n
  
  res_jac__p(var_AR,var_Te )= - tauIC_ARAZ * tau_IC*F0/rho0_corr/BB2 * BR0 * BgradPe_Te__p + eta_T_T * JR0 
  res_jac__n(var_AR,var_Te )= - tauIC_ARAZ * tau_IC*F0/rho0_corr/BB2 * BR0 * BgradPe_Te__n
  
  ! AZ equation:
  res_jac__p(var_AZ, var_AR) = - Up0 * BR0_AR + UR0 * Bp0_AR &
                               - tauIC_ARAZ * tau_IC*F0/rho0_corr/BB2    * BZ0 * BgradPe_AR__p &
                               + tauIC_ARAZ * tau_IC*F0/rho0_corr/BB2**2 * BZ0 * BgradPe * BB2_AR__p &
                               + eta_T * JZ0_AR__p
  res_jac__n(var_AZ, var_AR) = - tauIC_ARAZ * tau_IC*F0/rho0_corr/BB2    * BZ0_AR__n * BgradPe &
                               - tauIC_ARAZ * tau_IC*F0/rho0_corr/BB2    * BZ0       * BgradPe_AR__n &
                               + tauIC_ARAZ * tau_IC*F0/rho0_corr/BB2**2 * BZ0       * BgradPe * BB2_AR__n &
                               + eta_T * JZ0_AZ__n
  res_jac__nn(var_AR, var_AR) = + eta_T * JZ0_AR__nn
  
  res_jac__p(var_AZ, var_AZ) = + UR0 * Bp0_AZ &
                               - tauIC_ARAZ * tau_IC*F0/rho0_corr/BB2    * BZ0 * BgradPe_AZ__p &
                               + tauIC_ARAZ * tau_IC*F0/rho0_corr/BB2**2 * BZ0 * BgradPe * BB2_AZ__p &
                               + eta_T * JZ0_AZ__p
  res_jac__n(var_AZ, var_AZ) = - Up0 * BR0_AZ__n &
                               - tauIC_ARAZ * tau_IC*F0/rho0_corr/BB2    * BZ0 * BgradPe_AZ__n &
                               + tauIC_ARAZ * tau_IC*F0/rho0_corr/BB2**2 * BZ0 * BgradPe * BB2_AZ__n &
                               + eta_T * JZ0_AZ__n
  res_jac__nn(var_AZ, var_AZ) = + eta_T * JZ0_AZ__nn
  
  res_jac__p(var_AZ, var_A3) = - Up0 * BR0_A3 + UR0 * Bp0_A3 &
                               - tauIC_ARAZ * tau_IC*F0/rho0_corr/BB2    * BZ0_A3 * BgradPe &
                               - tauIC_ARAZ * tau_IC*F0/rho0_corr/BB2    * BZ0    * BgradPe_A3 &
                               + tauIC_ARAZ * tau_IC*F0/rho0_corr/BB2**2 * BZ0    * BgradPe * BB2_A3 &
                               + eta_T * JZ0_A3__p
  res_jac__n (var_AZ, var_A3) = + eta_T * JZ0_A3__n
  res_jac__nn(var_AZ, var_A3) = + eta_T * JZ0_A3__nn
  
  res_jac__p (var_AZ,var_UR) =   UR * Bp0
  res_jac__p (var_AZ,var_Up) = - Up * BR0
  
  res_jac__p (var_AZ,var_rho)= - tauIC_ARAZ * tau_IC*F0/rho0_corr   /BB2 * BZ0 * BgradPe_rho__p &
                               + tauIC_ARAZ * tau_IC*F0/rho0_corr**2/BB2 * BZ0 * BgradPe * rho
  res_jac__n (var_AZ,var_rho)= - tauIC_ARAZ * tau_IC*F0/rho0_corr   /BB2 * BZ0 * BgradPe_rho__n
  
  res_jac__p (var_AZ,var_Te )= - tauIC_ARAZ * tau_IC*F0/rho0_corr/BB2 * BZ0 * BgradPe_Te__p + eta_T_T * JZ0
  res_jac__n (var_AZ,var_Te )= - tauIC_ARAZ * tau_IC*F0/rho0_corr/BB2 * BZ0 * BgradPe_Te__n
  
  ! A3 equation:
  res_jac__p(var_A3, var_AR) = + UZ0 * BR0_AR * R &
                               - tau_IC*F0/rho0_corr/BB2    * R*Bp0_AR * BgradPe &
                               - tau_IC*F0/rho0_corr/BB2    * R*Bp0    * BgradPe_AR__p &
                               + tau_IC*F0/rho0_corr/BB2**2 * R*Bp0    * BgradPe * BB2_AR__p &
                               + R * eta_T * Jp0_AR__p
  res_jac__n(var_A3, var_AR) = - UR0 * BZ0_AR__n * R &
                               - tau_IC*F0/rho0_corr/BB2    * R*Bp0    * BgradPe_AR__n &
                               + tau_IC*F0/rho0_corr/BB2**2 * R*Bp0    * BgradPe * BB2_AR__n &
                               + R * eta_T * Jp0_AR__n
  res_jac__nn(var_A3, var_AR) = + R * eta_T * Jp0_AR__nn
  
  res_jac__p(var_A3, var_AZ) = - UR0 * BZ0_AZ * R &
                               - tau_IC*F0/rho0_corr/BB2    * R*Bp0_AZ * BgradPe &
                               - tau_IC*F0/rho0_corr/BB2    * R*Bp0    * BgradPe_AZ__p &
                               + tau_IC*F0/rho0_corr/BB2**2 * R*Bp0    * BgradPe * BB2_AZ__p &
                               + R * eta_T * Jp0_AZ__p
  res_jac__n(var_A3, var_AZ) = + UZ0 * BR0_AZ__n * R &
                               - tau_IC*F0/rho0_corr/BB2    * R*Bp0    * BgradPe_AZ__n &
                               + tau_IC*F0/rho0_corr/BB2**2 * R*Bp0    * BgradPe * BB2_AZ__n &
                               + R * eta_T * Jp0_AZ__n
  res_jac__nn(var_A3, var_AZ) = + R * eta_T * Jp0_AZ__nn
  
  res_jac__p(var_A3, var_A3) = - UR0 * BZ0_A3 * R + UZ0 * BR0_A3 * R &
                               - tau_IC*F0/rho0_corr/BB2    * R*Bp0    * BgradPe_A3 &
                               + tau_IC*F0/rho0_corr/BB2**2 * R*Bp0    * BgradPe * BB2_A3 &
                               + R * eta_T * Jp0_A3__p
  res_jac__n (var_A3, var_A3) = + R * eta_T * Jp0_A3__n
  res_jac__nn(var_A3, var_A3) = + R * eta_T * Jp0_A3__nn
  
  
  res_jac__p(var_A3,var_UR) = - UR * BZ0
  res_jac__n(var_A3,var_UZ) = + UZ * BR0
  
  res_jac__p(var_A3,var_rho) = - tau_IC*F0/rho0_corr   /BB2 * R*Bp0 * BgradPe_rho__p &
                               + tau_IC*F0/rho0_corr**2/BB2 * R*Bp0 * BgradPe * rho
  res_jac__n (var_A3,var_rho)= - tau_IC*F0/rho0_corr/BB2 * R*Bp0 * BgradPe_rho__n
  
  res_jac__p (var_A3,var_Te )= - tau_IC*F0/rho0_corr/BB2 * R*Bp0 * BgradPe_Te__p + R * eta_T_T * Jp0
  res_jac__n (var_A3,var_Te )= - tau_IC*F0/rho0_corr/BB2 * R*Bp0 * BgradPe_Te__n
  
  ! UR equation:
  res_jac__p(var_UR,var_AR ) =  rho0 * ( VdiaGradUR_AR__p  - VdiaP0_AR__p * Up0 / R ) - PneoR_AR__p
  res_jac__n(var_UR,var_AR ) =  rho0 * ( VdiaGradUR_AR__n  - VdiaP0_AR__n * Up0 / R ) - PneoR_AR__n
  
  res_jac__p(var_UR,var_AZ ) =  rho0 * ( VdiaGradUR_AZ__p  - VdiaP0_AZ__p * Up0 / R ) - PneoR_AZ__p
  res_jac__n(var_UR,var_AZ ) =  rho0 * ( VdiaGradUR_AZ__n  - VdiaP0_AZ__n * Up0 / R ) - PneoR_AZ__n
  
  res_jac__p(var_UR,var_A3 ) =  rho0 * ( VdiaGradUR_A3__p  - VdiaP0_A3__p * Up0 / R ) - PneoR_A3__p
  res_jac__n(var_UR,var_A3 ) =  rho0 * ( VdiaGradUR_A3__n  - VdiaP0_A3__n * Up0 / R ) - PneoR_A3__n
  
  res_jac__p(var_UR, var_UR) =  rho0 * UgradUR_UR__p + rho0 * VdiaGradUR_UR__p &
                               - visco_T * vsR_UR__p - PneoR_UR + particle_source(ms,mt) * UR 
  res_jac__n(var_UR, var_UR) =  rho0 * UgradUR_UR__n + rho0 * VdiaGradUR_UR__n &
                               - visco_T * vsR_UR__n
  res_jac__nn(var_UR, var_UR)= - visco_T * vsR_UR__nn
  
  res_jac__p(var_UR, var_UZ) =  rho0 * UgradUR_UZ - visco_T * vsR_UZ__p - PneoR_UZ
  res_jac__n(var_UR, var_UZ) =                    - visco_T * vsR_UZ__n           
  res_jac__nn(var_UR, var_UZ)=                    - visco_T * vsR_UZ__nn          
  
  res_jac__p(var_UR, var_Up) = rho0 * UgradUR_Up - 2.0*rho0*Up*Up0**2 / R - rho0 * VdiaP0 * Up / R &
                                             - visco_T * vsR_Up__p 
  res_jac__n(var_UR, var_Up) =               - visco_T * vsR_Up__n 
  res_jac__nn(var_UR, var_Up)=               - visco_T * vsR_Up__nn
  
  ! approximation
  res_jac__p (var_UR,var_rho) = rho * UgradUR - rho * Up0 * Up0 / R                     &
                              + rho0 * VdiaGradUR_rho__p + rho * VdiaGradUR             &
                              - rho0 * VdiaP0_rho__p * Up0 / R - rho * VdiaP0 * Up0 / R &
                              + rho_R*(Ti0+Te0) + rho*(Ti0_R+Te0_R) - PneoR_rho__p
  res_jac__n (var_UR,var_rho) = rho0 * VdiaGradUR_rho__n                                &
                              - rho0 * VdiaP0_rho__n * Up0 / R                          &
                                                                   - PneoR_rho__n   
  
  res_jac__p (var_UR,var_Ti ) = rho0 * VdiaGradUR_Ti__p  - rho0 * VdiaP0_Ti__p * Up0 / R  &
                              + (rho0_R*Ti     + rho0*Ti_R       ) - PneoR_Ti__p !- visco_T_T * Qvisc_UR__p ! this term should not be included ?
  res_jac__n (var_UR,var_Ti ) = rho0 * VdiaGradUR_Ti__n  - rho0 * VdiaP0_Ti__n * Up0 / R &
                                                                    - PneoR_Ti__n 
  
  res_jac__p (var_UR,var_Te ) =  (rho0_R*Te     + rho0*Te_R       ) + visco_T_T * vsR 
  
  ! UZ equation:
  res_jac__p(var_UZ,var_AR ) =  rho0 * VdiaGradUZ_AR__p - PneoZ_AR__p
  res_jac__n(var_UZ,var_AR ) =  rho0 * VdiaGradUZ_AR__n - PneoZ_AR__n
  
  res_jac__p(var_UZ,var_AZ ) =  rho0 * VdiaGradUZ_AZ__p - PneoZ_AZ__p
  res_jac__n(var_UZ,var_AZ ) =  rho0 * VdiaGradUZ_AZ__n - PneoZ_AZ__n
  
  res_jac__p(var_UZ,var_A3 ) =  rho0 * VdiaGradUZ_A3__p - PneoZ_A3__p
  res_jac__n(var_UZ,var_A3 ) =  rho0 * VdiaGradUZ_A3__n - PneoZ_A3__n
  
  res_jac__p(var_UZ, var_UR) =  rho0 * UgradUZ_UR - visco_T * vsZ_UR__p - PneoZ_UR
  res_jac__n(var_UZ, var_UR) =                    - visco_T * vsR_UR__n           
  res_jac__nn(var_UZ, var_UR)=                    - visco_T * vsR_UR__nn          
  
  res_jac__p(var_UZ, var_UZ) = rho0 * UgradUZ_UZ__p - visco_T * vsZ_UZ__p - PneoZ_UZ + particle_source(ms,mt) * UZ
  res_jac__n(var_UZ, var_UZ) = rho0 * UgradUZ_UZ__n - visco_T * vsZ_UZ__n
  res_jac__nn(var_UZ, var_UZ)=                      - visco_T * vsZ_UZ__nn
  
  res_jac__p(var_UZ, var_Up) = rho0 * UgradUZ_Up - visco_T * vsZ_Up__p 
  res_jac__n(var_UZ, var_Up) =                   - visco_T * vsZ_Up__n 
  res_jac__nn(var_UZ, var_Up)=                   - visco_T * vsZ_Up__nn
  
  ! approximation
  res_jac__p (var_UZ,var_rho) =  rho0 * VdiaGradUZ_rho__p + rho_Z*(Ti0+Te0) + rho*(Ti0_Z+Te0_Z) - PneoZ_rho__p
  res_jac__n (var_UZ,var_rho) =  rho0 * VdiaGradUZ_rho__n                                        - PneoZ_rho__n
  
  res_jac__p (var_UZ,var_Ti ) =  rho0 * VdiaGradUZ_Ti__p  + (rho0_Z*Ti       + rho0*Ti_Z       ) - PneoZ_Ti__p
  res_jac__n (var_UZ,var_Ti ) =  rho0 * VdiaGradUZ_Ti__n                                         - PneoZ_Ti__n
  
  res_jac__p (var_UZ,var_Te ) =  (rho0_Z*Te     + rho0*Te_Z       ) - visco_T_T * vsZ
  
  ! Up equation:
  res_jac__p(var_Up,var_AR ) =  rho0 * (VdiaGradUp_AR__p  + VdiaP0_AR__p * UR0 / R) - PneoR_AR__p
  res_jac__n(var_Up,var_AR ) =  rho0 * (VdiaGradUp_AR__n  + VdiaP0_AR__n * UR0 / R) - PneoR_AR__n
  
  res_jac__p(var_Up,var_AZ ) =  rho0 * (VdiaGradUp_AZ__p  + VdiaP0_AZ__p * UR0 / R) - PneoR_AZ__p
  res_jac__n(var_Up,var_AZ ) =  rho0 * (VdiaGradUp_AZ__n  + VdiaP0_AZ__n * UR0 / R) - PneoR_AZ__n
  
  res_jac__p(var_Up,var_A3 ) =  rho0 * (VdiaGradUp_A3__p  + VdiaP0_A3__p * UR0 / R) - PneoR_A3__p
  res_jac__n(var_Up,var_A3 ) =  rho0 * (VdiaGradUp_A3__n  + VdiaP0_A3__n * UR0 / R) - PneoR_A3__n
  
  res_jac__p(var_Up, var_UR) =  rho0 * UgradUp_UR     + rho0 * UR * Up0 / R      + rho0 * VdiaP0 * UR / R &
                               - visco_T * vsp_UR__p
  res_jac__n(var_Up, var_UR) = - visco_T * vsp_UR__n
  res_jac__nn(var_Up, var_UR)= - visco_T * vsp_UR__nn
  
  res_jac__p(var_Up, var_UZ) =  rho0 * UgradUp_UZ - visco_T * vsp_UZ__p 
  res_jac__n(var_Up, var_UZ) =                    - visco_T * vsp_UZ__n 
  res_jac__nn(var_Up, var_UZ)=                    - visco_T * vsp_UZ__nn
  
  res_jac__p(var_Up, var_Up) =  rho0 * UgradUp_Up__p + rho0 * UR0 * Up / R - rho0 * VdiaGradUp_Up__p &
                                             - visco_T * vsR_Up__p + particle_source(ms,mt) * Up
  res_jac__n(var_Up, var_Up) =  rho0 * UgradUp_Up__n                       - rho0 * VdiaGradUp_Up__n &
                                             - visco_T * vsR_Up__n 
  res_jac__nn(var_Up, var_Up)=               - visco_T * vsR_Up__nn
  
  res_jac__p (var_Up,var_rho) = rho0 *  ( VdiaGradUp_rho__p + VdiaP0_rho__p * UR0 / R ) &
                              + rho  *  ( VdiaGradUp        + VdiaP0 * UR0 / R )
  res_jac__n (var_Up,var_rho) = rho0 *  ( VdiaGradUp_rho__n + VdiaP0_rho__n * UR0 / R ) &
                              + (rho_p__n*(Ti0+Te0) + rho*(Ti0_p+Te0_p))/R
  
  res_jac__p (var_Up,var_Ti ) = rho0 * ( VdiaGradUp_Ti__p  + VdiaP0_Ti__p * UR0 / R ) 
  res_jac__n (var_Up,var_Ti ) = rho0 * ( VdiaGradUp_Ti__n  + VdiaP0_Ti__n * UR0 / R ) &
                              +  (rho0_p * Ti   + rho0 * Ti_p__n   ) / R
  
  res_jac__p (var_Up,var_Te ) =  - visco_T_T * vsp
  res_jac__n (var_Up,var_Te ) =  (rho0_p * Te   + rho0 * Te_p__n   ) / R - visco_T_T * vsp 
  
  ! rho equation
  res_jac__p(var_rho, var_UR) =  UR * rho0_R + rho0 * divU_UR
  
  res_jac__p(var_rho, var_UZ) =  UZ * rho0_Z + rho0 * divU_UZ
  
  res_jac__p(var_rho, var_Up) =  Up * rho0_p / R
  res_jac__n(var_rho, var_Up) =  rho0 * divU_Up__n
  
  res_jac__p(var_rho, var_rho) = UR0 * rho_R + UZ0 * rho_Z + rho * divU &
                               - D_prof * (rho_R / R + rho_RR + rho_ZZ)
  res_jac__n(var_rho, var_rho) = Up0 * rho_p__n / R
  res_jac__nn(var_rho, var_rho) = - D_prof * rho_pp / R**2
  
  ! Ti equation
  res_jac__p(var_Ti, var_UR) =  rho0 * UR * Ti0_R + Ti0 * UR * rho0_R + gamma * p0 * divU_UR

  res_jac__p(var_Ti, var_UZ) =  rho0 * UZ * Ti0_Z + Ti0 * UZ * rho0_Z + gamma * p0 * divU_UZ
  
  res_jac__p(var_Ti, var_Up) =  rho0 * Up * Ti0_p / R + Ti0 * Up * rho0_p / R
  res_jac__n(var_Ti, var_Up) =  gamma * p0 * divU_Up__n
  
  res_jac__p(var_Ti, var_rho)=  rho * UgradTi + Ti0 * (UR0*rho_R + UZ0*rho_Z) + gamma * rho * Ti0 * divU  
  res_jac__n(var_Ti, var_rho)=                  Ti0 * (UZ0*rho_p__n / R)
  
  res_jac__p(var_Ti, var_Ti) =   rho0 * (UR0 * Ti_R + UZ0 * Ti_Z)                      &
                               + Ti * (UR0 * rho0_R + UZ0 * rho0_Z + Up0 * rho0_p / R) &
                               + gamma * rho0 * Ti * divU                              &
                               - ZKi_prof * (Ti_R / R + Ti_RR + Ti_ZZ) 
  res_jac__n(var_Ti, var_Ti) =   rho0 * Up0 * Ti_p__n / R
  res_jac__nn(var_Ti, var_Ti)= - ZKi_prof * Ti_pp / R**2 
  
  res_jac__p(var_Ti, var_Te) =   gamma * rho0 * Te * divU
  
  ! Te equation
  res_jac__p(var_Te, var_UR) =  rho0 * UR * Te0_R + Te0 * UR * rho0_R + gamma * p0 * divU_UR
  
  res_jac__p(var_Te, var_UZ) =  rho0 * UZ * Te0_Z + Te0 * UZ * rho0_Z + gamma * p0 * divU_UZ
  
  res_jac__p(var_Te, var_Up) =  rho0 * Up * Te0_p / R + Te0 * Up * rho0_p / R
  res_jac__n(var_Te, var_Up) =  gamma * p0 * divU_Up__n
  
  res_jac__p(var_Te, var_rho)=  rho * UgradTe + Te0 * (UR0*rho_R + UZ0*rho_Z) + gamma * rho * Te0 * divU
  res_jac__n(var_Te, var_rho)=                  Te0 * (UZ0*rho_p__n / R)
  
  res_jac__p(var_Te, var_Ti) =  gamma * rho0 * Ti * divU
  
  res_jac__p(var_Te, var_Te) =   rho0 * ( UR0 * Te_R + UZ0 * Te_Z )                      &
                               + Te   * (UR0 * rho0_R + UZ0 * rho0_Z + Up0 * rho0_p / R) &
                               + gamma * rho0 * Te * divU                                &
                               - ZKe_prof * (Te_R / R + Te_RR + Te_ZZ )
  res_jac__n(var_Te, var_Te) =   rho0 * Up0 * Te_p__n / R
  res_jac__nn(var_Te, var_Te)= - ZKe_prof * Te_pp / R**2

  do jj =1, n_var
    Qjac_p (var_AR,jj)  =  Qjac_p (var_AR,jj) - vms_coeff_AR * tscale * dot_product(vms_AR__p(:) , res_jac__p(:, jj))
    Qjac_k (var_AR,jj)  =  Qjac_k (var_AR,jj) - vms_coeff_AR * tscale * dot_product(vms_AR__k(:) , res_jac__p(:, jj))
    Qjac_n (var_AR,jj)  =  Qjac_n (var_AR,jj) - vms_coeff_AR * tscale * dot_product(vms_AR__p(:) , res_jac__n(:, jj))
    Qjac_kn(var_AR,jj)  =  Qjac_kn(var_AR,jj) - vms_coeff_AR * tscale * dot_product(vms_AR__k(:) , res_jac__n(:, jj))
    Qjac_pnn(var_AR,jj) =  Qjac_pnn(var_AR,jj) - vms_coeff_AR * tscale * dot_product(vms_AR__p(:) , res_jac__nn(:, jj))
    
    Qjac_p (var_AZ,jj)  =  Qjac_p (var_AZ,jj) - vms_coeff_AZ * tscale * dot_product(vms_AZ__p(:) , res_jac__p(:, jj))
    Qjac_k (var_AZ,jj)  =  Qjac_k (var_AZ,jj) - vms_coeff_AZ * tscale * dot_product(vms_AZ__k(:) , res_jac__p(:, jj))
    Qjac_n (var_AZ,jj)  =  Qjac_n (var_AZ,jj) - vms_coeff_AZ * tscale * dot_product(vms_AZ__p(:) , res_jac__n(:, jj))
    Qjac_kn(var_AZ,jj)  =  Qjac_kn(var_AZ,jj) - vms_coeff_AZ * tscale * dot_product(vms_AZ__k(:) , res_jac__n(:, jj))
    Qjac_pnn(var_AZ,jj) =  Qjac_pnn(var_AZ,jj)- vms_coeff_AZ * tscale * dot_product(vms_AZ__p(:) , res_jac__nn(:, jj))
  
    Qjac_p (var_A3,jj)  =  Qjac_p (var_A3,jj) - vms_coeff_A3 * tscale * dot_product(vms_A3__p(:) , res_jac__p(:, jj))
    Qjac_k (var_A3,jj)  =  Qjac_k (var_A3,jj) - vms_coeff_A3 * tscale * dot_product(vms_A3__k(:) , res_jac__p(:, jj))
    Qjac_n (var_A3,jj)  =  Qjac_n (var_A3,jj) - vms_coeff_A3 * tscale * dot_product(vms_A3__p(:) , res_jac__n(:, jj))
    Qjac_kn(var_A3,jj)  =  Qjac_kn(var_A3,jj) - vms_coeff_A3 * tscale * dot_product(vms_A3__k(:) , res_jac__n(:, jj))
    Qjac_pnn(var_A3,jj) =  Qjac_pnn(var_A3,jj)- vms_coeff_A3 * tscale * dot_product(vms_A3__p(:) , res_jac__nn(:, jj))
    
    Qjac_p (var_UR,jj)  =  Qjac_p (var_UR,jj) - vms_coeff_UR * tscale * dot_product(vms_UR__p(:) , res_jac__p(:, jj))
    Qjac_k (var_UR,jj)  =  Qjac_k (var_UR,jj) - vms_coeff_UR * tscale * dot_product(vms_UR__k(:) , res_jac__p(:, jj))
    Qjac_n (var_UR,jj)  =  Qjac_n (var_UR,jj) - vms_coeff_UR * tscale * dot_product(vms_UR__p(:) , res_jac__n(:, jj))
    Qjac_kn(var_UR,jj)  =  Qjac_kn(var_UR,jj) - vms_coeff_UR * tscale * dot_product(vms_UR__k(:) , res_jac__n(:, jj))
    Qjac_pnn(var_UR,jj) =  Qjac_pnn(var_UR,jj)- vms_coeff_UR * tscale * dot_product(vms_UR__p(:) , res_jac__nn(:, jj))
  
    Qjac_p (var_UZ,jj)  =  Qjac_p (var_UZ,jj) - vms_coeff_UZ * tscale * dot_product(vms_UZ__p(:) , res_jac__p(:, jj))
    Qjac_k (var_UZ,jj)  =  Qjac_k (var_UZ,jj) - vms_coeff_UZ * tscale * dot_product(vms_UZ__k(:) , res_jac__p(:, jj))
    Qjac_n (var_UZ,jj)  =  Qjac_n (var_UZ,jj) - vms_coeff_UZ * tscale * dot_product(vms_UZ__p(:) , res_jac__n(:, jj))
    Qjac_kn(var_UZ,jj)  =  Qjac_kn(var_UZ,jj) - vms_coeff_UZ * tscale * dot_product(vms_UZ__k(:) , res_jac__n(:, jj))
    Qjac_pnn(var_UZ,jj) =  Qjac_pnn(var_UZ,jj)- vms_coeff_UZ * tscale * dot_product(vms_UZ__p(:) , res_jac__nn(:, jj))
    
    Qjac_p (var_Up,jj)  =  Qjac_p (var_Up,jj) - vms_coeff_Up * tscale * dot_product(vms_Up__p(:) , res_jac__p(:, jj))
    Qjac_k (var_Up,jj)  =  Qjac_k (var_Up,jj) - vms_coeff_Up * tscale * dot_product(vms_Up__k(:) , res_jac__p(:, jj))
    Qjac_n (var_Up,jj)  =  Qjac_n (var_Up,jj) - vms_coeff_Up * tscale * dot_product(vms_Up__p(:) , res_jac__n(:, jj))
    Qjac_kn(var_Up,jj)  =  Qjac_kn(var_Up,jj) - vms_coeff_Up * tscale * dot_product(vms_Up__k(:) , res_jac__n(:, jj))
    Qjac_pnn(var_Up,jj) =  Qjac_pnn(var_Up,jj)- vms_coeff_Up * tscale * dot_product(vms_Up__p(:) , res_jac__nn(:, jj))
    
    Qjac_p (var_rho,jj) =  Qjac_p (var_rho,jj)- vms_coeff_rho* tscale * dot_product(vms_rho__p(:), res_jac__p(:, jj))
    Qjac_k (var_rho,jj) =  Qjac_k (var_rho,jj)- vms_coeff_rho* tscale * dot_product(vms_rho__k(:), res_jac__p(:, jj))
    Qjac_n (var_rho,jj) =  Qjac_n (var_rho,jj)- vms_coeff_rho* tscale * dot_product(vms_rho__p(:), res_jac__n(:, jj))
    Qjac_kn(var_rho,jj) =  Qjac_kn(var_rho,jj)- vms_coeff_rho* tscale * dot_product(vms_rho__k(:), res_jac__n(:, jj))
    Qjac_pnn(var_rho,jj)=  Qjac_pnn(var_rho,jj)- vms_coeff_rho* tscale *dot_product(vms_rho__p(:), res_jac__nn(:, jj))
    

    Qjac_p (var_Ti,jj)   =  Qjac_p (var_Ti,jj)  - vms_coeff_Ti  * tscale * dot_product(vms_Ti__p(:)  , res_jac__p(:, jj))
    Qjac_k (var_Ti,jj)   =  Qjac_k (var_Ti,jj)  - vms_coeff_Ti  * tscale * dot_product(vms_Ti__k(:)  , res_jac__p(:, jj))
    Qjac_n (var_Ti,jj)   =  Qjac_n (var_Ti,jj)  - vms_coeff_Ti  * tscale * dot_product(vms_Ti__p(:)  , res_jac__n(:, jj))
    Qjac_kn(var_Ti,jj)   =  Qjac_kn(var_Ti,jj)  - vms_coeff_Ti  * tscale * dot_product(vms_Ti__k(:)  , res_jac__n(:, jj))
    Qjac_pnn(var_Ti,jj)  =  Qjac_pnn(var_Ti,jj)  - vms_coeff_Ti  * tscale * dot_product(vms_Ti__p(:) , res_jac__nn(:, jj))
  
    Qjac_p (var_Te,jj)   =  Qjac_p (var_Te,jj)  - vms_coeff_Te  * tscale * dot_product(vms_Te__p(:)  , res_jac__p(:, jj))
    Qjac_k (var_Te,jj)   =  Qjac_k (var_Te,jj)  - vms_coeff_Te  * tscale * dot_product(vms_Te__k(:)  , res_jac__p(:, jj))
    Qjac_n (var_Te,jj)   =  Qjac_n (var_Te,jj)  - vms_coeff_Te  * tscale * dot_product(vms_Te__p(:)  , res_jac__n(:, jj))
    Qjac_kn(var_Te,jj)   =  Qjac_kn(var_Te,jj)  - vms_coeff_Te  * tscale * dot_product(vms_Te__k(:)  , res_jac__n(:, jj))
    Qjac_pnn(var_Te,jj)  =  Qjac_pnn(var_Te,jj)  - vms_coeff_Te  * tscale * dot_product(vms_Te__p(:) , res_jac__nn(:, jj))
  
  enddo
       
  end subroutine add_vms_to_elm

subroutine my_fft(in_fft,out_fft,n)
  implicit none
  real*8     :: in_fft(*)
  complex*16 :: out_fft(*)
  integer    :: i, n
  real*8     :: tmp_fft(2*n+2)
  tmp_fft(1:n) = in_fft(1:n)
  call RFT2(tmp_fft,n,1)
  do i=1,n
    out_fft(i) = cmplx(tmp_fft(2*i-1),tmp_fft(2*i))
  enddo
  return
end subroutine my_fft


end subroutine element_matrix_fft

end module mod_elt_matrix_fft
