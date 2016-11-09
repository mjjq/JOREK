module mod_elt_matrix_fft
  implicit none
contains

subroutine element_matrix_fft(element, nodes, xpoint2, xcase2, minRad, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, ELM, RHS, tid)
! NOT YET IMPLEMENTED

use mod_parameters
use data_structure
use gauss
use basis_at_gaussian
use phys_module

implicit none

include 'mpif.h'

type (type_element)   :: element
type (type_node)      :: nodes(n_vertex_max)

real*8, dimension (:,:), pointer  :: ELM
real*8, dimension (:)  , pointer  :: RHS
integer, intent(in) :: tid

integer    :: i, j, k, l, index_ij, index_kl, index, xcase2
integer    :: in, im, ij, kl, ivar, kvar, ms, mt, mp
real*8     :: wst, xjac, xjac_s, xjac_t, xjac_R, xjac_Z, xjac3, BigR, phi
real*8     :: current_source(n_gauss,n_gauss),particle_source(n_gauss,n_gauss),heat_source(n_gauss,n_gauss), source_pellet
real*8     :: minRad, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint(2), Z_xpoint(2), dj_dpsi, dj_dz

real*8     :: uR0, uR0_R, uR0_Z, uR0_p, uR0_s, uR0_t, uR0_ss, uR0_st, uR0_tt
real*8     :: uZ0, uZ0_R, uZ0_Z, uZ0_p, uZ0_s, uZ0_t, uZ0_ss, uZ0_st, uZ0_tt
real*8     :: up0, up0_R, up0_Z, up0_p, up0_s, up0_t, up0_ss, up0_st, up0_tt
real*8     :: AR0, AR0_R, AR0_Z, AR0_p, AR0_s, AR0_t, AR0_ss, AR0_st, AR0_tt
real*8     :: AZ0, AZ0_R, AZ0_Z, AZ0_p, AZ0_s, AZ0_t, AZ0_ss, AZ0_st, AZ0_tt
real*8     :: A30, A30_R, A30_Z, A30_p, A30_s, A30_t, A30_ss, A30_st, A30_tt
real*8     :: r0, r0_R, r0_Z, r0_p, r0_s, r0_t, r0_ss, r0_st, r0_tt
real*8     :: T0, T0_R, T0_Z, T0_p, T0_s, T0_t, T0_ss, T0_st, T0_tt
real*8     :: p0, p0_R, p0_Z, p0_p, p0_s, p0_t, p0_ss, p0_st, p0_tt
real*8     :: uR, uR_R, uR_Z, uR_p, uR_s, uR_t, uR_ss, uR_st, uR_tt
real*8     :: uZ, uZ_R, uZ_Z, uZ_p, uZ_s, uZ_t, uZ_ss, uZ_st, uZ_tt
real*8     :: up, up_R, up_Z, up_p, up_s, up_t, up_ss, up_st, up_tt
real*8     :: AR, AR_R, AR_Z, AR_p, AR_s, AR_t, AR_ss, AR_st, AR_tt
real*8     :: AZ, AZ_R, AZ_Z, AZ_p, AZ_s, AZ_t, AZ_ss, AZ_st, AZ_tt
real*8     :: A3, A3_R, A3_Z, A3_p, A3_s, A3_t, A3_ss, A3_st, A3_tt
real*8     :: T, T_R, T_Z, T_p, T_s, T_t, T_ss, T_st, T_tt
real*8     :: r, r_R, r_Z, r_p, r_s, r_t, r_ss, r_st, r_tt

real*8     :: v, v_R, v_Z, v_s, v_t, v_p, v_ss, v_st, v_tt, v_RR, v_ZZ, v_RZ
real*8     :: bf, bf_R, bf_Z, bf_s, bf_t, bf_p, bf_ss, bf_st, bf_tt, bf_RR, bf_ZZ, bf_RZ
   
real*8     :: B0grad_T0,    B0grad_T0_AR,    B0grad_T0_AZ,    B0grad_T0_A3
real*8     :: B0grad_r0,    B0grad_r0_AR,    B0grad_r0_AZ,    B0grad_r0_A3
real*8     :: BB2,          BB2_AR,          BB2_AZ,          BB2_A3
real*8     :: BR0,          BR0_AR,          BR0_AZ,          BR0_A3 
real*8     :: BZ0,          BZ0_AR,          BZ0_AZ,          BZ0_A3 
real*8     :: Bp0,          Bp0_AR,          Bp0_AZ,          Bp0_A3 
real*8     :: B0grad_vstar,  B0grad_vstar_AR,  B0grad_vstar_AZ,  B0grad_vstar_A3
real*8     :: u0grad_vstar, u0grad_T0,    u0grad_r0,       gradr0grad_vstar, gradT0grad_vstar, gradbfgrad_vstar, u0grad_bf, B0grad_bf
real*8     :: divu,         divu_uR,         divu_uZ,          divu_up

real*8     :: ZK_prof, D_prof, psi_norm, theta, zeta

real*8     :: eta_T, visco_T, deta_dT, d2eta_d2T, dvisco_dT, visco_num_T, eta_num_T, eta_R, eta_Z, eta_p

real*8     :: Qvisc_uR, Qvisc_uZ, Qvisc_up, Qvisc_T

real*8     :: Fprofile(n_gauss,n_gauss), Fprof

real*8     :: amat(n_var,n_var), Pjac(n_var,n_var), Qjac(n_var,n_var), rhs_ij(n_var), Pvec_prev(n_var), Qvec(n_var)

logical    :: xpoint2

real*8, dimension(n_gauss,n_gauss)    :: x_g, x_s, x_t
real*8, dimension(n_gauss,n_gauss)    :: x_ss, x_st, x_tt
real*8, dimension(n_gauss,n_gauss)    :: y_g, y_s, y_t
real*8, dimension(n_gauss,n_gauss)    :: y_ss, y_st, y_tt

real*8, dimension(:,:,:,:) , pointer :: eq_g, eq_s, eq_t
real*8, dimension(:,:,:,:) , pointer :: eq_p
real*8, dimension(:,:,:,:) , pointer :: eq_ss, eq_st, eq_tt   
real*8, dimension(:,:,:,:) , pointer :: delta_g, delta_s, delta_t


integer*4  :: rank
integer    :: my_id, ierr

return
end subroutine element_matrix_fft
end module mod_elt_matrix_fft
