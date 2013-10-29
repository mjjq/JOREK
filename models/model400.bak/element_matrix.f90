module mod_elt_matrix
contains





subroutine element_matrix_fft(element, nodes, xpoint2, xcase2, minRad, R_axis, Z_axis, psi_axis, psi_bnd, Z_xpoint, ELM, RHS, tid)
!--------------------------------------------------------------------------
! This is just a wrapper to the real routine since I combined both into one
!--------------------------------------------------------------------------

  use data_structure

  implicit none

  type (type_element)               :: element
  type (type_node)                  :: nodes(n_vertex_max)

  integer    :: xcase2
  logical    :: xpoint2
  real*8     :: minRad, R_axis, Z_axis, psi_axis, psi_bnd, Z_xpoint(2)
  real*8, dimension (:,:), pointer  :: ELM
  real*8, dimension (:)  , pointer  :: RHS
  integer, intent(in)               :: tid

  call element_matrix(element, nodes, xpoint2, xcase2, minRad, R_axis, Z_axis, psi_axis, psi_bnd, Z_xpoint, ELM, RHS, tid)

  return

end subroutine element_matrix_fft






subroutine element_matrix(element, nodes, xpoint2, xcase2, minRad, R_axis, Z_axis, psi_axis, psi_bnd, Z_xpoint, ELM, RHS, tid)
!---------------------------------------------------------------
! calculates the matrix contribution of one element
!---------------------------------------------------------------

use constants
use parameters
use data_structure
use gauss
use basis_at_gaussian
use phys_module
use tr_module 
use profiles, only: interpolProf

implicit none

type (type_element)                         :: element
type (type_node)                            :: nodes(n_vertex_max)
type (type_surface_list)                    :: flux_list

integer, intent(in)                         :: tid
real*8, dimension (:,:)  , pointer          :: ELM
real*8, dimension (:)    , pointer          :: RHS
real*8, dimension(:,:,:) , pointer          :: ELM_p, ELM_n, ELM_k, ELM_kn
real*8, dimension(:,:)   , pointer          :: RHS_p, RHS_k 
real*8, dimension(n_tor,n_plane)            :: HHZ, HHZ_p
real*8                                      :: in_fft(1:n_plane)
complex*16                                  :: out_fft(1:n_plane)

integer    :: i,  j,  k,  l
integer    :: i2, j2, k2, l2, ms, mt, mp, m, in, im, ik, id
integer    :: ij1, ij2, ij3, ij4, ij5, ij6, ij7, ij8, ij9, kl1, kl2, kl3, kl4, kl5, kl6, kl7, kl8, kl9
integer    :: n_tor_loop, n_tor_loop2
integer    :: index, index_k, index_m, index_ij, index_kl
integer    :: i_ij, i_kl, ij(n_var), kl(n_var)

logical    :: xpoint2
integer    :: xcase2
real*8     :: minRad, R_axis, Z_axis, psi_axis, psi_bnd, Z_xpoint(2)
real*8     :: wst

real*8     :: zn, dn_dpsi, dn_dz, dn_dpsi2, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz
real*8     :: zTi,dTi_dpsi,dTi_dz,dTi_dpsi2,dTi_dz2,dTi_dpsi_dz,dTi_dpsi3,dTi_dpsi_dz2,dTi_dpsi2_dz
real*8     :: zTe,dTe_dpsi,dTe_dz,dTe_dpsi2,dTe_dz2,dTe_dpsi_dz,dTe_dpsi3,dTe_dpsi_dz2,dTe_dpsi2_dz
real*8     :: current_source, particle_source, heat_source_i, heat_source_e

real*8     :: eta_Te,	visco_Te,   deta_dTe, d2eta_d2Te, dvisco_dTe
real*8     :: eta_numm, visco_numm, visco_par_numm, D_perp_numm, K_perp_numm
real*8     :: ZK_i_prof, ZK_i_par, dZK_i_par
real*8     :: ZK_e_prof, ZK_e_par, dZK_e_par, D_prof 
real*8     :: atn_D, datn_D, atn_D_n, pol_D, dpol_D, D_min, psi_D
real*8     :: prof(1:3),Diff(1:3,1:10)
real*8     :: Q, dQ_drho, dQ_dTe

real*8     :: eps_cyl, theta, zeta
real*8     :: BB2, BB2_psi, psi_norm
real*8     :: x_g,	      x_s,	      x_t
real*8     :: y_g,	      y_s,	      y_t
real*8     :: delta_g(n_var), delta_s(n_var), delta_t(n_var)
real*8     :: xjac, BigR,  BigR_x
real*8     :: v,     v_x,     v_y,     v_p,	v_s,	 v_t,	  v_ss,     v_st,     v_tt,   v_xs,   v_ys,   v_xt,   v_yt,   v_xx,   v_yy
real*8     :: ps0,   ps0_x,   ps0_y,   ps0_p,	ps0_s,   ps0_t
real*8     :: psi,   psi_x,   psi_y,   psi_p,	psi_s,   psi_t,   psi_ss,   psi_st,   psi_tt, psi_xs, psi_ys, psi_xt, psi_yt, psi_xx, psi_yy
real*8     :: u0,    u0_x,    u0_y,    u0_p,	u0_s,	 u0_t, vv2
real*8     :: u,     u_x,     u_y,     u_p,	u_s,	 u_t,  delta_u_x, delta_u_y
real*8     :: zj0,   zj0_x,   zj0_y,   zj0_p,	zj0_s,   zj0_t
real*8     :: zj,    zj_x,    zj_y,    zj_p,	zj_s,	 zj_t
real*8     :: w0,    w0_x,    w0_y,    w0_p,	w0_s,	 w0_t,    w0_ss,    w0_st,    w0_tt
real*8     :: w,     w_x,     w_y,     w_p,	w_s,	 w_t,	  w_ss,     w_st,     w_tt
real*8     :: Vpar0, Vpar0_x, Vpar0_y, Vpar0_p, Vpar0_s, Vpar0_t, Vpar0_ss, Vpar0_st, Vpar0_tt
real*8     :: Vpar,  Vpar_x,  Vpar_y,  Vpar_p,  Vpar_s,  Vpar_t,  Vpar_ss,  Vpar_st,  Vpar_tt
real*8     :: r0,    r0_x,    r0_y,    r0_p,	r0_s,	 r0_t,    r0_ss,    r0_st,    r0_tt,  r0_hat,  r0_x_hat,  r0_y_hat
real*8     :: rho,   rho_x,   rho_y,   rho_p,	rho_s,   rho_t,   rho_ss,   rho_st,   rho_tt, rho_hat, rho_x_hat, rho_y_hat
real*8     :: Ti0,   Ti0_x,   Ti0_y,   Ti0_p,	Ti0_s,   Ti0_t,   Ti0_ss,   Ti0_st,   Ti0_tt
real*8     :: Ti,    Ti_x,    Ti_y,    Ti_p,	Ti_s,	 Ti_t,    Ti_ss,    Ti_st,    Ti_tt
real*8     :: Te0,   Te0_x,   Te0_y,   Te0_p,	Te0_s,   Te0_t,   Te0_ss,   Te0_st,   Te0_tt
real*8     :: Te,    Te_x,    Te_y,    Te_p,	Te_s,	 Te_t,    Te_ss,    Te_st,    Te_tt
real*8     :: P0,    P0_x,    P0_y,    P0_s,	P0_t,	 P0_p
real*8     :: KPe0
real*8     :: KPe
real*8     :: Bgrad_rho, Bgrad_rho_star, Bgrad_rho_psi, Bgrad_rho_star_psi, Bgrad_rho_rho, Bgrad_rho_rho_n, Bgrad_rho_k_star
real*8     :: Bgrad_Ti,  Bgrad_Ti_star,  Bgrad_Ti_psi,  Bgrad_Ti_star_psi,  Bgrad_Ti_Ti,   Bgrad_Ti_Ti_n,   Bgrad_Ti_k_star
real*8     :: Bgrad_Te,  Bgrad_Te_star,  Bgrad_Te_psi,  Bgrad_Te_star_psi,  Bgrad_Te_Te,   Bgrad_Te_Te_n,   Bgrad_Te_k_star
real*8     :: TG_num1, TG_num2, TG_num5, TG_num6, TG_num7, TG_num8
real*8     :: Jb, dJb_psi, dJb_rho, dJb_Ti, dJb_Te

real*8     :: rhs_ij_1,   rhs_ij_2,   rhs_ij_3,   rhs_ij_4, rhs_ij_5, rhs_ij_6, rhs_ij_7, rhs_ij_8, rhs_ij_9
real*8     :: rhs_ij_5_k, rhs_ij_6_k, rhs_ij_8_k, rhs_ij_9_k
real*8     :: amat_11, amat_12, amat_13, amat_15, amat_16, amat_18
real*8     :: amat_21, amat_22, amat_23, amat_24, amat_25, amat_26, amat_28
real*8     :: amat_31, amat_33
real*8     :: amat_42, amat_44
real*8     :: amat_51, amat_52, amat_55, amat_57
real*8     :: amat_61, amat_62, amat_65, amat_66, amat_67, amat_68
real*8     :: amat_71, amat_72, amat_75, amat_76, amat_77, amat_78
real*8     :: amat_81, amat_82, amat_83, amat_85, amat_86, amat_87, amat_88, amat_89
real*8     :: amat_91, amat_98, amat_99
real*8     :: amat_12_n
real*8     :: amat_23_n
real*8     :: amat_51_k, amat_55_k, amat_55_n, amat_55_kn, amat_57_k, amat_57_n
real*8     :: amat_61_k, amat_65_k, amat_66_k, amat_66_n,  amat_66_kn, amat_67_k,  amat_67_n
real*8     :: amat_75_k, amat_75_n, amat_76_n, amat_77_k,  amat_77_n, amat_77_kn, amat_78_n
real*8     :: amat_81_k, amat_81_n, amat_85_k, amat_85_n,  amat_87_k, amat_87_n,  amat_88_k, amat_88_n, amat_88_kn
real*8     :: amat_91_k, amat_98_k, amat_98_n, amat_98_kn

ELM_p  => thread_struct(tid)%ELM_p 
ELM_n  => thread_struct(tid)%ELM_n 
ELM_k  => thread_struct(tid)%ELM_k
ELM_kn => thread_struct(tid)%ELM_kn
RHS_p  => thread_struct(tid)%RHS_p 
RHS_k  => thread_struct(tid)%RHS_k 

ELM_p = 0.d0; ELM_n = 0.d0; ELM_k = 0.d0; ELM_kn = 0.d0; RHS_p = 0.d0; RHS_k = 0.d0; ELM = 0.d0; RHS = 0.d0

! --- Take time evolution parameters from phys_module
theta = time_evol_theta
zeta  = time_evol_zeta

! Taylor Galerkin (TG2) stabilisation switches
!TG_num1    = 0.d0; TG_num2    = 2.d0; TG_num5    = 2.d0; TG_num6    = 0.d0; TG_num7    = 0.d0; TG_num8    = 0.d0;
TG_num1    = 0.d0; TG_num2    = 0.d0; TG_num5    = 0.d0; TG_num6    = 0.d0; TG_num7    = 0.d0; TG_num8    = 0.d0;





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! Begin integration loop over Gaussian integration points !!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do ms=1, n_gauss

  do mt=1, n_gauss

    wst = wgauss(ms)*wgauss(mt)

    x_g = 0.d0 ; x_s = 0.d0 ; x_t = 0.d0
    y_g = 0.d0 ; y_s = 0.d0 ; y_t = 0.d0

    current_source = 0.d0 ; particle_source = 0.d0 ; heat_source_i = 0.d0 ; heat_source_e = 0.d0

    do i=1,n_vertex_max
      do j=1,n_order+1

        x_g = x_g + nodes(i)%x(j,1) * element%size(i,j) * H(i,j,ms,mt)
        x_s = x_s + nodes(i)%x(j,1) * element%size(i,j) * H_s(i,j,ms,mt)
        x_t = x_t + nodes(i)%x(j,1) * element%size(i,j) * H_t(i,j,ms,mt)

        y_g = y_g + nodes(i)%x(j,2) * element%size(i,j) * H(i,j,ms,mt)
        y_s = y_s + nodes(i)%x(j,2) * element%size(i,j) * H_s(i,j,ms,mt)
        y_t = y_t + nodes(i)%x(j,2) * element%size(i,j) * H_t(i,j,ms,mt)

      enddo
    enddo
    
    xjac    = x_s*y_t - x_t*y_s

    BigR    = x_g
    BigR_x  = 1.d0

    eps_cyl = 1.d0	    ! for cylinder geometry : epscyl = eps

    do mp = 1, n_plane

      ! -----------------------------------------
      ! --- Reconstruct variables from eq_ values
      ! -----------------------------------------
      ps0   = 0.d0; ps0_p   = 0.d0; ps0_s   = 0.d0; ps0_t   = 0.d0
      u0    = 0.d0; u0_p    = 0.d0; u0_s    = 0.d0; u0_t    = 0.d0
      zj0   = 0.d0; zj0_p   = 0.d0; zj0_s   = 0.d0; zj0_t   = 0.d0
      w0    = 0.d0; w0_p    = 0.d0; w0_s    = 0.d0; w0_t    = 0.d0; w0_ss    = 0.d0; w0_tt    = 0.d0; w0_st    = 0.d0
      r0    = 0.d0; r0_p    = 0.d0; r0_s    = 0.d0; r0_t    = 0.d0; r0_ss    = 0.d0; r0_tt    = 0.d0; r0_st    = 0.d0
      Ti0   = 0.d0; Ti0_p   = 0.d0; Ti0_s   = 0.d0; Ti0_t   = 0.d0; Ti0_ss   = 0.d0; Ti0_tt   = 0.d0; Ti0_st   = 0.d0
      Te0   = 0.d0; Te0_p   = 0.d0; Te0_s   = 0.d0; Te0_t   = 0.d0; Te0_ss   = 0.d0; Te0_tt   = 0.d0; Te0_st   = 0.d0
      Vpar0 = 0.d0; Vpar0_p = 0.d0; Vpar0_s = 0.d0; Vpar0_t = 0.d0; Vpar0_ss = 0.d0; Vpar0_tt = 0.d0; Vpar0_st = 0.d0
      Kpe0  = 0.d0
      delta_g = 0.d0 ; delta_s = 0.d0 ; delta_t = 0.d0

      do i=1,n_vertex_max
        do j=1,n_order+1
 	  do in=1,n_tor

            ps0      = ps0	+ nodes(i)%values(in,j,1) * element%size(i,j) * H(i,j,ms,mt)  * HZ(in,mp)
            ps0_p    = ps0_p	+ nodes(i)%values(in,j,1) * element%size(i,j) * H(i,j,ms,mt)  * HZ_p(in,mp)
            ps0_s    = ps0_s	+ nodes(i)%values(in,j,1) * element%size(i,j) * H_s(i,j,ms,mt)* HZ(in,mp)
            ps0_t    = ps0_t	+ nodes(i)%values(in,j,1) * element%size(i,j) * H_t(i,j,ms,mt)* HZ(in,mp)

            u0       = u0	+ nodes(i)%values(in,j,2) * element%size(i,j) * H(i,j,ms,mt)  * HZ(in,mp)
            u0_p     = u0_p	+ nodes(i)%values(in,j,2) * element%size(i,j) * H(i,j,ms,mt)  * HZ_p(in,mp)
            u0_s     = u0_s	+ nodes(i)%values(in,j,2) * element%size(i,j) * H_s(i,j,ms,mt)* HZ(in,mp)
            u0_t     = u0_t	+ nodes(i)%values(in,j,2) * element%size(i,j) * H_t(i,j,ms,mt)* HZ(in,mp)

            zj0      = zj0	+ nodes(i)%values(in,j,3) * element%size(i,j) * H(i,j,ms,mt)  * HZ(in,mp)
            zj0_p    = zj0_p	+ nodes(i)%values(in,j,3) * element%size(i,j) * H(i,j,ms,mt)  * HZ_p(in,mp)
            zj0_s    = zj0_s	+ nodes(i)%values(in,j,3) * element%size(i,j) * H_s(i,j,ms,mt)* HZ(in,mp)
            zj0_t    = zj0_t	+ nodes(i)%values(in,j,3) * element%size(i,j) * H_t(i,j,ms,mt)* HZ(in,mp)

            w0       = w0	+ nodes(i)%values(in,j,4) * element%size(i,j) * H(i,j,ms,mt)  * HZ(in,mp)
            w0_p     = w0_p	+ nodes(i)%values(in,j,4) * element%size(i,j) * H(i,j,ms,mt)  * HZ_p(in,mp)
            w0_s     = w0_s	+ nodes(i)%values(in,j,4) * element%size(i,j) * H_s(i,j,ms,mt)* HZ(in,mp)
            w0_t     = w0_t	+ nodes(i)%values(in,j,4) * element%size(i,j) * H_t(i,j,ms,mt)* HZ(in,mp)
            w0_ss    = w0_ss	+ nodes(i)%values(in,j,4) * element%size(i,j) * H_ss(i,j,ms,mt)*HZ(in,mp)
            w0_tt    = w0_tt	+ nodes(i)%values(in,j,4) * element%size(i,j) * H_tt(i,j,ms,mt)*HZ(in,mp)
            w0_st    = w0_st	+ nodes(i)%values(in,j,4) * element%size(i,j) * H_st(i,j,ms,mt)*HZ(in,mp)

            r0       = r0	+ nodes(i)%values(in,j,5) * element%size(i,j) * H(i,j,ms,mt)  * HZ(in,mp)
            r0_p     = r0_p	+ nodes(i)%values(in,j,5) * element%size(i,j) * H(i,j,ms,mt)  * HZ_p(in,mp)
            r0_s     = r0_s	+ nodes(i)%values(in,j,5) * element%size(i,j) * H_s(i,j,ms,mt)* HZ(in,mp)
            r0_t     = r0_t	+ nodes(i)%values(in,j,5) * element%size(i,j) * H_t(i,j,ms,mt)* HZ(in,mp)
            r0_ss    = r0_ss	+ nodes(i)%values(in,j,5) * element%size(i,j) * H_ss(i,j,ms,mt)*HZ(in,mp)
            r0_tt    = r0_tt	+ nodes(i)%values(in,j,5) * element%size(i,j) * H_tt(i,j,ms,mt)*HZ(in,mp)
            r0_st    = r0_st	+ nodes(i)%values(in,j,5) * element%size(i,j) * H_st(i,j,ms,mt)*HZ(in,mp)

            Ti0      = Ti0	+ nodes(i)%values(in,j,6) * element%size(i,j) * H(i,j,ms,mt)  * HZ(in,mp)
            Ti0_p    = Ti0_p	+ nodes(i)%values(in,j,6) * element%size(i,j) * H(i,j,ms,mt)  * HZ_p(in,mp)
            Ti0_s    = Ti0_s	+ nodes(i)%values(in,j,6) * element%size(i,j) * H_s(i,j,ms,mt)* HZ(in,mp)
            Ti0_t    = Ti0_t	+ nodes(i)%values(in,j,6) * element%size(i,j) * H_t(i,j,ms,mt)* HZ(in,mp)
            Ti0_ss   = Ti0_ss	+ nodes(i)%values(in,j,6) * element%size(i,j) * H_ss(i,j,ms,mt)*HZ(in,mp)
            Ti0_tt   = Ti0_tt	+ nodes(i)%values(in,j,6) * element%size(i,j) * H_tt(i,j,ms,mt)*HZ(in,mp)
            Ti0_st   = Ti0_st	+ nodes(i)%values(in,j,6) * element%size(i,j) * H_st(i,j,ms,mt)*HZ(in,mp)

            Te0      = Te0	+ nodes(i)%values(in,j,8) * element%size(i,j) * H(i,j,ms,mt)  * HZ(in,mp)
            Te0_p    = Te0_p	+ nodes(i)%values(in,j,8) * element%size(i,j) * H(i,j,ms,mt)  * HZ_p(in,mp)
            Te0_s    = Te0_s	+ nodes(i)%values(in,j,8) * element%size(i,j) * H_s(i,j,ms,mt)* HZ(in,mp)
            Te0_t    = Te0_t	+ nodes(i)%values(in,j,8) * element%size(i,j) * H_t(i,j,ms,mt)* HZ(in,mp)
            Te0_ss   = Te0_ss	+ nodes(i)%values(in,j,8) * element%size(i,j) * H_ss(i,j,ms,mt)*HZ(in,mp)
            Te0_tt   = Te0_tt	+ nodes(i)%values(in,j,8) * element%size(i,j) * H_tt(i,j,ms,mt)*HZ(in,mp)
            Te0_st   = Te0_st	+ nodes(i)%values(in,j,8) * element%size(i,j) * H_st(i,j,ms,mt)*HZ(in,mp)

            Vpar0    = Vpar0	+ nodes(i)%values(in,j,7) * element%size(i,j) * H(i,j,ms,mt)  * HZ(in,mp)
            Vpar0_p  = Vpar0_p  + nodes(i)%values(in,j,7) * element%size(i,j) * H(i,j,ms,mt)  * HZ_p(in,mp)
            Vpar0_s  = Vpar0_s  + nodes(i)%values(in,j,7) * element%size(i,j) * H_s(i,j,ms,mt)* HZ(in,mp)
            Vpar0_t  = Vpar0_t  + nodes(i)%values(in,j,7) * element%size(i,j) * H_t(i,j,ms,mt)* HZ(in,mp)
            Vpar0_ss = Vpar0_ss + nodes(i)%values(in,j,7) * element%size(i,j) * H_ss(i,j,ms,mt)*HZ(in,mp)
            Vpar0_tt = Vpar0_tt + nodes(i)%values(in,j,7) * element%size(i,j) * H_tt(i,j,ms,mt)*HZ(in,mp)
            Vpar0_st = Vpar0_st + nodes(i)%values(in,j,7) * element%size(i,j) * H_st(i,j,ms,mt)*HZ(in,mp)

!            Kpe0     = Kpe0	+ nodes(i)%values(in,j,9) * element%size(i,j) * H(i,j,ms,mt)  * HZ(in,mp)

 	    do k=1,n_var
              delta_g(k) = delta_g(k) + nodes(i)%deltas(in,j,k) * element%size(i,j) * H(i,j,ms,mt)   * HZ(in,mp)
              delta_s(k) = delta_s(k) + nodes(i)%deltas(in,j,k) * element%size(i,j) * H_s(i,j,ms,mt) * HZ(in,mp)
              delta_t(k) = delta_t(k) + nodes(i)%deltas(in,j,k) * element%size(i,j) * H_t(i,j,ms,mt) * HZ(in,mp)
	    enddo      		    
      
          enddo
        enddo
      enddo

      ps0_x    = (   y_t * ps0_s - y_s * ps0_t ) / xjac
      ps0_y    = ( - x_t * ps0_s + x_s * ps0_t ) / xjac
      u0_x     = (   y_t * u0_s - y_s * u0_t ) / xjac
      u0_y     = ( - x_t * u0_s + x_s * u0_t ) / xjac
      vv2      = BigR**2 *  ( u0_x * u0_x + u0_y *u0_y  )
      zj0_x    = (   y_t * zj0_s - y_s * zj0_t ) / xjac
      zj0_y    = ( - x_t * zj0_s + x_s * zj0_t ) / xjac
      w0_x     = (   y_t * w0_s - y_s * w0_t ) / xjac
      w0_y     = ( - x_t * w0_s + x_s * w0_t ) / xjac
      r0_x     = (   y_t * r0_s - y_s * r0_t ) / xjac
      r0_y     = ( - x_t * r0_s + x_s * r0_t ) / xjac
      r0_hat   = BigR**2 * r0
      r0_x_hat = 2.d0 * BigR * BigR_x  * r0 + BigR**2 * r0_x
      r0_y_hat = BigR**2 * r0_y
      Ti0_x    = (   y_t * Ti0_s - y_s * Ti0_t ) / xjac
      Ti0_y    = ( - x_t * Ti0_s + x_s * Ti0_t ) / xjac
      Te0_x    = (   y_t * Te0_s - y_s * Te0_t ) / xjac
      Te0_y    = ( - x_t * Te0_s + x_s * Te0_t ) / xjac
      Vpar0_x  = (   y_t * Vpar0_s - y_s * Vpar0_t ) / xjac
      Vpar0_y  = ( - x_t * Vpar0_s + x_s * Vpar0_t ) / xjac
      delta_u_x= (   y_t * delta_s(2) - y_s * delta_t(2) ) / xjac
      delta_u_y= ( - x_t * delta_s(2) + x_s * delta_t(2) ) / xjac
      P0       = r0   * (Ti0 + Te0)
      P0_x     = r0_x * (Ti0 + Te0) + r0 * (Ti0_x + Te0_x)
      P0_y     = r0_y * (Ti0 + Te0) + r0 * (Ti0_y + Te0_y)
      P0_s     = r0_s * (Ti0 + Te0) + r0 * (Ti0_s + Te0_s)
      P0_t     = r0_t * (Ti0 + Te0) + r0 * (Ti0_t + Te0_t)
      P0_p     = r0_p * (Ti0 + Te0) + r0 * (Ti0_p + Te0_p)

      
      ! -------------------------------------
      ! --- Temperature dependent resistivity
      ! -------------------------------------
      if ( eta_T_dependent ) then
 	eta_Te     = eta   * (abs(Te0)/Te_0)**(-1.5d0)
 	deta_dTe   = - eta   * (1.5d0)  * abs(Te0)**(-2.5d0) * Te_0**(1.5d0)
 	d2eta_d2Te =   eta   * (3.75d0) * abs(Te0)**(-3.5d0) * Te_0**(1.5d0)
      else
 	eta_Te     = eta
 	deta_dTe   = 0.d0
 	d2eta_d2Te = 0.d0
      end if
      
      ! -----------------------------------
      ! --- Temperature dependent viscosity
      ! -----------------------------------
      if ( visco_T_dependent ) then
 	visco_Te   = visco * (abs(Te0)/Te_0)**(-1.5d0)
 	dvisco_dTe = - visco * (1.5d0)  * abs(Te0)**(-2.5d0) * Te_0**(1.5d0)
      else
 	visco_Te   = visco
 	dvisco_dTe = 0.d0
      end if
      
      ! -------------------------------------------------------------
      ! --- D_perp and K_perp profiles (for fixed pedestal gradients)
      ! -------------------------------------------------------------
      ! --- First need psi_norm
      psi_norm = (ps0 - psi_axis)/(psi_bnd - psi_axis)
      if (xpoint2) then
       if ((psi_norm .lt. 1.d0) .and. (y_g .lt. Z_xpoint(1)) .and. (xcase2 .ne. 2)) then
         psi_norm = 2.d0 - psi_norm
       endif
       if ((psi_norm .lt. 1.d0) .and. (y_g .gt. Z_xpoint(2)) .and. (xcase2 .ne. 1)) then
         psi_norm = 2.d0 - psi_norm
       endif
      endif

      ! --- Take values from input file (rho_coef, T_coef...) 
      do id = 1,10
	if ((id .eq. 7) .or. (id .eq. 8) .or. (id .eq. 9)) then
	  Diff(1,id) = rho_coef(id-6) ; Diff(2,id) = Ti_coef(id-6) ; Diff(3,id) = Te_coef(id-6)
	else
          Diff(1,id) = D_perp(id) ; Diff(2,id) = ZK_i_perp(id) ; Diff(3,id) = ZK_e_perp(id)
	endif
      enddo      
      if (Diff(1,10) .eq. 1.d0) then 
        Diff(1,4) = rho_coef(4) ; Diff(1,5) = rho_coef(5) 
      endif
      if (Diff(2,10) .eq. 1.d0) then 
        Diff(2,4) = Ti_coef(4)  ; Diff(2,5) = Ti_coef(5) 
      endif
      if (Diff(3,10) .eq. 1.d0) then 
        Diff(3,4) = Te_coef(4)  ; Diff(3,5) = Te_coef(5) 
      endif
            
      ! --- Build profiles
      do id = 1,3
  	if (psi_norm .gt. Diff(id,5)) then
  	  if (id .eq. 1) then
	    psi_D = 2.d0*Diff(id,5) - psi_norm
	  else
	    psi_D = Diff(id,5)
	  endif
  	else 
  	  psi_D = psi_norm
  	endif
  	if (psi_norm .lt. 0.5d0) then
  	  psi_D = 0.5d0
  	endif 
  	if (Diff(id,7) .ge. 0.d0) then
  	  Diff(id,7) = -0.1d0
  	endif 
        if (xcase2 .ne. 2) psi_D = psi_D * (0.5d0 - 0.5d0 * tanh((Z_xpoint(1)-y_g)/0.1d0))
        if (xcase2 .ne. 1) psi_D = psi_D * (0.5d0 - 0.5d0 * tanh((y_g-Z_xpoint(2))/0.1d0))
  	atn_D    = 0.5d0 - 0.5d0 * tanh((psi_D-Diff(id,5))/Diff(id,4))
  	datn_D   =       - 0.5d0 / cosh((psi_D-Diff(id,5))/Diff(id,4))**2.d0 /(Diff(id,4)*(psi_bnd - psi_axis))
  	pol_D    = 1 + Diff(id,7)*psi_D	   + Diff(id,8)*psi_D**2.d0	 + Diff(id,9)*psi_D**3.d0
  	dpol_D   =    (Diff(id,7)       + 2.d0*Diff(id,8)*psi_D	    + 3.d0*Diff(id,9)*psi_D**2.d0)/(psi_bnd - psi_axis)
  	D_min    = 1.d0/( -(1+Diff(id,7)*Diff(id,5)+Diff(id,8)*Diff(id,5)**2.d0+Diff(id,9)*Diff(id,5)**3.d0) * 0.5d0/(Diff(id,4)*(psi_bnd - psi_axis))&
  	           + 0.5d0 * (Diff(id,7)     + 2.d0*Diff(id,8)*Diff(id,5)+ 3.d0*Diff(id,9)*Diff(id,5)**2.d0)/(psi_bnd - psi_axis) )

  	prof(id) = (1.d0-Diff(id,10)) * ( Diff(id,1) * (1.d0-Diff(id,2)+Diff(id,2)*(0.5d0 - 0.5d0 * tanh((psi_norm-Diff(id,5))/Diff(id,4)))) &
  		        	        + Diff(id,6) * (0.5d0 - 0.5d0 * tanh((-psi_norm+Diff(id,5)+Diff(id,3))/Diff(id,4)))) &
  	                +Diff(id,10)  * ( Diff(id,1) / (dpol_D*atn_D + pol_D*datn_D) / D_min ) &
			              * (1 + Diff(id,6) - Diff(id,6) * tanh(-(psi_norm-(1+4*Diff(id,4)))/Diff(id,4)))   !higher Kperp in SOL
      enddo
      
      ! --- Allocate profiles to corresponding names
      if ( num_d_perp ) then
        D_prof  = interpolProf(num_d_perp_x, num_d_perp_y, num_d_perp_len, psi_norm)
      else
        D_prof  = prof(1)
      end if
      
      if ( num_zk_i_perp ) then
        ZK_i_prof = interpolProf(num_zk_i_perp_x, num_zk_i_perp_y, num_zk_i_perp_len, psi_norm)
      else
        ZK_i_prof = prof(2) 
      end if
      
      if ( num_zk_e_perp ) then
        ZK_e_prof = interpolProf(num_zk_e_perp_x, num_zk_e_perp_y, num_zk_e_perp_len, psi_norm)
      else
        ZK_e_prof = prof(3) 
      end if

      ! --- Avoid negative density 
      if ( r0 .lt. rho_1 ) then
	D_prof = D_prof * 1.d0
      end if	  
      if ( Te0 .lt. Te_1 ) then
	ZK_e_prof = ZK_e_prof * 1.d0
      end if
      if ( Ti0 .lt. Ti_1 ) then
	ZK_i_prof = ZK_i_prof * 1.d0
      end if

      ! -----------------------------------------------------
      ! --- Parallel conductivity profiles (Braginskii model)
      ! -----------------------------------------------------
      ZK_i_par  = K_i_par * abs(Ti0)**(2.5d0) ! = K_i_par
      ZK_e_par  = K_e_par * abs(Te0)**(2.5d0) ! = K_e_par
      
      dZK_i_par  = K_i_par * (2.5d0) * abs(Ti0)**(1.5d0) ! = 0.d0
      dZK_e_par  = K_e_par * (2.5d0) * abs(Te0)**(1.5d0) ! = 0.d0
      
      Q = Q_bar * (abs(r0)**(2.d0)) * abs(Te0)**(-1.5d0)
      dQ_drho = (2.d0) * Q_bar * abs(r0) * abs(Te0)**(-1.5d0)
      dQ_dTe = (-1.5d0) * Q_bar * (abs(r0)**(2.d0)) * abs(Te0)**(-2.5d0)
     
      ! -------------------------
      ! --- Hyper diffusivitities
      ! -------------------------
      eta_numm       = eta_num  		     ! hyperresistivity
      visco_numm     = visco_num		     ! hyperviscosity
      visco_par_numm = visco_par_num		     ! hyperviscosity
      D_perp_numm    = D_perp_num		     ! hyperdiffusivity
      K_perp_numm    = ZK_perp_num		     ! hyperconductivity
      
      if (psi_norm .lt. 0.4d0) eta_numm   = eta_numm   * 1.d2
      if (psi_norm .lt. 0.4d0) visco_numm = visco_numm * 1.d2

      ! ---------------------------------------------------
      ! --- Bootstrap current coefficients (Wesson formula)
      ! ---------------------------------------------------
      if (bootstrap) then
        call bootstrap_current_rhs(BigR, minRad, R_axis, &
                                   psi_axis, psi_bnd,    &
                                   ps0, ps0_x, ps0_y,    &
                                   r0,  r0_x,  r0_y,     &
                                   Ti0, Ti0_x, Ti0_y,    &
                                   Te0, Te0_x, Te0_y,    &
                                   Jb)
      else
        Jb = 0.d0
      endif

      ! --------------------------------------------------------------
      ! --- Heating, current and particle source (the same for all mp)
      ! --------------------------------------------------------------
      if (mp .eq. 1) then
        call current(xpoint2, xcase2, x_g,y_g, Z_xpoint, ps0,psi_axis,psi_bnd,current_source)

        !call sources(xpoint2, xcase2, y_g, Z_xpoint, ps0,psi_axis,psi_bnd,particle_source,heat_source_i,heat_source_e)
        ! --- New source profile: source with exactly the same profile as the initial equilibirum profiles.
        call density(	   xpoint2, xcase2, y_g, Z_xpoint, ps0,psi_axis,psi_bnd, &
        		   zn,dn_dpsi,  dn_dz, dn_dpsi2, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz)

        call temperature_i(xpoint2, xcase2, y_g, Z_xpoint, ps0,psi_axis,psi_bnd, &
        		   zTi,dTi_dpsi,dTi_dz,dTi_dpsi2,dTi_dz2,dTi_dpsi_dz,dTi_dpsi3,dTi_dpsi_dz2,dTi_dpsi2_dz)
               
        call temperature_e(xpoint2, xcase2, y_g, Z_xpoint, ps0,psi_axis,psi_bnd, &
        		   zTe,dTe_dpsi,dTe_dz,dTe_dpsi2,dTe_dz2,dTe_dpsi_dz,dTe_dpsi3,dTe_dpsi_dz2,dTe_dpsi2_dz)

        particle_source = particlesource * ( zn -  r0 )
        heat_source_i	= heatsource_i   * ( zTi - Ti0 )
        heat_source_e	= heatsource_e   * ( zTe - Te0 )

        particle_source = particle_source * ( 0.5d0 - 0.5d0 * tanh((psi_norm-0.99)/0.005) )
        heat_source_i	= heat_source_i   * ( 0.5d0 - 0.5d0 * tanh((psi_norm-0.99)/0.005) )
        heat_source_e	= heat_source_e   * ( 0.5d0 - 0.5d0 * tanh((psi_norm-0.99)/0.005) )
      endif

      ! ------------------------------------
      ! --- Now the equations, first the RHS
      ! ------------------------------------
      do i=1,n_vertex_max

 	do j=1,n_order+1

          ! --- If we're doing the fft, don't loop...
          n_tor_loop = n_tor
          if (n_tor .gt. 3) n_tor_loop = 1
	  
	  do im=1,n_tor_loop

            ! --- Build up test function and parallel gradient terms            
	    if (n_tor .gt. 3) then
	      HHZ(im,mp)   = 1.d0
	      HHZ_p(im,mp) = 1.d0
 	      index_ij = n_var*(n_order+1)*(i-1) + n_var * (j-1) + 1   ! index in the ELM matrix
	    else
	      HHZ(im,mp)   = HZ(im,mp)
	      HHZ_p(im,mp) = HZ_p(im,mp)
	      index_ij = n_tor*n_var*(n_order+1)*(i-1) + n_tor * n_var * (j-1) + im   ! index in the ELM matrix
	    endif
            
            v	=  H(i,j,ms,mt) * element%size(i,j) * HHZ(im,mp)
            v_x = (  y_t * h_s(i,j,ms,mt) - y_s * h_t(i,j,ms,mt) ) * element%size(i,j) / xjac * HHZ(im,mp)
            v_y = (- x_t * h_s(i,j,ms,mt) + x_s * h_t(i,j,ms,mt) ) * element%size(i,j) / xjac * HHZ(im,mp)

            v_s = h_s(i,j,ms,mt) * element%size(i,j) * HHZ(im,mp)
            v_t = h_t(i,j,ms,mt) * element%size(i,j) * HHZ(im,mp)
            v_p = H(i,j,ms,mt)   * element%size(i,j) * HHZ_p(im,mp)

            v_ss = h_ss(i,j,ms,mt) * element%size(i,j) * HHZ(im,mp)
            v_tt = h_tt(i,j,ms,mt) * element%size(i,j) * HHZ(im,mp)
            v_st = h_st(i,j,ms,mt) * element%size(i,j) * HHZ(im,mp)

 	    Bgrad_rho_star   = ( v_x  * ps0_y - v_y  * ps0_x ) / BigR			      ! only the part without	 d/d(phi)
 	    Bgrad_rho_k_star = ( F0 / BigR * v_p )	     / BigR			      ! only the part containing d/d(phi)
 	    Bgrad_rho	     = ( F0 / BigR * r0_p +  r0_x * ps0_y - r0_y * ps0_x ) / BigR
 	    Bgrad_Ti_star     = ( v_x  * ps0_y - v_y  * ps0_x ) / BigR  		      ! only the part without	 d/d(phi)
 	    Bgrad_Ti_k_star   = ( F0 / BigR * v_p	    ) / BigR			      ! only the part containing d/d(phi)
 	    Bgrad_Ti	      = ( F0 / BigR * Ti0_p +  Ti0_x * ps0_y - Ti0_y * ps0_x ) / BigR
 	    Bgrad_Te_star     = ( v_x  * ps0_y - v_y  * ps0_x ) / BigR  		      ! only the part without	 d/d(phi)
 	    Bgrad_Te_k_star   = ( F0 / BigR * v_p	    ) / BigR			      ! only the part containing d/d(phi)
 	    Bgrad_Te	      = ( F0 / BigR * Te0_p +  Te0_x * ps0_y - Te0_y * ps0_x ) / BigR

 	    BB2 = (F0*F0 + ps0_x * ps0_x + ps0_y * ps0_y )/BigR**2


    	    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    	    !!!!!!!!!! RHS equation 1 (psi - induction) !!!!!!!!!!!!
    	    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    	    rhs_ij_1 = + v * eta_Te  * (zj0 - current_source + Jb)/ BigR   * xjac * tstep &
    		       + v * (ps0_x * u0_y - ps0_y * u0_x)		   * xjac * tstep &
    		       - v * eps_cyl * F0 / BigR  * u0_p		   * xjac * tstep &
    		       + eta_numm * (v_x * zj0_x + v_y * zj0_y) 	   * xjac * tstep &
    		       + zeta * v / BigR				   * xjac * delta_g(1)
    				 
    				 
    	    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    	    !!!!!!!!!! RHS equation 2 (U - momentum) !!!!!!!!!!!!
    	    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    	    rhs_ij_2 = - 0.5d0 * vv2 * (v_x * r0_y_hat - v_y * r0_x_hat)	     * xjac * tstep &
    		       - r0_hat * BigR**2 * w0 * (v_x * u0_y - v_y * u0_x)	     * xjac * tstep &
    		       + v * (ps0_x * zj0_y - ps0_y * zj0_x )			     * xjac * tstep &
    		       - visco_Te * BigR * (v_x * w0_x + v_y * w0_y)		     * xjac * tstep &
    		       - v * eps_cyl * F0 / BigR * zj0_p			     * xjac * tstep &
    		       + BigR**2 * (v_x * p0_y - v_y * p0_x)			     * xjac * tstep &
    		       - zeta * BigR * r0_hat * (v_x * delta_u_x + v_y * delta_u_y)  * xjac	    &		      
    		       - visco_numm  *  							    & 
    		         (    ( 	v_ss  * (x_t**2+y_t**2) 				    & 
    		     	       +	v_tt  * (x_s**2+y_s**2) 				    &
    		     	       - 2.d0 * v_st  * (x_s*x_t + y_s*y_t) )				    &
    		            * ( 	w0_ss * (x_t**2+y_t**2) 				    &
    		     	       +	w0_tt * (x_s**2+y_s**2) 				    &
    		     	       - 2.d0 * w0_st * (x_s*x_t + y_s*y_t) ) ) 			    &
    		         / xjac**4						     * xjac * tstep &
    		       - TG_num2 * tstep * 0.25d0 * r0_hat * BigR**3				    &
		        	 * (w0_x * u0_y - w0_y * u0_x)  				    &
    		     	     	 * ( v_x * u0_y - v_y  * u0_x)  		     * xjac * tstep

    	    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    	    !!!!!!!!!! RHS equation 3 (j - current) !!!!!!!!!!!!
    	    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    	    rhs_ij_3 = 0.d0 !- ( v_x * ps0_x  + v_y * ps0_y + v*zj0) / BigR * xjac * tstep

    	    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    	    !!!!!!!!!! RHS equation 4 (w - vorticity) !!!!!!!!!!
    	    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    	    rhs_ij_4 = 0.d0 !- ( v_x * u0_x   + v_y * u0_y  + v*w0)  * BigR * xjac * tstep

    	    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    	    !!!!!!!!!! RHS equation 5 (rho - continuity) !!!!!!!!!
    	    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    	    rhs_ij_5 =  v * BigR * particle_source				      * xjac * tstep &
    		      + v * BigR**2 * ( r0_x * u0_y - r0_y * u0_x)		      * xjac * tstep &
    		      + v * 2.d0 * BigR * r0 * u0_y				      * xjac * tstep &
    		      - (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star * Bgrad_rho      * xjac * tstep &
    		      - D_prof * BigR  * (v_x*r0_x + v_y*r0_y)  		      * xjac * tstep &
    		      - v * F0 / BigR * Vpar0 * r0_p				      * xjac * tstep &
    		      - v * Vpar0 * (r0_x * ps0_y - r0_y * ps0_x)		      * xjac * tstep &
    		      - v * F0 / BigR * r0 * vpar0_p				      * xjac * tstep &
    		      - v * r0 * (vpar0_x * ps0_y - vpar0_y * ps0_x)		      * xjac * tstep &
    		      + zeta * v * BigR 					      * xjac *delta_g(5)& 
    		      - D_perp_numm *								     &
    		        (    (         v_ss  * (x_t**2 + y_t**2)				     &
    		     	      +        v_tt  * (x_s**2 + y_s**2)				     &
    		     	      - 2.d0 * v_st  * (x_s*x_t + y_s*y_t) )				     &
    		           * (         r0_ss * (x_t**2 + y_t**2)				     &
    		     	      +        r0_tt * (x_s**2 + y_s**2)				     &
    		     	      - 2.d0 * r0_st * (x_s*x_t + y_s*y_t) ) )  			     &
    		        / xjac**3 *tstep							     &
    		      - TG_num5 * tstep * 0.25d0 * BigR**3 * (r0_x * u0_y - r0_y * u0_x)	     &
    		     	     	* ( v_x * u0_y - v_y  * u0_x)			      * xjac * tstep &
    		      - TG_num5 * tstep * 0.25d0 / BigR * vpar0**2				     &
    		     	     	* (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)		     &
    		     	     	* ( v_x * ps0_y -  v_y * ps0_x )		      * xjac * tstep

    	    rhs_ij_5_k = - (D_par-D_prof) * BigR / BB2 * Bgrad_rho_k_star * Bgrad_rho * xjac * tstep &
    		    	 - D_prof * BigR  * ( v_p*r0_p * eps_cyl**2 /BigR**2 )        * xjac * tstep &
    		    	 - TG_num5 * tstep * 0.25d0 / BigR * vpar0**2				     &
    		     	     	   * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)		     &
    		     	     	   * (  			    F0 / BigR * v_p ) * xjac * tstep

    	    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    	    !!!!!!!!!! RHS equation 6 (Ti - Ion energy) !!!!!!!!!
    	    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    	    rhs_ij_6 =   v * BigR * heat_source_i					   * xjac * tstep &
    		       + v * r0 * BigR**2 * ( Ti0_x * u0_y - Ti0_y * u0_x)		   * xjac * tstep &
    		       + v * r0 * 2.d0* (GAMMA-1.d0) * BigR * Ti0 * u0_y		   * xjac * tstep &
    		       - v * r0 * F0 / BigR * Vpar0 * Ti0_p				   * xjac * tstep &
    		       - v * r0 * Vpar0 * (Ti0_x * ps0_y - Ti0_y * ps0_x)		   * xjac * tstep &
    		       - v * r0 * (GAMMA-1.d0) * Ti0 * (vpar0_x * ps0_y - vpar0_y * ps0_x) * xjac * tstep &
    		       - v * r0 * (GAMMA-1.d0) * Ti0 * F0 / BigR * vpar0_p		   * xjac * tstep &
    		       -  (ZK_i_par-ZK_i_prof) * BigR / BB2 * Bgrad_Ti_star * Bgrad_Ti     * xjac * tstep &
    		       -  ZK_i_prof * BigR * (v_x*Ti0_x + v_y*Ti0_y )			   * xjac * tstep &
    		       + zeta * v * r0 * BigR						   * xjac *delta_g(6)&
    		       ! Temperature transfer (not yet tested!)
    		       !+ v * BigR * Q * (Te0 - Ti0)					   * xjac * tstep &
    		       - K_perp_numm  * 								  &
    		         (    ( 	v_ss   * (x_t**2 + y_t**2)					  &
    		     	       +	v_tt   * (x_s**2 + y_s**2)					  &
    		     	       - 2.d0 * v_st   * (x_s*x_t + y_s*y_t) )  				  &
    		            * ( 	Ti0_ss * (x_t**2 + y_t**2)					  &
    		     	       +	Ti0_tt * (x_s**2 + y_s**2)					  &
    		     	       - 2.d0 * Ti0_st * (x_s*x_t + y_s*y_t) ) )				  &
    		         / xjac**4							   * xjac * tstep &
    		       - TG_num6 * tstep * 0.25d0 * BigR**3 * Ti0 * (r0_x * u0_y - r0_y * u0_x) 	  &
    		     	     	 * ( v_x * u0_y - v_y * u0_x)				   * xjac * tstep &
    		       - TG_num6 * tstep * 0.25d0 * BigR**3 * r0 * (Ti0_x * u0_y - Ti0_y * u0_x)	  &
    		     	     	 * ( v_x * u0_y - v_y * u0_x)				   * xjac * tstep

    	    rhs_ij_6_k = - (ZK_i_par-ZK_i_prof) * BigR / BB2 * Bgrad_Ti_k_star * Bgrad_Ti  * xjac * tstep &
    		    	 - ZK_i_prof * BigR * ( v_p*Ti0_p /BigR**2 )			   * xjac * tstep &
    		    	 - TG_num6 * tstep * 0.25d0 / BigR * vpar0**2					  &
    		      	 	   * Ti0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)		  &
    		      	 	   * ( v_x * ps0_y - v_y * ps0_x + F0 / BigR * v_p)	   * xjac * tstep &
    		    	 - TG_num6 * tstep * 0.25d0 / BigR * vpar0**2					  &
    		      	 	   * r0 * (Ti0_x * ps0_y - Ti0_y * ps0_x + F0 / BigR * Ti0_p)		  &
    		      	 	   * ( v_x * ps0_y - v_y * ps0_x + F0 / BigR * v_p)	   * xjac * tstep

    	    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    	    !!!!!!!!!! RHS equation 7 (Vpar - parallel momentum) !!!!!!!!!
    	    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    	    rhs_ij_7 = - v * F0 / BigR * P0_p						      * xjac * tstep &
    		       - v * (P0_x * ps0_y - P0_y * ps0_x)				      * xjac * tstep &
    		       - visco_par * (v_x * vpar0_x + v_y * vpar0_y) * BigR		      * xjac * tstep &
    		       + zeta * v * r0 * F0**2 / BigR					      * xjac *delta_g(7)&  
    		       - visco_par_numm  *								     &
    		         (    ( 	v_ss	 * (x_t**2 + y_t**2)					     &
    		     	       +	v_tt	 * (x_s**2 + y_s**2)					     &
    		     	       - 2.d0 * v_st	 * (x_s*x_t + y_s*y_t) )				     &
    		            * ( 	vpar0_ss * (x_t**2 + y_t**2)					     &
    		     	       +	vpar0_tt * (x_s**2 + y_s**2)					     &
    		     	       - 2.d0 * vpar0_st * (x_s*x_t + y_s*y_t) ) )				     &
    		         / xjac**4							      * xjac * tstep &
    		       - TG_NUM7 * tstep * 0.25d0 * r0 * Vpar0**2 * BB2 				     &
    		     	     	 * (-(ps0_x * vpar0_y - ps0_y * vpar0_x) + F0 / BigR * vpar0_p) / BigR       &
    		     	     	 * (-(ps0_x * v_y     - ps0_y * v_x)	 + F0 / BigR * v_p)   * xjac * tstep &
    		       - TG_NUM7 * tstep * 0.25d0 * v  * Vpar0**2 * BB2 				     &
    		     	     	 * (-(ps0_x * vpar0_y - ps0_y * vpar0_x) + F0 / BigR * vpar0_p) / BigR       &
    		     	     	 * (-(ps0_x * r0_y    - ps0_y * r0_x)	 + F0 / BigR * r0_p)  * xjac * tstep

    	    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    	    !!!!!!!!!! RHS equation 8 (Te - electron energy) !!!!!!!!!
    	    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            rhs_ij_8 =   v * BigR * heat_source_e                                          * xjac * tstep &
                       + v * r0 * BigR**2 * ( Te0_s * u0_t - Te0_t * u0_s)                        * tstep &
                       + v * r0 * 2.d0* (GAMMA-1.d0) * BigR * Te0 * u0_y                   * xjac * tstep &
                       - v * r0 * F0 / BigR * Vpar0 * Te0_p                                * xjac * tstep &
                       - v * r0 * Vpar0 * (Te0_s * ps0_t - Te0_t * ps0_s)                         * tstep &
                       - v * r0 * (GAMMA-1.d0) * Te0 * (vpar0_s * ps0_t - vpar0_t * ps0_s)        * tstep &
                       - v * r0 * (GAMMA-1.d0) * Te0 * F0 / BigR * vpar0_p                 * xjac * tstep &
                       ! Current terms (not yet tested!)
                       !+ sigma * v /BigR * Te0_x * ps0_xp                                 * xjac * tstep &
                       !+ sigma * v /BigR * Te0_y * ps0_yp                                 * xjac * tstep &
                       !- sigma * v /BigR * zj0 * Te0_p                                    * xjac * tstep &
                       !- sigma * (GAMMA-1.d0) * v * Te0/(BigR*r0) * r0_x * ps0_xp         * xjac * tstep &
                       !- sigma * (GAMMA-1.d0) * v * Te0/(BigR*r0) * r0_y * ps0_yp         * xjac * tstep &
                       !+ sigma * (GAMMA-1.d0) * v * Te0/(BigR*r0) * zj0 * r0_p            * xjac * tstep &
                       ! Temperature transfer (not yet tested!)
                       !- v * BigR * Q * (Te0 - Ti0)                                       * xjac * tstep &
                       - (ZK_e_par-ZK_e_prof) * BigR / BB2 * Bgrad_Te_star * Bgrad_Te      * xjac * tstep &
                       - ZK_e_prof * BigR * (v_x*Te0_x + v_y*Te0_y )                       * xjac * tstep &
                       + zeta * v * r0 * delta_g(8) * BigR                                 * xjac         &
                       - K_perp_numm  *                                                                   &
                         (    (         v_ss   * (x_t**2 + y_t**2)			                  &
                               +        v_tt   * (x_s**2 + y_s**2)			                  &
                               - 2.d0 * v_st   * (x_s*x_t + y_s*y_t) )	                                  &
                            * (         Te0_ss * (x_t**2 + y_t**2)			                  &
                               +        Te0_tt * (x_s**2 + y_s**2)			                  &
                               - 2.d0 * Te0_st * (x_s*x_t + y_s*y_t) ) )	                          &
                         / xjac**4                                                         * xjac * tstep &
                       - TG_num8 * tstep * 0.25d0 * BigR**3 * Te0 * (r0_x * u0_y - r0_y * u0_x)           &
                                 * ( v_x * u0_y - v_y * u0_x)                              * xjac * tstep &
                       - TG_num8 * tstep * 0.25d0 * BigR**3 * r0 * (Te0_x * u0_y - Te0_y * u0_x)          &
                                 * ( v_x * u0_y - v_y * u0_x)                              * xjac * tstep

            rhs_ij_8_k = - (ZK_e_par-ZK_e_prof) * BigR / BB2 * Bgrad_Te_k_star * Bgrad_Te  * xjac * tstep &
                         - ZK_e_prof * BigR * ( v_p*Te0_p /BigR**2 )                       * xjac * tstep &
                         - TG_num8 * tstep * 0.25d0 / BigR * vpar0**2                                     &
                                   * Te0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)               &
                                   * ( v_x * ps0_y -  v_y * ps0_x + F0 / BigR * v_p)       * xjac * tstep &
                         - TG_num8 * tstep * 0.25d0 / BigR * vpar0**2                                     &
                                   * r0 * (Te0_x * ps0_y - Te0_y * ps0_x + F0 / BigR * Te0_p)             &
                                   * ( v_x * ps0_y -  v_y * ps0_x + F0 / BigR * v_p)       * xjac * tstep
    		       
    	    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    	    !!!!!!!!!! RHS equation 9 (Grad(Kpar.Grad(Te)) !!!!!!!!!
    	    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    	    rhs_ij_9   = -(ZK_e_par / BB2* Bgrad_Te_star   * Bgrad_Te + v*KPe0) * BigR * xjac * tstep	   

!    	    rhs_ij_9_k = -(ZK_e_par / BB2* Bgrad_Te_k_star * Bgrad_Te	    )	   * BigR * xjac * tstep   

	    ! --- Fill up the matrix
 	    ij1 = index_ij
 	    ij2 = index_ij + 1*n_tor_loop
 	    ij3 = index_ij + 2*n_tor_loop
 	    ij4 = index_ij + 3*n_tor_loop
 	    ij5 = index_ij + 4*n_tor_loop
 	    ij6 = index_ij + 5*n_tor_loop
 	    ij7 = index_ij + 6*n_tor_loop
 	    ij8 = index_ij + 7*n_tor_loop
! 	    ij9 = index_ij + 8*n_tor_loop

	    if (n_tor .gt. 3) then
 	      RHS_p(mp,ij1) = RHS_p(mp,ij1) + rhs_ij_1   * wst
 	      RHS_p(mp,ij2) = RHS_p(mp,ij2) + rhs_ij_2   * wst
 	      RHS_p(mp,ij3) = RHS_p(mp,ij3) + rhs_ij_3   * wst
 	      RHS_p(mp,ij4) = RHS_p(mp,ij4) + rhs_ij_4   * wst
 	      RHS_p(mp,ij5) = RHS_p(mp,ij5) + rhs_ij_5   * wst
 	      RHS_k(mp,ij5) = RHS_k(mp,ij5) + rhs_ij_5_k * wst
 	      RHS_p(mp,ij6) = RHS_p(mp,ij6) + rhs_ij_6   * wst
 	      RHS_k(mp,ij6) = RHS_k(mp,ij6) + rhs_ij_6_k * wst
 	      RHS_p(mp,ij7) = RHS_p(mp,ij7) + rhs_ij_7   * wst
 	      RHS_p(mp,ij8) = RHS_p(mp,ij8) + rhs_ij_8   * wst
 	      RHS_k(mp,ij8) = RHS_k(mp,ij8) + rhs_ij_8_k * wst
 	      RHS_p(mp,ij9) = RHS_p(mp,ij9) + rhs_ij_9   * wst
! 	      RHS_k(mp,ij9) = RHS_k(mp,ij9) + rhs_ij_9_k * wst
	    else
	      RHS(ij1) = RHS(ij1) + rhs_ij_1 * wst
              RHS(ij2) = RHS(ij2) + rhs_ij_2 * wst
              RHS(ij3) = RHS(ij3) + rhs_ij_3 * wst
              RHS(ij4) = RHS(ij4) + rhs_ij_4 * wst
              RHS(ij5) = RHS(ij5) + rhs_ij_5 * wst + rhs_ij_5_k * wst
              RHS(ij6) = RHS(ij6) + rhs_ij_6 * wst + rhs_ij_6_k * wst
              RHS(ij7) = RHS(ij7) + rhs_ij_7 * wst
              RHS(ij8) = RHS(ij8) + rhs_ij_8 * wst + rhs_ij_8_k * wst
!              RHS(ij9) = RHS(ij9) + rhs_ij_9 * wst + rhs_ij_9_k * wst
	    endif

            ! ---------------------------------
            ! --- Now the LHS (linearised part)
            ! ---------------------------------
 	    do k=1,n_vertex_max

 	      do l=1,n_order+1

                ! --- If we're doing the fft, don't loop...
                n_tor_loop2 = n_tor
                if (n_tor .gt. 3) n_tor_loop2 = 1
	  
                do in = 1, n_tor_loop2

            	  ! --- Build up test function and parallel gradient terms
	          if (n_tor .gt. 3) then
	            HHZ(in,mp)   = 1.d0
	            HHZ_p(in,mp) = 1.d0
 	    	    index_kl = n_var*(n_order+1)*(k-1) + n_var * (l-1) + 1   ! index in the ELM matrix
	          else
	            HHZ(in,mp)   = HZ(in,mp)
	            HHZ_p(in,mp) = HZ_p(in,mp)
                    index_kl = n_tor*n_var*(n_order+1)*(k-1) + n_tor * n_var * (l-1) + in   ! index in the ELM matrix
	          endif

                  psi    = h(k,l,ms,mt)    * element%size(k,l) * HHZ(in,mp)
                  psi_s  = h_s(k,l,ms,mt)  * element%size(k,l) * HHZ(in,mp)
                  psi_t  = h_t(k,l,ms,mt)  * element%size(k,l) * HHZ(in,mp)
                  psi_p  = h(k,l,ms,mt)    * element%size(k,l) * HHZ_p(in,mp)

                  psi_ss = h_ss(k,l,ms,mt) * element%size(k,l) * HHZ(in,mp)
                  psi_tt = h_tt(k,l,ms,mt) * element%size(k,l) * HHZ(in,mp)
                  psi_st = h_st(k,l,ms,mt) * element%size(k,l) * HHZ(in,mp)

                  psi_x  = (   y_t * psi_s - y_s * psi_t ) / xjac
                  psi_y  = ( - x_t * psi_s + x_s * psi_t ) / xjac

                  u    = psi  ; zj    = psi  ; w    = psi   ; rho    = psi   ; Ti    = psi   ; Te    = psi   ; Vpar    = psi   ; KPe = psi
                  u_x  = psi_x; zj_x  = psi_x; w_x  = psi_x ; rho_x  = psi_x ; Ti_x  = psi_x ; Te_x  = psi_x ; Vpar_x  = psi_x
                  u_y  = psi_y; zj_y  = psi_y; w_y  = psi_y ; rho_y  = psi_y ; Ti_y  = psi_y ; Te_y  = psi_y ; Vpar_y  = psi_y
                  u_p  = psi_p; zj_p  = psi_p; w_p  = psi_p ; rho_p  = psi_p ; Ti_p  = psi_p ; Te_p  = psi_p ; Vpar_p  = psi_p
                  u_s  = psi_s; zj_s  = psi_s; w_s  = psi_s ; rho_s  = psi_s ; Ti_s  = psi_s ; Te_s  = psi_s ; Vpar_s  = psi_s
                  u_t  = psi_t; zj_t  = psi_t; w_t  = psi_t ; rho_t  = psi_t ; Ti_t  = psi_t ; Te_t  = psi_t ; Vpar_t  = psi_t
                                               w_ss = psi_ss; rho_ss = psi_ss; Ti_ss = psi_ss; Te_ss = psi_ss; Vpar_ss = psi_ss
                                               w_tt = psi_tt; rho_tt = psi_tt; Ti_tt = psi_tt; Te_tt = psi_tt; Vpar_tt = psi_tt
                                               w_st = psi_st; rho_st = psi_st; Ti_st = psi_st; Te_st = psi_st; Vpar_st = psi_st

                  rho_hat   = BigR**2 * rho
                  rho_x_hat = 2.d0 * BigR * BigR_x  * rho + BigR**2 * rho_x
                  rho_y_hat = BigR**2 * rho_y

                  Bgrad_rho_star_psi = ( v_x   * psi_y - v_y   * psi_x ) / BigR
                  Bgrad_Ti_star_psi  = ( v_x   * psi_y - v_y   * psi_x ) / BigR
                  Bgrad_Te_star_psi  = ( v_x   * psi_y - v_y   * psi_x ) / BigR

                  Bgrad_rho_psi      = ( r0_x  * psi_y - r0_y  * psi_x ) / BigR
                  Bgrad_Ti_psi       = ( Ti0_x * psi_y - Ti0_y * psi_x ) / BigR
                  Bgrad_Te_psi       = ( Te0_x * psi_y - Te0_y * psi_x ) / BigR

                  Bgrad_rho_rho      = ( rho_x * ps0_y - rho_y * ps0_x ) / BigR
                  Bgrad_Ti_Ti        = ( Ti_x  * ps0_y - Ti_y  * ps0_x ) / BigR
                  Bgrad_Te_Te        = ( Te_x  * ps0_y - Te_y  * ps0_x ) / BigR       ! F0 due to absence of normalisation

                  Bgrad_rho_rho_n    = ( F0 / BigR * rho_p ) / BigR
                  Bgrad_Ti_Ti_n      = ( F0 / BigR * Ti_p  ) / BigR
                  Bgrad_Te_Te_n      = ( F0 / BigR * Te_p  ) / BigR

                  BB2_psi            = 2.d0 * (psi_x * ps0_x + psi_y * ps0_y ) /BigR**2

                  ! --- Bootstrap current coefficients (Wesson formula)
                  if (bootstrap) then
                    call bootstrap_current_lhs(BigR, minRad, R_axis, &
                                               psi_axis, psi_bnd,    &
                                               ps0, ps0_x, ps0_y,    &
                                               psi, psi_x, psi_y,    &
                                               r0,  r0_x,  r0_y,     &
                                               rho, rho_x, rho_y,    &
                                               Ti0, Ti0_x, Ti0_y,    &
                                               Ti,  Ti_x,  Ti_y,     &
                                               Te0, Te0_x, Te0_y,    &
                                               Te,  Te_x,  Te_y,     &
                                               dJb_psi, dJb_rho, dJb_Ti, dJb_Te)
                  else
                    dJb_psi = 0.d0
                    dJb_rho = 0.d0
                    dJb_Ti  = 0.d0
                    dJb_Te  = 0.d0
                  endif

    		  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    		  !!!!!!!!!! LHS equation 1 (psi - induction) !!!!!!!!!!!!
    		  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    		  amat_11 = + v * psi / BigR					       * xjac * (1.d0+zeta)   &
    			    - v * (psi_x * u0_y - psi_y * u0_x) 		       * xjac * theta * tstep &
                            - v * eta_Te * dJb_psi / BigR                              * xjac * theta * tstep

    		  amat_12 = -  v * (ps0_x * u_y - ps0_y * u_x)  		       * xjac * theta * tstep

    		  amat_12_n = +  eps_cyl * F0 / BigR * v * u_p  		       * xjac * theta * tstep

    		  amat_13 = - eta_numm * (v_x * zj_x + v_y * zj_y)		       * xjac * theta * tstep &
    			    - eta_Te * v * zj / BigR				       * xjac * theta * tstep

                  amat_15 = - v * eta_Te * dJb_rho / BigR                              * xjac * theta * tstep

                  amat_16 = - v * eta_Te * dJb_Ti  / BigR                              * xjac * theta * tstep

    		  amat_18 = - deta_dTe * v * Te * (zj0 - current_source + Jb) / BigR   * xjac * theta * tstep &
                            - v * eta_Te * dJb_Te  / BigR                              * xjac * theta * tstep


    		  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    		  !!!!!!!!!! LHS equation 2 (U - momentum) !!!!!!!!!!!!
    		  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    		  amat_21 = - v * (psi_x * zj0_y - psi_y * zj0_x )				      * xjac * theta * tstep

    		  amat_22 = - BigR**3 * abs(r0) * (v_x * u_x + v_y * u_y)			      * xjac * (1.d0 + zeta) &
    			    + r0_hat * BigR**2 * w0 * (v_x * u_y  - v_y  * u_x) 		      * xjac * theta * tstep &
    			    + BigR**2 * (u_x * u0_x + u_y * u0_y) * (v_x * r0_y_hat - v_y * r0_x_hat) * xjac * theta * tstep &
    			    + TG_num2 * tstep * 0.25d0 * r0_hat * BigR**3 * (w0_x * u_y - w0_y * u_x)			     &
    				      * ( v_x * u0_y - v_y * u0_x)				      * xjac * theta * tstep &
    			    + TG_num2 * tstep * 0.25d0 * r0_hat * BigR**3 * (w0_x * u0_y - w0_y * u0_x) 		     &
    				      * ( v_x * u_y - v_y * u_x)				      * xjac * theta * tstep

    		  amat_23 = - v * (ps0_x * zj_y  - ps0_y * zj_x)				      * xjac * theta * tstep

    		  amat_23_n = + eps_cyl * F0 / BigR * v * zj_p  				      * xjac * theta * tstep

    		  amat_24 = r0_hat * BigR**2 * w  * ( v_x * u0_y - v_y * u0_x)  		      * xjac * theta * tstep &
    			  + BigR * ( v_x * w_x + v_y * w_y) * visco_Te  			      * xjac * theta * tstep &
    			  + visco_numm  *										     &
    			    (	 (	   v_ss * (x_t**2 + y_t**2)							     & 
    				  +	   v_tt * (x_s**2 + y_s**2)							     &
    				  - 2.d0 * v_st * (x_s*x_t + y_s*y_t) ) 						     &
    			       * (	   w_ss * (x_t**2 + y_t**2)							     &
    				  +	   w_tt * (x_s**2 + y_s**2)							     &
    				  - 2.d0 * w_st * (x_s*x_t + y_s*y_t) ) )						     &
    			    / xjac**4								      * xjac * theta * tstep &
    			  + TG_num2 * tstep * 0.25d0 * r0_hat * BigR**3 * (w_x * u0_y - w_y * u0_x)			     &
    				    * ( v_x * u0_y - v_y * u0_x)				      * xjac * theta * tstep

    		  amat_25 = + 0.5d0 * vv2 * (v_x * rho_y_hat - v_y * rho_x_hat) 		      * xjac * theta * tstep &
    			    + rho_hat * BigR**2 * w0 * (v_x * u0_y - v_y * u0_x)		      * xjac * theta * tstep &
    			    - BigR**2 * (v_x * rho_y * (Ti0 + Te0)   - v_y * rho_x * (Ti0 + Te0)  )   * xjac * theta * tstep &
    			    - BigR**2 * (v_x * rho   * (Ti0_y + Te0_y) - v_y * rho * (Ti0_x + Te0_x) )* xjac * theta * tstep &
    			    + TG_num2 * tstep * 0.25d0 * rho_hat * BigR**3 * (w0_x * u0_y - w0_y * u0_x)		     &
    				      * ( v_x * u0_y - v_y * u0_x)				      * xjac * theta * tstep

    		  amat_26 = - BigR**2 * (v_x * r0_y * Ti   - v_y * r0_x * Ti)			      * xjac * theta * tstep &
    			    - BigR**2 * (v_x * r0   * Ti_y - v_y * r0	* Ti_x) 		      * xjac * theta * tstep  

    		  amat_28 = - BigR**2 * (v_x * r0_y * Te   - v_y * r0_x * Te)			      * xjac * theta * tstep &
    			    - BigR**2 * (v_x * r0   * Te_y - v_y * r0	* Te_x) 		      * xjac * theta * tstep &
    			    + dvisco_dTe * Te * ( v_x * w0_x + v_y * w0_y ) * BigR		      * xjac * theta * tstep

    		  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    		  !!!!!!!!!! LHS equation 3 (j - current) !!!!!!!!!!!!
    		  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    		  amat_31 = (v_x * psi_x + v_y * psi_y ) / BigR * xjac * tstep

    		  amat_33 = v * zj / BigR			* xjac * tstep 

    		  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    		  !!!!!!!!!! LHS equation 4 (w - vorticity) !!!!!!!!!!
    		  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    		  amat_42 = (v_x * u_x + v_y * u_y) * BigR * xjac * tstep

    		  amat_44 =  v * w * BigR		   * xjac * tstep 
    		  
    		  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    		  !!!!!!!!!! LHS equation 5 (rho - continuity) !!!!!!!!!
    		  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    		  amat_51 = - (D_par-D_prof) * BigR * BB2_psi/ BB2**2 * Bgrad_rho_star     * Bgrad_rho     * xjac * theta * tstep &
    			    + (D_par-D_prof) * BigR / BB2	      * Bgrad_rho_star_psi * Bgrad_rho     * xjac * theta * tstep &
    			    + (D_par-D_prof) * BigR / BB2	      * Bgrad_rho_star     * Bgrad_rho_psi * xjac * theta * tstep &
    			    + v * Vpar0 * (r0_x * psi_y - r0_y * psi_x) 				   * xjac * theta * tstep &
    			    + v * r0 * (vpar0_x * psi_y - vpar0_y * psi_x)				   * xjac * theta * tstep &
    			    + TG_num5 * tstep * 0.25d0 / BigR * vpar0**2							  &
    				      * (r0_x * psi_y - r0_y * psi_x)								  &
    				      * ( v_x * ps0_y -  v_y * ps0_x )  				   * xjac * theta * tstep &
    			    + TG_num5 * tstep * 0.25d0 / BigR * vpar0**2							  &
    				      * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)					  &
    				      * ( v_x * psi_y -  v_y * psi_x )  				   * xjac * theta * tstep

    		  amat_51_k = - (D_par-D_prof) * BigR * BB2_psi/ BB2**2 * Bgrad_rho_k_star * Bgrad_rho     * xjac * theta * tstep &
    			      + (D_par-D_prof) * BigR / BB2		* Bgrad_rho_k_star * Bgrad_rho_psi * xjac * theta * tstep &
    			      + TG_num5 * tstep * 0.25d0 / BigR * vpar0**2							  &
    					* (r0_x * psi_y - r0_y * psi_x) 							  &
    					* (			       + F0 / BigR * v_p)		   * xjac * theta * tstep

    		  amat_52 =  - v * BigR**2 * ( r0_x * u_y - r0_y * u_x) 				   * xjac * theta * tstep &
    			     - v * 2.d0 * BigR * r0 * u_y						   * xjac * theta * tstep &
    			     + TG_num5 * tstep * 0.25d0 * BigR**3 * (r0_x * u_y  - r0_y * u_x)  				  &
    				       * ( v_x * u0_y - v_y  * u0_x)					   * xjac * theta * tstep &
    			     + TG_num5 * tstep * 0.25d0 * BigR**3 * (r0_x * u0_y - r0_y * u0_x) 				  &
    				       * ( v_x * u_y  - v_y  * u_x)					   * xjac * theta * tstep

    		  amat_55 = v * rho * BigR								   * xjac * (1.d0 + zeta) &
    			    - v * BigR**2 * ( rho_x * u0_y - rho_y * u0_x)				   * xjac * theta * tstep &
    			    - v * 2.d0 * BigR * rho * u0_y						   * xjac * theta * tstep &
    			    + (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star * Bgrad_rho_rho		   * xjac * theta * tstep &
    			    + D_prof * BigR  * (v_x*rho_x + v_y*rho_y ) 				   * xjac * theta * tstep &
    			    + v * Vpar0 * (rho_x * ps0_y - rho_y * ps0_x)				   * xjac * theta * tstep &
    			    + v * rho * (vpar0_x * ps0_y - vpar0_y * ps0_x)				   * xjac * theta * tstep &
    			    + v * rho * F0 / BigR * vpar0_p						   * xjac * theta * tstep &
    			    + D_perp_num  *											  &
    			      (    (	     v_ss   * (x_t**2 + y_t**2) 							  &
    				    +	     v_tt   * (x_s**2 + y_s**2) 							  &
    				    - 2.d0 * v_st   * (x_s*x_t + y_s*y_t) )							  &
    				 * (	     rho_ss * (x_t**2 + y_t**2) 							  &
    				    +	     rho_tt * (x_s**2 + y_s**2) 							  &
    				    - 2.d0 * rho_st * (x_s*x_t + y_s*y_t) ) )							  &
    			      / xjac**4 								   * xjac * theta * tstep &
    			    + TG_num5 * tstep * 0.25d0 * BigR**3 * (rho_x * u0_y - rho_y * u0_x)				  &
    				      * ( v_x  * u0_y - v_y   * u0_x )  				   * xjac * theta * tstep &
    			    + TG_num5 * tstep * 0.25d0 / BigR * vpar0**2							  &
    				      * (rho_x * ps0_y - rho_y * ps0_x )							  &
    				      * ( v_x * ps0_y -  v_y * ps0_x )  				   * xjac * theta * tstep

    		  amat_55_k = + (D_par-D_prof) * BigR / BB2 * Bgrad_rho_k_star * Bgrad_rho_rho  	   * xjac * theta * tstep &
    			      + TG_num5 * tstep * 0.25d0 / BigR * vpar0**2							  &
    					* (rho_x * ps0_y - rho_y * ps0_x )							  &
    					* (				 + F0 / BigR * v_p)		   * xjac * theta * tstep

    		  amat_55_n = + (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star   * Bgrad_rho_rho_n	   * xjac * theta * tstep &
    			      + v * F0 / BigR * Vpar0 * rho_p						   * xjac * theta * tstep &
    			      + TG_num5 * tstep * 0.25d0 / BigR * vpar0**2							  &
    					* (				 + F0 / BigR * rho_p)					  &
    					* ( v_x * ps0_y -  v_y * ps0_x )				   * xjac * theta * tstep

    		  amat_55_kn = + (D_par-D_prof) * BigR / BB2 * Bgrad_rho_k_star * Bgrad_rho_rho_n	   * xjac * theta * tstep &
    			       + D_prof * BigR  * ( v_p*rho_p * eps_cyl**2 /BigR**2 )			   * xjac * theta * tstep &
    			       + TG_num5 * tstep * 0.25d0 / BigR * vpar0**2							  &
    					 * (				  + F0 / BigR * rho_p)  				  &
    					 * (				  + F0 / BigR * v_p)		   * xjac * theta * tstep

    		  amat_57   = + v * F0 / BigR * Vpar * r0_p						   * xjac * theta * tstep &
    			      + v * Vpar * (r0_x * ps0_y - r0_y * ps0_x)				   * xjac * theta * tstep &
    			      + v * r0 * (vpar_x * ps0_y - vpar_y * ps0_x)				   * xjac * theta * tstep &
    			      + TG_num5 * tstep * 0.25d0 / BigR * 2.d0*vpar0*vpar						  &
    					* (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)					  &
    					* ( v_x * ps0_y -  v_y * ps0_x  		 )		   * xjac * theta * tstep 

    		  amat_57_k = + TG_num5 * tstep * 0.25d0 / BigR * 2.d0*vpar0*vpar						  &
    					* (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)					  &
    					* (			       + F0 / BigR * v_p)		   * xjac * theta * tstep 

    		  amat_57_n = + v * r0 * F0 / BigR * vpar_p						   * xjac * theta * tstep
    		  
    		  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    		  !!!!!!!!!! LHS equation 6 (Ti - Ion energy) !!!!!!!!!
    		  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    		  amat_61 = + v * r0 * Vpar0 * (Ti0_x * psi_y - Ti0_y * psi_x)  				* xjac * theta * tstep &
    			    + v * (GAMMA-1.d0) * r0 * Ti0 * (vpar0_x * psi_y - vpar0_y * psi_x) 		* xjac * theta * tstep &
    			    - (ZK_i_par-ZK_i_prof) * BigR * BB2_psi / BB2**2 * Bgrad_Ti_star * Bgrad_Ti 	* xjac * theta * tstep &
    			    + (ZK_i_par-ZK_i_prof) * BigR / BB2 	     * Bgrad_Ti_star_psi * Bgrad_Ti	* xjac * theta * tstep &
    			    + (ZK_i_par-ZK_i_prof) * BigR / BB2 	     * Bgrad_Ti_star * Bgrad_Ti_psi	* xjac * theta * tstep & 
    			    + TG_num6 * tstep * 0.25d0 / BigR * vpar0**2							       &
    				      * Ti0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)  				       &
    				      * ( v_x * psi_y -  v_y * psi_x )  					* xjac * theta * tstep &
    			    + TG_num6 * tstep * 0.25d0 / BigR * vpar0**2							       &
    				      * r0 * (Ti0_x * ps0_y - Ti0_y * ps0_x + F0 / BigR * Ti0_p)				       &
    				      * ( v_x * psi_y -  v_y * ps0_x )  					* xjac * theta * tstep &
    			    + TG_num6 * tstep * 0.25d0 / BigR * vpar0**2							       &
    				      * Ti0 * (r0_x * psi_y - r0_y * psi_x)							       &
    				      * ( v_x * ps0_y -  v_y * ps0_x )  					* xjac * theta * tstep &
    			    + TG_num6 * tstep * 0.25d0 / BigR * vpar0**2							       &
    				      * r0 * (Ti0_x * psi_y - Ti0_y * psi_x)							       &
    				      * ( v_x * ps0_y -  v_y * ps0_x )  					* xjac * theta * tstep  	     

    		  amat_61_k = - (ZK_i_par-ZK_i_prof) * BigR * BB2_psi / BB2**2 * Bgrad_Ti_k_star * Bgrad_Ti	* xjac * theta * tstep &
    			      + (ZK_i_par-ZK_i_prof) * BigR / BB2	       * Bgrad_Ti_k_star * Bgrad_Ti_psi * xjac * theta * tstep &
    			      + TG_num6 * tstep * 0.25d0 / BigR * vpar0**2							       &
    					* Ti0 * (r0_x * psi_y - r0_y * psi_x)							       &
    					* ( F0 / BigR * v_p)							* xjac * theta * tstep &
    			      + TG_num6 * tstep * 0.25d0 / BigR * vpar0**2							       &
    					* r0 * (Ti0_x * psi_y - Ti0_y * psi_x)  						       &
    					* ( F0 / BigR * v_p)							* xjac * theta * tstep  	  

    		  amat_62 = - v * r0 * BigR**2 * ( Ti0_x * u_y - Ti0_y * u_x)					* xjac * theta * tstep &
    			    - v * 2.d0 * (GAMMA-1.d0) * r0 * BigR * Ti0 * u_y					* xjac * theta * tstep &
    			    + TG_num6 * tstep * 0.25d0 * BigR**2 * Ti0* (r0_x * u_y - r0_y * u_x)				       &
    				      * ( v_x * u0_y - v_y * u0_x)						* xjac * theta * tstep &
    			    + TG_num6 * tstep * 0.25d0 * BigR**2 * r0* (Ti0_x * u_y - Ti0_y * u_x)				       &
    				      * ( v_x * u0_y - v_y * u0_x)						* xjac * theta * tstep &
    			    + TG_num6 * tstep * 0.25d0 * BigR**2 * Ti0* (r0_x * u0_y - r0_y * u0_x)				       &
    				      * ( v_x * u_y - v_y * u_x)						* xjac * theta * tstep &
    			    + TG_num6 * tstep * 0.25d0 * BigR**2 * r0* (Ti0_x * u0_y - Ti0_y * u0_x)				       &
    				      * ( v_x * u_y - v_y * u_x)						* xjac * theta * tstep

    		  amat_65 = - v * rho * BigR**2 * (Ti0_x * u0_y - Ti0_y * u0_x) 				* xjac * theta * tstep &
    			    + v * rho * Vpar0 * F0/BigR * Ti0_p 						* xjac * theta * tstep &
    			    + v * rho * Vpar0 * (Ti0_x * ps0_y - Ti0_y * ps0_x) 				* xjac * theta * tstep &
    			    - v * 2.d0 * (GAMMA-1.d0) * rho * BigR * Ti0 * u0_y 				* xjac * theta * tstep &
    			    + v * (GAMMA-1.d0) * rho * Ti0 * F0/BigR * Vpar0_p  				* xjac * theta * tstep &
    			    + v * (GAMMA-1.d0) * rho * Ti0 * (Vpar0_x * ps0_y - Vpar0_y * ps0_x)		* xjac * theta * tstep &
    			    ! Temperature transfer (not yet tested!)
    			    !- v * BigR * dQ_drho * rho * (Te0 - Ti0)						* xjac * theta * tstep & 
    			    + TG_num6 * tstep * 0.25d0 * BigR**2 * Ti0* (rho_x * u0_y - rho_y * u0_x)				       &
    				      * ( v_x * u0_y - v_y * u0_x)						* xjac * theta * tstep &
    			    + TG_num6 * tstep * 0.25d0 * BigR**2 * rho * (Ti0_x * u0_y - Ti0_y * u0_x)  			       &
    				      * ( v_x * u0_y - v_y * u0_x)						* xjac * theta * tstep &
    			    + TG_num6 * tstep * 0.25d0 / BigR * vpar0**2							       &
    				      * Ti0 * (rho_x * ps0_y - rho_y * ps0_x + F0 / BigR * rho_p)				       &
    				      * ( v_x * ps0_y -  v_y * ps0_x )  					* xjac * theta * tstep &
    			    + TG_num6 * tstep * 0.25d0 / BigR * vpar0**2							       &
    				      * rho * (Ti0_x * ps0_y - Ti0_y * ps0_x + F0 / BigR * Ti0_p)				       &
    				      * ( v_x * ps0_y -  v_y * ps0_x )  					* xjac * theta * tstep

    		  amat_65_k = + TG_num6 * tstep * 0.25d0 / BigR * vpar0**2							       &
    					* Ti0 * (rho_x * ps0_y - rho_y * ps0_x + F0 / BigR * rho_p)				       &
    					* ( F0 / BigR * v_p)							* xjac * theta * tstep &
    			      + TG_num6 * tstep * 0.25d0 / BigR * vpar0**2							       &
    					* rho * (Ti0_x * ps0_y - Ti0_y * ps0_x + F0 / BigR * Ti0_p)				       &
    					* ( F0 / BigR * v_p)							* xjac * theta * tstep  

    		  amat_66 = + v * abs(r0) * Ti * BigR								* xjac * (1.d0 + zeta) &
    			    - v * r0 * BigR**2 * ( Ti_x * u0_y - Ti_y * u0_x)					* xjac * theta * tstep &
    			    + v * r0 * Vpar0 * (Ti_x * ps0_y - Ti_y * ps0_x)					* xjac * theta * tstep &
    			    - 2.d0 * v * r0 * (GAMMA-1.d0) * Ti * BigR * u0_y					* xjac * theta * tstep &
    			    + v * r0 * (GAMMA-1.d0) * Ti * F0/BigR * Vpar0_p					* xjac * theta * tstep &
    			    + v * r0 * (GAMMA-1.d0) * Ti * (Vpar0_x * ps0_y - Vpar0_y * ps0_x)  		* xjac * theta * tstep &
    			    + (ZK_i_par-ZK_i_prof) * BigR / BB2 * Bgrad_Ti_star * Bgrad_Ti_Ti			* xjac * theta * tstep &
    			    + dZK_i_par * Ti * BigR / BB2 * Bgrad_Ti_star * Bgrad_Ti				* xjac * theta * tstep &
    			    + ZK_i_prof * BigR * (v_x*Ti_x + v_y*Ti_y ) 					* xjac * theta * tstep & 
    			    + K_perp_numm  *											       &
    			      (    (	     v_ss  * (x_t**2 + y_t**2)  							       &
    				    +	     v_tt  * (x_s**2 + y_s**2)  							       &
    				    - 2.d0 * v_st  * (x_s*x_t + y_s*y_t) )							       &
    				 * (	     Ti_ss * (x_t**2 + y_t**2)  							       &
    				    +	     Ti_tt * (x_s**2 + y_s**2)  							       &
    				    - 2.d0 * Ti_st * (x_s*x_t + y_s*y_t) ) )	 / xjac**4			* xjac * theta * tstep &
    			    ! Temperature transfer (not yet tested!)
    			    !+ v * BigR * Q * Ti								* xjac * theta * tstep &
    			    + TG_num6 * tstep * 0.25d0 * BigR**2 * Ti * (r0_x * u0_y - r0_y * u0_x)				       &
    				      * ( v_x * u0_y - v_y * u0_x)						* xjac * theta * tstep &
    			    + TG_num6 * tstep * 0.25d0 * BigR**2 * r0* (Ti_x * u0_y - Ti_y * u0_x)				       &
    				      * ( v_x * u0_y - v_y * u0_x)						* xjac * theta * tstep &	    
    			    + TG_num6 * tstep * 0.25d0 / BigR * vpar0**2							       &
    				      * Ti * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)					       &
    				      * ( v_x * ps0_y -  v_y * ps0_x )  					* xjac * theta * tstep &
    			    + TG_num6 * tstep * 0.25d0 / BigR * vpar0**2							       &
    				      * r0 * (Ti_x * ps0_y - Ti_y * ps0_x )							       &
    				      * ( v_x * ps0_y -  v_y * ps0_x )  					* xjac * theta * tstep

    		  amat_66_k = + (ZK_i_par-ZK_i_prof) * BigR / BB2 * Bgrad_Ti_k_star * Bgrad_Ti_Ti		* xjac * theta * tstep &
    			      + dZK_i_par * Ti * BigR / BB2 * Bgrad_Ti_k_star * Bgrad_Ti			* xjac * theta * tstep &
    			      + TG_num6 * tstep * 0.25d0 / BigR * vpar0**2							       &
    					* Ti * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p) 				       &
    					* ( F0 / BigR * v_p)							* xjac * theta * tstep &
    			      + TG_num6 * tstep * 0.25d0 / BigR * vpar0**2							       &
    					* r0 * (Ti_x * ps0_y - Ti_y * ps0_x )							       &
    					* ( F0 / BigR * v_p)							* xjac * theta * tstep
    			      
    		  amat_66_n = + (ZK_i_par-ZK_i_prof) * BigR / BB2 * Bgrad_Ti_star   * Bgrad_Ti_Ti_n		* xjac * theta * tstep &
    			      + v * r0 * Vpar0 * F0/BigR * Ti_p 						* xjac * theta * tstep &
    			      + TG_num6 * tstep * 0.25d0 / BigR * vpar0**2							       &
    					* r0 * ( F0 / BigR * Ti_p)								       &
    					* ( v_x * ps0_y -  v_y * ps0_x )					* xjac * theta * tstep

    		  amat_66_kn = + (ZK_i_par-ZK_i_prof) * BigR / BB2 * Bgrad_Ti_k_star * Bgrad_Ti_Ti_n		* xjac * theta * tstep &
    			       + ZK_i_prof * BigR   * (v_p*Ti_p /BigR**2 )					* xjac * theta * tstep &
    			       + TG_num6 * tstep * 0.25d0 / BigR * vpar0**2							       &
    					 * r0 * ( F0 / BigR * Ti_p)								       &
    					 * ( F0 / BigR * v_p)							* xjac * theta * tstep  

    		  amat_67 = + v * r0 * F0/BigR * Vpar * Ti0_p							* xjac * theta * tstep &
    			    + v * r0 * Vpar * (Ti0_x * ps0_y - Ti0_y * ps0_x)					* xjac * theta * tstep &
    			    + v * (GAMMA-1.d0) * r0 * Ti0 * (Vpar_x * ps0_y - Vpar_y * ps0_x)			* xjac * theta * tstep &
    			    + TG_num6 * tstep * 0.25d0 / BigR * 2.d0 * vpar0*vpar						       &
    				      * Ti0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)  				       &
    				      * ( v_x * ps0_y -  v_y * ps0_x )  					* xjac * theta * tstep &
    			    + TG_num6 * tstep * 0.25d0 / BigR * 2.d0 * vpar0*vpar						       &
    				      * r0 * (Ti0_x * ps0_y - Ti0_y * ps0_x + F0 / BigR * Ti0_p)				       &
    				      * ( v_x * ps0_y -  v_y * ps0_x )  					* xjac * theta * tstep
    		     
    		  amat_67_k = + TG_num6 * tstep * 0.25d0 / BigR * 2.d0 * vpar0*vpar						       &
    					* Ti0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)				       &
    					* ( F0 / BigR * v_p)							* xjac * theta * tstep &
    			      + TG_num6 * tstep * 0.25d0 / BigR * 2.d0 * vpar0*vpar						       &
    					* r0 * (Ti0_x * ps0_y - Ti0_y * ps0_x + F0 / BigR * Ti0_p)				       &
    					* ( F0 / BigR * v_p)							* xjac * theta * tstep

    		  amat_67_n = v * (GAMMA-1.d0) * r0 * Ti0 * F0/BigR * Vpar_p					* xjac * theta * tstep       
    		 
    		  ! Temperature transfer (not yet tested!)
    		  amat_68 = 0.d0!- v * BigR * dQ_dTe * Te * (Te0 - Ti0) 					* xjac * theta * tstep & 
    				!- v * BigR * Q * Te								* xjac * theta * tstep  
    		  
    		  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    		  !!!!!!!!!! LHS equation 7 (Vpar - parallel momentum) !!!!!!!!!
    		  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    		  amat_71 = + v * (P0_x * psi_y - P0_y * psi_x) 				   * xjac * theta * tstep &
    			    !+ vpar0 * (F0/BigR)**2 * (vpar0_x * ps_y - vpar0_y * ps_x) 	   * xjac * theta * tstep &
    			    + TG_NUM7 * tstep * 0.25d0 * r0 * Vpar0**2 * BB2						  &
    				      * (-(psi_x * vpar0_y - psi_y * vpar0_x)) / BigR					  &
    				      * (-(psi_x * v_y     - psi_y * v_x)    )  		   * xjac * theta * tstep &
    			    + TG_NUM7 * tstep * 0.25d0 * v  * Vpar0**2 * BB2						  &
    				      * (-(psi_x * vpar0_y - psi_y * vpar0_x)) / BigR					  &
    				      * (-(psi_x * r0_y    - psi_y * r0_x)   )  		   * xjac * theta * tstep 

    		  amat_72 = 0.d0 ! &
    				 !+ F0 * (u_s * vpar0_t - u_t * vpar0_s) * theta * tstep 

    		  amat_75 = + v * (rho_x * (Ti0 + Te0) * ps0_y - rho_y * (Ti0 + Te0) * ps0_x)	   * xjac * theta * tstep &
    			    + v * (rho * (Ti0_x + Te0_x) * ps0_y - rho * (Ti0_y + Te0_y) * ps0_x)  * xjac * theta * tstep &
    			    + v * F0 / BigR * rho * (Ti0_p + Te0_p)				   * xjac * theta * tstep &
    			    + TG_NUM7 * tstep * 0.25d0 * rho * Vpar0**2 * BB2						  &
    				      * (-(ps0_x * vpar0_y - ps0_y * vpar0_x) + F0 / BigR * vpar0_p) / BigR		  &
    				      * (-(ps0_x * v_y     - ps0_y * v_x)     ) 		   * xjac * theta * tstep &
    			    + TG_NUM7 * tstep * 0.25d0 * v * Vpar0**2 * BB2						  &
    				      * (-(ps0_x * vpar0_y - ps0_y * vpar0_x) + F0 / BigR * vpar0_p) / BigR		  &
    				      * (-(ps0_x * rho_y   - ps0_y * rho_x)   ) 		   * xjac * theta * tstep

    		  amat_75_k = + TG_NUM7 * tstep * 0.25d0 * rho * Vpar0**2 * BB2 					  &
    					* (-(ps0_x * vpar0_y - ps0_y * vpar0_x) + F0 / BigR * vpar0_p) / BigR		  &
    					* ( F0 / BigR * v_p )					   * xjac * theta * tstep

    		  amat_75_n = v * F0 / BigR * rho_p * (Ti0 + Te0)				   * xjac * theta * tstep & 
    			      + TG_NUM7 * tstep * 0.25d0 * v * Vpar0**2 * BB2						  &
    					* (-(ps0_x * vpar0_y - ps0_y * vpar0_x) + F0 / BigR * vpar0_p) / BigR		  &
    					* ( F0 / BigR * rho_p ) 				   * xjac * theta * tstep

    		  amat_76 = + v * (Ti_x * r0 * ps0_y - Ti_y * r0 * ps0_x)			   * xjac * theta * tstep &
    			    + v * (Ti * r0_x * ps0_y   - Ti * r0_y * ps0_x)			   * xjac * theta * tstep &
    			    + v * F0 / BigR * Ti * r0_p 					   * xjac * theta * tstep 
    		  
    		  amat_76_n = v * F0 / BigR * Ti_p * r0 					   * xjac * theta * tstep

    		  amat_77 = v * Vpar * abs(r0) * F0**2 / BigR					   * xjac * (1.d0 + zeta) &
    			    + visco_par * (v_x * Vpar_x + v_y * Vpar_y) * BigR  		   * xjac * theta * tstep &
    			    + visco_par_numm  * 									  &
    			      (    (	    v_ss    * (x_t**2 + y_t**2) 						  &
    				    +	    v_tt    * (x_s**2 + y_s**2) 						  &
    				    -2.d0 * v_st    * (x_s*x_t + y_s*y_t) )						  &
    				 * (	    Vpar_ss * (x_t**2 + y_t**2) 						  &
    				    +	    Vpar_tt * (x_s**2 + y_s**2) 						  &
    				    -2.d0 * Vpar_st * (x_s*x_t + y_s*y_t) ) )						  &
    			      / xjac**4 							   * xjac * theta * tstep &
    			    !+ F0 * (u0_x * vpar_y - u0_y * vpar_x)				   * xjac * theta * tstep &
    			    !+ vpar * (F0/BigR)**2  * (vpar0_x * ps0_y - vpar0_y * ps0_x)	   * xjac * theta * tstep &
    			    !+ vpar0 * (F0/BigR)**2 * (vpar_x  * ps0_y - vpar_y * ps0_x)	   * xjac * theta * tstep &
    			    !+ (F0/BigR)**3 * vpar  * vpar0_p					   * xjac * theta * tstep &
    			    !+ (F0/BigR)**3 * vpar0 * vpar_p					   * xjac * theta * tstep &
    			    + TG_NUM7 * tstep * 0.5d0 * r0 * Vpar * Vpar0 * BB2 					  &
    				      * (-(ps0_x * vpar0_y - ps0_y * vpar0_x) + F0 / BigR * vpar0_p) / BigR		  &
    				      * (-(ps0_x * v_y     - ps0_y * v_x)    )  		   * xjac * theta * tstep &
    			    + TG_NUM7 * tstep * 0.5d0 * v * Vpar * Vpar0 * BB2  					  &
    				      * (-(ps0_x * vpar0_y - ps0_y * vpar0_x) + F0 / BigR * vpar0_p) / BigR		  &
    				      * (-(ps0_x * r0_y    - ps0_y * r0_x)    + F0 / BigR * r0_p)  * xjac * theta * tstep &
    			    + TG_NUM7 * tstep * 0.25d0 * r0 * Vpar0**2 * BB2						  &
    				      * (-(ps0_x * vpar_y  - ps0_y * vpar_x) ) / BigR					  &
    				      * (-(ps0_x * v_y     - ps0_y * v_x)  )			   * xjac * theta * tstep &
    			    + TG_NUM7 * tstep * 0.25d0 * v * Vpar0**2 * BB2						  &
    				      * (-(ps0_x * vpar_y  - ps0_y * vpar_x) ) / BigR					  &
    				      * (-(ps0_x * r0_y    - ps0_y * r0_x)    + F0 / BigR * r0_p)  * xjac * theta * tstep

    		  amat_77_k = + TG_NUM7 * tstep * 0.5d0 * r0 * Vpar * Vpar0 * BB2					  &
    					* (-(ps0_x * vpar0_y - ps0_y * vpar0_x) + F0 / BigR * vpar0_p) / BigR		  &
    					* ( F0 / BigR * v_p )					   * xjac * theta * tstep &
    			      + TG_NUM7 * tstep * 0.25d0 * r0 * Vpar0**2 * BB2  					  &
    					* (-(ps0_x * vpar_y - ps0_y * vpar_x) ) / BigR  				  &
    					* ( F0 / BigR * v_p)					   * xjac * theta * tstep 

    		  amat_77_n = + TG_NUM7 * tstep * 0.25d0 * r0 * Vpar0**2 * BB2  					  &
    					* ( F0 / BigR * vpar_p ) / BigR 						  &
    					* (-(ps0_x * v_y - ps0_y * v_x) )			   * xjac * theta * tstep &
    			      + TG_NUM7 * tstep * 0.25d0 * v * Vpar0**2 * BB2						  &
    					* ( F0 / BigR * vpar_p ) / BigR 						  &
    					* (-(ps0_x * r0_y - ps0_y * r0_x) + F0 / BigR * r0_p)	   * xjac * theta * tstep

    		  amat_77_kn = + TG_NUM7 * tstep * 0.25d0 * r0 * Vpar0**2 * BB2 					  &
    					 * ( F0 / BigR * vpar_p) / BigR 						  &
    					 * ( F0 / BigR * v_p )  				   * xjac * theta * tstep    

    		  amat_78 = + v * (Te_x * r0 * ps0_y   - Te_y * r0 * ps0_x)			   * xjac * theta * tstep &
    			    + v * (Te * r0_x * ps0_y   - Te * r0_y * ps0_x)			   * xjac * theta * tstep &
    			    + v * F0 / BigR * Te * r0_p 					   * xjac * theta * tstep 
    			    
    		  amat_78_n = v * F0 / BigR * Te_p * r0 					   * xjac * theta * tstep

    		  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    		  !!!!!!!!!! LHS equation 8 (Te - electron energy) !!!!!!!!!
    		  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    		  amat_81 = + v * r0 * Vpar0 * (Te0_x * psi_y - Te0_y * psi_x)  				* xjac * theta * tstep &
    			    + v * (GAMMA-1.d0) * r0 * Te0 * (Vpar0_x * psi_y - vpar0_y * psi_x) 		* xjac * theta * tstep &
    			    - (ZK_e_par-ZK_e_prof) * BigR * BB2_psi / BB2**2 * Bgrad_Te_star * Bgrad_Te			* xjac * theta * tstep &
    			    + (ZK_e_par-ZK_e_prof) * BigR / BB2 	* Bgrad_Te_star_psi	 * Bgrad_Te		       * xjac * theta * tstep &
    			    + (ZK_e_par-ZK_e_prof) * BigR / BB2 	* Bgrad_Te_star 	 * Bgrad_Te_psi 	       * xjac * theta * tstep &
    			    !- dZK_e_prof * psi * BigR / BB2 * Bgrad_Te_star * Bgrad_Te 			* xjac * theta * tstep & 
    			    + TG_num8 * tstep * 0.25d0 / BigR * vpar0**2							       &
    				      * Te0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)  				       &
    				      * ( v_x * psi_y -  v_y * psi_x )  					* xjac * theta * tstep &
    			    + TG_num8 * tstep * 0.25d0 / BigR * vpar0**2							       &
    				      * r0 * (Te0_x * ps0_y - Te0_y * ps0_x + F0 / BigR * Te0_p)				       &
    				      * ( v_x * psi_y -  v_y * ps0_x )  					* xjac * theta * tstep &
    			    + TG_num8 * tstep * 0.25d0 / BigR * vpar0**2							       &
    				      * Te0 * (r0_x * psi_y - r0_y * psi_x)							       &
    				      * ( v_x * ps0_y -  v_y * ps0_x )  					* xjac * theta * tstep &
    			    + TG_num8 * tstep * 0.25d0 / BigR * vpar0**2							       &
    				      * r0 * (Te0_x * psi_y - Te0_y * psi_x)							       &
    				      * ( v_x * ps0_y -  v_y * ps0_x )  					* xjac * theta * tstep

    		  amat_81_k = - (ZK_e_par-ZK_e_prof) * BigR * BB2_psi / BB2**2 * Bgrad_Te_k_star * Bgrad_Te	       * xjac * theta * tstep &
    			      + (ZK_e_par-ZK_e_prof) * BigR / BB2	   * Bgrad_Te_k_star * Bgrad_Te_psi	       * xjac * theta * tstep &
    			      !- dZK_e_prof * psi * BigR / BB2 * Bgrad_Te_k_star * Bgrad_Te			* xjac * theta * tstep &
    			      + TG_num8 * tstep * 0.25d0 / BigR * vpar0**2							       &
    					* Te0 * (r0_x * psi_y - r0_y * psi_x)							       &
    					* ( F0 / BigR * v_p)							* xjac * theta * tstep &
    			      + TG_num8 * tstep * 0.25d0 / BigR * vpar0**2							       &
    					* r0 * (Te0_x * psi_y - Te0_y * psi_x)  						       &
    					* ( F0 / BigR * v_p)							* xjac * theta * tstep    

    		  ! Current terms (not yet tested!)
    		  amat_81_n = 0.d0 !- sigma * v / BigR * Te0_x * psi_xp 					* xjac * theta * tstep & 
    				   !- sigma * v / BigR * Te0_y * psi_yp 					* xjac * theta * tstep &
    				   !+ v * sigma * (GAMMA-1.d0) * Te0 / (BigR*r0) * r0_x * psi_xp		* xjac * theta * tstep &
    				   !+ v * sigma * (GAMMA-1.d0) * Te0 / (BigR*r0) * r0_y * psi_yp		* xjac * theta * tstep 

    		  amat_82 = - v * r0 * BigR**2 * ( Te0_x * u_y - Te0_y * u_x)					* xjac * theta * tstep &
    			    - v * 2.d0 * (GAMMA-1.d0) * r0 * BigR * Te0 * u_y					* xjac * theta * tstep &
    			    + TG_num8 * tstep * 0.25d0 * BigR**2 * Te0* (r0_x * u_y - r0_y * u_x)				       &
    				      * ( v_x * u0_y - v_y * u0_x)						* xjac * theta * tstep &
    			    + TG_num8 * tstep * 0.25d0 * BigR**2 * r0* (Te0_x * u_y - Te0_y * u_x)				       &
    				      * ( v_x * u0_y - v_y * u0_x)						* xjac * theta * tstep &
    			    + TG_num8 * tstep * 0.25d0 * BigR**2 * Te0* (r0_x * u0_y - r0_y * u0_x)				       &
    				      * ( v_x * u_y - v_y * u_x)						* xjac * theta * tstep &
    			    + TG_num8 * tstep * 0.25d0 * BigR**2 * r0* (Te0_x * u0_y - Te0_y * u0_x)				       &
    				      * ( v_x * u_y - v_y * u_x)						* xjac * theta * tstep
    							 
    		  ! Current terms (not yet tested!)
    		  amat_83 = 0.d0 !+ v * sigma * zj / BigR * Te0_p						* xjac * theta * tstep &
    				 !- v * sigma * (GAMMA-1.d0) * zj * Te0 / (BigR*r0) * r0_p			* xjac * theta * tstep 
    			    
    			    
    		  amat_85 = - v * rho * BigR**2 * (Te0_x * u0_y - Te0_y * u0_x) 				* xjac * theta * tstep &
    			    + v * rho * Vpar0 * F0/BigR * Te0_p 						* xjac * theta * tstep &
    			    + v * rho * Vpar0 * (Te0_x * ps0_y - Te0_y * ps0_x) 				* xjac * theta * tstep &
    			    - v * 2.d0 * (GAMMA-1.d0) * rho * BigR * Te0 * u0_y 				* xjac * theta * tstep &
    			    + v * (GAMMA-1.d0) * rho * Te0 * F0/BigR * Vpar0_p  				* xjac * theta * tstep &
    			    + v * (GAMMA-1.d0) * rho * Te0 * (Vpar0_x * ps0_y - Vpar0_y * ps0_x)		* xjac * theta * tstep &  
    			    ! Current terms (not yet tested!)
    			    !+ v * sigma * (GAMMA-1.d0) * Te0 / (BigR*r0) * rho_x * ps0_xp			* xjac * theta * tstep &
    			    !+ v * sigma * (GAMMA-1.d0) * Te0 / (BigR*r0) * rho_y * ps0_yp			* xjac * theta * tstep &
    			    !- v * sigma * (GAMMA-1.d0) * Te0 * rho / (BigR*r0**2) * r0_x * ps0_xp		* xjac * theta * tstep &
    			    !- v * sigma * (GAMMA-1.d0) * Te0 * rho / (BigR*r0**2) * r0_y * ps0_yp		* xjac * theta * tstep &
    			    !+ v * sigma * (GAMMA-1.d0) * Te0 * rho / (BigR*r0**2) * zj0 * r0_p 		* xjac * theta * tstep &
    			    ! Temperature transfer (not yet tested!)
    			    !+ v * BigR * dQ_drho * rho * (Te0 - Ti0)						* xjac * theta * tstep & 
    			    + TG_num8 * tstep * 0.25d0 * BigR**2 * Te0* (rho_x * u0_y - rho_y * u0_x)				       &
    				      * ( v_x * u0_y - v_y * u0_x)						* xjac * theta * tstep &
    			    + TG_num8 * tstep * 0.25d0 * BigR**2 * rho * (Te0_x * u0_y - Te0_y * u0_x)  			       &
    				      * ( v_x * u0_y - v_y * u0_x)						* xjac * theta * tstep &
    			    + TG_num8 * tstep * 0.25d0 / BigR * vpar0**2							       &
    				      * Te0 * (rho_x * ps0_y - rho_y * ps0_x + F0 / BigR * rho_p)				       &
    				      * ( v_x * ps0_y -  v_y * ps0_x )  					* xjac * theta * tstep &
    			    + TG_num8 * tstep * 0.25d0 / BigR * vpar0**2							       &
    				      * rho * (Te0_x * ps0_y - Te0_y * ps0_x + F0 / BigR * Te0_p)				       &
    				      * ( v_x * ps0_y -  v_y * ps0_x )  					* xjac * theta * tstep  

    		  amat_85_k = + TG_num8 * tstep * 0.25d0 / BigR * vpar0**2							       &
    					* Te0 * (rho_x * ps0_y - rho_y * ps0_x + F0 / BigR * rho_p)				       &
    					* ( F0 / BigR * v_p)							* xjac * theta * tstep &
    			      + TG_num8 * tstep * 0.25d0 / BigR * vpar0**2							       &
    					* rho * (Te0_x * ps0_y - Te0_y * ps0_x + F0 / BigR * Te0_p)				       &
    					* ( F0 / BigR * v_p)							* xjac * theta * tstep  
    			  
    		  ! Current terms (not yet tested!)
    		  amat_85_n = 0.d0!- v * sigma * (GAMMA-1.d0) * Te0 / (BigR*r0) * zj0 * rho_p			* xjac * theta * tstep 

    		  ! Temperature transfer (not yet tested!)
    		  amat_86 = 0.d0!- v * BigR * Q * Ti								* xjac * theta * tstep      

    		  amat_87 = + v * r0 * F0/BigR * Vpar * Te0_p							* xjac * theta * tstep &
    			    + v * r0 * Vpar * (Te0_x * ps0_y - Te0_y * ps0_x)					* xjac * theta * tstep &
    			    + v * (GAMMA-1.d0) * r0 * Te0 * (Vpar_x * ps0_y - Vpar_y * ps0_x)			* xjac * theta * tstep &
    			    + TG_num8 * tstep * 0.25d0 / BigR * 2.d0 * vpar0*vpar						       &
    				      * Te0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)  				       &
    				      * ( v_x * ps0_y -  v_y * ps0_x )  					* xjac * theta * tstep &
    			    + TG_num8 * tstep * 0.25d0 / BigR * 2.d0 * vpar0*vpar						       &
    				      * r0 * (Te0_x * ps0_y - Te0_y * ps0_x + F0 / BigR * Te0_p)				       &
    				      * ( v_x * ps0_y -  v_y * ps0_x )  					* xjac * theta * tstep
    		 
    		  amat_87_k = + TG_num8 * tstep * 0.25d0 / BigR * 2.d0 * vpar0*vpar						       &
    					* Te0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)				       &
    					* ( F0 / BigR * v_p)							* xjac * theta * tstep &
    			      + TG_num8 * tstep * 0.25d0 / BigR * 2.d0 * vpar0*vpar						       &
    					* r0 * (Te0_x * ps0_y - Te0_y * ps0_x + F0 / BigR * Te0_p)				       &
    					* ( F0 / BigR * v_p)							* xjac * theta * tstep
    						    
    		  amat_87_n = + v * (GAMMA-1.d0) * r0 * Te0 * F0/BigR * Vpar_p  				* xjac * theta * tstep 
    						    
    		  amat_88 = v * abs(r0) * Te * BigR								* xjac * (1.d0 + zeta) &
    			    - v * r0 * BigR**2 * ( Te_x * u0_y - Te_y * u0_x)					* xjac * theta * tstep &
    			    + v * r0 * Vpar0 * (Te_x * ps0_y - Te_y * ps0_x)					* xjac * theta * tstep &
    			    - 2.d0 * v * r0 * (GAMMA-1.d0) * Te * BigR * u0_y					* xjac * theta * tstep &
    			    + v * r0 * (GAMMA-1.d0) * Te * F0/BigR * Vpar0_p					* xjac * theta * tstep &
    			    + v * r0 * (GAMMA-1.d0) * Te * (Vpar0_x * ps0_y - Vpar0_y * ps0_x)  		* xjac * theta * tstep &
    			    ! Current terms (not yet tested!)
    			    !- v * sigma / BigR * Te_x * ps0_xp 						* xjac * theta * tstep &
    			    !- v * sigma / BigR * Te_y * ps0_yp 						* xjac * theta * tstep &
    			    !+ v * sigma * (GAMMA-1.d0) * Te / (BigR*r0) * r0_x * ps0_xp			* xjac * theta * tstep &
    			    !+ v * sigma * (GAMMA-1.d0) * Te / (BigR*r0) * r0_y * ps0_yp			* xjac * theta * tstep &
    			    !- v * sigma * (GAMMA-1.d0) * Te / (BigR*r0) * zj0 * r0_p				* xjac * theta * tstep &
    			    ! Temperature transfer (not yet tested!)
    			    !+ v * BigR * dQ_dTe * Te * (Te0 - Ti0)						* xjac * theta * tstep & 
    			    !+ v * BigR * Q * Te								* xjac * theta * tstep &     
    			    + (ZK_e_par-ZK_e_prof) * BigR / BB2 * Bgrad_Te_star * Bgrad_Te_Te				* xjac * theta * tstep &
                            + dZK_e_par * Te * BigR / BB2 * Bgrad_Te_star * Bgrad_Te                            * xjac * theta * tstep &
    			    + ZK_e_prof * BigR * (v_x*Te_x + v_y*Te_y ) 					* xjac * theta * tstep & 
    			    + K_perp_numm  *											       &
    			      (    (	     v_ss  * (x_t**2 + y_t**2)  							       &
    				    +	     v_tt  * (x_s**2 + y_s**2)  							       &
    				    - 2.d0 * v_st  * (x_s*x_t + y_s*y_t) )							       &
    				 * (	     Te_ss * (x_t**2 + y_t**2)  							       &
    				    +	     Te_tt * (x_s**2 + y_s**2)  							       &
    				    - 2.d0 * Te_st * (x_s*x_t + y_s*y_t) ) )	 / xjac**4			* xjac * theta * tstep & 
    			    + TG_num8 * tstep * 0.25d0 * BigR**2 * Te * (r0_x * u0_y - r0_y * u0_x)				       &
    				      * ( v_x * u0_y - v_y * u0_x)						* xjac * theta * tstep &
    			    + TG_num8 * tstep * 0.25d0 * BigR**2 * r0* (Te_x * u0_y - Te_y * u0_x)				       &
    				      * ( v_x * u0_y - v_y * u0_x)						* xjac * theta * tstep &
    			    + TG_num8 * tstep * 0.25d0 / BigR * vpar0**2							       &
    				      * Te * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)					       &
    				      * ( v_x * ps0_y -  v_y * ps0_x )  					* xjac * theta * tstep &
    			    + TG_num8 * tstep * 0.25d0 / BigR * vpar0**2							       &
    				      * r0 * (Te_x * ps0_y - Te_y * ps0_x )							       &
    				      * ( v_x * ps0_y -  v_y * ps0_x )  					* xjac * theta * tstep 

    		  amat_88_k = + (ZK_e_par-ZK_e_prof) * BigR / BB2 * Bgrad_Te_k_star * Bgrad_Te_Te  			* xjac * theta * tstep &
                              + dZK_e_par * Te * BigR / BB2 * Bgrad_Te_k_star * Bgrad_Te                        * xjac * theta * tstep &
    			      + TG_num8 * tstep * 0.25d0 / BigR * vpar0**2							       &
    					* Te * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p) 				       &
    					* ( F0 / BigR * v_p)							* xjac * theta * tstep &
    			      + TG_num8 * tstep * 0.25d0 / BigR * vpar0**2							       &
    					* r0 * (Te_x * ps0_y - Te_y * ps0_x )							       &
    					* ( F0 / BigR * v_p)							* xjac * theta * tstep
    		     
    		  amat_88_n = + v * r0 * Vpar0 * F0/BigR * Te_p 						* xjac * theta * tstep &
    			      ! Current terms (not yet tested!)
    			      !+ v ** sigma / BigR * zj0 * Te_p 						* xjac * theta * tstep &
    			      + (ZK_e_par-ZK_e_prof) * BigR / BB2 * Bgrad_Te_star * Bgrad_Te_Te_n  			* xjac * theta * tstep &
    			      + TG_num8 * tstep * 0.25d0 / BigR * vpar0**2							       &
    					* r0 * ( F0 / BigR * Te_p)								       &
    					* ( v_x * ps0_y -  v_y * ps0_x )					* xjac * theta * tstep
    				     
    		  amat_88_kn = + (ZK_e_par-ZK_e_prof) * BigR / BB2 * Bgrad_Te_k_star * Bgrad_Te_Te_n			* xjac * theta * tstep &
    			       + ZK_e_prof * v_p*Te_p / BigR							* xjac * theta * tstep &
    			       + TG_num8 * tstep * 0.25d0 / BigR * vpar0**2							       &
    					 * r0 * ( F0 / BigR * Te_p)								       &
    					 * ( F0 / BigR * v_p)							* xjac * theta * tstep  

!    		  amat_89 = - v * KPe * BigR									* xjac * theta * tstep  

    		  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    		  !!!!!!!!!! LHS equation 9 (Grad(Kpar.Grad(Te)) !!!!!!!!!
    		  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    		  amat_91    = - ZK_e_par * BB2_psi / BB2**2 * Bgrad_Te_star	 * Bgrad_Te        * BigR * xjac * tstep &
!    			       + ZK_e_par           / BB2    * Bgrad_Te_star_psi * Bgrad_Te        * BigR * xjac * tstep &
!    			       + ZK_e_par           / BB2    * Bgrad_Te_star     * Bgrad_Te_psi    * BigR * xjac * tstep 

!    		  amat_91_k  = - ZK_e_par * BB2_psi / BB2**2 * Bgrad_Te_k_star   * Bgrad_Te        * BigR * xjac * tstep &
!    			       + ZK_e_par           / BB2    * Bgrad_Te_k_star   * Bgrad_Te_psi    * BigR * xjac * tstep 

!    		  amat_98    = + ZK_e_par           / BB2    * Bgrad_Te_star     * Bgrad_Te_Te     * BigR * xjac * tstep &
!    			       + dZK_e_par * Te     / BB2    * Bgrad_Te_star	 * Bgrad_Te	   * BigR * xjac * tstep 

!    		  amat_98_k  = + ZK_e_par           / BB2    * Bgrad_Te_k_star   * Bgrad_Te_Te     * BigR * xjac * tstep &
!    			       + dZK_e_par * Te     / BB2    * Bgrad_Te_k_star   * Bgrad_Te        * BigR * xjac * tstep 
    		     
!    		  amat_98_n  = + ZK_e_par	    / BB2    * Bgrad_Te_star	 * Bgrad_Te_Te_n   * BigR * xjac * tstep 
    				     
!    		  amat_98_kn = + ZK_e_par	    / BB2    * Bgrad_Te_k_star   * Bgrad_Te_Te_n   * BigR * xjac * tstep 

!    		  amat_99    = v * KPe  							   * BigR * xjac * tstep 

            	  ! --------------------------
	    	  ! --- Now fill up the matrix
            	  ! --------------------------

	    	  ! --- Indices for the matrix
 	    	  kl1 = index_kl
 	    	  kl2 = index_kl + 1*n_tor_loop2
 	    	  kl3 = index_kl + 2*n_tor_loop2
 	    	  kl4 = index_kl + 3*n_tor_loop2
 	    	  kl5 = index_kl + 4*n_tor_loop2
 	    	  kl6 = index_kl + 5*n_tor_loop2
 	    	  kl7 = index_kl + 6*n_tor_loop2
 	    	  kl8 = index_kl + 7*n_tor_loop2
! 	    	  kl9 = index_kl + 8*n_tor_loop2

            	  ! --- Equation 1 (psi - induction)
                  if (n_tor .gt. 3) then
 	    	    ELM_p(mp,ij1,kl1)  =  ELM_p(mp,ij1,kl1)  + wst * amat_11
 	    	    ELM_p(mp,ij1,kl2)  =  ELM_p(mp,ij1,kl2)  + wst * amat_12  
 	    	    ELM_n(mp,ij1,kl2)  =  ELM_n(mp,ij1,kl2)  + wst * amat_12_n
 	    	    ELM_p(mp,ij1,kl3)  =  ELM_p(mp,ij1,kl3)  + wst * amat_13
 	    	    ELM_p(mp,ij1,kl8)  =  ELM_p(mp,ij1,kl8)  + wst * amat_18
		  else
		    ELM(ij1,kl1) = ELM(ij1,kl1) + wst * amat_11
		    ELM(ij1,kl2) = ELM(ij1,kl2) + wst * amat_12 + wst * amat_12_n
		    ELM(ij1,kl3) = ELM(ij1,kl3) + wst * amat_13
		    ELM(ij1,kl8) = ELM(ij1,kl8) + wst * amat_18
		  endif

            	  ! --- Equation 2 (U - momentum)
                  if (n_tor .gt. 3) then
 	    	    ELM_p(mp,ij2,kl1)  =  ELM_p(mp,ij2,kl1)  + wst * amat_21
 	    	    ELM_p(mp,ij2,kl2)  =  ELM_p(mp,ij2,kl2)  + wst * amat_22
 	    	    ELM_p(mp,ij2,kl3)  =  ELM_p(mp,ij2,kl3)  + wst * amat_23
 	    	    ELM_n(mp,ij2,kl3)  =  ELM_n(mp,ij2,kl3)  + wst * amat_23_n
 	    	    ELM_p(mp,ij2,kl4)  =  ELM_p(mp,ij2,kl4)  + wst * amat_24
 	    	    ELM_p(mp,ij2,kl5)  =  ELM_p(mp,ij2,kl5)  + wst * amat_25
 	    	    ELM_p(mp,ij2,kl6)  =  ELM_p(mp,ij2,kl6)  + wst * amat_26
 	    	    ELM_p(mp,ij2,kl8)  =  ELM_p(mp,ij2,kl8)  + wst * amat_28
		  else
		    ELM(ij2,kl1) = ELM(ij2,kl1) + wst * amat_21
		    ELM(ij2,kl2) = ELM(ij2,kl2) + wst * amat_22
		    ELM(ij2,kl3) = ELM(ij2,kl3) + wst * amat_23 + wst * amat_23_n
		    ELM(ij2,kl4) = ELM(ij2,kl4) + wst * amat_24
		    ELM(ij2,kl5) = ELM(ij2,kl5) + wst * amat_25
		    ELM(ij2,kl6) = ELM(ij2,kl6) + wst * amat_26
		    ELM(ij2,kl8) = ELM(ij2,kl8) + wst * amat_28
		  endif

            	  ! --- Equation 3 (j - current)
                  if (n_tor .gt. 3) then
 	    	    ELM_p(mp,ij3,kl1)  =  ELM_p(mp,ij3,kl1)  + wst * amat_31
 	    	    ELM_p(mp,ij3,kl3)  =  ELM_p(mp,ij3,kl3)  + wst * amat_33
		  else
		    ELM(ij3,kl1) = ELM(ij3,kl1) + wst * amat_31
		    ELM(ij3,kl3) = ELM(ij3,kl3) + wst * amat_33
		  endif

	    	  ! --- Equation 4 (w - vorticity)
                  if (n_tor .gt. 3) then
 	    	    ELM_p(mp,ij4,kl2)  =  ELM_p(mp,ij4,kl2)  + wst * amat_42
 	    	    ELM_p(mp,ij4,kl4)  =  ELM_p(mp,ij4,kl4)  + wst * amat_44
		  else
		    ELM(ij4,kl2) = ELM(ij4,kl2) + wst * amat_42
		    ELM(ij4,kl4) = ELM(ij4,kl4) + wst * amat_44
		  endif

	    	  ! --- Equation 5 (rho - continuity)
                  if (n_tor .gt. 3) then
 	    	    ELM_p(mp,ij5,kl1)  =  ELM_p(mp,ij5,kl1)  + wst * amat_51
 	    	    ELM_k(mp,ij5,kl1)  =  ELM_k(mp,ij5,kl1)  + wst * amat_51_k
 	    	    ELM_p(mp,ij5,kl2)  =  ELM_p(mp,ij5,kl2)  + wst * amat_52
 	    	    ELM_p(mp,ij5,kl5)  =  ELM_p(mp,ij5,kl5)  + wst * amat_55
 	    	    ELM_k(mp,ij5,kl5)  =  ELM_k(mp,ij5,kl5)  + wst * amat_55_k
 	    	    ELM_n(mp,ij5,kl5)  =  ELM_n(mp,ij5,kl5)  + wst * amat_55_n
 	    	    ELM_kn(mp,ij5,kl5) =  ELM_kn(mp,ij5,kl5) + wst * amat_55_kn
 	    	    ELM_p(mp,ij5,kl7)  =  ELM_p(mp,ij5,kl7)  + wst * amat_57
 	    	    ELM_k(mp,ij5,kl7)  =  ELM_k(mp,ij5,kl7)  + wst * amat_57_k
 	    	    ELM_n(mp,ij5,kl7)  =  ELM_n(mp,ij5,kl7)  + wst * amat_57_n
		  else
		    ELM(ij5,kl1) = ELM(ij5,kl1) + wst * amat_51 + wst * amat_51_k
		    ELM(ij5,kl2) = ELM(ij5,kl2) + wst * amat_52
		    ELM(ij5,kl5) = ELM(ij5,kl5) + wst * amat_55 + wst * amat_55_k + wst * amat_55_n + wst * amat_55_kn
		    ELM(ij5,kl7) = ELM(ij5,kl7) + wst * amat_57 + wst * amat_57_k + wst * amat_57_n
		  endif

	    	  ! --- Equation 6 (Ti - Ion energy)
                  if (n_tor .gt. 3) then
 	    	    ELM_p(mp,ij6,kl1)  =  ELM_p(mp,ij6,kl1)  + wst * amat_61
 	    	    ELM_k(mp,ij6,kl1)  =  ELM_k(mp,ij6,kl1)  + wst * amat_61_k
 	    	    ELM_p(mp,ij6,kl2)  =  ELM_p(mp,ij6,kl2)  + wst * amat_62
 	    	    ELM_p(mp,ij6,kl5)  =  ELM_p(mp,ij6,kl5)  + wst * amat_65   
 	    	    ELM_k(mp,ij6,kl5)  =  ELM_k(mp,ij6,kl5)  + wst * amat_65_k
 	    	    ELM_p(mp,ij6,kl6)  =  ELM_p(mp,ij6,kl6)  + wst * amat_66
 	    	    ELM_k(mp,ij6,kl6)  =  ELM_k(mp,ij6,kl6)  + wst * amat_66_k
 	    	    ELM_n(mp,ij6,kl6)  =  ELM_n(mp,ij6,kl6)  + wst * amat_66_n
 	    	    ELM_kn(mp,ij6,kl6) =  ELM_kn(mp,ij6,kl6) + wst * amat_66_kn
 	    	    ELM_p(mp,ij6,kl7)  =  ELM_p(mp,ij6,kl7)  + wst * amat_67
 	    	    ELM_k(mp,ij6,kl7)  =  ELM_k(mp,ij6,kl7)  + wst * amat_67_k
 	    	    ELM_n(mp,ij6,kl7)  =  ELM_n(mp,ij6,kl7)  + wst * amat_67_n
 	    	    ELM_p(mp,ij6,kl8)  =  ELM_p(mp,ij6,kl8)  + wst * amat_68
		  else
		    ELM(ij6,kl1) = ELM(ij6,kl1) + wst * amat_61 + wst * amat_61_k
		    ELM(ij6,kl2) = ELM(ij6,kl2) + wst * amat_62
		    ELM(ij6,kl5) = ELM(ij6,kl5) + wst * amat_65 + wst * amat_65_k
		    ELM(ij6,kl6) = ELM(ij6,kl6) + wst * amat_66 + wst * amat_66_k + wst * amat_66_n + wst * amat_66_kn
		    ELM(ij6,kl7) = ELM(ij6,kl7) + wst * amat_67 + wst * amat_67_k + wst * amat_67_n
		    ELM(ij6,kl8) = ELM(ij6,kl8) + wst * amat_68
		  endif

	    	  ! --- Equation 7 (Vpar - parallel momentum)
                  if (n_tor .gt. 3) then
 	    	    ELM_p(mp,ij7,kl1)  =  ELM_p(mp,ij7,kl1)  + wst * amat_71
 	    	    ELM_p(mp,ij7,kl2)  =  ELM_p(mp,ij7,kl2)  + wst * amat_72
 	    	    ELM_p(mp,ij7,kl5)  =  ELM_p(mp,ij7,kl5)  + wst * amat_75
 	    	    ELM_k(mp,ij7,kl5)  =  ELM_k(mp,ij7,kl5)  + wst * amat_75_k
 	    	    ELM_n(mp,ij7,kl5)  =  ELM_n(mp,ij7,kl5)  + wst * amat_75_n
 	    	    ELM_p(mp,ij7,kl6)  =  ELM_p(mp,ij7,kl6)  + wst * amat_76
 	    	    ELM_n(mp,ij7,kl6)  =  ELM_n(mp,ij7,kl6)  + wst * amat_76_n
 	    	    ELM_p(mp,ij7,kl7)  =  ELM_p(mp,ij7,kl7)  + wst * amat_77
 	    	    ELM_k(mp,ij7,kl7)  =  ELM_k(mp,ij7,kl7)  + wst * amat_77_k
 	    	    ELM_n(mp,ij7,kl7)  =  ELM_n(mp,ij7,kl7)  + wst * amat_77_n
 	    	    ELM_kn(mp,ij7,kl7) =  ELM_kn(mp,ij7,kl7) + wst * amat_77_kn
 	    	    ELM_p(mp,ij7,kl8)  =  ELM_p(mp,ij7,kl8)  + wst * amat_78
 	    	    ELM_n(mp,ij7,kl8)  =  ELM_n(mp,ij7,kl8)  + wst * amat_78_n
		  else 	    	  
		    ELM(ij7,kl1) = ELM(ij7,kl1) + wst * amat_71
		    ELM(ij7,kl2) = ELM(ij7,kl2) + wst * amat_72
		    ELM(ij7,kl5) = ELM(ij7,kl5) + wst * amat_75 + wst * amat_75_k + wst * amat_75_n
		    ELM(ij7,kl6) = ELM(ij7,kl6) + wst * amat_76 + wst * amat_76_n
		    ELM(ij7,kl7) = ELM(ij7,kl7) + wst * amat_77 + wst * amat_77_k + wst * amat_77_n + wst * amat_77_kn
		    ELM(ij7,kl8) = ELM(ij7,kl8) + wst * amat_78 + wst * amat_78_n
		  endif

	    	  ! --- Equation 8 (Te - electron energy)
                  if (n_tor .gt. 3) then
 	    	    ELM_p(mp,ij8,kl1)  =  ELM_p(mp,ij8,kl1)  + wst * amat_81
 	    	    ELM_k(mp,ij8,kl1)  =  ELM_k(mp,ij8,kl1)  + wst * amat_81_k
 	    	    ELM_n(mp,ij8,kl1)  =  ELM_n(mp,ij8,kl1)  + wst * amat_81_n
 	    	    ELM_p(mp,ij8,kl2)  =  ELM_p(mp,ij8,kl2)  + wst * amat_82
 	    	    ELM_p(mp,ij8,kl3)  =  ELM_p(mp,ij8,kl3)  + wst * amat_83
 	    	    ELM_p(mp,ij8,kl5)  =  ELM_p(mp,ij8,kl5)  + wst * amat_85
 	    	    ELM_k(mp,ij8,kl5)  =  ELM_k(mp,ij8,kl5)  + wst * amat_85_k
 	    	    ELM_n(mp,ij8,kl5)  =  ELM_n(mp,ij8,kl5)  + wst * amat_85_n
 	    	    ELM_p(mp,ij8,kl6)  =  ELM_p(mp,ij8,kl6)  + wst * amat_86
 	    	    ELM_p(mp,ij8,kl7)  =  ELM_p(mp,ij8,kl7)  + wst * amat_87
 	    	    ELM_k(mp,ij8,kl7)  =  ELM_k(mp,ij8,kl7)  + wst * amat_87_k
 	    	    ELM_n(mp,ij8,kl7)  =  ELM_n(mp,ij8,kl7)  + wst * amat_87_n
 	    	    ELM_p(mp,ij8,kl8)  =  ELM_p(mp,ij8,kl8)  + wst * amat_88
 	    	    ELM_k(mp,ij8,kl8)  =  ELM_k(mp,ij8,kl8)  + wst * amat_88_k
 	    	    ELM_n(mp,ij8,kl8)  =  ELM_n(mp,ij8,kl8)  + wst * amat_88_n
 	    	    ELM_kn(mp,ij8,kl8) =  ELM_kn(mp,ij8,kl8) + wst * amat_88_kn
! 	    	    ELM_p(mp,ij8,kl9)  =  ELM_p(mp,ij8,kl9)  + wst * amat_89
		  else 	    	  		  
		    ELM(ij8,kl1) = ELM(ij8,kl1) + wst * amat_81 + wst * amat_81_k + wst * amat_81_n
		    ELM(ij8,kl2) = ELM(ij8,kl2) + wst * amat_82
		    ELM(ij8,kl3) = ELM(ij8,kl3) + wst * amat_83
		    ELM(ij8,kl5) = ELM(ij8,kl5) + wst * amat_85 + wst * amat_85_k + wst * amat_85_n
		    ELM(ij8,kl6) = ELM(ij8,kl6) + wst * amat_86
		    ELM(ij8,kl7) = ELM(ij8,kl7) + wst * amat_87 + wst * amat_87_k + wst * amat_87_n
		    ELM(ij8,kl8) = ELM(ij8,kl8) + wst * amat_88 + wst * amat_88_k + wst * amat_88_n + wst * amat_88_kn
!		    ELM(ij8,kl9) = ELM(ij8,kl9) + wst * amat_89
		  endif
 	    	
	    	  ! --- Equation 9 (Grad(Kpar.Grad(Te))
!                  if (n_tor .gt. 3) then
! 	    	    ELM_p(mp,ij9,kl1)  =  ELM_p(mp,ij9,kl1)  + wst * amat_91
! 	    	    ELM_k(mp,ij9,kl1)  =  ELM_k(mp,ij9,kl1)  + wst * amat_91_k
! 	    	    ELM_p(mp,ij9,kl8)  =  ELM_p(mp,ij9,kl8)  + wst * amat_98
! 	    	    ELM_k(mp,ij9,kl8)  =  ELM_k(mp,ij9,kl8)  + wst * amat_98_k
! 	    	    ELM_n(mp,ij9,kl8)  =  ELM_n(mp,ij9,kl8)  + wst * amat_98_n
! 	    	    ELM_kn(mp,ij9,kl8) =  ELM_kn(mp,ij9,kl8) + wst * amat_98_kn
! 	    	    ELM_p(mp,ij9,kl9)  =  ELM_p(mp,ij9,kl9)  + wst * amat_99
!		  else 	    	  		  
!		    ELM(ij9,kl1) = ELM(ij9,kl1) + wst * amat_91 + wst * amat_91_k
!		    ELM(ij9,kl8) = ELM(ij9,kl8) + wst * amat_98 + wst * amat_98_k + wst * amat_98_n + wst * amat_98_kn
!		    ELM(ij9,kl9) = ELM(ij9,kl9) + wst * amat_99
!		  endif
 	    	
 	        enddo ! inner n_tor_loop
 	      enddo   ! inner n_order+1
 	    enddo     ! inner n_vertex_max

 	  enddo       ! outer n_tor_loop
 	enddo	      ! outer n_order+1
      enddo	      ! outer n_vertex_max

    enddo	      ! n_plane

  enddo	              ! n_gauss
enddo	              ! n_gauss


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! Apply FFT !!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (n_tor .gt. 3) then

  do i=1,n_vertex_max*n_var*(n_order+1)

    do j=1, n_vertex_max*n_var*(n_order+1)

      in_fft =  ELM_p(1:n_plane,i,j)

      call my_fft(in_fft, out_fft, n_plane)

      do k=1,(n_tor+1)/2

    	index_k = n_tor*(i-1) + max(2*(k-1),1)

    	do m=1,(n_tor+1)/2

    	  index_m = n_tor*(j-1) + max(2*(m-1),1)

    	  l = (k-1) + (m-1)

    	  if ( (l .ge. 0) .and. (l .lt. n_plane/2) ) then

    	    ELM(index_k,  index_m  ) = ELM(index_k,  index_m)	+ real(out_fft(l+1))
    	    ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)	- imag(out_fft(l+1))
    	    ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) - imag(out_fft(l+1))
    	    ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) - real(out_fft(l+1))

    	  elseif ( (l .lt. 0) .and. (abs(l) .lt. n_plane/2) ) then

    	    ELM(index_k,  index_m  ) = ELM(index_k,  index_m)	+ real(out_fft(abs(l)+1))
    	    ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)	+ imag(out_fft(abs(l)+1))
    	    ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + imag(out_fft(abs(l)+1))
    	    ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) - real(out_fft(abs(l)+1))

    	  endif

    	  l = (k-1) - (m-1)

    	  if ( (l .ge. 0) .and. (l .lt. n_plane/2) ) then

    	    ELM(index_k,  index_m  ) = ELM(index_k,  index_m)	+ real(out_fft(l+1))
    	    ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)	- imag(out_fft(l+1))
    	    ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + imag(out_fft(l+1))
    	    ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + real(out_fft(l+1))

    	  elseif ( (l .lt. 0) .and. (abs(l) .lt. n_plane/2) ) then

    	    ELM(index_k,  index_m  ) = ELM(index_k,  index_m)	+ real(out_fft(abs(l)+1))
    	    ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)	+ imag(out_fft(abs(l)+1))
    	    ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) - imag(out_fft(abs(l)+1))
    	    ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + real(out_fft(abs(l)+1))

    	  endif

    	enddo
      enddo
    enddo
  enddo

  do i=1,n_vertex_max*n_var*(n_order+1)

    do j=1, n_vertex_max*n_var*(n_order+1)

      if (maxval(abs(ELM_n(1:n_plane,i,j))) .ne. 0.d0) then

    	in_fft =  ELM_n(1:n_plane,i,j)

    	call my_fft(in_fft, out_fft, n_plane)

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
      endif
    enddo
  enddo

  do i=1,n_vertex_max*n_var*(n_order+1)

    do j=1, n_vertex_max*n_var*(n_order+1)

      if (maxval(abs(ELM_k(1:n_plane,i,j))) .ne. 0.d0) then

    	in_fft =  ELM_k(1:n_plane,i,j)

    	call my_fft(in_fft, out_fft, n_plane)

    	do k=1,(n_tor+1)/2

    	  ik	  = max(2*(k-1),1)
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
      endif
    enddo
  enddo


  do i=1,n_vertex_max*n_var*(n_order+1)

    do j=1, n_vertex_max*n_var*(n_order+1)

      if (maxval(abs(ELM_kn(1:n_plane,i,j))) .ne. 0.d0) then

    	in_fft =  ELM_kn(1:n_plane,i,j)

    	call my_fft(in_fft, out_fft, n_plane)

    	do k=1,(n_tor+1)/2

    	  ik	  = max(2*(k-1),1)
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
      endif
    enddo
  enddo

  ELM = 0.5d0 * ELM

  do j=1, n_vertex_max*n_var*(n_order+1)

    in_fft = RHS_p(1:n_plane,j)

    call my_fft(in_fft, out_fft, n_plane)

    index = n_tor*(j-1) + 1

    RHS(index) = real(out_fft(1))

    do k=2,(n_tor+1)/2

      index = n_tor*(j-1) + 2*(k-1)

      RHS(index)   =   real(out_fft(k))
      RHS(index+1) = - imag(out_fft(k))

    enddo

  enddo

  do j=1, n_vertex_max*n_var*(n_order+1)

    in_fft = RHS_k(1:n_plane,j)

    call my_fft(in_fft, out_fft, n_plane)

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

endif

return
end subroutine element_matrix





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! The FFT routine !!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine my_fft(in_fft,out_fft,n)
      
  real*8     :: in_fft(*)
  complex*16 :: out_fft(*)      
  real*8     :: tmp_fft(2*n+2)
  integer    :: i
      
  tmp_fft(1:n) = in_fft(1:n)      
  call RFT2(tmp_fft,n,1)
      
  do i=1,n
    out_fft(i) = cmplx(tmp_fft(2*i-1),tmp_fft(2*i))
  enddo
      
  return
  
end subroutine my_fft




end module mod_elt_matrix
