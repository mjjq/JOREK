! Contains all the variables for the equations, to share with each equation routine
module equation_variables

  use parameters
  
  implicit none

    
  ! --- Profiles and sources
  real*8 	:: zn,  dn_dpsi,  dn_dz,  dn_dpsi2,  dn_dz2,  dn_dpsi_dz,  dn_dpsi3,  dn_dpsi_dz2,  dn_dpsi2_dz
  real*8 	:: zTi, dTi_dpsi, dTi_dz, dTi_dpsi2, dTi_dz2, dTi_dpsi_dz, dTi_dpsi3, dTi_dpsi_dz2, dTi_dpsi2_dz
  real*8 	:: zTe, dTe_dpsi, dTe_dz, dTe_dpsi2, dTe_dz2, dTe_dpsi_dz, dTe_dpsi3, dTe_dpsi_dz2, dTe_dpsi2_dz
  real*8 	:: current_source, particle_source, heat_source_i, heat_source_e, V_source
  
  ! --- Diffusivities
  real*8 	:: eta_Te,    deta_dTe, d2eta_d2Te
  real*8 	:: visco_Te,  dvisco_dTe
  real*8 	:: D_prof
  real*8 	:: Ki_prof, Ki_par, dKi_par
  real*8 	:: Ke_prof, Ke_par, dKe_par
  
  ! --- Hyper diffusivities
  real*8 	:: eta_numm, visco_numm, visco_par_numm, D_perp_numm, Ki_perp_numm, Ki_par_num, Ke_perp_numm, Ke_par_num 
  real*8 	:: v_ps0_x,   v_ps0_y
  real*8 	:: Ti0_ps0_x, Ti0_ps0_y
  real*8 	:: Te0_ps0_x, Te0_ps0_y
  
  ! --- Neoclassical coefficients
  real*8 	:: tau_IC
  
  ! --- R,Z-coords and jacobians
  real*8 	:: x_g, x_s, x_t, x_ss, x_st, x_tt
  real*8 	:: y_g, y_s, y_t, y_ss, y_st, y_tt
  real*8 	:: BigR,  BigR_x
  real*8 	:: xjac, xjac_x, xjac_y
  
  ! --- Various
  real*8 	:: eps_cyl, theta, zeta
  
  ! --- Deltas
  real*8 	:: delta_g(n_var), delta_s(n_var), delta_t(n_var)
  
  ! --- Test functions
  real*8 	:: v,     v_x,	v_y,	 v_p,	  v_s,     v_t,     v_ss,     v_st,	v_tt,	v_xx,	v_yy,	v_xy, v_pp
  ! --- Variable 1
  real*8 	:: ps0,   ps0_x,   ps0_y,   ps0_p,   ps0_s,   ps0_t,	ps0_ss,   ps0_st,   ps0_tt, ps0_xx, ps0_xy, ps0_yy
  real*8 	:: psi,   psi_x,   psi_y,   psi_p,   psi_s,   psi_t,	psi_ss,   psi_st,   psi_tt, psi_xx, psi_yy, psi_xy, psi_pp
  ! --- Variable 2
  real*8 	:: u0,    u0_x,    u0_y,    u0_p,    u0_s,    u0_t, vv2
  real*8 	:: u,     u_x,     u_y,     u_p,     u_s,     u_t,  delta_u_x, delta_u_y
  ! --- Variable 3
  real*8 	:: zj0,   zj0_x,   zj0_y,   zj0_p,   zj0_s,   zj0_t
  real*8 	:: zj,    zj_x,    zj_y,    zj_p,    zj_s,    zj_t
  ! --- Variable 4
  real*8 	:: w0,    w0_x,    w0_y,    w0_p,    w0_s,    w0_t,	w0_ss,    w0_st,    w0_tt
  real*8 	:: w,     w_x,     w_y,     w_p,     w_s,     w_t,	w_ss,	  w_st,     w_tt
  ! --- Variable 5
  real*8 	:: r0,    r0_x,    r0_y,    r0_p,    r0_s,    r0_t,	r0_ss,    r0_st,    r0_tt,  r0_hat,  r0_x_hat,  r0_y_hat
  real*8 	:: rho,   rho_x,   rho_y,   rho_p,   rho_s,   rho_t,	rho_ss,   rho_st,   rho_tt
  ! --- Variable 6
  real*8 	:: Ti0,   Ti0_x,   Ti0_y,   Ti0_p,   Ti0_s,   Ti0_t,	Ti0_ss,   Ti0_st,   Ti0_tt, Ti0_pp, Ti0_xx, Ti0_xy, Ti0_yy
  real*8 	:: Ti,    Ti_x,    Ti_y,    Ti_p,    Ti_s,    Ti_t,	Ti_ss,    Ti_st,    Ti_tt,  Ti_pp,  Ti_xx,  Ti_xy,  Ti_yy
  ! --- Variable 7
  real*8 	:: Vpar0, Vpar0_x, Vpar0_y, Vpar0_p, Vpar0_s, Vpar0_t, Vpar0_ss, Vpar0_st, Vpar0_tt
  real*8 	:: Vpar,  Vpar_x,  Vpar_y,  Vpar_p,  Vpar_s,  Vpar_t,  Vpar_ss,  Vpar_st,  Vpar_tt
  ! --- Variable 8
  real*8 	:: Te0,   Te0_x,   Te0_y,   Te0_p,   Te0_s,   Te0_t,   Te0_ss,   Te0_st,   Te0_tt,  Te0_pp, Te0_xx, Te0_xy, Te0_yy
  real*8 	:: Te,    Te_x,    Te_y,    Te_p,    Te_s,    Te_t,    Te_ss,	 Te_st,    Te_tt,   Te_pp,  Te_xx,  Te_xy,  Te_yy
  ! --- Pressure
  real*8 	:: P0,    P0_x,	P0_y,	 P0_s,    P0_t,    P0_p
  
  ! --- Parallel gradients
  real*8 	:: BB2, BB2_psi
  real*8 	:: Bgrad_rho, Bgrad_rho_star, Bgrad_rho_k_star
  real*8 	:: Bgrad_Ti,  Bgrad_Ti_star,  Bgrad_Ti_k_star
  real*8 	:: Bgrad_Te,  Bgrad_Te_star,  Bgrad_Te_k_star
  
  ! --- Linearized equation terms
  real*8 	:: rhs_tmp(n_var),        rhs_k_tmp(n_var)
  real*8 	:: amat_tmp(n_var,n_var), amat_k_tmp(n_var,n_var), amat_n_tmp(n_var,n_var), amat_kn_tmp(n_var,n_var)
        
  ! --- !

  ! --- Declare variables as private for each thread (one module for each call to element_matrix)
  !$omp threadprivate(														   &
  !$omp 	zn,  dn_dpsi,  dn_dz,  dn_dpsi2,  dn_dz2,  dn_dpsi_dz,  dn_dpsi3,  dn_dpsi_dz2,  dn_dpsi2_dz,                      & 
  !$omp 	zTi, dTi_dpsi, dTi_dz, dTi_dpsi2, dTi_dz2, dTi_dpsi_dz, dTi_dpsi3, dTi_dpsi_dz2, dTi_dpsi2_dz,			   &
  !$omp 	zTe, dTe_dpsi, dTe_dz, dTe_dpsi2, dTe_dz2, dTe_dpsi_dz, dTe_dpsi3, dTe_dpsi_dz2, dTe_dpsi2_dz,			   &
  !$omp 	current_source, particle_source, heat_source_i, heat_source_e, V_source,					   &
  !$omp 	eta_Te,    deta_dTe, d2eta_d2Te,										   &
  !$omp 	visco_Te,  dvisco_dTe,												   &
  !$omp 	D_prof, 													   &
  !$omp 	Ki_prof, Ki_par, dKi_par,										           &
  !$omp 	Ke_prof, Ke_par, dKe_par,										           &
  !$omp 	eta_numm, visco_numm, visco_par_numm, D_perp_numm, Ki_perp_numm, Ki_par_num, Ke_perp_numm, Ke_par_num,		   &
  !$omp 	tau_IC,												   		   &
  !$omp 	v_ps0_x,   v_ps0_y,												   &
  !$omp 	Ti0_ps0_x, Ti0_ps0_y,												   &
  !$omp 	Te0_ps0_x, Te0_ps0_y,												   &
  !$omp 	x_g, x_s, x_t, x_ss, x_st, x_tt,										   &
  !$omp 	y_g, y_s, y_t, y_ss, y_st, y_tt,										   &
  !$omp 	BigR,  BigR_x,													   &
  !$omp 	xjac, xjac_x, xjac_y,												   &
  !$omp 	eps_cyl, theta, zeta,												   &
  !$omp 	delta_g, delta_s, delta_t,											   &
  !$omp 	v,     v_x,	v_y,	 v_p,	  v_s,     v_t,      v_ss,    v_st,	v_tt,	v_xx,	v_yy,	v_xy,  v_pp,	   &
  !$omp 	ps0,   ps0_x,	ps0_y,   ps0_p,   ps0_s,   ps0_t,    ps0_ss,  ps0_st,	ps0_tt, ps0_xx, ps0_xy, ps0_yy, 	   &
  !$omp 	psi,   psi_x,	psi_y,   psi_p,   psi_s,   psi_t,    psi_ss,  psi_st,	psi_tt, psi_xx, psi_yy, psi_xy, psi_pp,    &
  !$omp 	u0,    u0_x,	u0_y,	 u0_p,    u0_s,    u0_t, vv2,								   &
  !$omp 	u,     u_x,	u_y,	 u_p,	  u_s,     u_t,  delta_u_x, delta_u_y,						   &
  !$omp 	zj0,   zj0_x,	zj0_y,   zj0_p,   zj0_s,   zj0_t,								   &
  !$omp 	zj,    zj_x,	zj_y,	 zj_p,    zj_s,    zj_t,								   &
  !$omp 	w0,    w0_x,	w0_y,	 w0_p,    w0_s,    w0_t,     w0_ss,   w0_st,	w0_tt,					   &
  !$omp 	w,     w_x,	w_y,	 w_p,	  w_s,     w_t,      w_ss,    w_st,	w_tt,					   &
  !$omp 	r0,    r0_x,	r0_y,	 r0_p,    r0_s,    r0_t,     r0_ss,   r0_st,	r0_tt,  r0_hat,  r0_x_hat,  r0_y_hat,	   &
  !$omp 	rho,   rho_x,	rho_y,   rho_p,   rho_s,   rho_t,    rho_ss,  rho_st,	rho_tt, 				   &
  !$omp 	Ti0,   Ti0_x,	Ti0_y,   Ti0_p,   Ti0_s,   Ti0_t,    Ti0_ss,  Ti0_st,	Ti0_tt, Ti0_pp, Ti0_xx, Ti0_xy, Ti0_yy,    &
  !$omp 	Ti,    Ti_x,	Ti_y,	 Ti_p,    Ti_s,    Ti_t,     Ti_ss,   Ti_st,	Ti_tt,  Ti_pp,  Ti_xx,  Ti_xy,  Ti_yy,	   &
  !$omp 	Vpar0, Vpar0_x, Vpar0_y, Vpar0_p, Vpar0_s, Vpar0_t, Vpar0_ss, Vpar0_st, Vpar0_tt,				   &
  !$omp 	Vpar,  Vpar_x,  Vpar_y,  Vpar_p,  Vpar_s,  Vpar_t,  Vpar_ss,  Vpar_st,  Vpar_tt,				   &
  !$omp 	Te0,   Te0_x,	Te0_y,   Te0_p,   Te0_s,   Te0_t,   Te0_ss,   Te0_st,	Te0_tt,  Te0_pp, Te0_xx, Te0_xy, Te0_yy,   &
  !$omp 	Te,    Te_x,	Te_y,	 Te_p,    Te_s,    Te_t,    Te_ss,    Te_st,	Te_tt,   Te_pp,  Te_xx,  Te_xy,  Te_yy,    &
  !$omp 	P0,    P0_x, P0_y,    P0_s,    P0_t,	P0_p,									   &
  !$omp 	BB2, BB2_psi,                                                                                                      &
  !$omp 	Bgrad_rho, Bgrad_rho_star, Bgrad_rho_k_star,									   &
  !$omp 	Bgrad_Ti,  Bgrad_Ti_star,  Bgrad_Ti_k_star,									   &
  !$omp 	Bgrad_Te,  Bgrad_Te_star,  Bgrad_Te_k_star,									   &
  !$omp 	rhs_tmp, rhs_k_tmp,												   &
  !$omp 	amat_tmp, amat_k_tmp, amat_n_tmp, amat_kn_tmp)


end module equation_variables
