! Contains all the variables for the equations, to share with each equation routine
module equation_variables

  use parameters
  
  implicit none

    
  ! --- Profiles and sources
  real*8 	:: zn, dn_dpsi, dn_dz, dn_dpsi2, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz
  real*8 	:: zT, dT_dpsi, dT_dz, dT_dpsi2, dT_dz2, dT_dpsi_dz, dT_dpsi3, dT_dpsi_dz2, dT_dpsi2_dz
  real*8 	:: current_source, particle_source, heat_source
  
  ! --- Diffusivities
  real*8 	:: eta_T,    deta_dT, d2eta_d2T
  real*8 	:: visco_T,  dvisco_dT
  real*8 	:: D_prof, ZK_prof
  real*8 	:: K_par, dZK_par
  
  ! --- Hyper diffusivities
  real*8 	:: eta_numm, visco_numm, visco_par_numm, D_perp_numm, K_perp_numm, K_par_num 
  real*8 	:: T0_ps0_x, T0_ps0_y, v_ps0_x, v_ps0_y
  
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
  real*8 	:: T0,    T0_x,    T0_y,    T0_p,    T0_s,    T0_t,	T0_ss,    T0_st,    T0_tt, T0_pp, T0_xx, T0_xy, T0_yy
  real*8 	:: T,     T_x,     T_y,     T_p,     T_s,     T_t,	T_ss,	  T_st,     T_tt,  T_pp,  T_xx,  T_xy,  T_yy
  ! --- Variable 7
  real*8 	:: Vpar0, Vpar0_x, Vpar0_y, Vpar0_p, Vpar0_s, Vpar0_t, Vpar0_ss, Vpar0_st, Vpar0_tt
  real*8 	:: Vpar,  Vpar_x,  Vpar_y,  Vpar_p,  Vpar_s,  Vpar_t,  Vpar_ss,  Vpar_st,  Vpar_tt
  ! --- Pressure
  real*8 	:: P0,    P0_x,	P0_y,	 P0_s,    P0_t,    P0_p
  
  ! --- Parallel gradients
  real*8 	:: BB2, BB2_psi
  real*8 	:: Bgrad_rho, Bgrad_rho_star, Bgrad_rho_k_star
  real*8 	:: Bgrad_T,   Bgrad_T_star,   Bgrad_T_k_star
  
  ! --- Toroidally localised hyper-diffusivity
  real*8 	:: K_perp_numm2, K_par_num2  
  real*8 	:: T3, T3_ss, T3_st, T3_tt
 
  ! --- Linearized equation terms
  real*8 	:: rhs_tmp(n_var),        rhs_k_tmp(n_var)
  real*8 	:: amat_tmp(n_var,n_var), amat_k_tmp(n_var,n_var), amat_n_tmp(n_var,n_var), amat_kn_tmp(n_var,n_var)
        
  ! --- !

  ! --- Save variables for omp copies
  save 		zn, dn_dpsi, dn_dz, dn_dpsi2, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz
  save 		zT, dT_dpsi, dT_dz, dT_dpsi2, dT_dz2, dT_dpsi_dz, dT_dpsi3, dT_dpsi_dz2, dT_dpsi2_dz
  save 		current_source, particle_source, heat_source
  save 		eta_T,    deta_dT, d2eta_d2T
  save 		visco_T,  dvisco_dT
  save 		D_prof, ZK_prof
  save 		K_par, dZK_par
  save 		eta_numm, visco_numm, visco_par_numm, D_perp_numm, K_perp_numm, K_par_num
  save 		T0_ps0_x, T0_ps0_y, v_ps0_x, v_ps0_y
  save 		x_g, x_s, x_t, x_ss, x_st, x_tt
  save 		y_g, y_s, y_t, y_ss, y_st, y_tt
  save 		BigR,  BigR_x
  save 		xjac, xjac_x, xjac_y
  save 		eps_cyl, theta, zeta
  save 		delta_g, delta_s, delta_t
  save 		v,     v_x,     v_y,	 v_p,	  v_s,	   v_t,     v_ss,     v_st,     v_tt,   v_xx,   v_yy,   v_xy, v_pp
  save 		ps0,   ps0_x,   ps0_y,   ps0_p,   ps0_s,   ps0_t,   ps0_ss,   ps0_st,   ps0_tt, ps0_xx, ps0_xy, ps0_yy
  save 		psi,   psi_x,   psi_y,   psi_p,   psi_s,   psi_t,   psi_ss,   psi_st,   psi_tt, psi_xx, psi_yy, psi_xy, psi_pp
  save 		u0,    u0_x,    u0_y,	 u0_p,	  u0_s,    u0_t, vv2
  save 		u,     u_x,     u_y,	 u_p,	  u_s,	   u_t,  delta_u_x, delta_u_y
  save 		zj0,   zj0_x,   zj0_y,   zj0_p,   zj0_s,   zj0_t
  save 		zj,    zj_x,    zj_y,	 zj_p,	  zj_s,    zj_t
  save 		w0,    w0_x,    w0_y,	 w0_p,	  w0_s,    w0_t,    w0_ss,    w0_st,    w0_tt
  save 		w,     w_x,     w_y,	 w_p,	  w_s,	   w_t,     w_ss,     w_st,     w_tt
  save 		r0,    r0_x,    r0_y,	 r0_p,	  r0_s,    r0_t,    r0_ss,    r0_st,    r0_tt,  r0_hat,  r0_x_hat,  r0_y_hat
  save 		rho,   rho_x,   rho_y,   rho_p,   rho_s,   rho_t,   rho_ss,   rho_st,   rho_tt
  save 		T0,    T0_x,    T0_y,	 T0_p,	  T0_s,    T0_t,    T0_ss,    T0_st,    T0_tt, T0_pp, T0_xx, T0_xy, T0_yy
  save 		T,     T_x,     T_y,	 T_p,	  T_s,	   T_t,     T_ss,     T_st,     T_tt,  T_pp,  T_xx,  T_xy,  T_yy
  save 		Vpar0, Vpar0_x, Vpar0_y, Vpar0_p, Vpar0_s, Vpar0_t, Vpar0_ss, Vpar0_st, Vpar0_tt
  save 		Vpar,  Vpar_x,  Vpar_y,  Vpar_p,  Vpar_s,  Vpar_t,  Vpar_ss,  Vpar_st,  Vpar_tt
  save 		P0,    P0_x,    P0_y,	 P0_s,	  P0_t,    P0_p
  save 		BB2, BB2_psi
  save 		Bgrad_rho, Bgrad_rho_star, Bgrad_rho_k_star
  save 		Bgrad_T,   Bgrad_T_star,   Bgrad_T_k_star
  save 		K_perp_numm2, K_par_num2
  save 		T3, T3_ss, T3_st, T3_tt
  save 		rhs_tmp,  rhs_k_tmp
  save 		amat_tmp, amat_k_tmp, amat_n_tmp, amat_kn_tmp

  ! --- !

  ! --- Declare variables as private for each thread (one module for each call to element_matrix)
  !$omp threadprivate(														&
  !$omp 	zn, dn_dpsi, dn_dz, dn_dpsi2, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz,				&
  !$omp 	zT, dT_dpsi, dT_dz, dT_dpsi2, dT_dz2, dT_dpsi_dz, dT_dpsi3, dT_dpsi_dz2, dT_dpsi2_dz,				&
  !$omp 	current_source, particle_source, heat_source,									&
  !$omp 	eta_T,    deta_dT, d2eta_d2T,											&
  !$omp 	visco_T,  dvisco_dT,												&
  !$omp 	D_prof, ZK_prof,												&
  !$omp 	K_par, dZK_par,													&
  !$omp 	eta_numm, visco_numm, visco_par_numm, D_perp_numm, K_perp_numm, K_par_num,					&
  !$omp 	T0_ps0_x, T0_ps0_y, v_ps0_x, v_ps0_y,										&
  !$omp 	x_g, x_s, x_t, x_ss, x_st, x_tt,										&
  !$omp 	y_g, y_s, y_t, y_ss, y_st, y_tt,										&
  !$omp 	BigR,  BigR_x,													&
  !$omp 	xjac, xjac_x, xjac_y,												&
  !$omp 	eps_cyl, theta, zeta,												&
  !$omp 	delta_g, delta_s, delta_t,											&
  !$omp 	v,     v_x,	v_y,	 v_p,	  v_s,     v_t,     v_ss,     v_st,	v_tt,	v_xx,	v_yy,	v_xy, v_pp,	&
  !$omp 	ps0,   ps0_x,	ps0_y,   ps0_p,   ps0_s,   ps0_t,   ps0_ss,   ps0_st,	ps0_tt, ps0_xx, ps0_xy, ps0_yy,		&
  !$omp 	psi,   psi_x,	psi_y,   psi_p,   psi_s,   psi_t,   psi_ss,   psi_st,	psi_tt, psi_xx, psi_yy, psi_xy, psi_pp,	&
  !$omp 	u0,    u0_x,	u0_y,	 u0_p,    u0_s,    u0_t, vv2,								&
  !$omp 	u,     u_x,	u_y,	 u_p,	  u_s,     u_t,  delta_u_x, delta_u_y,						&
  !$omp 	zj0,   zj0_x,	zj0_y,   zj0_p,   zj0_s,   zj0_t,								&
  !$omp 	zj,    zj_x,	zj_y,	 zj_p,    zj_s,    zj_t,								&
  !$omp 	w0,    w0_x,	w0_y,	 w0_p,    w0_s,    w0_t,    w0_ss,    w0_st,	w0_tt,					&
  !$omp 	w,     w_x,	w_y,	 w_p,	  w_s,     w_t,     w_ss,     w_st,	w_tt,					&
  !$omp 	r0,    r0_x,	r0_y,	 r0_p,    r0_s,    r0_t,    r0_ss,    r0_st,	r0_tt,  r0_hat,  r0_x_hat,  r0_y_hat,	&
  !$omp 	rho,   rho_x,	rho_y,   rho_p,   rho_s,   rho_t,   rho_ss,   rho_st,	rho_tt,					&
  !$omp 	T0,    T0_x,	T0_y,	 T0_p,    T0_s,    T0_t,    T0_ss,    T0_st,	T0_tt, T0_pp, T0_xx, T0_xy, T0_yy,	&
  !$omp 	T,     T_x,	T_y,	 T_p,	  T_s,     T_t,     T_ss,     T_st,	T_tt,  T_pp,  T_xx,  T_xy,  T_yy,	&
  !$omp 	Vpar0, Vpar0_x, Vpar0_y, Vpar0_p, Vpar0_s, Vpar0_t, Vpar0_ss, Vpar0_st, Vpar0_tt,				&
  !$omp 	Vpar,  Vpar_x,  Vpar_y,  Vpar_p,  Vpar_s,  Vpar_t,  Vpar_ss,  Vpar_st,  Vpar_tt,				&
  !$omp 	P0,    P0_x,	P0_y,	 P0_s,    P0_t,    P0_p,								&
  !$omp 	BB2, BB2_psi,													&
  !$omp 	Bgrad_rho, Bgrad_rho_star, Bgrad_rho_k_star,									&
  !$omp 	Bgrad_T,   Bgrad_T_star,   Bgrad_T_k_star,									&
  !$omp 	K_perp_numm2, K_par_num2,											&
  !$omp 	T3, T3_ss, T3_st, T3_tt,											&
  !$omp 	rhs_tmp,  rhs_k_tmp,												&
  !$omp 	amat_tmp, amat_k_tmp, amat_n_tmp, amat_kn_tmp)


end module equation_variables
