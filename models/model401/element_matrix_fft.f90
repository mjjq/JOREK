! Contains all the variables for the equations, to share with each equation routine
module equation_variables

  use parameters
  
  implicit none

    
  ! --- Profiles and sources
  real*8     :: zn, dn_dpsi, dn_dz, dn_dpsi2, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz
  real*8     :: zT, dT_dpsi, dT_dz, dT_dpsi2, dT_dz2, dT_dpsi_dz, dT_dpsi3, dT_dpsi_dz2, dT_dpsi2_dz
  real*8     :: current_source, particle_source, heat_source
  
  ! --- Diffusivities
  real*8     :: eta_T,    deta_dT, d2eta_d2T
  real*8     :: visco_T,  dvisco_dT
  real*8     :: D_prof, ZK_prof
  real*8     :: K_par, dZK_par
  
  ! --- Hyper diffusivities
  real*8     :: eta_numm, visco_numm, visco_par_numm, D_perp_numm, K_perp_numm, K_par_num 
  real*8     :: T0_ps0_x, T0_ps0_y, v_ps0_x, v_ps0_y
  
  ! --- R,Z-coords and jacobians
  real*8     :: x_g, x_s, x_t, x_ss, x_st, x_tt
  real*8     :: y_g, y_s, y_t, y_ss, y_st, y_tt
  real*8     :: BigR,  BigR_x
  real*8     :: xjac, xjac_x, xjac_y
  
  ! --- Various
  real*8     :: eps_cyl, theta, zeta
  
  ! --- Deltas
  real*8     :: delta_g(n_var), delta_s(n_var), delta_t(n_var)
  
  ! --- Test functions
  real*8     :: v,     v_x,	v_y,	 v_p,	  v_s,     v_t,     v_ss,     v_st,	v_tt,	v_xx,	v_yy,	v_xy, v_pp
  ! --- Variable 1
  real*8     :: ps0,   ps0_x,	ps0_y,   ps0_p,   ps0_s,   ps0_t,   ps0_ss,   ps0_st,	ps0_tt, ps0_xx, ps0_xy, ps0_yy
  real*8     :: psi,   psi_x,	psi_y,   psi_p,   psi_s,   psi_t,   psi_ss,   psi_st,	psi_tt, psi_xx, psi_yy, psi_xy, psi_pp
  ! --- Variable 2
  real*8     :: u0,    u0_x,	u0_y,	 u0_p,    u0_s,    u0_t, vv2
  real*8     :: u,     u_x,	u_y,	 u_p,	  u_s,     u_t,  delta_u_x, delta_u_y
  ! --- Variable 3
  real*8     :: zj0,   zj0_x,	zj0_y,   zj0_p,   zj0_s,   zj0_t
  real*8     :: zj,    zj_x,	zj_y,	 zj_p,    zj_s,    zj_t
  ! --- Variable 4
  real*8     :: w0,    w0_x,	w0_y,	 w0_p,    w0_s,    w0_t,    w0_ss,    w0_st,	w0_tt
  real*8     :: w,     w_x,	w_y,	 w_p,	  w_s,     w_t,     w_ss,     w_st,	w_tt
  ! --- Variable 5
  real*8     :: r0,    r0_x,	r0_y,	 r0_p,    r0_s,    r0_t,    r0_ss,    r0_st,	r0_tt,  r0_hat,  r0_x_hat,  r0_y_hat
  real*8     :: rho,   rho_x,	rho_y,   rho_p,   rho_s,   rho_t,   rho_ss,   rho_st,	rho_tt
  ! --- Variable 6
  real*8     :: T0,    T0_x,	T0_y,	 T0_p,    T0_s,    T0_t,    T0_ss,    T0_st,	T0_tt, T0_pp, T0_xx, T0_xy, T0_yy
  real*8     :: T,     T_x,	T_y,	 T_p,	  T_s,     T_t,     T_ss,     T_st,	T_tt,  T_pp,  T_xx,  T_xy,  T_yy
  ! --- Variable 7
  real*8     :: Vpar0, Vpar0_x, Vpar0_y, Vpar0_p, Vpar0_s, Vpar0_t, Vpar0_ss, Vpar0_st, Vpar0_tt
  real*8     :: Vpar,  Vpar_x,  Vpar_y,  Vpar_p,  Vpar_s,  Vpar_t,  Vpar_ss,  Vpar_st,  Vpar_tt
  ! --- Pressure
  real*8     :: P0,    P0_x,	P0_y,	 P0_s,    P0_t,    P0_p
  
  ! --- Parallel gradients
  real*8     :: BB2, BB2_psi
  real*8     :: Bgrad_rho, Bgrad_rho_star, Bgrad_rho_k_star
  real*8     :: Bgrad_T,   Bgrad_T_star,   Bgrad_T_k_star
  
  ! --- Toroidally localised hyper-diffusivity
  real*8     :: K_perp_numm2, K_par_num2  
  real*8     :: T3, T3_ss, T3_st, T3_tt
 

  save zn, dn_dpsi, dn_dz, dn_dpsi2, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz
  save zT, dT_dpsi, dT_dz, dT_dpsi2, dT_dz2, dT_dpsi_dz, dT_dpsi3, dT_dpsi_dz2, dT_dpsi2_dz
  save current_source, particle_source, heat_source
  save eta_T,    deta_dT, d2eta_d2T
  save visco_T,  dvisco_dT
  save D_prof, ZK_prof
  save K_par, dZK_par
  save eta_numm, visco_numm, visco_par_numm, D_perp_numm, K_perp_numm, K_par_num
  save T0_ps0_x, T0_ps0_y, v_ps0_x, v_ps0_y
  save x_g, x_s, x_t, x_ss, x_st, x_tt
  save y_g, y_s, y_t, y_ss, y_st, y_tt
  save BigR,  BigR_x
  save xjac, xjac_x, xjac_y
  save eps_cyl, theta, zeta
  save delta_g, delta_s, delta_t
  save v,     v_x,     v_y,	v_p,	 v_s,	  v_t,     v_ss,     v_st,     v_tt,   v_xx,   v_yy,   v_xy, v_pp
  save ps0,   ps0_x,   ps0_y,   ps0_p,   ps0_s,   ps0_t,   ps0_ss,   ps0_st,   ps0_tt, ps0_xx, ps0_xy, ps0_yy
  save psi,   psi_x,   psi_y,   psi_p,   psi_s,   psi_t,   psi_ss,   psi_st,   psi_tt, psi_xx, psi_yy, psi_xy, psi_pp
  save u0,    u0_x,    u0_y,	u0_p,	 u0_s,    u0_t, vv2
  save u,     u_x,     u_y,	u_p,	 u_s,	  u_t,  delta_u_x, delta_u_y
  save zj0,   zj0_x,   zj0_y,   zj0_p,   zj0_s,   zj0_t
  save zj,    zj_x,    zj_y,	zj_p,	 zj_s,    zj_t
  save w0,    w0_x,    w0_y,	w0_p,	 w0_s,    w0_t,    w0_ss,    w0_st,    w0_tt
  save w,     w_x,     w_y,	w_p,	 w_s,	  w_t,     w_ss,     w_st,     w_tt
  save r0,    r0_x,    r0_y,	r0_p,	 r0_s,    r0_t,    r0_ss,    r0_st,    r0_tt,  r0_hat,  r0_x_hat,  r0_y_hat
  save rho,   rho_x,   rho_y,   rho_p,   rho_s,   rho_t,   rho_ss,   rho_st,   rho_tt
  save T0,    T0_x,    T0_y,	T0_p,	 T0_s,    T0_t,    T0_ss,    T0_st,    T0_tt, T0_pp, T0_xx, T0_xy, T0_yy
  save T,     T_x,     T_y,	T_p,	 T_s,	  T_t,     T_ss,     T_st,     T_tt,  T_pp,  T_xx,  T_xy,  T_yy
  save Vpar0, Vpar0_x, Vpar0_y, Vpar0_p, Vpar0_s, Vpar0_t, Vpar0_ss, Vpar0_st, Vpar0_tt
  save Vpar,  Vpar_x,  Vpar_y,  Vpar_p,  Vpar_s,  Vpar_t,  Vpar_ss,  Vpar_st,  Vpar_tt
  save P0,    P0_x,    P0_y,	P0_s,	 P0_t,    P0_p
  save BB2, BB2_psi
  save Bgrad_rho, Bgrad_rho_star, Bgrad_rho_k_star
  save Bgrad_T,   Bgrad_T_star,   Bgrad_T_k_star
  save K_perp_numm2, K_par_num2
  save T3, T3_ss, T3_st, T3_tt

  !$omp threadprivate(zn, dn_dpsi, dn_dz, dn_dpsi2, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz,&
  !$omp   	   zT, dT_dpsi, dT_dz, dT_dpsi2, dT_dz2, dT_dpsi_dz, dT_dpsi3, dT_dpsi_dz2, dT_dpsi2_dz,&
  !$omp   	   current_source, particle_source, heat_source,&
  !$omp   	   eta_T,    deta_dT, d2eta_d2T,&
  !$omp   	   visco_T,  dvisco_dT,&
  !$omp   	   D_prof, ZK_prof,&
  !$omp   	   K_par, dZK_par,&
  !$omp   	   eta_numm, visco_numm, visco_par_numm, D_perp_numm, K_perp_numm, K_par_num,&
  !$omp   	   T0_ps0_x, T0_ps0_y, v_ps0_x, v_ps0_y,&
  !$omp   	   x_g, x_s, x_t, x_ss, x_st, x_tt,&
  !$omp   	   y_g, y_s, y_t, y_ss, y_st, y_tt,&
  !$omp   	   BigR,  BigR_x,&
  !$omp   	   xjac, xjac_x, xjac_y,&
  !$omp   	   eps_cyl, theta, zeta,&
  !$omp   	   delta_g, delta_s, delta_t,&
  !$omp   	   v,	  v_x,     v_y,     v_p,     v_s,     v_t,     v_ss,	 v_st,     v_tt,   v_xx,   v_yy,   v_xy, v_pp,&
  !$omp   	   ps0,   ps0_x,   ps0_y,   ps0_p,   ps0_s,   ps0_t,   ps0_ss,   ps0_st,   ps0_tt, ps0_xx, ps0_xy, ps0_yy,&
  !$omp   	   psi,   psi_x,   psi_y,   psi_p,   psi_s,   psi_t,   psi_ss,   psi_st,   psi_tt, psi_xx, psi_yy, psi_xy, psi_pp,&
  !$omp   	   u0,    u0_x,    u0_y,    u0_p,    u0_s,    u0_t, vv2,&
  !$omp   	   u,	  u_x,     u_y,     u_p,     u_s,     u_t,  delta_u_x, delta_u_y,&
  !$omp   	   zj0,   zj0_x,   zj0_y,   zj0_p,   zj0_s,   zj0_t,&
  !$omp   	   zj,    zj_x,    zj_y,    zj_p,    zj_s,    zj_t,&
  !$omp   	   w0,    w0_x,    w0_y,    w0_p,    w0_s,    w0_t,    w0_ss,	 w0_st,    w0_tt,&
  !$omp   	   w,	  w_x,     w_y,     w_p,     w_s,     w_t,     w_ss,	 w_st,     w_tt,&
  !$omp   	   r0,    r0_x,    r0_y,    r0_p,    r0_s,    r0_t,    r0_ss,	 r0_st,    r0_tt,  r0_hat,  r0_x_hat,  r0_y_hat,&
  !$omp   	   rho,   rho_x,   rho_y,   rho_p,   rho_s,   rho_t,   rho_ss,   rho_st,   rho_tt,&
  !$omp   	   T0,    T0_x,    T0_y,    T0_p,    T0_s,    T0_t,    T0_ss,	 T0_st,    T0_tt, T0_pp, T0_xx, T0_xy, T0_yy,&
  !$omp   	   T,	  T_x,     T_y,     T_p,     T_s,     T_t,     T_ss,	 T_st,     T_tt,  T_pp,  T_xx,  T_xy,  T_yy,&
  !$omp   	   Vpar0, Vpar0_x, Vpar0_y, Vpar0_p, Vpar0_s, Vpar0_t, Vpar0_ss, Vpar0_st, Vpar0_tt,&
  !$omp   	   Vpar,  Vpar_x,  Vpar_y,  Vpar_p,  Vpar_s,  Vpar_t,  Vpar_ss,  Vpar_st,  Vpar_tt,&
  !$omp   	   P0,    P0_x,    P0_y,    P0_s,    P0_t,    P0_p,&
  !$omp   	   BB2, BB2_psi,&
  !$omp   	   Bgrad_rho, Bgrad_rho_star, Bgrad_rho_k_star,&
  !$omp   	   Bgrad_T,   Bgrad_T_star,   Bgrad_T_k_star,&
  !$omp   	   K_perp_numm2, K_par_num2,&
  !$omp   	   T3, T3_ss, T3_st, T3_tt)

end module equation_variables





module mod_elt_matrix_fft
contains







  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------ Calculates the matrix contribution of one element ---------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  subroutine element_matrix_fft(element, nodes, xpoint2, xcase2, minRad, R_axis, Z_axis, psi_axis, psi_bnd, Z_xpoint, ELM, RHS, tid)

    ! --- Modules
    use equation_variables
    use constants
    use parameters
    use data_structure
    use gauss
    use basis_at_gaussian
    use phys_module
    use tr_module 
    use profiles, only: interpolProf

    implicit none
    
    ! --- Structures
    type (type_element) 			:: element
    type (type_node)				:: nodes(n_vertex_max)
    type (type_surface_list)			:: flux_list

    ! --- Matrix elements and toroidal functions
    integer, intent(in) 			:: tid
    real*8, dimension (:,:)	     		:: ELM
    real*8, dimension (:)	     		:: RHS
    real*8, dimension(:,:,:) , pointer  	:: ELM_p, ELM_n, ELM_k, ELM_kn
    real*8, dimension(:,:)   , pointer  	:: RHS_p, RHS_k 
    
    ! --- Indexes
    integer    :: i_ij, ij_tmp
    integer    :: i_kl, kl_tmp
    integer    :: ms, mt
    integer    :: i_plane
    integer    :: i_order,  i_vertex, i_tor
    integer    :: j_order,  j_vertex, j_tor
    integer    :: n_tor_loop, n_tor_loop2
    integer    :: index_ij, index_kl
    
    ! --- Routine variables (Xpoint and axis)
    logical    :: xpoint2
    integer    :: xcase2
    real*8     :: minRad, R_axis, Z_axis, psi_axis, psi_bnd, Z_xpoint(2)
    
    ! --- Linearized equation terms
    real*8     :: wst
    real*8     :: rhs_tmp    (n_var)
    real*8     :: rhs_k_tmp  (n_var)
    real*8     :: amat_tmp   (n_var,n_var)
    real*8     :: amat_k_tmp (n_var,n_var)
    real*8     :: amat_n_tmp (n_var,n_var)
    real*8     :: amat_kn_tmp(n_var,n_var)
        
    ! --- Initialise rhs and lhs terms
    rhs_tmp  = 0.d0; rhs_k_tmp  = 0.d0
    amat_tmp = 0.d0; amat_k_tmp = 0.d0; amat_n_tmp = 0.d0; amat_kn_tmp = 0.d0
    
    ! --- Matrix elements pointers
    ELM_p  => thread_struct(tid)%ELM_p  ; ELM_p = 0.d0
    ELM_n  => thread_struct(tid)%ELM_n  ; ELM_n = 0.d0
    ELM_k  => thread_struct(tid)%ELM_k  ; ELM_k = 0.d0
    ELM_kn => thread_struct(tid)%ELM_kn ; ELM_kn = 0.d0
    RHS_p  => thread_struct(tid)%RHS_p  ; RHS_p = 0.d0
    RHS_k  => thread_struct(tid)%RHS_k  ; RHS_k = 0.d0
    
    ELM = 0.d0; RHS = 0.d0
        
    ! --- Take time evolution parameters from phys_module
    theta = time_evol_theta
    zeta  = time_evol_zeta

    ! for cylinder geometry : epscyl = eps
    eps_cyl = 1.d0

    ! --- If we're doing the fft, don't loop...
    if (n_tor .gt. 3) then
      n_tor_loop  = 1
      n_tor_loop2 = 1
    else
      n_tor_loop  = n_tor
      n_tor_loop2 = n_tor
    endif
    	      



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!! Begin integration loop over Gaussian integration points !!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do ms =1,n_gauss
      do mt =1,n_gauss
      
    	wst = wgauss(ms)*wgauss(mt)

    	call ELM_build_RZ_and_Jacobians(element, nodes, ms, mt)

    	do i_plane =1,n_plane

    	  call ELM_build_variables(element, nodes, ms, mt, i_plane)
          
    	  call ELM_build_diffusivities_and_sources(element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, Z_xpoint, i_plane)

    	  ! --- Now the equations, first the RHS
    	  do i_vertex =1,n_vertex_max

    	    do i_order =1,n_order+1

    	      do i_tor =1,n_tor_loop

    		! --- Index in the ELM matrix	    
    		if (n_tor .gt. 3) then
    		  index_ij =       n_var*(n_order+1)*(i_vertex-1) +       n_var*(i_order-1) + 1
    		else
    		  index_ij = n_tor*n_var*(n_order+1)*(i_vertex-1) + n_tor*n_var*(i_order-1) + i_tor
    		endif
		
	        call ELM_build_test_functions(element, nodes, ms, mt, i_plane, i_vertex, i_order, i_tor, &
					      v, v_s,  v_t,	   v_p,  v_x,  v_y,			 &
						 v_ss, v_tt, v_st, v_pp, v_xx, v_yy, v_xy		 )
		
		call ELM_main_rhs_1(rhs_tmp(1))
		call ELM_main_rhs_2(rhs_tmp(2))
		call ELM_main_rhs_3(rhs_tmp(3))
		call ELM_main_rhs_4(rhs_tmp(4))
		call ELM_main_rhs_5(rhs_tmp(5), rhs_k_tmp(5))
		call ELM_main_rhs_6(rhs_tmp(6), rhs_k_tmp(6))
		call ELM_main_rhs_7(rhs_tmp(7))
    		

    		! --- Fill up the matrix
    		if (n_tor .gt. 3) then
    		  do i_ij =1,n_var
		    ij_tmp = index_ij + (i_ij-1)*n_tor_loop
		    RHS_p(i_plane,ij_tmp) = RHS_p(i_plane,ij_tmp) + rhs_tmp  (i_ij) * wst
		    RHS_k(i_plane,ij_tmp) = RHS_k(i_plane,ij_tmp) + rhs_k_tmp(i_ij) * wst
		  enddo
    		else
    		  do i_ij =1,n_var
		    ij_tmp = index_ij + (i_ij-1)*n_tor_loop
    		    RHS(ij_tmp) = RHS(ij_tmp) + (rhs_tmp(i_ij) + rhs_k_tmp(i_ij)) * wst
		  enddo
    		endif

    		! --- And the LHS (linearised part)
    		do j_vertex =1,n_vertex_max

    		  do j_order =1,n_order+1

    		    do j_tor =1,n_tor_loop2

    		      ! --- Index in the ELM matrix
    		      if (n_tor .gt. 3) then
    			index_kl =       n_var*(n_order+1)*(j_vertex-1) +       n_var*(j_order-1) + 1
    		      else
    			index_kl = n_tor*n_var*(n_order+1)*(j_vertex-1) + n_tor*n_var*(j_order-1) + j_tor
    		      endif

		      call ELM_build_test_functions(element, nodes, ms, mt, i_plane, j_vertex, j_order, j_tor,  &
                                                    psi, psi_s,  psi_t,          psi_p,  psi_x,  psi_y,         &
				                         psi_ss, psi_tt, psi_st, psi_pp, psi_xx, psi_yy, psi_xy )
    		      
		      u    = psi; u_x    = psi_x; u_y    = psi_y; u_p    = psi_p; u_s    = psi_s; u_t    = psi_t
		      zj   = psi; zj_x   = psi_x; zj_y   = psi_y; zj_p   = psi_p; zj_s   = psi_s; zj_t   = psi_t
		      w    = psi; w_x    = psi_x; w_y    = psi_y; w_p    = psi_p; w_s    = psi_s; w_t    = psi_t
		      rho  = psi; rho_x  = psi_x; rho_y  = psi_y; rho_p  = psi_p; rho_s  = psi_s; rho_t  = psi_t
		      T    = psi; T_x    = psi_x; T_y    = psi_y; T_p    = psi_p; T_s    = psi_s; T_t    = psi_t
		      Vpar = psi; Vpar_x = psi_x; Vpar_y = psi_y; Vpar_p = psi_p; Vpar_s = psi_s; Vpar_t = psi_t
		      
		      w_ss    = psi_ss; w_tt	= psi_tt; w_st    = psi_st
   		      rho_ss  = psi_ss; rho_tt  = psi_tt; rho_st  = psi_st
   		      T_ss    = psi_ss; T_tt	= psi_tt; T_st    = psi_st
   		      Vpar_ss = psi_ss; Vpar_tt = psi_tt; Vpar_st = psi_st
                      
		      T_xx = psi_xx; T_yy = psi_yy; T_xy = psi_xy; T_pp = psi_pp
    		      
		      BB2_psi = 2.d0 * (psi_x * ps0_x + psi_y * ps0_y ) /BigR**2
   
    		      call ELM_main_lhs_1(amat_tmp(1,1), amat_tmp(1,2), amat_n_tmp(1,2), amat_tmp(1,3), amat_tmp(1,6))
    		      call ELM_main_lhs_2(amat_tmp(2,1), amat_tmp(2,2), amat_tmp(2,3), amat_n_tmp(2,3), amat_tmp(2,4), amat_tmp(2,5), amat_tmp(2,6))
    		      call ELM_main_lhs_3(amat_tmp(3,1), amat_tmp(3,3))
    		      call ELM_main_lhs_4(amat_tmp(4,2), amat_tmp(4,4))
    		      call ELM_main_lhs_5(amat_tmp(5,1), amat_k_tmp(5,1), amat_tmp(5,2), amat_tmp(5,5), amat_k_tmp(5,5), amat_n_tmp(5,5), amat_kn_tmp(5,5), amat_tmp(5,7), amat_n_tmp(5,7))
    		      call ELM_main_lhs_6(amat_tmp(6,1), amat_k_tmp(6,1), amat_tmp(6,2), amat_tmp(6,5), amat_tmp(6,6), amat_k_tmp(6,6), amat_n_tmp(6,6), amat_kn_tmp(6,6), amat_tmp(6,7), amat_n_tmp(6,7))
    		      call ELM_main_lhs_7(amat_tmp(7,1), amat_tmp(7,2), amat_tmp(7,5), amat_n_tmp(7,5), amat_tmp(7,6), amat_n_tmp(7,6), amat_tmp(7,7))
		      
    		      ! --- Fill up the matrix
    		      if (n_tor .gt. 3) then
    		  	do i_ij =1,n_var
		  	  ij_tmp = index_ij + (i_ij-1)*n_tor_loop
    		  	  do i_kl =1,n_var
		  	    kl_tmp = index_kl + (i_kl-1)*n_tor_loop
    			    ELM_p (i_plane,ij_tmp,kl_tmp) = ELM_p (i_plane,ij_tmp,kl_tmp) + wst * amat_tmp   (i_ij,i_kl)
    			    ELM_k (i_plane,ij_tmp,kl_tmp) = ELM_k (i_plane,ij_tmp,kl_tmp) + wst * amat_k_tmp (i_ij,i_kl)
    			    ELM_n (i_plane,ij_tmp,kl_tmp) = ELM_n (i_plane,ij_tmp,kl_tmp) + wst * amat_n_tmp (i_ij,i_kl)
    			    ELM_kn(i_plane,ij_tmp,kl_tmp) = ELM_kn(i_plane,ij_tmp,kl_tmp) + wst * amat_kn_tmp(i_ij,i_kl)
		  	  enddo
		  	enddo
    		      else
    		  	do i_ij =1,n_var
		  	  ij_tmp = index_ij + (i_ij-1)*n_tor_loop
    		  	  do i_kl =1,n_var
		  	    kl_tmp = index_kl + (i_kl-1)*n_tor_loop
    			    ELM(ij_tmp,kl_tmp) = ELM(ij_tmp,kl_tmp) + (amat_tmp(i_ij,i_kl) + amat_k_tmp(i_ij,i_kl) + amat_n_tmp(i_ij,i_kl) + amat_kn_tmp(i_ij,i_kl)) * wst
		  	  enddo
		  	enddo
    		      endif


    		    
    		    enddo ! inner n_tor_loop
    		  enddo   ! inner n_order+1
    		enddo	  ! inner n_vertex_max

    	      enddo	  ! outer n_tor_loop
    	    enddo	  ! outer n_order+1
    	  enddo 	  ! outer n_vertex_max

    	enddo		  ! n_plane

      enddo		  ! n_gauss
    enddo		  ! n_gauss




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!! Apply FFT !!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (n_tor .gt. 3) then
      call ELM_apply_fft(RHS, RHS_p, RHS_k, ELM, ELM_p, ELM_n, ELM_k, ELM_kn, tid)
    endif
    
    return
  end subroutine element_matrix_fft
















  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !-------------------------------------- Compute the RZ-coordinates and the Jacobians ------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  subroutine ELM_build_RZ_and_Jacobians(element, nodes, ms, mt)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_build_RZ_and_Jacobians

    ! --- Modules
    use parameters    
    use basis_at_gaussian
    use equation_variables
    use data_structure
    
    implicit none
    
    ! --- Routine Variables
    type (type_element) 	:: element
    type (type_node)    	:: nodes(n_vertex_max)
    integer            		:: ms, mt
    
    ! --- Internal Variables
    integer             	:: i, j
        
    ! --- Empty before integration
    x_g = 0.d0 ; x_s = 0.d0 ; x_t = 0.d0; x_st  = 0.d0; x_ss  = 0.d0; x_tt  = 0.d0
    y_g = 0.d0 ; y_s = 0.d0 ; y_t = 0.d0; y_st  = 0.d0; y_ss  = 0.d0; y_tt  = 0.d0

    ! --- Integrate
    do i=1,n_vertex_max
      do j=1,n_order+1

    	x_g  = x_g  + nodes(i)%x(j,1) * element%size(i,j) * H   (i,j,ms,mt)
    	x_s  = x_s  + nodes(i)%x(j,1) * element%size(i,j) * H_s (i,j,ms,mt)
    	x_t  = x_t  + nodes(i)%x(j,1) * element%size(i,j) * H_t (i,j,ms,mt)

    	x_ss = x_ss + nodes(i)%x(j,1) * element%size(i,j) * H_ss(i,j,ms,mt)
    	x_st = x_st + nodes(i)%x(j,1) * element%size(i,j) * H_st(i,j,ms,mt)
    	x_tt = x_tt + nodes(i)%x(j,1) * element%size(i,j) * H_tt(i,j,ms,mt)

    	y_g  = y_g  + nodes(i)%x(j,2) * element%size(i,j) * H   (i,j,ms,mt)
    	y_s  = y_s  + nodes(i)%x(j,2) * element%size(i,j) * H_s (i,j,ms,mt)
    	y_t  = y_t  + nodes(i)%x(j,2) * element%size(i,j) * H_t (i,j,ms,mt)

    	y_ss = y_ss + nodes(i)%x(j,2) * element%size(i,j) * H_ss(i,j,ms,mt)
    	y_st = y_st + nodes(i)%x(j,2) * element%size(i,j) * H_st(i,j,ms,mt)
    	y_tt = y_tt + nodes(i)%x(j,2) * element%size(i,j) * H_tt(i,j,ms,mt)
      
      enddo
    enddo
    
    BigR    = x_g
    BigR_x  = 1.d0
    
    ! --- Jacobians
    xjac    = x_s*y_t - x_t*y_s

    xjac_x  = (x_ss* y_t**2 - y_ss*x_t*y_t - 2.d0*x_st*y_s*y_t   &	 
    	     + y_st*(x_s*y_t + x_t*y_s) 			 &
    	     + x_tt* y_s**2 - y_tt*x_s*y_s) / xjac
    	   
    xjac_y  = (y_tt* x_s**2 - x_tt*y_s*x_s - 2.d0*y_st*x_t*x_s   &	 
    	     + x_st*(y_t*x_s + y_s*x_t) 			 &
    	     + y_ss* x_t**2 - x_ss*y_t*x_t) / xjac

    return

  end subroutine ELM_build_RZ_and_Jacobians








  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !----------------------------------------- Compute the variables for the equations --------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  subroutine ELM_build_variables(element, nodes, ms, mt, i_plane)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_build_variables

    ! --- Modules
    use parameters    
    use basis_at_gaussian
    use equation_variables
    use data_structure
    use phys_module
    
    implicit none
    
    ! --- Routine variables
    type (type_element) 	:: element
    type (type_node)    	:: nodes(n_vertex_max)
    integer             	:: ms, mt, i_plane
    
    ! --- Internal variables
    integer             	:: i, j, k, i_tor
        
    ! --- Empty before integration
    ps0   = 0.d0; ps0_p   = 0.d0; ps0_s   = 0.d0; ps0_t   = 0.d0
    u0    = 0.d0; u0_p    = 0.d0; u0_s    = 0.d0; u0_t    = 0.d0
    zj0   = 0.d0; zj0_p   = 0.d0; zj0_s   = 0.d0; zj0_t   = 0.d0
    w0    = 0.d0; w0_p    = 0.d0; w0_s    = 0.d0; w0_t    = 0.d0; w0_ss    = 0.d0; w0_tt    = 0.d0; w0_st    = 0.d0
    r0    = 0.d0; r0_p    = 0.d0; r0_s    = 0.d0; r0_t    = 0.d0; r0_ss    = 0.d0; r0_tt    = 0.d0; r0_st    = 0.d0
    T0    = 0.d0; T0_p    = 0.d0; T0_s    = 0.d0; T0_t    = 0.d0; T0_ss    = 0.d0; T0_tt    = 0.d0; T0_st    = 0.d0; T0_pp = 0.d0
    T3    = 0.d0; T3_ss   = 0.d0; T3_tt   = 0.d0; T3_st   = 0.d0
    Vpar0 = 0.d0; Vpar0_p = 0.d0; Vpar0_s = 0.d0; Vpar0_t = 0.d0; Vpar0_ss = 0.d0; Vpar0_tt = 0.d0; Vpar0_st = 0.d0
    delta_g = 0.d0 ; delta_s = 0.d0 ; delta_t = 0.d0

    ! --- Integrate
    do i =1,n_vertex_max
      do j=1,n_order+1
    	do i_tor =1,n_tor

          ! --- Variable 1
    	  ps0	       = ps0	    + nodes(i)%values(i_tor,j,1) * element%size(i,j) * H   (i,j,ms,mt) * HZ   (i_tor,i_plane)
    	  ps0_p        = ps0_p      + nodes(i)%values(i_tor,j,1) * element%size(i,j) * H   (i,j,ms,mt) * HZ_p (i_tor,i_plane)
    	  ps0_s        = ps0_s      + nodes(i)%values(i_tor,j,1) * element%size(i,j) * H_s (i,j,ms,mt) * HZ   (i_tor,i_plane)
    	  ps0_t        = ps0_t      + nodes(i)%values(i_tor,j,1) * element%size(i,j) * H_t (i,j,ms,mt) * HZ   (i_tor,i_plane)

          ! --- Variable 2
    	  u0	       = u0	    + nodes(i)%values(i_tor,j,2) * element%size(i,j) * H   (i,j,ms,mt) * HZ   (i_tor,i_plane)
    	  u0_p         = u0_p	    + nodes(i)%values(i_tor,j,2) * element%size(i,j) * H   (i,j,ms,mt) * HZ_p (i_tor,i_plane)
    	  u0_s         = u0_s	    + nodes(i)%values(i_tor,j,2) * element%size(i,j) * H_s (i,j,ms,mt) * HZ   (i_tor,i_plane)
    	  u0_t         = u0_t	    + nodes(i)%values(i_tor,j,2) * element%size(i,j) * H_t (i,j,ms,mt) * HZ   (i_tor,i_plane)

          ! --- Variable 3
    	  zj0	       = zj0	    + nodes(i)%values(i_tor,j,3) * element%size(i,j) * H   (i,j,ms,mt) * HZ   (i_tor,i_plane)
    	  zj0_p        = zj0_p      + nodes(i)%values(i_tor,j,3) * element%size(i,j) * H   (i,j,ms,mt) * HZ_p (i_tor,i_plane)
    	  zj0_s        = zj0_s      + nodes(i)%values(i_tor,j,3) * element%size(i,j) * H_s (i,j,ms,mt) * HZ   (i_tor,i_plane)
    	  zj0_t        = zj0_t      + nodes(i)%values(i_tor,j,3) * element%size(i,j) * H_t (i,j,ms,mt) * HZ   (i_tor,i_plane)

          ! --- Variable 4
    	  w0	       = w0	    + nodes(i)%values(i_tor,j,4) * element%size(i,j) * H   (i,j,ms,mt) * HZ   (i_tor,i_plane)
    	  w0_p         = w0_p	    + nodes(i)%values(i_tor,j,4) * element%size(i,j) * H   (i,j,ms,mt) * HZ_p (i_tor,i_plane)
    	  w0_s         = w0_s	    + nodes(i)%values(i_tor,j,4) * element%size(i,j) * H_s (i,j,ms,mt) * HZ   (i_tor,i_plane)
    	  w0_t         = w0_t	    + nodes(i)%values(i_tor,j,4) * element%size(i,j) * H_t (i,j,ms,mt) * HZ   (i_tor,i_plane)
    	  w0_ss        = w0_ss      + nodes(i)%values(i_tor,j,4) * element%size(i,j) * H_ss(i,j,ms,mt) * HZ   (i_tor,i_plane)
    	  w0_tt        = w0_tt      + nodes(i)%values(i_tor,j,4) * element%size(i,j) * H_tt(i,j,ms,mt) * HZ   (i_tor,i_plane)
    	  w0_st        = w0_st      + nodes(i)%values(i_tor,j,4) * element%size(i,j) * H_st(i,j,ms,mt) * HZ   (i_tor,i_plane)

          ! --- Variable 5
    	  r0	       = r0	    + nodes(i)%values(i_tor,j,5) * element%size(i,j) * H   (i,j,ms,mt) * HZ   (i_tor,i_plane)
    	  r0_p         = r0_p	    + nodes(i)%values(i_tor,j,5) * element%size(i,j) * H   (i,j,ms,mt) * HZ_p (i_tor,i_plane)
    	  r0_s         = r0_s	    + nodes(i)%values(i_tor,j,5) * element%size(i,j) * H_s (i,j,ms,mt) * HZ   (i_tor,i_plane)
    	  r0_t         = r0_t	    + nodes(i)%values(i_tor,j,5) * element%size(i,j) * H_t (i,j,ms,mt) * HZ   (i_tor,i_plane)
    	  r0_ss        = r0_ss      + nodes(i)%values(i_tor,j,5) * element%size(i,j) * H_ss(i,j,ms,mt) * HZ   (i_tor,i_plane)
    	  r0_tt        = r0_tt      + nodes(i)%values(i_tor,j,5) * element%size(i,j) * H_tt(i,j,ms,mt) * HZ   (i_tor,i_plane)
    	  r0_st        = r0_st      + nodes(i)%values(i_tor,j,5) * element%size(i,j) * H_st(i,j,ms,mt) * HZ   (i_tor,i_plane)

          ! --- Variable 6
    	  T0	       = T0	    + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H   (i,j,ms,mt) * HZ   (i_tor,i_plane)
    	  T0_p         = T0_p	    + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H   (i,j,ms,mt) * HZ_p (i_tor,i_plane)
    	  T0_s         = T0_s	    + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H_s (i,j,ms,mt) * HZ   (i_tor,i_plane)
    	  T0_t         = T0_t	    + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H_t (i,j,ms,mt) * HZ   (i_tor,i_plane)
    	  T0_ss        = T0_ss      + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H_ss(i,j,ms,mt) * HZ   (i_tor,i_plane)
    	  T0_tt        = T0_tt      + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H_tt(i,j,ms,mt) * HZ   (i_tor,i_plane)
    	  T0_st        = T0_st      + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H_st(i,j,ms,mt) * HZ   (i_tor,i_plane)
    	  T0_pp        = T0_pp      + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H   (i,j,ms,mt) * HZ_pp(i_tor,i_plane)
    
          ! --- Toroidally localised hyper-parallel-conductivity
    	  if(i_tor .ne. 1) then
    	    T3         = T3	    + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H   (i,j,ms,mt) * HZ   (i_tor,i_plane)
    	    T3_ss      = T3_ss      + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H_ss(i,j,ms,mt) * HZ   (i_tor,i_plane)
    	    T3_tt      = T3_tt      + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H_tt(i,j,ms,mt) * HZ   (i_tor,i_plane)
    	    T3_st      = T3_st      + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H_st(i,j,ms,mt) * HZ   (i_tor,i_plane)
    	  endif

          ! --- Variable 7
    	  Vpar0        = Vpar0      + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H   (i,j,ms,mt) * HZ   (i_tor,i_plane)
    	  Vpar0_p      = Vpar0_p    + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H   (i,j,ms,mt) * HZ_p (i_tor,i_plane)
    	  Vpar0_s      = Vpar0_s    + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H_s (i,j,ms,mt) * HZ   (i_tor,i_plane)
    	  Vpar0_t      = Vpar0_t    + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H_t (i,j,ms,mt) * HZ   (i_tor,i_plane)
    	  Vpar0_ss     = Vpar0_ss   + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H_ss(i,j,ms,mt) * HZ   (i_tor,i_plane)
    	  Vpar0_tt     = Vpar0_tt   + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H_tt(i,j,ms,mt) * HZ   (i_tor,i_plane)
    	  Vpar0_st     = Vpar0_st   + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H_st(i,j,ms,mt) * HZ   (i_tor,i_plane)

          ! --- Deltas
    	  do k=1,n_var
    	    delta_g(k) = delta_g(k) + nodes(i)%deltas(i_tor,j,k) * element%size(i,j) * H   (i,j,ms,mt) * HZ   (i_tor,i_plane)
    	    delta_s(k) = delta_s(k) + nodes(i)%deltas(i_tor,j,k) * element%size(i,j) * H_s (i,j,ms,mt) * HZ   (i_tor,i_plane)
    	    delta_t(k) = delta_t(k) + nodes(i)%deltas(i_tor,j,k) * element%size(i,j) * H_t (i,j,ms,mt) * HZ   (i_tor,i_plane)
    	  enddo 		  
    
    	enddo
      enddo
    enddo

    ! --- Variable 1
    ps0_x    = (   y_t * ps0_s - y_s * ps0_t ) / xjac
    ps0_y    = ( - x_t * ps0_s + x_s * ps0_t ) / xjac
    
    ! --- Variable 2
    u0_x     = (   y_t * u0_s - y_s * u0_t ) / xjac
    u0_y     = ( - x_t * u0_s + x_s * u0_t ) / xjac
    vv2      = BigR**2 *  ( u0_x * u0_x + u0_y *u0_y  )
    
    ! --- Variable 3
    zj0_x    = (   y_t * zj0_s - y_s * zj0_t ) / xjac
    zj0_y    = ( - x_t * zj0_s + x_s * zj0_t ) / xjac
    
    ! --- Variable 4
    w0_x     = (   y_t * w0_s - y_s * w0_t ) / xjac
    w0_y     = ( - x_t * w0_s + x_s * w0_t ) / xjac
    
    ! --- Variable 5
    r0_x     = (   y_t * r0_s - y_s * r0_t ) / xjac
    r0_y     = ( - x_t * r0_s + x_s * r0_t ) / xjac
    r0_hat   = BigR**2 * r0
    r0_x_hat = 2.d0 * BigR * BigR_x  * r0 + BigR**2 * r0_x
    r0_y_hat = BigR**2 * r0_y
    
    ! --- Variable 6
    T0_x     = (   y_t * T0_s  - y_s * T0_t ) / xjac
    T0_y     = ( - x_t * T0_s  + x_s * T0_t ) / xjac
    T0_xx    = (T0_ss * y_t**2 - 2.d0*T0_st * y_s*y_t + T0_tt * y_s**2  &	      
    	  	+ T0_s * (y_st*y_t - y_tt*y_s ) 			&    
    	  	+ T0_t * (y_st*y_s - y_ss*y_t ) )      / xjac**2	&	  
    	      - xjac_x * (T0_s * y_t - T0_t * y_s)     / xjac**2
    T0_yy    = (T0_ss * x_t**2 - 2.d0*T0_st * x_s*x_t + T0_tt * x_s**2  &	      
    	  	+ T0_s * (x_st*x_t - x_tt*x_s ) 			&    
    	  	+ T0_t * (x_st*x_s - x_ss*x_t ) )      / xjac**2	&	  
    	      - xjac_y * (- T0_s * x_t + T0_t * x_s )  / xjac**2
    T0_xy    = (- T0_ss * y_t*x_t - T0_tt * x_s*y_s			&
    	  	+ T0_st * (y_s*x_t  + y_t*x_s  )			&	 
    	  	- T0_s  * (x_st*y_t - x_tt*y_s )			&    
    	  	- T0_t  * (x_st*y_s - x_ss*y_t )  )    / xjac**2	&	  
    	      - xjac_x * (- T0_s * x_t + T0_t * x_s )  / xjac**2
    T0_ps0_x = T0_xx * ps0_y - T0_xy * ps0_x + T0_x * ps0_xy - T0_y * ps0_xx
    T0_ps0_y = T0_xy * ps0_y - T0_yy * ps0_x + T0_x * ps0_yy - T0_y * ps0_xy
    
    ! --- Variable 7
    Vpar0_x  = (   y_t * Vpar0_s - y_s * Vpar0_t ) / xjac
    Vpar0_y  = ( - x_t * Vpar0_s + x_s * Vpar0_t ) / xjac
    
    ! --- Deltas
    delta_u_x= (   y_t * delta_s(2) - y_s * delta_t(2) ) / xjac
    delta_u_y= ( - x_t * delta_s(2) + x_s * delta_t(2) ) / xjac
    
    ! --- Pressure
    P0       = r0   * T0
    P0_x     = r0_x * T0 + r0 * T0_x
    P0_y     = r0_y * T0 + r0 * T0_y
    P0_s     = r0_s * T0 + r0 * T0_s
    P0_t     = r0_t * T0 + r0 * T0_t
    P0_p     = r0_p * T0 + r0 * T0_p
    
    ! --- Magnetic field amplitude (squared)
    BB2       = (F0*F0 + ps0_x * ps0_x + ps0_y * ps0_y )/BigR**2
    
    return

  end subroutine ELM_build_variables







  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------- Compute the diffusivities and source ---------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  subroutine ELM_build_diffusivities_and_sources(element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, Z_xpoint, i_plane)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_build_diffusivities_and_sources

    ! --- Modules
    use parameters    
    use basis_at_gaussian
    use phys_module
    use equation_variables
    use data_structure
    
    implicit none
    
    ! --- Routine variables
    type (type_element) 	:: element
    type (type_node)    	:: nodes(n_vertex_max)
    logical             	:: xpoint2
    integer             	:: xcase2
    integer             	:: i_plane
    real*8              	:: R_axis, Z_axis, psi_axis, psi_bnd, Z_xpoint(2)
    
    ! --- Internal variables
    integer             	:: id
    real*8              	:: psi_norm, psi_D
    real*8              	:: atn_D, datn_D, atn_D_n, pol_D, dpol_D, D_min
    real*8              	:: prof(1:2),Diff(1:2,1:10)
        
    ! -------------------------------------
    ! --- Temperature dependent resistivity
    ! -------------------------------------
    if ( eta_T_dependent ) then
      if ( T0 .lt. T_1 ) then
    	eta_T	  =   eta   * (T_1/T_0)**(-1.5d0)
    	deta_dT   = - eta   * (1.5d0)  * T_1**(-2.5d0) * T_0**(1.5d0)
    	d2eta_d2T =   eta   * (3.75d0) * T_1**(-3.5d0) * T_0**(1.5d0)
      else
    	eta_T	  =   eta   * (T0 /T_0)**(-1.5d0)
    	deta_dT   = - eta   * (1.5d0)  * T0 **(-2.5d0) * T_0**(1.5d0)
    	d2eta_d2T =   eta   * (3.75d0) * T0 **(-3.5d0) * T_0**(1.5d0)
      endif
    else
      eta_T	= eta
      deta_dT	= 0.d0
      d2eta_d2T = 0.d0
    end if
    
    ! -----------------------------------
    ! --- Temperature dependent viscosity
    ! -----------------------------------
    if ( visco_T_dependent ) then
      if ( T0 .lt. T_1 ) then
    	visco_T   =   visco * (T_1/T_0)**(-1.5d0)
    	dvisco_dT = - visco * (1.5d0)  * T_1**(-2.5d0) * T_0**(1.5d0)
      else
    	visco_T   =   visco * (T0 /T_0)**(-1.5d0)
    	dvisco_dT = - visco * (1.5d0)  * T0 **(-2.5d0) * T_0**(1.5d0)
      endif
    else
      visco_T	= visco
      dvisco_dT = 0.d0
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
    	Diff(1,id) = rho_coef(id-6) ; Diff(2,id) = T_coef(id-6)
      else
    	Diff(1,id) = D_perp(id)     ; Diff(2,id) = ZK_perp(id)
      endif
    enddo      
    if (Diff(1,10) .eq. 1.d0) then 
      Diff(1,4) = rho_coef(4) ; Diff(1,5) = rho_coef(5) 
    endif
    if (Diff(2,10) .eq. 1.d0) then 
      Diff(2,4) = T_coef(4)   ; Diff(2,5) = T_coef(5) 
    endif
    
    ! --- Build profiles
    do id = 1,2
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
      if ((xcase2 .ne. 2) .and. (id .eq. 1))  psi_D = psi_D * (0.5d0 - 0.5d0 * tanh((Z_xpoint(1)-y_g)/0.1d0))
      if ((xcase2 .ne. 1) .and. (id .eq. 1))  psi_D = psi_D * (0.5d0 - 0.5d0 * tanh((y_g-Z_xpoint(2))/0.1d0))
      atn_D    = 0.5d0 - 0.5d0 * tanh((psi_D-Diff(id,5))/Diff(id,4))
      datn_D   =       - 0.5d0 / cosh((psi_D-Diff(id,5))/Diff(id,4))**2.d0 /(Diff(id,4)*(psi_bnd - psi_axis))
      pol_D    = 1 + Diff(id,7)*psi_D	 + Diff(id,8)*psi_D**2.d0      + Diff(id,9)*psi_D**3.d0
      dpol_D   =    (Diff(id,7)       + 2.d0*Diff(id,8)*psi_D	  + 3.d0*Diff(id,9)*psi_D**2.d0)/(psi_bnd - psi_axis)
      D_min    = 1.d0/( -(1+Diff(id,7)*Diff(id,5)+Diff(id,8)*Diff(id,5)**2.d0+Diff(id,9)*Diff(id,5)**3.d0) * 0.5d0/(Diff(id,4)*(psi_bnd - psi_axis))&
    	  	 + 0.5d0 * (Diff(id,7)     + 2.d0*Diff(id,8)*Diff(id,5)+ 3.d0*Diff(id,9)*Diff(id,5)**2.d0)/(psi_bnd - psi_axis) )

      prof(id) = (1.d0-Diff(id,10)) * ( Diff(id,1) * (1.d0-Diff(id,2)+Diff(id,2)*(0.5d0 - 0.5d0 * tanh((psi_norm-Diff(id,5))/Diff(id,4)))) &
    	  			      + Diff(id,6) * (0.5d0 - 0.5d0 * tanh((-psi_norm+Diff(id,5)+Diff(id,3))/Diff(id,4)))) &
    	  	      +Diff(id,10)  * ( Diff(id,1) / (dpol_D*atn_D + pol_D*datn_D) / D_min ) &
    	  			    * (1 + Diff(id,6) - Diff(id,6) * tanh(-(psi_norm-(1+4*Diff(id,4)))/Diff(id,4)))   !higher Kperp in SOL
    enddo
    
    ! --- Allocate profiles to corresponding names
    D_prof  = prof(1)
    ZK_prof = prof(2) 
    
    ! --- Avoid negative density 
    if ( r0 .lt. rho_1 ) then
      D_prof = D_prof * 1.d2
    endif	
    if ( T0 .lt. T_1 ) then
      ZK_prof = ZK_prof * 1.d0
    endif

    ! -----------------------------------------------------
    ! --- Parallel conductivity profiles (Braginskii model)
    ! -----------------------------------------------------
    if ( T0 .lt. T_1 ) then
      K_par   = ZK_par  	 * T_1**(2.5d0) 
      dZK_par = ZK_par * (2.5d0) * T_1**(1.5d0)
    else
      K_par   = ZK_par  	 * T0 **(2.5d0) 
      dZK_par = ZK_par * (2.5d0) * T0 **(1.5d0)
    endif
    
    ! -------------------------
    ! --- Hyper diffusivitities
    ! -------------------------
    eta_numm	   = eta_num			   ! hyper-resistivity
    visco_numm     = visco_num  		   ! hyper-viscosity
    visco_par_numm = visco_par_num		   ! hyper-viscosity
    D_perp_numm    = D_perp_num 		   ! hyper-diffusivity
    K_perp_numm    = 0.d0!ZK_perp_num		   ! hyper-conductivity
    K_perp_numm2   = ZK_perp_num		   ! hyper-conductivity
    K_par_num      = 0.d0!1.d-10		   ! hyper-parallel-conductivity
    K_par_num2     = 0.d0!1.d-10		   ! hyper-parallel-conductivity
    
    if (psi_norm .lt. 0.4d0) eta_numm	= eta_numm   * 1.d2
    if (psi_norm .lt. 0.4d0) visco_numm = visco_numm * 1.d2
    if (psi_norm .lt. 0.2d0) visco_numm = visco_numm * 1.d2
    if (psi_norm .lt. 0.2d0) eta_numm	= eta_numm   * 1.d2

    ! -------------------------------------------------------------------
    ! --- Heating, current and particle source (the same for all i_plane)
    ! -------------------------------------------------------------------
    if (i_plane .eq. 1) then
      call current(xpoint2, xcase2, x_g,y_g, Z_xpoint, ps0,psi_axis,psi_bnd,current_source)

      !call sources(xpoint2, xcase2, y_g, Z_xpoint, ps0,psi_axis,psi_bnd,particle_source,heat_source_i,heat_source_e)
      ! --- New source profile: source with exactly the same profile as the initial equilibirum profiles.
      !call density(    xpoint2, xcase2, y_g, Z_xpoint, ps0,psi_axis,psi_bnd, &
      !	  	       zn,dn_dpsi,  dn_dz, dn_dpsi2, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz)

      !call temperature(xpoint2, xcase2, y_g, Z_xpoint, ps0,psi_axis,psi_bnd, &
      !	  	       zT,dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2,dT_dpsi2_dz)

      !particle_source = particlesource * ( zn - r0 )
      !heat_source    = heatsource     * ( zT - T0 )

      !particle_source = particle_source * ( 0.5d0 - 0.5d0 * tanh((psi_norm-0.99)/0.005) )
      !heat_source    = heat_source	* ( 0.5d0 - 0.5d0 * tanh((psi_norm-0.99)/0.005) )
      particle_source = particlesource * (0.5d0 - 0.5d0 * tanh((psi_norm-0.5)/0.005) )
      heat_source     = heatsource     * (0.5d0 - 0.5d0 * tanh((psi_norm-0.5)/0.005) ) 
    endif
    
    return

  end subroutine ELM_build_diffusivities_and_sources










  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !-------------------------------- Compute the test functions (ie. the Bezier polynomials) -------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  subroutine ELM_build_test_functions(element, nodes, ms, mt, i_plane, i_vertex, i_order, i_tor, &
                                      vv, vv_s,  vv_t,         vv_p,  vv_x,  vv_y,               &
				          vv_ss, vv_tt, vv_st, vv_pp, vv_xx, vv_yy, vv_xy        )
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_build_test_functions

    ! --- Modules
    use parameters    
    use basis_at_gaussian
    use equation_variables
    use data_structure
    
    implicit none
    
    ! --- Routine variables
    type (type_element) 	:: element
    type (type_node)    	:: nodes(n_vertex_max)
    integer             	:: ms, mt, i_plane, i_vertex, i_order, i_tor
    real*8              	:: vv, vv_s,  vv_t,         vv_p,  vv_x,  vv_y
    real*8              	::     vv_ss, vv_tt, vv_st, vv_pp, vv_xx, vv_yy, vv_xy
    
    ! --- Internal variables
    real*8, dimension(n_tor,n_plane) :: HHZ, HHZ_p, HHZ_pp
    
    ! --- Toroidal functions		
    if (n_tor .gt. 3) then
      HHZ   (i_tor,i_plane) = 1.d0
      HHZ_p (i_tor,i_plane) = 1.d0
      HHZ_pp(i_tor,i_plane) = 1.d0
    else
      HHZ   (i_tor,i_plane) = HZ   (i_tor,i_plane)
      HHZ_p (i_tor,i_plane) = HZ_p (i_tor,i_plane)
      HHZ_pp(i_tor,i_plane) = HZ_pp(i_tor,i_plane)
    endif
    
    ! --- Test functions		
    vv	  = H   (i_vertex,i_order,ms,mt) * element%size(i_vertex,i_order) * HHZ   (i_tor,i_plane)
    vv_s  = H_s (i_vertex,i_order,ms,mt) * element%size(i_vertex,i_order) * HHZ   (i_tor,i_plane)
    vv_t  = H_t (i_vertex,i_order,ms,mt) * element%size(i_vertex,i_order) * HHZ   (i_tor,i_plane)
    vv_p  = H   (i_vertex,i_order,ms,mt) * element%size(i_vertex,i_order) * HHZ_p (i_tor,i_plane)
    vv_pp = H   (i_vertex,i_order,ms,mt) * element%size(i_vertex,i_order) * HHZ_pp(i_tor,i_plane)

    vv_ss = H_ss(i_vertex,i_order,ms,mt) * element%size(i_vertex,i_order) * HHZ   (i_tor,i_plane)
    vv_tt = H_tt(i_vertex,i_order,ms,mt) * element%size(i_vertex,i_order) * HHZ   (i_tor,i_plane)
    vv_st = H_st(i_vertex,i_order,ms,mt) * element%size(i_vertex,i_order) * HHZ   (i_tor,i_plane)

    vv_x = (  y_t * vv_s - y_s * vv_t ) / xjac
    vv_y = (- x_t * vv_s + x_s * vv_t ) / xjac

    vv_xx = (vv_ss * y_t**2 - 2.d0*vv_st * y_s*y_t + vv_tt * y_s**2   &		
    	   + vv_s * (y_st*y_t - y_tt*y_s )			      &    
    	   + vv_t * (y_st*y_s - y_ss*y_t ) )    / xjac**2	      &	    
    	   - xjac_x * (vv_s * y_t - vv_t * y_s) / xjac**2

    vv_yy = (vv_ss * x_t**2 - 2.d0*vv_st * x_s*x_t + vv_tt * x_s**2   &		
    	   + vv_s * (x_st*x_t - x_tt*x_s )			      &    
    	   + vv_t * (x_st*x_s - x_ss*x_t ) )      / xjac**2	      &	    
    	   - xjac_y * (- vv_s * x_t + vv_t * x_s ) / xjac**2

    vv_xy = (- vv_ss * y_t*x_t - vv_tt * x_s*y_s		  &
    	     + vv_st * (y_s*x_t  + y_t*x_s  )			  &	   
    	     - vv_s  * (x_st*y_t - x_tt*y_s )			  &	     
    	     - vv_t  * (x_st*y_s - x_ss*y_t )  )       / xjac**2  &		  
    	     - xjac_x * (- vv_s * x_t + vv_t * x_s )   / xjac**2		
    
    return

  end subroutine ELM_build_test_functions












  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------ Equation 1 (psi - induction) ------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!! RHS !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ELM_main_rhs_1(rhs_1)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_rhs_1

    ! --- Modules
    use phys_module
    use equation_variables
    
    implicit none
    
    ! --- Routine variables
    real*8 :: rhs_1
    
    ! --- The RHS term		
    rhs_1 = + v * eta_T  * (zj0 - current_source)/ BigR   * xjac * tstep &
    	    + v * (ps0_x * u0_y - ps0_y * u0_x) 	  * xjac * tstep &
    	    - v * eps_cyl * F0 / BigR  * u0_p		  * xjac * tstep &
    	    + eta_numm * (v_x * zj0_x + v_y * zj0_y)	  * xjac * tstep &
    	    + zeta * v / BigR				  * xjac * delta_g(1)
    
    return

  end subroutine ELM_main_rhs_1

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!! LHS !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ELM_main_lhs_1(amat_11, amat_12, amat_12_n, amat_13, amat_16)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_lhs_1

    ! --- Modules
    use phys_module
    use equation_variables
    
    implicit none
    
    ! --- Routine variables
    real*8 :: amat_11, amat_12, amat_12_n, amat_13, amat_16
    
    ! --- The LHS terms		
    amat_11   = + v * psi / BigR  				   * xjac * (1.d0+zeta)   &
                - v * (psi_x * u0_y - psi_y * u0_x)		   * xjac * theta * tstep

    amat_12   = -  v * (ps0_x * u_y - ps0_y * u_x)		   * xjac * theta * tstep

    amat_12_n = +  eps_cyl * F0 / BigR * v * u_p		   * xjac * theta * tstep

    amat_13   = - eta_numm * (v_x * zj_x + v_y * zj_y)		   * xjac * theta * tstep &
                - eta_T * v * zj / BigR				   * xjac * theta * tstep

    amat_16   = - deta_dT * v * T * (zj0 - current_source) / BigR  * xjac * theta * tstep

    
    return

  end subroutine ELM_main_lhs_1








  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------ Equation 2 (U - momentum) ---------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!! RHS !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ELM_main_rhs_2(rhs_2)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_rhs_2

    ! --- Modules
    use phys_module
    use equation_variables
    
    implicit none
    
    ! --- Routine variables
    real*8 :: rhs_2
    
    ! --- The RHS term		
    rhs_2 = - 0.5d0 * vv2 * (v_x * r0_y_hat - v_y * r0_x_hat)	          * xjac * tstep &
    	    - r0_hat * BigR**2 * w0 * (v_x * u0_y - v_y * u0_x) 	  * xjac * tstep &
    	    + v * (ps0_x * zj0_y - ps0_y * zj0_x )			  * xjac * tstep &
    	    - visco_T * BigR * (v_x * w0_x + v_y * w0_y)		  * xjac * tstep &
    	    - v * eps_cyl * F0 / BigR * zj0_p				  * xjac * tstep &
    	    + BigR**2 * (v_x * p0_y - v_y * p0_x)			  * xjac * tstep &
    	    - zeta * BigR * r0_hat * (v_x * delta_u_x + v_y * delta_u_y)  * xjac	 &		   
    	    - visco_numm  *								 & 
    	      (    (	     v_ss  * (x_t**2+y_t**2)					 & 
    	     	    +	     v_tt  * (x_s**2+y_s**2)					 &
    	     	    - 2.d0 * v_st  * (x_s*x_t + y_s*y_t) )				 &
    	     	 * (	     w0_ss * (x_t**2+y_t**2)					 &
    	     	    +	     w0_tt * (x_s**2+y_s**2)					 &
    	     	    - 2.d0 * w0_st * (x_s*x_t + y_s*y_t) ) )				 &
    	      / xjac**4 						  * xjac * tstep 

    
    return

  end subroutine ELM_main_rhs_2
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!! LHS !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ELM_main_lhs_2(amat_21, amat_22, amat_23, amat_23_n, amat_24, amat_25, amat_26)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_lhs_2

    ! --- Modules
    use phys_module
    use equation_variables
    
    implicit none
    
    ! --- Routine variables
    real*8 :: amat_21, amat_22, amat_23, amat_23_n, amat_24, amat_25, amat_26
    
    ! --- Internal variables
    real*8 :: rho_hat, rho_x_hat, rho_y_hat
    
    rho_hat   = BigR**2 * rho
    rho_x_hat = 2.d0 * BigR * BigR_x  * rho + BigR**2 * rho_x
    rho_y_hat = BigR**2 * rho_y

    ! --- The LHS terms
    amat_21   = - v * (psi_x * zj0_y - psi_y * zj0_x )					* xjac * theta * tstep

    amat_22   = + r0_hat * BigR**2 * w0 * (v_x * u_y  - v_y  * u_x)			* xjac * theta * tstep &
                + BigR**2 * (u_x*u0_x + u_y*u0_y) * (v_x*r0_y_hat - v_y*r0_x_hat)       * xjac * theta * tstep 
    
    if (r0 .lt. rho_1) then
      amat_22 = amat_22 - BigR**3 * rho_1 * (v_x * u_x + v_y * u_y)			* xjac * (1.d0 + zeta) 
    else
      amat_22 = amat_22 - BigR**3 * r0    * (v_x * u_x + v_y * u_y)			* xjac * (1.d0 + zeta) 
    endif

    amat_23   = - v * (ps0_x * zj_y  - ps0_y * zj_x)					* xjac * theta * tstep

    amat_23_n = + eps_cyl * F0 / BigR * v * zj_p					* xjac * theta * tstep

    amat_24   = r0_hat * BigR**2 * w  * ( v_x * u0_y - v_y * u0_x)			* xjac * theta * tstep &
                + BigR * ( v_x * w_x + v_y * w_y) * visco_T 				* xjac * theta * tstep &
                + visco_numm  *										       &
                  (    (	 v_ss * (x_t**2 + y_t**2)						       & 
        	        +	 v_tt * (x_s**2 + y_s**2)						       &
        	        - 2.d0 * v_st * (x_s*x_t + y_s*y_t) )						       &
        	     * (	 w_ss * (x_t**2 + y_t**2)						       &
        	        +	 w_tt * (x_s**2 + y_s**2)						       &
        	        - 2.d0 * w_st * (x_s*x_t + y_s*y_t) ) ) / xjac**4 		* xjac * theta * tstep 

    amat_25   = + 0.5d0 * vv2 * (v_x * rho_y_hat - v_y * rho_x_hat)			* xjac * theta * tstep &
                + rho_hat * BigR**2 * w0 * (v_x * u0_y - v_y * u0_x)			* xjac * theta * tstep &
                - BigR**2 * (v_x * rho_y * T0   - v_y * rho_x * T0  )			* xjac * theta * tstep &
                - BigR**2 * (v_x * rho   * T0_y - v_y * rho   * T0_x )			* xjac * theta * tstep 

    amat_26   = - BigR**2 * (v_x * r0_y * T   - v_y * r0_x * T)				* xjac * theta * tstep &
                - BigR**2 * (v_x * r0   * T_y - v_y * r0   * T_x) 			* xjac * theta * tstep &
                + dvisco_dT * T * ( v_x * w0_x + v_y * w0_y ) * BigR			* xjac * theta * tstep
    
    return

  end subroutine ELM_main_lhs_2




  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------ Equation 3 (j - current) ----------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!! RHS !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ELM_main_rhs_3(rhs_3)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_rhs_3

    ! --- Modules
    use phys_module
    use equation_variables
    
    implicit none
    
    ! --- Routine variables
    real*8 :: rhs_3
    
    ! --- The RHS term		
    !rhs_3 = 0.d0
    rhs_3 = - ( v_x * ps0_x  + v_y * ps0_y + v*zj0) / BigR * xjac * tstep

    
    return

  end subroutine ELM_main_rhs_3

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!! LHS !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ELM_main_lhs_3(amat_31, amat_33)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_lhs_3

    ! --- Modules
    use phys_module
    use equation_variables
    
    implicit none
    
    ! --- Routine variables
    real*8 :: amat_31, amat_33
    
    ! --- The LHS terms
    amat_31 = (v_x * psi_x + v_y * psi_y ) / BigR * xjac * tstep

    amat_33 = v * zj / BigR			  * xjac * tstep 

    
    return

  end subroutine ELM_main_lhs_3







  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------ Equation 4 (w - vorticity) --------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!! RHS !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ELM_main_rhs_4(rhs_4)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_rhs_4

    ! --- Modules
    use phys_module
    use equation_variables
    
    implicit none
    
    ! --- Routine variables
    real*8 :: rhs_4
    
    ! --- The RHS term		
    rhs_4 = 0.d0
    !rhs_4 = - ( v_x * u0_x   + v_y * u0_y  + v*w0)  * BigR * xjac * tstep

    return

  end subroutine ELM_main_rhs_4

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!! LHS !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ELM_main_lhs_4(amat_42, amat_44)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_lhs_4

    ! --- Modules
    use phys_module
    use equation_variables
    
    implicit none
    
    ! --- Routine variables
    real*8 :: amat_42, amat_44
    
    ! --- The LHS terms
    amat_42 = (v_x * u_x + v_y * u_y) * BigR * xjac * tstep

    amat_44 =  v * w * BigR		     * xjac * tstep 
    		      
    
    return

  end subroutine ELM_main_lhs_4







  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------ Equation 5 (rho - continuity) -----------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!! RHS !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ELM_main_rhs_5(rhs_5, rhs_5_k)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_rhs_5

    ! --- Modules
    use phys_module
    use equation_variables
    
    implicit none
    
    ! --- Routine variables
    real*8 :: rhs_5, rhs_5_k
    
    ! --- Parallel gradient terms	
    Bgrad_rho_star   = ( v_x  * ps0_y - v_y  * ps0_x ) / BigR
    Bgrad_rho_k_star = ( F0 / BigR * v_p	     ) / BigR
    		
    ! --- The RHS term		
    rhs_5 =  v * BigR * particle_source				           * xjac * tstep &
    	   + v * BigR**2 * ( r0_x * u0_y - r0_y * u0_x) 		   * xjac * tstep &
    	   + v * 2.d0 * BigR * r0 * u0_y				   * xjac * tstep &
    	   - (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star * Bgrad_rho	   * xjac * tstep &
    	   - D_prof * BigR  * (v_x*r0_x + v_y*r0_y)			   * xjac * tstep &
    	   - v * F0 / BigR * Vpar0 * r0_p				   * xjac * tstep &
    	   - v * Vpar0 * (r0_x * ps0_y - r0_y * ps0_x)  		   * xjac * tstep &
    	   - v * F0 / BigR * r0 * vpar0_p				   * xjac * tstep &
    	   - v * r0 * (vpar0_x * ps0_y - vpar0_y * ps0_x)		   * xjac * tstep &
    	   + zeta * v * BigR						   * xjac *delta_g(5)& 
    	   - D_perp_numm *								  &
    	     (    (	    v_ss  * (x_t**2 + y_t**2)					  &
    	     	   +	    v_tt  * (x_s**2 + y_s**2)					  &
    	     	   - 2.d0 * v_st  * (x_s*x_t + y_s*y_t) )				  &
    	     	* (	    r0_ss * (x_t**2 + y_t**2)					  &
    	     	   +	    r0_tt * (x_s**2 + y_s**2)					  &
    	     	   - 2.d0 * r0_st * (x_s*x_t + y_s*y_t) ) )				  &
    	     / xjac**3 *tstep								  

    rhs_5_k = - (D_par-D_prof) * BigR / BB2 * Bgrad_rho_k_star * Bgrad_rho * xjac * tstep &
    	      - D_prof * BigR  * ( v_p*r0_p * eps_cyl**2 /BigR**2 )	   * xjac * tstep 
    
    return

  end subroutine ELM_main_rhs_5

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!! LHS !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ELM_main_lhs_5(amat_51, amat_51_k, amat_52, amat_55, amat_55_k, amat_55_n, amat_55_kn, amat_57, amat_57_n)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_lhs_5

    ! --- Modules
    use phys_module
    use equation_variables
    
    implicit none
    
    ! --- Routine variables
    real*8 :: amat_51, amat_51_k, amat_52, amat_55, amat_55_k, amat_55_n, amat_55_kn, amat_57, amat_57_n
    
    ! --- Internal variables
    real*8 :: Bgrad_rho_star_psi, Bgrad_rho_psi, Bgrad_rho_rho, Bgrad_rho_rho_n
    
    ! --- Parallel gradient terms	
    Bgrad_rho_star_psi = ( v_x   * psi_y - v_y   * psi_x ) / BigR
    Bgrad_rho_psi      = ( r0_x  * psi_y - r0_y  * psi_x ) / BigR
    Bgrad_rho_rho      = ( rho_x * ps0_y - rho_y * ps0_x ) / BigR
    Bgrad_rho_rho_n    = ( F0 / BigR * rho_p ) / BigR

    ! --- The LHS terms
    amat_51    = - (D_par-D_prof) * BigR * BB2_psi/ BB2**2 * Bgrad_rho_star     * Bgrad_rho     * xjac * theta * tstep &
                 + (D_par-D_prof) * BigR / BB2		   * Bgrad_rho_star_psi * Bgrad_rho     * xjac * theta * tstep &
                 + (D_par-D_prof) * BigR / BB2		   * Bgrad_rho_star     * Bgrad_rho_psi * xjac * theta * tstep &
                 + v * Vpar0 * (r0_x * psi_y - r0_y * psi_x)				        * xjac * theta * tstep &
                 + v * r0 * (vpar0_x * psi_y - vpar0_y * psi_x)				        * xjac * theta * tstep 

    amat_51_k  = - (D_par-D_prof) * BigR * BB2_psi/ BB2**2 * Bgrad_rho_k_star * Bgrad_rho       * xjac * theta * tstep &
        	 + (D_par-D_prof) * BigR / BB2		   * Bgrad_rho_k_star * Bgrad_rho_psi   * xjac * theta * tstep 

    amat_52    = - v * BigR**2 * ( r0_x * u_y - r0_y * u_x)				        * xjac * theta * tstep &
                 - v * 2.d0 * BigR * r0 * u_y						        * xjac * theta * tstep 

    amat_55    = + v * rho * BigR							        * xjac * (1.d0 + zeta) &
                 - v * BigR**2 * ( rho_x * u0_y - rho_y * u0_x)				        * xjac * theta * tstep &
                 - v * 2.d0 * BigR * rho * u0_y						        * xjac * theta * tstep &
                 + (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star * Bgrad_rho_rho		        * xjac * theta * tstep &
                 + D_prof * BigR  * (v_x*rho_x + v_y*rho_y )				        * xjac * theta * tstep &
                 + v * Vpar0 * (rho_x * ps0_y - rho_y * ps0_x)				        * xjac * theta * tstep &
                 + v * rho * (vpar0_x * ps0_y - vpar0_y * ps0_x)			        * xjac * theta * tstep &
                 + v * rho * F0 / BigR * vpar0_p					        * xjac * theta * tstep &
                 + D_perp_num  *							        		       &
        	   (    (         v_ss   * (x_t**2 + y_t**2)				        		       &
        	         +        v_tt   * (x_s**2 + y_s**2)				        		       &
        	         - 2.d0 * v_st   * (x_s*x_t + y_s*y_t) )			        		       &
        	      * (         rho_ss * (x_t**2 + y_t**2)				        		       &
        	         +        rho_tt * (x_s**2 + y_s**2)				        		       &
        	         - 2.d0 * rho_st * (x_s*x_t + y_s*y_t) ) ) / xjac**4		        * xjac * theta * tstep 

    amat_55_k  = + (D_par-D_prof) * BigR / BB2 * Bgrad_rho_k_star * Bgrad_rho_rho	        * xjac * theta * tstep 

    amat_55_n  = + (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star   * Bgrad_rho_rho_n	        * xjac * theta * tstep &
        	 + v * F0 / BigR * Vpar0 * rho_p 					        * xjac * theta * tstep 

    amat_55_kn = + (D_par-D_prof) * BigR / BB2 * Bgrad_rho_k_star * Bgrad_rho_rho_n	        * xjac * theta * tstep &
        	 + D_prof * BigR  * ( v_p*rho_p * eps_cyl**2 /BigR**2 ) 		        * xjac * theta * tstep 

    amat_57    = + v * F0 / BigR * Vpar * r0_p						        * xjac * theta * tstep &
        	 + v * Vpar * (r0_x * ps0_y - r0_y * ps0_x)				        * xjac * theta * tstep &
        	 + v * r0 * (vpar_x * ps0_y - vpar_y * ps0_x)				        * xjac * theta * tstep 

    amat_57_n  = + v * r0 * F0 / BigR * vpar_p						        * xjac * theta * tstep
    		      
    
    return

  end subroutine ELM_main_lhs_5


  
  
  
  
  
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------ Equation 6 (Ti - Ion energy) ------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!! RHS !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ELM_main_rhs_6(rhs_6, rhs_6_k)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_rhs_6

    ! --- Modules
    use phys_module
    use equation_variables
    
    implicit none
    
    ! --- Routine variables
    real*8 :: rhs_6, rhs_6_k
    
    ! --- Parallel gradient terms	
    v_ps0_x  = v_xx  * ps0_y - v_xy  * ps0_x + v_x  * ps0_xy - v_y * ps0_xx
    v_ps0_y  = v_xy  * ps0_y - v_yy  * ps0_x + v_x  * ps0_yy - v_y * ps0_xy    
    Bgrad_T_star     = ( v_x  * ps0_y - v_y  * ps0_x ) / BigR
    Bgrad_T_k_star   = ( F0 / BigR * v_p	     ) / BigR
    		
    ! --- The RHS term		
    rhs_6 =   v * BigR * heat_source						* xjac * tstep &
            + v * r0 * BigR**2 * ( T0_x * u0_y - T0_y * u0_x)			* xjac * tstep &
            + v * r0 * 2.d0* (GAMMA-1.d0) * BigR * T0 * u0_y		        * xjac * tstep &
            - v * r0 * F0 / BigR * Vpar0 * T0_p 				* xjac * tstep &
            - v * r0 * Vpar0 * (T0_x * ps0_y - T0_y * ps0_x)			* xjac * tstep &
            - v * r0 * (GAMMA-1.d0) * T0 * (vpar0_x * ps0_y - vpar0_y * ps0_x)  * xjac * tstep &
            - v * r0 * (GAMMA-1.d0) * T0 * F0 / BigR * vpar0_p  		* xjac * tstep &
            -  (K_par-ZK_prof) * BigR / BB2 * Bgrad_T_star * Bgrad_T	        * xjac * tstep &
            -  ZK_prof * BigR * (v_x*T0_x + v_y*T0_y )  			* xjac * tstep &
            -  K_par_num * (v_ps0_x  * ps0_y - v_ps0_y  * ps0_x) &
    	     		 * (T0_ps0_x * ps0_y - T0_ps0_y * ps0_x)		* xjac * tstep &
            + zeta * v * r0 * BigR						* xjac *delta_g(6)&
            - K_perp_numm2  *								       &
    	      (    (	     v_ss  * (x_t**2 + y_t**2)  				       &
    	     	    +	     v_tt  * (x_s**2 + y_s**2)  				       &
    	     	    - 2.d0 * v_st  * (x_s*x_t + y_s*y_t) )				       &
    	     	 * (	     T3_ss * (x_t**2 + y_t**2)  				       &
    	     	    +	     T3_tt * (x_s**2 + y_s**2)  				       &
    	     	    - 2.d0 * T3_st * (x_s*x_t + y_s*y_t) ) )				       &
    	      / xjac**4 							* xjac * tstep 

    rhs_6_k = - (K_par-ZK_prof) * BigR / BB2 * Bgrad_T_k_star * Bgrad_T 	* xjac * tstep &
    	      -  K_par_num2 * (T0_pp  * v_pp)					* xjac * tstep &
    	      - ZK_prof * BigR * ( v_p*T0_p /BigR**2 )  			* xjac * tstep 

    
    return

  end subroutine ELM_main_rhs_6


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!! LHS !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ELM_main_lhs_6(amat_61, amat_61_k, amat_62, amat_65, amat_66, amat_66_k, amat_66_n, amat_66_kn, amat_67, amat_67_n)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_lhs_6

    ! --- Modules
    use phys_module
    use equation_variables
    
    implicit none
    
    ! --- Routine variables
    real*8 :: amat_61, amat_61_k, amat_62, amat_65, amat_66, amat_66_k, amat_66_n, amat_66_kn, amat_67, amat_67_n
    
    ! --- Internal variables
    real*8 :: Bgrad_T_star_psi, Bgrad_T_psi, Bgrad_T_T, Bgrad_T_T_n
    real*8 :: T_ps0_x,  T_ps0_y
    real*8 :: T0_psi_x, T0_psi_y
    real*8 :: v_psi_x,  v_psi_y
    
    ! --- Parallel gradient terms	
    Bgrad_T_star_psi   = ( v_x   * psi_y - v_y   * psi_x ) / BigR
    Bgrad_T_psi        = ( T0_x  * psi_y - T0_y  * psi_x ) / BigR
    Bgrad_T_T	       = ( T_x   * ps0_y - T_y   * ps0_x ) / BigR
    Bgrad_T_T_n        = ( F0 / BigR * T_p   ) / BigR

    T_ps0_x = T_xx * ps0_y - T_xy * ps0_x + T_x * ps0_xy - T_y * ps0_xx
    T_ps0_y = T_xy * ps0_y - T_yy * ps0_x + T_x * ps0_yy - T_y * ps0_xy
    
    T0_psi_x = T0_xx * psi_y - T0_xy * psi_x + T0_x * psi_xy - T0_y * psi_xx
    T0_psi_y = T0_xy * psi_y - T0_yy * psi_x + T0_x * psi_yy - T0_y * psi_xy
    
    v_psi_x = v_xx * psi_y - v_xy * psi_x + v_x * psi_xy - v_y * psi_xx
    v_psi_y = v_xy * psi_y - v_yy * psi_x + v_x * psi_yy - v_y * psi_xy

    ! --- The LHS terms
    amat_61    = + v * r0 * Vpar0 * (T0_x * psi_y - T0_y * psi_x)  				  * xjac * theta * tstep &
                 + v * (GAMMA-1.d0) * r0 * T0 * (vpar0_x * psi_y - vpar0_y * psi_x)		  * xjac * theta * tstep &
                 + K_par_num * (v_psi_x *ps0_y - v_psi_y *ps0_x + v_ps0_x *psi_y - v_ps0_y *psi_x)			 &
              	   	     * (T0_ps0_x*ps0_y - T0_ps0_y*ps0_x)				  * xjac * theta * tstep &
                 + K_par_num * (T0_psi_x*ps0_y - T0_psi_y*ps0_x + T0_ps0_x*psi_y - T0_ps0_y*psi_x)			 &
              	   	     * (v_ps0_x *ps0_y - v_ps0_y *ps0_x)				  * xjac * theta * tstep &
                 - (K_par-ZK_prof) * BigR * BB2_psi / BB2**2 * Bgrad_T_star	* Bgrad_T	  * xjac * theta * tstep &
                 + (K_par-ZK_prof) * BigR / BB2 	     * Bgrad_T_star_psi * Bgrad_T	  * xjac * theta * tstep &
                 + (K_par-ZK_prof) * BigR / BB2 	     * Bgrad_T_star	* Bgrad_T_psi	  * xjac * theta * tstep 

    amat_61_k  = - (K_par-ZK_prof) * BigR * BB2_psi / BB2**2 * Bgrad_T_k_star   * Bgrad_T  	  * xjac * theta * tstep &
        	 + (K_par-ZK_prof) * BigR / BB2  	     * Bgrad_T_k_star   * Bgrad_T_psi	  * xjac * theta * tstep 

    amat_62    = - v * r0 * BigR**2 * ( T0_x * u_y - T0_y * u_x)				  * xjac * theta * tstep &
                 - v * 2.d0 * (GAMMA-1.d0) * r0 * BigR * T0 * u_y  				  * xjac * theta * tstep 

    amat_65    = - v * rho * BigR**2 * (T0_x * u0_y - T0_y * u0_x) 				  * xjac * theta * tstep &
                 + v * rho * Vpar0 * F0/BigR * T0_p						  * xjac * theta * tstep &
                 + v * rho * Vpar0 * (T0_x * ps0_y - T0_y * ps0_x) 				  * xjac * theta * tstep &
                 - v * 2.d0 * (GAMMA-1.d0) * rho * BigR * T0 * u0_y				  * xjac * theta * tstep &
                 + v * (GAMMA-1.d0) * rho * T0 * F0/BigR * Vpar0_p 				  * xjac * theta * tstep &
                 + v * (GAMMA-1.d0) * rho * T0 * (Vpar0_x * ps0_y - Vpar0_y * ps0_x)		  * xjac * theta * tstep 

    amat_66    = - v * r0 * BigR**2 * (T_x * u0_y  - T_y * u0_x)				  * xjac * theta * tstep &
                 + v * r0 * Vpar0   * (T_x * ps0_y - T_y * ps0_x) 				  * xjac * theta * tstep &
                 - 2.d0 * v * r0 * (GAMMA-1.d0) * T * BigR * u0_y				  * xjac * theta * tstep &
                 + v * r0 * (GAMMA-1.d0) * T * F0/BigR * Vpar0_p				  * xjac * theta * tstep &
                 + v * r0 * (GAMMA-1.d0) * T * (Vpar0_x * ps0_y - Vpar0_y * ps0_x)		  * xjac * theta * tstep &
                 + K_par_num * (v_ps0_x * ps0_y - v_ps0_y * ps0_x) &
              	   	     * (T_ps0_x * ps0_y - T_ps0_y * ps0_x)				  * xjac * theta * tstep &
                 + (K_par-ZK_prof) * BigR / BB2 * Bgrad_T_star * Bgrad_T_T			  * xjac * theta * tstep &
                 + dZK_par * T     * BigR / BB2 * Bgrad_T_star * Bgrad_T			  * xjac * theta * tstep &
                 + ZK_prof * BigR * (v_x*T_x + v_y*T_y )					  * xjac * theta * tstep & 
                 + K_perp_numm  *								     			 &
              	   (	(	  v_ss * (x_t**2 + y_t**2)					     			 &
              	   	 +	  v_tt * (x_s**2 + y_s**2)					     			 &
              	   	 - 2.d0 * v_st * (x_s*x_t + y_s*y_t) )  				     			 &
              	      * (	  T_ss * (x_t**2 + y_t**2)					     			 &
              	   	 +	  T_tt * (x_s**2 + y_s**2)					     			 &
              	   	 - 2.d0 * T_st * (x_s*x_t + y_s*y_t) ) )     / xjac**4  		  * xjac * theta * tstep 
    
    if (r0 .lt. rho_1) then
      amat_66  = amat_66 + v * rho_1 * T * BigR  						  * xjac * (1.d0 + zeta)
    else
      amat_66  = amat_66 + v * r0    * T * BigR  						  * xjac * (1.d0 + zeta)
    endif
    
    !if ((im .ne. 1) .and. (in .ne. 1)) then
    !  amat_66  = amat_66 + K_perp_numm2  *										 &
    !    		  (    (       v_ss * (x_t**2 + y_t**2) 							 &
    !    		      +        v_tt * (x_s**2 + y_s**2) 							 &
    !    		      - 2.d0 * v_st * (x_s*x_t + y_s*y_t) )							 &
    !    		     * (       T_ss * (x_t**2 + y_t**2) 							 &
    !    		      +        T_tt * (x_s**2 + y_s**2) 							 &
    !    		      - 2.d0 * T_st * (x_s*x_t + y_s*y_t) ) )	  / xjac**4		  * xjac * theta * tstep 
    !endif

    amat_66_k  = + (K_par-ZK_prof) * BigR / BB2 * Bgrad_T_k_star * Bgrad_T_T			  * xjac * theta * tstep &
        	 + dZK_par * T	   * BigR / BB2 * Bgrad_T_k_star * Bgrad_T			  * xjac * theta * tstep 
        	
    amat_66_n  = + (K_par-ZK_prof) * BigR / BB2 * Bgrad_T_star   * Bgrad_T_T_n  		  * xjac * theta * tstep &
        	 + v * r0 * Vpar0  * F0/BigR * T_p						  * xjac * theta * tstep 

    amat_66_kn = + (K_par-ZK_prof) * BigR / BB2 * Bgrad_T_k_star * Bgrad_T_T_n  		  * xjac * theta * tstep &
        	 +  K_par_num2 * (T_pp  * v_pp) 						  * xjac * theta * tstep &
        	 + ZK_prof	   * BigR	* (v_p*T_p /BigR**2 )				  * xjac * theta * tstep 

    amat_67    = + v * r0 * F0/BigR * Vpar * T0_p  						  * xjac * theta * tstep &
                 + v * r0 * Vpar * (T0_x * ps0_y - T0_y * ps0_x)				  * xjac * theta * tstep &
                 + v * (GAMMA-1.d0) * r0 * T0 * (Vpar_x * ps0_y - Vpar_y * ps0_x)  		  * xjac * theta * tstep 
       
    amat_67_n  = v * (GAMMA-1.d0) * r0 * T0 * F0/BigR * Vpar_p					  * xjac * theta * tstep       
    		     

    
    return

  end subroutine ELM_main_lhs_6


  
  
  
  
  
  
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------ Equation 7 (Vpar - parallel momentum) ---------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!! RHS !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ELM_main_rhs_7(rhs_7)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_rhs_7

    ! --- Modules
    use phys_module
    use equation_variables
    
    implicit none
    
    ! --- Routine variables
    real*8 :: rhs_7
    
    ! --- The RHS term		
    rhs_7 = - v * F0 / BigR * P0_p				   * xjac * tstep &
            - v * (P0_x * ps0_y - P0_y * ps0_x) 		   * xjac * tstep &
            - visco_par * (v_x * vpar0_x + v_y * vpar0_y) * BigR   * xjac * tstep &
            + zeta * v * r0 * F0**2 / BigR			   * xjac *delta_g(7)&  
            - visco_par_numm  * 						  &
    	      (    (	     v_ss     * (x_t**2 + y_t**2)			  &
    	     	    +	     v_tt     * (x_s**2 + y_s**2)			  &
    	     	    - 2.d0 * v_st     * (x_s*x_t + y_s*y_t) )			  &
    	     	 * (	     vpar0_ss * (x_t**2 + y_t**2)			  &
    	     	    +	     vpar0_tt * (x_s**2 + y_s**2)			  &
    	     	    - 2.d0 * vpar0_st * (x_s*x_t + y_s*y_t) ) ) 		  &
    	      / xjac**4 					   * xjac * tstep 
    
    return

  end subroutine ELM_main_rhs_7

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!! LHS !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ELM_main_lhs_7(amat_71, amat_72, amat_75, amat_75_n, amat_76, amat_76_n, amat_77)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_lhs_7

    ! --- Modules
    use phys_module
    use equation_variables
    
    implicit none
    
    ! --- Routine variables
    real*8 :: amat_71, amat_72, amat_75, amat_75_n, amat_76, amat_76_n, amat_77
    
    ! --- The LHS terms
    amat_71   = + v * (P0_x * psi_y - P0_y * psi_x)				     * xjac * theta * tstep !&
                !+ vpar0 * (F0/BigR)**2 * (vpar0_x * ps_y - vpar0_y * ps_x)	     * xjac * theta * tstep &

    amat_72   = 0.d0!+ F0 * (u_s * vpar0_t - u_t * vpar0_s) * theta * tstep 

    amat_75   = + v * (rho_x * T0   * ps0_y - rho_y * T0   * ps0_x)		     * xjac * theta * tstep &
                + v * (rho   * T0_x * ps0_y - rho   * T0_y * ps0_x)		     * xjac * theta * tstep &
                + v * F0 / BigR * rho * T0_p					     * xjac * theta * tstep 

    amat_75_n = v * F0 / BigR * rho_p * T0					     * xjac * theta * tstep  

    amat_76   = + v * (T_x * r0	  * ps0_y - T_y * r0   * ps0_x)			     * xjac * theta * tstep &
                + v * (T   * r0_x * ps0_y - T   * r0_y * ps0_x)			     * xjac * theta * tstep &
                + v * F0 / BigR * T * r0_p					     * xjac * theta * tstep 
    
    amat_76_n = v * F0 / BigR * T_p * r0					     * xjac * theta * tstep

    amat_77   = + visco_par * (v_x * Vpar_x + v_y * Vpar_y) * BigR		     * xjac * theta * tstep &
                + visco_par_numm  *						        		    &
              	  (    (	v_ss	* (x_t**2 + y_t**2)			        		    &
              	  	+	v_tt	* (x_s**2 + y_s**2)			        		    &
              	  	-2.d0 * v_st	* (x_s*x_t + y_s*y_t) ) 		        		    &
              	     * (	Vpar_ss * (x_t**2 + y_t**2)			        		    &
              	  	+	Vpar_tt * (x_s**2 + y_s**2)			        		    &
              	  	-2.d0 * Vpar_st * (x_s*x_t + y_s*y_t) ) ) / xjac**4	     * xjac * theta * tstep !&
                !+ F0 * (u0_x * vpar_y - u0_y * vpar_x) 			     * xjac * theta * tstep &
                !+ vpar * (F0/BigR)**2  * (vpar0_x * ps0_y - vpar0_y * ps0_x)	     * xjac * theta * tstep &
                !+ vpar0 * (F0/BigR)**2 * (vpar_x  * ps0_y - vpar_y * ps0_x)	     * xjac * theta * tstep &
                !+ (F0/BigR)**3 * vpar  * vpar0_p				     * xjac * theta * tstep &
                !+ (F0/BigR)**3 * vpar0 * vpar_p				     * xjac * theta * tstep &
    
    if (r0 .lt. rho_1) then
      amat_77 = amat_77 + v * Vpar * rho_1 * F0**2 / BigR		             * xjac * (1.d0 + zeta)
    else
      amat_77 = amat_77 + v * Vpar * r0    * F0**2 / BigR		             * xjac * (1.d0 + zeta)
    endif

    
    return

  end subroutine ELM_main_lhs_7

















  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !-------------------------------------- Compute the FFT for the matrix elements -----------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  subroutine ELM_apply_fft(RHS, RHS_p, RHS_k, ELM, ELM_p, ELM_n, ELM_k, ELM_kn, tid)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_apply_fft

    ! --- Modules
    use phys_module
    use parameters    
    use data_structure
    
    implicit none
    
    ! --- Matrix elements and toroidal functions
    integer, intent(in)     	:: tid
    real*8, dimension (:,:) 	:: ELM
    real*8, dimension (:)   	:: RHS
    real*8, dimension(:,:,:)	:: ELM_p, ELM_n, ELM_k, ELM_kn
    real*8, dimension(:,:)  	:: RHS_p, RHS_k 
    real*8		    	:: in_fft(1:n_plane)
    complex*16  	    	:: out_fft(1:n_plane)
    
    ! --- FFT Indexes
    integer    :: i, j, k, l, m, ik, im
    integer    :: index, index_k, index_m
            
    ! --- RHS_p
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

    ! --- RHS_k
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

    ! --- ELM_p
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
      enddo
    enddo

    ! --- ELM_n
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

    ! --- ELM_k
    do i=1,n_vertex_max*n_var*(n_order+1)
      do j=1, n_vertex_max*n_var*(n_order+1)
        if (maxval(abs(ELM_k(1:n_plane,i,j))) .ne. 0.d0) then

          in_fft =  ELM_k(1:n_plane,i,j)

          call my_fft(in_fft, out_fft, n_plane)

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
        endif
      enddo
    enddo


    ! --- ELM_kn
    do i=1,n_vertex_max*n_var*(n_order+1)
      do j=1, n_vertex_max*n_var*(n_order+1)
        if (maxval(abs(ELM_kn(1:n_plane,i,j))) .ne. 0.d0) then

          in_fft =  ELM_kn(1:n_plane,i,j)

          call my_fft(in_fft, out_fft, n_plane)

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
        endif
      enddo
    enddo

    ELM = 0.5d0 * ELM

    return
  end subroutine ELM_apply_fft
















  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------ The FFT routine -------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  subroutine my_fft(in_fft,out_fft,n)
!DEC$ ATTRIBUTES FORCEINLINE :: my_fft
      
    real*8     :: in_fft(*)
    complex*16 :: out_fft(*)	  
    real*8     :: tmp_fft(2*n+2)
    integer    :: i, n
    	
    tmp_fft(1:n) = in_fft(1:n)      
    call RFT2(tmp_fft,n,1)
    	
    do i=1,n
      out_fft(i) = cmplx(tmp_fft(2*i-1),tmp_fft(2*i))
    enddo
    	
    return
  
  end subroutine my_fft





end module mod_elt_matrix_fft





! This module contains nothing but it is needed by construct_matrix.
! Can be removed once the other models have also combined element_matrix and element_matrix_fft.
module mod_elt_matrix
contains

  subroutine element_matrix(element, nodes, xpoint2, xcase2, minRad, R_axis, Z_axis, psi_axis, psi_bnd, Z_xpoint, ELM, RHS, tid)
  !--------------------------------------------------------------------------
  ! This is just a wrapper to the real routine since I combined both into one
  !--------------------------------------------------------------------------

    use data_structure
    use mod_elt_matrix_fft

    implicit none

    type (type_element) 	      :: element
    type (type_node)		      :: nodes(n_vertex_max)

    integer    :: xcase2
    logical    :: xpoint2
    real*8     :: minRad, R_axis, Z_axis, psi_axis, psi_bnd, Z_xpoint(2)
    real*8, dimension (:,:), pointer  :: ELM
    real*8, dimension (:)  , pointer  :: RHS
    integer, intent(in) 	      :: tid

    call element_matrix_fft(element, nodes, xpoint2, xcase2, minRad, R_axis, Z_axis, psi_axis, psi_bnd, Z_xpoint, ELM, RHS, tid)

    return

  end subroutine element_matrix

end module mod_elt_matrix

