module mod_elt_matrix_fft
contains

  ! --- Include all the routines directly for runtime efficiency
  INCLUDE "../model333/construct_variables.f90"
  INCLUDE "../model333/equations.f90"
  INCLUDE "../model333/equations_numm.f90"
  INCLUDE "construct_neutral_variables.f90"
  INCLUDE "equations_neutral_terms.f90"


  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------ Calculates the matrix contribution of one element ---------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  subroutine element_matrix_fft(element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, ELM, RHS, tid)

    ! --- Modules
    use equation_variables
    use constants
    use mod_parameters
    use data_structure
    use gauss
    use basis_at_gaussian
    use phys_module
    use tr_module 
    use profiles, only: interpolProf
    use diffusivities, only: get_dperp, get_zkperp
    use corr_neg
    use pellet_module
    use mod_elm_apply_fft

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
    real*8     :: R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint(2), Z_xpoint(2)
    
    ! --- Integration weight
    real*8     :: wst
        
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

    	  call ELM_build_neutral_variables(element, nodes, ms, mt, i_plane)
          
    	  call ELM_build_diffusivities_and_sources        (element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, i_plane)
          call ELM_build_neutral_diffusivities_and_sources(element, nodes, xpoint2, xcase2,         R_axis, Z_axis, psi_axis, psi_bnd,           Z_xpoint, i_plane)

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
		
	        ! --- Build test functions (which we choose to be the basis functions)
		call ELM_build_basis_functions(element, nodes, ms, mt, i_plane, i_vertex, i_order, i_tor, &
					       v, v_s,  v_t,	    v_p,  v_x,  v_y,			  &
						  v_ss, v_tt, v_st, v_pp, v_xx, v_yy, v_xy		  )
		
		rhs_tmp   = 0.d0
		rhs_k_tmp = 0.d0
		call ELM_main_rhs_1  	(rhs_tmp, rhs_k_tmp)
		call ELM_main_rhs_2  	(rhs_tmp, rhs_k_tmp)
		call ELM_main_rhs_3  	(rhs_tmp, rhs_k_tmp)
		call ELM_main_rhs_4  	(rhs_tmp, rhs_k_tmp)
		call ELM_main_rhs_5  	(rhs_tmp, rhs_k_tmp)
		call ELM_main_rhs_6  	(rhs_tmp, rhs_k_tmp)
		call ELM_main_rhs_7  	(rhs_tmp, rhs_k_tmp)
		call ELM_neutral_rhs_2  (rhs_tmp, rhs_k_tmp)
		call ELM_neutral_rhs_5  (rhs_tmp, rhs_k_tmp)
		call ELM_neutral_rhs_6  (rhs_tmp, rhs_k_tmp)
		call ELM_neutral_rhs_7  (rhs_tmp, rhs_k_tmp)
		call ELM_neutral_rhs_8  (rhs_tmp, rhs_k_tmp)
		call ELM_main_rhs_2_numm(rhs_tmp, rhs_k_tmp)
		call ELM_main_rhs_5_numm(rhs_tmp, rhs_k_tmp)
		call ELM_main_rhs_6_numm(rhs_tmp, rhs_k_tmp)
		call ELM_main_rhs_7_numm(rhs_tmp, rhs_k_tmp)
    		

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

		      ! --- Build basis functions
		      call ELM_build_basis_functions(element, nodes, ms, mt, i_plane, j_vertex, j_order, j_tor,  &
                                                     psi, psi_s,  psi_t,          psi_p,  psi_x,  psi_y,         &
				                          psi_ss, psi_tt, psi_st, psi_pp, psi_xx, psi_yy, psi_xy )
    		      
		      u    = psi; u_x    = psi_x; u_y    = psi_y; u_p    = psi_p; u_s    = psi_s; u_t    = psi_t
		      zj   = psi; zj_x   = psi_x; zj_y   = psi_y; zj_p   = psi_p; zj_s   = psi_s; zj_t   = psi_t
		      w    = psi; w_x    = psi_x; w_y    = psi_y; w_p    = psi_p; w_s    = psi_s; w_t    = psi_t
		      rho  = psi; rho_x  = psi_x; rho_y  = psi_y; rho_p  = psi_p; rho_s  = psi_s; rho_t  = psi_t
		      T    = psi; T_x    = psi_x; T_y    = psi_y; T_p    = psi_p; T_s    = psi_s; T_t    = psi_t
		      Vpar = psi; Vpar_x = psi_x; Vpar_y = psi_y; Vpar_p = psi_p; Vpar_s = psi_s; Vpar_t = psi_t
		      rhon = psi; rhon_x = psi_x; rhon_y = psi_y; rhon_p = psi_p; rhon_s = psi_s; rhon_t = psi_t
		      
		      u_ss    = psi_ss; u_tt	= psi_tt; u_st    = psi_st
		      zj_ss   = psi_ss; zj_tt	= psi_tt; zj_st   = psi_st
		      w_ss    = psi_ss; w_tt	= psi_tt; w_st    = psi_st
   		      rho_ss  = psi_ss; rho_tt  = psi_tt; rho_st  = psi_st
   		      T_ss    = psi_ss; T_tt	= psi_tt; T_st    = psi_st
   		      Vpar_ss = psi_ss; Vpar_tt = psi_tt; Vpar_st = psi_st
   		      rhon_ss = psi_ss; rhon_tt = psi_tt; rhon_st = psi_st
                      
		      u_xx    = psi_xx; u_yy    = psi_yy; u_xy    = psi_xy; u_pp    = psi_pp
		      zj_xx   = psi_xx; zj_yy   = psi_yy; zj_xy   = psi_xy; zj_pp   = psi_pp
		      w_xx    = psi_xx; w_yy    = psi_yy; w_xy    = psi_xy; w_pp    = psi_pp
		      rho_xx  = psi_xx; rho_yy  = psi_yy; rho_xy  = psi_xy; rho_pp  = psi_pp
		      T_xx    = psi_xx; T_yy    = psi_yy; T_xy    = psi_xy; T_pp    = psi_pp
		      Vpar_xx = psi_xx; Vpar_yy = psi_yy; Vpar_xy = psi_xy; Vpar_pp = psi_pp
		      rhon_xx = psi_xx; rhon_yy = psi_yy; rhon_xy = psi_xy; rhon_pp = psi_pp
                      
		      BB2_psi = 2.d0 * (psi_x * ps0_x + psi_y * ps0_y ) /R**2
   
    		      amat_tmp    = 0.d0
    		      amat_k_tmp  = 0.d0
    		      amat_n_tmp  = 0.d0
    		      amat_kn_tmp = 0.d0
		      call ELM_main_lhs_1     (amat_tmp, amat_k_tmp, amat_n_tmp, amat_kn_tmp)
    		      call ELM_main_lhs_2     (amat_tmp, amat_k_tmp, amat_n_tmp, amat_kn_tmp)
    		      call ELM_main_lhs_3     (amat_tmp, amat_k_tmp, amat_n_tmp, amat_kn_tmp)
    		      call ELM_main_lhs_4     (amat_tmp, amat_k_tmp, amat_n_tmp, amat_kn_tmp)
    		      call ELM_main_lhs_5     (amat_tmp, amat_k_tmp, amat_n_tmp, amat_kn_tmp)
    		      call ELM_main_lhs_6     (amat_tmp, amat_k_tmp, amat_n_tmp, amat_kn_tmp)
    		      call ELM_main_lhs_7     (amat_tmp, amat_k_tmp, amat_n_tmp, amat_kn_tmp)
    		      call ELM_neutral_lhs_2  (amat_tmp, amat_k_tmp, amat_n_tmp, amat_kn_tmp)
    		      call ELM_neutral_lhs_5  (amat_tmp, amat_k_tmp, amat_n_tmp, amat_kn_tmp)
    		      call ELM_neutral_lhs_6  (amat_tmp, amat_k_tmp, amat_n_tmp, amat_kn_tmp)
    		      call ELM_neutral_lhs_7  (amat_tmp, amat_k_tmp, amat_n_tmp, amat_kn_tmp)
    		      call ELM_neutral_lhs_8  (amat_tmp, amat_k_tmp, amat_n_tmp, amat_kn_tmp)
    		      call ELM_main_lhs_2_numm(amat_tmp, amat_k_tmp, amat_n_tmp, amat_kn_tmp)
    		      call ELM_main_lhs_5_numm(amat_tmp, amat_k_tmp, amat_n_tmp, amat_kn_tmp)
    		      call ELM_main_lhs_6_numm(amat_tmp, amat_k_tmp, amat_n_tmp, amat_kn_tmp)
    		      call ELM_main_lhs_7_numm(amat_tmp, amat_k_tmp, amat_n_tmp, amat_kn_tmp)
		      
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



end module mod_elt_matrix_fft





! This module contains nothing (just a wrapper) but it is needed by construct_matrix for the other models.
! Can be removed once the other models have also combined element_matrix and element_matrix_fft.
module mod_elt_matrix
contains

  subroutine element_matrix(element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, ELM, RHS, tid)
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
    real*8     :: R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint(2), Z_xpoint(2)
    real*8, dimension (:,:), pointer  :: ELM
    real*8, dimension (:)  , pointer  :: RHS
    integer, intent(in) 	      :: tid

    call element_matrix_fft(element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, ELM, RHS, tid)

    return

  end subroutine element_matrix

end module mod_elt_matrix

