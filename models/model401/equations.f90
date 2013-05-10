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
subroutine ELM_main_rhs_1(rhs,rhs_k)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_rhs_1

  ! --- Modules
  use parameters
  use phys_module
  use equation_variables
  
  implicit none
  
  ! --- Routine variables
  real*8 :: rhs(n_var),rhs_k(n_var)
  
  ! --- The RHS term	      
  rhs(1) = + v * eta_T  * (zj0 - current_source)/ BigR  * xjac * tstep &
	   + v * (ps0_x * u0_y - ps0_y * u0_x)  	* xjac * tstep &
	   - v * eps_cyl * F0 / BigR  * u0_p		* xjac * tstep &
	   + eta_numm * (v_x * zj0_x + v_y * zj0_y)	* xjac * tstep &
	   + zeta * v / BigR				* xjac * delta_g(1)
  
  return

end subroutine ELM_main_rhs_1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! LHS !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ELM_main_lhs_1(amat, amat_k, amat_n, amat_kn)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_lhs_1

  ! --- Modules
  use phys_module
  use equation_variables
  
  implicit none
  
  ! --- Routine variables
  real*8 :: amat(n_var,n_var), amat_k(n_var,n_var), amat_n(n_var,n_var), amat_kn(n_var,n_var)
  
  ! --- The LHS terms	      
  amat(1,1)   = + v * psi / BigR				   * xjac * (1.d0+zeta)   &
		- v * (psi_x * u0_y - psi_y * u0_x)		   * xjac * theta * tstep

  amat(1,2)   = -  v * (ps0_x * u_y - ps0_y * u_x)		   * xjac * theta * tstep

  amat_n(1,2) = +  eps_cyl * F0 / BigR * v * u_p		   * xjac * theta * tstep

  amat(1,3)   = - eta_numm * (v_x * zj_x + v_y * zj_y)  	   * xjac * theta * tstep &
	      - eta_T * v * zj / BigR				   * xjac * theta * tstep

  amat(1,6)   = - deta_dT * v * T * (zj0 - current_source) / BigR  * xjac * theta * tstep

  
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
subroutine ELM_main_rhs_2(rhs,rhs_k)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_rhs_2

  ! --- Modules
  use phys_module
  use equation_variables
  
  implicit none
  
  ! --- Routine variables
  real*8 :: rhs(n_var),rhs_k(n_var)
  
  ! --- The RHS term	      
  rhs(2) = - 0.5d0 * vv2 * (v_x * r0_y_hat - v_y * r0_x_hat)		 * xjac * tstep &
	   - r0_hat * BigR**2 * w0 * (v_x * u0_y - v_y * u0_x)  	 * xjac * tstep &
	   + v * (ps0_x * zj0_y - ps0_y * zj0_x )			 * xjac * tstep &
	   - visco_T * BigR * (v_x * w0_x + v_y * w0_y) 		 * xjac * tstep &
	   - v * eps_cyl * F0 / BigR * zj0_p				 * xjac * tstep &
	   + BigR**2 * (v_x * p0_y - v_y * p0_x)			 * xjac * tstep &
	   - zeta * BigR * r0_hat * (v_x * delta_u_x + v_y * delta_u_y)  * xjac 	&		  
	   - visco_numm  *								& 
	     (    (	    v_ss  * (x_t**2+y_t**2)					& 
		   +	    v_tt  * (x_s**2+y_s**2)					&
		   - 2.d0 * v_st  * (x_s*x_t + y_s*y_t) )				&
		* (	    w0_ss * (x_t**2+y_t**2)					&
		   +	    w0_tt * (x_s**2+y_s**2)					&
		   - 2.d0 * w0_st * (x_s*x_t + y_s*y_t) ) )				&
	     / xjac**4  						 * xjac * tstep 

  
  return

end subroutine ELM_main_rhs_2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! LHS !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ELM_main_lhs_2(amat, amat_k, amat_n, amat_kn)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_lhs_2

  ! --- Modules
  use phys_module
  use equation_variables
  
  implicit none
  
  ! --- Routine variables
  real*8 :: amat(n_var,n_var), amat_k(n_var,n_var), amat_n(n_var,n_var), amat_kn(n_var,n_var)
  
  ! --- Internal variables
  real*8 :: rho_hat, rho_x_hat, rho_y_hat
  
  rho_hat   = BigR**2 * rho
  rho_x_hat = 2.d0 * BigR * BigR_x  * rho + BigR**2 * rho_x
  rho_y_hat = BigR**2 * rho_y

  ! --- The LHS terms
  amat(2,1)   = - v * (psi_x * zj0_y - psi_y * zj0_x )  				      * xjac * theta * tstep

  amat(2,2)   = + r0_hat * BigR**2 * w0 * (v_x * u_y  - v_y  * u_x)			      * xjac * theta * tstep &
		+ BigR**2 * (u_x*u0_x + u_y*u0_y) * (v_x*r0_y_hat - v_y*r0_x_hat)	      * xjac * theta * tstep 
  
  if (r0 .lt. rho_1) then
    amat(2,2) = amat(2,2) - BigR**3 * rho_1 * (v_x * u_x + v_y * u_y)			      * xjac * (1.d0 + zeta) 
  else
    amat(2,2) = amat(2,2) - BigR**3 * r0    * (v_x * u_x + v_y * u_y)			      * xjac * (1.d0 + zeta) 
  endif

  amat(2,3)   = - v * (ps0_x * zj_y  - ps0_y * zj_x)					      * xjac * theta * tstep

  amat_n(2,3) = + eps_cyl * F0 / BigR * v * zj_p					      * xjac * theta * tstep

  amat(2,4)   = r0_hat * BigR**2 * w  * ( v_x * u0_y - v_y * u0_x)			      * xjac * theta * tstep &
		+ BigR * ( v_x * w_x + v_y * w_y) * visco_T				      * xjac * theta * tstep &
		+ visco_numm  * 										     &
		  (    (	 v_ss * (x_t**2 + y_t**2)							     & 
			+	 v_tt * (x_s**2 + y_s**2)							     &
			- 2.d0 * v_st * (x_s*x_t + y_s*y_t) )							     &
		     * (	 w_ss * (x_t**2 + y_t**2)							     &
			+	 w_tt * (x_s**2 + y_s**2)							     &
			- 2.d0 * w_st * (x_s*x_t + y_s*y_t) ) ) / xjac**4		      * xjac * theta * tstep 

  amat(2,5)   = + 0.5d0 * vv2 * (v_x * rho_y_hat - v_y * rho_x_hat)			      * xjac * theta * tstep &
		+ rho_hat * BigR**2 * w0 * (v_x * u0_y - v_y * u0_x)			      * xjac * theta * tstep &
		- BigR**2 * (v_x * rho_y * T0	- v_y * rho_x * T0  )			      * xjac * theta * tstep &
		- BigR**2 * (v_x * rho   * T0_y - v_y * rho   * T0_x )  		      * xjac * theta * tstep 

  amat(2,6)   = - BigR**2 * (v_x * r0_y * T   - v_y * r0_x * T) 			      * xjac * theta * tstep &
		- BigR**2 * (v_x * r0	* T_y - v_y * r0   * T_x)			      * xjac * theta * tstep &
		+ dvisco_dT * T * ( v_x * w0_x + v_y * w0_y ) * BigR			      * xjac * theta * tstep
  
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
subroutine ELM_main_rhs_3(rhs,rhs_k)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_rhs_3

  ! --- Modules
  use phys_module
  use equation_variables
  
  implicit none
  
  ! --- Routine variables
  real*8 :: rhs(n_var),rhs_k(n_var)
  
  ! --- The RHS term	      
  !rhs(3) = 0.d0
  rhs(3) = - ( v_x * ps0_x  + v_y * ps0_y + v*zj0) / BigR * xjac * tstep

  
  return

end subroutine ELM_main_rhs_3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! LHS !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ELM_main_lhs_3(amat, amat_k, amat_n, amat_kn)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_lhs_3

  ! --- Modules
  use phys_module
  use equation_variables
  
  implicit none
  
  ! --- Routine variables
  real*8 :: amat(n_var,n_var), amat_k(n_var,n_var), amat_n(n_var,n_var), amat_kn(n_var,n_var)
  
  ! --- The LHS terms
  amat(3,1) = (v_x * psi_x + v_y * psi_y ) / BigR     * xjac * tstep

  amat(3,3) = v * zj / BigR			      * xjac * tstep 

  
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
subroutine ELM_main_rhs_4(rhs,rhs_k)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_rhs_4

  ! --- Modules
  use phys_module
  use equation_variables
  
  implicit none
  
  ! --- Routine variables
  real*8 :: rhs(n_var),rhs_k(n_var)
  
  ! --- The RHS term	      
  rhs(4) = 0.d0
  !rhs(4) = - ( v_x * u0_x   + v_y * u0_y  + v*w0)  * BigR * xjac * tstep

  return

end subroutine ELM_main_rhs_4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! LHS !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ELM_main_lhs_4(amat, amat_k, amat_n, amat_kn)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_lhs_4

  ! --- Modules
  use phys_module
  use equation_variables
  
  implicit none
  
  ! --- Routine variables
  real*8 :: amat(n_var,n_var), amat_k(n_var,n_var), amat_n(n_var,n_var), amat_kn(n_var,n_var)
  
  ! --- The LHS terms
  amat(4,2) = (v_x * u_x + v_y * u_y) * BigR  * xjac * tstep

  amat(4,4) =  v * w * BigR		      * xjac * tstep 
		    
  
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
subroutine ELM_main_rhs_5(rhs,rhs_k)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_rhs_5

  ! --- Modules
  use phys_module
  use equation_variables
  
  implicit none
  
  ! --- Routine variables
  real*8 :: rhs(n_var),rhs_k(n_var)
  
  ! --- Parallel gradient terms       
  Bgrad_rho_star   = ( v_x  * ps0_y - v_y  * ps0_x ) / BigR
  Bgrad_rho_k_star = ( F0 / BigR * v_p  	   ) / BigR
	      
  ! --- The RHS term	      
  rhs(5) =  v * BigR * particle_source  				  * xjac * tstep &
	  + v * BigR**2 * ( r0_x * u0_y - r0_y * u0_x)  		  * xjac * tstep &
	  + v * 2.d0 * BigR * r0 * u0_y 				  * xjac * tstep &
	  - (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star * Bgrad_rho	  * xjac * tstep &
	  - D_prof * BigR  * (v_x*r0_x + v_y*r0_y)			  * xjac * tstep &
	  - v * F0 / BigR * Vpar0 * r0_p				  * xjac * tstep &
	  - v * Vpar0 * (r0_x * ps0_y - r0_y * ps0_x)			  * xjac * tstep &
	  - v * F0 / BigR * r0 * vpar0_p				  * xjac * tstep &
	  - v * r0 * (vpar0_x * ps0_y - vpar0_y * ps0_x)		  * xjac * tstep &
	  + zeta * v * BigR						  * xjac *delta_g(5)& 
	  - D_perp_numm *								 &
	    (	 (	   v_ss  * (x_t**2 + y_t**2)					 &
		  +	   v_tt  * (x_s**2 + y_s**2)					 &
		  - 2.d0 * v_st  * (x_s*x_t + y_s*y_t) )				 &
	       * (	   r0_ss * (x_t**2 + y_t**2)					 &
		  +	   r0_tt * (x_s**2 + y_s**2)					 &
		  - 2.d0 * r0_st * (x_s*x_t + y_s*y_t) ) )				 &
	    / xjac**3 *tstep								 

  rhs_k(5) = - (D_par-D_prof) * BigR / BB2 * Bgrad_rho_k_star * Bgrad_rho * xjac * tstep &
	     - D_prof * BigR  * ( v_p*r0_p * eps_cyl**2 /BigR**2 )	  * xjac * tstep 
  
  return

end subroutine ELM_main_rhs_5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! LHS !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ELM_main_lhs_5(amat, amat_k, amat_n, amat_kn)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_lhs_5

  ! --- Modules
  use phys_module
  use equation_variables
  
  implicit none
  
  ! --- Routine variables
  real*8 :: amat(n_var,n_var), amat_k(n_var,n_var), amat_n(n_var,n_var), amat_kn(n_var,n_var)
  
  ! --- Internal variables
  real*8 :: Bgrad_rho_star_psi, Bgrad_rho_psi, Bgrad_rho_rho, Bgrad_rho_rho_n
  
  ! --- Parallel gradient terms       
  Bgrad_rho_star_psi = ( v_x   * psi_y - v_y   * psi_x ) / BigR
  Bgrad_rho_psi      = ( r0_x  * psi_y - r0_y  * psi_x ) / BigR
  Bgrad_rho_rho      = ( rho_x * ps0_y - rho_y * ps0_x ) / BigR
  Bgrad_rho_rho_n    = ( F0 / BigR * rho_p ) / BigR

  ! --- The LHS terms
  amat(5,1)    = - (D_par-D_prof) * BigR * BB2_psi/ BB2**2 * Bgrad_rho_star	* Bgrad_rho	* xjac * theta * tstep &
		 + (D_par-D_prof) * BigR / BB2  	   * Bgrad_rho_star_psi * Bgrad_rho	* xjac * theta * tstep &
		 + (D_par-D_prof) * BigR / BB2  	   * Bgrad_rho_star	* Bgrad_rho_psi * xjac * theta * tstep &
		 + v * Vpar0 * (r0_x * psi_y - r0_y * psi_x)					* xjac * theta * tstep &
		 + v * r0 * (vpar0_x * psi_y - vpar0_y * psi_x) 				* xjac * theta * tstep 

  amat_k(5,1)  = - (D_par-D_prof) * BigR * BB2_psi/ BB2**2 * Bgrad_rho_k_star * Bgrad_rho	* xjac * theta * tstep &
		 + (D_par-D_prof) * BigR / BB2  	   * Bgrad_rho_k_star * Bgrad_rho_psi	* xjac * theta * tstep 

  amat(5,2)    = - v * BigR**2 * ( r0_x * u_y - r0_y * u_x)					* xjac * theta * tstep &
		 - v * 2.d0 * BigR * r0 * u_y							* xjac * theta * tstep 

  amat(5,5)    = + v * rho * BigR								* xjac * (1.d0 + zeta) &
		 - v * BigR**2 * ( rho_x * u0_y - rho_y * u0_x) 				* xjac * theta * tstep &
		 - v * 2.d0 * BigR * rho * u0_y 						* xjac * theta * tstep &
		 + (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star * Bgrad_rho_rho 		* xjac * theta * tstep &
		 + D_prof * BigR  * (v_x*rho_x + v_y*rho_y )					* xjac * theta * tstep &
		 + v * Vpar0 * (rho_x * ps0_y - rho_y * ps0_x)  				* xjac * theta * tstep &
		 + v * rho * (vpar0_x * ps0_y - vpar0_y * ps0_x)				* xjac * theta * tstep &
		 + v * rho * F0 / BigR * vpar0_p						* xjac * theta * tstep &
		 + D_perp_num  *										       &
		   (	(	  v_ss   * (x_t**2 + y_t**2)							       &
			 +	  v_tt   * (x_s**2 + y_s**2)							       &
			 - 2.d0 * v_st   * (x_s*x_t + y_s*y_t) )						       &
		      * (	  rho_ss * (x_t**2 + y_t**2)							       &
			 +	  rho_tt * (x_s**2 + y_s**2)							       &
			 - 2.d0 * rho_st * (x_s*x_t + y_s*y_t) ) ) / xjac**4			* xjac * theta * tstep 

  amat_k(5,5)  = + (D_par-D_prof) * BigR / BB2 * Bgrad_rho_k_star * Bgrad_rho_rho		* xjac * theta * tstep 

  amat_n(5,5)  = + (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star   * Bgrad_rho_rho_n		* xjac * theta * tstep &
		 + v * F0 / BigR * Vpar0 * rho_p						* xjac * theta * tstep 

  amat_kn(5,5) = + (D_par-D_prof) * BigR / BB2 * Bgrad_rho_k_star * Bgrad_rho_rho_n		* xjac * theta * tstep &
		 + D_prof * BigR  * ( v_p*rho_p * eps_cyl**2 /BigR**2 ) 			* xjac * theta * tstep 

  amat(5,7)    = + v * F0 / BigR * Vpar * r0_p  						* xjac * theta * tstep &
		 + v * Vpar * (r0_x * ps0_y - r0_y * ps0_x)					* xjac * theta * tstep &
		 + v * r0 * (vpar_x * ps0_y - vpar_y * ps0_x)					* xjac * theta * tstep 

  amat_n(5,7)  = + v * r0 * F0 / BigR * vpar_p  						* xjac * theta * tstep
		    
  
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
subroutine ELM_main_rhs_6(rhs,rhs_k)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_rhs_6

  ! --- Modules
  use phys_module
  use equation_variables
  
  implicit none
  
  ! --- Routine variables
  real*8 :: rhs(n_var),rhs_k(n_var)
  
  ! --- Parallel gradient terms       
  v_ps0_x  = v_xx  * ps0_y - v_xy  * ps0_x + v_x  * ps0_xy - v_y * ps0_xx
  v_ps0_y  = v_xy  * ps0_y - v_yy  * ps0_x + v_x  * ps0_yy - v_y * ps0_xy    
  Bgrad_T_star     = ( v_x  * ps0_y - v_y  * ps0_x ) / BigR
  Bgrad_T_k_star   = ( F0 / BigR * v_p  	   ) / BigR
	      
  ! --- The RHS term	      
  rhs(6) =   v * BigR * heat_source					       * xjac * tstep &
	   + v * r0 * BigR**2 * ( T0_x * u0_y - T0_y * u0_x)		       * xjac * tstep &
	   + v * r0 * 2.d0* (GAMMA-1.d0) * BigR * T0 * u0_y		       * xjac * tstep &
	   - v * r0 * F0 / BigR * Vpar0 * T0_p  			       * xjac * tstep &
	   - v * r0 * Vpar0 * (T0_x * ps0_y - T0_y * ps0_x)		       * xjac * tstep &
	   - v * r0 * (GAMMA-1.d0) * T0 * (vpar0_x * ps0_y - vpar0_y * ps0_x)  * xjac * tstep &
	   - v * r0 * (GAMMA-1.d0) * T0 * F0 / BigR * vpar0_p		       * xjac * tstep &
	   -  (K_par-ZK_prof) * BigR / BB2 * Bgrad_T_star * Bgrad_T	       * xjac * tstep &
	   -  ZK_prof * BigR * (v_x*T0_x + v_y*T0_y )			       * xjac * tstep &
	   -  K_par_num * (v_ps0_x  * ps0_y - v_ps0_y  * ps0_x) &
			* (T0_ps0_x * ps0_y - T0_ps0_y * ps0_x) 	       * xjac * tstep &
	   + zeta * v * r0 * BigR					       * xjac *delta_g(6)&
	   - K_perp_numm2  *								      &
	     (    (	    v_ss  * (x_t**2 + y_t**2)					      &
		   +	    v_tt  * (x_s**2 + y_s**2)					      &
		   - 2.d0 * v_st  * (x_s*x_t + y_s*y_t) )				      &
		* (	    T3_ss * (x_t**2 + y_t**2)					      &
		   +	    T3_tt * (x_s**2 + y_s**2)					      &
		   - 2.d0 * T3_st * (x_s*x_t + y_s*y_t) ) )				      &
	     / xjac**4  						       * xjac * tstep 

  rhs_k(6) = - (K_par-ZK_prof) * BigR / BB2 * Bgrad_T_k_star * Bgrad_T         * xjac * tstep &
	     -  K_par_num2 * (T0_pp  * v_pp)				       * xjac * tstep &
	     - ZK_prof * BigR * ( v_p*T0_p /BigR**2 )			       * xjac * tstep 

  
  return

end subroutine ELM_main_rhs_6


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! LHS !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ELM_main_lhs_6(amat, amat_k, amat_n, amat_kn)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_lhs_6

  ! --- Modules
  use phys_module
  use equation_variables
  
  implicit none
  
  ! --- Routine variables
  real*8 :: amat(n_var,n_var), amat_k(n_var,n_var), amat_n(n_var,n_var), amat_kn(n_var,n_var)
  
  ! --- Internal variables
  real*8 :: Bgrad_T_star_psi, Bgrad_T_psi, Bgrad_T_T, Bgrad_T_T_n
  real*8 :: T_ps0_x,  T_ps0_y
  real*8 :: T0_psi_x, T0_psi_y
  real*8 :: v_psi_x,  v_psi_y
  
  ! --- Parallel gradient terms       
  Bgrad_T_star_psi   = ( v_x   * psi_y - v_y   * psi_x ) / BigR
  Bgrad_T_psi	     = ( T0_x  * psi_y - T0_y  * psi_x ) / BigR
  Bgrad_T_T	     = ( T_x   * ps0_y - T_y   * ps0_x ) / BigR
  Bgrad_T_T_n	     = ( F0 / BigR * T_p   ) / BigR

  T_ps0_x = T_xx * ps0_y - T_xy * ps0_x + T_x * ps0_xy - T_y * ps0_xx
  T_ps0_y = T_xy * ps0_y - T_yy * ps0_x + T_x * ps0_yy - T_y * ps0_xy
  
  T0_psi_x = T0_xx * psi_y - T0_xy * psi_x + T0_x * psi_xy - T0_y * psi_xx
  T0_psi_y = T0_xy * psi_y - T0_yy * psi_x + T0_x * psi_yy - T0_y * psi_xy
  
  v_psi_x = v_xx * psi_y - v_xy * psi_x + v_x * psi_xy - v_y * psi_xx
  v_psi_y = v_xy * psi_y - v_yy * psi_x + v_x * psi_yy - v_y * psi_xy

  ! --- The LHS terms
  amat(6,1)    = + v * r0 * Vpar0 * (T0_x * psi_y - T0_y * psi_x)				  * xjac * theta * tstep &
		 + v * (GAMMA-1.d0) * r0 * T0 * (vpar0_x * psi_y - vpar0_y * psi_x)		  * xjac * theta * tstep &
		 + K_par_num * (v_psi_x *ps0_y - v_psi_y *ps0_x + v_ps0_x *psi_y - v_ps0_y *psi_x)			 &
			     * (T0_ps0_x*ps0_y - T0_ps0_y*ps0_x)				  * xjac * theta * tstep &
		 + K_par_num * (T0_psi_x*ps0_y - T0_psi_y*ps0_x + T0_ps0_x*psi_y - T0_ps0_y*psi_x)			 &
			     * (v_ps0_x *ps0_y - v_ps0_y *ps0_x)				  * xjac * theta * tstep &
		 - (K_par-ZK_prof) * BigR * BB2_psi / BB2**2 * Bgrad_T_star	* Bgrad_T	  * xjac * theta * tstep &
		 + (K_par-ZK_prof) * BigR / BB2 	     * Bgrad_T_star_psi * Bgrad_T	  * xjac * theta * tstep &
		 + (K_par-ZK_prof) * BigR / BB2 	     * Bgrad_T_star	* Bgrad_T_psi	  * xjac * theta * tstep 

  amat_k(6,1)  = - (K_par-ZK_prof) * BigR * BB2_psi / BB2**2 * Bgrad_T_k_star	* Bgrad_T	  * xjac * theta * tstep &
		 + (K_par-ZK_prof) * BigR / BB2 	     * Bgrad_T_k_star	* Bgrad_T_psi	  * xjac * theta * tstep 

  amat(6,2)    = - v * r0 * BigR**2 * ( T0_x * u_y - T0_y * u_x)				  * xjac * theta * tstep &
		 - v * 2.d0 * (GAMMA-1.d0) * r0 * BigR * T0 * u_y				  * xjac * theta * tstep 

  amat(6,5)    = - v * rho * BigR**2 * (T0_x * u0_y - T0_y * u0_x)				  * xjac * theta * tstep &
		 + v * rho * Vpar0 * F0/BigR * T0_p						  * xjac * theta * tstep &
		 + v * rho * Vpar0 * (T0_x * ps0_y - T0_y * ps0_x)				  * xjac * theta * tstep &
		 - v * 2.d0 * (GAMMA-1.d0) * rho * BigR * T0 * u0_y				  * xjac * theta * tstep &
		 + v * (GAMMA-1.d0) * rho * T0 * F0/BigR * Vpar0_p				  * xjac * theta * tstep &
		 + v * (GAMMA-1.d0) * rho * T0 * (Vpar0_x * ps0_y - Vpar0_y * ps0_x)		  * xjac * theta * tstep 

  amat(6,6)    = - v * r0 * BigR**2 * (T_x * u0_y  - T_y * u0_x)				  * xjac * theta * tstep &
		 + v * r0 * Vpar0   * (T_x * ps0_y - T_y * ps0_x)				  * xjac * theta * tstep &
		 - 2.d0 * v * r0 * (GAMMA-1.d0) * T * BigR * u0_y				  * xjac * theta * tstep &
		 + v * r0 * (GAMMA-1.d0) * T * F0/BigR * Vpar0_p				  * xjac * theta * tstep &
		 + v * r0 * (GAMMA-1.d0) * T * (Vpar0_x * ps0_y - Vpar0_y * ps0_x)		  * xjac * theta * tstep &
		 + K_par_num * (v_ps0_x * ps0_y - v_ps0_y * ps0_x) &
			     * (T_ps0_x * ps0_y - T_ps0_y * ps0_x)				  * xjac * theta * tstep &
		 + (K_par-ZK_prof) * BigR / BB2 * Bgrad_T_star * Bgrad_T_T			  * xjac * theta * tstep &
		 + dZK_par * T     * BigR / BB2 * Bgrad_T_star * Bgrad_T			  * xjac * theta * tstep &
		 + ZK_prof * BigR * (v_x*T_x + v_y*T_y )					  * xjac * theta * tstep & 
		 + K_perp_numm  *											 &
		   (	(	  v_ss * (x_t**2 + y_t**2)								 &
			 +	  v_tt * (x_s**2 + y_s**2)								 &
			 - 2.d0 * v_st * (x_s*x_t + y_s*y_t) )  							 &
		      * (	  T_ss * (x_t**2 + y_t**2)								 &
			 +	  T_tt * (x_s**2 + y_s**2)								 &
			 - 2.d0 * T_st * (x_s*x_t + y_s*y_t) ) )     / xjac**4  		  * xjac * theta * tstep 
  
  if (r0 .lt. rho_1) then
    amat(6,6)  = amat(6,6) + v * rho_1 * T * BigR						  * xjac * (1.d0 + zeta)
  else
    amat(6,6)  = amat(6,6) + v * r0    * T * BigR						  * xjac * (1.d0 + zeta)
  endif
  
  !if ((im .ne. 1) .and. (in .ne. 1)) then
  !  amat(6,6)  = amat(6,6) + K_perp_numm2  *									       &
  !			(    (       v_ss * (x_t**2 + y_t**2)							       &
  !			    +	     v_tt * (x_s**2 + y_s**2)							       &
  !			    - 2.d0 * v_st * (x_s*x_t + y_s*y_t) )						       &
  !			   * (       T_ss * (x_t**2 + y_t**2)							       &
  !			    +	     T_tt * (x_s**2 + y_s**2)							       &
  !			    - 2.d0 * T_st * (x_s*x_t + y_s*y_t) ) )	/ xjac**4		* xjac * theta * tstep 
  !endif

  amat_k(6,6)  = + (K_par-ZK_prof) * BigR / BB2 * Bgrad_T_k_star * Bgrad_T_T			  * xjac * theta * tstep &
		 + dZK_par * T     * BigR / BB2 * Bgrad_T_k_star * Bgrad_T			  * xjac * theta * tstep 
	      
  amat_n(6,6)  = + (K_par-ZK_prof) * BigR / BB2 * Bgrad_T_star   * Bgrad_T_T_n  		  * xjac * theta * tstep &
		 + v * r0 * Vpar0  * F0/BigR * T_p						  * xjac * theta * tstep 

  amat_kn(6,6) = + (K_par-ZK_prof) * BigR / BB2 * Bgrad_T_k_star * Bgrad_T_T_n  		  * xjac * theta * tstep &
		 +  K_par_num2 * (T_pp  * v_pp) 						  * xjac * theta * tstep &
		 + ZK_prof	   * BigR	* (v_p*T_p /BigR**2 )				  * xjac * theta * tstep 

  amat(6,7)    = + v * r0 * F0/BigR * Vpar * T0_p						  * xjac * theta * tstep &
		 + v * r0 * Vpar * (T0_x * ps0_y - T0_y * ps0_x)				  * xjac * theta * tstep &
		 + v * (GAMMA-1.d0) * r0 * T0 * (Vpar_x * ps0_y - Vpar_y * ps0_x)		  * xjac * theta * tstep 
     
  amat_n(6,7)  = v * (GAMMA-1.d0) * r0 * T0 * F0/BigR * Vpar_p  				  * xjac * theta * tstep       
		   

  
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
subroutine ELM_main_rhs_7(rhs,rhs_k)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_rhs_7

  ! --- Modules
  use phys_module
  use equation_variables
  
  implicit none
  
  ! --- Routine variables
  real*8 :: rhs(n_var),rhs_k(n_var)
  
  ! --- The RHS term	      
  rhs(7) = - v * F0 / BigR * P0_p				  * xjac * tstep &
	   - v * (P0_x * ps0_y - P0_y * ps0_x)  		  * xjac * tstep &
	   - visco_par * (v_x * vpar0_x + v_y * vpar0_y) * BigR   * xjac * tstep &
	   + zeta * v * r0 * F0**2 / BigR			  * xjac *delta_g(7)&  
	   - visco_par_numm  *  						 &
	     (    (	    v_ss     * (x_t**2 + y_t**2)			 &
		   +	    v_tt     * (x_s**2 + y_s**2)			 &
		   - 2.d0 * v_st     * (x_s*x_t + y_s*y_t) )			 &
		* (	    vpar0_ss * (x_t**2 + y_t**2)			 &
		   +	    vpar0_tt * (x_s**2 + y_s**2)			 &
		   - 2.d0 * vpar0_st * (x_s*x_t + y_s*y_t) ) )  		 &
	     / xjac**4  					  * xjac * tstep 
  
  return

end subroutine ELM_main_rhs_7

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! LHS !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ELM_main_lhs_7(amat, amat_k, amat_n, amat_kn)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_lhs_7

  ! --- Modules
  use phys_module
  use equation_variables
  
  implicit none
  
  ! --- Routine variables
  real*8 :: amat(n_var,n_var), amat_k(n_var,n_var), amat_n(n_var,n_var), amat_kn(n_var,n_var)
  
  ! --- The LHS terms
  amat(7,1)   = + v * (P0_x * psi_y - P0_y * psi_x)				   * xjac * theta * tstep !&
	      !+ vpar0 * (F0/BigR)**2 * (vpar0_x * ps_y - vpar0_y * ps_x)	   * xjac * theta * tstep &

  amat(7,2)   = 0.d0!+ F0 * (u_s * vpar0_t - u_t * vpar0_s) * theta * tstep 

  amat(7,5)   = + v * (rho_x * T0   * ps0_y - rho_y * T0   * ps0_x)		   * xjac * theta * tstep &
		+ v * (rho   * T0_x * ps0_y - rho   * T0_y * ps0_x)		   * xjac * theta * tstep &
		+ v * F0 / BigR * rho * T0_p					   * xjac * theta * tstep 

  amat_n(7,5) = v * F0 / BigR * rho_p * T0					   * xjac * theta * tstep  

  amat(7,6)   = + v * (T_x * r0 	* ps0_y - T_y * r0   * ps0_x)		   * xjac * theta * tstep &
		+ v * (T   * r0_x * ps0_y - T	* r0_y * ps0_x) 		   * xjac * theta * tstep &
		+ v * F0 / BigR * T * r0_p					   * xjac * theta * tstep 
  
  amat_n(7,6) = v * F0 / BigR * T_p * r0					   * xjac * theta * tstep

  amat(7,7)   = + visco_par * (v_x * Vpar_x + v_y * Vpar_y) * BigR		   * xjac * theta * tstep &
		+ visco_par_numm  *									  &
		  (    (	v_ss	* (x_t**2 + y_t**2)						  &
			+	v_tt	* (x_s**2 + y_s**2)						  &
			-2.d0 * v_st	* (x_s*x_t + y_s*y_t) ) 					  &
		     * (	Vpar_ss * (x_t**2 + y_t**2)						  &
			+	Vpar_tt * (x_s**2 + y_s**2)						  &
			-2.d0 * Vpar_st * (x_s*x_t + y_s*y_t) ) ) / xjac**4	   * xjac * theta * tstep !&
		!+ F0 * (u0_x * vpar_y - u0_y * vpar_x) 			   * xjac * theta * tstep &
		!+ vpar * (F0/BigR)**2  * (vpar0_x * ps0_y - vpar0_y * ps0_x)	   * xjac * theta * tstep &
		!+ vpar0 * (F0/BigR)**2 * (vpar_x  * ps0_y - vpar_y * ps0_x)	   * xjac * theta * tstep &
		!+ (F0/BigR)**3 * vpar  * vpar0_p				   * xjac * theta * tstep &
		!+ (F0/BigR)**3 * vpar0 * vpar_p				   * xjac * theta * tstep &
  
  if (r0 .lt. rho_1) then
    amat(7,7) = amat(7,7) + v * Vpar * rho_1 * F0**2 / BigR			   * xjac * (1.d0 + zeta)
  else
    amat(7,7) = amat(7,7) + v * Vpar * r0    * F0**2 / BigR			   * xjac * (1.d0 + zeta)
  endif

  
  return

end subroutine ELM_main_lhs_7



