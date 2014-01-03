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
  rhs(1) = + v * eta_Te  * (zj0 - current_source)/ BigR				* xjac * tstep &
	   + v * (ps0_x * u0_y - ps0_y * u0_x)  				* xjac * tstep &
           - v * tau_IC/(r0*BB2) * F0**2/BigR**2 * (ps0_x*pe0_y - ps0_y*pe0_x)	* xjac * tstep & 		  
           + v * tau_IC/(r0*BB2) * F0**3/BigR**3 * eps_cyl * pe0_p		* xjac * tstep &
	   - v * eps_cyl * F0 / BigR  * u0_p					* xjac * tstep &
	   + eta_numm * (v_x * zj0_x + v_y * zj0_y)				* xjac * tstep &
	   + zeta * v / BigR							* xjac * delta_g(1)
  
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
  amat(1,1)   = + v * psi / BigR									* xjac * (1.d0+zeta)   &
		- v * (psi_x * u0_y - psi_y * u0_x)							* xjac * theta * tstep &
                - v * tau_IC/(r0*BB2**2) * BB2_psi * F0**2/BigR**2 * (ps0_x*pe0_y - ps0_y*pe0_x)	* xjac * theta * tstep &
                + v * tau_IC/(r0*BB2**2) * BB2_psi * F0**3/BigR**3 * eps_cyl * pe0_p			* xjac * theta * tstep &
                + v * tau_IC/(r0*BB2) * F0**2/BigR**2 * (psi_x * pe0_y - psi_y * pe0_x)			* xjac * theta * tstep

  amat(1,2)   = -  v * (ps0_x * u_y - ps0_y * u_x)							* xjac * theta * tstep

  amat_n(1,2) = +  eps_cyl * F0 / BigR * v * u_p							* xjac * theta * tstep

  amat(1,3)   = - eta_numm * (v_x * zj_x + v_y * zj_y)							* xjac * theta * tstep &
	      - eta_Te * v * zj / BigR									* xjac * theta * tstep

  amat(1,5)   = + v * tau_IC/(r0*BB2) * F0**2/BigR**2 * Te0 * (ps0_x*rho_y - ps0_y*rho_x)		* xjac * theta * tstep &
	        + v * tau_IC/(r0*BB2) * F0**2/BigR**2 * rho * (ps0_x*Te0_y - ps0_y*Te0_x)		* xjac * theta * tstep &
	        - v * tau_IC/(r0*BB2) * F0**3/BigR**3 * eps_cyl * rho * Te0_p				* xjac * theta * tstep &
		- v * tau_IC * rho /(r0**2 * BB2) * F0**2/BigR**2 * (ps0_x*pe0_y - ps0_y*pe0_x)		* xjac * theta * tstep &		    
		+ v * tau_IC * rho /(r0**2 * BB2) * F0**3/BigR**3 * eps_cyl * pe0_p			* xjac * theta * tstep 
  
  amat_n(1,5) = - v * tau_IC/(r0*BB2) * F0**3/BigR**3 * eps_cyl * Te0 * rho_p				* xjac * theta * tstep
  
  amat(1,8)   = - deta_dTe * v * Te * (zj0 - current_source) / BigR					* xjac * theta * tstep &
                + v * tau_IC/(r0*BB2) * F0**2/BigR**2 * r0 * (ps0_x*Te_y - ps0_y*Te_x) 			* xjac * theta * tstep &
	        + v * tau_IC/(r0*BB2) * F0**2/BigR**2 * Te * (ps0_x*r0_y - ps0_y*r0_x)			* xjac * theta * tstep &
		- v * tau_IC/(r0*BB2) * F0**3/BigR**3 * eps_cyl * Te * r0_p				* xjac * theta * tstep 

  amat_n(1,8) = - v * tau_IC/(r0*BB2) * F0**3/BigR**3 * eps_cyl * r0 * Te_p				* xjac * theta * tstep 

  
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
  real*8 :: rhs(n_var),rhs_k(n_var), mi_e
  
  ! ----------------------------	      
  ! --- The RHS term (main part)	      
  rhs(2) = 											       &
           ! --- Time derivative
	   - zeta * BigR * r0_hat * (v_x * delta_u_x + v_y * delta_u_y)  		* xjac 	       &		  
           ! --- Convective terms
	   - 0.5d0 * vv2 * (v_x * r0_y_hat - v_y * r0_x_hat)		 		* xjac * tstep &
	   - r0_hat * BigR**2 * w0 * (v_x * u0_y - v_y * u0_x)  	 		* xjac * tstep &
           ! --- [psi,j]
	   + v * (ps0_x * zj0_y - ps0_y * zj0_x )			 		* xjac * tstep &
	   - v * eps_cyl * F0 / BigR * zj0_p				 		* xjac * tstep &
           ! --- Grad(p)
	   + BigR**2 * (v_x * p0_y     - v_y * p0_x)			 		* xjac * tstep &
           ! --- Source and viscosity
	   + BigR**3 * particle_source * (v_x * u0_x + v_y * u0_y)			* xjac * tstep &
	   - visco_Te * BigR * (v_x * w0_x + v_y * w0_y) 		 		* xjac * tstep &
	   - visco_numm  *									       & 
	     (    (	    v_ss  * (x_t**2+y_t**2)						       & 
		   +	    v_tt  * (x_s**2+y_s**2)						       &
		   - 2.d0 * v_st  * (x_s*x_t + y_s*y_t) )					       &
		* (	    w0_ss * (x_t**2+y_t**2)						       &
		   +	    w0_tt * (x_s**2+y_s**2)						       &
		   - 2.d0 * w0_st * (x_s*x_t + y_s*y_t) ) )					       &
	     / xjac**4  						 		* xjac * tstep 
  
  ! --------------------------------------------------      
  ! --- The RHS term (diamagnetic and neoclassic part)	      
  mi_e = MASS_PROTON / EL_CHG
  rhs(2) = rhs(2)										       &
	   ! --- Main diamagnetic terms
	   - tau_IC * v * BigR**4        * (Pi0_x * w0_y - Pi0_y * w0_x)		* xjac * tstep &
	   - tau_IC     * BigR**3 * Pi0_y * (v_x  * u0_x + v_y  * u0_y)			* xjac * tstep &
	   - tau_IC * v * BigR**4 * (u0_xy * (Pi0_xx-Pi0_yy) - Pi0_xy * (u0_xx-u0_yy) )	* xjac * tstep &
	   - visco_Te * BigR * (v_x * Wdia0_x + v_y * Wdia0_y) 		 		* xjac * tstep &
	   ! --- Second diamagnetic terms (RHS of gyro-viscous cancellation)
	   - mi_e * F0 / BB2 * BigR**2 * (v_x * W0*Pi0_y - v_y * W0*Pi0_x)		* xjac * tstep &
	   - mi_e * F0 / BB2 * BigR**2 * (v_x * W0_y*Pi0 - v_y * W0_x*Pi0)		* xjac * tstep &
	   ! --- Neoclassic term
           + amu_neo_prof * BB2 / (Btheta2+epsil)**2.d0 * (ps0_x*v_x + ps0_y*v_y) * BigR	       &
                    * (  r0                         * (ps0_x*u0_x  + ps0_y*u0_y)		       &
		       + tau_IC                     * (ps0_x*Pi0_x + ps0_y*Pi0_y)		       &
                       + aki_neo_prof * tau_IC * r0 * (ps0_x*Ti0_x + ps0_y*Ti0_y)		       &
                       - r0 * Vpar0 * Btheta2					     )	* xjac * tstep 

  
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
  real*8 :: Btheta2_psi, mi_e
  
  rho_hat     = BigR**2 * rho
  rho_x_hat   = 2.d0 * BigR * BigR_x  * rho + BigR**2 * rho_x
  rho_y_hat   = BigR**2 * rho_y
  Btheta2_psi = 2.d0 * (psi_x * ps0_x + psi_y * ps0_y ) / BigR**2

  ! -----------------------------	      
  ! --- The LHS terms (main part)
  amat(2,1)   = - v * (psi_x * zj0_y - psi_y * zj0_x )								* xjac * theta * tstep

  amat(2,2)   = + r0_hat * BigR**2 * w0 * (v_x * u_y  - v_y  * u_x)						* xjac * theta * tstep &
		+ BigR**2 * (u_x*u0_x + u_y*u0_y) * (v_x*r0_y_hat - v_y*r0_x_hat)				* xjac * theta * tstep &
		- BigR**3 * particle_source * (v_x * u_x + v_y * u_y)						* xjac * theta * tstep
  
  if (r0 .lt. rho_1) then
    amat(2,2) = amat(2,2) - BigR**3 * rho_1 * (v_x * u_x + v_y * u_y)						* xjac * (1.d0 + zeta) 
  else
    amat(2,2) = amat(2,2) - BigR**3 * r0    * (v_x * u_x + v_y * u_y)						* xjac * (1.d0 + zeta) 
  endif

  amat(2,3)   = - v * (ps0_x * zj_y  - ps0_y * zj_x)								* xjac * theta * tstep

  amat_n(2,3) = + eps_cyl * F0 / BigR * v * zj_p								* xjac * theta * tstep

  amat(2,4)   = r0_hat * BigR**2 * w  * ( v_x * u0_y - v_y * u0_x)						* xjac * theta * tstep &
		+ visco_Te * BigR * ( v_x * w_x + v_y * w_y)							* xjac * theta * tstep &
		+ visco_numm  * 												       &
		  (    (	 v_ss * (x_t**2 + y_t**2)									       & 
			+	 v_tt * (x_s**2 + y_s**2)									       &
			- 2.d0 * v_st * (x_s*x_t + y_s*y_t) )									       &
		     * (	 w_ss * (x_t**2 + y_t**2)									       &
			+	 w_tt * (x_s**2 + y_s**2)									       &
			- 2.d0 * w_st * (x_s*x_t + y_s*y_t) ) ) / xjac**4					* xjac * theta * tstep 

  amat(2,5)   = + 0.5d0 * vv2 * (v_x * rho_y_hat - v_y * rho_x_hat)						* xjac * theta * tstep &
		+ rho_hat * BigR**2 * w0 * (v_x * u0_y - v_y * u0_x)						* xjac * theta * tstep &
		- BigR**2 * (v_x * rho_y * T0   - v_y * rho_x * T0  )						* xjac * theta * tstep &
		- BigR**2 * (v_x * rho   * T0_y - v_y * rho   * T0_x)						* xjac * theta * tstep

  amat(2,6)   = - BigR**2 * (v_x * r0_y * Ti   - v_y * r0_x * Ti)						* xjac * theta * tstep &
		- BigR**2 * (v_x * r0	* Ti_y - v_y * r0   * Ti_x)						* xjac * theta * tstep 
  
  amat(2,8)   = - BigR**2 * (v_x * r0_y * Te   - v_y * r0_x * Te)						* xjac * theta * tstep &
		- BigR**2 * (v_x * r0	* Te_y - v_y * r0   * Te_x)						* xjac * theta * tstep &
		+ dvisco_dTe * Te * ( v_x * w0_x    + v_y * w0_y    ) * BigR					* xjac * theta * tstep &
		+ dvisco_dTe * Te * ( v_x * Wdia0_x + v_y * Wdia0_y ) * BigR					* xjac * theta * tstep 
  
  
  ! ---------------------------------------------------    
  ! --- The LHS terms (diamagnetic and neoclassic part)
  mi_e = MASS_PROTON / EL_CHG
  amat(2,1)   = amat(2,1)													       &
	        - mi_e * F0 * BB2_psi/BB2**2.d0 * BigR**2 * (v_x * W0*Pi0_y - v_y * W0*Pi0_x)			* xjac * theta * tstep &
	        - mi_e * F0 * BB2_psi/BB2**2.d0 * BigR**2 * (v_x * W0_y*Pi0 - v_y * W0_x*Pi0)			* xjac * theta * tstep &
                - amu_neo_prof * BB2 / (Btheta2+epsil)**2 * (psi_x*v_x+psi_y*v_y) * BigR					       &
                               * (  r0                         * (ps0_x*u0_x  + ps0_y*u0_y)					       &
                                  + tau_IC                     * (ps0_x*Pi0_x + ps0_y*Pi0_y)					       &
                                  + aki_neo_prof * tau_IC * r0 * (ps0_x*Ti0_x + ps0_y*Ti0_y)					       &
                                  - r0 * Vpar0 * Btheta2)	 						* xjac * theta * tstep &
                - amu_neo_prof * BB2 / (Btheta2+epsil)**2 * (ps0_x*v_x+ps0_y*v_y) * BigR					       &
		               * (  r0                         * (psi_x*u0_x  + psi_y*u0_y)					       &
                                  + tau_IC                     * (psi_x*Pi0_x + psi_y*Pi0_y)					       &
                                  + aki_neo_prof * tau_IC * r0 * (psi_x*Ti0_x + psi_y*Ti0_y)	)		* xjac * theta * tstep &
                + amu_neo_prof * BB2 * 2.d0*Btheta2_psi / (Btheta2+epsil)**3 * (ps0_x*v_x+ps0_y*v_y) * BigR			       &
                               * (  r0                         * (ps0_x*u0_x  + ps0_y*u0_y)					       &
                                  + tau_IC                     * (ps0_x*Pi0_x + ps0_y*Pi0_y)					       &
                                  + aki_neo_prof * tau_IC * r0 * (ps0_x*Ti0_x + ps0_y*Ti0_y)	)		* xjac * theta * tstep &
                - amu_neo_prof * BB2 * Btheta2_psi / (Btheta2+epsil)**2								       &
		               * r0 * vpar0 * (ps0_x*v_x + ps0_y*v_y) * BigR					* xjac * tstep * theta

  amat(2,2)   = amat(2,2)													       &
	        + tau_IC * BigR**3 * Pi0_y * (v_x* u_x + v_y * u_y)						* xjac * theta * tstep &
	        + tau_IC * v * BigR**4 * (u_xy * (Pi0_xx-Pi0_yy) - Pi0_xy * (u_xx-u_yy))			* xjac * theta * tstep &
                - amu_neo_prof * BB2 / (Btheta2+epsil)**2 * (ps0_x*v_x + ps0_y*v_y) * BigR					       &
		               * r0 * (ps0_x*u_x + ps0_y*u_y)							* xjac * theta * tstep
  
  amat(2,4)   = amat(2,4)													       &
	        + mi_e * F0 / BB2 * BigR**2 * (v_x * W*Pi0_y - v_y * W*Pi0_x)			 		* xjac * theta * tstep &
	        + mi_e * F0 / BB2 * BigR**2 * (v_x * W_y*Pi0 - v_y * W_x*Pi0)			 		* xjac * theta * tstep &
                + tau_IC * v * BigR**4 * (Pi0_x * w_y - Pi0_y * w_x)              				* xjac * theta * tstep 

  amat(2,5)   = amat(2,5)													       &
	        + mi_e * F0 / BB2 * BigR**2 * (v_x * W0*rho*Ti0_y - v_y * W0*rho*Ti0_x)			 	* xjac * theta * tstep &
	        + mi_e * F0 / BB2 * BigR**2 * (v_x * W0*rho_y*Ti0 - v_y * W0*rho_x*Ti0)			 	* xjac * theta * tstep &
	        + mi_e * F0 / BB2 * BigR**2 * (v_x * W0_y*rho*Ti0 - v_y * W0_x*rho*Ti0)			 	* xjac * theta * tstep &
                + tau_IC * v * BigR**4 * Ti0  * (rho_x  * w0_y - rho_y  * w0_x)  				* xjac * theta * tstep &
                + tau_IC * v * BigR**4 * rho  * (Ti0_x  * w0_y - Ti0_y  * w0_x)  				* xjac * theta * tstep &
		+ tau_IC     * BigR**3 * (Ti0_y*rho + Ti0*rho_y) * (v_x*u0_x + v_y*u0_y) 			* xjac * theta * tstep &
		+ tau_IC * v * BigR**4 * ( u0_xy        * (rho_xx*Ti0 + 2.d0*rho_x*Ti0_x + rho*Ti0_xx  				       &
			    	                          -rho_yy*Ti0 - 2.d0*rho_y*Ti0_y - rho*Ti0_yy)	   			       &					   
			                 -(u0_xx-u0_yy) * (rho_xy*Ti0 + rho_x*Ti0_y + rho_y*Ti0_x + rho*Ti0_xy))* xjac * theta * tstep &
                - amu_neo_prof * BB2 / (Btheta2+epsil)**2 * (ps0_x*v_x + ps0_y*v_y) * BigR					       &
                               * (  rho                         * (ps0_x*u0_x                  + ps0_y*u0_y)			       &
                                  + tau_IC                      * (ps0_x*rho_x*Ti0   + ps0_y*rho_y*Ti0  )			       &
                                  + tau_IC                      * (ps0_x*rho  *Ti0_x + ps0_y*rho  *Ti0_y)			       &
                                  + aki_neo_prof * tau_IC * rho * (ps0_x*Ti0_x                 + ps0_y*Ti0_y)			       &
                                  -rho * Vpar0 * Btheta2						      )	* xjac * tstep * theta

  amat(2,6)   = amat(2,6)													       &
	        + mi_e * F0 / BB2 * BigR**2 * (v_x * W0*r0*Ti_y - v_y * W0*r0*Ti_x)			 	* xjac * theta * tstep &
	        + mi_e * F0 / BB2 * BigR**2 * (v_x * W0*r0_y*Ti - v_y * W0*r0_x*Ti)			 	* xjac * theta * tstep &
	        + mi_e * F0 / BB2 * BigR**2 * (v_x * W0_y*r0*Ti - v_y * W0_x*r0*Ti)			 	* xjac * theta * tstep &
                + tau_IC * v * BigR**4 * r0 * (Ti_x *w0_y - Ti_y *w0_x)  					* xjac * theta * tstep &
                + tau_IC * v * BigR**4 * Ti  * (r0_x*w0_y - r0_y*w0_x)  					* xjac * theta * tstep &
		+ tau_IC     * BigR**3 * (r0_y*Ti + r0*Ti_y) * (v_x*u0_x + v_y*u0_y) 				* xjac * theta * tstep &
		+ tau_IC * v * BigR**4 * ( u0_xy        * (r0_xx*Ti + 2.d0*r0_x*Ti_x + r0*Ti_xx      				       &
				                          -r0_yy*Ti - 2.d0*r0_y*Ti_y - r0*Ti_yy)				       &					 
			                 -(u0_xx-u0_yy) * (r0_xy*Ti + r0_x*Ti_y + r0_y*Ti_x + r0*Ti_xy)	)	* xjac * theta * tstep &
                - amu_neo_prof * BB2 / (Btheta2+epsil)**2 * (ps0_x*v_x + ps0_y*v_y) * BigR					       &
                               * (  tau_IC                     * (ps0_x*r0_x*Ti   + ps0_y*r0_y*Ti  )				       &
                                  + tau_IC                     * (ps0_x*r0  *Ti_x + ps0_y*r0  *Ti_y)				       &
                                  + aki_neo_prof * tau_IC * r0 * (ps0_x*Ti_x+ps0_y*Ti_y)		   )	* xjac * tstep * theta
  
  amat(2,7)   = amat(2,7)													       &
                + amu_neo_prof * BB2 / (Btheta2+epsil) * r0 * vpar * (ps0_x*v_x + ps0_y*v_y) * BigR		* xjac * tstep * theta 
  
  amat(2,9)   = amat(2,9)													       &
                + visco_Te * BigR * ( v_x * Wdia_x + v_y * Wdia_y)						* xjac * theta * tstep 
  
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
  rhs(5) =  v * BigR * particle_source						* xjac * tstep &
	  + v * BigR**2 * ( r0_x * u0_y - r0_y * u0_x)				* xjac * tstep &
	  + v * 2.d0 * BigR * r0 * u0_y						* xjac * tstep &
	  - (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star * Bgrad_rho		* xjac * tstep &
	  - D_prof * BigR  * (v_x*r0_x + v_y*r0_y)				* xjac * tstep &
	  - v * F0 / BigR * Vpar0 * r0_p					* xjac * tstep &
	  - v * Vpar0 * (r0_x * ps0_y - r0_y * ps0_x)				* xjac * tstep &
	  - v * F0 / BigR * r0 * vpar0_p					* xjac * tstep &
	  - v * r0 * (vpar0_x * ps0_y - vpar0_y * ps0_x)			* xjac * tstep &
	  + tau_IC * v * 2.d0 * pi0_y * BigR					* xjac * tstep &
	  + zeta * v * BigR							* xjac *delta_g(5)& 
	  - D_perp_numm *								       &
	    (	 (	   v_ss  * (x_t**2 + y_t**2)					       &
		  +	   v_tt  * (x_s**2 + y_s**2)					       &
		  - 2.d0 * v_st  * (x_s*x_t + y_s*y_t) )				       &
	       * (	   r0_ss * (x_t**2 + y_t**2)					       &
		  +	   r0_tt * (x_s**2 + y_s**2)					       &
		  - 2.d0 * r0_st * (x_s*x_t + y_s*y_t) ) )				       &
	    / xjac**3 *tstep								 

  rhs_k(5) = - (D_par-D_prof) * BigR / BB2 * Bgrad_rho_k_star * Bgrad_rho	* xjac * tstep &
	     - D_prof * BigR  * ( v_p*r0_p * eps_cyl**2 /BigR**2 )		* xjac * tstep 
  
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
                 - tau_IC * v * 2.d0 * (rho_y*Ti0 + rho*Ti0_y) * BigR				* xjac * theta * tstep &
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

  amat(5,6)    = - tau_IC * v * 2.d0 * (Ti_y*r0 + Ti*r0_y) * BigR				* xjac * theta * tstep 
  
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
  Bgrad_Ti_star     = ( v_x  * ps0_y - v_y  * ps0_x ) / BigR
  Bgrad_Ti_k_star   = ( F0 / BigR * v_p  	    ) / BigR
	      
  ! --- The RHS term	      
  rhs(6) =   v * BigR * heat_source_i						* xjac * tstep &
	   + v * r0 * BigR**2 * ( Ti0_x * u0_y - Ti0_y * u0_x)			* xjac * tstep &
	   - v * r0 * F0 / BigR * Vpar0 * Ti0_p					* xjac * tstep &
	   - v * r0 * Vpar0 * (Ti0_x * ps0_y - Ti0_y * ps0_x)			* xjac * tstep &
	   + v * r0 * (GAMMA-1.d0) * Ti0 * u0_y * 2.d0 * BigR			* xjac * tstep &
	   - v * r0 * (GAMMA-1.d0) * Ti0 * (vpar0_x * ps0_y - vpar0_y * ps0_x)	* xjac * tstep &
	   - v * r0 * (GAMMA-1.d0) * Ti0 * F0 / BigR * vpar0_p			* xjac * tstep &
	   -  (Ki_par-Ki_prof) * BigR / BB2 * Bgrad_Ti_star * Bgrad_Ti		* xjac * tstep &
	   -  Ki_prof * BigR * (v_x*Ti0_x + v_y*Ti0_y )				* xjac * tstep &
	   -  Ki_par_num * (v_ps0_x   * ps0_y - v_ps0_y   * ps0_x)			       &
			 * (Ti0_ps0_x * ps0_y - Ti0_ps0_y * ps0_x)		* xjac * tstep &
	   + zeta * v * r0 * BigR						* xjac *delta_g(6)&
	   - Ki_perp_numm *								       &
	     (    (	    v_ss  *  (x_t**2 + y_t**2)					       &
		   +	    v_tt  *  (x_s**2 + y_s**2)					       &
		   - 2.d0 * v_st  *  (x_s*x_t + y_s*y_t) )				       &
		* (	    Ti0_ss * (x_t**2 + y_t**2)					       &
		   +	    Ti0_tt * (x_s**2 + y_s**2)					       &
		   - 2.d0 * Ti0_st * (x_s*x_t + y_s*y_t) ) )				       &
	     / xjac**4								* xjac * tstep 

  rhs_k(6) = - (Ki_par-Ki_prof) * BigR / BB2 * Bgrad_Ti_k_star * Bgrad_Ti	* xjac * tstep &
	     - Ki_par_num * (Ti0_pp  * v_pp)					* xjac * tstep &
	     - Ki_prof * BigR * ( v_p*Ti0_p /BigR**2 )				* xjac * tstep 

  
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
  real*8 :: Bgrad_Ti_star_psi, Bgrad_Ti_psi, Bgrad_Ti_Ti, Bgrad_Ti_Ti_n
  real*8 :: Ti_ps0_x,  Ti_ps0_y
  real*8 :: Ti0_psi_x, Ti0_psi_y
  real*8 :: v_psi_x,   v_psi_y
  
  ! --- Parallel gradient terms       
  Bgrad_Ti_star_psi  = ( v_x   * psi_y - v_y   * psi_x ) / BigR
  Bgrad_Ti_psi	     = ( Ti0_x * psi_y - Ti0_y * psi_x ) / BigR
  Bgrad_Ti_Ti	     = ( Ti_x  * ps0_y - Ti_y  * ps0_x ) / BigR
  Bgrad_Ti_Ti_n	     = ( F0 / BigR * Ti_p   ) / BigR

  Ti_ps0_x  = Ti_xx  * ps0_y - Ti_xy  * ps0_x + Ti_x  * ps0_xy - Ti_y  * ps0_xx
  Ti_ps0_y  = Ti_xy  * ps0_y - Ti_yy  * ps0_x + Ti_x  * ps0_yy - Ti_y  * ps0_xy
  
  Ti0_psi_x = Ti0_xx * psi_y - Ti0_xy * psi_x + Ti0_x * psi_xy - Ti0_y * psi_xx
  Ti0_psi_y = Ti0_xy * psi_y - Ti0_yy * psi_x + Ti0_x * psi_yy - Ti0_y * psi_xy
  
  v_psi_x   = v_xx   * psi_y - v_xy   * psi_x + v_x   * psi_xy - v_y   * psi_xx
  v_psi_y   = v_xy   * psi_y - v_yy   * psi_x + v_x   * psi_yy - v_y   * psi_xy

  ! --- The LHS terms
  amat(6,1)    = + v * r0 * Vpar0 * (Ti0_x * psi_y - Ti0_y * psi_x)					* xjac * theta * tstep &
		 + v * (GAMMA-1.d0) * r0 * Ti0 * (vpar0_x * psi_y - vpar0_y * psi_x)			* xjac * theta * tstep &
		 + Ki_par_num * (v_psi_x  *ps0_y - v_psi_y  *ps0_x + v_ps0_x  *psi_y - v_ps0_y  *psi_x)			       &
			      * (Ti0_ps0_x*ps0_y - Ti0_ps0_y*ps0_x)					* xjac * theta * tstep &
		 + Ki_par_num * (Ti0_psi_x*ps0_y - Ti0_psi_y*ps0_x + Ti0_ps0_x*psi_y - Ti0_ps0_y*psi_x)			       &
			      * (v_ps0_x *ps0_y - v_ps0_y *ps0_x)					* xjac * theta * tstep &
		 - (Ki_par-Ki_prof) * BigR * BB2_psi / BB2**2 * Bgrad_Ti_star     * Bgrad_Ti		* xjac * theta * tstep &
		 + (Ki_par-Ki_prof) * BigR / BB2 	      * Bgrad_Ti_star_psi * Bgrad_Ti		* xjac * theta * tstep &
		 + (Ki_par-Ki_prof) * BigR / BB2 	      * Bgrad_Ti_star     * Bgrad_Ti_psi	* xjac * theta * tstep 

  amat_k(6,1)  = - (Ki_par-Ki_prof) * BigR * BB2_psi / BB2**2 * Bgrad_Ti_k_star   * Bgrad_Ti		* xjac * theta * tstep &
		 + (Ki_par-Ki_prof) * BigR / BB2 	      * Bgrad_Ti_k_star   * Bgrad_Ti_psi	* xjac * theta * tstep 

  amat(6,2)    = - v * r0 * BigR**2 * ( Ti0_x * u_y - Ti0_y * u_x)					* xjac * theta * tstep &
		 - v * 2.d0 * (GAMMA-1.d0) * r0 * BigR * Ti0 * u_y					* xjac * theta * tstep 

  amat(6,5)    = - v * rho * BigR**2 * (Ti0_x * u0_y  - Ti0_y * u0_x)					* xjac * theta * tstep &
		 + v * rho * Vpar0   * (Ti0_x * ps0_y - Ti0_y * ps0_x)					* xjac * theta * tstep &
		 + v * rho * Vpar0   * F0/BigR * Ti0_p							* xjac * theta * tstep &
		 - v * (GAMMA-1.d0) * rho * Ti0 * u0_y * 2.d0 * BigR					* xjac * theta * tstep &
		 + v * (GAMMA-1.d0) * rho * Ti0 * F0/BigR * Vpar0_p					* xjac * theta * tstep &
		 + v * (GAMMA-1.d0) * rho * Ti0 * (Vpar0_x * ps0_y - Vpar0_y * ps0_x)			* xjac * theta * tstep 

  amat(6,6)    = - v * r0 * BigR**2 * (Ti_x * u0_y  - Ti_y * u0_x)					* xjac * theta * tstep &
		 + v * r0 * Vpar0   * (Ti_x * ps0_y - Ti_y * ps0_x)					* xjac * theta * tstep &
		 - v * r0 * (GAMMA-1.d0) * Ti * BigR * u0_y * 2.d0					* xjac * theta * tstep &
		 + v * r0 * (GAMMA-1.d0) * Ti * F0/BigR * Vpar0_p					* xjac * theta * tstep &
		 + v * r0 * (GAMMA-1.d0) * Ti * (Vpar0_x * ps0_y - Vpar0_y * ps0_x)			* xjac * theta * tstep &
		 + Ki_par_num * (v_ps0_x  * ps0_y - v_ps0_y  * ps0_x) 							       &
			      * (Ti_ps0_x * ps0_y - Ti_ps0_y * ps0_x)					* xjac * theta * tstep &
		 + (Ki_par-Ki_prof) * BigR / BB2 * Bgrad_Ti_star * Bgrad_Ti_Ti				* xjac * theta * tstep &
		 + dKi_par * Ti     * BigR / BB2 * Bgrad_Ti_star * Bgrad_Ti				* xjac * theta * tstep &
		 + Ki_prof * BigR * (v_x*Ti_x + v_y*Ti_y )						* xjac * theta * tstep & 
		 + Ki_perp_numm  *											       &
		   (	(	  v_ss  * (x_t**2 + y_t**2)								       &
			 +	  v_tt  * (x_s**2 + y_s**2)								       &
			 - 2.d0 * v_st  * (x_s*x_t + y_s*y_t) )  							       &
		      * (	  Ti_ss * (x_t**2 + y_t**2)								       &
			 +	  Ti_tt * (x_s**2 + y_s**2)								       &
			 - 2.d0 * Ti_st * (x_s*x_t + y_s*y_t) ) )     / xjac**4  			* xjac * theta * tstep 
  
  if (r0 .lt. rho_1) then
    amat(6,6)  = amat(6,6) + v * rho_1 * Ti * BigR							* xjac * (1.d0 + zeta)
  else
    amat(6,6)  = amat(6,6) + v * r0    * Ti * BigR							* xjac * (1.d0 + zeta)
  endif
  
  amat_k(6,6)  = + (Ki_par-Ki_prof) * BigR / BB2 * Bgrad_Ti_k_star * Bgrad_Ti_Ti			* xjac * theta * tstep &
		 + dKi_par * Ti     * BigR / BB2 * Bgrad_Ti_k_star * Bgrad_Ti				* xjac * theta * tstep 
	      
  amat_n(6,6)  = + (Ki_par-Ki_prof) * BigR / BB2 * Bgrad_Ti_star   * Bgrad_Ti_Ti_n  			* xjac * theta * tstep &
		 + v * r0 * Vpar0  * F0/BigR * Ti_p							* xjac * theta * tstep 

  amat_kn(6,6) = + (Ki_par-Ki_prof) * BigR / BB2 * Bgrad_Ti_k_star * Bgrad_Ti_Ti_n  			* xjac * theta * tstep &
		 +  Ki_par_num * (Ti_pp  * v_pp) 							* xjac * theta * tstep &
		 +  Ki_prof	    * BigR	* (v_p*Ti_p /BigR**2 )					* xjac * theta * tstep 

  amat(6,7)    = + v * r0 * F0/BigR * Vpar * Ti0_p							* xjac * theta * tstep &
		 + v * r0 * Vpar * (Ti0_x * ps0_y - Ti0_y * ps0_x)					* xjac * theta * tstep &
		 + v * (GAMMA-1.d0) * r0 * Ti0 * (Vpar_x * ps0_y - Vpar_y * ps0_x)			* xjac * theta * tstep 
     
  amat_n(6,7)  = v * (GAMMA-1.d0) * r0 * Ti0 * F0/BigR * Vpar_p  					* xjac * theta * tstep       
		   

  
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
  rhs(7) = 											       &
           ! --- Time derivative terms (including dB/dt)
	   + zeta * v * r0 * F0**2 / BigR						* xjac *delta_g(7)&  
           + zeta * v * r0 * vpar0 * (ps0_x * delta_ps_x + ps0_y * delta_ps_y) / BigR 	* xjac 	       &
           ! --- Convection terms
	   - 0.5d0 * r0 * vpar0**2 * BB2 * (ps0_x * v_y  - ps0_y * v_x)			* xjac * tstep &
           + 0.5d0 * r0 * vpar0**2 * BB2 * F0 / BigR * v_p				* xjac * tstep &
           - 0.5d0 * v  * vpar0**2 * BB2 * (ps0_x * r0_y - ps0_y * r0_x)		* xjac * tstep &
           + 0.5d0 * v  * vpar0**2 * BB2 * F0 / BigR * r0_p				* xjac * tstep &
           ! --- Parallel pressure gradient
           - v * F0 / BigR * P0_p							* xjac * tstep &
	   - v * (P0_x * ps0_y - P0_y * ps0_x)						* xjac * tstep &
           ! --- Viscosity and toroidal source
           + visco_par * (v_x * Vt0_x   + v_y * Vt0_y)   * BigR 			* xjac * tstep &
	   - visco_par * (v_x * vpar0_x + v_y * vpar0_y) * BigR				* xjac * tstep &
	   - visco_par_numm  *  								       &
	     (    (	    v_ss     * (x_t**2 + y_t**2)					       &
		   +	    v_tt     * (x_s**2 + y_s**2)					       &
		   - 2.d0 * v_st     * (x_s*x_t + y_s*y_t) )					       &
		* (	    vpar0_ss * (x_t**2 + y_t**2)					       &
		   +	    vpar0_tt * (x_s**2 + y_s**2)					       &
		   - 2.d0 * vpar0_st * (x_s*x_t + y_s*y_t) ) )  				       &
	     / xjac**4									* xjac * tstep &
           ! --- Diamagnetic term (RHS of gyroviscous cancellation)
           !+ v * F0 / BigR * Pi0_p*w0							* xjac * tstep &
           !+ v * F0 / BigR * Pi0  *w0_p   						* xjac * tstep &
	   !+ v * (Pi0_x*w0   * ps0_y - Pi0_y*w0   * ps0_x)				* xjac * tstep &
	   !+ v * (Pi0  *w0_x * ps0_y - Pi0  *w0_y * ps0_x)				* xjac * tstep &
           ! --- Neoclassic term
           + v * amu_neo_prof * BB2 / (Btheta2+epsil) * BigR					       &
	       * (  r0                         * (ps0_x*u0_x  + ps0_y*u0_y)			       &
                  + tau_IC                     * (ps0_x*Pi0_x + ps0_y*Pi0_y)			       &
                  + aki_neo_prof * tau_IC * r0 * (ps0_x*Ti0_x + ps0_y*Ti0_y)			       &
                  - r0 * Vpar0 * Btheta2					)	* xjac * tstep 
  
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
  
  ! --- Internal variables
  real*8 :: Btheta2_psi
  
  Btheta2_psi = 2.d0 * (psi_x * ps0_x + psi_y * ps0_y ) / BigR**2
  
  ! -----------------------------
  ! --- The LHS terms (main part)
  amat(7,1)   = + v * (P0_x * psi_y - P0_y * psi_x)								* xjac * theta * tstep &
                + 0.5d0 * r0 * vpar0**2 * BB2     * (psi_x * v_y  - psi_y * v_x)				* xjac * theta * tstep &
                + 0.5d0 * v  * vpar0**2 * BB2     * (psi_x * r0_y - psi_y * r0_x)				* xjac * theta * tstep &
                + 0.5d0 * r0 * vpar0**2 * BB2_psi * (ps0_x * v_y  - ps0_y * v_x)				* xjac * theta * tstep &
                + 0.5d0 * v  * vpar0**2 * BB2_psi * (ps0_x * r0_y - ps0_y * r0_x)				* xjac * theta * tstep &
                + v * r0 * vpar0 / BigR * (ps0_x * psi_x + ps0_y * psi_y)					* xjac * (1.d0 + zeta) 
           
  amat(7,5)   = + v * (rho_x * (Ti0 + Te0)     * ps0_y - rho_y * (Ti0 + Te0)     * ps0_x)			* xjac * theta * tstep &
		+ v * (rho   * (Ti0_x + Te0_x) * ps0_y - rho   * (Ti0_y + Te0_y) * ps0_x)			* xjac * theta * tstep &
		+ v * F0 / BigR * rho * (Ti0_p + Te0_p)								* xjac * theta * tstep & 
		+ 0.5d0 * rho * vpar0**2 * BB2 * (ps0_x * v_y   - ps0_y * v_x)					* xjac * theta * tstep &
		- 0.5d0 * rho * vpar0**2 * BB2 * F0 / BigR * v_p						* xjac * theta * tstep &
                + 0.5d0 * v   * vpar0**2 * BB2 * (ps0_x * rho_y - ps0_y * rho_x)				* xjac * theta * tstep &
                - 0.5d0 * v   * vpar0**2 * BB2 * F0 / BigR * rho_p						* xjac * theta * tstep 

  amat_n(7,5) = v * F0 / BigR * rho_p * (Ti0 + Te0)								* xjac * theta * tstep  

  amat(7,6)   = + v * (Ti_x * r0   * ps0_y - Ti_y * r0   * ps0_x)						* xjac * theta * tstep &
		+ v * (Ti   * r0_x * ps0_y - Ti	  * r0_y * ps0_x)						* xjac * theta * tstep &
		+ v * F0 / BigR * Ti * r0_p									* xjac * theta * tstep 
  
  amat_n(7,6) = v * F0 / BigR * Ti_p * r0									* xjac * theta * tstep 

  amat(7,7)   = + r0 * vpar0 * vpar * BB2 * (ps0_x * v_y  - ps0_y * v_x)					* xjac * theta * tstep &
                - r0 * vpar0 * vpar * BB2 * F0 / BigR * v_p							* xjac * theta * tstep &
                + v  * vpar0 * vpar * BB2 * (ps0_x * r0_y - ps0_y * r0_x)					* xjac * theta * tstep &
                - v  * vpar0 * vpar * BB2 * F0 / BigR * r0_p							* xjac * theta * tstep &
                + visco_par * (v_x * Vpar_x + v_y * Vpar_y) * BigR						* xjac * theta * tstep &
		+ visco_par_numm  *											  	       &
		  (    (	v_ss	* (x_t**2 + y_t**2)								  	       &
			+	v_tt	* (x_s**2 + y_s**2)								  	       &
			-2.d0 * v_st	* (x_s*x_t + y_s*y_t) ) 							  	       &
		     * (	Vpar_ss * (x_t**2 + y_t**2)								  	       &
			+	Vpar_tt * (x_s**2 + y_s**2)								  	       &
			-2.d0 * Vpar_st * (x_s*x_t + y_s*y_t) ) ) / xjac**4					* xjac * theta * tstep
  
  if (r0 .lt. rho_1) then
    amat(7,7) = amat(7,7) + v * Vpar * rho_1 * F0**2 / BigR							* xjac * (1.d0 + zeta)
  else
    amat(7,7) = amat(7,7) + v * Vpar * r0    * F0**2 / BigR							* xjac * (1.d0 + zeta)
  endif

  amat(7,8)   = + v * (Te_x * r0   * ps0_y - Te_y * r0   * ps0_x)						* xjac * theta * tstep &
		+ v * (Te   * r0_x * ps0_y - Te	  * r0_y * ps0_x)						* xjac * theta * tstep &
		+ v * F0 / BigR * Te * r0_p									* xjac * theta * tstep 
  
  amat_n(7,8) = v * F0 / BigR * Te_p * r0									* xjac * theta * tstep

  
  ! ---------------------------------------------------
  ! --- The LHS terms (diamagnetic and neoclassic part)
  amat(7,1)   = amat(7,1)													       &
                !- v * (Pi0_x*w0   * psi_y - Pi0_y*w0   * psi_x)   						* xjac * theta * tstep &
	        !- v * (Pi0  *w0_x * psi_y - Pi0  *w0_y * psi_x)   						* xjac * theta * tstep &
		- v * amu_neo_prof * BB2 / (Btheta2+epsil) * BigR								       &
		                   * (  r0                         * (psi_x*u0_x  + psi_y*u0_y)					       &
                                      + tau_IC                     * (psi_x*Pi0_x + psi_y*Pi0_y)				       &
                                      + aki_neo_prof * tau_IC * r0 * (psi_x*Ti0_x + psi_y*Ti0_y)  )		* xjac * theta * tstep &
                + v * amu_neo_prof * Btheta2_psi * BB2 / Btheta2**2.d0 * BigR							       &
                                   * (  r0                         * (ps0_x*u0_x  + ps0_y*u0_y)					       &
                                      + tau_IC                     * (ps0_x*Pi0_x + ps0_y*Pi0_y)				       &
                                      + aki_neo_prof * tau_IC * r0 * (ps0_x*Ti0_x + ps0_y*Ti0_y)  )		* xjac * theta * tstep
           
  amat(7,2)   = amat(7,2)													       &
                - v * amu_neo_prof * BB2 / (Btheta2+epsil) * r0 * (ps0_x*u_x + ps0_y*u_y) * BigR		* xjac * theta * tstep 
  
  amat(7,4)   = amat(7,4)!													       &
                !- v * F0 / BigR * Pi0_p*w									* xjac * theta * tstep &
	        !- v * (Pi0_x*w   * ps0_y - Pi0_y*w   * ps0_x)							* xjac * theta * tstep &
	        !- v * (Pi0  *w_x * ps0_y - Pi0  *w_y * ps0_x)							* xjac * theta * tstep 
  
  amat_n(7,4) = amat_n(7,4)!													       &
                !- v * F0 / BigR * Pi0*w_p									* xjac * theta * tstep 

  amat(7,5)   = amat(7,5)													       &
                !- v * F0 / BigR * rho*Ti0_p*w0									* xjac * theta * tstep &
                !- v * F0 / BigR * rho*Ti0*w0_p									* xjac * theta * tstep &
	        !- v * (rho_x*Ti0  *w0   * ps0_y - rho_y*Ti0  *w0   * ps0_x)					* xjac * theta * tstep &
	        !- v * (rho  *Ti0_x*w0   * ps0_y - rho  *Ti0_y*w0   * ps0_x)					* xjac * theta * tstep &
	        !- v * (rho  *Ti0  *w0_x * ps0_y - rho  *Ti0  *w0_y * ps0_x)					* xjac * theta * tstep &
                - v * amu_neo_prof * BB2 / (Btheta2+epsil) * BigR								       &
                                   * (  rho * (ps0_x*u0_x + ps0_y*u0_y)								       &
                                      + tau_IC                      * (ps0_x*rho_x*Ti0   + ps0_y*rho_y*Ti0  ) 			       &
                                      + tau_IC                      * (ps0_x*rho  *Ti0_x + ps0_y*rho  *Ti0_y) 			       &
                                      + aki_neo_prof * tau_IC * rho * (ps0_x*Ti0_x       + ps0_y*Ti0_y)				       &
                                      - rho * Vpar0 * Btheta2						      )	* xjac * tstep * theta

  amat_n(7,5) = amat_n(7,5)!													       & 
                !- v * F0 / BigR * rho_p*Ti0*w0									* xjac * theta * tstep 

  amat(7,6)   = amat(7,6)													       &
                !- v * F0 / BigR * r0_p*Ti*w0									* xjac * theta * tstep &
                !- v * F0 / BigR * r0*Ti*w0_p									* xjac * theta * tstep &
	        !- v * (r0_x*Ti  *w0   * ps0_y - r0_y*Ti  *w0   * ps0_x)   					* xjac * theta * tstep &
	        !- v * (r0  *Ti_x*w0   * ps0_y - r0  *Ti_y*w0   * ps0_x)   					* xjac * theta * tstep &
	        !- v * (r0  *Ti  *w0_x * ps0_y - r0  *Ti  *w0_y * ps0_x)   					* xjac * theta * tstep &
                - v * amu_neo_prof * BB2 / (Btheta2+epsil) * BigR								       &
                                   * (  tau_IC                     * (ps0_x*r0_x*Ti   + ps0_y*r0_y*Ti  )			       &
                                      + tau_IC                     * (ps0_x*r0  *Ti_x + ps0_y*r0  *Ti_y)			       &
                                      + aki_neo_prof * tau_IC * r0 * (ps0_x*Ti_x      + ps0_y*Ti_y     )   )	* xjac * tstep * theta
  
  amat_n(7,6) = amat_n(7,6)!													       &
                !- v * F0 / BigR * r0*Ti_p*w0									* xjac * theta * tstep 

  amat(7,7)   = amat(7,7)													       &
                + v * amu_neo_prof * BB2 * r0 * vpar * BigR 							* xjac * tstep * theta 

  
  return

end subroutine ELM_main_lhs_7







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
subroutine ELM_main_rhs_8(rhs,rhs_k)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_rhs_8

  ! --- Modules
  use phys_module
  use equation_variables
  
  implicit none
  
  ! --- Routine variables
  real*8 :: rhs(n_var),rhs_k(n_var)
  
  ! --- Parallel gradient terms       
  v_ps0_x  = v_xx  * ps0_y - v_xy  * ps0_x + v_x  * ps0_xy - v_y * ps0_xx
  v_ps0_y  = v_xy  * ps0_y - v_yy  * ps0_x + v_x  * ps0_yy - v_y * ps0_xy    
  Bgrad_Te_star     = ( v_x  * ps0_y - v_y  * ps0_x ) / BigR
  Bgrad_Te_k_star   = ( F0 / BigR * v_p  	    ) / BigR
	      
  ! --- The RHS term	      
  rhs(8) =   v * BigR * heat_source_e						* xjac * tstep &
	   + v * r0 * BigR**2 * ( Te0_x * u0_y - Te0_y * u0_x)			* xjac * tstep &
	   - v * r0 * F0 / BigR * Vpar0 * Te0_p					* xjac * tstep &
	   - v * r0 * Vpar0 * (Te0_x * ps0_y - Te0_y * ps0_x)			* xjac * tstep &
	   + v * r0 * (GAMMA-1.d0) * Te0 * u0_y * 2.d0 * BigR			* xjac * tstep &
	   - v * r0 * (GAMMA-1.d0) * Te0 * (vpar0_x * ps0_y - vpar0_y * ps0_x)	* xjac * tstep &
	   - v * r0 * (GAMMA-1.d0) * Te0 * F0 / BigR * vpar0_p			* xjac * tstep &
	   -  (Ke_par-Ke_prof) * BigR / BB2 * Bgrad_Te_star * Bgrad_Te		* xjac * tstep &
	   -  Ke_prof * BigR * (v_x*Te0_x + v_y*Te0_y )				* xjac * tstep &
	   -  Ke_par_num * (v_ps0_x   * ps0_y - v_ps0_y   * ps0_x)			       &
			 * (Te0_ps0_x * ps0_y - Te0_ps0_y * ps0_x)		* xjac * tstep &
	   + zeta * v * r0 * BigR						* xjac *delta_g(8)&
	   - Ke_perp_numm *								       &
	     (    (	    v_ss  *  (x_t**2 + y_t**2)					       &
		   +	    v_tt  *  (x_s**2 + y_s**2)					       &
		   - 2.d0 * v_st  *  (x_s*x_t + y_s*y_t) )				       &
		* (	    Te0_ss * (x_t**2 + y_t**2)					       &
		   +	    Te0_tt * (x_s**2 + y_s**2)					       &
		   - 2.d0 * Te0_st * (x_s*x_t + y_s*y_t) ) )				       &
	     / xjac**4								* xjac * tstep 

  rhs_k(8) = - (Ke_par-Ke_prof) * BigR / BB2 * Bgrad_Te_k_star * Bgrad_Te	* xjac * tstep &
	     - Ke_par_num * (Te0_pp  * v_pp)					* xjac * tstep &
	     - Ke_prof * BigR * ( v_p*Te0_p /BigR**2 )				* xjac * tstep 

  
  return

end subroutine ELM_main_rhs_8


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! LHS !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ELM_main_lhs_8(amat, amat_k, amat_n, amat_kn)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_lhs_8

  ! --- Modules
  use phys_module
  use equation_variables
  
  implicit none
  
  ! --- Routine variables
  real*8 :: amat(n_var,n_var), amat_k(n_var,n_var), amat_n(n_var,n_var), amat_kn(n_var,n_var)
  
  ! --- Internal variables
  real*8 :: Bgrad_Te_star_psi, Bgrad_Te_psi, Bgrad_Te_Te, Bgrad_Te_Te_n
  real*8 :: Te_ps0_x,  Te_ps0_y
  real*8 :: Te0_psi_x, Te0_psi_y
  real*8 :: v_psi_x,   v_psi_y
  
  ! --- Parallel gradient terms       
  Bgrad_Te_star_psi  = ( v_x   * psi_y - v_y   * psi_x ) / BigR
  Bgrad_Te_psi	     = ( Te0_x * psi_y - Te0_y * psi_x ) / BigR
  Bgrad_Te_Te	     = ( Te_x  * ps0_y - Te_y  * ps0_x ) / BigR
  Bgrad_Te_Te_n	     = ( F0 / BigR * Te_p   ) / BigR

  Te_ps0_x  = Te_xx  * ps0_y - Te_xy  * ps0_x + Te_x  * ps0_xy - Te_y  * ps0_xx
  Te_ps0_y  = Te_xy  * ps0_y - Te_yy  * ps0_x + Te_x  * ps0_yy - Te_y  * ps0_xy
  
  Te0_psi_x = Te0_xx * psi_y - Te0_xy * psi_x + Te0_x * psi_xy - Te0_y * psi_xx
  Te0_psi_y = Te0_xy * psi_y - Te0_yy * psi_x + Te0_x * psi_yy - Te0_y * psi_xy
  
  v_psi_x   = v_xx   * psi_y - v_xy   * psi_x + v_x   * psi_xy - v_y   * psi_xx
  v_psi_y   = v_xy   * psi_y - v_yy   * psi_x + v_x   * psi_yy - v_y   * psi_xy

  ! --- The LHS terms
  amat(8,1)    = + v * r0 * Vpar0 * (Te0_x * psi_y - Te0_y * psi_x)					* xjac * theta * tstep &
		 + v * (GAMMA-1.d0) * r0 * Te0 * (vpar0_x * psi_y - vpar0_y * psi_x)			* xjac * theta * tstep &
		 + Ke_par_num * (v_psi_x  *ps0_y - v_psi_y  *ps0_x + v_ps0_x  *psi_y - v_ps0_y  *psi_x)			       &
			      * (Te0_ps0_x*ps0_y - Te0_ps0_y*ps0_x)					* xjac * theta * tstep &
		 + Ke_par_num * (Te0_psi_x*ps0_y - Te0_psi_y*ps0_x + Te0_ps0_x*psi_y - Te0_ps0_y*psi_x)			       &
			      * (v_ps0_x *ps0_y - v_ps0_y *ps0_x)					* xjac * theta * tstep &
		 - (Ke_par-Ke_prof) * BigR * BB2_psi / BB2**2 * Bgrad_Te_star     * Bgrad_Te		* xjac * theta * tstep &
		 + (Ke_par-Ke_prof) * BigR / BB2 	      * Bgrad_Te_star_psi * Bgrad_Te		* xjac * theta * tstep &
		 + (Ke_par-Ke_prof) * BigR / BB2 	      * Bgrad_Te_star     * Bgrad_Te_psi	* xjac * theta * tstep 

  amat_k(8,1)  = - (Ke_par-Ke_prof) * BigR * BB2_psi / BB2**2 * Bgrad_Te_k_star   * Bgrad_Te		* xjac * theta * tstep &
		 + (Ke_par-Ke_prof) * BigR / BB2 	      * Bgrad_Te_k_star   * Bgrad_Te_psi	* xjac * theta * tstep 

  amat(8,2)    = - v * r0 * BigR**2 * ( Te0_x * u_y - Te0_y * u_x)					* xjac * theta * tstep &
		 - v * 2.d0 * (GAMMA-1.d0) * r0 * BigR * Te0 * u_y					* xjac * theta * tstep 

  amat(8,5)    = - v * rho * BigR**2 * (Te0_x * u0_y  - Te0_y * u0_x)					* xjac * theta * tstep &
		 + v * rho * Vpar0   * (Te0_x * ps0_y - Te0_y * ps0_x)					* xjac * theta * tstep &
		 + v * rho * Vpar0   * F0/BigR * Te0_p							* xjac * theta * tstep &
		 - v * (GAMMA-1.d0) * rho * Te0 * u0_y * 2.d0 * BigR					* xjac * theta * tstep &
		 + v * (GAMMA-1.d0) * rho * Te0 * F0/BigR * Vpar0_p					* xjac * theta * tstep &
		 + v * (GAMMA-1.d0) * rho * Te0 * (Vpar0_x * ps0_y - Vpar0_y * ps0_x)			* xjac * theta * tstep 

  amat(8,7)    = + v * r0 * F0/BigR * Vpar * Te0_p							* xjac * theta * tstep &
		 + v * r0 * Vpar * (Te0_x * ps0_y - Te0_y * ps0_x)					* xjac * theta * tstep &
		 + v * (GAMMA-1.d0) * r0 * Te0 * (Vpar_x * ps0_y - Vpar_y * ps0_x)			* xjac * theta * tstep 
     
  amat_n(8,7)  = v * (GAMMA-1.d0) * r0 * Te0 * F0/BigR * Vpar_p  					* xjac * theta * tstep       
		   
  amat(8,8)    = - v * r0 * BigR**2 * (Te_x * u0_y  - Te_y * u0_x)					* xjac * theta * tstep &
		 + v * r0 * Vpar0   * (Te_x * ps0_y - Te_y * ps0_x)					* xjac * theta * tstep &
		 - v * r0 * (GAMMA-1.d0) * Te * BigR * u0_y * 2.d0					* xjac * theta * tstep &
		 + v * r0 * (GAMMA-1.d0) * Te * F0/BigR * Vpar0_p					* xjac * theta * tstep &
		 + v * r0 * (GAMMA-1.d0) * Te * (Vpar0_x * ps0_y - Vpar0_y * ps0_x)			* xjac * theta * tstep &
		 + Ke_par_num * (v_ps0_x  * ps0_y - v_ps0_y  * ps0_x) 							       &
			      * (Te_ps0_x * ps0_y - Te_ps0_y * ps0_x)					* xjac * theta * tstep &
		 + (Ke_par-Ke_prof) * BigR / BB2 * Bgrad_Te_star * Bgrad_Te_Te				* xjac * theta * tstep &
		 + dKe_par * Te     * BigR / BB2 * Bgrad_Te_star * Bgrad_Te				* xjac * theta * tstep &
		 + Ke_prof * BigR * (v_x*Te_x + v_y*Te_y )						* xjac * theta * tstep & 
		 + Ke_perp_numm  *											       &
		   (	(	  v_ss  * (x_t**2 + y_t**2)								       &
			 +	  v_tt  * (x_s**2 + y_s**2)								       &
			 - 2.d0 * v_st  * (x_s*x_t + y_s*y_t) )  							       &
		      * (	  Te_ss * (x_t**2 + y_t**2)								       &
			 +	  Te_tt * (x_s**2 + y_s**2)								       &
			 - 2.d0 * Te_st * (x_s*x_t + y_s*y_t) ) )     / xjac**4  			* xjac * theta * tstep 
  
  if (r0 .lt. rho_1) then
    amat(8,8)  = amat(8,8) + v * rho_1 * Te * BigR							* xjac * (1.d0 + zeta)
  else
    amat(8,8)  = amat(8,8) + v * r0    * Te * BigR							* xjac * (1.d0 + zeta)
  endif
  
  amat_k(8,8)  = + (Ke_par-Ke_prof) * BigR / BB2 * Bgrad_Te_k_star * Bgrad_Te_Te			* xjac * theta * tstep &
		 + dKe_par * Te     * BigR / BB2 * Bgrad_Te_k_star * Bgrad_Te				* xjac * theta * tstep 
	      
  amat_n(8,8)  = + (Ke_par-Ke_prof) * BigR / BB2 * Bgrad_Te_star   * Bgrad_Te_Te_n  			* xjac * theta * tstep &
		 + v * r0 * Vpar0  * F0/BigR * Te_p							* xjac * theta * tstep 

  amat_kn(8,8) = + (Ke_par-Ke_prof) * BigR / BB2 * Bgrad_Te_k_star * Bgrad_Te_Te_n  			* xjac * theta * tstep &
		 +  Ke_par_num * (Te_pp  * v_pp) 							* xjac * theta * tstep &
		 +  Ke_prof	    * BigR	* (v_p*Te_p /BigR**2 )					* xjac * theta * tstep 


  
  return

end subroutine ELM_main_lhs_8







!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------ Equation 9 (w* - Diamagnetic vorticity) -------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! RHS !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ELM_main_rhs_9(rhs,rhs_k)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_rhs_4

  ! --- Modules
  use phys_module
  use equation_variables
  
  implicit none
  
  ! --- Routine variables
  real*8 :: rhs(n_var),rhs_k(n_var)
  
  ! --- The RHS term	      
  rhs(9) = 0.d0
  !rhs(9) = - tau_IC / r0 * ( v_x * Pi0_x   + v_y * Pi0_y  + v*w0)  * BigR * xjac * tstep

  return

end subroutine ELM_main_rhs_9

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! LHS !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ELM_main_lhs_9(amat, amat_k, amat_n, amat_kn)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_lhs_4

  ! --- Modules
  use phys_module
  use equation_variables
  
  implicit none
  
  ! --- Routine variables
  real*8 :: amat(n_var,n_var), amat_k(n_var,n_var), amat_n(n_var,n_var), amat_kn(n_var,n_var)
  
  ! --- The LHS terms
  amat(9,5) = + tau_IC / r0             * (v_x * rho  *Ti0_x + v_y * rho  *Ti0_y) * BigR  		* xjac * tstep &
              + tau_IC / r0             * (v_x * rho_x*Ti0   + v_y * rho_y*Ti0  ) * BigR  		* xjac * tstep &
              - tau_IC * rho / r0**2.d0 * (v_x * Pi0_x       + v_y * Pi0_y      ) * BigR  		* xjac * tstep 

  amat(9,6) = + tau_IC / r0             * (v_x * r0   *Ti_x  + v_y * r0   *Ti_y ) * BigR  		* xjac * tstep &
              + tau_IC / r0             * (v_x * r0_x *Ti    + v_y * r0_y *Ti   ) * BigR  		* xjac * tstep 

  amat(9,9) =  v * Wdia * BigR		      								* xjac * tstep 
		    
  
  return

end subroutine ELM_main_lhs_9








