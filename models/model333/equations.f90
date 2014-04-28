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
  
  ! ----------------------------	      
  ! --- The RHS term (main part)	      
  rhs(1) = 											&
           ! --- Time derivative
	   + zeta * v / BigR							* xjac * delta_g(1)&
           ! --- Resistive term
           + v * eta_T  * (zj0 - current_source)/ BigR				* xjac * tstep  &
           ! --- VxB
	   + v * (ps0_x * u0_y - ps0_y * u0_x)  				* xjac * tstep  &
           ! --- Integration term
	   - v * eps_cyl * F0 / BigR  * u0_p					* xjac * tstep  &
           ! --- Numerical resistivity
	   + eta_numm * (v_x * zj0_x + v_y * zj0_y)				* xjac * tstep  
  
  ! -----------------------------------    
  ! --- The RHS term (diamagnetic part)	      
  rhs(1) = rhs(1)										&
           - v * tau_IC/(r0*BB2) * F0**2/BigR**2 * (ps0_x*p0_y - ps0_y*p0_x)	* xjac * tstep  & 		  
           + v * tau_IC/(r0*BB2) * F0**3/BigR**3 * eps_cyl * p0_p		* xjac * tstep  
  
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
  
  ! -----------------------------	      
  ! --- The LHS terms (main part)
  amat(1,1)   = + v * psi / BigR									* xjac * (1.d0+zeta)	&
		- v * (psi_x * u0_y - psi_y * u0_x)							* xjac * theta * tstep	

  amat(1,2)   = -  v * (ps0_x * u_y - ps0_y * u_x)							* xjac * theta * tstep

  amat_n(1,2) = +  eps_cyl * F0 / BigR * v * u_p							* xjac * theta * tstep

  amat(1,3)   = - eta_numm * (v_x * zj_x + v_y * zj_y)							* xjac * theta * tstep	&
	        - eta_T * v * zj / BigR									* xjac * theta * tstep

  amat(1,6)   = - deta_dT * v * T * (zj0 - current_source) / BigR					* xjac * theta * tstep	

  ! ------------------------------------
  ! --- The LHS terms (diamagnetic part)
  amat(1,1)   = amat(1,1)													&
                - v * tau_IC/(r0*BB2**2) * BB2_psi * F0**2/BigR**2 * (ps0_x*p0_y - ps0_y*p0_x)		* xjac * theta * tstep	&
                + v * tau_IC/(r0*BB2**2) * BB2_psi * F0**3/BigR**3 * eps_cyl * p0_p			* xjac * theta * tstep	&
                + v * tau_IC/(r0*BB2) * F0**2/BigR**2 * (psi_x * p0_y - psi_y * p0_x)			* xjac * theta * tstep
  
  amat(1,5)   = amat(1,5)													&
                + v * tau_IC/(r0*BB2) * F0**2/BigR**2 * T0  * (ps0_x*rho_y - ps0_y*rho_x)		* xjac * theta * tstep	&
	        + v * tau_IC/(r0*BB2) * F0**2/BigR**2 * rho * (ps0_x*T0_y  - ps0_y*T0_x )		* xjac * theta * tstep	&
	        - v * tau_IC/(r0*BB2) * F0**3/BigR**3 * eps_cyl * rho * T0_p				* xjac * theta * tstep	&
		- v * tau_IC * rho /(r0**2 * BB2) * F0**2/BigR**2 * (ps0_x*p0_y - ps0_y*p0_x)		* xjac * theta * tstep	&		    
		+ v * tau_IC * rho /(r0**2 * BB2) * F0**3/BigR**3 * eps_cyl * p0_p			* xjac * theta * tstep 
  
  amat_n(1,5) = amat_n(1,5)													&
                - v * tau_IC/(r0*BB2) * F0**3/BigR**3 * eps_cyl * T0 * rho_p				* xjac * theta * tstep
  
  amat(1,6)   = amat(1,6)													&
                + v * tau_IC/(r0*BB2) * F0**2/BigR**2 * r0 * (ps0_x*T_y  - ps0_y*T_x ) 			* xjac * theta * tstep	&
	        + v * tau_IC/(r0*BB2) * F0**2/BigR**2 * T  * (ps0_x*r0_y - ps0_y*r0_x)			* xjac * theta * tstep	&
		- v * tau_IC/(r0*BB2) * F0**3/BigR**3 * eps_cyl * T * r0_p				* xjac * theta * tstep 

  amat_n(1,6) = amat_n(1,6)													&
                - v * tau_IC/(r0*BB2) * F0**3/BigR**3 * eps_cyl * r0 * T_p				* xjac * theta * tstep 
  
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
  
  ! ----------------------------	      
  ! --- The RHS term (main part)	      
  rhs(2) = 													&
           ! --- Time derivative
	   - zeta * BigR * r0_hat * (v_x * delta_u_x + v_y * delta_u_y)  			* xjac		&		  
           ! --- Convective terms
	   - 0.5d0 * vv2 * (v_x * r0_y_hat - v_y * r0_x_hat)		 			* xjac * tstep	&
	   - r0_hat * BigR**2 * w0 * (v_x * u0_y - v_y * u0_x)  	 			* xjac * tstep	&
           ! --- [psi,j]
	   + v * (ps0_x * zj0_y - ps0_y * zj0_x )			 			* xjac * tstep	&
	   - v * eps_cyl * F0 / BigR * zj0_p				 			* xjac * tstep	&
           ! --- Grad(p)
	   + BigR**2 * (v_x * p0_y     - v_y * p0_x)			 			* xjac * tstep	&
           ! --- Source interaction
	   + BigR**3 * total_rho_source * (v_x * u0_x + v_y * u0_y)				* xjac * tstep	&
           ! --- Viscosity
	   - visco_T * BigR * (v_x * w0_x + v_y * w0_y) 		 			* xjac * tstep	&
           ! --- Numerical iscosity
           - visco_numm * (v_xx  + v_x /BigR + v_yy ) 								&
	                * (w0_xx + w0_x/Bigr + w0_yy)						* xjac * tstep 
  
  ! --------------------------------------------------      
  ! --- The RHS term (diamagnetic and neoclassic part)	      
  rhs(2) = rhs(2)												&
	   ! --- Main diamagnetic terms
	   - tau_IC * v * BigR**4        * (P0_x * w0_y - P0_y * w0_x)				* xjac * tstep	&
	   - tau_IC     * BigR**3 * P0_y * (v_x  * u0_x + v_y  * u0_y)				* xjac * tstep	&
	   - tau_IC * v * BigR**4 * (u0_xy * (P0_xx-P0_yy) - P0_xy * (u0_xx-u0_yy) )		* xjac * tstep	&
	   ! --- Inverse diamagnetic terms (needed when including diamagnetic vorticity directly into W - equation4)
	   + tau_IC * BigR**4 *      (P0_xx + P0_x/BigR + P0_yy) * (v_x * u0_y - v_y * u0_x)	* xjac * tstep	&
	   - tau_IC * BigR**4 / r0 * (r0_x * P0_x + r0_y * P0_y) * (v_x * u0_y - v_y * u0_x)	* xjac * tstep	&
	   ! --- Neoclassic term
           + amu_neo_prof * BB2 / (Btheta2+epsil)**2.d0 * (ps0_x*v_x + ps0_y*v_y) * BigR			&
                    * (  r0                         * (ps0_x*u0_x + ps0_y*u0_y)					&
		       + tau_IC                     * (ps0_x*P0_x + ps0_y*P0_y)					&
                       + aki_neo_prof * tau_IC * r0 * (ps0_x*T0_x + ps0_y*T0_y)					&
                       - r0 * Vpar0 * Btheta2					     )		* xjac * tstep 

  
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
  real*8 :: Btheta2_psi
  
  rho_hat     = BigR**2 * rho
  rho_x_hat   = 2.d0 * BigR * BigR_x  * rho + BigR**2 * rho_x
  rho_y_hat   = BigR**2 * rho_y
  Btheta2_psi = 2.d0 * (psi_x * ps0_x + psi_y * ps0_y ) / BigR**2

  ! -----------------------------	      
  ! --- The LHS terms (main part)
  amat(2,1)   = - v * (psi_x * zj0_y - psi_y * zj0_x )								* xjac * theta * tstep

  amat(2,2)   = + r0_hat * BigR**2 * w0 * (v_x * u_y  - v_y  * u_x)						* xjac * theta * tstep	&
		+ BigR**2 * (u_x*u0_x + u_y*u0_y) * (v_x*r0_y_hat - v_y*r0_x_hat)				* xjac * theta * tstep	&
		- BigR**3 * total_rho_source * (v_x * u_x + v_y * u_y)						* xjac * theta * tstep
  
  if (r0 .lt. rho_1) then
    amat(2,2) = amat(2,2) - BigR**3 * rho_1 * (v_x * u_x + v_y * u_y)						* xjac * (1.d0 + zeta) 
  else
    amat(2,2) = amat(2,2) - BigR**3 * r0    * (v_x * u_x + v_y * u_y)						* xjac * (1.d0 + zeta) 
  endif

  amat(2,3)   = - v * (ps0_x * zj_y  - ps0_y * zj_x)								* xjac * theta * tstep

  amat_n(2,3) = + eps_cyl * F0 / BigR * v * zj_p								* xjac * theta * tstep

  amat(2,4)   = r0_hat * BigR**2 * w  * ( v_x * u0_y - v_y * u0_x)						* xjac * theta * tstep	&
		+ visco_T * BigR * ( v_x * w_x + v_y * w_y)							* xjac * theta * tstep	&
                + visco_numm * (v_xx + v_x/BigR + v_yy) 										&
		             * (w_xx + w_x/BigR + w_yy) 							* xjac * theta * tstep    

  amat(2,5)   = + 0.5d0 * vv2 * (v_x * rho_y_hat - v_y * rho_x_hat)						* xjac * theta * tstep	&
		+ rho_hat * BigR**2 * w0 * (v_x * u0_y - v_y * u0_x)						* xjac * theta * tstep	&
		- BigR**2 * (v_x * rho_y * T0   - v_y * rho_x * T0  )						* xjac * theta * tstep	&
		- BigR**2 * (v_x * rho   * T0_y - v_y * rho   * T0_x)						* xjac * theta * tstep

  amat(2,6)   = - BigR**2 * (v_x * r0_y * T   - v_y * r0_x * T)							* xjac * theta * tstep	&
		- BigR**2 * (v_x * r0	* T_y - v_y * r0   * T_x)						* xjac * theta * tstep	&
		+ dvisco_dT * T * ( v_x * w0_x    + v_y * w0_y    ) * BigR					* xjac * theta * tstep
    
  
  ! ---------------------------------------------------    
  ! --- The LHS terms (diamagnetic and neoclassic part)
  amat(2,1)   = amat(2,1)														&
	        ! --- Neoclassical term
                - amu_neo_prof * BB2 / (Btheta2+epsil)**2 * (psi_x*v_x+psi_y*v_y) * BigR						&
                               * (  r0                         * (ps0_x*u0_x + ps0_y*u0_y)						&
                                  + tau_IC                     * (ps0_x*P0_x + ps0_y*P0_y)						&
                                  + aki_neo_prof * tau_IC * r0 * (ps0_x*T0_x + ps0_y*T0_y)						&
                                  - r0 * Vpar0 * Btheta2)	 						* xjac * theta * tstep	&
                - amu_neo_prof * BB2 / (Btheta2+epsil)**2 * (ps0_x*v_x+ps0_y*v_y) * BigR						&
		               * (  r0                         * (psi_x*u0_x + psi_y*u0_y)						&
                                  + tau_IC                     * (psi_x*P0_x + psi_y*P0_y)						&
                                  + aki_neo_prof * tau_IC * r0 * (psi_x*T0_x + psi_y*T0_y)	)		* xjac * theta * tstep	&
                + amu_neo_prof * BB2 * 2.d0*Btheta2_psi / (Btheta2+epsil)**3 * (ps0_x*v_x+ps0_y*v_y) * BigR				&
                               * (  r0                         * (ps0_x*u0_x + ps0_y*u0_y)						&
                                  + tau_IC                     * (ps0_x*P0_x + ps0_y*P0_y)						&
                                  + aki_neo_prof * tau_IC * r0 * (ps0_x*T0_x + ps0_y*T0_y)	)		* xjac * theta * tstep	&
                - amu_neo_prof * BB2 * Btheta2_psi / (Btheta2+epsil)**2									&
		               * r0 * vpar0 * (ps0_x*v_x + ps0_y*v_y) * BigR					* xjac * tstep * theta

  amat(2,2)   = amat(2,2)														&
	        ! --- Main diamagnetic terms
	        + tau_IC * BigR**3 * P0_y * (v_x* u_x + v_y * u_y)						* xjac * theta * tstep	&
	        + tau_IC * v * BigR**4 * (u_xy * (P0_xx-P0_yy) - P0_xy * (u_xx-u_yy))				* xjac * theta * tstep	&
	        ! --- Inverse diamagnetic terms
	        - tau_IC * BigR**4      * (P0_xx + P0_x/BigR + P0_yy)  * (v_x * u_y - v_y * u_x)		* xjac * theta * tstep	&
	        + tau_IC * BigR**4 / r0 * (r0_x * P0_x + r0_y * P0_y ) * (v_x * u_y - v_y * u_x)		* xjac * theta * tstep	&
	        ! --- Neoclassical term
                - amu_neo_prof * BB2 / (Btheta2+epsil)**2 * (ps0_x*v_x + ps0_y*v_y) * BigR						&
		               * r0 * (ps0_x*u_x + ps0_y*u_y)							* xjac * theta * tstep
  
  amat(2,4)   = amat(2,4)														&
	        ! --- Main diamagnetic terms
                + tau_IC * v * BigR**4 * (P0_x * w_y - P0_y * w_x)              				* xjac * theta * tstep 

  amat(2,5)   = amat(2,5)														&
	        ! --- Main diamagnetic terms
                + tau_IC * v * BigR**4 * T0  * (rho_x  * w0_y - rho_y  * w0_x)  				* xjac * theta * tstep	&
                + tau_IC * v * BigR**4 * rho * (T0_x   * w0_y - T0_y   * w0_x)  				* xjac * theta * tstep	&
		+ tau_IC     * BigR**3 * (T0_y*rho + T0*rho_y) * (v_x*u0_x + v_y*u0_y) 				* xjac * theta * tstep	&
		+ tau_IC * v * BigR**4 * ( u0_xy        * (rho_xx*T0 + 2.d0*rho_x*T0_x + rho*T0_xx  					&
			    	                          -rho_yy*T0 - 2.d0*rho_y*T0_y - rho*T0_yy)	   				&					   
			                 -(u0_xx-u0_yy) * (rho_xy*T0 + rho_x*T0_y + rho_y*T0_x + rho*T0_xy))	* xjac * theta * tstep	&
	        ! --- Inverse diamagnetic terms
	        - tau_IC * BigR**4 * ( rho_xx*T0 + rho*T0_xx + 2.d0*rho_x*T0_x								&
		                     + rho_x*T0/BigR + rho*T0_x/BigR									&
		                     + rho_yy*T0 + rho*T0_yy + 2.d0*rho_y*T0_y )  * (v_x * u0_y - v_y * u0_x)	* xjac * theta * tstep	&
	        + tau_IC * BigR**4 / r0 * (rho_x * P0_x + rho_y * P0_y)           * (v_x * u0_y - v_y * u0_x)	* xjac * theta * tstep	&
	        + tau_IC * BigR**4 / r0 * ( r0_x * (rho_x*T0 + rho*T0_x) 								&
		                          + r0_y * (rho_y*T0 + rho*T0_y) )        * (v_x * u0_y - v_y * u0_x)	* xjac * theta * tstep	&
	        - tau_IC * BigR**4 * rho/r0**2.d0 * (r0_x * P0_x + r0_y * P0_y)   * (v_x * u0_y - v_y * u0_x)	* xjac * theta * tstep	&
	        ! --- Neoclassical term
                - amu_neo_prof * BB2 / (Btheta2+epsil)**2 * (ps0_x*v_x + ps0_y*v_y) * BigR						&
                               * (  rho                         * (ps0_x*u0_x       + ps0_y*u0_y      )					&
                                  + tau_IC                      * (ps0_x*rho_x*T0   + ps0_y*rho_y*T0  )					&
                                  + tau_IC                      * (ps0_x*rho  *T0_x + ps0_y*rho  *T0_y)					&
                                  + aki_neo_prof * tau_IC * rho * (ps0_x*T0_x                 + ps0_y*T0_y)				&
                                  -rho * Vpar0 * Btheta2						      )	* xjac * tstep * theta

  amat(2,6)   = amat(2,6)														&
	        ! --- Main diamagnetic terms
                + tau_IC * v * BigR**4 * r0  * (T_x*w0_y - T_y*w0_x)  						* xjac * theta * tstep	&
                + tau_IC * v * BigR**4 * T  * (r0_x*w0_y - r0_y*w0_x)  						* xjac * theta * tstep	&
		+ tau_IC     * BigR**3 * (r0_y*T + r0*T_y) * (v_x*u0_x + v_y*u0_y) 				* xjac * theta * tstep	&
		+ tau_IC * v * BigR**4 * ( u0_xy        * (r0_xx*T + 2.d0*r0_x*T_x + r0*T_xx      					&
				                          -r0_yy*T - 2.d0*r0_y*T_y - r0*T_yy)						&					 
			                 -(u0_xx-u0_yy) * (r0_xy*T + r0_x*T_y + r0_y*T_x + r0*T_xy)	)	* xjac * theta * tstep	&
	        ! --- Inverse diamagnetic terms
	        - tau_IC * BigR**4 * ( r0_xx*T + r0*T_xx + 2.d0*r0_x*T_x       								&
		                     + r0_x*T/BigR + r0*T_x/BigR									&
		                     + r0_yy*T + r0*T_yy + 2.d0*r0_y*T_y ) * (v_x * u0_y - v_y * u0_x)		* xjac * theta * tstep	&
	        + tau_IC * BigR**4 / r0 * ( r0_x * (r0_x*T + r0*T_x)									&
		                          + r0_y * (r0_y*T + r0*T_y) )     * (v_x * u0_y - v_y * u0_x)		* xjac * theta * tstep	&
	        ! --- Neoclassical term
                - amu_neo_prof * BB2 / (Btheta2+epsil)**2 * (ps0_x*v_x + ps0_y*v_y) * BigR						&
                               * (  tau_IC                     * (ps0_x*r0_x*T   + ps0_y*r0_y*T  )					&
                                  + tau_IC                     * (ps0_x*r0  *T_x + ps0_y*r0  *T_y)					&
                                  + aki_neo_prof * tau_IC * r0 * (ps0_x*T_x+ps0_y*T_y)		   )		* xjac * tstep * theta
  
  amat(2,7)   = amat(2,7)														&
	        ! --- Neoclassical term
                + amu_neo_prof * BB2 / (Btheta2+epsil) * r0 * vpar * (ps0_x*v_x + ps0_y*v_y) * BigR		* xjac * tstep * theta 
  
  
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
  rhs(3) = 0.d0
  !rhs(3) = - ( v_x * ps0_x  + v_y * ps0_y + v*zj0) / BigR * xjac * tstep

  
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
  
  ! ----------------	      
  ! --- The RHS term	      
  !rhs(4) = 0.d0
  
  ! ----------------------------	      
  ! --- The RHS term (main part)	      
  rhs(4) = - ( v_x * u0_x   + v_y * u0_y  + v*w0)  			* BigR * xjac * tstep	&
  
  ! -----------------------------------    
  ! --- The RHS term (diamagnetic part)	      
  rhs(4) = rhs(4)										&
  	   - tau_IC / r0 * ( v_x  * P0_x + v_y  * P0_y)         	* BigR * xjac * tstep	

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
  
  ! -----------------------------	      
  ! --- The LHS terms (main part)
  amat(4,2) = (v_x * u_x + v_y * u_y) * BigR  									* xjac * tstep

  amat(4,4) =  v * w * BigR		      									* xjac * tstep 
		    
  ! ------------------------------------
  ! --- The LHS terms (diamagnetic part)
  amat(4,5) = amat(4,5)									        				&
              + tau_IC / r0                        * (v_x   * rho  *T0_x + v_y   * rho  *T0_y) * BigR  		* xjac * tstep	&
              + tau_IC / r0                        * (v_x   * rho_x*T0   + v_y   * rho_y*T0  ) * BigR  		* xjac * tstep	&
              - tau_IC * rho / r0**2.d0            * (v_x   * P0_x       + v_y   * P0_y      ) * BigR  		* xjac * tstep	

  amat(4,6) = amat(4,6)									        	    			&
              + tau_IC / r0                        * (v_x   * r0   *T_x  + v_y   * r0   *T_y ) * BigR    	* xjac * tstep	&
              + tau_IC / r0                        * (v_x   * r0_x *T    + v_y   * r0_y *T   ) * BigR    	* xjac * tstep	

  
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
  Bgrad_rho        = ( r0_x * ps0_y - r0_y * ps0_x &
                     + F0 / BigR * r0_p            ) / BigR
  Bgrad_rho_star   = ( v_x  * ps0_y - v_y  * ps0_x ) / BigR
  Bgrad_rho_k_star = ( F0 / BigR * v_p  	   ) / BigR
	      
  ! ----------------------------	      
  ! --- The RHS term (main part)	      
  rhs(5)   =   											&
             ! --- Time derivative
	     + zeta * v * BigR							* xjac *delta_g(5)& 
             ! --- Div(rho.V)
	     + v * BigR**2 * ( r0_x * u0_y - r0_y * u0_x)			* xjac * tstep	&
	     + v * 2.d0 * BigR * r0 * u0_y					* xjac * tstep	&
	     - v * F0 / BigR * Vpar0 * r0_p					* xjac * tstep	&
	     - v * F0 / BigR * r0    * vpar0_p					* xjac * tstep	&
	     - v * Vpar0 * (r0_x    * ps0_y - r0_y    * ps0_x)			* xjac * tstep	&
	     - v * r0    * (vpar0_x * ps0_y - vpar0_y * ps0_x)			* xjac * tstep	&
             ! --- Source
             + v * BigR * total_rho_source					* xjac * tstep	&
             ! --- Diffusivity
	     - (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star * Bgrad_rho		* xjac * tstep	&
	     - D_prof * BigR  * (v_x*r0_x + v_y*r0_y)				* xjac * tstep	&
             ! --- Numerical diffusivity
             - D_perp_numm * (v_xx  + v_x /Bigr + v_yy )					&
	                   * (r0_xx + r0_x/Bigr + r0_yy) * BigR 		* xjac * tstep 

  rhs_k(5) = - (D_par-D_prof) * BigR / BB2 * Bgrad_rho_k_star * Bgrad_rho	* xjac * tstep	&
	     - D_prof * BigR  * ( v_p*r0_p * eps_cyl**2 /BigR**2 )		* xjac * tstep 
  
  ! -----------------------------------    
  ! --- The RHS term (diamagnetic part)	      
  rhs(5)   = rhs(5)									 	&
	     + tau_IC * v * 2.d0 * p0_y * BigR					* xjac * tstep	
  
  
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

  ! -----------------------------	      
  ! --- The LHS terms (main part)
  amat(5,1)    = - (D_par-D_prof) * BigR * BB2_psi/ BB2**2 * Bgrad_rho_star	* Bgrad_rho	* xjac * theta * tstep	&
		 + (D_par-D_prof) * BigR / BB2  	   * Bgrad_rho_star_psi * Bgrad_rho	* xjac * theta * tstep	&
		 + (D_par-D_prof) * BigR / BB2  	   * Bgrad_rho_star	* Bgrad_rho_psi * xjac * theta * tstep	&
		 + v * Vpar0 * (r0_x * psi_y - r0_y * psi_x)					* xjac * theta * tstep	&
		 + v * r0 * (vpar0_x * psi_y - vpar0_y * psi_x) 				* xjac * theta * tstep 

  amat_k(5,1)  = - (D_par-D_prof) * BigR * BB2_psi/ BB2**2 * Bgrad_rho_k_star * Bgrad_rho	* xjac * theta * tstep	&
		 + (D_par-D_prof) * BigR / BB2  	   * Bgrad_rho_k_star * Bgrad_rho_psi	* xjac * theta * tstep 

  amat(5,2)    = - v * BigR**2 * ( r0_x * u_y - r0_y * u_x)					* xjac * theta * tstep	&
		 - v * 2.d0 * BigR * r0 * u_y							* xjac * theta * tstep 

  amat(5,5)    = + v * rho * BigR								* xjac * (1.d0 + zeta)	&
		 - v * BigR**2 * ( rho_x * u0_y - rho_y * u0_x) 				* xjac * theta * tstep	&
		 - v * 2.d0 * BigR * rho * u0_y 						* xjac * theta * tstep	&
		 + (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star * Bgrad_rho_rho 		* xjac * theta * tstep	&
		 + D_prof * BigR  * (v_x*rho_x + v_y*rho_y )					* xjac * theta * tstep	&
		 + v * Vpar0 * (rho_x * ps0_y - rho_y * ps0_x)  				* xjac * theta * tstep	&
		 + v * rho * (vpar0_x * ps0_y - vpar0_y * ps0_x)				* xjac * theta * tstep	&
		 + v * rho * F0 / BigR * vpar0_p						* xjac * theta * tstep	&
                 + D_perp_numm * (v_xx   + v_x  /BigR + v_yy  )								&
		               * (rho_xx + rho_x/BigR + rho_yy) * BigR 				* xjac * theta * tstep	

  amat_k(5,5)  = + (D_par-D_prof) * BigR / BB2 * Bgrad_rho_k_star * Bgrad_rho_rho		* xjac * theta * tstep 

  amat_n(5,5)  = + (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star   * Bgrad_rho_rho_n		* xjac * theta * tstep	&
		 + v * F0 / BigR * Vpar0 * rho_p						* xjac * theta * tstep 

  amat_kn(5,5) = + (D_par-D_prof) * BigR / BB2 * Bgrad_rho_k_star * Bgrad_rho_rho_n		* xjac * theta * tstep	&
		 + D_prof * BigR  * ( v_p*rho_p * eps_cyl**2 /BigR**2 ) 			* xjac * theta * tstep 

  amat(5,7)    = + v * F0 / BigR * Vpar * r0_p  						* xjac * theta * tstep	&
		 + v * Vpar * (r0_x * ps0_y - r0_y * ps0_x)					* xjac * theta * tstep	&
		 + v * r0 * (vpar_x * ps0_y - vpar_y * ps0_x)					* xjac * theta * tstep 

  amat_n(5,7)  = + v * r0 * F0 / BigR * vpar_p  						* xjac * theta * tstep
		    
  ! ------------------------------------
  ! --- The LHS terms (diamagnetic part)
  amat(5,5)    = amat(5,5)												&
                 - tau_IC * v * 2.d0 * (rho_y*T0 + rho*T0_y) * BigR				* xjac * theta * tstep	
  
  amat(5,6)    = amat(5,6)												&
                 - tau_IC * v * 2.d0 * (T_y*r0 + T*r0_y) * BigR					* xjac * theta * tstep 

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
  Bgrad_T          = ( T0_x * ps0_y - T0_y * ps0_x &
                     + F0 / BigR * T0_p              ) / BigR
  Bgrad_T_star     = ( v_x   * ps0_y - v_y   * ps0_x ) / BigR
  Bgrad_T_k_star   = ( F0 / BigR * v_p  	     ) / BigR
	      
  ! -----------------------------	      
  ! --- The RHS terms (main part)
  rhs(6) =   												&
             ! --- Time derivative
	     + zeta * v * r0 * BigR							* xjac *delta_g(6)&
             + zeta * v * T0 * BigR 							* xjac *delta_g(5)&
             ! --- Convective terms
	     + v * r0 * BigR**2 * ( T0_x * u0_y - T0_y * u0_x)				* xjac * tstep	&
	     + v * T0 * BigR**2 * ( r0_x * u0_y - r0_y * u0_x)				* xjac * tstep	&
	     + v * r0 * GAMMA * T0 * u0_y * 2.d0 * BigR					* xjac * tstep	&
	     - v * r0 *         F0 / BigR * Vpar0 * T0_p				* xjac * tstep	&
	     - v * T0 *         F0 / BigR * Vpar0 * r0_p				* xjac * tstep	&
	     - v * r0 * GAMMA * F0 / BigR * T0    * vpar0_p				* xjac * tstep	&
	     - v * r0 *      Vpar0 * (T0_x    * ps0_y - T0_y    * ps0_x)		* xjac * tstep	&
	     - v * T0 *      Vpar0 * (r0_x    * ps0_y - r0_y    * ps0_x)		* xjac * tstep	&
	     - v * r0 * GAMMA * T0 * (vpar0_x * ps0_y - vpar0_y * ps0_x)		* xjac * tstep	&
             ! --- Source
             + v * BigR * heat_source							* xjac * tstep	&
             ! --- Conductivity
	     - (K_par-K_prof) * BigR / BB2 * Bgrad_T_star * Bgrad_T			* xjac * tstep	&
	     - K_prof * BigR * (v_x*T0_x + v_y*T0_y )					* xjac * tstep	&
             ! --- Numerical conductivity
             - K_perp_numm * (v_xx  + v_x /Bigr + v_yy )						&
	                   * (T0_xx + T0_x/Bigr + T0_yy) * BigR 			* xjac * tstep  

  rhs_k(6) = - (K_par-K_prof) * BigR / BB2 * Bgrad_T_k_star * Bgrad_T			* xjac * tstep	&
	     - K_prof * BigR * ( v_p*T0_p /BigR**2 )					* xjac * tstep 

  
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
  
  ! --- Parallel gradient terms       
  Bgrad_T_star_psi  = ( v_x  * psi_y - v_y  * psi_x ) / BigR
  Bgrad_T_psi	    = ( T0_x * psi_y - T0_y * psi_x ) / BigR
  Bgrad_T_T	    = ( T_x  * ps0_y - T_y  * ps0_x ) / BigR
  Bgrad_T_T_n	    = ( F0 / BigR * T_p             ) / BigR

  ! -----------------------------	      
  ! --- The LHS terms (main part)
  amat(6,1)    = + v         * r0 * Vpar0 * (T0_x    * psi_y - T0_y    * psi_x)				* xjac * theta * tstep	&
	         + v         * T0 * Vpar0 * (r0_x    * psi_y - r0_y    * psi_x)				* xjac * theta * tstep	&
		 + v * GAMMA * r0 * T0    * (vpar0_x * psi_y - vpar0_y * psi_x)				* xjac * theta * tstep	&
		 - (K_par-K_prof) * BigR * BB2_psi / BB2**2 * Bgrad_T_star     * Bgrad_T		* xjac * theta * tstep	&
		 + (K_par-K_prof) * BigR / BB2 	            * Bgrad_T_star_psi * Bgrad_T		* xjac * theta * tstep	&
		 + (K_par-K_prof) * BigR / BB2 	            * Bgrad_T_star     * Bgrad_T_psi		* xjac * theta * tstep 

  amat_k(6,1)  = - (K_par-K_prof) * BigR * BB2_psi / BB2**2 * Bgrad_T_k_star   * Bgrad_T		* xjac * theta * tstep	&
		 + (K_par-K_prof) * BigR / BB2 	            * Bgrad_T_k_star   * Bgrad_T_psi		* xjac * theta * tstep 

  amat(6,2)    = - v * r0 * BigR**2 * ( T0_x * u_y - T0_y * u_x)					* xjac * theta * tstep	&
	         - v * T0 * BigR**2 * ( r0_x * u_y - r0_y * u_x)					* xjac * theta * tstep	&
		 - v * 2.d0 * GAMMA * r0 * BigR * T0 * u_y						* xjac * theta * tstep 

  amat(6,5)    = + v * rho * T0 * BigR 									* xjac * (1.d0 + zeta)  &
                 - v * rho * BigR**2 * (T0_x  * u0_y  - T0_y  * u0_x )					* xjac * theta * tstep	&
	         - v * T0  * BigR**2 * (rho_x * u0_y  - rho_y * u0_x )					* xjac * theta * tstep	&
		 + v * rho * Vpar0   * (T0_x  * ps0_y - T0_y  * ps0_x)					* xjac * theta * tstep	&
	         + v * T0  * Vpar0   * (rho_x * ps0_y - rho_y * ps0_x)					* xjac * theta * tstep	&
		 + v * rho * Vpar0   * F0/BigR * T0_p							* xjac * theta * tstep	&
		 - v * GAMMA * rho * T0 * u0_y * 2.d0 * BigR						* xjac * theta * tstep	&
		 + v * GAMMA * rho * T0 * F0/BigR * Vpar0_p						* xjac * theta * tstep	&
		 + v * GAMMA * rho * T0 * (Vpar0_x * ps0_y - Vpar0_y * ps0_x)				* xjac * theta * tstep

  amat_n(6,5)  = + v * T0 *                F0 / BigR * Vpar0 * rho_p					* xjac * theta * tstep
  
  amat(6,6)    = - v * r0 * BigR**2 * (T_x  * u0_y - T_y  * u0_x)					* xjac * theta * tstep	&
	         - v * T  * BigR**2 * (r0_x * u0_y - r0_y * u0_x)					* xjac * theta * tstep	&
	         + v * T  *                F0 / BigR * Vpar0 * r0_p					* xjac * theta * tstep	&
		 + v * r0 * Vpar0   * (T_x  * ps0_y - T_y  * ps0_x)					* xjac * theta * tstep	&
	         + v * T  * Vpar0   * (r0_x * ps0_y - r0_y * ps0_x)					* xjac * theta * tstep	&
		 - v * r0 * GAMMA * T * BigR * u0_y * 2.d0						* xjac * theta * tstep	&
		 + v * r0 * GAMMA * T * F0/BigR * Vpar0_p						* xjac * theta * tstep	&
		 + v * r0 * GAMMA * T * (Vpar0_x * ps0_y - Vpar0_y * ps0_x)				* xjac * theta * tstep	&
		 + (K_par-K_prof) * BigR / BB2 * Bgrad_T_star * Bgrad_T_T				* xjac * theta * tstep	&
		 + dK_par * T     * BigR / BB2 * Bgrad_T_star * Bgrad_T					* xjac * theta * tstep	&
		 + K_prof * BigR * (v_x*T_x + v_y*T_y )							* xjac * theta * tstep	& 
                 + K_perp_numm * (v_xx + v_x/BigR + v_yy)									&
		               * (T_xx + T_x/BigR + T_yy) * BigR 					* xjac * theta * tstep
  
  if (r0 .lt. rho_1) then
    amat(6,6)  = amat(6,6) + v * rho_1 * T * BigR							* xjac * (1.d0 + zeta)
  else
    amat(6,6)  = amat(6,6) + v * r0    * T * BigR							* xjac * (1.d0 + zeta)
  endif
  
  amat_k(6,6)  = + (K_par-K_prof) * BigR / BB2 * Bgrad_T_k_star * Bgrad_T_T				* xjac * theta * tstep	&
		 + dK_par * T     * BigR / BB2 * Bgrad_T_k_star * Bgrad_T				* xjac * theta * tstep 
	      
  amat_n(6,6)  = + (K_par-K_prof) * BigR / BB2 * Bgrad_T_star   * Bgrad_T_T_n  				* xjac * theta * tstep	&
		 + v * r0 * Vpar0  * F0/BigR * T_p							* xjac * theta * tstep 

  amat_kn(6,6) = + (K_par-K_prof) * BigR / BB2 * Bgrad_T_k_star * Bgrad_T_T_n  				* xjac * theta * tstep	&
		 +  K_prof	    * BigR	* (v_p*T_p /BigR**2 )					* xjac * theta * tstep 

  amat(6,7)    = + v * r0 * F0/BigR * Vpar * T0_p							* xjac * theta * tstep	&
	         + v * T0 *                F0 / BigR * Vpar * r0_p					* xjac * theta * tstep	&
		 + v * r0 * Vpar * (T0_x * ps0_y - T0_y * ps0_x)					* xjac * theta * tstep	&
	         + v * T0 * Vpar * (r0_x * ps0_y - r0_y * ps0_x)					* xjac * theta * tstep	&
		 + v * GAMMA * r0 * T0 * (Vpar_x * ps0_y - Vpar_y * ps0_x)				* xjac * theta * tstep 
     
  amat_n(6,7)  = + v * GAMMA * r0 * T0 * F0/BigR * Vpar_p  						* xjac * theta * tstep       

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
  
  ! -----------------------------	      
  ! --- The RHS terms (main part)
  rhs(7) = 													&
             ! --- Time derivative terms (including dB/dt)
	     + zeta * v * r0 * F0**2 / BigR							* xjac *delta_g(7)&  
             + zeta * v * r0 * vpar0 * (ps0_x * delta_ps_x + ps0_y * delta_ps_y) / BigR 	* xjac		&
             ! --- Convection terms
	     - 0.5d0 * r0 * vpar0**2 * BB2 * (ps0_x * v_y  - ps0_y * v_x)			* xjac * tstep	&
             - 0.5d0 * v  * vpar0**2 * BB2 * (ps0_x * r0_y - ps0_y * r0_x)			* xjac * tstep	&
             + 0.5d0 * v  * vpar0**2 * BB2 * F0 / BigR * r0_p					* xjac * tstep	&
             ! --- Parallel pressure gradient
             - v * F0 / BigR * P0_p								* xjac * tstep	&
	     - v * (P0_x * ps0_y - P0_y * ps0_x)						* xjac * tstep	&
             ! --- Density source term
             - v * total_rho_source * vpar0 * BB2 * BigR 					* xjac * tstep	&
             ! --- Viscosity and toroidal source
             + visco_par * (v_x * Vt0_x   + v_y * Vt0_y)   * BigR 				* xjac * tstep	&
	     - visco_par * (v_x * vpar0_x + v_y * vpar0_y) * BigR				* xjac * tstep	&
             ! --- Numerical viscosity
             - visco_par_numm * (v_xx     + v_x    /Bigr + v_yy    )						&
	                      * (vpar0_xx + vpar0_x/Bigr + vpar0_yy) * BigR 			* xjac * tstep 
  
  rhs_k(7) = + 0.5d0 * r0 * vpar0**2 * BB2 * F0 / BigR * v_p					* xjac * tstep 

  ! -----------------------------------    
  ! --- The RHS term (Neoclassical part)	      
  rhs(7)   = rhs(7)											 	&
             ! --- Neoclassic term
             + v * amu_neo_prof * BB2 / (Btheta2+epsil) * BigR							&
	         * (  r0                         * (ps0_x*u0_x + ps0_y*u0_y)					&
                    + tau_IC                     * (ps0_x*P0_x + ps0_y*P0_y)					&
                    + aki_neo_prof * tau_IC * r0 * (ps0_x*T0_x + ps0_y*T0_y)					&
                    - r0 * Vpar0 * Btheta2					)		* xjac * tstep 
  
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
  real*8 :: Vt_x_psi, Vt_y_psi
    
  Btheta2_psi = 2.d0 * (psi_x * ps0_x + psi_y * ps0_y ) / BigR**2
  Vt_x_psi    = dV_dpsi_source * psi_x
  Vt_y_psi    = dV_dpsi_source * psi_y
  
  ! -----------------------------
  ! --- The LHS terms (main part)
  amat(7,1)   = + v * (P0_x * psi_y - P0_y * psi_x)								* xjac * theta * tstep	&
                + 0.5d0 * r0 * vpar0**2 * BB2     * (psi_x * v_y  - psi_y * v_x)				* xjac * theta * tstep	&
                + 0.5d0 * r0 * vpar0**2 * BB2_psi * (ps0_x * v_y  - ps0_y * v_x)				* xjac * theta * tstep	&
                + 0.5d0 * v  * vpar0**2 * BB2     * (psi_x * r0_y - psi_y * r0_x)				* xjac * theta * tstep	&
                + 0.5d0 * v  * vpar0**2 * BB2_psi * (ps0_x * r0_y - ps0_y * r0_x)				* xjac * theta * tstep	&
                - 0.5d0 * v  * vpar0**2 * BB2_psi * F0 / BigR * r0_p						* xjac * theta * tstep	&
                + v * total_rho_source * vpar0 * BB2_psi * BigR 						* xjac * theta * tstep	&
                - visco_par * (v_x * Vt_x_psi   + v_y * Vt_y_psi)   * BigR 					* xjac * theta * tstep	&
                + v * r0 * vpar0 / BigR * (ps0_x * psi_x + ps0_y * psi_y)					* xjac * (1.d0 + zeta) 
  
  amat_k(7,1) = - 0.5d0 * r0 * vpar0**2 * BB2_psi * F0 / BigR * v_p						* xjac * theta * tstep 
           
  amat(7,5)   = + v * (rho_x * T0   * ps0_y - rho_y * T0   * ps0_x)						* xjac * theta * tstep	&
		+ v * (rho   * T0_x * ps0_y - rho   * T0_y * ps0_x)						* xjac * theta * tstep	&
		+ v * F0 / BigR * rho * T0_p									* xjac * theta * tstep	& 
		+ 0.5d0 * rho * vpar0**2 * BB2 * (ps0_x * v_y   - ps0_y * v_x)					* xjac * theta * tstep	&
                + 0.5d0 * v   * vpar0**2 * BB2 * (ps0_x * rho_y - ps0_y * rho_x)				* xjac * theta * tstep

  amat_k(7,5) = - 0.5d0 * rho * vpar0**2 * BB2 * F0 / BigR * v_p						* xjac * theta * tstep 

  amat_n(7,5) = + v * F0 / BigR * rho_p * T0									* xjac * theta * tstep	& 
                - 0.5d0 * v   * vpar0**2 * BB2 * F0 / BigR * rho_p						* xjac * theta * tstep 

  amat(7,6)   = + v * (T_x * r0   * ps0_y - T_y * r0   * ps0_x)							* xjac * theta * tstep	&
		+ v * (T   * r0_x * ps0_y - T	* r0_y * ps0_x)							* xjac * theta * tstep	&
		+ v * F0 / BigR * T * r0_p									* xjac * theta * tstep
  
  amat_n(7,6) = + v * F0 / BigR * T_p * r0									* xjac * theta * tstep 

  amat(7,7)   = + r0 * vpar0 * vpar * BB2 * (ps0_x * v_y  - ps0_y * v_x)					* xjac * theta * tstep	&
                + v  * vpar0 * vpar * BB2 * (ps0_x * r0_y - ps0_y * r0_x)					* xjac * theta * tstep	&
                - v  * vpar0 * vpar * BB2 * F0 / BigR * r0_p							* xjac * theta * tstep	&
                + v * total_rho_source  * vpar * BB2 * BigR  							* xjac * theta * tstep	&
                + visco_par * (v_x * Vpar_x + v_y * Vpar_y) * BigR						* xjac * theta * tstep	&
                + visco_par_numm * (v_xx    + v_x   /BigR + v_yy   )									&
		                 * (vpar_xx + vpar_x/BigR + vpar_yy) * BigR 					* xjac * theta * tstep 
  
  if (r0 .lt. rho_1) then
    amat(7,7) = amat(7,7) + v * Vpar * rho_1 * F0**2 / BigR							* xjac * (1.d0 + zeta)
  else
    amat(7,7) = amat(7,7) + v * Vpar * r0    * F0**2 / BigR							* xjac * (1.d0 + zeta)
  endif

  amat_k(7,7) = - r0 * vpar0 * vpar * BB2 * F0 / BigR * v_p							* xjac * theta * tstep 

  
  ! -----------------------------------
  ! --- The LHS terms (Neoclassic part)
  amat(7,1)   = amat(7,1)														&
		- v * amu_neo_prof * BB2 / (Btheta2+epsil) * BigR									&
		                   * (  r0                         * (psi_x*u0_x + psi_y*u0_y)						&
                                      + tau_IC                     * (psi_x*P0_x + psi_y*P0_y)						&
                                      + aki_neo_prof * tau_IC * r0 * (psi_x*T0_x + psi_y*T0_y)  )		* xjac * theta * tstep	&
                + v * amu_neo_prof * Btheta2_psi * BB2 / Btheta2**2.d0 * BigR								&
                                   * (  r0                         * (ps0_x*u0_x + ps0_y*u0_y)						&
                                      + tau_IC                     * (ps0_x*P0_x + ps0_y*P0_y)						&
                                      + aki_neo_prof * tau_IC * r0 * (ps0_x*T0_x + ps0_y*T0_y)  )		* xjac * theta * tstep
           
  amat(7,2)   = amat(7,2)														&
                - v * amu_neo_prof * BB2 / (Btheta2+epsil) * r0 * (ps0_x*u_x + ps0_y*u_y) * BigR		* xjac * theta * tstep 
  
  amat(7,5)   = amat(7,5)														&
                - v * amu_neo_prof * BB2 / (Btheta2+epsil) * BigR									&
                                   * (  rho * (ps0_x*u0_x + ps0_y*u0_y)									&
                                      + tau_IC                      * (ps0_x*rho_x*T0   + ps0_y*rho_y*T0  ) 				&
                                      + tau_IC                      * (ps0_x*rho  *T0_x + ps0_y*rho  *T0_y) 				&
                                      + aki_neo_prof * tau_IC * rho * (ps0_x*T0_x       + ps0_y*T0_y)					&
                                      - rho * Vpar0 * Btheta2						      )	* xjac * tstep * theta

  amat(7,6)   = amat(7,6)														&
                - v * amu_neo_prof * BB2 / (Btheta2+epsil) * BigR									&
                                   * (  tau_IC                     * (ps0_x*r0_x*T   + ps0_y*r0_y*T  )					&
                                      + tau_IC                     * (ps0_x*r0  *T_x + ps0_y*r0  *T_y)					&
                                      + aki_neo_prof * tau_IC * r0 * (ps0_x*T_x      + ps0_y*T_y     )   )	* xjac * tstep * theta
  
  amat(7,7)   = amat(7,7)														&
                + v * amu_neo_prof * BB2 * r0 * vpar * BigR 							* xjac * tstep * theta 


  return

end subroutine ELM_main_lhs_7










