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
  rhs(1) = 												&
           ! --- Time derivative
	   + zeta * v / R								* xjac * delta_g(1)&
           ! --- Resistive term
           + v * eta_Te  * (zj0 - current_source - jb)/ R				* xjac * tstep  &
           ! --- VxB
	   + v * (ps0_x * u0_y - ps0_y * u0_x)  					* xjac * tstep  &
           ! --- Integration term
	   - v * eps_cyl * F0 / R  * u0_p						* xjac * tstep  &
           ! --- Numerical resistivity
	   + eta_numm * (v_x * zj0_x + v_y * zj0_y)					* xjac * tstep  
  
  ! -----------------------------------    
  ! --- The RHS term (diamagnetic part)	      
  rhs(1) = rhs(1)											&
           - v * tau_IC/(rho_corr*BB2) * F0**2/R**2 * (ps0_x*pe0_y - ps0_y*pe0_x)	* xjac * tstep  & 		  
           + v * tau_IC/(rho_corr*BB2) * F0**3/R**3 * eps_cyl * pe0_p			* xjac * tstep  
  
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
  amat(1,1)   = + v * psi / R										* xjac * (1.d0+zeta)	&
		- v * (psi_x * u0_y - psi_y * u0_x)							* xjac * theta * tstep	

  amat(1,2)   = -  v * (ps0_x * u_y - ps0_y * u_x)							* xjac * theta * tstep

  amat_n(1,2) = +  eps_cyl * F0 / R * v * u_p								* xjac * theta * tstep

  amat(1,3)   = - eta_numm * (v_x * zj_x + v_y * zj_y)							* xjac * theta * tstep	&
	        - eta_Te * v * zj / R									* xjac * theta * tstep

  amat(1,8)   = - deta_dTe * v * Te * (zj0 - current_source - jb) / R					* xjac * theta * tstep	

  ! ------------------------------------
  ! --- The LHS terms (diamagnetic part)
  amat(1,1)   = amat(1,1)													&
                - v * tau_IC/(rho_corr*BB2**2) * BB2_psi * F0**2/R**2 * (ps0_x*pe0_y - ps0_y*pe0_x)	* xjac * theta * tstep	&
                + v * tau_IC/(rho_corr*BB2**2) * BB2_psi * F0**3/R**3 * eps_cyl * pe0_p			* xjac * theta * tstep	&
                + v * tau_IC/(rho_corr*BB2)              * F0**2/R**2 * (psi_x*pe0_y - psi_y*pe0_x)	* xjac * theta * tstep
  
  amat(1,5)   = amat(1,5)													&
                + v * tau_IC/(rho_corr*BB2)             * F0**2/R**2 * Te0 * (ps0_x*rho_y - ps0_y*rho_x)* xjac * theta * tstep	&
	        + v * tau_IC/(rho_corr*BB2)             * F0**2/R**2 * rho * (ps0_x*Te0_y - ps0_y*Te0_x)* xjac * theta * tstep	&
	        - v * tau_IC/(rho_corr*BB2)             * F0**3/R**3 * eps_cyl * rho * Te0_p		* xjac * theta * tstep	&
		- v * tau_IC * rho /(rho_corr**2 * BB2) * F0**2/R**2 * (ps0_x*pe0_y - ps0_y*pe0_x)	* xjac * theta * tstep	&		    
		+ v * tau_IC * rho /(rho_corr**2 * BB2) * F0**3/R**3 * eps_cyl * pe0_p			* xjac * theta * tstep 
  
  amat_n(1,5) = amat_n(1,5)													&
                - v * tau_IC/(rho_corr*BB2) * F0**3/R**3 * eps_cyl * Te0 * rho_p			* xjac * theta * tstep
  
  amat(1,8)   = amat(1,8)													&
                + v * tau_IC/(rho_corr*BB2) * F0**2/R**2 * r0 * (ps0_x*Te_y - ps0_y*Te_x) 		* xjac * theta * tstep	&
	        + v * tau_IC/(rho_corr*BB2) * F0**2/R**2 * Te * (ps0_x*r0_y - ps0_y*r0_x)		* xjac * theta * tstep	&
		- v * tau_IC/(rho_corr*BB2) * F0**3/R**3 * eps_cyl * Te * r0_p				* xjac * theta * tstep 

  amat_n(1,8) = amat_n(1,8)													&
                - v * tau_IC/(rho_corr*BB2) * F0**3/R**3 * eps_cyl * r0 * Te_p				* xjac * theta * tstep 
  
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
  rhs(2) = 															&
           ! --- Time derivative
	   - zeta * R * r0_hat * (v_x * delta_u_x + v_y * delta_u_y)  						* xjac		&		  
           ! --- Convective terms
	   - 0.5d0 * vv2 * (v_x * r0_y_hat - v_y * r0_x_hat)		 					* xjac * tstep	&
	   - r0_hat * R**2 * w0 * (v_x * u0_y - v_y * u0_x)  	 						* xjac * tstep	&
           ! --- [psi,j]
	   + v * (ps0_x * zj0_y - ps0_y * zj0_x )			 					* xjac * tstep	&
	   - v * eps_cyl * F0 / R * zj0_p				 					* xjac * tstep	&
           ! --- Grad(p)
	   + R**2 * (v_x * p0_y     - v_y * p0_x)			 					* xjac * tstep	&
           ! --- Source interaction
	   + R**3 * particle_source * (v_x * u0_x + v_y * u0_y)							* xjac * tstep	&
           ! --- Viscosity
	   - visco_Te * R * (v_x * w0_x + v_y * w0_y) 		 						* xjac * tstep	&
           ! --- Numerical iscosity
           - visco_numm * (v_xx  + v_x /R + v_yy ) 										&
	                * (w0_xx + w0_x/R + w0_yy)								* xjac * tstep 
  
  ! --------------------------------------------------      
  ! --- The RHS term (diamagnetic and neoclassic part)	      
  rhs(2) = rhs(2)														&
           ! --- Diamagnetic terms	
           - tau_IC * R**3 * Pi0_y * (v_x*u0_x  + v_y*u0_y )                                    		* xjac * tstep  &
           - tau_IC * R**4 * Pi0_y * (v_x*u0_xy + v_y*u0_xy)                                    		* xjac * tstep  &
           + tau_IC * R**4 * Pi0_x * (v_x*u0_xy + v_y*u0_yy)                                    		* xjac * tstep  &
           ! --- Inverse diamagnetic terms (needed when including diamagnetic vorticity directly into W - equation4)
           + tau_IC * W_dia * R**4 *            (Pi0_xx + Pi0_x/R + Pi0_yy)    * (v_x * u0_y - v_y * u0_x)	* xjac * tstep  &
           - tau_IC * W_dia * R**4 / rho_corr * (r0_x * Pi0_x + r0_y * Pi0_y)  * (v_x * u0_y - v_y * u0_x)	* xjac * tstep  &
	   ! --- Neoclassic term
           + amu_neo_prof * BB2 / (Btheta2+epsil)**2.d0 * (ps0_x*v_x + ps0_y*v_y) * R						&
                    * (  r0                         * (ps0_x*u0_x  + ps0_y*u0_y)						&
		       + tau_IC                     * (ps0_x*Pi0_x + ps0_y*Pi0_y)						&
                       + aki_neo_prof * tau_IC * r0 * (ps0_x*Ti0_x + ps0_y*Ti0_y)						&
                       - r0 * Vpar0 * Btheta2					     )				* xjac * tstep 

  
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
  
  rho_hat     = R**2 * rho
  rho_x_hat   = 2.d0 * R * R_x  * rho + R**2 * rho_x
  rho_y_hat   = R**2 * rho_y
  Btheta2_psi = 2.d0 * (psi_x * ps0_x + psi_y * ps0_y ) / R**2

  ! -----------------------------	      
  ! --- The LHS terms (main part)
  amat(2,1)   = - v * (psi_x * zj0_y - psi_y * zj0_x )								* xjac * theta * tstep

  amat(2,2)   = + r0_hat * R**2 * w0 * (v_x * u_y  - v_y  * u_x)						* xjac * theta * tstep	&
		+ R**2 * (u_x*u0_x + u_y*u0_y) * (v_x*r0_y_hat - v_y*r0_x_hat)					* xjac * theta * tstep	&
		- R**3 * particle_source * (v_x * u_x + v_y * u_y)						* xjac * theta * tstep
  
  if (r0 .lt. rho_1) then
    amat(2,2) = amat(2,2) - R**3 * rho_1 * (v_x * u_x + v_y * u_y)						* xjac * (1.d0 + zeta) 
  else
    amat(2,2) = amat(2,2) - R**3 * r0    * (v_x * u_x + v_y * u_y)						* xjac * (1.d0 + zeta) 
  endif

  amat(2,3)   = - v * (ps0_x * zj_y  - ps0_y * zj_x)								* xjac * theta * tstep

  amat_n(2,3) = + eps_cyl * F0 / R * v * zj_p									* xjac * theta * tstep

  amat(2,4)   = r0_hat * R**2 * w  * ( v_x * u0_y - v_y * u0_x)							* xjac * theta * tstep	&
		+ visco_Te * R * ( v_x * w_x + v_y * w_y)							* xjac * theta * tstep	&
                + visco_numm * (v_xx + v_x/R + v_yy) 										&
		             * (w_xx + w_x/R + w_yy) 								* xjac * theta * tstep    

  amat(2,5)   = + 0.5d0 * vv2 * (v_x * rho_y_hat - v_y * rho_x_hat)						* xjac * theta * tstep	&
		+ rho_hat * R**2 * w0 * (v_x * u0_y - v_y * u0_x)						* xjac * theta * tstep	&
		- R**2 * (v_x * rho_y * T0   - v_y * rho_x * T0  )						* xjac * theta * tstep	&
		- R**2 * (v_x * rho   * T0_y - v_y * rho   * T0_x)						* xjac * theta * tstep

  amat(2,6)   = - R**2 * (v_x * r0_y * Ti   - v_y * r0_x * Ti)							* xjac * theta * tstep	&
		- R**2 * (v_x * r0	* Ti_y - v_y * r0   * Ti_x)						* xjac * theta * tstep 
  
  amat(2,8)   = - R**2 * (v_x * r0_y * Te   - v_y * r0_x * Te)							* xjac * theta * tstep	&
		- R**2 * (v_x * r0	* Te_y - v_y * r0   * Te_x)						* xjac * theta * tstep	&
		+ dvisco_dTe * Te * ( v_x * w0_x    + v_y * w0_y    ) * R					* xjac * theta * tstep 
  
  
  ! ---------------------------------------------------    
  ! --- The LHS terms (diamagnetic and neoclassic part)
  amat(2,1)   = amat(2,1)														&
	        ! --- Neoclassical term
                - amu_neo_prof * BB2 / (Btheta2+epsil)**2 * (psi_x*v_x+psi_y*v_y) * R							&
                               * (  r0                         * (ps0_x*u0_x  + ps0_y*u0_y)						&
                                  + tau_IC                     * (ps0_x*Pi0_x + ps0_y*Pi0_y)						&
                                  + aki_neo_prof * tau_IC * r0 * (ps0_x*Ti0_x + ps0_y*Ti0_y)						&
                                  - r0 * Vpar0 * Btheta2)	 						* xjac * theta * tstep	&
                - amu_neo_prof * BB2 / (Btheta2+epsil)**2 * (ps0_x*v_x+ps0_y*v_y) * R							&
		               * (  r0                         * (psi_x*u0_x  + psi_y*u0_y)						&
                                  + tau_IC                     * (psi_x*Pi0_x + psi_y*Pi0_y)						&
                                  + aki_neo_prof * tau_IC * r0 * (psi_x*Ti0_x + psi_y*Ti0_y)	)		* xjac * theta * tstep	&
                + amu_neo_prof * BB2 * 2.d0*Btheta2_psi / (Btheta2+epsil)**3 * (ps0_x*v_x+ps0_y*v_y) * R				&
                               * (  r0                         * (ps0_x*u0_x  + ps0_y*u0_y)						&
                                  + tau_IC                     * (ps0_x*Pi0_x + ps0_y*Pi0_y)						&
                                  + aki_neo_prof * tau_IC * r0 * (ps0_x*Ti0_x + ps0_y*Ti0_y)	)		* xjac * theta * tstep	&
                - amu_neo_prof * BB2 * Btheta2_psi / (Btheta2+epsil)**2									&
		               * r0 * vpar0 * (ps0_x*v_x + ps0_y*v_y) * R					* xjac * tstep * theta

  amat(2,2)   = amat(2,2)														&
                ! --- Diamagnetic terms
                + tau_IC * R**3 * Pi0_y * (v_x*u_x  + v_y*u_y )                                    		* xjac * theta * tstep  &
                + tau_IC * R**4 * Pi0_y * (v_x*u_xy + v_y*u_xy)                                    		* xjac * theta * tstep  &
                - tau_IC * R**4 * Pi0_x * (v_x*u_xy + v_y*u_yy)                                    		* xjac * theta * tstep  &
                ! --- Inverse diamagnetic terms
                - tau_IC * W_dia * R**4            * (Pi0_xx + Pi0_x/R + Pi0_yy) * (v_x*u_y - v_y*u_x)		* xjac * theta * tstep  &
                + tau_IC * W_dia * R**4 / rho_corr * (r0_x*Pi0_x + r0_y*Pi0_y )  * (v_x*u_y - v_y*u_x)		* xjac * theta * tstep  &
	        ! --- Neoclassical term
                - amu_neo_prof * BB2 / (Btheta2+epsil)**2 * (ps0_x*v_x + ps0_y*v_y) * R							&
		               * r0 * (ps0_x*u_x + ps0_y*u_y)							* xjac * theta * tstep
  
  amat(2,5)   = amat(2,5)														&
                ! --- Diamagnetic terms
                + tau_IC * R**3 * (rho*Ti0_y + rho_y*Ti0) * (v_x*u0_x  + v_y*u0_y )				* xjac * theta * tstep  &
                + tau_IC * R**4 * (rho*Ti0_y + rho_y*Ti0) * (v_x*u0_xy + v_y*u0_xy)				* xjac * theta * tstep  &
                - tau_IC * R**4 * (rho*Ti0_x + rho_x*Ti0) * (v_x*u0_xy + v_y*u0_yy)				* xjac * theta * tstep  &
                ! --- Inverse diamagnetic terms
                - tau_IC * W_dia * R**4 * ( rho_xx*Ti0 + rho*Ti0_xx + 2.d0*rho_x*Ti0_x                                                  &
                                          + rho_x*Ti0/R + rho*Ti0_x/R                                                                   &
                                          + rho_yy*Ti0 + rho*Ti0_yy + 2.d0*rho_y*Ti0_y )* (v_x*u0_y - v_y*u0_x) * xjac * theta * tstep  &
                + tau_IC * W_dia * R**4 / rho_corr * (rho_x*Pi0_x + rho_y*Pi0_y)        * (v_x*u0_y - v_y*u0_x) * xjac * theta * tstep  &
                + tau_IC * W_dia * R**4 / rho_corr * ( r0_x * (rho_x*Ti0 + rho*Ti0_x)							&
                                                     + r0_y * (rho_y*Ti0 + rho*Ti0_y) ) * (v_x*u0_y - v_y*u0_x) * xjac * theta * tstep  &
                - tau_IC * W_dia * R**4 * rho/rho_corr**2.d0 * (r0_x*Pi0_x + r0_y*Pi0_y)* (v_x*u0_y - v_y*u0_x) * xjac * theta * tstep  &
	        ! --- Neoclassical term
                - amu_neo_prof * BB2 / (Btheta2+epsil)**2 * (ps0_x*v_x + ps0_y*v_y) * R							&
                               * (  rho                         * (ps0_x*u0_x                  + ps0_y*u0_y)				&
                                  + tau_IC                      * (ps0_x*rho_x*Ti0   + ps0_y*rho_y*Ti0  )				&
                                  + tau_IC                      * (ps0_x*rho  *Ti0_x + ps0_y*rho  *Ti0_y)				&
                                  + aki_neo_prof * tau_IC * rho * (ps0_x*Ti0_x                 + ps0_y*Ti0_y)				&
                                  -rho * Vpar0 * Btheta2						      )	* xjac * tstep * theta

  amat(2,6)   = amat(2,6)														&
                ! --- Diamagnetic terms
                + tau_IC * R**3 * (r0*Ti_y + r0_y*Ti) * (v_x*u0_x  + v_y*u0_y )					* xjac * theta * tstep  &
                + tau_IC * R**4 * (r0*Ti_y + r0_y*Ti) * (v_x*u0_xy + v_y*u0_xy)					* xjac * theta * tstep  &
                - tau_IC * R**4 * (r0*Ti_x + r0_x*Ti) * (v_x*u0_xy + v_y*u0_yy)					* xjac * theta * tstep  &
                ! --- Inverse diamagnetic terms
                - tau_IC * W_dia * R**4 * ( r0_xx*Ti + r0*Ti_xx + 2.d0*r0_x*Ti_x                                                        &
                                          + r0_x*Ti/R + r0*Ti_x/R                                                                       &
                                          + r0_yy*Ti + r0*Ti_yy + 2.d0*r0_y*Ti_y  ) * (v_x*u0_y - v_y*u0_x)	* xjac * theta * tstep  &
                + tau_IC * W_dia * R**4 / rho_corr * ( r0_x * (r0_x*Ti + r0*Ti_x)							&
                                                     + r0_y * (r0_y*Ti + r0*Ti_y) ) * (v_x*u0_y - v_y*u0_x)	* xjac * theta * tstep  &
	        ! --- Neoclassical term
                - amu_neo_prof * BB2 / (Btheta2+epsil)**2 * (ps0_x*v_x + ps0_y*v_y) * R							&
                               * (  tau_IC                     * (ps0_x*r0_x*Ti   + ps0_y*r0_y*Ti  )					&
                                  + tau_IC                     * (ps0_x*r0  *Ti_x + ps0_y*r0  *Ti_y)					&
                                  + aki_neo_prof * tau_IC * r0 * (ps0_x*Ti_x+ps0_y*Ti_y)		   )	* xjac * tstep * theta
  
  amat(2,7)   = amat(2,7)														&
	        ! --- Neoclassical term
                + amu_neo_prof * BB2 / (Btheta2+epsil) * r0 * vpar * (ps0_x*v_x + ps0_y*v_y) * R		* xjac * tstep * theta 
  
  
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
  !rhs(3) = - ( v_x * ps0_x  + v_y * ps0_y + v*zj0) / R * xjac * tstep

  
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
  amat(3,1) = (v_x * psi_x + v_y * psi_y ) / R		* xjac * tstep

  amat(3,3) = v * zj / R				* xjac * tstep 

  
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
  rhs(4) = 0.d0
  
  ! ----------------------------	      
  ! --- The RHS term (main part)	      
  !rhs(4) = - ( v_x * u0_x   + v_y * u0_y  + v*w0) * R 				* xjac * tstep	&
  
  ! -----------------------------------    
  ! --- The RHS term (diamagnetic part)	      
  !rhs(4) = rhs(4)										&
  !	    - tau_IC / rho_corr * ( v_x  * Pi0_x + v_y  * Pi0_y) * R 		* xjac * tstep	

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
  amat(4,2) = (v_x * u_x + v_y * u_y) * R  								* xjac * tstep

  amat(4,4) =  v * w * R		      								* xjac * tstep 
		    
  ! ------------------------------------
  ! --- The LHS terms (diamagnetic part)
  amat(4,5) = amat(4,5)									        			&
              + tau_IC * W_dia / rho_corr              * (v_x * rho  *Ti0_x + v_y * rho  *Ti0_y) * R	* xjac * tstep	&
              + tau_IC * W_dia / rho_corr              * (v_x * rho_x*Ti0   + v_y * rho_y*Ti0  ) * R	* xjac * tstep	&
              - tau_IC * W_dia * rho / rho_corr**2.d0  * (v_x * Pi0_x       + v_y * Pi0_y      ) * R	* xjac * tstep	

  amat(4,6) = amat(4,6)									        		        &
              + tau_IC * W_dia / rho_corr              * (v_x * r0   *Ti_x  + v_y * r0   *Ti_y ) * R	* xjac * tstep	&
              + tau_IC * W_dia / rho_corr              * (v_x * r0_x *Ti    + v_y * r0_y *Ti   ) * R	* xjac * tstep	

  
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
                     + F0 / R * r0_p               ) / R
  Bgrad_rho_star   = ( v_x  * ps0_y - v_y  * ps0_x ) / R
  Bgrad_rho_k_star = ( F0 / R * v_p  	           ) / R
	      
  ! ----------------------------	      
  ! --- The RHS term (main part)	      
  rhs(5)   =   											&
             ! --- Time derivative
	     + zeta * v * R							* xjac *delta_g(5)& 
             ! --- Div(rho.V)
	     + v * R**2 * ( r0_x * u0_y - r0_y * u0_x)				* xjac * tstep	&
	     + v * 2.d0 * R * r0 * u0_y						* xjac * tstep	&
	     - v * F0 / R * Vpar0 * r0_p					* xjac * tstep	&
	     - v * F0 / R * r0    * vpar0_p					* xjac * tstep	&
	     - v * Vpar0 * (r0_x    * ps0_y - r0_y    * ps0_x)			* xjac * tstep	&
	     - v * r0    * (vpar0_x * ps0_y - vpar0_y * ps0_x)			* xjac * tstep	&
             ! --- Source
             + v * R * particle_source						* xjac * tstep	&
             ! --- Diffusivity
	     - (D_par-D_prof) * R / BB2 * Bgrad_rho_star * Bgrad_rho		* xjac * tstep	&
	     - D_prof * R  * (v_x*r0_x + v_y*r0_y)				* xjac * tstep	&
             ! --- Numerical diffusivity
             - D_perp_numm * (v_xx  + v_x /R + v_yy )						&
	                   * (r0_xx + r0_x/R + r0_yy) * R 			* xjac * tstep 

  rhs_k(5) = - (D_par-D_prof) * R / BB2 * Bgrad_rho_k_star * Bgrad_rho		* xjac * tstep	&
	     - D_prof * R  * ( v_p*r0_p * eps_cyl**2 /R**2 )			* xjac * tstep 
  
  ! -----------------------------------    
  ! --- The RHS term (diamagnetic part)	      
  rhs(5)   = rhs(5)									 	&
	     + tau_IC * v * 2.d0 * pi0_y * R					* xjac * tstep	
  
  
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
  Bgrad_rho_star_psi = ( v_x   * psi_y - v_y   * psi_x ) / R
  Bgrad_rho_psi      = ( r0_x  * psi_y - r0_y  * psi_x ) / R
  Bgrad_rho_rho      = ( rho_x * ps0_y - rho_y * ps0_x ) / R
  Bgrad_rho_rho_n    = ( F0 / R * rho_p ) / R

  ! -----------------------------	      
  ! --- The LHS terms (main part)
  amat(5,1)    = - (D_par-D_prof) * R * BB2_psi/ BB2**2 * Bgrad_rho_star	* Bgrad_rho	* xjac * theta * tstep	&
		 + (D_par-D_prof) * R / BB2  	        * Bgrad_rho_star_psi    * Bgrad_rho	* xjac * theta * tstep	&
		 + (D_par-D_prof) * R / BB2  	        * Bgrad_rho_star	* Bgrad_rho_psi * xjac * theta * tstep	&
		 + v * Vpar0 * (r0_x * psi_y - r0_y * psi_x)					* xjac * theta * tstep	&
		 + v * r0 * (vpar0_x * psi_y - vpar0_y * psi_x) 				* xjac * theta * tstep 

  amat_k(5,1)  = - (D_par-D_prof) * R * BB2_psi/ BB2**2 * Bgrad_rho_k_star * Bgrad_rho		* xjac * theta * tstep	&
		 + (D_par-D_prof) * R / BB2  	        * Bgrad_rho_k_star * Bgrad_rho_psi	* xjac * theta * tstep 

  amat(5,2)    = - v * R**2 * ( r0_x * u_y - r0_y * u_x)					* xjac * theta * tstep	&
		 - v * 2.d0 * R * r0 * u_y							* xjac * theta * tstep 

  amat(5,5)    = + v * rho * R									* xjac * (1.d0 + zeta)	&
		 - v * R**2 * ( rho_x * u0_y - rho_y * u0_x) 					* xjac * theta * tstep	&
		 - v * 2.d0 * R * rho * u0_y 							* xjac * theta * tstep	&
		 + (D_par-D_prof) * R / BB2 * Bgrad_rho_star * Bgrad_rho_rho 			* xjac * theta * tstep	&
		 + D_prof * R  * (v_x*rho_x + v_y*rho_y )					* xjac * theta * tstep	&
		 + v * Vpar0 * (  rho_x * ps0_y -   rho_y * ps0_x)  				* xjac * theta * tstep	&
		 + v * rho   * (vpar0_x * ps0_y - vpar0_y * ps0_x)				* xjac * theta * tstep	&
		 + v * rho * F0 / R * vpar0_p							* xjac * theta * tstep	&
                 + D_perp_numm * (v_xx   + v_x  /R + v_yy  )								&
		               * (rho_xx + rho_x/R + rho_yy) * R 				* xjac * theta * tstep	

  amat_k(5,5)  = + (D_par-D_prof) * R / BB2 * Bgrad_rho_k_star * Bgrad_rho_rho			* xjac * theta * tstep 

  amat_n(5,5)  = + (D_par-D_prof) * R / BB2 * Bgrad_rho_star   * Bgrad_rho_rho_n		* xjac * theta * tstep	&
		 + v * F0 / R * Vpar0 * rho_p							* xjac * theta * tstep 

  amat_kn(5,5) = + (D_par-D_prof) * R / BB2 * Bgrad_rho_k_star * Bgrad_rho_rho_n		* xjac * theta * tstep	&
		 + D_prof * R  * ( v_p*rho_p * eps_cyl**2 /R**2 ) 				* xjac * theta * tstep 

  amat(5,7)    = + v * F0 / R * Vpar * r0_p  							* xjac * theta * tstep	&
		 + v * Vpar * (  r0_x * ps0_y -   r0_y * ps0_x)					* xjac * theta * tstep	&
		 + v * r0   * (vpar_x * ps0_y - vpar_y * ps0_x)					* xjac * theta * tstep 

  amat_n(5,7)  = + v * r0 * F0 / R * vpar_p  							* xjac * theta * tstep
		    
  ! ------------------------------------
  ! --- The LHS terms (diamagnetic part)
  amat(5,5)    = amat(5,5)												&
                 - tau_IC * v * 2.d0 * (rho_y*Ti0 + rho*Ti0_y) * R				* xjac * theta * tstep	
  
  amat(5,6)    = amat(5,6)												&
                 - tau_IC * v * 2.d0 * (Ti_y*r0 + Ti*r0_y) * R					* xjac * theta * tstep 
  
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
  Bgrad_Ti          = ( Ti0_x * ps0_y - Ti0_y * ps0_x &
                      + F0 / R * Ti0_p                ) / R
  Bgrad_Ti_star     = ( v_x   * ps0_y - v_y   * ps0_x ) / R
  Bgrad_Ti_k_star   = ( F0 / R * v_p  	              ) / R
	      
  ! -----------------------------	      
  ! --- The RHS terms (main part)
  rhs(6) =   												&
             ! --- Time derivative
	     + zeta * v * r0 * R							* xjac *delta_g(6)&
             ! --- Convective terms
	     + v * r0 * R**2 * ( Ti0_x * u0_y - Ti0_y * u0_x)				* xjac * tstep	&
	     + v * r0 * (GAMMA-1.d0) * Ti0 * u0_y * 2.d0 * R				* xjac * tstep	&
	     - v * r0 *                F0 / R * Vpar0 * Ti0_p				* xjac * tstep	&
	     - v * r0 * (GAMMA-1.d0) * F0 / R * Ti0   * vpar0_p				* xjac * tstep	&
	     - v * r0 *              Vpar0 * (Ti0_x   * ps0_y - Ti0_y   * ps0_x)	* xjac * tstep	&
	     - v * r0 * (GAMMA-1.d0) * Ti0 * (vpar0_x * ps0_y - vpar0_y * ps0_x)	* xjac * tstep	&
             ! --- Source
             + v * R * heat_source_i							* xjac * tstep	&
             ! --- Conductivity
	     - (Ki_par-Ki_prof) * R / BB2 * Bgrad_Ti_star * Bgrad_Ti			* xjac * tstep	&
	     - Ki_prof * R * (v_x*Ti0_x + v_y*Ti0_y )					* xjac * tstep	&
             ! --- Numerical conductivity
	     - Ki_par_num * (v_ps0_x   * ps0_y - v_ps0_y   * ps0_x)					&
		 	  * (Ti0_ps0_x * ps0_y - Ti0_ps0_y * ps0_x)			* xjac * tstep	&
             - Ki_perp_numm * (v_xx   + v_x  /R + v_yy  )						&
	                    * (Ti0_xx + Ti0_x/R + Ti0_yy) * R 				* xjac * tstep  

  rhs_k(6) = - (Ki_par-Ki_prof) * R / BB2 * Bgrad_Ti_k_star * Bgrad_Ti			* xjac * tstep	&
	     - Ki_par_num * (Ti0_pp  * v_pp)						* xjac * tstep	&
	     - Ki_prof * R * ( v_p*Ti0_p /R**2 )					* xjac * tstep 
  
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
  Bgrad_Ti_star_psi  = ( v_x   * psi_y - v_y   * psi_x ) / R
  Bgrad_Ti_psi	     = ( Ti0_x * psi_y - Ti0_y * psi_x ) / R
  Bgrad_Ti_Ti	     = ( Ti_x  * ps0_y - Ti_y  * ps0_x ) / R
  Bgrad_Ti_Ti_n	     = ( F0 / R * Ti_p   ) / R

  Ti_ps0_x  = Ti_xx  * ps0_y - Ti_xy  * ps0_x + Ti_x  * ps0_xy - Ti_y  * ps0_xx
  Ti_ps0_y  = Ti_xy  * ps0_y - Ti_yy  * ps0_x + Ti_x  * ps0_yy - Ti_y  * ps0_xy
  
  Ti0_psi_x = Ti0_xx * psi_y - Ti0_xy * psi_x + Ti0_x * psi_xy - Ti0_y * psi_xx
  Ti0_psi_y = Ti0_xy * psi_y - Ti0_yy * psi_x + Ti0_x * psi_yy - Ti0_y * psi_xy
  
  v_psi_x   = v_xx   * psi_y - v_xy   * psi_x + v_x   * psi_xy - v_y   * psi_xx
  v_psi_y   = v_xy   * psi_y - v_yy   * psi_x + v_x   * psi_yy - v_y   * psi_xy

  ! -----------------------------	      
  ! --- The LHS terms (main part)
  amat(6,1)    = + v * r0 * Vpar0 * (Ti0_x * psi_y - Ti0_y * psi_x)					* xjac * theta * tstep	&
		 + v * (GAMMA-1.d0) * r0 * Ti0 * (vpar0_x * psi_y - vpar0_y * psi_x)			* xjac * theta * tstep	&
		 + Ki_par_num * (v_psi_x  *ps0_y - v_psi_y  *ps0_x + v_ps0_x  *psi_y - v_ps0_y  *psi_x)				&
			      * (Ti0_ps0_x*ps0_y - Ti0_ps0_y*ps0_x)					* xjac * theta * tstep	&
		 + Ki_par_num * (Ti0_psi_x*ps0_y - Ti0_psi_y*ps0_x + Ti0_ps0_x*psi_y - Ti0_ps0_y*psi_x)				&
			      * (v_ps0_x *ps0_y - v_ps0_y *ps0_x)					* xjac * theta * tstep	&
		 - (Ki_par-Ki_prof) * R * BB2_psi / BB2**2 * Bgrad_Ti_star     * Bgrad_Ti		* xjac * theta * tstep	&
		 + (Ki_par-Ki_prof) * R / BB2 	           * Bgrad_Ti_star_psi * Bgrad_Ti		* xjac * theta * tstep	&
		 + (Ki_par-Ki_prof) * R / BB2 	           * Bgrad_Ti_star     * Bgrad_Ti_psi		* xjac * theta * tstep 

  amat_k(6,1)  = - (Ki_par-Ki_prof) * R * BB2_psi / BB2**2 * Bgrad_Ti_k_star   * Bgrad_Ti		* xjac * theta * tstep	&
		 + (Ki_par-Ki_prof) * R / BB2 	           * Bgrad_Ti_k_star   * Bgrad_Ti_psi		* xjac * theta * tstep 

  amat(6,2)    = - v * r0 * R**2 * ( Ti0_x * u_y - Ti0_y * u_x)						* xjac * theta * tstep	&
		 - v * 2.d0 * (GAMMA-1.d0) * r0 * R * Ti0 * u_y						* xjac * theta * tstep 

  amat(6,5)    = - v * rho * R**2 * (Ti0_x * u0_y  - Ti0_y * u0_x)					* xjac * theta * tstep	&
		 + v * rho * Vpar0   * (Ti0_x * ps0_y - Ti0_y * ps0_x)					* xjac * theta * tstep	&
		 + v * rho * Vpar0   * F0/R * Ti0_p							* xjac * theta * tstep	&
		 - v * (GAMMA-1.d0) * rho * Ti0 * u0_y * 2.d0 * R					* xjac * theta * tstep	&
		 + v * (GAMMA-1.d0) * rho * Ti0 * F0/R * Vpar0_p					* xjac * theta * tstep	&
		 + v * (GAMMA-1.d0) * rho * Ti0 * (Vpar0_x * ps0_y - Vpar0_y * ps0_x)			* xjac * theta * tstep 

  amat(6,6)    = - v * r0 * R**2    * (Ti_x * u0_y  - Ti_y * u0_x)					* xjac * theta * tstep	&
		 + v * r0 * Vpar0   * (Ti_x * ps0_y - Ti_y * ps0_x)					* xjac * theta * tstep	&
		 - v * r0 * (GAMMA-1.d0) * Ti * R * u0_y * 2.d0						* xjac * theta * tstep	&
		 + v * r0 * (GAMMA-1.d0) * Ti * F0/R * Vpar0_p						* xjac * theta * tstep	&
		 + v * r0 * (GAMMA-1.d0) * Ti * (Vpar0_x * ps0_y - Vpar0_y * ps0_x)			* xjac * theta * tstep	&
		 + Ki_par_num * (v_ps0_x  * ps0_y - v_ps0_y  * ps0_x) 								&
			      * (Ti_ps0_x * ps0_y - Ti_ps0_y * ps0_x)					* xjac * theta * tstep	&
		 + (Ki_par-Ki_prof) * R / BB2 * Bgrad_Ti_star * Bgrad_Ti_Ti				* xjac * theta * tstep	&
		 + dKi_par * Ti     * R / BB2 * Bgrad_Ti_star * Bgrad_Ti				* xjac * theta * tstep	&
		 + Ki_prof * R * (v_x*Ti_x + v_y*Ti_y )							* xjac * theta * tstep	& 
                 + Ki_perp_numm * (v_xx  + v_x /R + v_yy )									&
		                * (Ti_xx + Ti_x/R + Ti_yy) * R 						* xjac * theta * tstep	
  
  if (r0 .lt. rho_1) then
    amat(6,6)  = amat(6,6) + v * rho_1 * Ti * R								* xjac * (1.d0 + zeta)
  else
    amat(6,6)  = amat(6,6) + v * r0    * Ti * R								* xjac * (1.d0 + zeta)
  endif
  
  amat_k(6,6)  = + (Ki_par-Ki_prof) * R / BB2 * Bgrad_Ti_k_star * Bgrad_Ti_Ti				* xjac * theta * tstep	&
		 + dKi_par * Ti     * R / BB2 * Bgrad_Ti_k_star * Bgrad_Ti				* xjac * theta * tstep 
	      
  amat_n(6,6)  = + (Ki_par-Ki_prof) * R / BB2 * Bgrad_Ti_star   * Bgrad_Ti_Ti_n  			* xjac * theta * tstep	&
		 + v * r0 * Vpar0  * F0/R * Ti_p							* xjac * theta * tstep 

  amat_kn(6,6) = + (Ki_par-Ki_prof) * R / BB2 * Bgrad_Ti_k_star * Bgrad_Ti_Ti_n  			* xjac * theta * tstep	&
		 +  Ki_par_num * (Ti_pp  * v_pp) 							* xjac * theta * tstep	&
		 +  Ki_prof	    * R	* (v_p*Ti_p /R**2 )						* xjac * theta * tstep 

  amat(6,7)    = + v * r0 * F0/R * Vpar * Ti0_p								* xjac * theta * tstep	&
		 + v                * r0 * Vpar * (Ti0_x  * ps0_y - Ti0_y  * ps0_x)			* xjac * theta * tstep	&
		 + v * (GAMMA-1.d0) * r0 * Ti0  * (Vpar_x * ps0_y - Vpar_y * ps0_x)			* xjac * theta * tstep 
     
  amat_n(6,7)  = v * (GAMMA-1.d0) * r0 * Ti0 * F0/R * Vpar_p  						* xjac * theta * tstep       
		   

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
	     + zeta * v * r0 * F0**2 / R							* xjac *delta_g(7)&  
             + zeta * v * r0 * vpar0 * (ps0_x * delta_ps_x + ps0_y * delta_ps_y) / R	 	* xjac		&
             ! --- Convection terms
	     - 0.5d0 * r0 * vpar0**2 * BB2 * (ps0_x * v_y  - ps0_y * v_x)			* xjac * tstep	&
             - 0.5d0 * v  * vpar0**2 * BB2 * (ps0_x * r0_y - ps0_y * r0_x)			* xjac * tstep	&
             + 0.5d0 * v  * vpar0**2 * BB2 * F0 / R * r0_p					* xjac * tstep	&
             ! --- Parallel pressure gradient
             - v * F0 / R * P0_p								* xjac * tstep	&
	     - v * (P0_x * ps0_y - P0_y * ps0_x)						* xjac * tstep	&
             ! --- Viscosity and toroidal source
             + visco_par * (v_x * Vt0_x   + v_y * Vt0_y)   * R	 				* xjac * tstep	&
	     - visco_par * (v_x * vpar0_x + v_y * vpar0_y) * R					* xjac * tstep	&
             ! --- Numerical viscosity
             - visco_par_numm * (v_xx     + v_x    /R + v_yy    )						&
	                      * (vpar0_xx + vpar0_x/R + vpar0_yy) * R	 			* xjac * tstep 
  
  rhs_k(7) = + 0.5d0 * r0 * vpar0**2 * BB2 * F0 / R * v_p					* xjac * tstep 

  ! -----------------------------------    
  ! --- The RHS term (Neoclassical part)	      
  rhs(7)   = rhs(7)											 	&
             ! --- Neoclassic term
             + v * amu_neo_prof * BB2 / (Btheta2+epsil) * R							&
	         * (  r0                         * (ps0_x*u0_x  + ps0_y*u0_y)					&
                    + tau_IC                     * (ps0_x*Pi0_x + ps0_y*Pi0_y)					&
                    + aki_neo_prof * tau_IC * r0 * (ps0_x*Ti0_x + ps0_y*Ti0_y)					&
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
    
  Btheta2_psi = 2.d0 * (psi_x * ps0_x + psi_y * ps0_y ) / R**2
  Vt_x_psi    = dV_dpsi_source * psi_x
  Vt_y_psi    = dV_dpsi_source * psi_y
  
  ! -----------------------------
  ! --- The LHS terms (main part)
  amat(7,1)   = + v * (P0_x * psi_y - P0_y * psi_x)								* xjac * theta * tstep	&
                + 0.5d0 * r0 * vpar0**2 * BB2     * (psi_x * v_y  - psi_y * v_x)				* xjac * theta * tstep	&
                + 0.5d0 * r0 * vpar0**2 * BB2_psi * (ps0_x * v_y  - ps0_y * v_x)				* xjac * theta * tstep	&
                + 0.5d0 * v  * vpar0**2 * BB2     * (psi_x * r0_y - psi_y * r0_x)				* xjac * theta * tstep	&
                + 0.5d0 * v  * vpar0**2 * BB2_psi * (ps0_x * r0_y - ps0_y * r0_x)				* xjac * theta * tstep	&
                - 0.5d0 * v  * vpar0**2 * BB2_psi * F0 / R * r0_p						* xjac * theta * tstep	&
                - visco_par * (v_x * Vt_x_psi   + v_y * Vt_y_psi)   * R 					* xjac * theta * tstep	&
                + v * r0 * vpar0 / R * (ps0_x * psi_x + ps0_y * psi_y)						* xjac * (1.d0 + zeta) 
  
  amat_k(7,1) = - 0.5d0 * r0 * vpar0**2 * BB2_psi * F0 / R * v_p						* xjac * theta * tstep 
           
  amat(7,5)   = + v * (rho_x * (Ti0 + Te0)     * ps0_y - rho_y * (Ti0 + Te0)     * ps0_x)			* xjac * theta * tstep	&
		+ v * (rho   * (Ti0_x + Te0_x) * ps0_y - rho   * (Ti0_y + Te0_y) * ps0_x)			* xjac * theta * tstep	&
		+ v * F0 / R * rho * (Ti0_p + Te0_p)								* xjac * theta * tstep	& 
		+ 0.5d0 * rho * vpar0**2 * BB2 * (ps0_x * v_y   - ps0_y * v_x)					* xjac * theta * tstep	&
                + 0.5d0 * v   * vpar0**2 * BB2 * (ps0_x * rho_y - ps0_y * rho_x)				* xjac * theta * tstep 

  amat_k(7,5) = - 0.5d0 * rho * vpar0**2 * BB2 * F0 / R * v_p							* xjac * theta * tstep 

  amat_n(7,5) = + v * F0 / R * rho_p * (Ti0 + Te0)								* xjac * theta * tstep	& 
                - 0.5d0 * v   * vpar0**2 * BB2 * F0 / R * rho_p							* xjac * theta * tstep 

  amat(7,6)   = + v * (Ti_x * r0   * ps0_y - Ti_y * r0   * ps0_x)						* xjac * theta * tstep	&
		+ v * (Ti   * r0_x * ps0_y - Ti	  * r0_y * ps0_x)						* xjac * theta * tstep	&
		+ v * F0 / R * Ti * r0_p									* xjac * theta * tstep 
  
  amat_n(7,6) = + v * F0 / R * Ti_p * r0									* xjac * theta * tstep 

  amat(7,7)   = + r0 * vpar0 * vpar * BB2 * (ps0_x * v_y  - ps0_y * v_x)					* xjac * theta * tstep	&
                + v  * vpar0 * vpar * BB2 * (ps0_x * r0_y - ps0_y * r0_x)					* xjac * theta * tstep	&
                - v  * vpar0 * vpar * BB2 * F0 / R * r0_p							* xjac * theta * tstep	&
                + visco_par * (v_x * Vpar_x + v_y * Vpar_y) * R							* xjac * theta * tstep	&
                + visco_par_numm * (v_xx    + v_x   /R + v_yy   )									&
		                 * (vpar_xx + vpar_x/R + vpar_yy) * R	 					* xjac * theta * tstep 
  
  if (r0 .lt. rho_1) then
    amat(7,7) = amat(7,7) + v * Vpar * rho_1 * F0**2 / R							* xjac * (1.d0 + zeta)
  else
    amat(7,7) = amat(7,7) + v * Vpar * r0    * F0**2 / R							* xjac * (1.d0 + zeta)
  endif

  amat_k(7,7) = - r0 * vpar0 * vpar * BB2 * F0 / R * v_p							* xjac * theta * tstep 

  amat(7,8)   = + v * (Te_x * r0   * ps0_y - Te_y * r0   * ps0_x)						* xjac * theta * tstep	&
		+ v * (Te   * r0_x * ps0_y - Te	  * r0_y * ps0_x)						* xjac * theta * tstep	&
		+ v * F0 / R * Te * r0_p									* xjac * theta * tstep 
  
  amat_n(7,8) = + v * F0 / R * Te_p * r0									* xjac * theta * tstep

  
  ! -----------------------------------
  ! --- The LHS terms (Neoclassic part)
  amat(7,1)   = amat(7,1)														&
		- v * amu_neo_prof * BB2 / (Btheta2+epsil) * R									&
		                   * (  r0                         * (psi_x*u0_x  + psi_y*u0_y)						&
                                      + tau_IC                     * (psi_x*Pi0_x + psi_y*Pi0_y)					&
                                      + aki_neo_prof * tau_IC * r0 * (psi_x*Ti0_x + psi_y*Ti0_y)  )		* xjac * theta * tstep	&
                + v * amu_neo_prof * Btheta2_psi * BB2 / Btheta2**2.d0 * R								&
                                   * (  r0                         * (ps0_x*u0_x  + ps0_y*u0_y)						&
                                      + tau_IC                     * (ps0_x*Pi0_x + ps0_y*Pi0_y)					&
                                      + aki_neo_prof * tau_IC * r0 * (ps0_x*Ti0_x + ps0_y*Ti0_y)  )		* xjac * theta * tstep
           
  amat(7,2)   = amat(7,2)														&
                - v * amu_neo_prof * BB2 / (Btheta2+epsil) * r0 * (ps0_x*u_x + ps0_y*u_y) * R			* xjac * theta * tstep 
  
  amat(7,5)   = amat(7,5)														&
                - v * amu_neo_prof * BB2 / (Btheta2+epsil) * R										&
                                   * (  rho * (ps0_x*u0_x + ps0_y*u0_y)									&
                                      + tau_IC                      * (ps0_x*rho_x*Ti0   + ps0_y*rho_y*Ti0  ) 				&
                                      + tau_IC                      * (ps0_x*rho  *Ti0_x + ps0_y*rho  *Ti0_y) 				&
                                      + aki_neo_prof * tau_IC * rho * (ps0_x*Ti0_x       + ps0_y*Ti0_y)					&
                                      - rho * Vpar0 * Btheta2						      )	* xjac * tstep * theta

  amat(7,6)   = amat(7,6)														&
                - v * amu_neo_prof * BB2 / (Btheta2+epsil) * R										&
                                   * (  tau_IC                     * (ps0_x*r0_x*Ti   + ps0_y*r0_y*Ti  )				&
                                      + tau_IC                     * (ps0_x*r0  *Ti_x + ps0_y*r0  *Ti_y)				&
                                      + aki_neo_prof * tau_IC * r0 * (ps0_x*Ti_x      + ps0_y*Ti_y     )   )	* xjac * tstep * theta
  
  amat(7,7)   = amat(7,7)														&
                + v * amu_neo_prof * BB2 * r0 * vpar * R 							* xjac * tstep * theta 


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
  Bgrad_Te          = ( Te0_x * ps0_y - Te0_y * ps0_x &
                      + F0 / R * Te0_p              ) / R
  Bgrad_Te_star     = ( v_x  * ps0_y - v_y  * ps0_x ) / R
  Bgrad_Te_k_star   = ( F0 / R * v_p  	            ) / R
	      
  ! -----------------------------	      
  ! --- The RHS terms (main part)
  rhs(8) =   rhs(8) 											&
             ! --- Time derivative
	     + zeta * v * r0 * R							* xjac *delta_g(8)&
             ! --- Convective terms
	     + v * r0 * R**2 * ( Te0_x * u0_y - Te0_y * u0_x)				* xjac * tstep	&
	     + v * r0 * (GAMMA-1.d0) * Te0 * u0_y * 2.d0 * R				* xjac * tstep	&
	     - v * r0 *                F0 / R * Vpar0 * Te0_p				* xjac * tstep	&
	     - v * r0 * (GAMMA-1.d0) * F0 / R * Te0   * vpar0_p				* xjac * tstep	&
	     - v * r0 *              Vpar0 * (Te0_x   * ps0_y - Te0_y   * ps0_x)	* xjac * tstep	&
	     - v * r0 * (GAMMA-1.d0) * Te0 * (vpar0_x * ps0_y - vpar0_y * ps0_x)	* xjac * tstep	&
             ! --- Source
             + v * R * heat_source_e							* xjac * tstep	&
             ! --- Conductivity
	     - (Ke_par-Ke_prof) * R / BB2 * Bgrad_Te_star * Bgrad_Te			* xjac * tstep	&
	     - Ke_prof * R * (v_x*Te0_x + v_y*Te0_y )					* xjac * tstep	&
             ! --- Numerical conductivity
	     - Ke_par_num * (v_ps0_x   * ps0_y - v_ps0_y   * ps0_x)					&
			  * (Te0_ps0_x * ps0_y - Te0_ps0_y * ps0_x)			* xjac * tstep	&
             - Ke_perp_numm * (v_xx   + v_x  /R + v_yy  )						&
	                    * (Te0_xx + Te0_x/R + Te0_yy) * R 				* xjac * tstep  

  rhs_k(8) = - (Ke_par-Ke_prof) * R / BB2 * Bgrad_Te_k_star * Bgrad_Te			* xjac * tstep	&
	     - Ke_par_num * (Te0_pp  * v_pp)						* xjac * tstep	&
	     - Ke_prof * R * ( v_p*Te0_p /R**2 )					* xjac * tstep 

  
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
  Bgrad_Te_star_psi  = ( v_x   * psi_y - v_y   * psi_x ) / R
  Bgrad_Te_psi	     = ( Te0_x * psi_y - Te0_y * psi_x ) / R
  Bgrad_Te_Te	     = ( Te_x  * ps0_y - Te_y  * ps0_x ) / R
  Bgrad_Te_Te_n	     = ( F0 / R * Te_p   ) / R

  Te_ps0_x  = Te_xx  * ps0_y - Te_xy  * ps0_x + Te_x  * ps0_xy - Te_y  * ps0_xx
  Te_ps0_y  = Te_xy  * ps0_y - Te_yy  * ps0_x + Te_x  * ps0_yy - Te_y  * ps0_xy
  
  Te0_psi_x = Te0_xx * psi_y - Te0_xy * psi_x + Te0_x * psi_xy - Te0_y * psi_xx
  Te0_psi_y = Te0_xy * psi_y - Te0_yy * psi_x + Te0_x * psi_yy - Te0_y * psi_xy
  
  v_psi_x   = v_xx   * psi_y - v_xy   * psi_x + v_x   * psi_xy - v_y   * psi_xx
  v_psi_y   = v_xy   * psi_y - v_yy   * psi_x + v_x   * psi_yy - v_y   * psi_xy

  ! -----------------------------
  ! --- The LHS terms (main part)
  amat(8,1)    = + v * r0 * Vpar0 * (Te0_x * psi_y - Te0_y * psi_x)					* xjac * theta * tstep	&
		 + v * (GAMMA-1.d0) * r0 * Te0 * (vpar0_x * psi_y - vpar0_y * psi_x)			* xjac * theta * tstep	&
		 + Ke_par_num * (v_psi_x  *ps0_y - v_psi_y  *ps0_x + v_ps0_x  *psi_y - v_ps0_y  *psi_x)				&
			      * (Te0_ps0_x*ps0_y - Te0_ps0_y*ps0_x)					* xjac * theta * tstep	&
		 + Ke_par_num * (Te0_psi_x*ps0_y - Te0_psi_y*ps0_x + Te0_ps0_x*psi_y - Te0_ps0_y*psi_x)				&
			      * (v_ps0_x *ps0_y - v_ps0_y *ps0_x)					* xjac * theta * tstep	&
		 - (Ke_par-Ke_prof) * R * BB2_psi / BB2**2 * Bgrad_Te_star     * Bgrad_Te		* xjac * theta * tstep	&
		 + (Ke_par-Ke_prof) * R / BB2 	           * Bgrad_Te_star_psi * Bgrad_Te		* xjac * theta * tstep	&
		 + (Ke_par-Ke_prof) * R / BB2 	           * Bgrad_Te_star     * Bgrad_Te_psi		* xjac * theta * tstep 

  amat_k(8,1)  = - (Ke_par-Ke_prof) * R * BB2_psi / BB2**2 * Bgrad_Te_k_star   * Bgrad_Te		* xjac * theta * tstep	&
		 + (Ke_par-Ke_prof) * R / BB2 	           * Bgrad_Te_k_star   * Bgrad_Te_psi		* xjac * theta * tstep 

  amat(8,2)    = - v * r0 * R**2 * ( Te0_x * u_y - Te0_y * u_x)						* xjac * theta * tstep	&
		 - v * 2.d0 * (GAMMA-1.d0) * r0 * R * Te0 * u_y						* xjac * theta * tstep 

  amat(8,5)    = - v * rho * R**2 * (Te0_x * u0_y  - Te0_y * u0_x)					* xjac * theta * tstep	&
		 + v * rho * Vpar0   * (Te0_x * ps0_y - Te0_y * ps0_x)					* xjac * theta * tstep	&
		 + v * rho * Vpar0   * F0/R * Te0_p							* xjac * theta * tstep	&
		 - v * (GAMMA-1.d0) * rho * Te0 * u0_y * 2.d0 * R					* xjac * theta * tstep	&
		 + v * (GAMMA-1.d0) * rho * Te0 * F0/R * Vpar0_p					* xjac * theta * tstep	&
		 + v * (GAMMA-1.d0) * rho * Te0 * (Vpar0_x * ps0_y - Vpar0_y * ps0_x)			* xjac * theta * tstep 

  amat(8,7)    = + v * r0 * F0/R * Vpar * Te0_p								* xjac * theta * tstep	&
		 + v * r0 * Vpar * (Te0_x * ps0_y - Te0_y * ps0_x)					* xjac * theta * tstep	&
		 + v * (GAMMA-1.d0) * r0 * Te0 * (Vpar_x * ps0_y - Vpar_y * ps0_x)			* xjac * theta * tstep 
     
  amat_n(8,7)  = v * (GAMMA-1.d0) * r0 * Te0 * F0/R * Vpar_p  						* xjac * theta * tstep       
		   
  amat(8,8)    = - v * r0 * R**2 * (Te_x * u0_y  - Te_y * u0_x)						* xjac * theta * tstep	&
		 + v * r0 * Vpar0   * (Te_x * ps0_y - Te_y * ps0_x)					* xjac * theta * tstep	&
		 - v * r0 * (GAMMA-1.d0) * Te * R * u0_y * 2.d0						* xjac * theta * tstep	&
		 + v * r0 * (GAMMA-1.d0) * Te * F0/R * Vpar0_p						* xjac * theta * tstep	&
		 + v * r0 * (GAMMA-1.d0) * Te * (Vpar0_x * ps0_y - Vpar0_y * ps0_x)			* xjac * theta * tstep	&
		 + Ke_par_num * (v_ps0_x  * ps0_y - v_ps0_y  * ps0_x) 								&
			      * (Te_ps0_x * ps0_y - Te_ps0_y * ps0_x)					* xjac * theta * tstep	&
		 + (Ke_par-Ke_prof) * R / BB2 * Bgrad_Te_star * Bgrad_Te_Te				* xjac * theta * tstep	&
		 + dKe_par * Te     * R / BB2 * Bgrad_Te_star * Bgrad_Te				* xjac * theta * tstep	&
		 + Ke_prof * R * (v_x*Te_x + v_y*Te_y )							* xjac * theta * tstep	& 
                 + Ke_perp_numm * (v_xx  + v_x /R + v_yy )									&
		                * (Te_xx + Te_x/R + Te_yy) * R	 					* xjac * theta * tstep	
  
  if (r0 .lt. rho_1) then
    amat(8,8)  = amat(8,8) + v * rho_1 * Te * R								* xjac * (1.d0 + zeta)
  else
    amat(8,8)  = amat(8,8) + v * r0    * Te * R								* xjac * (1.d0 + zeta)
  endif
  
  amat_k(8,8)  = + (Ke_par-Ke_prof) * R / BB2 * Bgrad_Te_k_star * Bgrad_Te_Te				* xjac * theta * tstep	&
		 + dKe_par * Te     * R / BB2 * Bgrad_Te_k_star * Bgrad_Te				* xjac * theta * tstep 
	      
  amat_n(8,8)  = + (Ke_par-Ke_prof) * R / BB2 * Bgrad_Te_star   * Bgrad_Te_Te_n  			* xjac * theta * tstep	&
		 + v * r0 * Vpar0  * F0/R * Te_p							* xjac * theta * tstep 

  amat_kn(8,8) = + (Ke_par-Ke_prof) * R / BB2 * Bgrad_Te_k_star * Bgrad_Te_Te_n  			* xjac * theta * tstep	&
		 +  Ke_par_num * (Te_pp  * v_pp) 							* xjac * theta * tstep	&
		 +  Ke_prof	    * R	* (v_p*Te_p /R**2 )						* xjac * theta * tstep 


  return

end subroutine ELM_main_lhs_8









