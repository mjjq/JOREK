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
subroutine ELM_main_rhs_2_numm(rhs,rhs_k)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_rhs_2_numm

  ! --- Modules
  use phys_module
  use equation_variables
  
  implicit none
  
  ! --- Routine variables
  real*8 :: rhs(n_var),rhs_k(n_var)
  
  ! ------------------------------------------------------
  ! --- The RHS term (Taylor Galerkin (TG2) stabilisation)	      
  rhs(2) = rhs(2)												&
           - TG_num2 * tstep * 0.25d0 * r0_hat * BigR**3 * (w0_x * u0_y - w0_y * u0_x)  			&
                                                         * (v_x  * u0_y - v_y  * u0_x)		* xjac * tstep
  return

end subroutine ELM_main_rhs_2_numm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! LHS !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ELM_main_lhs_2_numm(amat, amat_k, amat_n, amat_kn)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_lhs_2_numm

  ! --- Modules
  use phys_module
  use equation_variables
  
  implicit none
  
  ! --- Routine variables
  real*8 :: amat(n_var,n_var), amat_k(n_var,n_var), amat_n(n_var,n_var), amat_kn(n_var,n_var)
  
  ! --- Internal variables
  real*8 :: rho_hat
  
  rho_hat = BigR**2 * rho
  
  ! ------------------------------------------------------
  ! --- The LHS term (Taylor Galerkin (TG2) stabilisation)	      
  amat(2,2)   = amat(2,2)														&
                + TG_num2 * tstep * 0.25d0 * r0_hat * BigR**3 * (w0_x * u_y  - w0_y * u_x )                        			&
                                                              * (v_x  * u0_y - v_y  * u0_x)			* xjac * theta * tstep	&
                + TG_num2 * tstep * 0.25d0 * r0_hat * BigR**3 * (w0_x * u0_y - w0_y * u0_x)						&
                                                              * (v_x  * u_y  - v_y  * u_x )			* xjac * theta * tstep

  amat(2,4)   = amat(2,4)														&
                + TG_num2 * tstep * 0.25d0 * r0_hat * BigR**3 * (w_x * u0_y - w_y * u0_x)						&
                                                              * (v_x * u0_y - v_y * u0_x)			* xjac * theta * tstep

  amat(2,5)   = amat(2,5)														&
                + TG_num2 * tstep * 0.25d0 * rho_hat * BigR**3 * (w0_x * u0_y - w0_y * u0_x)						&
                                                               * (v_x  * u0_y - v_y  * u0_x) 			* xjac * theta * tstep
  
  return

end subroutine ELM_main_lhs_2_numm





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
subroutine ELM_main_rhs_5_numm(rhs,rhs_k)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_rhs_5_numm

  ! --- Modules
  use phys_module
  use equation_variables
  
  implicit none
  
  ! --- Routine variables
  real*8 :: rhs(n_var),rhs_k(n_var)
  
  ! ------------------------------------------------------
  ! --- The RHS term (Taylor Galerkin (TG2) stabilisation)	      
  rhs(5)   = rhs(5)										&
             - TG_num5 * tstep * 0.25d0 * BigR**3 * (r0_x * u0_y - r0_y * u0_x)			&
                                                  * (v_x  * u0_y - v_y  * u0_x)	* xjac * tstep	&
             - TG_num5 * tstep * 0.25d0 / BigR * vpar0**2					&
                       * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)			&
                       * ( v_x * ps0_y -  v_y * ps0_x )				* xjac * tstep
  
  rhs_k(5) = rhs_k(5)										&
             - TG_num5 * tstep * 0.25d0 / BigR * vpar0**2					&
                       * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)			&
                       * (                              F0 / BigR * v_p )	* xjac * tstep
  
  
  return

end subroutine ELM_main_rhs_5_numm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! LHS !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ELM_main_lhs_5_numm(amat, amat_k, amat_n, amat_kn)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_lhs_5_numm

  ! --- Modules
  use phys_module
  use equation_variables
  
  implicit none
  
  ! --- Routine variables
  real*8 :: amat(n_var,n_var), amat_k(n_var,n_var), amat_n(n_var,n_var), amat_kn(n_var,n_var)
  
  ! ------------------------------------------------------
  ! --- The LHS term (Taylor Galerkin (TG2) stabilisation)	      
  amat(5,1)    = amat(5,1)												&
                 + TG_num5 * tstep * 0.25d0 / BigR * vpar0**2								&
                           * (r0_x * psi_y - r0_y * psi_x)								&
                           * ( v_x * ps0_y -  v_y * ps0_x )					* xjac * theta * tstep	&
                 + TG_num5 * tstep * 0.25d0 / BigR * vpar0**2								&
                           * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)						&
                           * ( v_x * psi_y -  v_y * psi_x )					* xjac * theta * tstep

  amat_k(5,1)  = amat_k(5,1)												&
                 + TG_num5 * tstep * 0.25d0 / BigR * vpar0**2								&
                           * (r0_x * psi_y - r0_y * psi_x)								&
                           * (                            + F0 / BigR * v_p)			* xjac * theta * tstep

  amat(5,2)    = amat(5,2)												&
                 + TG_num5 * tstep * 0.25d0 * BigR**3 * (r0_x * u_y  - r0_y * u_x )					&
                                                      * (v_x  * u0_y - v_y  * u0_x)		* xjac * theta * tstep	&
                 + TG_num5 * tstep * 0.25d0 * BigR**3 * (r0_x * u0_y - r0_y * u0_x)					&
                                                      * (v_x  * u_y  - v_y  * u_x )		* xjac * theta * tstep

  amat(5,5)    = amat(5,5)												&
                 + TG_num5 * tstep * 0.25d0 * BigR**3 * (rho_x * u0_y - rho_y * u0_x)					&
                                                      * ( v_x  * u0_y - v_y   * u0_x )		* xjac * theta * tstep	&
                 + TG_num5 * tstep * 0.25d0 / BigR * vpar0**2								&
                           * (rho_x * ps0_y - rho_y * ps0_x )								&
                           * ( v_x * ps0_y -  v_y * ps0_x )					* xjac * theta * tstep

  amat_k(5,5)  = amat_k(5,5)												&
                 + TG_num5 * tstep * 0.25d0 / BigR * vpar0**2								&
                           * (rho_x * ps0_y - rho_y * ps0_x )								&
                           * (                              + F0 / BigR * v_p)			* xjac * theta * tstep

  amat_n(5,5)  = amat_n(5,5)												&
                 + TG_num5 * tstep * 0.25d0 / BigR * vpar0**2								&
                           * (                              + F0 / BigR * rho_p)					&
                           * ( v_x * ps0_y -  v_y * ps0_x )					* xjac * theta * tstep

  amat_kn(5,5) = amat_kn(5,5)												&
                 + TG_num5 * tstep * 0.25d0 / BigR * vpar0**2								&
                           * (                              + F0 / BigR * rho_p)					&
                           * (                              + F0 / BigR * v_p)			* xjac * theta * tstep

  amat(5,7)    = amat(5,7)												&
                 + TG_num5 * tstep * 0.25d0 / BigR * 2.d0*vpar0*vpar							&
                           * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)						&
                           * ( v_x * ps0_y -  v_y * ps0_x                   )			* xjac * theta * tstep 

  amat_k(5,7)  = amat_k(5,7)												&
                 + TG_num5 * tstep * 0.25d0 / BigR * 2.d0*vpar0*vpar							&
                           * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)						&
                           * (                            + F0 / BigR * v_p)			* xjac * theta * tstep 

  
  return

end subroutine ELM_main_lhs_5_numm







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
subroutine ELM_main_rhs_6_numm(rhs,rhs_k)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_rhs_6_numm

  ! --- Modules
  use phys_module
  use equation_variables
  
  implicit none
  
  ! --- Routine variables
  real*8 :: rhs(n_var),rhs_k(n_var)
  
  ! ------------------------------------------------------
  ! --- The RHS term (Taylor Galerkin (TG2) stabilisation)	      
  rhs(6)   = rhs(6)											&
             - TG_num6 * tstep * 0.25d0 * BigR**3 * Ti0 * (r0_x * u0_y  - r0_y  * u0_x)			&
                                                        * (v_x  * u0_y  - v_y   * u0_x)	* xjac * tstep	&
             - TG_num6 * tstep * 0.25d0 * BigR**3 * r0  * (Ti0_x * u0_y - Ti0_y * u0_x)			&
                                                        * (v_x   * u0_y - v_y   * u0_x)	* xjac * tstep	&
             - TG_num6 * tstep * 0.25d0 / BigR * vpar0**2 * Ti0						&
                       * (r0_x  * ps0_y - r0_y  * ps0_x + F0 / BigR * r0_p)				&
                       * (v_x   * ps0_y - v_y   * ps0_x)					* xjac * tstep	&
             - TG_num6 * tstep * 0.25d0 / BigR * vpar0**2 * r0						&
                       * (Ti0_x * ps0_y - Ti0_y * ps0_x + F0 / BigR * Ti0_p)				&
                       * (v_x   * ps0_y - v_y   * ps0_x)				* xjac * tstep
  
  rhs_k(6) = rhs_k(6)											&
             - TG_num6 * tstep * 0.25d0 / BigR * vpar0**2 * Ti0						&
                       * (r0_x  * ps0_y - r0_y  * ps0_x + F0 / BigR * r0_p)				&
                       * (                                F0 / BigR * v_p)		* xjac * tstep	&
             - TG_num6 * tstep * 0.25d0 / BigR * vpar0**2 * r0						&
                       * (Ti0_x * ps0_y - Ti0_y * ps0_x + F0 / BigR * Ti0_p)				&
                       * (                                F0 / BigR * v_p)		* xjac * tstep
  
  return

end subroutine ELM_main_rhs_6_numm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! LHS !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ELM_main_lhs_6_numm(amat, amat_k, amat_n, amat_kn)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_lhs_6_numm

  ! --- Modules
  use phys_module
  use equation_variables
  
  implicit none
  
  ! --- Routine variables
  real*8 :: amat(n_var,n_var), amat_k(n_var,n_var), amat_n(n_var,n_var), amat_kn(n_var,n_var)
  
  ! ------------------------------------------------------
  ! --- The RHS term (Taylor Galerkin (TG2) stabilisation)	      
  amat(6,1)    = amat(6,1)													&
             	 + TG_num6 * tstep * 0.25d0 / BigR * vpar0**2 * Ti0								&
             		   * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)							&
             		   * ( v_x * psi_y -  v_y * psi_x )						* xjac * theta * tstep	&
             	 + TG_num6 * tstep * 0.25d0 / BigR * vpar0**2 * r0								&
             		   * (Ti0_x * ps0_y - Ti0_y * ps0_x + F0 / BigR * Ti0_p)						&
             		   * ( v_x * psi_y -  v_y * ps0_x )						* xjac * theta * tstep	&
             	 + TG_num6 * tstep * 0.25d0 / BigR * vpar0**2 * Ti0								&
             		   * (r0_x * psi_y - r0_y * psi_x)									&
             		   * ( v_x * ps0_y -  v_y * ps0_x )						* xjac * theta * tstep	&
             	 + TG_num6 * tstep * 0.25d0 / BigR * vpar0**2									&
             		   * r0 * (Ti0_x * psi_y - Ti0_y * psi_x)								&
             		   * ( v_x * ps0_y -  v_y * ps0_x )						* xjac * theta * tstep		  

  amat_k(6,1)  = amat_k(6,1)													&
                 + TG_num6 * tstep * 0.25d0 / BigR * vpar0**2 * Ti0								&
                	   * (r0_x * psi_y - r0_y * psi_x)									&
                	   * ( F0 / BigR * v_p)								* xjac * theta * tstep	&
                 + TG_num6 * tstep * 0.25d0 / BigR * vpar0**2 * r0								&
                	   * (Ti0_x * psi_y - Ti0_y * psi_x)									&
                	   * ( F0 / BigR * v_p)								* xjac * theta * tstep	     

  amat(6,2)    = amat(6,2)													&
                 + TG_num6 * tstep * 0.25d0 * BigR**2 * Ti0 * (r0_x  * u_y  - r0_y  * u_x)					&
        	                                            * ( v_x  * u0_y - v_y   * u0_x)		* xjac * theta * tstep	&
                 + TG_num6 * tstep * 0.25d0 * BigR**2 * r0  * (Ti0_x * u_y  - Ti0_y * u_x)					&
        	                                            * ( v_x  * u0_y - v_y   * u0_x)		* xjac * theta * tstep	&
                 + TG_num6 * tstep * 0.25d0 * BigR**2 * Ti0 * (r0_x  * u0_y - r0_y  * u0_x)					&
        	                                            * ( v_x  * u_y  - v_y   * u_x)		* xjac * theta * tstep	&
                 + TG_num6 * tstep * 0.25d0 * BigR**2 * r0  * (Ti0_x * u0_y - Ti0_y * u0_x)					&
        	                                            * ( v_x  * u_y  - v_y   * u_x)		* xjac * theta * tstep

  amat(6,5)    = amat(6,5)													&
             	 + TG_num6 * tstep * 0.25d0 * BigR**2 * Ti0 * (rho_x * u0_y - rho_y * u0_x)					&
             		                                    * ( v_x * u0_y  - v_y   * u0_x)		* xjac * theta * tstep	&
             	 + TG_num6 * tstep * 0.25d0 * BigR**2 * rho * (Ti0_x * u0_y - Ti0_y * u0_x)					&
             		                                    * ( v_x * u0_y  - v_y   * u0_x)		* xjac * theta * tstep	&
             	 + TG_num6 * tstep * 0.25d0 / BigR * vpar0**2 * Ti0								&
             		   * (rho_x * ps0_y - rho_y * ps0_x + F0 / BigR * rho_p)						&
             		   * ( v_x  * ps0_y -   v_y * ps0_x )						* xjac * theta * tstep	&
             	 + TG_num6 * tstep * 0.25d0 / BigR * vpar0**2 * rho								&
             		   * (Ti0_x * ps0_y - Ti0_y * ps0_x + F0 / BigR * Ti0_p)						&
             		   * (  v_x * ps0_y -   v_y * ps0_x )						* xjac * theta * tstep

  amat_k(6,5)  = amat_k(6,5)													&
                 + TG_num6 * tstep * 0.25d0 / BigR * vpar0**2 * Ti0                                                             &
                           * (rho_x * ps0_y - rho_y * ps0_x + F0 / BigR * rho_p)						&
                           * ( F0 / BigR * v_p)								* xjac * theta * tstep	&
                 + TG_num6 * tstep * 0.25d0 / BigR * vpar0**2 * rho								&
                           * (Ti0_x * ps0_y - Ti0_y * ps0_x + F0 / BigR * Ti0_p)						&
                           * ( F0 / BigR * v_p)								* xjac * theta * tstep  

  amat(6,6)    = amat(6,6)													&
             	 + TG_num6 * tstep * 0.25d0 * BigR**2 * Ti * (r0_x * u0_y - r0_y * u0_x)					&
             		                                   * ( v_x * u0_y - v_y * u0_x)			* xjac * theta * tstep	&
             	 + TG_num6 * tstep * 0.25d0 * BigR**2 * r0 * (Ti_x * u0_y - Ti_y * u0_x)					&
             		                                   * ( v_x * u0_y - v_y * u0_x)			* xjac * theta * tstep	&		 
             	 + TG_num6 * tstep * 0.25d0 / BigR * vpar0**2 * Ti								&
             		   * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)							&
             		   * ( v_x * ps0_y -  v_y * ps0_x )						* xjac * theta * tstep	&
             	 + TG_num6 * tstep * 0.25d0 / BigR * vpar0**2 * r0								&
             		   * (Ti_x * ps0_y - Ti_y * ps0_x )									&
             		   * ( v_x * ps0_y -  v_y * ps0_x )						* xjac * theta * tstep

  amat_k(6,6)  = amat_k(6,6)													&
                 + TG_num6 * tstep * 0.25d0 / BigR * vpar0**2 * Ti								&
                	   * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)							&
                	   * ( F0 / BigR * v_p)								* xjac * theta * tstep	&
                 + TG_num6 * tstep * 0.25d0 / BigR * vpar0**2 * r0								&
                	   * (Ti_x * ps0_y - Ti_y * ps0_x )									&
                	   * ( F0 / BigR * v_p)								* xjac * theta * tstep
  	      
  amat_n(6,6)  = amat_n(6,6)													&
                 + TG_num6 * tstep * 0.25d0 / BigR * vpar0**2 * r0								&
                           * ( F0 / BigR * Ti_p)										&
                           * ( v_x * ps0_y -  v_y * ps0_x )						* xjac * theta * tstep

  amat_kn(6,6) = amat_kn(6,6)													&
                 + TG_num6 * tstep * 0.25d0 / BigR * vpar0**2 * r0								&
                           * ( F0 / BigR * Ti_p)										&
                           * ( F0 / BigR * v_p)								* xjac * theta * tstep  

  amat(6,7)    = amat(6,7)													&
             	 + TG_num6 * tstep * 0.25d0 / BigR * 2.d0 * vpar0 * vpar * Ti0							&
             		   * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)							&
             		   * ( v_x * ps0_y -  v_y * ps0_x )						* xjac * theta * tstep	&
             	 + TG_num6 * tstep * 0.25d0 / BigR * 2.d0 * vpar0 * vpar * r0							&
             		   * (Ti0_x * ps0_y - Ti0_y * ps0_x + F0 / BigR * Ti0_p)						&
             		   * ( v_x * ps0_y -  v_y * ps0_x )						* xjac * theta * tstep
     
  amat_k(6,7)  = amat_k(6,7)													&
                 + TG_num6 * tstep * 0.25d0 / BigR * 2.d0 * vpar0 * vpar * Ti0							&
                           * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)							&
                           * ( F0 / BigR * v_p)								* xjac * theta * tstep	&
                 + TG_num6 * tstep * 0.25d0 / BigR * 2.d0 * vpar0 * vpar * r0							&
                           * (Ti0_x * ps0_y - Ti0_y * ps0_x + F0 / BigR * Ti0_p)						&
                           * ( F0 / BigR * v_p)								* xjac * theta * tstep
  
  return

end subroutine ELM_main_lhs_6_numm








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
subroutine ELM_main_rhs_7_numm(rhs,rhs_k)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_rhs_7_numm

  ! --- Modules
  use phys_module
  use equation_variables
  
  implicit none
  
  ! --- Routine variables
  real*8 :: rhs(n_var),rhs_k(n_var)
  
  ! ------------------------------------------------------
  ! --- The RHS term (Taylor Galerkin (TG2) stabilisation)	      
  rhs(7)   = rhs(7)												&
             - TG_NUM7 * tstep * 0.25d0 * r0 * Vpar0**2 * BB2 / BigR						&
                       * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p)			&
                       * (-(ps0_s * v_t     - ps0_t * v_s)    /xjac + F0 / BigR * v_p    )	* xjac * tstep	&
             - TG_NUM7 * tstep * 0.25d0 * v  * Vpar0**2 * BB2 / BigR						&
                       * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p)			&
                       * (-(ps0_s * r0_t    - ps0_t * r0_s)   /xjac + F0 / BigR * r0_p   )	* xjac * tstep
  
  return

end subroutine ELM_main_rhs_7_numm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! LHS !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ELM_main_lhs_7_numm(amat, amat_k, amat_n, amat_kn)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_lhs_7_numm

  ! --- Modules
  use phys_module
  use equation_variables
  
  implicit none
  
  ! --- Routine variables
  real*8 :: amat(n_var,n_var), amat_k(n_var,n_var), amat_n(n_var,n_var), amat_kn(n_var,n_var)
  
  ! ------------------------------------------------------
  ! --- The RHS term (Taylor Galerkin (TG2) stabilisation)	      
  amat(7,1)    = amat(7,1)														&
                 + TG_NUM7 * tstep * 0.25d0 * r0 * Vpar0**2 * BB2									&
                 	   * (-(psi_s * vpar0_t - psi_t * vpar0_s)/xjac) / BigR								&
                 	   * (-(psi_s * v_t	- psi_t * v_s)    /xjac)					* xjac * theta * tstep	&
                 + TG_NUM7 * tstep * 0.25d0 * v  * Vpar0**2 * BB2									&
                 	   * (-(psi_s * vpar0_t - psi_t * vpar0_s)/xjac) / BigR								&
                 	   * (-(psi_s * r0_t	- psi_t * r0_s)   /xjac)					* xjac * theta * tstep 

  amat(7,5)    = amat(7,5)														&
                 + TG_NUM7 * tstep * 0.25d0 * rho * Vpar0**2 * BB2									&
                 	   * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR					&
                 	   * (-(ps0_s * v_t	- ps0_t * v_s)    /xjac )					* xjac * theta * tstep	&
                 + TG_NUM7 * tstep * 0.25d0 * v * Vpar0**2 * BB2									&
                 	   * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR					&
                 	   * (-(ps0_s * rho_t	- ps0_t * rho_s)  /xjac )					* xjac * theta * tstep

  amat_k(7,5)  = amat_k(7,5)														&
                 + TG_NUM7 * tstep * 0.25d0 * rho * Vpar0**2 * BB2									&
                 	   * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR					&
                 	   * ( F0 / BigR * v_p )								* xjac * theta * tstep

  amat_n(7,5)  = amat_n(7,5)														&
                 + TG_NUM7 * tstep * 0.25d0 * v * Vpar0**2 * BB2									&
                 	   * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR					&
                 	   * ( F0 / BigR * rho_p )								* xjac * theta * tstep

  amat(7,7)    = amat(7,7)														&
            	 + TG_NUM7 * tstep * 0.5d0 * r0 * Vpar * Vpar0 * BB2									&
            	 	   * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR					&
            	 	   * (-(ps0_s * v_t	- ps0_t * v_s)    /xjac )					* xjac * theta * tstep	&
            	 + TG_NUM7 * tstep * 0.5d0 * v * Vpar * Vpar0 * BB2									&
            	 	   * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR					&
            	 	   * (-(ps0_s * r0_t	- ps0_t * r0_s)   /xjac + F0 / BigR * r0_p)			* xjac * theta * tstep	&
            	 + TG_NUM7 * tstep * 0.25d0 * r0 * Vpar0**2 * BB2									&
            	 	   * (-(ps0_s * vpar_t - ps0_t * vpar_s)/xjac ) / BigR								&
            	 	   * (-(ps0_s * v_t    - ps0_t * v_s)  /xjac )						* xjac * theta * tstep	&
            	 + TG_NUM7 * tstep * 0.25d0 * v * Vpar0**2 * BB2									&
            	 	   * (-(ps0_s * vpar_t - ps0_t * vpar_s)/xjac ) / BigR								&
            	 	   * (-(ps0_s * r0_t   - ps0_t * r0_s) /xjac + F0 / BigR * r0_p)			* xjac * theta * tstep

  amat_k(7,7)  = amat_k(7,7)														&
                 + TG_NUM7 * tstep * 0.5d0 * r0 * Vpar * Vpar0 * BB2									&
            	 	   * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR					&
            	 	   * ( F0 / BigR * v_p )								* xjac * theta * tstep	&
            	 + TG_NUM7 * tstep * 0.25d0 * r0 * Vpar0**2 * BB2									&
            	 	   * (-(ps0_s * vpar_t - ps0_t * vpar_s)/xjac ) / BigR								&
            	 	   * ( F0 / BigR * v_p)									* xjac * theta * tstep 

  amat_n(7,7)  = amat_n(7,7)														&
                 + TG_NUM7 * tstep * 0.25d0 * r0 * Vpar0**2 * BB2									&
            	 	   * ( F0 / BigR * vpar_p ) / BigR										&
            	 	   * (-(ps0_s * v_t - ps0_t * v_s)  /xjac )						* xjac * theta * tstep	&
            	 + TG_NUM7 * tstep * 0.25d0 * v * Vpar0**2 * BB2									&
            	 	   * ( F0 / BigR * vpar_p ) / BigR										&
            	 	   * (-(ps0_s * r0_t - ps0_t * r0_s)/xjac + F0 / BigR * r0_p)				* xjac * theta * tstep

  amat_kn(7,7) = amat_kn(7,7)														&
                 + TG_NUM7 * tstep * 0.25d0 * r0 * Vpar0**2 * BB2									&
            	 	   * ( F0 / BigR * vpar_p) / BigR										&
            	 	   * ( F0 / BigR * v_p )								* xjac * theta * tstep    

  
  return

end subroutine ELM_main_lhs_7_numm







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
subroutine ELM_main_rhs_8_numm(rhs,rhs_k)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_rhs_8_numm

  ! --- Modules
  use phys_module
  use equation_variables
  
  implicit none
  
  ! --- Routine variables
  real*8 :: rhs(n_var),rhs_k(n_var)
  
  ! ------------------------------------------------------
  ! --- The RHS term (Taylor Galerkin (TG2) stabilisation)	      
  rhs(8)   = rhs(8)											&
             - TG_num8 * tstep * 0.25d0 * BigR**3 * Te0 * ( r0_x * u0_y -  r0_y * u0_x)			&
             	                                        * (  v_x * u0_y -   v_y * u0_x)	* xjac * tstep	&
             - TG_num8 * tstep * 0.25d0 * BigR**3 * r0  * (Te0_x * u0_y - Te0_y * u0_x)			&
             	                                        * (  v_x * u0_y -   v_y * u0_x)	* xjac * tstep	&
             - TG_num8 * tstep * 0.25d0 / BigR * vpar0**2 * Te0						&
             	       * (r0_x  * ps0_y - r0_y  * ps0_x + F0 / BigR * r0_p)				&
             	       * ( v_x  * ps0_y -  v_y  * ps0_x)				* xjac * tstep	&
             - TG_num8 * tstep * 0.25d0 / BigR * vpar0**2 * r0						&
             	       * (Te0_x * ps0_y - Te0_y * ps0_x + F0 / BigR * Te0_p)				&
             	       * (  v_x * ps0_y -   v_y * ps0_x)				* xjac * tstep 
  
  rhs_k(8) = rhs_k(8)											&
             - TG_num8 * tstep * 0.25d0 / BigR * vpar0**2 * Te0						&
             	       * (r0_x  * ps0_y - r0_y  * ps0_x + F0 / BigR * r0_p)				&
             	       * (                                F0 / BigR * v_p)		* xjac * tstep	&
             - TG_num8 * tstep * 0.25d0 / BigR * vpar0**2 * r0						&
             	       * (Te0_x * ps0_y - Te0_y * ps0_x + F0 / BigR * Te0_p)				&
             	       * (                                F0 / BigR * v_p)		* xjac * tstep 
  
  return

end subroutine ELM_main_rhs_8_numm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! LHS !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ELM_main_lhs_8_numm(amat, amat_k, amat_n, amat_kn)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_main_lhs_8_numm

  ! --- Modules
  use phys_module
  use equation_variables
  
  implicit none
  
  ! --- Routine variables
  real*8 :: amat(n_var,n_var), amat_k(n_var,n_var), amat_n(n_var,n_var), amat_kn(n_var,n_var)
  
  ! ------------------------------------------------------
  ! --- The RHS term (Taylor Galerkin (TG2) stabilisation)	      
  amat(8,1)    = amat(8,1)													&
                 + TG_num8 * tstep * 0.25d0 / BigR * vpar0**2 * Te0								&
                     	   * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)							&
                     	   * ( v_x * psi_y -  v_y * psi_x )						* xjac * theta * tstep	&
                 + TG_num8 * tstep * 0.25d0 / BigR * vpar0**2 * r0								&
                     	   * (Te0_x * ps0_y - Te0_y * ps0_x + F0 / BigR * Te0_p)						&
                     	   * (  v_x * psi_y -   v_y * ps0_x )						* xjac * theta * tstep	&
                 + TG_num8 * tstep * 0.25d0 / BigR * vpar0**2 * Te0								&
                     	   * (r0_x * psi_y - r0_y * psi_x)									&
                     	   * ( v_x * ps0_y -  v_y * ps0_x)						* xjac * theta * tstep	&
                 + TG_num8 * tstep * 0.25d0 / BigR * vpar0**2 * r0								&
                     	   * (Te0_x * psi_y - Te0_y * psi_x)									&
                     	   * (  v_x * ps0_y -   v_y * ps0_x)						* xjac * theta * tstep

  amat_k(8,1)  = amat_k(8,1)													&
                 + TG_num8 * tstep * 0.25d0 / BigR * vpar0**2 * Te0								&
                   	   * (r0_x * psi_y - r0_y * psi_x)									&
                   	   * ( F0 / BigR * v_p)								* xjac * theta * tstep	&
                 + TG_num8 * tstep * 0.25d0 / BigR * vpar0**2 * r0								&
                   	   * (Te0_x * psi_y - Te0_y * psi_x)									&
                   	   * ( F0 / BigR * v_p)								* xjac * theta * tstep    

  amat(8,2)    = amat(8,2)													&
                 + TG_num8 * tstep * 0.25d0 * BigR**2 * Te0 * ( r0_x * u_y  -  r0_y * u_x)					&
                     	                                    * (  v_x * u0_y -   v_y * u0_x)		* xjac * theta * tstep	&
                 + TG_num8 * tstep * 0.25d0 * BigR**2 * r0  * (Te0_x * u_y  - Te0_y * u_x)					&
                     	   * ( v_x * u0_y - v_y * u0_x)							* xjac * theta * tstep	&
                 + TG_num8 * tstep * 0.25d0 * BigR**2 * Te0 * ( r0_x * u0_y -  r0_y * u0_x)					&
                     	                                    * (  v_x * u_y  -   v_y * u_x)		* xjac * theta * tstep	&
                 + TG_num8 * tstep * 0.25d0 * BigR**2 * r0  * (Te0_x * u0_y - Te0_y * u0_x)					&
                     	                                    * (  v_x * u_y  -   v_y * u_x)		* xjac * theta * tstep
 	  	  	    
  amat(8,5)    = amat(8,5)													&
                 + TG_num8 * tstep * 0.25d0 * BigR**2 * Te0 * (rho_x * u0_y - rho_y * u0_x)					&
                     	                                    * ( v_x * u0_y - v_y * u0_x)		* xjac * theta * tstep	&
                 + TG_num8 * tstep * 0.25d0 * BigR**2 * rho * (Te0_x * u0_y - Te0_y * u0_x)					&
                     	                                    * ( v_x * u0_y - v_y * u0_x)		* xjac * theta * tstep	&
                 + TG_num8 * tstep * 0.25d0 / BigR * vpar0**2 * Te0								&
                 	   * (rho_x * ps0_y - rho_y * ps0_x + F0 / BigR * rho_p)						&
                 	   * ( v_x * ps0_y -  v_y * ps0_x )						* xjac * theta * tstep	&
                 + TG_num8 * tstep * 0.25d0 / BigR * vpar0**2 * rho								&
                 	   * (Te0_x * ps0_y - Te0_y * ps0_x + F0 / BigR * Te0_p)						&
                 	   * ( v_x * ps0_y -  v_y * ps0_x )						* xjac * theta * tstep  

  amat_k(8,5)  = amat_k(8,5)													&
  	         + TG_num8 * tstep * 0.25d0 / BigR * vpar0**2 * Te0								&
                 	   * (rho_x * ps0_y - rho_y * ps0_x + F0 / BigR * rho_p)						&
                 	   * ( F0 / BigR * v_p)								* xjac * theta * tstep	&
                 + TG_num8 * tstep * 0.25d0 / BigR * vpar0**2 * rho								&
                 	   * (Te0_x * ps0_y - Te0_y * ps0_x + F0 / BigR * Te0_p)						&
                 	   * ( F0 / BigR * v_p)								* xjac * theta * tstep  
			  
  amat(8,7)    = amat(8,7)													&
                 + TG_num8 * tstep * 0.25d0 / BigR * 2.d0 * vpar0 * vpar * Te0							&
                 	   * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)							&
                 	   * ( v_x * ps0_y -  v_y * ps0_x )						* xjac * theta * tstep	&
                 + TG_num8 * tstep * 0.25d0 / BigR * 2.d0 * vpar0 * vpar * r0							&
                 	   * (Te0_x * ps0_y - Te0_y * ps0_x + F0 / BigR * Te0_p)						&
                 	   * ( v_x * ps0_y -  v_y * ps0_x )						* xjac * theta * tstep
 	   
  amat_k(8,7)  = amat_k(8,7)													&
 	    	 + TG_num8 * tstep * 0.25d0 / BigR * 2.d0 * vpar0 * vpar * Te0							&
                 	   * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)							&
                 	   * ( F0 / BigR * v_p)								* xjac * theta * tstep	&
                 + TG_num8 * tstep * 0.25d0 / BigR * 2.d0 * vpar0 * vpar * r0							&
                 	   * (Te0_x * ps0_y - Te0_y * ps0_x + F0 / BigR * Te0_p)						&
                 	   * ( F0 / BigR * v_p)								* xjac * theta * tstep
 	  	  				    
  amat(8,8)    = amat(8,8)													&
                 + TG_num8 * tstep * 0.25d0 * BigR**2 * Te * (r0_x * u0_y - r0_y * u0_x)					&
                     	   * ( v_x * u0_y - v_y * u0_x)							* xjac * theta * tstep	&
                 + TG_num8 * tstep * 0.25d0 * BigR**2 * r0* (Te_x * u0_y - Te_y * u0_x)						&
                     	   * ( v_x * u0_y - v_y * u0_x)							* xjac * theta * tstep	&
                 + TG_num8 * tstep * 0.25d0 / BigR * vpar0**2 * Te								&
                     	   * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)							&
                     	   * ( v_x * ps0_y -  v_y * ps0_x )						* xjac * theta * tstep	&
                 + TG_num8 * tstep * 0.25d0 / BigR * vpar0**2 * r0								&
                     	   * (Te_x * ps0_y - Te_y * ps0_x )									&
                     	   * ( v_x * ps0_y -  v_y * ps0_x )						* xjac * theta * tstep 

  amat_k(8,8)  = amat_k(8,8)													&
                 + TG_num8 * tstep * 0.25d0 / BigR * vpar0**2 * Te								&
                   	   * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)							&
                   	   * ( F0 / BigR * v_p)								* xjac * theta * tstep	&
                 + TG_num8 * tstep * 0.25d0 / BigR * vpar0**2 * r0								&
                   	   * (Te_x * ps0_y - Te_y * ps0_x )									&
                   	   * ( F0 / BigR * v_p)								* xjac * theta * tstep
 	    	     
  amat_n(8,8)  = amat_n(8,8)													&
                 + TG_num8 * tstep * 0.25d0 / BigR * vpar0**2 * r0								&
                 	   * ( F0 / BigR * Te_p)										&
                 	   * ( v_x * ps0_y -  v_y * ps0_x )						* xjac * theta * tstep
 	  	  		     
  amat_kn(8,8) = amat_kn(8,8)													&
                 + TG_num8 * tstep * 0.25d0 / BigR * vpar0**2 * r0								&
                 	   * ( F0 / BigR * Te_p)										&
                 	   * ( F0 / BigR * v_p)								* xjac * theta * tstep  

  
  return

end subroutine ELM_main_lhs_8_numm









