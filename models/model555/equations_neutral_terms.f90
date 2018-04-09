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
subroutine ELM_neutral_rhs_2(rhs,rhs_k)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_neutral_rhs_2

  ! --- Modules
  use phys_module
  use equation_variables
  
  implicit none
  
  ! --- Routine variables
  real*8 :: rhs(n_var),rhs_k(n_var)
  
  ! ----------------	       
  ! --- The RHS term	       
  rhs(2) = rhs(2)												&
           + R**3 * (r0*rn0*S_ion)   * (v_x * u0_x + v_y * u0_y)				* xjac * tstep

  return

end subroutine ELM_neutral_rhs_2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! LHS !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ELM_neutral_lhs_2(amat, amat_k, amat_n, amat_kn)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_neutral_lhs_2

  ! --- Modules
  use phys_module
  use equation_variables
  
  implicit none
  
  ! --- Routine variables
  real*8 :: amat(n_var,n_var), amat_k(n_var,n_var), amat_n(n_var,n_var), amat_kn(n_var,n_var)
  
  ! -----------------	   	 
  ! --- The LHS terms
  amat(2,2)   = amat(2,2)														&
                - R**3 * (r0*rn0*S_ion)   * (v_x * u_x + v_y * u_y) 						* xjac * theta * tstep 
  
  amat(2,5)   = amat(2,5)														&
                - R**3 * (rho*rn0*S_ion) * (v_x * u0_x + v_y * u0_y) 						* xjac * theta * tstep

  amat(2,6)   = amat(2,6)														&
                - R**3 * (r0*rn0*S_ion_T*T) * (v_x * u0_x + v_y * u0_y) 					* xjac * theta * tstep      

  amat(2,8)   = amat(2,8)														&
                - R**3 * (r0*rhon*S_ion) * (v_x * u0_x + v_y * u0_y) 						* xjac * theta * tstep
  
  return

end subroutine ELM_neutral_lhs_2




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
subroutine ELM_neutral_rhs_5(rhs,rhs_k)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_neutral_rhs_5

  ! --- Modules
  use phys_module
  use equation_variables
  
  implicit none
  
  ! --- Routine variables
  real*8 :: rhs(n_var),rhs_k(n_var)
  
  ! ----------------------------	      
  ! --- The RHS term (main part)	      
  rhs(5)   = rhs(5)									 	&
             + v * r0 * rn0 * R * S_ion						* xjac * tstep
  
  return

end subroutine ELM_neutral_rhs_5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! LHS !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ELM_neutral_lhs_5(amat, amat_k, amat_n, amat_kn)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_neutral_lhs_5

  ! --- Modules
  use phys_module
  use equation_variables
  
  implicit none
  
  ! --- Routine variables
  real*8 :: amat(n_var,n_var), amat_k(n_var,n_var), amat_n(n_var,n_var), amat_kn(n_var,n_var)
  
  ! -----------------------------	      
  ! --- The LHS terms (main part)
  amat(5,5)    = amat(5,5)												&
		 - v * R * rho * rn0 * S_ion							* xjac * theta * tstep

  amat(5,6)    = amat(5,6)												&
                 - v * R * r0 * rn0 * S_ion_T * T						* xjac * theta * tstep
  
  amat(5,8)    = amat(5,8)												&
                 - R * v * r0 * S_ion * rhon							* xjac * theta * tstep
  
  return

end subroutine ELM_neutral_lhs_5







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
subroutine ELM_neutral_rhs_6(rhs,rhs_k)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_neutral_rhs_6

  ! --- Modules
  use phys_module
  use equation_variables
  
  implicit none
  
  ! --- Routine variables
  real*8 :: rhs(n_var),rhs_k(n_var)
  
  ! -----------------------------	      
  ! --- The RHS terms (main part)
  rhs(6) =   rhs(6)											&
             - v * R * ksiion * r0 * rn0 * S_ion                        	     	* xjac * tstep
  
  return

end subroutine ELM_neutral_rhs_6


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! LHS !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ELM_neutral_lhs_6(amat, amat_k, amat_n, amat_kn)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_neutral_lhs_6

  ! --- Modules
  use phys_module
  use equation_variables
  
  implicit none
  
  ! --- Routine variables
  real*8 :: amat(n_var,n_var), amat_k(n_var,n_var), amat_n(n_var,n_var), amat_kn(n_var,n_var)
  
  ! -----------------------------	      
  ! --- The LHS terms (main part)
  amat(6,5)    = amat(6,5)													&
		 + v * R * rho * rn0 * ksiion * S_ion                             			* xjac * theta * tstep
  
  amat(6,6)    = amat(6,6)													&
		 + v * R * r0 * rn0 * ksiion * S_ion_T * T               	     			* xjac * theta * tstep
  
  amat(6,8)    = amat(6,8)													&
                 + v * R * r0 * rhon * ksiion * S_ion							* xjac * theta * tstep			   

  return

end subroutine ELM_neutral_lhs_6








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
subroutine ELM_neutral_rhs_7(rhs,rhs_k)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_neutral_rhs_7

  ! --- Modules
  use phys_module
  use equation_variables
  
  implicit none
  
  ! --- Routine variables
  real*8 :: rhs(n_var),rhs_k(n_var)
  
  ! -----------------------------	      
  ! --- The RHS terms (main part)
  rhs(7)   = rhs(7)											 	&
             - v *(r0 * rn0 * S_ion) * vpar0 * BB2 * R	                        		* xjac * tstep 
  
  return

end subroutine ELM_neutral_rhs_7

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! LHS !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ELM_neutral_lhs_7(amat, amat_k, amat_n, amat_kn)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_neutral_lhs_7

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
  amat(7,1)   = amat(7,1)														&
	        + v *(r0 * rn0 * S_ion) * vpar0 * BB2_psi * R	         					* xjac * theta * tstep
  
  amat(7,5)   = amat(7,5)														&
		+ v *(rho * rn0 * S_ion) * vpar0 * BB2 * R							* xjac * theta * tstep

  amat(7,6)   = amat(7,6)														&
		+ v *(r0 * rn0 * S_ion_T * T) * vpar0 * BB2 * R							* xjac * theta * tstep
  
  amat(7,7)   = amat(7,7)														&
                + v *(r0 * rn0 * S_ion) * vpar * BB2 * R							* xjac * theta * tstep

  amat(7,8)   = amat(7,8)														&
                + v *(r0 * rhon * S_ion) * vpar0 * BB2 * R							* xjac * theta * tstep

  
  return

end subroutine ELM_neutral_lhs_7







!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------ Equation 8 (Rho_n - Neutral Density) ----------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! RHS !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ELM_neutral_rhs_8(rhs,rhs_k)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_neutral_rhs_8

  ! --- Modules
  use phys_module
  use equation_variables
  
  implicit none
  
  ! --- Routine variables
  real*8 :: rhs(n_var),rhs_k(n_var)
  
  ! -----------------------------	      
  ! --- The RHS terms (main part)
  rhs(8)   = 															&
             ! --- Time derivative terms
	     + zeta * v * R											* xjac *delta_g(8)&  
             ! --- The diffusion term
             + R* (- Dn0x * rn0_x * v_x - Dn0y * rn0_y * v_y)							* xjac * tstep	&         
             ! --- Ionisation term
	     - R * v * r0 * rn0 * S_ion										* xjac * tstep	&  
             ! --- Source
	     + R * v * source_neutral										* xjac * tstep

  rhs_k(8) = + R* (Dn0p * rn0_p * v_p*eps_cyl**2/R**2)								* xjac * tstep        
  
  return

end subroutine ELM_neutral_rhs_8


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! LHS !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ELM_neutral_lhs_8(amat, amat_k, amat_n, amat_kn)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_neutral_lhs_8

  ! --- Modules
  use phys_module
  use equation_variables
  
  implicit none
  
  ! --- Routine variables
  real*8 :: amat(n_var,n_var), amat_k(n_var,n_var), amat_n(n_var,n_var), amat_kn(n_var,n_var)
  
  ! -----------------------------
  ! --- The LHS terms (main part)
  amat(8,5)    = + R * v * rn0 * S_ion * rho									* xjac * theta * tstep     
                 
  amat(8,6)    = + R * v * r0 * rn0 * S_ion_T * T								* xjac * theta * tstep      

                 
  amat(8,8)    = + v * rhon * R 										* xjac * (1.d0 + zeta)	&
		 + R * (Dn0x * rhon_x * v_x + Dn0y * rhon_y * v_y)						* xjac * theta * tstep	&   
		 + R * v * r0 * rhon* S_ion                                                                	* xjac * theta * tstep     
  
  amat_kn(8,8) = + R * (Dn0p * rhon_p * v_p*eps_cyl**2/R**2)							* xjac * theta * tstep 
  

  return

end subroutine ELM_neutral_lhs_8









