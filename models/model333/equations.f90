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
  use mod_parameters
  use phys_module
  use equation_variables
  
  implicit none
  
  ! --- Routine variables
  real*8 :: rhs(n_var),rhs_k(n_var)
  
  ! ----------------------------              
  ! --- The RHS term (main part)              
  rhs(1) =                                                                                              &
           ! --- Time derivative
           + zeta * v / R                                                               * xjac * delta_g(1)&
           ! --- Resistive term
           + v * eta_T  * (zj0 - current_source - jb)/ R                                * xjac * tstep  &
           ! --- VxB
           + v * (ps0_x * u0_y - ps0_y * u0_x)                                          * xjac * tstep  &
           ! --- Integration term
           - v * eps_cyl * F0 / R  * u0_p                                               * xjac * tstep  &
           ! --- Numerical resistivity
           + eta_numm * (v_x * zj0_x + v_y * zj0_y)                                     * xjac * tstep  
  
  ! -----------------------------------    
  ! --- The RHS term (diamagnetic part)       
  rhs(1) = rhs(1)                                                                                       &
           - v * tau_IC/(r0_corr2*BB2) * F0**2/R**2 * (ps0_x*p0_y - ps0_y*p0_x)         * xjac * tstep  &                 
           + v * tau_IC/(r0_corr2*BB2) * F0**3/R**3 * eps_cyl * p0_p                    * xjac * tstep  
  
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
  amat(1,1)   = + v * psi / R                                                                           * xjac * (1.d0+zeta)    &
                - v * (psi_x * u0_y - psi_y * u0_x)                                                     * xjac * theta * tstep  

  amat(1,2)   = -  v * (ps0_x * u_y - ps0_y * u_x)                                                      * xjac * theta * tstep

  amat_n(1,2) = +  eps_cyl * F0 / R * v * u_p                                                           * xjac * theta * tstep

  amat(1,3)   = - eta_numm * (v_x * zj_x + v_y * zj_y)                                                  * xjac * theta * tstep  &
                - eta_T * v * zj / R                                                                    * xjac * theta * tstep

  amat(1,6)   = - deta_dT * v * T * (zj0 - current_source - jb) / R                                     * xjac * theta * tstep  

  ! ------------------------------------
  ! --- The LHS terms (diamagnetic part)
  amat(1,1)   = amat(1,1)                                                                                                       &
                - v * tau_IC/(r0_corr2*BB2**2) * BB2_psi * F0**2/R**2 * (ps0_x*p0_y - ps0_y*p0_x)       * xjac * theta * tstep  &
                + v * tau_IC/(r0_corr2*BB2**2) * BB2_psi * F0**3/R**3 * eps_cyl * p0_p                  * xjac * theta * tstep  &
                + v * tau_IC/(r0_corr2*BB2) * F0**2/R**2 * (psi_x * p0_y - psi_y * p0_x)                * xjac * theta * tstep
  
  amat(1,5)   = amat(1,5)                                                                                                       &
                + v * tau_IC/(r0_corr2*BB2) * F0**2/R**2 * T0  * (ps0_x*rho_y - ps0_y*rho_x)            * xjac * theta * tstep  &
                + v * tau_IC/(r0_corr2*BB2) * F0**2/R**2 * rho * (ps0_x*T0_y  - ps0_y*T0_x )            * xjac * theta * tstep  &
                - v * tau_IC/(r0_corr2*BB2) * F0**3/R**3 * eps_cyl * rho * T0_p                         * xjac * theta * tstep  &
                - v * tau_IC * rho /(r0_corr2**2 * BB2) * F0**2/R**2 * (ps0_x*p0_y - ps0_y*p0_x)        * xjac * theta * tstep  &                   
                + v * tau_IC * rho /(r0_corr2**2 * BB2) * F0**3/R**3 * eps_cyl * p0_p                   * xjac * theta * tstep 
  
  amat_n(1,5) = amat_n(1,5)                                                                                                     &
                - v * tau_IC/(r0_corr2*BB2) * F0**3/R**3 * eps_cyl * T0 * rho_p                         * xjac * theta * tstep
  
  amat(1,6)   = amat(1,6)                                                                                                       &
                + v * tau_IC/(r0_corr2*BB2) * F0**2/R**2 * r0 * (ps0_x*T_y  - ps0_y*T_x )               * xjac * theta * tstep  &
                + v * tau_IC/(r0_corr2*BB2) * F0**2/R**2 * T  * (ps0_x*r0_y - ps0_y*r0_x)               * xjac * theta * tstep  &
                - v * tau_IC/(r0_corr2*BB2) * F0**3/R**3 * eps_cyl * T * r0_p                           * xjac * theta * tstep 

  amat_n(1,6) = amat_n(1,6)                                                                                                     &
                - v * tau_IC/(r0_corr2*BB2) * F0**3/R**3 * eps_cyl * r0 * T_p                           * xjac * theta * tstep 
  
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
  rhs(2) =                                                                                                              &
           ! --- Time derivative
           - zeta * R**3 * r0_corr * (v_x * delta_u_x + v_y * delta_u_y)                                * xjac          &                 
           ! --- Convective terms
           - 0.5d0 * vv2 * (v_x * r0_y_hat - v_y * r0_x_hat)                                            * xjac * tstep  &
           - r0_corr * R**4 * w0 * (v_x * u0_y - v_y * u0_x)                                            * xjac * tstep  &
           ! --- [psi,j]
           + v * (ps0_x * zj0_y - ps0_y * zj0_x )                                                       * xjac * tstep  &
           - v * eps_cyl * F0 / R * zj0_p                                                               * xjac * tstep  &
           ! --- Grad(p)
           + R**2 * (v_x * p0_y     - v_y * p0_x)                                                       * xjac * tstep  &
           ! --- Source interaction
           + R**3 * total_rho_source * (v_x * u0_x + v_y * u0_y)                                        * xjac * tstep  &
           ! --- Viscosity
           - visco_T * R * (v_x * w0_x + v_y * w0_y)                                                    * xjac * tstep  &
           ! --- Numerical iscosity
           - visco_numm * (v_xx  + v_x /R + v_yy )                                                                      &
                        * (w0_xx + w0_x/R + w0_yy)                                                      * xjac * tstep 
  
  ! --------------------------------------------------      
  ! --- The RHS term (diamagnetic and neoclassic part)        
  rhs(2) = rhs(2)                                                                                                       &
           ! --- Diamagnetic terms
           - tau_IC * v * R**4         * (P0_x* w0_y - P0_y* w0_x)                                      * xjac * tstep  &
           - tau_IC     * R**3 * P0_y  * (v_x * u0_x + v_y * u0_y)                                      * xjac * tstep  &
           - tau_IC * v * R**4 * u0_xy * (P0_xx-P0_yy)                                                  * xjac * tstep  &
           + tau_IC * v * R**4 * P0_xy * (u0_xx-u0_yy)                                                  * xjac * tstep  &
           ! --- Diamagnetic viscosity
           + dvisco_dT * R * W_dia * (v_x*T0_x + v_y*T0_y)                                              * xjac * tstep  &
           + visco_T   * R * W_dia * (v_xx + v_x/R + v_yy)                                              * xjac * tstep  &
           ! --- Neoclassic term
           + amu_neo_prof * BB2 / (Btheta2+epsil)**2 * (ps0_x*v_x + ps0_y*v_y) * R                                      &
                    * (  r0                         * (ps0_x*u0_x + ps0_y*u0_y)                                         &
                       + tau_IC                     * (ps0_x*P0_x + ps0_y*P0_y)                                         &
                       + aki_neo_prof * tau_IC * r0 * (ps0_x*T0_x + ps0_y*T0_y)                                         &
                       - r0 * Vpar0 * Btheta2                                        )                  * xjac * tstep 

  
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
  real*8 :: P0_x_rho, P0_xx_rho, P0_y_rho, P0_yy_rho, P0_xy_rho
  real*8 :: P0_x_T,   P0_xx_T,   P0_y_T,   P0_yy_T,   P0_xy_T
  
  rho_hat     = R**2 * rho
  rho_x_hat   = 2.d0 * R * R_x  * rho + R**2 * rho_x
  rho_y_hat   = R**2 * rho_y
  Btheta2_psi = 2.d0 * (psi_x * ps0_x + psi_y * ps0_y ) / R**2

  P0_x_rho  = rho_x*T0 + rho*T0_x
  P0_xx_rho = rho_xx*T0 + 2.0*rho_x*T0_x + rho*T0_xx
  P0_y_rho  = rho_y*T0 + rho*T0_y
  P0_yy_rho = rho_yy*T0 + 2.0*rho_y*T0_y + rho*T0_yy
  P0_xy_rho = rho_xy*T0 + rho*T0_xy + rho_x*T0_y + rho_y*T0_x
  
  P0_x_T  = r0_x*T + r0*T_x
  P0_xx_T = r0_xx*T + 2.0*r0_x*T_x + r0*T_xx
  P0_y_T  = r0_y*T + r0*T_y
  P0_yy_T = r0_yy*T + 2.0*r0_y*T_y + r0*T_yy
  P0_xy_T = r0_xy*T + r0*T_xy + r0_x*T_y + r0_y*T_x
  
  if (Wdia) then
    W_dia_rho = - tau_IC *     rho/r0_corr2**2 * (p0_xx     + p0_x    /R + p0_yy    ) &
                + tau_IC          /r0_corr2    * (p0_xx_rho + p0_x_rho/R + p0_yy_rho) &
                + tau_IC * 2.0*rho/r0_corr2**3 * (r0_x *p0_x     + r0_y *p0_y    )    &
                - tau_IC          /r0_corr2**2 * (rho_x*p0_x     + rho_y*p0_y    )    &
                - tau_IC          /r0_corr2**2 * (r0_x *p0_x_rho + r0_y *p0_y_rho)
    W_dia_T   = + tau_IC          /r0_corr2    * (p0_xx_T + p0_x_T/R + p0_yy_T) &
                - tau_IC          /r0_corr2**2 * (r0_x  *p0_x_T + r0_y  *p0_y_T)
  else
    W_dia_rho = 0.d0
    W_dia_T   = 0.d0
  endif
  
  ! -----------------------------             
  ! --- The LHS terms (main part)
  amat(2,1)   = - v * (psi_x * zj0_y - psi_y * zj0_x )                                                                 * xjac * theta * tstep

  amat(2,2)   = - R**3 * r0_corr * (v_x * u_x + v_y * u_y)                                                             * xjac * (1.d0 + zeta)  &
                + r0_hat * R**2 * w0 * (v_x * u_y  - v_y  * u_x)                                                       * xjac * theta * tstep  &
                + R**2 * (u_x*u0_x + u_y*u0_y) * (v_x*r0_y_hat - v_y*r0_x_hat)                                         * xjac * theta * tstep  &
                - R**3 * total_rho_source * (v_x * u_x + v_y * u_y)                                                    * xjac * theta * tstep
  
  amat(2,3)   = - v * (ps0_x * zj_y  - ps0_y * zj_x)                                                                   * xjac * theta * tstep

  amat_n(2,3) = + eps_cyl * F0 / R * v * zj_p                                                                          * xjac * theta * tstep

  amat(2,4)   = r0_corr * R**4 * w  * ( v_x * u0_y - v_y * u0_x)                                                       * xjac * theta * tstep  &
                + visco_T * R * ( v_x * w_x + v_y * w_y)                                                               * xjac * theta * tstep  &
                + visco_numm * (v_xx + v_x/R + v_yy)                                                                                           &
                             * (w_xx + w_x/R + w_yy)                                                                   * xjac * theta * tstep    

  amat(2,5)   = + 0.5d0 * vv2 * (v_x * rho_y_hat - v_y * rho_x_hat)                                                    * xjac * theta * tstep  &
                + rho_hat * R**2 * w0 * (v_x * u0_y - v_y * u0_x)                                                      * xjac * theta * tstep  &
                - R**2 * (v_x * rho_y * T0   - v_y * rho_x * T0  )                                                     * xjac * theta * tstep  &
                - R**2 * (v_x * rho   * T0_y - v_y * rho   * T0_x)                                                     * xjac * theta * tstep

  amat(2,6)   = - R**2 * (v_x * r0_y * T   - v_y * r0_x * T)                                                           * xjac * theta * tstep  &
                - R**2 * (v_x * r0      * T_y - v_y * r0   * T_x)                                                      * xjac * theta * tstep  &
                + dvisco_dT * T * ( v_x * w0_x    + v_y * w0_y    ) * R                                                * xjac * theta * tstep
    
  
  ! ---------------------------------------------------    
  ! --- The LHS terms (diamagnetic and neoclassic part)
  amat(2,1)   = amat(2,1)                                                                                                                      &
                ! --- Neoclassical term
                - amu_neo_prof * BB2 / (Btheta2+epsil)**2 * (psi_x*v_x+psi_y*v_y) * R                                                          &
                               * (  r0                         * (ps0_x*u0_x + ps0_y*u0_y)                                                     &
                                  + tau_IC                     * (ps0_x*P0_x + ps0_y*P0_y)                                                     &
                                  + aki_neo_prof * tau_IC * r0 * (ps0_x*T0_x + ps0_y*T0_y)                                                     &
                                  - r0 * Vpar0 * Btheta2)                                                              * xjac * theta * tstep  &
                - amu_neo_prof * BB2 / (Btheta2+epsil)**2 * (ps0_x*v_x+ps0_y*v_y) * R                                                          &
                               * (  r0                         * (psi_x*u0_x + psi_y*u0_y)                                                     &
                                  + tau_IC                     * (psi_x*P0_x + psi_y*P0_y)                                                     &
                                  + aki_neo_prof * tau_IC * r0 * (psi_x*T0_x + psi_y*T0_y)      )                      * xjac * theta * tstep  &
                + amu_neo_prof * BB2 * 2.d0*Btheta2_psi / (Btheta2+epsil)**3 * (ps0_x*v_x+ps0_y*v_y) * R                                       &
                               * (  r0                         * (ps0_x*u0_x + ps0_y*u0_y)                                                     &
                                  + tau_IC                     * (ps0_x*P0_x + ps0_y*P0_y)                                                     &
                                  + aki_neo_prof * tau_IC * r0 * (ps0_x*T0_x + ps0_y*T0_y)      )                      * xjac * theta * tstep  &
                - amu_neo_prof * BB2 * Btheta2_psi / (Btheta2+epsil)**2                                                                        &
                               * r0 * vpar0 * (ps0_x*v_x + ps0_y*v_y) * R                                              * xjac * tstep * theta

  amat(2,2)   = amat(2,2)                                                                                                                      &
                ! --- Diamagnetic terms
                + tau_IC     * R**3 * P0_y  * (v_x  * u_x + v_y  * u_y)                                                * xjac * theta * tstep  &
                + tau_IC * v * R**4 * u_xy  * (P0_xx-P0_yy)                                                            * xjac * theta * tstep  &
                - tau_IC * v * R**4 * P0_xy * (u_xx -u_yy )                                                            * xjac * theta * tstep  &
                ! --- Neoclassical term
                - amu_neo_prof * BB2 / (Btheta2+epsil)**2 * (ps0_x*v_x + ps0_y*v_y) * R                                                        &
                               * r0 * (ps0_x*u_x + ps0_y*u_y)                                                          * xjac * theta * tstep
  
  amat(2,4)   = amat(2,4)                                                                                                                      &
                ! --- Main diamagnetic terms
                + tau_IC * v * R**4         * (P0_x* w_y - P0_y* w_x)                                                  * xjac * theta * tstep  

  amat(2,5)   = amat(2,5)                                                                                                                      &
                ! --- Diamagnetic terms
                + tau_IC * v * R**4 * (P0_x_rho * w0_y - P0_y_rho * w0_x)                                              * xjac * theta * tstep  &
                + tau_IC     * R**3 * P0_y_rho  * (v_x * u0_x + v_y * u0_y)                                            * xjac * theta * tstep  &
                + tau_IC * v * R**4 * u0_xy * (P0_xx_rho-P0_yy_rho)                                                    * xjac * theta * tstep  &
                - tau_IC * v * R**4 * P0_xy_rho * (u0_xx-u0_yy)                                                        * xjac * theta * tstep  &
                ! --- Diamagnetic viscosity
                - dvisco_dT * R * W_dia_rho * (v_x*T0_x + v_y*T0_y)                                                    * xjac * theta * tstep  &
                - visco_T   * R * W_dia_rho * (v_xx + v_x/R + v_yy)                                                    * xjac * theta * tstep  &
                ! --- Neoclassical term
                - amu_neo_prof * BB2 / (Btheta2+epsil)**2 * (ps0_x*v_x + ps0_y*v_y) * R                                                        &
                               * (  rho                         * (ps0_x*u0_x       + ps0_y*u0_y      )                                        &
                                  + tau_IC                      * (ps0_x*rho_x*T0   + ps0_y*rho_y*T0  )                                        &
                                  + tau_IC                      * (ps0_x*rho  *T0_x + ps0_y*rho  *T0_y)                                        &
                                  + aki_neo_prof * tau_IC * rho * (ps0_x*T0_x                 + ps0_y*T0_y)                                    &
                                  -rho * Vpar0 * Btheta2                                                      )        * xjac * tstep * theta

  amat(2,6)   = amat(2,6)                                                                                                                      &
                ! --- Diamagnetic terms
                + tau_IC * v * R**4 * (P0_x_T * w0_y - P0_y_T * w0_x)                                                  * xjac * theta * tstep  &
                + tau_IC     * R**3 * P0_y_T  * (v_x * u0_x + v_y * u0_y)                                              * xjac * theta * tstep  &
                + tau_IC * v * R**4 * u0_xy * (P0_xx_T-P0_yy_T)                                                        * xjac * theta * tstep  &
                - tau_IC * v * R**4 * P0_xy_T * (u0_xx-u0_yy)                                                          * xjac * theta * tstep  &
                ! --- Diamagnetic viscosity
                - d2visco_dT2*T * R * W_dia   * (v_x*T0_x + v_y*T0_y)                                                  * xjac * theta * tstep  &
                - dvisco_dT*T   * R * W_dia   * (v_xx + v_x/R + v_yy)                                                  * xjac * theta * tstep  &
                - dvisco_dT     * R * W_dia_T * (v_x*T0_x + v_y*T0_y)                                                  * xjac * theta * tstep  &
                - visco_T       * R * W_dia_T * (v_xx + v_x/R + v_yy)                                                  * xjac * theta * tstep  &
                - dvisco_dT     * R * W_dia   * (v_x*T_x  + v_y*T_y )                                                  * xjac * theta * tstep  &
                ! --- Neoclassical term
                - amu_neo_prof * BB2 / (Btheta2+epsil)**2 * (ps0_x*v_x + ps0_y*v_y) * R                                                        &
                               * (  tau_IC                     * (ps0_x*r0_x*T   + ps0_y*r0_y*T  )                                             &
                                  + tau_IC                     * (ps0_x*r0  *T_x + ps0_y*r0  *T_y)                                             &
                                  + aki_neo_prof * tau_IC * r0 * (ps0_x*T_x+ps0_y*T_y)             )                   * xjac * tstep * theta
  
  amat(2,7)   = amat(2,7)                                                                                                                      &
                ! --- Neoclassical term
                + amu_neo_prof * BB2 / (Btheta2+epsil) * r0 * vpar * (ps0_x*v_x + ps0_y*v_y) * R                       * xjac * tstep * theta 
  
  
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
  rhs(3) = - ( v_x * ps0_x  + v_y * ps0_y + v*zj0) / R * xjac

  
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
  amat(3,1) = (v_x * psi_x + v_y * psi_y ) / R          * xjac 

  amat(3,3) = v * zj / R                                * xjac  

  
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
  rhs(4) = rhs(4)                                                                                  &
           - ( v_x * u0_x   + v_y * u0_y  + v*w0)                              * R * xjac
  
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
  amat(4,2) = (v_x * u_x + v_y * u_y) * R               * xjac 

  amat(4,4) =  v * w * R                                                                * xjac 
                    
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
  Bgrad_rho_k_star = ( F0 / R * v_p                ) / R
              
  ! ----------------------------              
  ! --- The RHS term (main part)              
  rhs(5)   =                                                                                    &
             ! --- Time derivative
             + zeta * v * R                                                     * xjac *delta_g(5)& 
             ! --- Div(rho.V)
             + v * R**2 * ( r0_x * u0_y - r0_y * u0_x)                          * xjac * tstep  &
             + v * 2.d0 * R * r0 * u0_y                                         * xjac * tstep  &
             - v * F0 / R * Vpar0 * r0_p                                        * xjac * tstep  &
             - v * F0 / R * r0    * vpar0_p                                     * xjac * tstep  &
             - v * Vpar0 * (r0_x    * ps0_y - r0_y    * ps0_x)                  * xjac * tstep  &
             - v * r0    * (vpar0_x * ps0_y - vpar0_y * ps0_x)                  * xjac * tstep  &
             ! --- Source
             + v * R * total_rho_source                                         * xjac * tstep  &
             ! --- Diffusivity
             - (D_par-D_prof) * R / BB2 * Bgrad_rho_star * Bgrad_rho            * xjac * tstep  &
             - D_prof * R  * (v_x*r0_x + v_y*r0_y)                              * xjac * tstep  &
             ! --- Numerical diffusivity
             - D_perp_numm * (v_xx  + v_x /R + v_yy )                                           &
                           * (r0_xx + r0_x/R + r0_yy) * R                       * xjac * tstep 

  rhs_k(5) = - (D_par-D_prof) * R / BB2 * Bgrad_rho_k_star * Bgrad_rho          * xjac * tstep  &
             - D_prof * R  * ( v_p*r0_p * eps_cyl**2 /R**2 )                    * xjac * tstep 
  
  ! -----------------------------------    
  ! --- The RHS term (diamagnetic part)       
  rhs(5)   = rhs(5)                                                                             &
             + tau_IC * v * 2.d0 * p0_y * R                                     * xjac * tstep  
  
  
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
  amat(5,1)    = - (D_par-D_prof) * R * BB2_psi/ BB2**2 * Bgrad_rho_star        * Bgrad_rho     * xjac * theta * tstep  &
                 + (D_par-D_prof) * R / BB2             * Bgrad_rho_star_psi    * Bgrad_rho     * xjac * theta * tstep  &
                 + (D_par-D_prof) * R / BB2             * Bgrad_rho_star        * Bgrad_rho_psi * xjac * theta * tstep  &
                 + v * Vpar0 * (r0_x * psi_y - r0_y * psi_x)                                    * xjac * theta * tstep  &
                 + v * r0 * (vpar0_x * psi_y - vpar0_y * psi_x)                                 * xjac * theta * tstep 

  amat_k(5,1)  = - (D_par-D_prof) * R * BB2_psi/ BB2**2 * Bgrad_rho_k_star * Bgrad_rho          * xjac * theta * tstep  &
                 + (D_par-D_prof) * R / BB2             * Bgrad_rho_k_star * Bgrad_rho_psi      * xjac * theta * tstep 

  amat(5,2)    = - v * R**2 * ( r0_x * u_y - r0_y * u_x)                                        * xjac * theta * tstep  &
                 - v * 2.d0 * R * r0 * u_y                                                      * xjac * theta * tstep 

  amat(5,5)    = + v * rho * R                                                                  * xjac * (1.d0 + zeta)  &
                 - v * R**2 * ( rho_x * u0_y - rho_y * u0_x)                                    * xjac * theta * tstep  &
                 - v * 2.d0 * R * rho * u0_y                                                    * xjac * theta * tstep  &
                 + (D_par-D_prof) * R / BB2 * Bgrad_rho_star * Bgrad_rho_rho                    * xjac * theta * tstep  &
                 + D_prof * R  * (v_x*rho_x + v_y*rho_y )                                       * xjac * theta * tstep  &
                 + v * Vpar0 * (rho_x * ps0_y - rho_y * ps0_x)                                  * xjac * theta * tstep  &
                 + v * rho * (vpar0_x * ps0_y - vpar0_y * ps0_x)                                * xjac * theta * tstep  &
                 + v * rho * F0 / R * vpar0_p                                                   * xjac * theta * tstep  &
                 + D_perp_numm * (v_xx   + v_x  /R + v_yy  )                                                            &
                               * (rho_xx + rho_x/R + rho_yy) * R                                * xjac * theta * tstep  

  amat_k(5,5)  = + (D_par-D_prof) * R / BB2 * Bgrad_rho_k_star * Bgrad_rho_rho                  * xjac * theta * tstep 

  amat_n(5,5)  = + (D_par-D_prof) * R / BB2 * Bgrad_rho_star   * Bgrad_rho_rho_n                * xjac * theta * tstep  &
                 + v * F0 / R * Vpar0 * rho_p                                                   * xjac * theta * tstep 

  amat_kn(5,5) = + (D_par-D_prof) * R / BB2 * Bgrad_rho_k_star * Bgrad_rho_rho_n                * xjac * theta * tstep  &
                 + D_prof * R  * ( v_p*rho_p * eps_cyl**2 /R**2 )                               * xjac * theta * tstep 

  amat(5,7)    = + v * F0 / R * Vpar * r0_p                                                     * xjac * theta * tstep  &
                 + v * Vpar * (r0_x * ps0_y - r0_y * ps0_x)                                     * xjac * theta * tstep  &
                 + v * r0 * (vpar_x * ps0_y - vpar_y * ps0_x)                                   * xjac * theta * tstep 

  amat_n(5,7)  = + v * r0 * F0 / R * vpar_p                                                     * xjac * theta * tstep
                    
  ! ------------------------------------
  ! --- The LHS terms (diamagnetic part)
  amat(5,5)    = amat(5,5)                                                                                              &
                 - tau_IC * v * 2.d0 * (rho_y*T0 + rho*T0_y) * R                                * xjac * theta * tstep  
  
  amat(5,6)    = amat(5,6)                                                                                              &
                 - tau_IC * v * 2.d0 * (T_y*r0 + T*r0_y) * R                                    * xjac * theta * tstep 

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
                     + F0 / R * T0_p                 ) / R
  Bgrad_T_star     = ( v_x   * ps0_y - v_y   * ps0_x ) / R
  Bgrad_T_k_star   = ( F0 / R * v_p                  ) / R
              
  ! -----------------------------             
  ! --- The RHS terms (main part)
  rhs(6) =                                                                                              &
             ! --- Time derivative
             + zeta * v * r0_corr * R                                                   * xjac *delta_g(6)&
             + zeta * v * T0_corr * R                                                   * xjac *delta_g(5)&
             ! --- Convective terms
             + v * r0 * R**2 * ( T0_x * u0_y - T0_y * u0_x)                             * xjac * tstep  &
             + v * T0 * R**2 * ( r0_x * u0_y - r0_y * u0_x)                             * xjac * tstep  &
             + v * r0 * GAMMA * T0 * u0_y * 2.d0 * R                                    * xjac * tstep  &
             - v * r0 *         F0 / R * Vpar0 * T0_p                                   * xjac * tstep  &
             - v * T0 *         F0 / R * Vpar0 * r0_p                                   * xjac * tstep  &
             - v * r0 * GAMMA * F0 / R * T0    * vpar0_p                                * xjac * tstep  &
             - v * r0 *      Vpar0 * (T0_x    * ps0_y - T0_y    * ps0_x)                * xjac * tstep  &
             - v * T0 *      Vpar0 * (r0_x    * ps0_y - r0_y    * ps0_x)                * xjac * tstep  &
             - v * r0 * GAMMA * T0 * (vpar0_x * ps0_y - vpar0_y * ps0_x)                * xjac * tstep  &
             ! --- Source
             + v * R * heat_source                                                      * xjac * tstep  &
             ! --- Conductivity
             - (K_par-K_prof) * R / BB2 * Bgrad_T_star * Bgrad_T                        * xjac * tstep  &
             - K_prof * R * (v_x*T0_x + v_y*T0_y )                                      * xjac * tstep  &
             ! --- Numerical conductivity
             - K_perp_numm * (v_xx  + v_x /R + v_yy )                                                   &
                           * (T0_xx + T0_x/R + T0_yy) * R                               * xjac * tstep  &
             ! --- Ohmic heating
             + v * (gamma-1.d0) * eta_T_ohm * (zj0 / R)**2.d0 * R                       * xjac  * tstep
             

  rhs_k(6) = - (K_par-K_prof) * R / BB2 * Bgrad_T_k_star * Bgrad_T                      * xjac * tstep  &
             - K_prof * R * ( v_p*T0_p /R**2 )                                          * xjac * tstep 

  
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
  Bgrad_T_star_psi  = ( v_x  * psi_y - v_y  * psi_x ) / R
  Bgrad_T_psi       = ( T0_x * psi_y - T0_y * psi_x ) / R
  Bgrad_T_T         = ( T_x  * ps0_y - T_y  * ps0_x ) / R
  Bgrad_T_T_n       = ( F0 / R * T_p                ) / R

  ! -----------------------------             
  ! --- The LHS terms (main part)
  amat(6,1)    = + v         * r0 * Vpar0 * (T0_x    * psi_y - T0_y    * psi_x)                         * xjac * theta * tstep  &
                 + v         * T0 * Vpar0 * (r0_x    * psi_y - r0_y    * psi_x)                         * xjac * theta * tstep  &
                 + v * GAMMA * r0 * T0    * (vpar0_x * psi_y - vpar0_y * psi_x)                         * xjac * theta * tstep  &
                 - (K_par-K_prof) * R * BB2_psi / BB2**2    * Bgrad_T_star     * Bgrad_T                * xjac * theta * tstep  &
                 + (K_par-K_prof) * R / BB2                 * Bgrad_T_star_psi * Bgrad_T                * xjac * theta * tstep  &
                 + (K_par-K_prof) * R / BB2                 * Bgrad_T_star     * Bgrad_T_psi            * xjac * theta * tstep 

  amat_k(6,1)  = - (K_par-K_prof) * R * BB2_psi / BB2**2    * Bgrad_T_k_star   * Bgrad_T                * xjac * theta * tstep  &
                 + (K_par-K_prof) * R / BB2                 * Bgrad_T_k_star   * Bgrad_T_psi            * xjac * theta * tstep 

  amat(6,2)    = - v * r0 * R**2 * ( T0_x * u_y - T0_y * u_x)                                           * xjac * theta * tstep  &
                 - v * T0 * R**2 * ( r0_x * u_y - r0_y * u_x)                                           * xjac * theta * tstep  &
                 - v * 2.d0 * GAMMA * r0 * R * T0 * u_y                                                 * xjac * theta * tstep 

  amat(6,3)   = - v * (gamma-1.d0) * eta_T_ohm * 2.d0 * zj * zj0/(R**2.d0) * R                          * xjac * theta * tstep

  amat(6,5)    = + v * rho * T0_corr * R                                                                * xjac * (1.d0 + zeta)  &
                 - v * rho * R**2    * (T0_x  * u0_y  - T0_y  * u0_x )                                  * xjac * theta * tstep  &
                 - v * T0  * R**2    * (rho_x * u0_y  - rho_y * u0_x )                                  * xjac * theta * tstep  &
                 + v * rho * Vpar0   * (T0_x  * ps0_y - T0_y  * ps0_x)                                  * xjac * theta * tstep  &
                 + v * T0  * Vpar0   * (rho_x * ps0_y - rho_y * ps0_x)                                  * xjac * theta * tstep  &
                 + v * rho * Vpar0   * F0/R * T0_p                                                      * xjac * theta * tstep  &
                 - v * GAMMA * rho * T0 * u0_y * 2.d0 * R                                               * xjac * theta * tstep  &
                 + v * GAMMA * rho * T0 * F0/R * Vpar0_p                                                * xjac * theta * tstep  &
                 + v * GAMMA * rho * T0 * (Vpar0_x * ps0_y - Vpar0_y * ps0_x)                           * xjac * theta * tstep

  amat_n(6,5)  = + v * T0 *                F0 / R * Vpar0 * rho_p                                       * xjac * theta * tstep
  
  amat(6,6)    = + v * r0_corr * T * R                                                                  * xjac * (1.d0 + zeta)  &
                 - v * r0 * R**2    * (T_x  * u0_y - T_y  * u0_x)                                       * xjac * theta * tstep  &
                 - v * T  * R**2    * (r0_x * u0_y - r0_y * u0_x)                                       * xjac * theta * tstep  &
                 + v * T  *                F0 / R * Vpar0 * r0_p                                        * xjac * theta * tstep  &
                 + v * r0 * Vpar0   * (T_x  * ps0_y - T_y  * ps0_x)                                     * xjac * theta * tstep  &
                 + v * T  * Vpar0   * (r0_x * ps0_y - r0_y * ps0_x)                                     * xjac * theta * tstep  &
                 - v * r0 * GAMMA * T * R * u0_y * 2.d0                                                 * xjac * theta * tstep  &
                 + v * r0 * GAMMA * T * F0/R * Vpar0_p                                                  * xjac * theta * tstep  &
                 + v * r0 * GAMMA * T * (Vpar0_x * ps0_y - Vpar0_y * ps0_x)                             * xjac * theta * tstep  &
                 + (K_par-K_prof) * R / BB2 * Bgrad_T_star * Bgrad_T_T                                  * xjac * theta * tstep  &
                 + dK_par * T     * R / BB2 * Bgrad_T_star * Bgrad_T                                    * xjac * theta * tstep  &
                 + K_prof * R * (v_x*T_x + v_y*T_y )                                                    * xjac * theta * tstep  & 
                 + K_perp_numm * (v_xx + v_x/R + v_yy)                                                                          &
                               * (T_xx + T_x/R + T_yy) * R                                              * xjac * theta * tstep  &
                 - v * T * (gamma-1.d0) * deta_dT_ohm * (zj0 / R)**2.d0 * R                             * xjac * theta * tstep
  
  amat_k(6,6)  = + (K_par-K_prof) * R / BB2 * Bgrad_T_k_star * Bgrad_T_T                                * xjac * theta * tstep  &
                 + dK_par * T     * R / BB2 * Bgrad_T_k_star * Bgrad_T                                  * xjac * theta * tstep 
              
  amat_n(6,6)  = + (K_par-K_prof) * R / BB2 * Bgrad_T_star   * Bgrad_T_T_n                              * xjac * theta * tstep  &
                 + v * r0 * Vpar0  * F0/R * T_p                                                         * xjac * theta * tstep 

  amat_kn(6,6) = + (K_par-K_prof) * R / BB2 * Bgrad_T_k_star * Bgrad_T_T_n                              * xjac * theta * tstep  &
                 +  K_prof          * R * (v_p*T_p /R**2 )                                              * xjac * theta * tstep 

  amat(6,7)    = + v * r0 * F0/R * Vpar * T0_p                                                          * xjac * theta * tstep  &
                 + v * T0 *                F0 / R * Vpar * r0_p                                         * xjac * theta * tstep  &
                 + v * r0 * Vpar * (T0_x * ps0_y - T0_y * ps0_x)                                        * xjac * theta * tstep  &
                 + v * T0 * Vpar * (r0_x * ps0_y - r0_y * ps0_x)                                        * xjac * theta * tstep  &
                 + v * GAMMA * r0 * T0 * (Vpar_x * ps0_y - Vpar_y * ps0_x)                              * xjac * theta * tstep 
     
  amat_n(6,7)  = + v * GAMMA * r0 * T0 * F0/R * Vpar_p                                                  * xjac * theta * tstep       

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
  rhs(7) =                                                                                                      &
             ! --- Time derivative terms (including dB/dt)
             + zeta * v * r0_corr * F0**2 / R                                                   * xjac *delta_g(7)&  
             + zeta * v * r0_corr * vpar0 * (ps0_x * delta_ps_x + ps0_y * delta_ps_y) / R       * xjac          &
             ! --- Convection terms
             - 0.5d0 * r0 * vpar0**2 * BB2 * (ps0_x * v_y  - ps0_y * v_x)                       * xjac * tstep  &
             - 0.5d0 * v  * vpar0**2 * BB2 * (ps0_x * r0_y - ps0_y * r0_x)                      * xjac * tstep  &
             + 0.5d0 * v  * vpar0**2 * BB2 * F0 / R * r0_p                                      * xjac * tstep  &
             ! --- Parallel pressure gradient
             - v * F0 / R * P0_p                                                                * xjac * tstep  &
             - v * (P0_x * ps0_y - P0_y * ps0_x)                                                * xjac * tstep  &
             ! --- Density source term
             - v * total_rho_source * vpar0 * BB2 * R                                           * xjac * tstep  &
             ! --- Numerical viscosity
             - visco_par_numm * (v_xx     + v_x    /R + v_yy    )                                               &
                              * (vpar0_xx + vpar0_x/R + vpar0_yy) * R                           * xjac * tstep 
  
  rhs_k(7) = + 0.5d0 * r0 * vpar0**2 * BB2 * F0 / R * v_p                                       * xjac * tstep 

  ! --- The RHS term (Viscosity and source of parallel rotation)
  if (normalized_velocity_profile) then
    rhs(7)   = rhs(7)                                                                                              &
      + visco_parr * (v_x * Vt0_x   + v_y * Vt0_y)   * R                                           * xjac * tstep  &
      - visco_parr * (v_x * vpar0_x + v_y * vpar0_y) * R                                           * xjac * tstep  
  else
    rhs(7)   = rhs(7)                                                                                                                &
      - visco_parr * v_x * ( vpar0_x * F0**2 / R**2 -2.d0 * vpar0 * F0**2 / R**3 - 2.d0 * PI * F0 * Omega_tor0_x ) * R * xjac * tstep &
      - visco_parr * v_y * ( vpar0_y * F0**2 / R**2 - 2.d0 * PI * F0 * Omega_tor0_y ) * R                              * xjac * tstep            
  endif

  ! -----------------------------------    
  ! --- The RHS term (Neoclassical part)              
  rhs(7)   = rhs(7)                                                                                             &
             ! --- Neoclassic term
             + v * amu_neo_prof * BB2 / (Btheta2+epsil) * R                                                     &
                 * (  r0                         * (ps0_x*u0_x + ps0_y*u0_y)                                    &
                    + tau_IC                     * (ps0_x*P0_x + ps0_y*P0_y)                                    &
                    + aki_neo_prof * tau_IC * r0 * (ps0_x*T0_x + ps0_y*T0_y)                                    &
                    - r0 * Vpar0 * Btheta2                                      )               * xjac * tstep 
  
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
  real*8 :: Vt_x_psi, Vt_y_psi, Omega_tor_x_psi, Omega_tor_y_psi
    
  Btheta2_psi = 2.d0 * (psi_x * ps0_x + psi_y * ps0_y ) / R**2
  if (normalized_velocity_profile) then
    Vt_x_psi    = dV_dpsi_source * psi_x
    Vt_y_psi    = dV_dpsi_source * psi_y
  else
    Omega_tor_x_psi    = dV_dpsi_source * psi_x
    Omega_tor_y_psi    = dV_dpsi_source * psi_y
  endif
  
  ! -----------------------------
  ! --- The LHS terms (main part)
  amat(7,1)   = + v * (P0_x * psi_y - P0_y * psi_x)                                                             * xjac * theta * tstep  &
                + 0.5d0 * r0 * vpar0**2 * BB2     * (psi_x * v_y  - psi_y * v_x)                                * xjac * theta * tstep  &
                + 0.5d0 * r0 * vpar0**2 * BB2_psi * (ps0_x * v_y  - ps0_y * v_x)                                * xjac * theta * tstep  &
                + 0.5d0 * v  * vpar0**2 * BB2     * (psi_x * r0_y - psi_y * r0_x)                               * xjac * theta * tstep  &
                + 0.5d0 * v  * vpar0**2 * BB2_psi * (ps0_x * r0_y - ps0_y * r0_x)                               * xjac * theta * tstep  &
                - 0.5d0 * v  * vpar0**2 * BB2_psi * F0 / R * r0_p                                               * xjac * theta * tstep  &
                + v * total_rho_source * vpar0 * BB2_psi * R                                                    * xjac * theta * tstep  &
                + v * r0 * vpar0 / R * (ps0_x * psi_x + ps0_y * psi_y)                                          * xjac * (1.d0 + zeta) 
  if (normalized_velocity_profile) then
    amat(7,1) = amat(7,1)                                                                                                               &
                - visco_parr * (v_x * Vt_x_psi   + v_y * Vt_y_psi)   * R                                        * xjac * theta * tstep
  else
    amat(7,1) = amat(7,1)                                                                                                               &
                - visco_parr * 2.d0 * PI * F0 * (v_x * Omega_tor_x_psi   + v_y * Omega_tor_y_psi)   * R         * xjac * theta * tstep
  endif
  
  amat_k(7,1) = - 0.5d0 * r0 * vpar0**2 * BB2_psi * F0 / R * v_p                                                * xjac * theta * tstep 
           
  amat(7,5)   = + v * (rho_x * T0   * ps0_y - rho_y * T0   * ps0_x)                                             * xjac * theta * tstep  &
                + v * (rho   * T0_x * ps0_y - rho   * T0_y * ps0_x)                                             * xjac * theta * tstep  &
                + v * F0 / R * rho * T0_p                                                                       * xjac * theta * tstep  & 
                + 0.5d0 * rho * vpar0**2 * BB2 * (ps0_x * v_y   - ps0_y * v_x)                                  * xjac * theta * tstep  &
                + 0.5d0 * v   * vpar0**2 * BB2 * (ps0_x * rho_y - ps0_y * rho_x)                                * xjac * theta * tstep

  amat_k(7,5) = - 0.5d0 * rho * vpar0**2 * BB2 * F0 / R * v_p                                                   * xjac * theta * tstep 

  amat_n(7,5) = + v * F0 / R * rho_p * T0                                                                       * xjac * theta * tstep  & 
                - 0.5d0 * v   * vpar0**2 * BB2 * F0 / R * rho_p                                                 * xjac * theta * tstep 

  amat(7,6)   = + v * (T_x * r0   * ps0_y - T_y * r0   * ps0_x)                                                 * xjac * theta * tstep  &
                + v * (T   * r0_x * ps0_y - T   * r0_y * ps0_x)                                                 * xjac * theta * tstep  &
                + v * F0 / R * T * r0_p                                                                         * xjac * theta * tstep
  
  amat_n(7,6) = + v * F0 / R * T_p * r0                                                                         * xjac * theta * tstep 

  amat(7,7)   = + v * Vpar * r0_corr * F0**2 / R                                                                * xjac * (1.d0 + zeta)  &
                + r0 * vpar0 * vpar * BB2 * (ps0_x * v_y  - ps0_y * v_x)                                        * xjac * theta * tstep  &
                + v  * vpar0 * vpar * BB2 * (ps0_x * r0_y - ps0_y * r0_x)                                       * xjac * theta * tstep  &
                - v  * vpar0 * vpar * BB2 * F0 / R * r0_p                                                       * xjac * theta * tstep  &
                + v * total_rho_source  * vpar * BB2 * R                                                        * xjac * theta * tstep  &
                + visco_par_numm * (v_xx    + v_x   /R + v_yy   )                                                                       &
                                 * (vpar_xx + vpar_x/R + vpar_yy) * R                                           * xjac * theta * tstep 
  
  amat_k(7,7) = - r0 * vpar0 * vpar * BB2 * F0 / R * v_p                                                        * xjac * theta * tstep 

  ! -----------------------------------    
  ! --- The LHS term (Viscosity and source of parallel rotation)
  if (normalized_velocity_profile) then
    amat(7,7) = amat(7,7)                                                                          &
      + visco_parr * (v_x * Vpar_x + v_y * Vpar_y) * R * xjac * theta * tstep
  else
    amat(7,7) = amat(7,7)                                                                          &
      + visco_parr * F0**2 / R**2 * (v_x * (Vpar_x -2.d0 * vpar / R) + v_y * Vpar_y) * R * xjac * theta * tstep
  endif
  
  ! -----------------------------------
  ! --- The LHS terms (Neoclassic part)
  amat(7,1)   = amat(7,1)                                                                                                               &
                - v * amu_neo_prof * BB2 / (Btheta2+epsil) * R                                                                          &
                                   * (  r0                         * (psi_x*u0_x + psi_y*u0_y)                                          &
                                      + tau_IC                     * (psi_x*P0_x + psi_y*P0_y)                                          &
                                      + aki_neo_prof * tau_IC * r0 * (psi_x*T0_x + psi_y*T0_y)  )               * xjac * theta * tstep  &
                + v * amu_neo_prof * Btheta2_psi * BB2 / Btheta2**2 * R                                                         &
                                   * (  r0                         * (ps0_x*u0_x + ps0_y*u0_y)                                          &
                                      + tau_IC                     * (ps0_x*P0_x + ps0_y*P0_y)                                          &
                                      + aki_neo_prof * tau_IC * r0 * (ps0_x*T0_x + ps0_y*T0_y)  )               * xjac * theta * tstep
           
  amat(7,2)   = amat(7,2)                                                                                                               &
                - v * amu_neo_prof * BB2 / (Btheta2+epsil) * r0 * (ps0_x*u_x + ps0_y*u_y) * R                   * xjac * theta * tstep 
  
  amat(7,5)   = amat(7,5)                                                                                                               &
                - v * amu_neo_prof * BB2 / (Btheta2+epsil) * R                                                                          &
                                   * (  rho * (ps0_x*u0_x + ps0_y*u0_y)                                                                 &
                                      + tau_IC                      * (ps0_x*rho_x*T0   + ps0_y*rho_y*T0  )                             &
                                      + tau_IC                      * (ps0_x*rho  *T0_x + ps0_y*rho  *T0_y)                             &
                                      + aki_neo_prof * tau_IC * rho * (ps0_x*T0_x       + ps0_y*T0_y)                                   &
                                      - rho * Vpar0 * Btheta2                                                 ) * xjac * tstep * theta

  amat(7,6)   = amat(7,6)                                                                                                               &
                - v * amu_neo_prof * BB2 / (Btheta2+epsil) * R                                                                          &
                                   * (  tau_IC                     * (ps0_x*r0_x*T   + ps0_y*r0_y*T  )                                  &
                                      + tau_IC                     * (ps0_x*r0  *T_x + ps0_y*r0  *T_y)                                  &
                                      + aki_neo_prof * tau_IC * r0 * (ps0_x*T_x      + ps0_y*T_y     )   )      * xjac * tstep * theta
  
  amat(7,7)   = amat(7,7)                                                                                                               &
                + v * amu_neo_prof * BB2 * r0 * vpar * R                                                        * xjac * tstep * theta 


  return

end subroutine ELM_main_lhs_7










