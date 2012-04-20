subroutine bootstrap_current_rhs(BigR, minRad, R_axis, &
                                 psi_axis, psi_bnd,    &
                                 ps0, ps0_x, ps0_y,    &
				 r0,  r0_x,  r0_y,     &
				 Ti0, Ti0_x, Ti0_y,    &
                                 Te0, Te0_x, Te0_y,    &
				 Jb)
!----------------------------------------------------------------------------------------------------------
! calculates the bootstrap current for the RHS of the matrix, based on the Wesson formula (Appendix 14.12)
! Note that we do not linearise all the coeffs C1-4, they are considered as local constants...
!----------------------------------------------------------------------------------------------------------

  use constants
  use phys_module

  implicit none
  ! --- Routine parameters
  real*8, intent(in)  :: BigR, minRad, R_axis
  real*8, intent(in)  :: psi_axis, psi_bnd
  real*8, intent(in)  :: ps0, ps0_x, ps0_y
  real*8, intent(in)  :: r0,  r0_x,  r0_y
  real*8              :: rho, rho_x, rho_y, drho
  real*8, intent(in)  :: Ti0, Ti0_x, Ti0_y
  real*8              :: Ti,  Ti_x,  Ti_y,  dTi
  real*8, intent(in)  :: Te0, Te0_x, Te0_y
  real*8              :: Te,  Te_x,  Te_y,  dTe
  real*8              :: grad_psi, psi_n, X
  real*8              :: Nue, Nui
  real*8              :: DD, C1, C2, C3, C4, AC4, BC4
  real*8, intent(out) :: Jb

  ! --- Need central_density
  if (central_density .lt. 1.d-6) then
    write(*,*)'**********************!!WARNING!!*************************'
    write(*,*)'Element_matrix asks for the bootstrap current,'
    write(*,*)'but the central_density is not defined in the input file.'
    write(*,*)'Returning zero bootstrap current...'
    write(*,*)'**********************!!WARNING!!*************************'
    Jb = 0.d0
    return
  endif
  
  ! --- Renormalise temperature and density 
  ! --- Note for us density is n, not mi*n,
  ! --- but that's ok since formula is with p = n*T
  ! --- ie. rho*T for us...
  ! --- Temperature in Joules
  ! --- Density in 1/(cubic meters)
  rho   = r0   * central_density
  rho_x = r0_x * central_density
  rho_y = r0_y * central_density
  Ti    = Ti0   / (MU_ZERO*central_density)
  Ti_x  = Ti0_x / (MU_ZERO*central_density)
  Ti_y  = Ti0_y / (MU_ZERO*central_density)
  Te    = Te0   / (MU_ZERO*central_density)
  Te_x  = Te0_x / (MU_ZERO*central_density)
  Te_y  = Te0_y / (MU_ZERO*central_density)
        
  ! --- Derivatives with respect to psi
  drho = rho_x/ps0_x + rho_y/ps0_y
  dTi  = Ti_x /ps0_x + Ti_y /ps0_y
  dTe  = Te_x /ps0_x + Te_y /ps0_y
        
  ! --- Psi variables, including r~a*sqrt(psi) and X=sqrt(2*r/R0)
  psi_n    = (ps0 - psi_axis) / (psi_bnd - psi_axis)
  grad_psi = (ps0_x*ps0_x + ps0_y*ps0_y)**0.5d0
  X        = sqrt(2.d0*minRad*sqrt(psi_n)/R_axis)

  ! --- Nue* formula from Wesson : Nue* = R*q / ( eps**(3/2) * (Te/me)**(1/2) * Taue )
  ! --- where Taue is the electron collision time (formula for Nui* is very similar)
  Nue = 4.0d-28 * F0 * rho / ( Te**2.d0 * X**3.d0 * grad_psi )
  Nui = 2.7d-28 * F0 * rho / ( Ti**2.d0 * X**3.d0 * grad_psi )

  ! --- Coefficients
  DD  = 2.4d0 + 5.4d0*X  + 2.6d0*X*X

  C1  = (4.d0 + 2.6d0*X) / ( (1.d0 + 1.02d0*sqrt(Nue) + 1.07d0*Nue) * (1.d0 + 0.38d0*Nue*X**3.d0) )

  C2  = C1*Ti/Te
  
  C3  = (7.d0 + 6.5d0*X) / ( (1.d0 + 0.57d0*sqrt(Nue) + 0.61d0*Nue) * (1.d0 + 0.22d0*Nue*X**3.d0) ) - 2.5d0*C1

  AC4 = (  -1.17d0/(1.d0+0.46d0*X) + 0.35d0*sqrt(Nui)) / (1.d0 + 0.7d0*sqrt(Nui)) + 0.26d0*Nui*Nui*X**6.d0
  BC4 = (1.d0 - 0.125d0*Nui*Nui*X**6.d0) * (1.d0 + 0.125d0*Nue*Nue*X**6.d0)
  C4  = AC4*C2 / BC4

  ! --- Bootstrap Current
  Jb = BigR*X*rho*Te/DD * ( C1*(dTe/Te+drho/rho) + C2*(dTi/Ti+drho/rho) + C3*dTe/Te + C4*dTi/Ti)

  ! --- Current with denormalisation
  Jb = Jb / MU_ZERO


return
end subroutine bootstrap_current_rhs








subroutine bootstrap_current_lhs(BigR, minRad, R_axis, &
                                 psi_axis, psi_bnd,    &
                                 ps0, ps0_x, ps0_y,    &
                                 psi, psi_x, psi_y,    &
				 r0,  r0_x,  r0_y,     &
				 rrho,rrho_x,rrho_y,   &
				 Ti0, Ti0_x, Ti0_y,    &
				 TTi, TTi_x, TTi_y,    &
                                 Te0, Te0_x, Te0_y,    &
                                 TTe, TTe_x, TTe_y,    &
				 dJb_psi, dJb_rho, dJb_Ti, dJb_Te)
!----------------------------------------------------------------------------------------------------------
! calculates the bootstrap current for the LHS of the matrix, based on the Wesson formula (Appendix 14.12)
! Note that we do not linearise all the coeffs C1-4, they are considered as local constants...
! Note Te0 is the variable, Te is the renormalised variable, and TTe is the test function
!----------------------------------------------------------------------------------------------------------

  use constants
  use phys_module

  implicit none
  ! --- Routine parameters
  real*8, intent(in)  :: BigR, minRad, R_axis
  real*8, intent(in)  :: psi_axis, psi_bnd
  real*8, intent(in)  :: ps0, ps0_x, ps0_y
  real*8, intent(in)  :: psi, psi_x, psi_y
  real*8, intent(in)  :: r0,  r0_x,  r0_y
  real*8              :: rho, rho_x, rho_y, drho
  real*8, intent(in)  :: rrho,rrho_x,rrho_y
  real*8              :: vrho,vrho_x,vrho_y
  real*8              :: drho_rho, drho_psi
  real*8, intent(in)  :: Ti0, Ti0_x, Ti0_y
  real*8              :: Ti,  Ti_x,  Ti_y,  dTi
  real*8, intent(in)  :: TTi, TTi_x, TTi_y
  real*8              :: vTi, vTi_x, vTi_y
  real*8              :: dTi_ti, dTi_psi
  real*8, intent(in)  :: Te0, Te0_x, Te0_y
  real*8              :: Te,  Te_x,  Te_y,  dTe
  real*8, intent(in)  :: TTe, TTe_x, TTe_y
  real*8              :: vTe, vTe_x, vTe_y
  real*8              :: dTe_Te, dTe_psi
  real*8              :: grad_psi,  psi_n,   X
  real*8              :: dgrad_psi, dpsi_n, dX
  real*8              :: Nue, dNue_Te, dNue_psi, dNue_rho
  real*8              :: Nui, dNui_Ti, dNui_psi, dNui_rho
  real*8              :: DD, dDD
  real*8              :: C1,  dC1_psi,  dC1_rho,  dC1_Te
  real*8              :: C2,  dC2_psi,  dC2_rho,  dC2_Te,  dC2_Ti
  real*8              :: C3,  dC3_psi,  dC3_rho,  dC3_Te
  real*8              :: AC4, dAC4_psi, dAC4_rho,          dAC4_Ti
  real*8              :: BC4, dBC4_psi, dBC4_rho, dBC4_Te, dBC4_Ti
  real*8              :: C4,  dC4_psi,  dC4_rho,  dC4_Te,  dC4_Ti
  real*8              :: Jb1, dJb1_psi, dJb1_rho, dJb1_Te
  real*8              :: Jb2, dJb2_psi, dJb2_rho, dJb2_Te, dJb2_Ti
  real*8              :: Jb3, dJb3_psi, dJb3_rho, dJb3_Te
  real*8              :: Jb4, dJb4_psi, dJb4_rho, dJb4_Te, dJb4_Ti
  real*8              :: Jb
  real*8, intent(out) ::      dJb_psi, dJb_rho, dJb_Ti,  dJb_Te

  ! --- Need central_density
  if (central_density .lt. 1.d-6) then
    dJb_psi = 0.d0
    dJb_rho = 0.d0
    dJb_Ti  = 0.d0
    dJb_Te  = 0.d0
    return
  endif
  
  ! --- Renormalise temperature and density 
  ! --- Note for us density is n, not mi*n,
  ! --- but that's ok since formula is with p = n*T
  ! --- ie. rho*T for us...
  ! --- Temperature in Joules
  ! --- Density in 1/(cubic meters)
  rho    = r0     * central_density
  rho_x  = r0_x   * central_density
  rho_y  = r0_y   * central_density
  vrho   = rrho   * central_density
  vrho_x = rrho_x * central_density
  vrho_y = rrho_y * central_density
  Ti     = Ti0   / (MU_ZERO*central_density)
  Ti_x   = Ti0_x / (MU_ZERO*central_density)
  Ti_y   = Ti0_y / (MU_ZERO*central_density)
  vTi    = TTi   / (MU_ZERO*central_density)
  vTi_x  = TTi_x / (MU_ZERO*central_density)
  vTi_y  = TTi_y / (MU_ZERO*central_density)
  Te     = Te0   / (MU_ZERO*central_density)
  Te_x   = Te0_x / (MU_ZERO*central_density)
  Te_y   = Te0_y / (MU_ZERO*central_density)
  vTe    = TTe   / (MU_ZERO*central_density)
  vTe_x  = TTe_x / (MU_ZERO*central_density)
  vTe_y  = TTe_y / (MU_ZERO*central_density)
        
  ! --- Derivatives with respect to psi
  drho     = rho_x /ps0_x + rho_y /ps0_y
  drho_rho = vrho_x/ps0_x + vrho_y/ps0_y
  drho_psi = rho_x /psi_x + rho_y /psi_y
  dTi      = Ti_x  /ps0_x + Ti_y  /ps0_y
  dTi_Ti   = vTi_x /ps0_x + vTi_y /ps0_y
  dTi_psi  = Ti_x  /psi_x + Ti_y  /psi_y
  dTe      = Te_x  /ps0_x + Te_y  /ps0_y
  dTe_Te   = vTe_x /ps0_x + vTe_y /ps0_y
  dTe_psi  = Te_x  /psi_x + Te_y  /psi_y
        
  ! --- Psi variables, including r~a*sqrt(psi) and X=sqrt(2*r/R0)
  psi_n     = (ps0 - psi_axis) / (psi_bnd - psi_axis)
  dpsi_n    =  psi             / (psi_bnd - psi_axis)
  grad_psi  = (ps0_x*ps0_x + ps0_y*ps0_y)**0.5d0
  dgrad_psi = (psi_x*ps0_x + psi_y*ps0_y)/(ps0_x*ps0_x + ps0_y*ps0_y)**0.5d0
  X         = sqrt(2.d0*minRad/R_axis) * psi_n**0.25d0
  dX        = sqrt(2.d0*minRad/R_axis) * 0.25d0 * dpsi_n / psi_n**0.75d0

  ! --- Nue* formula from Wesson : Nue* = R*q / ( eps**(3/2) * (Te/me)**(1/2) * Taue )
  ! --- where Taue is the electron collision time (formula for Nui* is very similar)
  Nue      = 4.0d-28 * F0 * rho / ( Te**2.d0 * X**3.d0 * grad_psi       )
  dNue_Te  = 4.0d-28 * F0 * rho / ( Te**3.d0 * X**3.d0 * grad_psi       ) * (-2.d0 * vTe)
  dNue_rho = 4.0d-28 * F0 * vrho/ ( Te**2.d0 * X**3.d0 * grad_psi       )
  dNue_psi = 4.0d-28 * F0 * rho / ( Te**2.d0 * X**4.d0 * grad_psi**2.d0 ) * (-3.d0*dX*grad_psi - X*dgrad_psi)
  Nui      = 2.7d-28 * F0 * rho / ( Ti**2.d0 * X**3.d0 * grad_psi       )
  dNui_Ti  = 2.7d-28 * F0 * rho / ( Ti**3.d0 * X**3.d0 * grad_psi       ) * (-2.d0 * vTi)
  dNui_rho = 2.7d-28 * F0 * vrho/ ( Ti**2.d0 * X**3.d0 * grad_psi       )
  dNui_psi = 2.7d-28 * F0 * rho / ( Ti**2.d0 * X**4.d0 * grad_psi**2.d0 ) * (-3.d0*dX*grad_psi - X*dgrad_psi)

  ! --- **************************************************** ---
  ! --- Coefficient DD
  DD       = 2.4d0 + 5.4d0*X  + 2.6d0*X*X
  dDD      =         5.4d0*dX + 5.2d0*X*dX


  ! --- **************************************************** ---
  ! --- Coefficient C1
  C1       = (4.d0 + 2.6d0*X) / ( (1.d0 + 1.02d0*sqrt(Nue)          + 1.07d0* Nue    ) * (1.d0 + 0.38d0* Nue    *X**3.d0)   )
  
  dC1_psi  =         2.6d0*dX / ( (1.d0 + 1.02d0*sqrt(Nue)          + 1.07d0* Nue    ) * (1.d0 + 0.38d0* Nue    *X**3.d0)   ) &
            -(4.d0 + 2.6d0*X) * ( (1.d0 + 1.02d0*sqrt(Nue)	    + 1.07d0* Nue    ) *	 1.14d0* Nue	*X**2.d0*dX   &
           			+ (1.d0 + 1.02d0*sqrt(Nue)	    + 1.07d0* Nue    ) *	 0.38d0*dNue_psi*X**3.d0      &
           			+ (1.d0 + 0.51d0/sqrt(Nue)*dNue_psi + 1.07d0*dNue_psi) * (1.d0 + 0.38d0* Nue	*X**3.d0)   ) &
           		      / ( (1.d0 + 1.02d0*sqrt(Nue)	    + 1.07d0* Nue    ) * (1.d0 + 0.38d0* Nue	*X**3.d0)   )**2.d0 
  
  dC1_rho  =-(4.d0 + 2.6d0*X) * ( (       0.51d0/sqrt(Nue)*dNue_rho + 1.07d0*dNue_rho) * (1.d0 + 0.38d0* Nue    *X**3.d0)     &
           			+ (1.d0 + 1.02d0*sqrt(Nue)	    + 1.07d0* Nue    ) * (	 0.38d0*dNue_rho*X**3.d0)   ) &
           		      / ( (1.d0 + 1.02d0*sqrt(Nue)	    + 1.07d0* Nue    ) * (1.d0 + 0.38d0* Nue	*X**3.d0)   )**2.d0 
  
  dC1_Te   =-(4.d0 + 2.6d0*X) * ( (       0.51d0/sqrt(Nue)*dNue_Te  + 1.07d0*dNue_Te ) * (1.d0 + 0.38d0* Nue    *X**3.d0)     &
           			+ (1.d0 + 1.02d0*sqrt(Nue)	    + 1.07d0* Nue    ) * (	 0.38d0*dNue_Te *X**3.d0)   ) &
           		      / ( (1.d0 + 1.02d0*sqrt(Nue)	    + 1.07d0* Nue    ) * (1.d0 + 0.38d0* Nue	*X**3.d0)   )**2.d0 


  ! --- **************************************************** ---
  ! --- Coefficient C2
  C2       = C1     * Ti/Te
  dC2_psi  = dC1_psi* Ti/Te
  dC2_rho  = dC1_rho* Ti/Te
  dC2_Ti   = C1     *vTi/Te
  dC2_Te   = dC1_Te * Ti/Te - C1*Ti*vTe/Te**2.d0
  

  ! --- **************************************************** ---
  ! --- Coefficient C3
  C3       = (7.d0 + 6.5d0*X) / ( (1.d0 + 0.57d0*sqrt(Nue)          + 0.61d0* Nue    ) * (1.d0 + 0.22d0* Nue    *X**3.d0)   ) - 2.5d0*C1
  
  dC3_psi  =         6.5d0*dX / ( (1.d0 + 0.57d0*sqrt(Nue)          + 0.61d0* Nue    ) * (1.d0 + 0.22d0* Nue    *X**3.d0)   ) &
            -(7.d0 + 6.5d0*X) * ( (1.d0 + 0.57d0*sqrt(Nue)	    + 0.61d0* Nue    ) *	 0.66d0* Nue	*X**2.d0*dX   &
           			+ (1.d0 + 0.57d0*sqrt(Nue)	    + 0.61d0* Nue    ) *	 0.22d0*dNue_psi*X**3.d0      &
           			+ (1.d0 + 0.28d0/sqrt(Nue)*dNue_psi + 0.61d0*dNue_psi) * (1.d0 + 0.22d0* Nue	*X**3.d0)   ) &
           		      / ( (1.d0 + 0.57d0*sqrt(Nue)	    + 0.61d0* Nue    ) * (1.d0 + 0.22d0* Nue	*X**3.d0)   )**2.d0 &
	    - 2.5d0*dC1_psi

  dC3_rho  =-(7.d0 + 6.5d0*X) * ( (       0.28d0/sqrt(Nue)*dNue_rho + 0.61d0*dNue_rho) * (1.d0 + 0.22d0* Nue    *X**3.d0)     & 
           			+ (1.d0 + 0.57d0*sqrt(Nue)	    + 0.61d0* Nue    ) * (	 0.22d0*dNue_rho*X**3.d0)   ) &
           		      / ( (1.d0 + 0.57d0*sqrt(Nue)	    + 0.61d0* Nue    ) * (1.d0 + 0.22d0* Nue	*X**3.d0)   )**2.d0 &
	    - 2.5d0*dC1_rho
  
  dC3_Te   =-(7.d0 + 6.5d0*X) * ( (       0.28d0/sqrt(Nue)*dNue_Te  + 0.61d0*dNue_Te ) * (1.d0 + 0.22d0* Nue    *X**3.d0)     &
           			+ (1.d0 + 0.57d0*sqrt(Nue)	    + 0.61d0* Nue    ) * (	 0.22d0*dNue_Te *X**3.d0)   ) &
           		      / ( (1.d0 + 0.57d0*sqrt(Nue)	    + 0.61d0* Nue    ) * (1.d0 + 0.22d0* Nue	*X**3.d0)   )**2.d0 &
	    - 2.5d0*dC1_Te
  

  ! --- **************************************************** ---
  ! --- Coefficient AC4
  AC4      =     (  -1.17d0/(1.d0+0.46d0*X)       + 0.35d0*sqrt(Nui)         ) / (1.d0 + 0.7d0 *sqrt(Nui)         ) &
            + 0.26d0*Nui*Nui*X**6.d0

  dAC4_psi = (   (0.54d0*dX/(1.d0+0.46d0*X)**2.d0 + 0.18d0/sqrt(Nui)*dNui_psi) * (1.d0 + 0.7d0 *sqrt(Nui)         )   &
               - (  -1.17d0/(1.d0+0.46d0*X)	  + 0.35d0*sqrt(Nui)	     ) * (	 0.25d0/sqrt(Nui)*dNui_psi) ) &
	     / ( 1.d0 + 0.7d0 *sqrt(Nui) )**2.d0								      &
            + 1.56d0*Nui*Nui*X**5.d0*dX + 0.52d0*Nui*dNui_psi*X**6.d0

  dAC4_rho = (   (                                  0.18d0/sqrt(Nui)*dNui_rho) * (1.d0 + 0.7d0 *sqrt(Nui)         )   &
               - (  -1.17d0/(1.d0+0.46d0*X)	  + 0.35d0*sqrt(Nui)	     ) * (	 0.35d0/sqrt(Nui)*dNui_rho) ) &
	     / ( 1.d0 + 0.7d0 *sqrt(Nui) )**2.d0								      &
            + 0.52d0*Nui*dNui_rho*X**6.d0
  
  dAC4_Ti  = (   (                                  0.18d0/sqrt(Nui)*dNui_Ti ) * (1.d0 + 0.7d0 *sqrt(Nui)         )   &
               - (  -1.17d0/(1.d0+0.46d0*X)	  + 0.35d0*sqrt(Nui)	     ) * (	 0.35d0/sqrt(Nui)*dNui_Ti ) ) &
	     / ( 1.d0 + 0.7d0 *sqrt(Nui) )**2.d0								      &
            + 0.52d0*Nui*dNui_Ti*X**6.d0
  

  ! --- **************************************************** ---
  ! --- Coefficient BC4
  BC4      =  (1.d0 - 0.125d0*Nui*Nui     *X**6.d0   ) * (1.d0 + 0.125d0*Nue*Nue     *X**6.d0)

  dBC4_psi =  (       0.25d0 *Nui*dNui_psi*X**6.d0   ) * (1.d0 + 0.125d0*Nue*Nue     *X**6.d0   ) &
            + (       0.75d0 *Nui*Nui	  *X**5.d0*dX) * (1.d0 + 0.125d0*Nue*Nue     *X**6.d0	) &
	    + (1.d0 - 0.125d0*Nui*Nui	  *X**6.d0   ) * (	 0.25d0 *Nue*dNue_psi*X**6.d0	) &
	    + (1.d0 - 0.125d0*Nui*Nui	  *X**6.d0   ) * (	 0.75d0 *Nue*Nue     *X**5.d0*dX)

  dBC4_rho =  (       0.25d0 *Nui*dNui_rho*X**6.d0   ) * (1.d0 + 0.125d0*Nue*Nue     *X**6.d0) &
            + (1.d0 - 0.125d0*Nui*Nui     *X**6.d0   ) * (       0.25d0 *Nue*dNue_rho*X**6.d0)

  dBC4_Ti  =  (       0.25d0 *Nui*dNui_Ti *X**6.d0   ) * (1.d0 + 0.125d0*Nue*Nue     *X**6.d0)

  dBC4_Te  =  (1.d0 - 0.125d0*Nui*Nui     *X**6.d0   ) * (       0.25d0 *Nue*dNue_Te *X**6.d0)

  ! --- **************************************************** ---
  ! --- Coefficient C4
  C4       = AC4*C2 / BC4
  dC4_psi  = ( (dAC4_psi*C2 + AC4*dC2_psi)*BC4 - AC4*C2*dBC4_psi) / BC4**2.d0
  dC4_rho  = ( (dAC4_rho*C2 + AC4*dC2_rho)*BC4 - AC4*C2*dBC4_rho) / BC4**2.d0
  dC4_Ti   = ( (dAC4_Ti *C2 + AC4*dC2_Ti )*BC4 - AC4*C2*dBC4_Ti ) / BC4**2.d0
  dC4_Te   = ( (              AC4*dC2_Te )*BC4 - AC4*C2*dBC4_Te ) / BC4**2.d0


  ! --- **************************************************** ---
  ! --- Jb's
  Jb1       =  C1     *(dTe/Te + drho/rho)
  dJb1_psi  =  dC1_psi*(dTe/Te + drho/rho) + C1*(dTe_psi/Te + drho_psi/rho)
  dJb1_rho  =  dC1_rho*(dTe/Te + drho/rho) + C1*(             drho_rho/rho - drho*vrho/rho**2.d0)
  dJb1_Te   =  dC1_Te *(dTe/Te + drho/rho) + C1*(dTe_Te/Te - dTe*vTe/Te**2.d0)

  Jb2       =  C2     *(dTi/Ti + drho/rho)
  dJb2_psi  =  dC2_psi*(dTi/Ti + drho/rho) + C2*(dTi_psi/Ti + drho_psi/rho)
  dJb2_rho  =  dC2_rho*(dTi/Ti + drho/rho) + C2*(             drho_rho/rho - drho*vrho/rho**2.d0)
  dJb2_Ti   =  dC2_Ti *(dTi/Ti + drho/rho) + C2*(dTi_Ti/Ti - dTi*vTi/Ti**2.d0)
  dJb2_Te   =  dC2_Te *(dTi/Ti + drho/rho)

  Jb3       =  C3     *dTe/Te
  dJb3_psi  =  dC3_psi*dTe/Te + C3*dTe_psi/Te
  dJb3_rho  =  dC3_rho*dTe/Te
  dJb3_Te   =  dC3_Te *dTe/Te + C3*dTe_Te/Te - C3*dTe*vTe/Te**2.d0

  Jb4       =  C4     *dTi/Ti
  dJb4_psi  =  dC4_psi*dTi/Ti + C4*dTi_psi/Ti
  dJb4_rho  =  dC4_rho*dTi/Ti
  dJb4_Ti   =  dC4_Ti *dTi/Ti + C4*dTi_Ti/Ti - C4*dTi*vTi/Ti**2.d0
  dJb4_Te   =  dC4_Te *dTi/Ti


  ! --- **************************************************** ---
  ! --- Bootstrap Current
  Jb       =  BigR* X*rho*Te     /DD       * ( Jb1      + Jb2      + Jb3      + Jb4     )

  dJb_psi  =  BigR*dX*rho *Te    /DD       * ( Jb1      + Jb2      + Jb3      + Jb4     ) &
            - BigR* X*rho *Te*dDD/DD**2.d0 * ( Jb1	+ Jb2	   + Jb3      + Jb4	) &
            + BigR* X*rho *Te	 /DD	   * ( dJb1_psi + dJb2_psi + dJb3_psi + dJb4_psi)
 
  dJb_rho  =  BigR* X*vrho*Te    /DD       * ( Jb1      + Jb2      + Jb3      + Jb4     ) &
            + BigR* X*rho *Te    /DD       * ( dJb1_rho + dJb2_rho + dJb3_rho + dJb4_rho)

  dJb_Ti   =  BigR* X*rho *Te    /DD       * (            dJb2_Ti             + dJb4_Ti)

  dJb_Te   =  BigR* X*rho *vTe   /DD       * ( Jb1      + Jb2      + Jb3      + Jb4     ) &
            + BigR* X*rho *Te    /DD       * ( dJb1_Te  + dJb2_Te  + dJb3_Te  + dJb4_Te)


  ! --- **************************************************** ---
  ! --- Current with denormalisation
  dJb_psi = dJb_psi / MU_ZERO
  dJb_rho = dJb_rho / MU_ZERO
  dJb_Ti  = dJb_Ti  / MU_ZERO
  dJb_Te  = dJb_Te  / MU_ZERO




return
end subroutine bootstrap_current_lhs
