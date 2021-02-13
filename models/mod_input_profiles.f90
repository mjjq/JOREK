!> module with generic functions to build input profiles, in order to avoid duplication
module mod_input_profiles

contains


subroutine input_profiles_psi_component(psi,psi_axis,psi_bnd, &
                                        CORE, SOL, coef, &
                                        prof1,           &
                                        dprof1_dpsi,     &
                                        dprof1_dpsi2,    &
                                        dprof1_dpsi3,    &
                                        dprof1_dpsi4,    &
                                        dprof1_dpsi5)
  !-----------------------------------------------------------------------
  ! Analytical profile in the core of the plasma, 3rd order polynomial with tanh
  !-----------------------------------------------------------------------
  use phys_module
  
  implicit none
  
  ! --- Routine parameters
  real*8,  intent(in)  :: psi, psi_axis, psi_bnd
  real*8,  intent(in)  :: CORE, SOL, coef(10)
  real*8,  intent(out) :: prof1
  real*8,  intent(out) :: dprof1_dpsi
  real*8,  intent(out) :: dprof1_dpsi2
  real*8,  intent(out) :: dprof1_dpsi3
  real*8,  intent(out) :: dprof1_dpsi4
  real*8,  intent(out) :: dprof1_dpsi5
  
  ! --- Internal variables.
  real*8  :: psi_n, psi_star, delta_psi, sigg, psi_barrier
  real*8  :: atn, datn, d2atn, d3atn, d4atn, d5atn
  real*8  :: cosh1, sinh1, tanh1
  real*8  :: prof0       
  real*8  :: dprof0_dpsi 
  real*8  :: dprof0_dpsi2
  real*8  :: dprof0_dpsi3
  real*8  :: dprof0_dpsi4
  real*8  :: dprof0_dpsi5
  
  delta_psi = psi_bnd - psi_axis
  psi_n     = (psi - psi_axis) / delta_psi
  psi_n     = max( min(psi_n, 2.), 0. )

  prof0        = (CORE-SOL)*(1.d0 + coef(1) * psi_n        + coef(2) * psi_n**2        + coef(3) * psi_n**3)
  dprof0_dpsi  = (CORE-SOL)*(       coef(1)         + 2.d0 * coef(2) * psi_n    + 3.d0 * coef(3) * psi_n**2) / delta_psi
  dprof0_dpsi2 = (CORE-SOL)*(                         2.d0 * coef(2)            + 6.d0 * coef(3) * psi_n   ) / delta_psi**2
  dprof0_dpsi3 = (CORE-SOL)*(                                                     6.d0 * coef(3)           ) / delta_psi**3
  dprof0_dpsi4 = 0.d0
  dprof0_dpsi5 = 0.d0
  
  sigg        = coef(4)
  psi_barrier = coef(5)
  
  psi_star = (psi_n - psi_barrier)/sigg
  psi_star = min( max( psi_star, -40.d0), 40.d0) ! avoid floating-point exceptions
  
  tanh1 = tanh(psi_star)
  cosh1 = cosh(psi_star)
  sinh1 = sinh(psi_star)
  
  atn   = (0.5d0 - 0.5d0*tanh1)
  datn  = - 0.5d0 * cosh1**(-2)            / (sigg*delta_psi)
  d2atn = + 1.0d0 * cosh1**(-3) * sinh1    / (sigg*delta_psi)**2
  d3atn = - 3.0d0 * cosh1**(-4) * sinh1**2 / (sigg*delta_psi)**3 &
          + 1.0d0 * cosh1**(-2)            / (sigg*delta_psi)**3
  d4atn = +12.0d0 * cosh1**(-5) * sinh1**3 / (sigg*delta_psi)**4 &
          - 8.0d0 * cosh1**(-3) * sinh1    / (sigg*delta_psi)**4
  d5atn = -60.0d0 * cosh1**(-6) * sinh1**4 / (sigg*delta_psi)**5 &
          +60.0d0 * cosh1**(-4) * sinh1**2 / (sigg*delta_psi)**5 &
          - 8.0d0 * cosh1**(-2)            / (sigg*delta_psi)**5 
  
  prof1        = prof0        * atn
  dprof1_dpsi  = dprof0_dpsi  * atn +         prof0       * datn
  dprof1_dpsi2 = dprof0_dpsi2 * atn + 2.d0 * dprof0_dpsi  * datn + prof0              * d2atn
  dprof1_dpsi3 = dprof0_dpsi3 * atn + 3.d0 * dprof0_dpsi2 * datn + 3.d0 * dprof0_dpsi * d2atn + prof0 * d3atn
  dprof1_dpsi4 = dprof0_dpsi4 * atn + 4.d0 * dprof0_dpsi3 * datn + 6.d0 * dprof0_dpsi2 * d2atn &
                 + 4.d0 * dprof0_dpsi * d3atn + prof0 * d4atn
  dprof1_dpsi5 = dprof0_dpsi5 * atn + 5.d0 * dprof0_dpsi4 * datn + 10.d0 * dprof0_dpsi3 * d2atn &
                + 10.d0 * dprof0_dpsi2 * d3atn + 5.d0 * dprof0_dpsi * d4atn + prof0 * d5atn
  
  return
end subroutine input_profiles_psi_component








subroutine input_profiles_edge_perturbation(psi,psi_axis,psi_bnd, &
                                            coef,    &
                                            d_pert,  &
                                            d2_pert, &
                                            d3_pert, &
                                            d4_pert, &
                                            d5_pert)
  !-----------------------------------------------------------------------
  ! Analytical profile for FFprime perturbation at the plasma edge (basically the gradient of a tanh)
  !-----------------------------------------------------------------------
  use phys_module
  
  implicit none
  
  ! --- Routine parameters
  real*8,  intent(in)  :: psi, psi_axis, psi_bnd
  real*8,  intent(in)  :: coef(10)
  real*8,  intent(out) :: d_pert
  real*8,  intent(out) :: d2_pert
  real*8,  intent(out) :: d3_pert
  real*8,  intent(out) :: d4_pert
  real*8,  intent(out) :: d5_pert
  
  ! --- Internal variables.
  real*8  :: psi_n, psi_star, delta_psi, no_delta_psi, sigg, psi_barrier
  real*8  :: atn, datn, d2atn, d3atn, d4atn
  real*8  :: cosh1, sinh1, tanh1
  real*8  :: prof0       
  real*8  :: dprof0_dpsi 
  real*8  :: dprof0_dpsi2
  real*8  :: dprof0_dpsi3
  real*8  :: dprof0_dpsi4
  
  delta_psi = psi_bnd - psi_axis
  psi_n     = (psi - psi_axis) / delta_psi
  psi_n     = max( min(psi_n, 2.), 0. )
  no_delta_psi = 1.d0
  if (coef(9) .eq. 1.d0) no_delta_psi = delta_psi

  psi_star = (psi_n - coef(7))/coef(8)
  cosh1    = cosh(psi_star)
  sinh1    = sinh(psi_star)
  
  prof0         = + 0.5d0 * coef(6) * cosh1**(-2)            / (coef(8)*delta_psi)    * no_delta_psi
  dprof0_dpsi   = - 1.0d0 * coef(6) * cosh1**(-3) * sinh1    / (coef(8)*delta_psi)**2 * no_delta_psi
  dprof0_dpsi2  = + 3.0d0 * coef(6) * cosh1**(-4) * sinh1**2 / (coef(8)*delta_psi)**3 * no_delta_psi &
                  - 1.0d0 * coef(6) * cosh1**(-2)            / (coef(8)*delta_psi)**3 * no_delta_psi
  dprof0_dpsi3  = -12.0d0 * coef(6) * cosh1**(-5) * sinh1**3 / (coef(8)*delta_psi)**4 * no_delta_psi &
                  + 8.0d0 * coef(6) * cosh1**(-3) * sinh1    / (coef(8)*delta_psi)**4 * no_delta_psi
  dprof0_dpsi4  = +60.0d0 * coef(6) * cosh1**(-6) * sinh1**4 / (coef(8)*delta_psi)**5 * no_delta_psi &
                  -60.0d0 * coef(6) * cosh1**(-4) * sinh1**2 / (coef(8)*delta_psi)**5 * no_delta_psi &
                  + 8.0d0 * coef(6) * cosh1**(-2)            / (coef(8)*delta_psi)**5 * no_delta_psi
  
  sigg        = coef(4)
  psi_barrier = coef(5)
  
  psi_star = (psi_n - psi_barrier)/sigg
  psi_star = min( max( psi_star, -40.d0), 40.d0) ! avoid floating-point exceptions
  
  tanh1 = tanh(psi_star)
  cosh1 = cosh(psi_star)
  sinh1 = sinh(psi_star)
  
  atn   = (0.5d0 - 0.5d0*tanh1)
  datn  = - 0.5d0 * cosh1**(-2)            / (sigg*delta_psi)
  d2atn = + 1.0d0 * cosh1**(-3) * sinh1    / (sigg*delta_psi)**2
  d3atn = - 3.0d0 * cosh1**(-4) * sinh1**2 / (sigg*delta_psi)**3 &
          + 1.0d0 * cosh1**(-2)            / (sigg*delta_psi)**3
  d4atn = +12.0d0 * cosh1**(-5) * sinh1**3 / (sigg*delta_psi)**4 &
          - 8.0d0 * cosh1**(-3) * sinh1    / (sigg*delta_psi)**4
  
  d_pert  = prof0        * atn
  d2_pert = dprof0_dpsi  * atn +         prof0       * datn
  d3_pert = dprof0_dpsi2 * atn + 2.d0 * dprof0_dpsi  * datn + prof0              * d2atn
  d4_pert = dprof0_dpsi3 * atn + 3.d0 * dprof0_dpsi2 * datn + 3.d0 * dprof0_dpsi * d2atn + prof0 * d3atn
  d5_pert = dprof0_dpsi4 * atn + 4.d0 * dprof0_dpsi3 * datn + 6.d0 * dprof0_dpsi2 * d2atn &
                 + 4.d0 * dprof0_dpsi * d3atn + prof0 * d4atn
  
  return
end subroutine input_profiles_edge_perturbation








subroutine input_profiles_Z_component(xpoint2,xcase2,Z,Z_xpoint,&
                                      atn_both,   &
                                      datn_both,  &
                                      d2atn_both, &
                                      d3atn_both, &
                                      d4atn_both)  
  !-----------------------------------------------------------------------
  ! tanh profile in Z for the X-point private regions
  !-----------------------------------------------------------------------
  use phys_module
 
  implicit none
  
  ! --- Routine parameters
  logical, intent(in)  :: xpoint2
  integer, intent(in)  :: xcase2
  real*8,  intent(in)  :: Z, Z_xpoint(2)
  real*8,  intent(out) :: atn_both
  real*8,  intent(out) :: datn_both
  real*8,  intent(out) :: d2atn_both
  real*8,  intent(out) :: d3atn_both
  real*8,  intent(out) :: d4atn_both
 
  ! --- Internal variables.
  real*8  :: sigz
  real*8  :: atn,     datn,     d2atn,     d3atn,     d4atn
  real*8  :: atn_z,   datn_z,   d2atn_z,   d3atn_z,   d4atn_z
  real*8  :: atn_z_u, datn_z_u, d2atn_z_u, d3atn_z_u, d4atn_z_u
  real*8  :: cosh1, sinh1, tanh1
  real*8  :: Z_star, Z_star_u

  sigz            = 0.1d0

  if (xpoint2) then

    if (xcase2 .eq. 1) then
      atn_z_u   = 1.d0
      datn_z_u  = 0.d0
      d2atn_z_u = 0.d0
      d3atn_z_u = 0.d0
      d4atn_z_u = 0.d0
    else
      Z_star_u  = (Z-Z_xpoint(2))/sigz
      Z_star_u  = min( max( Z_star_u, -40.d0), 40.d0) ! avoid floating-point exceptions
      
      tanh1     = tanh(Z_star_u)
      cosh1     = cosh(Z_star_u)
      sinh1     = sinh(Z_star_u)
 
      atn_z_u   = (0.5d0 - 0.5d0*tanh1)
      datn_z_u  = - 0.5d0 * cosh1**(-2)            / sigz
      d2atn_z_u = + 1.0d0 * cosh1**(-3) * sinh1    / sigz**2
      d3atn_z_u = - 3.0d0 * cosh1**(-4) * sinh1**2 / sigz**3 &
                  + 1.0d0 * cosh1**(-2)            / sigz**3
      d4atn_z_u = +12.0d0 * cosh1**(-5) * sinh1**3 / sigz**4 &
                  - 8.0d0 * cosh1**(-3) * sinh1    / sigz**4
     
    endif
    
    if (xcase2 .eq. 2) then
      atn_z   = 1.d0
      datn_z  = 0.d0
      d2atn_z = 0.d0
      d3atn_z = 0.d0
      d4atn_z = 0.d0
    else
      Z_star  = (Z_xpoint(1)-Z)/sigz
      Z_star  = min( max( Z_star, -40.d0), 40.d0) ! avoid floating-point exceptions
 
      tanh1   = tanh(Z_star)
      cosh1   = cosh(Z_star)
      sinh1   = sinh(Z_star)
  
      atn_z   = (0.5d0 - 0.5d0*tanh1)
      datn_z  = - 0.5d0 * cosh1**(-2)            / sigz
      d2atn_z = + 1.0d0 * cosh1**(-3) * sinh1    / sigz**2
      d3atn_z = - 3.0d0 * cosh1**(-4) * sinh1**2 / sigz**3 &
                + 1.0d0 * cosh1**(-2)            / sigz**3
      d4atn_z = +12.0d0 * cosh1**(-5) * sinh1**3 / sigz**4 &
                - 8.0d0 * cosh1**(-3) * sinh1    / sigz**4
       
    endif
    
    atn_both   = atn_z * atn_z_u
    datn_both  = datn_z * atn_z_u + atn_z * datn_z_u
    d2atn_both = d2atn_z * atn_z_u + 2.d0 * datn_z * datn_z_u + atn_z * d2atn_z_u
    d3atn_both = d3atn_z * atn_z_u + 3.d0 * d2atn_z * datn_z_u + 3.d0 *  datn_z * d2atn_z_u + atn_z * d3atn_z_u
    d4atn_both = d4atn_z * atn_z_u + 4.d0 * d3atn_z * datn_z_u + 6.d0 * d2atn_z * d2atn_z_u + 4.d0 * datn_z * d3atn_z_u + atn_z * d4atn_z_u

  else

    atn_both   = 1.d0
    datn_both  = 0.d0
    d2atn_both = 0.d0
    d3atn_both = 0.d0
    d4atn_both = 0.d0

  endif


return
end subroutine input_profiles_Z_component








end module mod_input_profiles
