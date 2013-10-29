!> Determine the heat and particle sources at a given position.
subroutine sources(xpoint2, xcase2, Z, Z_xpoint, psi, psi_axis, psi_bnd, particle_source, heat_source_i,   &
  heat_source_e)

  use phys_module

  implicit none

  ! --- Routine parameters.
  logical, intent(in)	:: xpoint2
  integer, intent(in)	:: xcase2
  real*8,  intent(in)	:: Z
  real*8,  intent(in)	:: Z_xpoint(2)
  real*8,  intent(in)	:: psi
  real*8,  intent(in)	:: psi_axis
  real*8,  intent(in)	:: psi_bnd
  real*8,  intent(out)  :: particle_source
  real*8,  intent(out)  :: heat_source_i
  real*8,  intent(out)  :: heat_source_e

  ! --- Local variables
  real*8 :: psi_n, sig_Ti, sig_Te
  
  sig_Ti = 0.01
  sig_Te = 0.01

  psi_n = (psi - psi_axis) / (psi_bnd - psi_axis)
  
  particle_source = particlesource * (0.5d0 - 0.5d0*tanh((psi_n - particlesource_psin)/particlesource_sig))
  if(xcase2 .eq. 1) then
    heat_source_i   = heatsource_i   * (0.5d0 - 0.5d0*tanh((psi_n - 0.8d0)/sig_Ti )) * (0.5d0 + 0.5d0*tanh((Z - Z_xpoint(1))/0.01))
    heat_source_e   = heatsource_e   * (0.5d0 - 0.5d0*tanh((psi_n - 0.8d0)/sig_Te )) * (0.5d0 + 0.5d0*tanh((Z - Z_xpoint(1))/0.01))
  endif
  if(xcase2 .eq. 2) then
    heat_source_i   = heatsource_i   * (0.5d0 - 0.5d0*tanh((psi_n - 0.8d0)/sig_Ti )) * (0.5d0 + 0.5d0*tanh((Z_xpoint(2) - Z)/0.01))
    heat_source_e   = heatsource_e   * (0.5d0 - 0.5d0*tanh((psi_n - 0.8d0)/sig_Te )) * (0.5d0 + 0.5d0*tanh((Z_xpoint(2) - Z)/0.01))
  endif
  if(xcase2 .eq. 3) then
    heat_source_i   = heatsource_i   * (0.5d0 - 0.5d0*tanh((psi_n - 0.8d0)/sig_Ti )) * (0.5d0 + 0.5d0*tanh((Z - Z_xpoint(1))/0.01)) * (0.5d0 + 0.5d0*tanh((Z_xpoint(2) - Z)/0.01))
    heat_source_e   = heatsource_e   * (0.5d0 - 0.5d0*tanh((psi_n - 0.8d0)/sig_Te )) * (0.5d0 + 0.5d0*tanh((Z - Z_xpoint(1))/0.01)) * (0.5d0 + 0.5d0*tanh((Z_xpoint(2) - Z)/0.01))
  endif

  return
end subroutine sources














subroutine velocity(xpoint2,xcase2,Z,Z_xpoint,psi,psi_axis,psi_bnd,velocity_profile)
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  use phys_module
  
  implicit none
  
  ! --- Routine parameters.
  logical, intent(in)	:: xpoint2
  integer, intent(in)	:: xcase2
  real*8,  intent(in)	:: Z
  real*8,  intent(in)	:: Z_xpoint(2)
  real*8,  intent(in)	:: psi
  real*8,  intent(in)	:: psi_axis
  real*8,  intent(in)	:: psi_bnd
  real*8,  intent(out)  :: velocity_profile
  
  real*8  :: prof0, prof1, dprof0_dpsi, dprof0_dpsi2, dprof0_dpsi3, psi_barrier
  real*8  :: psi_n, sig_n, sigz, dprof1_dpsi, dprof1_dpsi2, dprof1_dpsi3
  real*8  :: atn, datn, d2atn, d3atn, atn_z, datn_z, d2atn_z
  real*8  :: atn_z_u, datn_z_u, d2atn_z_u, Z_star, Z_star_u
  real*8  :: cosh3, cosh3_u, tanh2, tanh2_u
  
  sig_n       = V_coef(4)
  psi_barrier = V_coef(5)
  psi_n       = (psi - psi_axis) / (psi_bnd - psi_axis)
  
  prof0 = (V_0-V_1)*(1.d0 +V_coef(1)*psi_n+ V_coef(2)*psi_n**2+ V_coef(3) * psi_n**2)
  atn	= (0.5d0 - 0.5d0*tanh((psi_n - psi_barrier)/sig_n))
  prof1 = prof0        * atn

  velocity_profile = prof1

  if (xpoint2) then

    sigz	    = 0.05d0
  
    Z_star   = (Z_xpoint(1)-Z)/sigz
    Z_star   = min( max( Z_star, -40.d0), 40.d0) ! avoid floating-point exceptions
    Z_star_u = (Z-Z_xpoint(2))/sigz
    Z_star_u = min( max( Z_star_u, -40.d0), 40.d0) ! avoid floating-point exceptions
  
    tanh2   = tanh(Z_star)
    tanh2_u = tanh(Z_star_u)
      
    atn_z	     = (0.5d0 - 0.5d0*tanh2)
    atn_z_u	     = (0.5d0 - 0.5d0*tanh2_u)
    
    if(xcase2 .eq. 1) then
      atn_z_u	       = 1.d0
    endif
    if(xcase2 .eq. 2) then
      atn_z	       = 1.d0
    endif
  
    velocity_profile =   prof1 * atn_z * atn_z_u
  
  else
    
    velocity_profile = prof1
  
  endif
  
  velocity_profile = velocity_profile + V_1

  return
end
