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
  
!  sig_Ti = 0.01
!  sig_Te = 0.01

  psi_n = (psi - psi_axis) / (psi_bnd - psi_axis)
  
  if (xpoint2) then
    if ((Z .lt. Z_xpoint(1)) .and. (psi_n .lt. 1.d0) ) then
      psi_n = 2.d0 - psi_n
    endif
  endif
  
  particle_source = particlesource * (0.5d0 - 0.5d0*tanh((psi_n - particlesource_psin)/particlesource_sig))   &
      + edgeparticlesource * (0.5d0 + 0.5d0*tanh((psi_n - edgeparticlesource_psin)/edgeparticlesource_sig))   &
      + particlesource_gauss * exp(-(psi_n - particlesource_gauss_psin)**2/(particlesource_gauss_sig**2))
  heat_source_i     = heatsource_i     * (0.5d0 - 0.5d0*tanh((psi_n - heatsource_psin    )/heatsource_sig    ))  &
      + heatsource_gauss * exp(-(psi_n - heatsource_gauss_psin)**2/(heatsource_gauss_sig**2))
  heat_source_e     = heatsource_e     * (0.5d0 - 0.5d0*tanh((psi_n - heatsource_psin    )/heatsource_sig    ))  &
      + heatsource_gauss * exp(-(psi_n - heatsource_gauss_psin)**2/(heatsource_gauss_sig**2))

!  if(xcase2 .eq. 1) then
!    heat_source_i   = heatsource_i   * (0.5d0 - 0.5d0*tanh((psi_n - 0.8d0)/sig_Ti )) * (0.5d0 + 0.5d0*tanh((Z - Z_xpoint(1))/0.01))
!    heat_source_e   = heatsource_e   * (0.5d0 - 0.5d0*tanh((psi_n - 0.8d0)/sig_Te )) * (0.5d0 + 0.5d0*tanh((Z - Z_xpoint(1))/0.01))
!  endif
!  if(xcase2 .eq. 2) then
!    heat_source_i   = heatsource_i   * (0.5d0 - 0.5d0*tanh((psi_n - 0.8d0)/sig_Ti )) * (0.5d0 + 0.5d0*tanh((Z_xpoint(2) - Z)/0.01))
!    heat_source_e   = heatsource_e   * (0.5d0 - 0.5d0*tanh((psi_n - 0.8d0)/sig_Te )) * (0.5d0 + 0.5d0*tanh((Z_xpoint(2) - Z)/0.01))
!  endif
!  if(xcase2 .eq. 3) then
!    heat_source_i   = heatsource_i   * (0.5d0 - 0.5d0*tanh((psi_n - 0.8d0)/sig_Ti )) * (0.5d0 + 0.5d0*tanh((Z - Z_xpoint(1))/0.01)) * (0.5d0 + 0.5d0*tanh((Z_xpoint(2) - Z)/0.01))
!    heat_source_e   = heatsource_e   * (0.5d0 - 0.5d0*tanh((psi_n - 0.8d0)/sig_Te )) * (0.5d0 + 0.5d0*tanh((Z - Z_xpoint(1))/0.01)) * (0.5d0 + 0.5d0*tanh((Z_xpoint(2) - Z)/0.01))
!  endif

  return
end subroutine sources












subroutine velocity(xpoint2, xcase2,Z,Z_xpoint,psi,psi_axis,psi_bnd,velocity_profile,dV_dpsi,dV_dz, &
                   dV_dpsi2,dV_dz2,dV_dpsi_dz,dV_dpsi3,dV_dpsi_dz2, dV_dpsi2_dz)
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  use phys_module
  
  implicit none
  
  logical :: xpoint2
  integer, intent(in)	:: xcase2
  real*8  :: prof0, prof1, dprof0_dpsi, dprof0_dpsi2, dprof0_dpsi3, velocity_profile, psi_barrier
  real*8  :: dV_dpsi, dV_dz, dV_dpsi2, dV_dz2, dV_dpsi_dz, dV_dpsi3, dV_dpsi_dz2, dV_dpsi2_dz
  real*8  :: Z, Z_xpoint(2), psi, psi_axis, psi_bnd,  psi_n, sig_n, sigz, dprof1_dpsi, dprof1_dpsi2, dprof1_dpsi3
  real*8  :: atn, datn, d2atn, d3atn
  real*8  :: atn_z,   datn_z,	d2atn_z
  real*8  :: atn_z_u, datn_z_u, d2atn_z_u, factor
  real*8  :: cosh1, cosh2, cosh3, cosh3_u
  real*8  :: tanh1, tanh2, tanh2_u
  ! for interpolating numerical profiles
  integer :: left, right, mid
  real*8  :: aux1, aux2, Z_star, Z_star_u
  
  sig_n       = V_coef(4)
  psi_barrier = V_coef(5)

  psi_n = (psi - psi_axis) / (psi_bnd - psi_axis)
  
  factor = 1.d0
  
  
  prof0        = (V_0-V_1)*(1.d0 +V_coef(1)*psi_n+ V_coef(2)*psi_n**2+ V_coef(3) * psi_n**2)
  dprof0_dpsi  = (V_0-V_1)*(V_coef(1) + 2.d0 * V_coef(2) * psi_n + 3.d0 * V_coef(3) * psi_n**2) / (psi_bnd - psi_axis)
  dprof0_dpsi2 = (V_0-V_1)*(2.d0 * V_coef(2) + 6.d0 * V_coef(3) * psi_n)			  / (psi_bnd - psi_axis)**2
  dprof0_dpsi3 = (V_0-V_1)*(6.d0 * V_coef(3))							    / (psi_bnd - psi_axis)**3
  
  atn	= (0.5d0 - 0.5d0*tanh((psi_n - psi_barrier)/sig_n))
  
  datn  = - 1.d0/cosh((psi_n - psi_barrier)/sig_n)**2 / (2.d0 * sig_n) / (psi_bnd - psi_axis)
  
  d2atn =   1.d0/cosh((psi_n - psi_barrier)/sig_n)**2 / (sig_n**2)  &
  	* tanh((psi_n - psi_barrier)/sig_n) / (psi_bnd - psi_axis)**2
  
  d3atn = - 1.d0/cosh((psi_n - psi_barrier)/sig_n)**4 / (sig_n**3)  &
  	* (-2.d0 + cosh(2.d0*(psi_n-psi_barrier)/sig_n) ) / (psi_bnd - psi_axis)**3
  
  prof1        = prof0        * atn
  dprof1_dpsi  = dprof0_dpsi  * atn +	      prof0	  * datn
  dprof1_dpsi2 = dprof0_dpsi2 * atn + 2.d0 * dprof0_dpsi  * datn + prof0	      * d2atn
  dprof1_dpsi3 = dprof0_dpsi3 * atn + 3.d0 * dprof0_dpsi2 * datn + 3.d0 * dprof0_dpsi * d2atn + prof0 * d3atn
  
  dV_dpsi     = dprof1_dpsi   * factor
  dV_dpsi2    = dprof1_dpsi2
  dV_dpsi3    = dprof1_dpsi3  * factor
  
  dV_dz       = 0.d0
  dV_dz2      = 0.d0
  dV_dpsi_dz  = 0.d0
  dV_dpsi2_dz = 0.d0
  dV_dpsi_dz2 = 0.d0

  velocity_profile = prof1
  
  if (xpoint2) then
    sigz	    = 0.1d0
    
    Z_star   = (Z_xpoint(1)-Z)/sigz
    Z_star   = min( max( Z_star, -40.d0), 40.d0) ! avoid floating-point exceptions
    Z_star_u = (Z-Z_xpoint(2))/sigz
    Z_star_u = min( max( Z_star_u, -40.d0), 40.d0) ! avoid floating-point exceptions
  
    tanh2   = tanh(Z_star)
    cosh3   = cosh(Z_star)
    tanh2_u = tanh(Z_star_u)
    cosh3_u = cosh(Z_star_u)
      
    atn_z	     = (0.5d0 - 0.5d0*tanh2)
    datn_z	     =  0.5d0/cosh3**2   / sigz
    d2atn_z	     =  1.0d0/cosh3**2   / sigz**2 * tanh2
    atn_z_u	     = (0.5d0 - 0.5d0*tanh2_u)
    datn_z_u	     = -0.5d0/cosh3_u**2 / sigz
    d2atn_z_u	     =  1.0d0/cosh3_u**2 / sigz**2 * tanh2_u
    
    if(xcase2 .eq. 1) then
      atn_z_u	       = 1.d0
      datn_z_u         = 0.d0
      d2atn_z_u        = 0.d0
    endif
    if(xcase2 .eq. 2) then
      atn_z	       = 1.d0
      datn_z	       = 0.d0
      d2atn_z	       = 0.d0
    endif
    
    velocity_profile= prof1	     *    atn_z * atn_z_u
    dV_dpsi	    = dprof1_dpsi    *    atn_z * atn_z_u
    dV_dpsi2	    = dprof1_dpsi2   *    atn_z * atn_z_u
    dV_dpsi3	    = dprof1_dpsi3   *    atn_z * atn_z_u  
    dV_dz	    = prof1	     * ( datn_z * atn_z_u  +	      atn_z * datn_z_u)
    dV_dz2	    = prof1	     * (d2atn_z * atn_z_u  +  2.d0 * datn_z * datn_z_u  +  atn_z * d2atn_z_u)  
    dV_dpsi_dz      = dprof1_dpsi    * ( datn_z * atn_z_u  +	      atn_z * datn_z_u)
    dV_dpsi2_dz     = dprof1_dpsi    * (d2atn_z * atn_z_u  +  2.d0 * datn_z * datn_z_u  +  atn_z * d2atn_z_u)  
    dV_dpsi_dz2     = dprof1_dpsi2   * ( datn_z * atn_z_u  +	      atn_z * datn_z_u)
  
  endif
  
  velocity_profile = velocity_profile + V_1

return
end
  
