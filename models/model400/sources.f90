!> Determine the heat and particle sources at a given position.
subroutine sources(xpoint2, Z, Z_xpoint, psi, psi_axis, psi_bnd, particle_source, heat_source_i,   &
  heat_source_e)

use phys_module

implicit none

! --- Routine parameters.
logical, intent(in)   :: xpoint2
real*8,  intent(in)   :: Z
real*8,  intent(in)   :: Z_xpoint
real*8,  intent(in)   :: psi
real*8,  intent(in)   :: psi_axis
real*8,  intent(in)   :: psi_bnd
real*8,  intent(out)  :: particle_source
real*8,  intent(out)  :: heat_source_i
real*8,  intent(out)  :: heat_source_e

! --- Local variables
real*8 :: psi_n, sig_Ti, sig_Te

sig_Ti = 0.01
sig_Te = 0.01

psi_n = (psi - psi_axis) / (psi_bnd - psi_axis)

particle_source = particlesource * (0.5d0 - 0.5d0*tanh((psi_n - particlesource_psin)/particlesource_sig))
heat_source_i   = heatsource_i   * (0.5d0 - 0.5d0*tanh((psi_n - 0.8d0)/sig_Ti )) * (0.5d0 + 0.5d0*tanh((Z - Z_xpoint)/0.01))
heat_source_e   = heatsource_e   * (0.5d0 - 0.5d0*tanh((psi_n - 0.8d0)/sig_Te )) * (0.5d0 + 0.5d0*tanh((Z - Z_xpoint)/0.01))

return
end subroutine sources
