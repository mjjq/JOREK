subroutine sources(xpoint2,Z,Z_xpoint,psi,psi_axis,psi_bnd,particle_source,heat_source_i,heat_source_e)
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
use phys_module

implicit none

logical :: xpoint2
real*8  :: Z, Z_xpoint, psi, psi_axis, psi_bnd, psi_n, sig_n, sig_Ti, sig_Te
real*8  :: particle_source, heat_source_i, heat_source_e

sig_n  = 0.01
sig_Ti = 0.01
sig_Te = 0.01

psi_n = (psi - psi_axis) / (psi_bnd - psi_axis)

particle_source = particlesource * (0.5d0 - 0.5d0*tanh((psi_n - 0.8d0)/sig_n  )) * (0.5d0 + 0.5d0*tanh((Z - Z_xpoint)/0.01))
heat_source_i   = heatsource_i   * (0.5d0 - 0.5d0*tanh((psi_n - 0.8d0)/sig_Ti )) * (0.5d0 + 0.5d0*tanh((Z - Z_xpoint)/0.01))
heat_source_e   = heatsource_e   * (0.5d0 - 0.5d0*tanh((psi_n - 0.8d0)/sig_Te )) * (0.5d0 + 0.5d0*tanh((Z - Z_xpoint)/0.01))

return
end
