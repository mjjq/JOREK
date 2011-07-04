subroutine sources(xpoint2,Z,Z_xpoint,psi,psi_axis,psi_bnd,particle_source,heat_source)
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
use phys_module

implicit none

logical :: xpoint2
real*8  :: Z, Z_xpoint, psi, psi_axis, psi_bnd, psi_n, sig_n, sig_T
real*8  :: particle_source, heat_source

sig_n = 0.1
sig_T = 0.1

psi_n = (psi - psi_axis) / (psi_bnd - psi_axis)

particle_source = particlesource !* (0.5d0 - 0.5d0*tanh((psi_n - 1.d0)/sig_n))
heat_source     = heatsource     * (0.5d0 - 0.5d0*tanh((psi_n - 0.5d0)/sig_T))

return
end
