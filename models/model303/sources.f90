!> Determine the heat and particle sources at a given position.
subroutine sources(xpoint2, xcase2, Z, Z_xpoint, psi, psi_axis, psi_bnd, particle_source, heat_source)

use phys_module

implicit none

! --- Routine parameters.
logical, intent(in)   :: xpoint2
integer, intent(in)   :: xcase2
real*8,  intent(in)   :: Z
real*8,  intent(in)   :: Z_xpoint(2)
real*8,  intent(in)   :: psi
real*8,  intent(in)   :: psi_axis
real*8,  intent(in)   :: psi_bnd
real*8,  intent(out)  :: particle_source
real*8,  intent(out)  :: heat_source

! --- Local variables
real*8 :: psi_n

psi_n = (psi - psi_axis) / (psi_bnd - psi_axis)

particle_source = particlesource * (0.5d0 - 0.5d0*tanh((psi_n - particlesource_psin)/particlesource_sig))
heat_source     = heatsource     * (0.5d0 - 0.5d0*tanh((psi_n - heatsource_psin    )/heatsource_sig    ))

return
end subroutine sources

!====MB===============parallel velocity profile which is kept by the // velocity source implemented in element_matrix.f90

subroutine velocity(xpoint2,Z,Z_xpoint,psi,psi_axis,psi_bnd,velocity_profile,dV_dpsi,dV_dz, &
                   dV_dpsi2,dV_dz2,dV_dpsi_dz,dV_dpsi3,dV_dpsi_dz2, dV_dpsi2_dz)
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
use phys_module

implicit none

logical :: xpoint2
real*8  :: prof0, prof1, dprof0_dpsi, dprof0_dpsi2, dprof0_dpsi3, velocity_profile, psi_barrier
real*8  :: dV_dpsi, dV_dz, dV_dpsi2, dV_dz2, dV_dpsi_dz, dV_dpsi3, dV_dpsi_dz2, dV_dpsi2_dz
real*8  :: Z, Z_xpoint, psi, psi_axis, psi_bnd,  psi_n, sig_n, sigz, dprof1_dpsi, dprof1_dpsi2, dprof1_dpsi3
real*8  :: atn, datn, d2atn, d3atn, atn_z, datn_z, d2atn_z, factor

sig_n       = V_coef(4)
psi_barrier = V_coef(5)

psi_n = (psi - psi_axis) / (psi_bnd - psi_axis)

factor = 1.d0


prof0        = (V_0-V_1)*(1.d0 +V_coef(1)*psi_n+ V_coef(2)*psi_n**2+ V_coef(3) * psi_n**2)
dprof0_dpsi  = (V_0-V_1)*(V_coef(1) + 2.d0 * V_coef(2) * psi_n + 3.d0 * V_coef(3) * psi_n**2) / (psi_bnd - psi_axis)
dprof0_dpsi2 = (V_0-V_1)*(2.d0 * V_coef(2) + 6.d0 * V_coef(3) * psi_n)                          / (psi_bnd - psi_axis)**2
dprof0_dpsi3 = (V_0-V_1)*(6.d0 * V_coef(3))                                                       / (psi_bnd - psi_axis)**3

atn   = (0.5d0 - 0.5d0*tanh((psi_n - psi_barrier)/sig_n))

datn  = - 1.d0/cosh((psi_n - psi_barrier)/sig_n)**2 / (2.d0 * sig_n) / (psi_bnd - psi_axis)

d2atn =   1.d0/cosh((psi_n - psi_barrier)/sig_n)**2 / (sig_n**2)  &
      * tanh((psi_n - psi_barrier)/sig_n) / (psi_bnd - psi_axis)**2

d3atn = - 1.d0/cosh((psi_n - psi_barrier)/sig_n)**4 / (sig_n**3)  &
      * (-2.d0 + cosh(2.d0*(psi_n-psi_barrier)/sig_n) ) / (psi_bnd - psi_axis)**3

prof1        = prof0        * atn
dprof1_dpsi  = dprof0_dpsi  * atn +         prof0       * datn
dprof1_dpsi2 = dprof0_dpsi2 * atn + 2.d0 * dprof0_dpsi  * datn + prof0              * d2atn
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
  sigz    = 0.1d0

  atn_z           =  (0.5d0 - 0.5d0*tanh((Z_xpoint-Z)/sigz))
  datn_z          =  0.5d0/cosh((Z_xpoint-Z)/sigz)**2 / sigz
  d2atn_z         = 1.d0/cosh((Z_xpoint-Z)/sigz)**2 /  sigz**2 * tanh((Z_xpoint-Z)/sigz)

  velocity_profile =   prof1        * atn_z
  dV_dpsi         =   dprof1_dpsi  * atn_z
  dV_dpsi2        =   dprof1_dpsi2 * atn_z
  dV_dpsi3        =   dprof1_dpsi3 * atn_z
  dV_dz           = + prof1        * datn_z
  dV_dz2          = + prof1        * d2atn_z
  dV_dpsi_dz      =   dprof1_dpsi  * datn_z
  dV_dpsi2_dz     =   dprof1_dpsi2 * datn_z
  dV_dpsi_dz2     =   dprof1_dpsi  * d2atn_z

endif

velocity_profile = velocity_profile + V_1

return
end
!============================================Marina 14.02.2011================
