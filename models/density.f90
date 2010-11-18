subroutine density(xpoint2,Z,Z_xpoint,psi,psi_axis,psi_bnd,density_profile,dn_dpsi,dn_dz, &
                   dn_dpsi2,dn_dz2,dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz)
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
use phys_module

implicit none

logical :: xpoint2
real*8  :: prof0, prof1, dprof0_dpsi, dprof0_dpsi2, dprof0_dpsi3, density_profile, psi_barrier
real*8  :: dn_dpsi, dn_dz, dn_dpsi2, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz
real*8  :: Z, Z_xpoint, psi, psi_axis, psi_bnd,  psi_n, sig_n, sigz, dprof1_dpsi, dprof1_dpsi2, dprof1_dpsi3
real*8  :: atn, datn, d2atn, d3atn, atn_z, datn_z, d2atn_z, factor

sig_n       = rho_coef(4)
psi_barrier = rho_coef(5)

psi_n = (psi - psi_axis) / (psi_bnd - psi_axis)

factor = 1.d0
!if (xpoint2) then
!  if ((Z .lt. Z_xpoint) .and. (psi_n .lt. 1.d0) ) then
!    psi_n = 2.d0 - psi_n
!    factor = -1.d0
!  endif
!endif

prof0        = (rho_0-rho_1)*(1.d0 + rho_coef(1) * psi_n + rho_coef(2) * psi_n**2 + rho_coef(3) * psi_n**3)
dprof0_dpsi  = (rho_0-rho_1)*(rho_coef(1) + 2.d0 * rho_coef(2) * psi_n + 3.d0 * rho_coef(3) * psi_n**2) / (psi_bnd - psi_axis)
dprof0_dpsi2 = (rho_0-rho_1)*(2.d0 * rho_coef(2) + 6.d0 * rho_coef(3) * psi_n)                          / (psi_bnd - psi_axis)**2
dprof0_dpsi3 = (rho_0-rho_1)*(6.d0 * rho_coef(3))                                                       / (psi_bnd - psi_axis)**3

atn   = (0.5d0 - 0.5d0*tanh((psi_n - psi_barrier)/sig_n))

datn  = - 1.d0/cosh((psi_n - psi_barrier)/sig_n)**2 / (2.d0 * sig_n) / (psi_bnd - psi_axis)

d2atn =   1.d0/cosh((psi_n - psi_barrier)/sig_n)**2 / (sig_n**2)  &
      * tanh((psi_n - psi_barrier)/sig_n) / (psi_bnd - psi_axis)**2

!d3atn = - 1.d0/cosh((psi_n - psi_barrier)/sig_n)**4 / (sig_n**3)  &
!      * (-2.d0 + cosh(2.d0*(psi_n-psi_barrier)/sig_n) ) / (psi_bnd - psi_axis)**3
d3atn = - 1.d0 / sig_n**3 / (psi_bnd - psi_axis)**3 &
      * ( -2.d0/cosh((psi_n-psi_barrier)/sig_n)**4 + 1.d0/cosh(2.d0*(psi_n-psi_barrier)/sig_n)**3 )

prof1        = prof0        * atn
dprof1_dpsi  = dprof0_dpsi  * atn +         prof0       * datn
dprof1_dpsi2 = dprof0_dpsi2 * atn + 2.d0 * dprof0_dpsi  * datn + prof0              * d2atn
dprof1_dpsi3 = dprof0_dpsi3 * atn + 3.d0 * dprof0_dpsi2 * datn + 3.d0 * dprof0_dpsi * d2atn + prof0 * d3atn

dn_dpsi     = dprof1_dpsi   * factor
dn_dpsi2    = dprof1_dpsi2
dn_dpsi3    = dprof1_dpsi3  * factor

dn_dz       = 0.d0
dn_dz2      = 0.d0
dn_dpsi_dz  = 0.d0
dn_dpsi2_dz = 0.d0
dn_dpsi_dz2 = 0.d0

density_profile = prof1

if (xpoint2) then
  sigz    = 0.1d0

  atn_z           =  (0.5d0 - 0.5d0*tanh((Z_xpoint-Z)/sigz))
  datn_z          =  0.5d0/cosh((Z_xpoint-Z)/sigz)**2 / sigz
  d2atn_z         = 1.d0/cosh((Z_xpoint-Z)/sigz)**2 /  sigz**2 * tanh((Z_xpoint-Z)/sigz)

  density_profile =   prof1        * atn_z
  dn_dpsi         =   dprof1_dpsi  * atn_z
  dn_dpsi2        =   dprof1_dpsi2 * atn_z
  dn_dpsi3        =   dprof1_dpsi3 * atn_z
  dn_dz           = + prof1        * datn_z
  dn_dz2          = + prof1        * d2atn_z
  dn_dpsi_dz      =   dprof1_dpsi  * datn_z
  dn_dpsi2_dz     =   dprof1_dpsi2 * datn_z
  dn_dpsi_dz2     =   dprof1_dpsi  * d2atn_z

endif

density_profile = density_profile + rho_1

return
end
