subroutine temperature(xpoint2,Z,Z_xpoint,psi,psi_axis,psi_bnd,temperature_profile, &
                       dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2, dT_dpsi2_dz)
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
use phys_module

implicit none

logical :: xpoint2
real*8  :: prof0, prof1, dprof0_dpsi, dprof0_dpsi2, dprof0_dpsi3, temperature_profile, psi_barrier
real*8  :: dT_dpsi, dT_dz, dT_dpsi2, dT_dz2, dT_dpsi_dz, dT_dpsi3, dT_dpsi_dz2, dT_dpsi2_dz
real*8  :: Z, Z_xpoint, psi, psi_axis, psi_bnd,  psi_n, sig_T, sigz, dprof1_dpsi, dprof1_dpsi2, dprof1_dpsi3
real*8  :: atn, datn, d2atn, d3atn, atn_z, datn_z, d2atn_z, factor

sig_T       = T_coef(4)
psi_barrier = T_coef(5)

psi_n = (psi - psi_axis) / (psi_bnd - psi_axis)

factor = 1.d0
!if (xpoint2) then
!  if ((Z .lt. Z_xpoint) .and. (psi_n .lt. 1.d0) ) then
!    psi_n  = 2.d0 - psi_n
!    factor = -1.d0
!  endif
!endif

prof0        = (T_0 - T_1) * (1.d0 + T_coef(1) * psi_n + T_coef(2) * psi_n**2 + T_coef(3) * psi_n**3)
dprof0_dpsi  = (T_0 - T_1) * (T_coef(1) + 2.d0 * T_coef(2) * psi_n + 3.d0 * T_coef(3) * psi_n**2) / (psi_bnd - psi_axis)
dprof0_dpsi2 = (T_0 - T_1) * (2.d0 * T_coef(2) + 6.d0 * T_coef(3) * psi_n)                        / (psi_bnd - psi_axis)**2
dprof0_dpsi3 = (T_0 - T_1) * (6.d0 * T_coef(3))                                                   / (psi_bnd - psi_axis)**3

atn   = (0.5d0 - 0.5d0*tanh((psi_n - psi_barrier)/sig_T))

datn  = - 1.d0/cosh((psi_n - psi_barrier)/sig_T)**2 / (2.d0 * sig_T) / (psi_bnd - psi_axis)

d2atn =   1.d0/cosh((psi_n - psi_barrier)/sig_T)**2 / (sig_T**2)  &
      * tanh((psi_n - psi_barrier)/sig_T) / (psi_bnd - psi_axis)**2

!d3atn = - 1.d0/cosh((psi_n - psi_barrier)/sig_T)**4 / (sig_T**3)  &
!      * (-2.d0 + cosh(2.d0*(psi_n-psi_barrier)/sig_T) ) / (psi_bnd - psi_axis)**3
d3atn = - 1.d0 / sig_T**3 / (psi_bnd - psi_axis)**3 &
      * ( -2.d0/cosh((psi_n-psi_barrier)/sig_T)**4 + 1.d0/cosh(2.d0*(psi_n-psi_barrier)/sig_T)**3 )

prof1        = prof0        * atn
dprof1_dpsi  = dprof0_dpsi  * atn +         prof0       * datn
dprof1_dpsi2 = dprof0_dpsi2 * atn + 2.d0 * dprof0_dpsi  * datn + prof0       * d2atn
dprof1_dpsi3 = dprof0_dpsi3 * atn + 3.d0 * dprof0_dpsi2 * datn + 3.d0 * dprof0_dpsi * d2atn + prof0 * d3atn

dT_dpsi     = dprof1_dpsi     * factor
dT_dpsi2    = dprof1_dpsi2
dT_dpsi3    = dprof1_dpsi3    * factor

dT_dz       = 0.d0
dT_dz2      = 0.d0
dT_dpsi_dz  = 0.d0
dT_dpsi2_dz = 0.d0
dT_dpsi_dz2 = 0.d0

temperature_profile = prof1

if (xpoint2) then
  sigz    = 0.05d0

  atn_z           =  (0.5d0 - 0.5d0*tanh((Z_xpoint-Z)/sigz))
  datn_z          =  0.5d0/cosh((Z_xpoint-Z)/sigz)**2 / sigz
  d2atn_z         = 1.d0/cosh((Z_xpoint-Z)/sigz)**2 /  sigz**2 * tanh((Z_xpoint-Z)/sigz)

  temperature_profile =   prof1        * atn_z
  dT_dpsi         =   dprof1_dpsi  * atn_z
  dT_dpsi2        =   dprof1_dpsi2 * atn_z
  dT_dpsi3        =   dprof1_dpsi3 * atn_z
  dT_dz           = + prof1        * datn_z
  dT_dz2          = + prof1        * d2atn_z
  dT_dpsi_dz      =   dprof1_dpsi  * datn_z
  dT_dpsi2_dz     =   dprof1_dpsi2 * datn_z
  dT_dpsi_dz2     =   dprof1_dpsi  * d2atn_z

endif

temperature_profile = temperature_profile + T_1

return
end

