subroutine FFprime(xpoint2,Z,Z_xpoint,psi,psi_axis,psi_bnd,FFprime_profile,dFF_dpsi,dFF_dz, &
                       dFF_dpsi2,dFF_dz2,dFF_dpsi_dz)
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
use phys_module

implicit none

logical :: xpoint2
real*8  :: prof0, prof1, dprof0_dpsi, dprof0_dpsi2, FFprime_profile, psi_barrier
real*8  :: dFF_dpsi, dFF_dz, dFF_dpsi2, dFF_dz2, dFF_dpsi_dz
real*8  :: Z, Z_xpoint, psi, psi_axis, psi_bnd,  psi_n, sig_F, sigz, dprof1_dpsi, dprof1_dpsi2
real*8  :: atn, datn, d2atn, d3atn, atn_z, datn_z, d2atn_z, factor, d_pert, d2_pert, d3_pert

sig_F       = FF_coef(4)
psi_barrier = FF_coef(5)

psi_n = (psi - psi_axis) / (psi_bnd - psi_axis)

factor = 1.d0
!if (xpoint2) then
!  if ((Z .lt. Z_xpoint) .and. (psi_n .lt. 1.d0) ) then
!    psi_n = 2.d0 - psi_n
!    factor = -1.d0
!  endif
!endif

d_pert  = + FF_coef(6)/cosh((psi_n - FF_coef(7))/FF_coef(8))**2 / (2.d0 * FF_coef(8)) / (psi_bnd - psi_axis)
d2_pert = - FF_coef(6)/cosh((psi_n - FF_coef(7))/FF_coef(8))**2 / (FF_coef(8)**2)  &
        * tanh((psi_n - FF_coef(7))/FF_coef(8)) / (psi_bnd - psi_axis)**2
d3_pert = + FF_coef(6)/cosh((psi_n - FF_coef(7))/FF_coef(8))**4 / (FF_coef(8)**3)  &
        * (-2.d0 + cosh(2.d0*(psi_n-FF_coef(7))/FF_coef(8)) ) / (psi_bnd - psi_axis)**3

prof0        = (FF_0 - FF_1) * ( 1.d0 + FF_coef(1) * psi_n + FF_coef(2) * psi_n**2 + FF_coef(3) * psi_n**3)
dprof0_dpsi  = (FF_0 - FF_1) * ( FF_coef(1) + 2.d0 * FF_coef(2) * psi_n + 3.d0 * FF_coef(3) * psi_n**2)    / (psi_bnd - psi_axis)
dprof0_dpsi2 = (FF_0 - FF_1) * (2.d0 * FF_coef(2) + 6.d0 * FF_coef(3) * psi_n) / (psi_bnd - psi_axis)**2

prof0        = prof0        + d_pert
dprof0_dpsi  = dprof0_dpsi  + d2_pert
dprof0_dpsi2 = dprof0_dpsi2 + d3_pert

atn   = (0.5d0 - 0.5d0*tanh((psi_n - psi_barrier)/sig_F))

datn  = - 1.d0/cosh((psi_n - psi_barrier)/sig_F)**2 / (2.d0 * sig_F) / (psi_bnd - psi_axis)

d2atn =   1.d0/cosh((psi_n - psi_barrier)/sig_F)**2 / (sig_F**2)  &
      * tanh((psi_n - psi_barrier)/sig_F) / (psi_bnd - psi_axis)**2

!d3atn = - 1.d0/cosh((psi_n - psi_barrier)/sig_F)**4 / (sig_F**3)  &
!      * (-2.d0 + cosh(2.d0*(psi_n-psi_barrier)/sig_F) ) / (psi_bnd - psi_axis)**3
d3atn = - 1.d0 / sig_F**3 / (psi_bnd - psi_axis)**3 &
      * ( -2.d0/cosh((psi_n-psi_barrier)/sig_F)**4 + 1.d0/cosh(2.d0*(psi_n-psi_barrier)/sig_F)**3 ) 

prof1        = prof0        * atn
dprof1_dpsi  = dprof0_dpsi  * atn +         prof0       * datn
dprof1_dpsi2 = dprof0_dpsi2 * atn + 2.d0 * dprof0_dpsi  * datn + prof0       * d2atn

dFF_dpsi     = dprof1_dpsi   * factor
dFF_dpsi2    = dprof1_dpsi2

dFF_dz       = 0.d0
dFF_dz2      = 0.d0
dFF_dpsi_dz  = 0.d0

FFprime_profile = prof1

if (xpoint2) then
  sigz    = 0.1d0

  atn_z           =  (0.5d0 - 0.5d0*tanh((Z_xpoint-Z)/sigz))
  datn_z          =  0.5d0/cosh((Z_xpoint-Z)/sigz)**2 / sigz
  d2atn_z         = 1.d0/cosh((Z_xpoint-Z)/sigz)**2 /  sigz**2 * tanh((Z_xpoint-Z)/sigz)

  FFprime_profile  =   prof1        * atn_z
  dFF_dpsi         =   dprof1_dpsi  * atn_z
  dFF_dpsi2        =   dprof1_dpsi2 * atn_z
  dFF_dz           = + prof1        * datn_z
  dFF_dz2          = + prof1        * d2atn_z
  dFF_dpsi_dz      =   dprof1_dpsi  * datn_z

endif

FFprime_profile = FFprime_profile + FF_1

return
end
