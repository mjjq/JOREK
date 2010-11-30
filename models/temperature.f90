subroutine temperature(xpoint2,Z,Z_xpoint,psi,psi_axis,psi_bnd,temperature_profile, &
                       dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2, dT_dpsi2_dz)
!-----------------------------------------------------------------------
! Determines the temperature value and its derivatives at the given
! position (Z, psi) from the analytical or numerical input profile.
!-----------------------------------------------------------------------
use phys_module

implicit none

! --- Routine parameters
logical, intent(in)  :: xpoint2
real*8,  intent(in)  :: Z, Z_xpoint, psi, psi_axis, psi_bnd
real*8,  intent(out) :: temperature_profile, dT_dpsi, dT_dz, dT_dpsi2, dT_dz2, &
                        dT_dpsi_dz, dT_dpsi3, dT_dpsi_dz2, dT_dpsi2_dz

! --- Internal variables.
real*8  :: prof0, prof1, dprof0_dpsi, dprof0_dpsi2, dprof0_dpsi3, psi_barrier
real*8  :: psi_n, delta_psi, sig_T, sigz, dprof1_dpsi, dprof1_dpsi2, dprof1_dpsi3
real*8  :: atn, datn, d2atn, d3atn, atn_z, datn_z, d2atn_z, factor
real*8  :: cosh1, cosh2, cosh3, tanh1, tanh2
! for interpolating numerical profiles
integer :: left, right, mid
real*8  :: aux1, aux2

delta_psi = psi_bnd - psi_axis
psi_n     = (psi - psi_axis) / delta_psi

!factor = 1.d0
!if (xpoint2) then
!  if ((Z .lt. Z_xpoint) .and. (psi_n .lt. 1.d0) ) then
!    psi_n = 2.d0 - psi_n
!    factor = -1.d0
!  endif
!endif

! --- Profile as a function of Psi_N.
if ( .not. num_T ) then ! use analytical representation
  
  prof0        = (T_0-T_1)*(1.d0 + T_coef(1) * psi_n + T_coef(2) * psi_n**2 + T_coef(3) * psi_n**3)
  dprof0_dpsi  = (T_0-T_1)*(T_coef(1) + 2.d0 * T_coef(2) * psi_n + 3.d0 * T_coef(3) * psi_n**2) / delta_psi
  dprof0_dpsi2 = (T_0-T_1)*(2.d0 * T_coef(2) + 6.d0 * T_coef(3) * psi_n)                        / delta_psi**2
  dprof0_dpsi3 = (T_0-T_1)*(6.d0 * T_coef(3))                                                   / delta_psi**3
  
  sig_T       = T_coef(4)
  psi_barrier = T_coef(5)
  
  tanh1 = tanh((psi_n - psi_barrier)/sig_T)
  cosh1 = cosh((psi_n - psi_barrier)/sig_T)
  cosh2 = cosh(2.d0*(psi_n - psi_barrier)/sig_T)
  
  atn   = (0.5d0 - 0.5d0*tanh1)
  datn  = - 1.d0/cosh1**2 / (2.d0 * sig_T) / delta_psi
  d2atn =   1.d0/cosh1**2 / sig_T**2 * tanh1 / delta_psi**2
  d3atn = - 1.d0/cosh1**4 / sig_T**3 * (-2.d0 + cosh2) / delta_psi**3

  prof1        = prof0        * atn
  dprof1_dpsi  = dprof0_dpsi  * atn +         prof0       * datn
  dprof1_dpsi2 = dprof0_dpsi2 * atn + 2.d0 * dprof0_dpsi  * datn + prof0              * d2atn
  dprof1_dpsi3 = dprof0_dpsi3 * atn + 3.d0 * dprof0_dpsi2 * datn + 3.d0 * dprof0_dpsi * d2atn + prof0 * d3atn
  
else ! use numerical representation.
  
  ! --- Interpolate profile and derivatives to position psi_n by bisections.
  left  = 1
  right = num_T_len
  do
    if ( right == left + 1 ) exit
    mid = (left + right) / 2
    if ( num_T_x(mid) >= psi_n ) then
      right = mid
    else
      left = mid
    end if
  end do
  aux1 = (psi_n - num_T_x(left)) / (num_T_x(right) - num_T_x(left))
  aux2 = (1. - aux1)
  prof1        = num_T_y0(left)   * aux2 + num_T_y0(right) * aux1
  dprof1_dpsi  = ( num_T_y1(left) * aux2 + num_T_y1(right) * aux1 ) / delta_psi
  dprof1_dpsi2 = ( num_T_y2(left) * aux2 + num_T_y2(right) * aux1 ) / delta_psi**2
  dprof1_dpsi3 = ( num_T_y3(left) * aux2 + num_T_y3(right) * aux1 ) / delta_psi**3
  
end if

! --- Additional explicit dependence of the profile on Z to ensure that the profile is
!     approximately zero in the private flux region below the x-point.
if (xpoint2) then
  
  sigz            = 0.05d0

  tanh2 = tanh((Z_xpoint-Z)/sigz)
  cosh3 = cosh((Z_xpoint-Z)/sigz)
  
  atn_z            = (0.5d0 - 0.5d0*tanh2)
  datn_z           = 0.5d0/cosh3**2 / sigz
  d2atn_z          = 1.0d0/cosh3**2 / sigz**2 * tanh2

  temperature_profile = prof1      * atn_z
  dT_dpsi         =   dprof1_dpsi  * atn_z
  dT_dpsi2        =   dprof1_dpsi2 * atn_z
  dT_dpsi3        =   dprof1_dpsi3 * atn_z
  dT_dz           = + prof1        * datn_z
  dT_dz2          = + prof1        * d2atn_z
  dT_dpsi_dz      =   dprof1_dpsi  * datn_z
  dT_dpsi2_dz     =   dprof1_dpsi2 * datn_z
  dT_dpsi_dz2     =   dprof1_dpsi  * d2atn_z
  
else
  
  temperature_profile = prof1
  dT_dpsi     = dprof1_dpsi!   * factor
  dT_dpsi2    = dprof1_dpsi2
  dT_dpsi3    = dprof1_dpsi3!  * factor
  dT_dz       = 0.d0
  dT_dz2      = 0.d0
  dT_dpsi_dz  = 0.d0
  dT_dpsi2_dz = 0.d0
  dT_dpsi_dz2 = 0.d0

end if

temperature_profile = temperature_profile + T_1

return
end subroutine temperature
