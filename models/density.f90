subroutine density(xpoint2,xcase2,Z,Z_xpoint,psi,psi_axis,psi_bnd,density_profile,&
                   dn_dpsi,  dn_dz, &                                             ! 1st order derivatives
                   dn_dpsi2, dn_dz2,      dn_dpsi_dz, &                           ! 2nd order derivatives
                   dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz,  dn_dz3, &                 ! 2rd order derivatives
                   dn_dpsi4, dn_dpsi_dz3, dn_dpsi2_dz2, dn_dpsi3_dz,  dn_dz4, &   ! 4th order derivatives
                   dn_dpsi5, dn_dpsi_dz4, dn_dpsi2_dz3, dn_dpsi3_dz2, dn_dpsi4_dz)! 5th order derivatives (z5 not needed)
!-----------------------------------------------------------------------
! Determines the density value and its derivatives at the given
! position (Z, psi) from the analytical or numerical input profile.
!-----------------------------------------------------------------------
use phys_module
use mod_input_profiles

implicit none

! --- Routine parameters
logical, intent(in)  :: xpoint2
integer, intent(in)  :: xcase2
real*8,  intent(in)  :: Z, Z_xpoint(2), psi, psi_axis, psi_bnd
real*8,  intent(out) :: density_profile
real*8,  intent(out) :: dn_dpsi,  dn_dz
real*8,  intent(out) :: dn_dpsi2, dn_dz2,      dn_dpsi_dz
real*8,  intent(out) :: dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz,  dn_dz3
real*8,  intent(out) :: dn_dpsi4, dn_dpsi_dz3, dn_dpsi2_dz2, dn_dpsi3_dz,  dn_dz4
real*8,  intent(out) :: dn_dpsi5, dn_dpsi_dz4, dn_dpsi2_dz3, dn_dpsi3_dz2, dn_dpsi4_dz

! --- Internal variables.
real*8  :: prof1, dprof1_dpsi, dprof1_dpsi2, dprof1_dpsi3, dprof1_dpsi4, dprof1_dpsi5
real*8  :: atn, datn, d2atn, d3atn, d4atn
real*8  :: delta_psi, psi_n
! for interpolating numerical profiles
integer :: left, right, mid
real*8  :: aux1, aux2

delta_psi = psi_bnd - psi_axis
psi_n     = (psi - psi_axis) / delta_psi
psi_n     = max( min(psi_n, 2.), 0. )

! --- Profile as a function of Psi_N.
if ( .not. num_rho ) then ! use analytical representation

  call input_profiles_psi_component(psi,psi_axis,psi_bnd, rho_0, rho_1, rho_coef, &
                                    prof1,dprof1_dpsi,dprof1_dpsi2,dprof1_dpsi3,dprof1_dpsi4,dprof1_dpsi5)

else ! use numerical representation.
  
  ! --- Interpolate profile and derivatives to position psi_n by bisections.
  left  = 1
  right = num_rho_len
  do
    if ( right == left + 1 ) exit
    mid = (left + right) / 2
    if ( num_rho_x(mid) >= psi_n ) then
      right = mid
    else
      left = mid
    end if
  end do
  aux1 = (psi_n - num_rho_x(left)) / (num_rho_x(right) - num_rho_x(left))
  aux2 = (1. - aux1)
  prof1        =   num_rho_y0(left) * aux2 + num_rho_y0(right) * aux1
  dprof1_dpsi  = ( num_rho_y1(left) * aux2 + num_rho_y1(right) * aux1 ) / delta_psi
  dprof1_dpsi2 = ( num_rho_y2(left) * aux2 + num_rho_y2(right) * aux1 ) / delta_psi**2
  dprof1_dpsi3 = ( num_rho_y3(left) * aux2 + num_rho_y3(right) * aux1 ) / delta_psi**3
  dprof1_dpsi4 = ( num_rho_y4(left) * aux2 + num_rho_y4(right) * aux1 ) / delta_psi**4
  dprof1_dpsi5 = ( num_rho_y5(left) * aux2 + num_rho_y5(right) * aux1 ) / delta_psi**5
  
end if

! compute z-tanh (will be 1 with zero deivatives if xpoint=.f.)
call input_profiles_Z_component(xpoint2,xcase2,Z,Z_xpoint, atn,datn,d2atn,d3atn,d4atn)  

density_profile = prof1        *   atn
dn_dpsi         = dprof1_dpsi  *   atn
dn_dpsi2        = dprof1_dpsi2 *   atn
dn_dpsi3        = dprof1_dpsi3 *   atn
dn_dpsi4        = dprof1_dpsi4 *   atn
dn_dpsi5        = dprof1_dpsi5 *   atn
dn_dz           = prof1        *  datn
dn_dz2          = prof1        * d2atn
dn_dz3          = prof1        * d3atn
dn_dz4          = prof1        * d4atn
dn_dpsi_dz      = dprof1_dpsi  *  datn
dn_dpsi_dz2     = dprof1_dpsi  * d2atn
dn_dpsi_dz3     = dprof1_dpsi  * d3atn
dn_dpsi_dz4     = dprof1_dpsi  * d4atn
dn_dpsi2_dz     = dprof1_dpsi2 *  datn
dn_dpsi2_dz2    = dprof1_dpsi2 * d2atn
dn_dpsi2_dz3    = dprof1_dpsi2 * d3atn
dn_dpsi3_dz     = dprof1_dpsi3 *  datn
dn_dpsi3_dz2    = dprof1_dpsi3 * d2atn
dn_dpsi4_dz     = dprof1_dpsi4 *  datn

density_profile = density_profile + rho_1

return
end subroutine density
