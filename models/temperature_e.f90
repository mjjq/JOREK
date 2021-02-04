subroutine temperature_e(xpoint2,xcase2,Z,Z_xpoint,psi,psi_axis,psi_bnd,temperature_e_profile, &
                         dTe_dpsi,  dTe_dz, &                                                ! 1st order derivatives
                         dTe_dpsi2, dTe_dz2,      dTe_dpsi_dz, &                             ! 2nd order derivatives
                         dTe_dpsi3, dTe_dpsi_dz2, dTe_dpsi2_dz,  dTe_dz3, &                  ! 2rd order derivatives
                         dTe_dpsi4, dTe_dpsi_dz3, dTe_dpsi2_dz2, dTe_dpsi3_dz,  dTe_dz4, &   ! 4th order derivatives
                         dTe_dpsi5, dTe_dpsi_dz4, dTe_dpsi2_dz3, dTe_dpsi3_dz2, dTe_dpsi4_dz)! 5th order derivatives (z5 not needed)
!-----------------------------------------------------------------------
! Determines the temperature value and its derivatives at the given
! position (Z, psi) from the analytical or numerical input profile.
!-----------------------------------------------------------------------
use phys_module
use mod_input_profiles

implicit none

! --- Routine parameters
logical, intent(in)  :: xpoint2
integer, intent(in)  :: xcase2
real*8,  intent(in)  :: Z, Z_xpoint(2), psi, psi_axis, psi_bnd
real*8,  intent(out) :: temperature_e_profile
real*8,  intent(out) :: dTe_dpsi,  dTe_dz
real*8,  intent(out) :: dTe_dpsi2, dTe_dz2,      dTe_dpsi_dz
real*8,  intent(out) :: dTe_dpsi3, dTe_dpsi_dz2, dTe_dpsi2_dz,  dTe_dz3
real*8,  intent(out) :: dTe_dpsi4, dTe_dpsi_dz3, dTe_dpsi2_dz2, dTe_dpsi3_dz,  dTe_dz4
real*8,  intent(out) :: dTe_dpsi5, dTe_dpsi_dz4, dTe_dpsi2_dz3, dTe_dpsi3_dz2, dTe_dpsi4_dz

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
if ( .not. num_Te ) then ! use analytical representation

  call input_profiles_psi_component(psi,psi_axis,psi_bnd, Te_0, Te_1, Te_coef, &
                                    prof1,dprof1_dpsi,dprof1_dpsi2,dprof1_dpsi3,dprof1_dpsi4,dprof1_dpsi5)

else ! use numerical representation.
  
  ! --- Interpolate profile and derivatives to position psi_n by bisections.
  left  = 1
  right = num_Te_len
  do
    if ( right == left + 1 ) exit
    mid = (left + right) / 2
    if ( num_Te_x(mid) >= psi_n ) then
      right = mid
    else
      left = mid
    end if
  end do
  aux1 = (psi_n - num_Te_x(left)) / (num_Te_x(right) - num_Te_x(left))
  aux2 = (1. - aux1)
  prof1        =   num_Te_y0(left) * aux2 + num_Te_y0(right) * aux1
  dprof1_dpsi  = ( num_Te_y1(left) * aux2 + num_Te_y1(right) * aux1 ) / delta_psi
  dprof1_dpsi2 = ( num_Te_y2(left) * aux2 + num_Te_y2(right) * aux1 ) / delta_psi**2
  dprof1_dpsi3 = ( num_Te_y3(left) * aux2 + num_Te_y3(right) * aux1 ) / delta_psi**3
  dprof1_dpsi4 = ( num_Te_y4(left) * aux2 + num_Te_y4(right) * aux1 ) / delta_psi**4
  dprof1_dpsi5 = ( num_Te_y5(left) * aux2 + num_Te_y5(right) * aux1 ) / delta_psi**5
  
end if

! compute z-tanh (will be 1 with zero deivatives if xpoint=.f.)
call input_profiles_Z_component(xpoint2,xcase2,Z,Z_xpoint, atn,datn,d2atn,d3atn,d4atn)  

temperature_e_profile = prof1        *   atn
dTe_dpsi              = dprof1_dpsi  *   atn
dTe_dpsi2             = dprof1_dpsi2 *   atn
dTe_dpsi3             = dprof1_dpsi3 *   atn
dTe_dpsi4             = dprof1_dpsi4 *   atn
dTe_dpsi5             = dprof1_dpsi5 *   atn
dTe_dz                = prof1        *  datn
dTe_dz2               = prof1        * d2atn
dTe_dz3               = prof1        * d3atn
dTe_dz4               = prof1        * d4atn
dTe_dpsi_dz           = dprof1_dpsi  *  datn
dTe_dpsi_dz2          = dprof1_dpsi  * d2atn
dTe_dpsi_dz3          = dprof1_dpsi  * d3atn
dTe_dpsi_dz4          = dprof1_dpsi  * d4atn
dTe_dpsi2_dz          = dprof1_dpsi2 *  datn
dTe_dpsi2_dz2         = dprof1_dpsi2 * d2atn
dTe_dpsi2_dz3         = dprof1_dpsi2 * d3atn
dTe_dpsi3_dz          = dprof1_dpsi3 *  datn
dTe_dpsi3_dz2         = dprof1_dpsi3 * d2atn
dTe_dpsi4_dz          = dprof1_dpsi4 *  datn



if (freeboundary_equil .and. num_Te) then                        !if the temperature profile is given in a file and there is freeboundary equilibrium
                                                                !the full profile is multiplied by a facto in order to iterate to a given current
  temperature_e_profile = temperature_e_profile * current_FB_fact
  dTe_dpsi              = dprof1_dpsi  * current_FB_fact
  dTe_dpsi2             = dprof1_dpsi2 * current_FB_fact
  dTe_dpsi3             = dprof1_dpsi3 * current_FB_fact
  dTe_dpsi4             = dprof1_dpsi4 * current_FB_fact
  dTe_dpsi5             = dprof1_dpsi5 * current_FB_fact
  dTe_dz                = prof1        * current_FB_fact
  dTe_dz2               = prof1        * current_FB_fact
  dTe_dz3               = prof1        * current_FB_fact
  dTe_dz4               = prof1        * current_FB_fact
  dTe_dpsi_dz           = dprof1_dpsi  * current_FB_fact
  dTe_dpsi_dz2          = dprof1_dpsi  * current_FB_fact
  dTe_dpsi_dz3          = dprof1_dpsi  * current_FB_fact
  dTe_dpsi_dz4          = dprof1_dpsi  * current_FB_fact
  dTe_dpsi2_dz          = dprof1_dpsi2 * current_FB_fact
  dTe_dpsi2_dz2         = dprof1_dpsi2 * current_FB_fact
  dTe_dpsi2_dz3         = dprof1_dpsi2 * current_FB_fact
  dTe_dpsi3_dz          = dprof1_dpsi3 * current_FB_fact
  dTe_dpsi3_dz2         = dprof1_dpsi3 * current_FB_fact
  dTe_dpsi4_dz          = dprof1_dpsi4 * current_FB_fact

end if

temperature_e_profile = temperature_e_profile + Te_1

return
end subroutine temperature_e
