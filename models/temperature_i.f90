subroutine temperature_i(xpoint2,xcase2,Z,Z_xpoint,psi,psi_axis,psi_bnd,temperature_i_profile, &
                         dTi_dpsi,  dTi_dz, &                                                ! 1st order derivatives
                         dTi_dpsi2, dTi_dz2,      dTi_dpsi_dz, &                             ! 2nd order derivatives
                         dTi_dpsi3, dTi_dpsi_dz2, dTi_dpsi2_dz,  dTi_dz3, &                  ! 2rd order derivatives
                         dTi_dpsi4, dTi_dpsi_dz3, dTi_dpsi2_dz2, dTi_dpsi3_dz,  dTi_dz4, &   ! 4th order derivatives
                         dTi_dpsi5, dTi_dpsi_dz4, dTi_dpsi2_dz3, dTi_dpsi3_dz2, dTi_dpsi4_dz)! 5th order derivatives (z5 not needed)
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
real*8,  intent(out) :: temperature_i_profile
real*8,  intent(out) :: dTi_dpsi,  dTi_dz
real*8,  intent(out) :: dTi_dpsi2, dTi_dz2,      dTi_dpsi_dz
real*8,  intent(out) :: dTi_dpsi3, dTi_dpsi_dz2, dTi_dpsi2_dz,  dTi_dz3
real*8,  intent(out) :: dTi_dpsi4, dTi_dpsi_dz3, dTi_dpsi2_dz2, dTi_dpsi3_dz,  dTi_dz4
real*8,  intent(out) :: dTi_dpsi5, dTi_dpsi_dz4, dTi_dpsi2_dz3, dTi_dpsi3_dz2, dTi_dpsi4_dz

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
if ( .not. num_Ti ) then ! use analytical representation

  call input_profiles_psi_component(psi,psi_axis,psi_bnd, Ti_0, Ti_1, Ti_coef, &
                                    prof1,dprof1_dpsi,dprof1_dpsi2,dprof1_dpsi3,dprof1_dpsi4,dprof1_dpsi5)

else ! use numerical representation.
  
  ! --- Interpolate profile and derivatives to position psi_n by bisections.
  left  = 1
  right = num_Ti_len
  do
    if ( right == left + 1 ) exit
    mid = (left + right) / 2
    if ( num_Ti_x(mid) >= psi_n ) then
      right = mid
    else
      left = mid
    end if
  end do
  aux1 = (psi_n - num_Ti_x(left)) / (num_Ti_x(right) - num_Ti_x(left))
  aux2 = (1. - aux1)
  prof1        =   num_Ti_y0(left) * aux2 + num_Ti_y0(right) * aux1
  dprof1_dpsi  = ( num_Ti_y1(left) * aux2 + num_Ti_y1(right) * aux1 ) / delta_psi
  dprof1_dpsi2 = ( num_Ti_y2(left) * aux2 + num_Ti_y2(right) * aux1 ) / delta_psi**2
  dprof1_dpsi3 = ( num_Ti_y3(left) * aux2 + num_Ti_y3(right) * aux1 ) / delta_psi**3
  dprof1_dpsi4 = ( num_Ti_y4(left) * aux2 + num_Ti_y4(right) * aux1 ) / delta_psi**4
  dprof1_dpsi5 = ( num_Ti_y5(left) * aux2 + num_Ti_y5(right) * aux1 ) / delta_psi**5
  
end if

! compute z-tanh (will be 1 with zero deivatives if xpoint=.f.)
call input_profiles_Z_component(xpoint2,xcase2,Z,Z_xpoint, atn,datn,d2atn,d3atn,d4atn)  

temperature_i_profile = prof1        *   atn
dTi_dpsi              = dprof1_dpsi  *   atn
dTi_dpsi2             = dprof1_dpsi2 *   atn
dTi_dpsi3             = dprof1_dpsi3 *   atn
dTi_dpsi4             = dprof1_dpsi4 *   atn
dTi_dpsi5             = dprof1_dpsi5 *   atn
dTi_dz                = prof1        *  datn
dTi_dz2               = prof1        * d2atn
dTi_dz3               = prof1        * d3atn
dTi_dz4               = prof1        * d4atn
dTi_dpsi_dz           = dprof1_dpsi  *  datn
dTi_dpsi_dz2          = dprof1_dpsi  * d2atn
dTi_dpsi_dz3          = dprof1_dpsi  * d3atn
dTi_dpsi_dz4          = dprof1_dpsi  * d4atn
dTi_dpsi2_dz          = dprof1_dpsi2 *  datn
dTi_dpsi2_dz2         = dprof1_dpsi2 * d2atn
dTi_dpsi2_dz3         = dprof1_dpsi2 * d3atn
dTi_dpsi3_dz          = dprof1_dpsi3 *  datn
dTi_dpsi3_dz2         = dprof1_dpsi3 * d2atn
dTi_dpsi4_dz          = dprof1_dpsi4 *  datn


if (freeboundary_equil .and. num_Ti) then                        !if the temperature profile is given in a file and there is freeboundary equilibrium
                                                                !the full profile is multiplied by a facto in order to iterate to a given current
  temperature_i_profile = temperature_i_profile * current_FB_fact
  dTi_dpsi              = dprof1_dpsi  * current_FB_fact
  dTi_dpsi2             = dprof1_dpsi2 * current_FB_fact
  dTi_dpsi3             = dprof1_dpsi3 * current_FB_fact
  dTi_dpsi4             = dprof1_dpsi4 * current_FB_fact
  dTi_dpsi5             = dprof1_dpsi5 * current_FB_fact
  dTi_dz                = prof1        * current_FB_fact
  dTi_dz2               = prof1        * current_FB_fact
  dTi_dz3               = prof1        * current_FB_fact
  dTi_dz4               = prof1        * current_FB_fact
  dTi_dpsi_dz           = dprof1_dpsi  * current_FB_fact
  dTi_dpsi_dz2          = dprof1_dpsi  * current_FB_fact
  dTi_dpsi_dz3          = dprof1_dpsi  * current_FB_fact
  dTi_dpsi_dz4          = dprof1_dpsi  * current_FB_fact
  dTi_dpsi2_dz          = dprof1_dpsi2 * current_FB_fact
  dTi_dpsi2_dz2         = dprof1_dpsi2 * current_FB_fact
  dTi_dpsi2_dz3         = dprof1_dpsi2 * current_FB_fact
  dTi_dpsi3_dz          = dprof1_dpsi3 * current_FB_fact
  dTi_dpsi3_dz2         = dprof1_dpsi3 * current_FB_fact
  dTi_dpsi4_dz          = dprof1_dpsi4 * current_FB_fact

end if

temperature_i_profile = temperature_i_profile + Ti_1

return
end subroutine temperature_i
