!====MB===============parallel velocity profile which is kept by the // velocity source implemented in element_matrix.f90
subroutine velocity(xpoint2,xcase2,Z,Z_xpoint,psi,psi_axis,psi_bnd,velocity_profile, &
                   dV_dpsi, dV_dz, &                                        ! 1st order derivatives
                   dV_dpsi2, dV_dz2, dV_dpsi_dz, &                          ! 2nd order derivatives
                   dV_dpsi3, dV_dpsi_dz2, dV_dpsi2_dz,  dV_dz3, &           ! 2rd order derivatives
                   dV_dpsi4, dV_dpsi_dz3, dV_dpsi2_dz2, dV_dpsi3_dz, dV_dz4)! 4th order derivatives
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
use phys_module
use mod_input_profiles

implicit none

! --- Routine parameters.
logical, intent(in)   :: xpoint2
integer, intent(in)   :: xcase2
real*8,  intent(in)   :: Z
real*8,  intent(in)   :: Z_xpoint(2)
real*8,  intent(in)   :: psi
real*8,  intent(in)   :: psi_axis
real*8,  intent(in)   :: psi_bnd
real*8,  intent(out)  :: velocity_profile
real*8,  intent(out)  :: dV_dpsi, dV_dz
real*8,  intent(out)  :: dV_dpsi2, dV_dz2, dV_dpsi_dz
real*8,  intent(out)  :: dV_dpsi3, dV_dpsi_dz2, dV_dpsi2_dz,  dV_dz3
real*8,  intent(out)  :: dV_dpsi4, dV_dpsi_dz3, dV_dpsi2_dz2, dV_dpsi3_dz, dV_dz4

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
if ( .not. num_rot ) then ! use analytical representation

  call input_profiles_psi_component(psi,psi_axis,psi_bnd, V_0, V_1, V_coef, &
                                    prof1,dprof1_dpsi,dprof1_dpsi2,dprof1_dpsi3,dprof1_dpsi4,dprof1_dpsi5)

else ! use numerical respresentation
!---------------------------------------------------------------------------------------------------------------

  left  = 1
  right = num_rot_len
  do
    if ( right == left + 1 ) exit
    mid = (left + right) / 2
    if ( num_rot_x(mid) >= psi_n ) then
      right = mid
    else
      left = mid
    end if
  end do
  aux1 = (psi_n - num_rot_x(left)) / (num_rot_x(right) - num_rot_x(left))
  aux2 = (1. - aux1)
  prof1        = num_rot_y0(left)   * aux2 + num_rot_y0(right) * aux1
  dprof1_dpsi  = ( num_rot_y1(left) * aux2 + num_rot_y1(right) * aux1 ) / delta_psi
  dprof1_dpsi2 = ( num_rot_y2(left) * aux2 + num_rot_y2(right) * aux1 ) / delta_psi**2
  dprof1_dpsi3 = ( num_rot_y3(left) * aux2 + num_rot_y3(right) * aux1 ) / delta_psi**3
  dprof1_dpsi4 = ( num_rot_y4(left) * aux2 + num_rot_y4(right) * aux1 ) / delta_psi**4
  
end if


! compute z-tanh (will be 1 with zero deivatives if xpoint=.f.)
call input_profiles_Z_component(xpoint2,xcase2,Z,Z_xpoint, atn,datn,d2atn,d3atn,d4atn)  

velocity_profile    = prof1        *   atn
dV_dpsi             = dprof1_dpsi  *   atn
dV_dpsi2            = dprof1_dpsi2 *   atn
dV_dpsi3            = dprof1_dpsi3 *   atn
dV_dpsi4            = dprof1_dpsi4 *   atn
dV_dpsi5            = dprof1_dpsi5 *   atn
dV_dz               = prof1        *  datn
dV_dz2              = prof1        * d2atn
dV_dz3              = prof1        * d3atn
dV_dz4              = prof1        * d4atn
dV_dpsi_dz          = dprof1_dpsi  *  datn
dV_dpsi_dz2         = dprof1_dpsi  * d2atn
dV_dpsi_dz3         = dprof1_dpsi  * d3atn
dV_dpsi2_dz         = dprof1_dpsi2 *  datn
dV_dpsi2_dz2        = dprof1_dpsi2 * d2atn
dV_dpsi3_dz         = dprof1_dpsi3 *  datn


if ( .not. num_rot ) then 
  velocity_profile = velocity_profile + V_1
end if

return
end subroutine velocity
!============================================Marina 14.02.2011================

