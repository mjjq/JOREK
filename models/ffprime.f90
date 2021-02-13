recursive subroutine FFprime(xpoint2,xcase2,Z,Z_xpoint,psi,psi_axis,psi_bnd,FFprime_profile, &
                             dFF_dpsi, dFF_dz, &                                             ! 1st order derivatives
                             dFF_dpsi2, dFF_dz2, dFF_dpsi_dz, &                              ! 2nd order derivatives
                             dFF_dpsi3, dFF_dpsi_dz2, dFF_dpsi2_dz,  dFF_dz3, &              ! 2rd order derivatives
                             dFF_dpsi4, dFF_dpsi_dz3, dFF_dpsi2_dz2, dFF_dpsi3_dz, dFF_dz4, &! 4th order derivatives
                             from_F_profile_710)
!-----------------------------------------------------------------------
! Determines the F*F' value and its derivatives at the given
! position (Z, psi) from the analytical or numerical input profile.
!-----------------------------------------------------------------------
use phys_module
use vacuum, only: current_FB_fact
use mod_F_profile
use mod_input_profiles

implicit none

! --- Routine parameters
logical, intent(in)  :: xpoint2
integer, intent(in)  :: xcase2
real*8,  intent(in)  :: Z, Z_xpoint(2), psi, psi_axis, psi_bnd
real*8,  intent(out) :: FFprime_profile
real*8,  intent(out) ::    dFF_dpsi, dFF_dz                                                ! 1st order derivatives
real*8,  intent(out) ::    dFF_dpsi2, dFF_dz2, dFF_dpsi_dz                                 ! 2nd order derivatives
real*8,  intent(out) ::    dFF_dpsi3, dFF_dpsi_dz2, dFF_dpsi2_dz,  dFF_dz3                 ! 2rd order derivatives
real*8,  intent(out) ::    dFF_dpsi4, dFF_dpsi_dz3, dFF_dpsi2_dz2, dFF_dpsi3_dz, dFF_dz4   ! 4th order derivatives
logical, intent(in)  :: from_F_profile_710

! --- Internal variables.
real*8  :: prof1,dprof1_dpsi,dprof1_dpsi2,dprof1_dpsi3,dprof1_dpsi4,dprof1_dpsi5
real*8  :: psi_n, psi_star, delta_psi
real*8  :: atn, datn, d2atn, d3atn, d4atn
real*8  :: d_pert, d2_pert, d3_pert, d4_pert, d5_pert
! for interpolating numerical profiles
integer :: left, right, mid
real*8  :: aux1, aux2
real*8  :: F_prof   
real*8  ::   dF_dpsi, dF_dz                                             ! 1st order derivatives
real*8  ::   dF_dpsi2, dF_dz2, dF_dpsi_dz                               ! 2nd order derivatives
real*8  ::   dF_dpsi3, dF_dpsi_dz2, dF_dpsi2_dz,  dF_dz3                ! 2rd order derivatives
real*8  ::   dF_dpsi4, dF_dpsi_dz3, dF_dpsi2_dz2, dF_dpsi3_dz, dF_dz4   ! 4th order derivatives


! --- the F-profile and FFprime need to be coherent. Always!!!
#ifdef fullmhd
  if (from_F_profile_710) then
    ! --- Call function
    call F_profile(xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd, &
                   F_prof, &
                     dF_dpsi, dF_dz, &                                          ! 1st order derivatives
                     dF_dpsi2, dF_dz2, dF_dpsi_dz, &                            ! 2nd order derivatives
                     dF_dpsi3, dF_dpsi_dz2, dF_dpsi2_dz,  dF_dz3, &             ! 2rd order derivatives
                     dF_dpsi4, dF_dpsi_dz3, dF_dpsi2_dz2, dF_dpsi3_dz, dF_dz4, &! 4th order derivatives
                   FFprime_profile, &
                     dFF_dpsi, dFF_dz, &                                             ! 1st order derivatives
                     dFF_dpsi2, dFF_dz2, dFF_dpsi_dz, &                              ! 2nd order derivatives
                     dFF_dpsi3, dFF_dpsi_dz2, dFF_dpsi2_dz,  dFF_dz3, &              ! 2rd order derivatives
                     dFF_dpsi4, dFF_dpsi_dz3, dFF_dpsi2_dz2, dFF_dpsi3_dz, dFF_dz4)  ! 4th order derivatives

    ! --- Because JOREK uses a negative FF' in the GS-equation and the current-routines
    ! --- But in Full-MHD, because we need to integrate FF', we can't do this, so we use the real F-profile and FF',
    ! --- and then reverse it for all the routines that use it.
    FFprime_profile = - FFprime_profile
    dFF_dpsi        = - dFF_dpsi
    dFF_dpsi2       = - dFF_dpsi2
    dFF_dpsi3       = - dFF_dpsi3
    dFF_dpsi4       = - dFF_dpsi4
    dFF_dz          = - dFF_dz
    dFF_dz2         = - dFF_dz2
    dFF_dz3         = - dFF_dz3
    dFF_dz4         = - dFF_dz4
    dFF_dpsi_dz     = - dFF_dpsi_dz
    dFF_dpsi_dz2    = - dFF_dpsi_dz2
    dFF_dpsi_dz3    = - dFF_dpsi_dz3
    dFF_dpsi2_dz    = - dFF_dpsi2_dz
    dFF_dpsi2_dz2   = - dFF_dpsi2_dz2
    dFF_dpsi3_dz    = - dFF_dpsi3_dz

    return
  endif
#endif


delta_psi = psi_bnd - psi_axis
psi_n     = (psi - psi_axis) / delta_psi
psi_n = max( min(psi_n, 2.), 0. )

! --- Profile as a function of Psi_N.
if ( .not. num_ffprime ) then ! use analytical representation
  
  call input_profiles_psi_component(psi,psi_axis,psi_bnd, FF_0, FF_1, FF_coef, &
                                    prof1,dprof1_dpsi,dprof1_dpsi2,dprof1_dpsi3,dprof1_dpsi4,dprof1_dpsi5)
  
  call input_profiles_edge_perturbation(psi,psi_axis,psi_bnd, FF_coef,    &
                                        d_pert,d2_pert,d3_pert,d4_pert,d5_pert)
  
  ! --- note: tanh is included inside both prof1 and d_pert
  prof1        = prof1        + d_pert
  dprof1_dpsi  = dprof1_dpsi  + d2_pert
  dprof1_dpsi2 = dprof1_dpsi2 + d3_pert
  dprof1_dpsi3 = dprof1_dpsi3 + d4_pert
  dprof1_dpsi4 = dprof1_dpsi4 + d5_pert
  
else ! use numerical representation.
  
  ! --- Interpolate profile and derivatives to position psi_n by bisections.
  left  = 1
  right = num_ffprime_len
  do
    if ( right == left + 1 ) exit
    mid = (left + right) / 2
    if ( num_ffprime_x(mid) >= psi_n ) then
      right = mid
    else
      left = mid
    end if
  end do
  aux1 = (psi_n - num_ffprime_x(left)) / (num_ffprime_x(right) - num_ffprime_x(left))
  aux2 = (1. - aux1)
  prof1        = num_ffprime_y0(left)   * aux2 + num_ffprime_y0(right) * aux1
  dprof1_dpsi  = ( num_ffprime_y1(left) * aux2 + num_ffprime_y1(right) * aux1 ) / delta_psi
  dprof1_dpsi2 = ( num_ffprime_y2(left) * aux2 + num_ffprime_y2(right) * aux1 ) / delta_psi**2
  dprof1_dpsi3 = ( num_ffprime_y3(left) * aux2 + num_ffprime_y3(right) * aux1 ) / delta_psi**3
  dprof1_dpsi4 = ( num_ffprime_y4(left) * aux2 + num_ffprime_y4(right) * aux1 ) / delta_psi**4
  
end if

! compute z-tanh (will be 1 with zero deivatives if xpoint=.f.)
call input_profiles_Z_component(xpoint2,xcase2,Z,Z_xpoint, atn,datn,d2atn,d3atn,d4atn)  

FFprime_profile  = prof1        * atn
dFF_dpsi         = dprof1_dpsi  * atn
dFF_dpsi2        = dprof1_dpsi2 * atn
dFF_dpsi3        = dprof1_dpsi3 * atn
dFF_dpsi4        = dprof1_dpsi4 * atn
dFF_dz           = prof1        * datn
dFF_dz2          = prof1        * d2atn
dFF_dz3          = prof1        * d3atn
dFF_dz4          = prof1        * d4atn
dFF_dpsi_dz      = dprof1_dpsi  * datn
dFF_dpsi_dz2     = dprof1_dpsi  * d2atn
dFF_dpsi_dz3     = dprof1_dpsi  * d3atn
dFF_dpsi2_dz     = dprof1_dpsi2 * datn
dFF_dpsi2_dz2    = dprof1_dpsi2 * d2atn
dFF_dpsi3_dz     = dprof1_dpsi3 * datn


if (freeboundary_equil .and. num_ffprime) then            !if the ffprime profile is given in a file and freeboundary equilibrium is on,
                                                         !the full profile is multiplied by a factor in order to iterate to a given current   
   FFprime_profile  = FFprime_profile  * current_FB_fact
   dFF_dpsi         = dFF_dpsi         * current_FB_fact
   dFF_dpsi2        = dFF_dpsi2        * current_FB_fact
   dFF_dpsi3        = dFF_dpsi3        * current_FB_fact
   dFF_dpsi4        = dFF_dpsi4        * current_FB_fact
   dFF_dz           = dFF_dz           * current_FB_fact
   dFF_dz2          = dFF_dz2          * current_FB_fact
   dFF_dz3          = dFF_dz3          * current_FB_fact
   dFF_dz4          = dFF_dz4          * current_FB_fact
   dFF_dpsi_dz      = dFF_dpsi_dz      * current_FB_fact
   dFF_dpsi_dz2     = dFF_dpsi_dz2     * current_FB_fact
   dFF_dpsi_dz3     = dFF_dpsi_dz3     * current_FB_fact
   dFF_dpsi2_dz     = dFF_dpsi2_dz     * current_FB_fact
   dFF_dpsi2_dz2    = dFF_dpsi2_dz2    * current_FB_fact
   dFF_dpsi3_dz     = dFF_dpsi3_dz     * current_FB_fact

end if

FFprime_profile = FFprime_profile + FF_1

return
end subroutine FFprime
