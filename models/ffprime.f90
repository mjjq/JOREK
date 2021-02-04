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
real*8  :: prof0, prof1, dprof0_dpsi, dprof0_dpsi2, dprof0_dpsi3, dprof0_dpsi4, psi_barrier
real*8  :: psi_n, psi_star, delta_psi, sig_F, sigz, dprof1_dpsi, dprof1_dpsi2, dprof1_dpsi3, dprof1_dpsi4
real*8  :: atn, datn, d2atn, d3atn, d4atn
real*8  :: atn_z,   datn_z,   d2atn_z, d3atn_z, d4atn_z
real*8  :: atn_z_u, datn_z_u, d2atn_z_u, d3atn_z_u, d4atn_z_u
real*8  :: cosh1, sinh1, cosh2, sinh2, cosh3, sinh3, cosh4, sinh4, cosh3_u, sinh3_u, cosh4_u, sinh4_u
real*8  :: tanh1, tanh2, tanh2_u
real*8  :: d_pert, d2_pert, d3_pert, d4_pert, d5_pert
! for interpolating numerical profiles
integer :: left, right, mid
real*8  :: aux1, aux2, Z_star, Z_star_u
real*8  :: F_prof   
real*8  ::   dF_dpsi, dF_dz                                             ! 1st order derivatives
real*8  ::   dF_dpsi2, dF_dz2, dF_dpsi_dz                               ! 2nd order derivatives
real*8  ::   dF_dpsi3, dF_dpsi_dz2, dF_dpsi2_dz,  dF_dz3                ! 2rd order derivatives
real*8  ::   dF_dpsi4, dF_dpsi_dz3, dF_dpsi2_dz2, dF_dpsi3_dz, dF_dz4   ! 4th order derivatives
real*8  :: no_delta_psi, ffprime_out, FFprime_profile2, F_prof2


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
no_delta_psi = 1.d0
if (FF_coef(9) .eq. 1.d0) no_delta_psi = delta_psi

psi_n = max( min(psi_n, 2.), 0. )

! --- Profile as a function of Psi_N.
if ( .not. num_ffprime ) then ! use analytical representation
  
  d_pert  = + FF_coef(6)/cosh((psi_n - FF_coef(7))/FF_coef(8))**2 / (2.d0 * FF_coef(8)) / delta_psi * no_delta_psi
  d2_pert = - FF_coef(6)/cosh((psi_n - FF_coef(7))/FF_coef(8))**2 / (FF_coef(8)**2)  &
            * tanh((psi_n - FF_coef(7))/FF_coef(8)) / delta_psi**2 * no_delta_psi
  d3_pert = + FF_coef(6)/cosh((psi_n - FF_coef(7))/FF_coef(8))**4 / (FF_coef(8)**3)  &
            * (-2.d0 + cosh(2.d0*(psi_n-FF_coef(7))/FF_coef(8)) ) / delta_psi**3 * no_delta_psi
  d4_pert = + 4.0 * FF_coef(6)/cosh((psi_n - FF_coef(7))/FF_coef(8))**5 / (FF_coef(8)**4) * sinh((psi_n - FF_coef(7))/FF_coef(8)) &
            * (-2.d0 + cosh(2.d0*(psi_n-FF_coef(7))/FF_coef(8)) ) / delta_psi**4 * no_delta_psi &
            + 4.0 * FF_coef(6)/cosh((psi_n - FF_coef(7))/FF_coef(8))**5 / (FF_coef(8)**4)  &
            * (-2.d0 * sinh(2.d0*(psi_n-FF_coef(7))/FF_coef(8)) ) / delta_psi**4 * no_delta_psi
  d5_pert = + 20.0 * FF_coef(6)/cosh((psi_n - FF_coef(7))/FF_coef(8))**6 / (FF_coef(8)**5) * sinh((psi_n - FF_coef(7))/FF_coef(8))**2 &
            * (-2.d0 + cosh(2.d0*(psi_n-FF_coef(7))/FF_coef(8)) ) / delta_psi**5 * no_delta_psi &
            + 4.0 * FF_coef(6)/cosh((psi_n - FF_coef(7))/FF_coef(8))**5 / (FF_coef(8)**5) * cosh((psi_n - FF_coef(7))/FF_coef(8)) &
            * (-2.d0 + cosh(2.d0*(psi_n-FF_coef(7))/FF_coef(8)) ) / delta_psi**5 * no_delta_psi &
            + 20.0 * FF_coef(6)/cosh((psi_n - FF_coef(7))/FF_coef(8))**6 / (FF_coef(8)**5) * sinh((psi_n - FF_coef(7))/FF_coef(8))**2 &
            * (-2.d0 * sinh(2.d0*(psi_n-FF_coef(7))/FF_coef(8)) ) / delta_psi**5 * no_delta_psi &
            + 4.0 * FF_coef(6)/cosh((psi_n - FF_coef(7))/FF_coef(8))**5 / (FF_coef(8)**5)  &
            * (-4.d0 * cosh(2.d0*(psi_n-FF_coef(7))/FF_coef(8)) ) / delta_psi**5 * no_delta_psi

  
  prof0        = (FF_0 - FF_1) * ( 1.d0 + FF_coef(1) * psi_n + FF_coef(2) * psi_n**2 + FF_coef(3) * psi_n**3)
  dprof0_dpsi  = (FF_0 - FF_1) * ( FF_coef(1) + 2.d0 * FF_coef(2) * psi_n + 3.d0 * FF_coef(3) * psi_n**2)    / delta_psi
  dprof0_dpsi2 = (FF_0 - FF_1) * (2.d0 * FF_coef(2) + 6.d0 * FF_coef(3) * psi_n) / delta_psi**2
  dprof0_dpsi3 = (FF_0 - FF_1) * (6.d0 * FF_coef(3))                             / delta_psi**3
  dprof0_dpsi4 = 0.d0
  
  prof0        = prof0        + d_pert
  dprof0_dpsi  = dprof0_dpsi  + d2_pert
  dprof0_dpsi2 = dprof0_dpsi2 + d3_pert
  dprof0_dpsi3 = dprof0_dpsi3 + d4_pert
  dprof0_dpsi4 = dprof0_dpsi4 + d5_pert
  
  sig_F        = FF_coef(4)
  psi_barrier  = FF_coef(5)
  
  psi_star = (psi_n - psi_barrier)/sig_F
  psi_star = min( max( psi_star, -40.d0), 40.d0) ! avoid floating-point exceptions
  
  tanh1 = tanh(psi_star)
  cosh1 = cosh(psi_star)
  sinh1 = sinh(psi_star)
  cosh2 = cosh(2.d0*psi_star)
  sinh2 = sinh(2.d0*psi_star)
  
  atn   = (0.5d0 - 0.5d0*tanh1)
  datn  = - 1.d0/cosh1**2 / (2.d0 * sig_F) / delta_psi
  d2atn =   1.d0/cosh1**2 / sig_F**2 * tanh1 / delta_psi**2
  d3atn = - 1.d0/cosh1**4 / sig_F**3 * (-2.d0 + cosh2) / delta_psi**3
  d4atn = - 4.d0/cosh1**5 / sig_F**4 * (-2.d0 + cosh2) / delta_psi**4 * sinh1 &
          - 1.d0/cosh1**4 / sig_F**4 * (-2.d0 * sinh2) / delta_psi**4
  
  prof1        = prof0        * atn
  dprof1_dpsi  = dprof0_dpsi  * atn +         prof0       * datn
  dprof1_dpsi2 = dprof0_dpsi2 * atn + 2.d0 * dprof0_dpsi  * datn + prof0       * d2atn
  dprof1_dpsi3 = dprof0_dpsi3 * atn + 3.d0 * dprof0_dpsi2 * datn + 3.d0 * dprof0_dpsi * d2atn + prof0 * d3atn
  dprof1_dpsi4 = dprof0_dpsi4 * atn + 4.d0 * dprof0_dpsi3 * datn + 6.d0 * dprof0_dpsi2 * d2atn & 
                 + 6.d0 * dprof0_dpsi * d3atn + prof0 * d4atn
  
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

! --- Additional explicit dependence of the profile on Z to ensure that the profile is
!     approximately zero in the private flux region below the x-point.
if ( xpoint2 ) then
  
  sigz = 0.1d0

  if (xcase2 .eq. 1) then
    atn_z_u   = 1.d0
    datn_z_u  = 0.d0
    d2atn_z_u = 0.d0
    d3atn_z_u = 0.d0
    d4atn_z_u = 0.d0
  else
    Z_star_u  = (Z-Z_xpoint(2))/sigz
    Z_star_u  = min( max( Z_star_u, -40.d0), 40.d0) ! avoid floating-point exceptions
    
    tanh2_u   = tanh(Z_star_u)
    cosh3_u   = cosh(Z_star_u)
    sinh3_u   = sinh(Z_star_u)
    cosh4_u   = cosh(2.0*Z_star_u)
    sinh4_u   = sinh(2.0*Z_star_u)

    atn_z_u   = (0.5d0 - 0.5d0*tanh2_u)
    datn_z_u  = -0.5d0/cosh3_u**2 / sigz
    d2atn_z_u =  1.0d0/cosh3_u**2 / sigz**2 * tanh2_u
    d3atn_z_u = -1.0d0/cosh3_u**4 / sigz**3 * (-2.d0 + cosh4_u) 
    d4atn_z_u = -4.0d0/cosh3_u**5 / sigz**4 * (-2.d0 + cosh4_u) * sinh3_u &
                -1.0d0/cosh3_u**4 / sigz**4 * (-2.d0 * sinh4_u) 
  endif
  if (xcase2 .eq. 2) then
    atn_z   = 1.d0
    datn_z  = 0.d0
    d2atn_z = 0.d0
    d3atn_z = 0.d0
    d4atn_z = 0.d0
  else
    Z_star  = (Z_xpoint(1)-Z)/sigz
    Z_star  = min( max( Z_star, -40.d0), 40.d0) ! avoid floating-point exceptions

    tanh2   = tanh(Z_star)
    cosh3   = cosh(Z_star)
    sinh3   = sinh(Z_star)
    cosh4   = cosh(2.0*Z_star)
    sinh4   = sinh(2.0*Z_star)
      
    atn_z   = (0.5d0 - 0.5d0*tanh2)
    datn_z  =  0.5d0/cosh3**2   / sigz
    d2atn_z =  1.0d0/cosh3**2   / sigz**2 * tanh2
    d3atn_z = -1.0d0/cosh3**4 / sigz**3 * (-2.d0 + cosh4) 
    d4atn_z = -4.0d0/cosh3**5 / sigz**4 * (-2.d0 + cosh4) * sinh3 &
              -1.0d0/cosh3**4 / sigz**4 * (-2.d0 * sinh4) 
  endif 
  
  FFprime_profile  = prof1        *    atn_z * atn_z_u

  dFF_dpsi         = dprof1_dpsi  *    atn_z * atn_z_u
  dFF_dpsi2        = dprof1_dpsi2 *    atn_z * atn_z_u
  dFF_dpsi3        = dprof1_dpsi3 *    atn_z * atn_z_u  
  dFF_dpsi4        = dprof1_dpsi4 *    atn_z * atn_z_u  
     
  dFF_dz           = prof1        * ( datn_z * atn_z_u +          atn_z * datn_z_u)
  dFF_dz2          = prof1        * (d2atn_z * atn_z_u + 2.d0 *  datn_z * datn_z_u        +   atn_z * d2atn_z_u)  
  dFF_dz3          = prof1        * (d3atn_z * atn_z_u + 4.d0 * d2atn_z * datn_z_u + 4.d0 *  datn_z * d2atn_z_u  + atn_z * d3atn_z_u)  
  dFF_dz4          = prof1        * (d4atn_z * atn_z_u + 4.d0 * d3atn_z * datn_z_u + 6.d0 * d2atn_z * d2atn_z_u &
                                     + 4.d0 * datn_z * d3atn_z_u + atn_z * d4atn_z_u)
     
  dFF_dpsi_dz      = dprof1_dpsi  * ( datn_z * atn_z_u +         atn_z * datn_z_u)
  dFF_dpsi_dz2     = dprof1_dpsi  * (d2atn_z * atn_z_u + 2.d0 * datn_z * datn_z_u  + atn_z * d2atn_z_u)  
  dFF_dpsi_dz3     = dprof1_dpsi  * (d3atn_z * atn_z_u + 4.d0 * d2atn_z * datn_z_u + 4.d0 *  datn_z * d2atn_z_u  + atn_z * d3atn_z_u)
     
  dFF_dpsi2_dz     = dprof1_dpsi2 * ( datn_z * atn_z_u +         atn_z * datn_z_u)
  dFF_dpsi2_dz2    = dprof1_dpsi2 * (d2atn_z * atn_z_u + 2.d0 * datn_z * datn_z_u  + atn_z * d2atn_z_u)
     
  dFF_dpsi3_dz     = dprof1_dpsi3 * ( datn_z * atn_z_u +         atn_z * datn_z_u)

else
 
  FFprime_profile  = prof1
  dFF_dpsi         = dprof1_dpsi
  dFF_dpsi2        = dprof1_dpsi2
  dFF_dpsi3        = dprof1_dpsi3
  dFF_dpsi4        = dprof1_dpsi4
  dFF_dz           = 0.d0
  dFF_dz2          = 0.d0
  dFF_dz3          = 0.d0
  dFF_dz4          = 0.d0
  dFF_dpsi_dz      = 0.d0
  dFF_dpsi_dz2     = 0.d0
  dFF_dpsi_dz3     = 0.d0
  dFF_dpsi2_dz     = 0.d0
  dFF_dpsi2_dz2    = 0.d0
  dFF_dpsi3_dz     = 0.d0
  
end if

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
