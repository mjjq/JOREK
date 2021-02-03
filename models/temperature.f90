subroutine temperature(xpoint2,xcase2,Z,Z_xpoint,psi,psi_axis,psi_bnd,temperature_profile, &
                       dT_dpsi,  dT_dz, &                                       ! 1st order derivatives
                       dT_dpsi2, dT_dz2, dT_dpsi_dz, &                          ! 2nd order derivatives
                       dT_dpsi3, dT_dpsi_dz2, dT_dpsi2_dz,  dT_dz3, &           ! 2rd order derivatives
                       dT_dpsi4, dT_dpsi_dz3, dT_dpsi2_dz2, dT_dpsi3_dz, dT_dz4)! 4th order derivatives
!-----------------------------------------------------------------------
! Determines the temperature value and its derivatives at the given
! position (Z, psi) from the analytical or numerical input profile.
!-----------------------------------------------------------------------
use phys_module
use vacuum, only: current_FB_fact

implicit none

! --- Routine parameters
logical, intent(in)  :: xpoint2
integer, intent(in)  :: xcase2
real*8,  intent(in)  :: Z, Z_xpoint(2), psi, psi_axis, psi_bnd
real*8,  intent(out) :: temperature_profile
real*8,  intent(out) :: dT_dpsi,  dT_dz                                          ! 1st order derivatives
real*8,  intent(out) :: dT_dpsi2, dT_dz2, dT_dpsi_dz                             ! 2nd order derivatives
real*8,  intent(out) :: dT_dpsi3, dT_dpsi_dz2, dT_dpsi2_dz,  dT_dz3              ! 2rd order derivatives
real*8,  intent(out) :: dT_dpsi4, dT_dpsi_dz3, dT_dpsi2_dz2, dT_dpsi3_dz, dT_dz4 ! 4th order derivatives

! --- Internal variables.
real*8  :: prof0, prof1, dprof0_dpsi, dprof0_dpsi2, dprof0_dpsi3, dprof0_dpsi4, psi_barrier
real*8  :: psi_n, psi_star, delta_psi, sig_T, sigz, dprof1_dpsi, dprof1_dpsi2, dprof1_dpsi3, dprof1_dpsi4
real*8  :: atn, datn, d2atn, d3atn, d4atn
real*8  :: atn_z,   datn_z,   d2atn_z, d3atn_z, d4atn_z
real*8  :: atn_z_u, datn_z_u, d2atn_z_u, d3atn_z_u, d4atn_z_u
real*8  :: cosh1, sinh1, cosh2, sinh2, cosh3, cosh4, sinh4, cosh3_u, cosh4_u, sinh4_u
real*8  :: tanh1, tanh2, tanh2_u
! for interpolating numerical profiles
integer :: left, right, mid
real*8  :: aux1, aux2, Z_star, Z_star_u

delta_psi = psi_bnd - psi_axis
psi_n     = (psi - psi_axis) / delta_psi

psi_n = max( min(psi_n, 2.), 0. )

! --- Profile as a function of Psi_N.
if ( .not. num_T ) then ! use analytical representation
  
  prof0        = (T_0-T_1)*(1.d0 + T_coef(1) * psi_n + T_coef(2) * psi_n**2 + T_coef(3) * psi_n**3)
  dprof0_dpsi  = (T_0-T_1)*(T_coef(1) + 2.d0 * T_coef(2) * psi_n + 3.d0 * T_coef(3) * psi_n**2) / delta_psi
  dprof0_dpsi2 = (T_0-T_1)*(2.d0 * T_coef(2) + 6.d0 * T_coef(3) * psi_n)                        / delta_psi**2
  dprof0_dpsi3 = (T_0-T_1)*(6.d0 * T_coef(3))                                                   / delta_psi**3
  dprof0_dpsi4 = 0.d0
  
  sig_T       = T_coef(4)
  psi_barrier = T_coef(5)
  
  psi_star = (psi_n - psi_barrier)/sig_T
  psi_star = min( max( psi_star, -40.d0), 40.d0) ! avoid floating-point exceptions
  
  tanh1 = tanh(psi_star)
  cosh1 = cosh(psi_star)
  sinh1 = sinh(psi_star)
  cosh2 = cosh(2.d0*psi_star)
  sinh2 = sinh(2.d0*psi_star)
  
  atn   = (0.5d0 - 0.5d0*tanh1)
  datn  = - 1.d0/cosh1**2 / (2.d0 * sig_T) / delta_psi
  d2atn =   1.d0/cosh1**2 / sig_T**2 * tanh1 / delta_psi**2
  d3atn = - 1.d0/cosh1**4 / sig_T**3 * (-2.d0 + cosh2) / delta_psi**3
  d4atn = - 4.d0/cosh1**5 / sig_T**4 * (-2.d0 + cosh2) / delta_psi**4 * sinh1 &
          - 1.d0/cosh1**4 / sig_T**4 * (-2.d0 * sinh2) / delta_psi**4
  
  prof1        = prof0        * atn
  dprof1_dpsi  = dprof0_dpsi  * atn +         prof0       * datn
  dprof1_dpsi2 = dprof0_dpsi2 * atn + 2.d0 * dprof0_dpsi  * datn + prof0              * d2atn
  dprof1_dpsi3 = dprof0_dpsi3 * atn + 3.d0 * dprof0_dpsi2 * datn + 3.d0 * dprof0_dpsi * d2atn + prof0 * d3atn
  dprof1_dpsi4 = dprof0_dpsi4 * atn + 4.d0 * dprof0_dpsi3 * datn + 6.d0 * dprof0_dpsi2 * d2atn & 
                 + 6.d0 * dprof0_dpsi * d3atn + prof0 * d4atn
  
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
  dprof1_dpsi4 = ( num_T_y4(left) * aux2 + num_T_y4(right) * aux1 ) / delta_psi**4
  
end if

! --- Additional explicit dependence of the profile on Z to ensure that the profile is
!     approximately zero in the private flux region below the x-point.
if (xpoint2) then
  
  sigz            = 0.05d0

  if(xcase2 .eq. 1) then
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
    cosh4_u   = cosh(2.0*Z_star_u)
    sinh4_u   = sinh(2.0*Z_star_u)

    atn_z_u   = (0.5d0 - 0.5d0*tanh2_u)
    datn_z_u  = -0.5d0/cosh3_u**2 / sigz
    d2atn_z_u =  1.0d0/cosh3_u**2 / sigz**2 * tanh2_u
    d3atn_z_u = -1.0d0/cosh3_u**4 / sigz**3 * (-2.d0 + cosh4_u) 
    d4atn_z_u =  4.0d0/cosh3_u**5 / sigz**4 * (-2.d0 + cosh4_u) &
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
    cosh4   = cosh(2.0*Z_star)
    sinh4   = sinh(2.0*Z_star)
      
    atn_z   = (0.5d0 - 0.5d0*tanh2)
    datn_z  =  0.5d0/cosh3**2   / sigz
    d2atn_z =  1.0d0/cosh3**2   / sigz**2 * tanh2
    d3atn_z = -1.0d0/cosh3**4 / sigz**3 * (-2.d0 + cosh4) 
    d4atn_z =  4.0d0/cosh3**5 / sigz**4 * (-2.d0 + cosh4) &
              -1.0d0/cosh3**4 / sigz**4 * (-2.d0 * sinh4) 
  endif  
  
  temperature_profile = prof1      *    atn_z * atn_z_u

  dT_dpsi         = dprof1_dpsi  *    atn_z * atn_z_u
  dT_dpsi2        = dprof1_dpsi2 *    atn_z * atn_z_u
  dT_dpsi3        = dprof1_dpsi3 *    atn_z * atn_z_u  
  dT_dpsi4        = dprof1_dpsi4 *    atn_z * atn_z_u  

  dT_dz           = prof1        * ( datn_z * atn_z_u +          atn_z * datn_z_u)
  dT_dz2          = prof1        * (d2atn_z * atn_z_u + 2.d0 *  datn_z * datn_z_u        +   atn_z * d2atn_z_u)  
  dT_dz3          = prof1        * (d3atn_z * atn_z_u + 4.d0 * d2atn_z * datn_z_u + 4.d0 *  datn_z * d2atn_z_u  + atn_z * d3atn_z_u)  
  dT_dz4          = prof1        * (d4atn_z * atn_z_u + 4.d0 * d3atn_z * datn_z_u + 6.d0 * d2atn_z * d2atn_z_u &
                                    + 4.d0 * datn_z * d3atn_z_u + atn_z * d4atn_z_u)

  dT_dpsi_dz      = dprof1_dpsi  * ( datn_z * atn_z_u +         atn_z * datn_z_u)
  dT_dpsi_dz2     = dprof1_dpsi  * (d2atn_z * atn_z_u + 2.d0 * datn_z * datn_z_u  + atn_z * d2atn_z_u)  
  dT_dpsi_dz3     = dprof1_dpsi  * (d3atn_z * atn_z_u + 4.d0 * d2atn_z * datn_z_u + 4.d0 *  datn_z * d2atn_z_u  + atn_z * d3atn_z_u)

  dT_dpsi2_dz     = dprof1_dpsi2 * ( datn_z * atn_z_u +         atn_z * datn_z_u)
  dT_dpsi2_dz2    = dprof1_dpsi2 * (d2atn_z * atn_z_u + 2.d0 * datn_z * datn_z_u  + atn_z * d2atn_z_u)

  dT_dpsi3_dz     = dprof1_dpsi3 * ( datn_z * atn_z_u +         atn_z * datn_z_u)

else
  
  temperature_profile = prof1
  dT_dpsi             = dprof1_dpsi
  dT_dpsi2            = dprof1_dpsi2
  dT_dpsi3            = dprof1_dpsi3
  dT_dpsi4            = dprof1_dpsi4
  dT_dz               = 0.d0
  dT_dz2              = 0.d0
  dT_dz3              = 0.d0
  dT_dz4              = 0.d0
  dT_dpsi_dz          = 0.d0
  dT_dpsi_dz2         = 0.d0
  dT_dpsi_dz3         = 0.d0
  dT_dpsi2_dz         = 0.d0
  dT_dpsi2_dz2        = 0.d0
  dT_dpsi3_dz         = 0.d0

end if

if (freeboundary_equil .and. num_T) then                        !if the temperature profile is given in a file and there is freeboundary equilibrium
                                                                !the full profile is multiplied by a facto in order to iterate to a given current
  temperature_profile = temperature_profile * current_FB_fact
  dT_dpsi             = dT_dpsi             * current_FB_fact
  dT_dpsi2            = dT_dpsi2            * current_FB_fact
  dT_dpsi3            = dT_dpsi3            * current_FB_fact
  dT_dpsi4            = dT_dpsi4            * current_FB_fact
  dT_dz               = dT_dz               * current_FB_fact
  dT_dz2              = dT_dz2              * current_FB_fact
  dT_dz3              = dT_dz3              * current_FB_fact
  dT_dz4              = dT_dz4              * current_FB_fact
  dT_dpsi_dz          = dT_dpsi_dz          * current_FB_fact
  dT_dpsi_dz2         = dT_dpsi_dz2         * current_FB_fact
  dT_dpsi_dz3         = dT_dpsi_dz3         * current_FB_fact
  dT_dpsi2_dz         = dT_dpsi2_dz         * current_FB_fact
  dT_dpsi2_dz2        = dT_dpsi2_dz2        * current_FB_fact
  dT_dpsi3_dz         = dT_dpsi3_dz         * current_FB_fact

end if

temperature_profile = temperature_profile + T_1

return
end subroutine temperature
