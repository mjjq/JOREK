subroutine neutral_density(xpoint2,xcase2,Z,Z_xpoint,psi,psi_axis,psi_bnd,density_profile, &
                           dn_dpsi, dn_dz, &                                        ! 1st order derivatives
                           dn_dpsi2, dn_dz2, dn_dpsi_dz, &                          ! 2nd order derivatives
                           dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz,  dn_dz3, &           ! 2rd order derivatives
                           dn_dpsi4, dn_dpsi_dz3, dn_dpsi2_dz2, dn_dpsi3_dz, dn_dz4)! 4th order derivatives
!-----------------------------------------------------------------------
! Determines the neutral density value and its derivatives at the given
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
real*8,  intent(out) :: dn_dpsi, dn_dz
real*8,  intent(out) :: dn_dpsi2, dn_dz2, dn_dpsi_dz
real*8,  intent(out) :: dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz,  dn_dz3
real*8,  intent(out) :: dn_dpsi4, dn_dpsi_dz3, dn_dpsi2_dz2, dn_dpsi3_dz, dn_dz4

! --- Internal variables.
real*8  :: prof0, prof1, dprof0_dpsi, dprof0_dpsi2, dprof0_dpsi3, dprof0_dpsi4, psi_barrier
real*8  :: psi_n, psi_star, delta_psi, dprof1_dpsi, dprof1_dpsi2, dprof1_dpsi3, dprof1_dpsi4, dprof1_dpsi5
real*8  :: atn_z,   datn_z,   d2atn_z, d3atn_z, d4atn_z
real*8  :: atn_z_u, datn_z_u, d2atn_z_u, d3atn_z_u, d4atn_z_u
real*8  :: Ztan_pos
! for interpolating numerical profiles
integer :: left, right, mid
real*8  :: aux1, aux2, Z_star, Z_star_u

delta_psi = psi_bnd - psi_axis
psi_n     = (psi - psi_axis) / delta_psi
psi_n     = max( min(psi_n, 2.), 0. )

! --- Profile as a function of Psi_N.
if ( .not. num_rho ) then ! use analytical representation
  
  call input_profiles_psi_component(psi,psi_axis,psi_bnd, rhon_0, rhon_1, rhon_coef, &
                                    prof1,dprof1_dpsi,dprof1_dpsi2,dprof1_dpsi3,dprof1_dpsi4,dprof1_dpsi5)
  
else ! use numerical representation.
  
  ! --- Interpolate profile and derivatives to position psi_n by bisections.
  left  = 1
  right = num_rhon_len
  do
    if ( right == left + 1 ) exit
    mid = (left + right) / 2
    if ( num_rhon_x(mid) >= psi_n ) then
      right = mid
    else
      left = mid
    end if
  end do
  aux1 = (psi_n - num_rhon_x(left)) / (num_rhon_x(right) - num_rhon_x(left))
  aux2 = (1. - aux1)
  prof1        = num_rhon_y0(left)   * aux2 + num_rhon_y0(right) * aux1
  dprof1_dpsi  = ( num_rhon_y1(left) * aux2 + num_rhon_y1(right) * aux1 ) / delta_psi
  dprof1_dpsi2 = ( num_rhon_y2(left) * aux2 + num_rhon_y2(right) * aux1 ) / delta_psi**2
  dprof1_dpsi3 = ( num_rhon_y3(left) * aux2 + num_rhon_y3(right) * aux1 ) / delta_psi**3
  dprof1_dpsi4 = ( num_rhon_y4(left) * aux2 + num_rhon_y4(right) * aux1 ) / delta_psi**4
  
end if

! --- Additional explicit dependence of the profile on Z to ensure that the profile is
!     approximately zero in the private flux region below the x-point.
!     ***************
!     ***************
!     ***************
!     IMPORTANT NOTE: Everything is the same as the usual density profile, except for
!                     this Z-tanh at the X-points. Instead of setting the profile to rhon_1
!                     in the divertor regions, we set it to the coefs
!                     rhon_coef(9) in the lower divertor
!                     and rhon_coef(10) in the upper divertor
!                     For the tanh width, rhon_coef(8) is used for both divertors
!                     For the tanh Z-position, rhon_coef(6) and rhon_coef(7) are used (lower and upper respectively)
!                     If rhon_coef(6) and/or rhon_coef(7) are 0.d0, then Z_xpoint(:) are used instead.
  
! compute z-tanh (will be 1 with zero deivatives if xpoint=.f.)
! we trick this to get both lower and upper tanh for different private regions levels of neutrals
call input_profiles_Z_component(xpoint2,xcase2,Z,(/Z_xpoint(1),+999./), atn_z,  datn_z,  d2atn_z,  d3atn_z,  d4atn_z)  
call input_profiles_Z_component(xpoint2,xcase2,Z,(/-999.,Z_xpoint(2)/), atn_z_u,datn_z_u,d2atn_z_u,d3atn_z_u,d4atn_z_u)  

density_profile = prof1        + (1.0 -   atn_z) * (rhon_coef(9)-rhon_1 -  prof1      ) + (1.0 -   atn_z_u) * (rhon_coef(10)-rhon_1 -  prof1      )
dn_dpsi         = dprof1_dpsi  + (1.0 -   atn_z) * (                    - dprof1_dpsi ) + (1.0 -   atn_z_u) * (                     - dprof1_dpsi )
dn_dpsi2        = dprof1_dpsi2 + (1.0 -   atn_z) * (                    - dprof1_dpsi2) + (1.0 -   atn_z_u) * (                     - dprof1_dpsi2)
dn_dpsi3        = dprof1_dpsi3 + (1.0 -   atn_z) * (                    - dprof1_dpsi3) + (1.0 -   atn_z_u) * (                     - dprof1_dpsi3)
dn_dpsi4        = dprof1_dpsi4 + (1.0 -   atn_z) * (                    - dprof1_dpsi4) + (1.0 -   atn_z_u) * (                     - dprof1_dpsi4)
dn_dz           = prof1        + (    -  datn_z) * (rhon_coef(9)-rhon_1 -  prof1      ) + (    -  datn_z_u) * (rhon_coef(10)-rhon_1 -  prof1      )
dn_dz2          = prof1        + (    - d2atn_z) * (rhon_coef(9)-rhon_1 -  prof1      ) + (    - d2atn_z_u) * (rhon_coef(10)-rhon_1 -  prof1      )
dn_dz3          = prof1        + (    - d3atn_z) * (rhon_coef(9)-rhon_1 -  prof1      ) + (    - d3atn_z_u) * (rhon_coef(10)-rhon_1 -  prof1      )
dn_dz4          = prof1        + (    - d4atn_z) * (rhon_coef(9)-rhon_1 -  prof1      ) + (    - d4atn_z_u) * (rhon_coef(10)-rhon_1 -  prof1      )
dn_dpsi_dz      = dprof1_dpsi  + (    -  datn_z) * (                    - dprof1_dpsi ) + (    -  datn_z_u) * (                     - dprof1_dpsi )
dn_dpsi_dz2     = dprof1_dpsi  + (    - d2atn_z) * (                    - dprof1_dpsi ) + (    - d2atn_z_u) * (                     - dprof1_dpsi )
dn_dpsi_dz3     = dprof1_dpsi  + (    - d3atn_z) * (                    - dprof1_dpsi ) + (    - d3atn_z_u) * (                     - dprof1_dpsi )
dn_dpsi2_dz     = dprof1_dpsi2 + (    -  datn_z) * (                    - dprof1_dpsi2) + (    -  datn_z_u) * (                     - dprof1_dpsi2)
dn_dpsi2_dz2    = dprof1_dpsi2 + (    - d2atn_z) * (                    - dprof1_dpsi2) + (    - d2atn_z_u) * (                     - dprof1_dpsi2)
dn_dpsi3_dz     = dprof1_dpsi3 + (    -  datn_z) * (                    - dprof1_dpsi3) + (    -  datn_z_u) * (                     - dprof1_dpsi3)


density_profile = density_profile + rhon_1

return
end subroutine neutral_density
