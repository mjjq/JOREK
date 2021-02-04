subroutine current(xpoint2,xcase2,R,Z,Z_xpoint,psi,psi_axis,psi_bnd,zjz)
!-----------------------------------------------------------------------
! Determine the current at a given position from the density,
! temperature, and FF' input profiles
!-----------------------------------------------------------------------

use mod_parameters
use phys_module

implicit none

! --- Routine parameters
logical, intent(in)    :: xpoint2
integer, intent(in)    :: xcase2
real*8,  intent(in)    :: R, Z
real*8,  intent(in)    :: Z_xpoint(2)
real*8,  intent(in)    :: psi
real*8,  intent(in)    :: psi_axis
real*8,  intent(in)    :: psi_bnd
real*8,  intent(out)   :: zjz        ! Current at the given position.

! --- local variables
real*8  :: psi_n
real*8  :: zn
real*8  ::    dn_dpsi, dn_dz                                           ! 1st order derivatives
real*8  ::    dn_dpsi2, dn_dz2, dn_dpsi_dz                             ! 2nd order derivatives
real*8  ::    dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz,  dn_dz3              ! 2rd order derivatives
real*8  ::    dn_dpsi4, dn_dpsi_dz3, dn_dpsi2_dz2, dn_dpsi3_dz, dn_dz4 ! 4th order derivatives
real*8  :: zT
real*8  ::    dT_dpsi,  dT_dz                                          ! 1st order derivatives
real*8  ::    dT_dpsi2, dT_dz2, dT_dpsi_dz                             ! 2nd order derivatives
real*8  ::    dT_dpsi3, dT_dpsi_dz2, dT_dpsi2_dz,  dT_dz3              ! 2rd order derivatives
real*8  ::    dT_dpsi4, dT_dpsi_dz3, dT_dpsi2_dz2, dT_dpsi3_dz, dT_dz4 ! 4th order derivatives
real*8  :: zTi
real*8  ::    dTi_dpsi,  dTi_dz                                             ! 1st order derivatives
real*8  ::    dTi_dpsi2, dTi_dz2, dTi_dpsi_dz                               ! 2nd order derivatives
real*8  ::    dTi_dpsi3, dTi_dpsi_dz2, dTi_dpsi2_dz,  dTi_dz3               ! 2rd order derivatives
real*8  ::    dTi_dpsi4, dTi_dpsi_dz3, dTi_dpsi2_dz2, dTi_dpsi3_dz, dTi_dz4 ! 4th order derivatives
real*8  :: zTe
real*8  ::    dTe_dpsi,  dTe_dz                                             ! 1st order derivatives
real*8  ::    dTe_dpsi2, dTe_dz2, dTe_dpsi_dz                               ! 2nd order derivatives
real*8  ::    dTe_dpsi3, dTe_dpsi_dz2, dTe_dpsi2_dz,  dTe_dz3               ! 2rd order derivatives
real*8  ::    dTe_dpsi4, dTe_dpsi_dz3, dTe_dpsi2_dz2, dTe_dpsi3_dz, dTe_dz4 ! 4th order derivatives
real*8  :: zFFprime
real*8  ::    dFF_dpsi, dFF_dz                                                ! 1st order derivatives
real*8  ::    dFF_dpsi2, dFF_dz2, dFF_dpsi_dz                                 ! 2nd order derivatives
real*8  ::    dFF_dpsi3, dFF_dpsi_dz2, dFF_dpsi2_dz,  dFF_dz3                 ! 2rd order derivatives
real*8  ::    dFF_dpsi4, dFF_dpsi_dz3, dFF_dpsi2_dz2, dFF_dpsi3_dz, dFF_dz4   ! 4th order derivatives

psi_n = (psi - psi_axis) / (psi_bnd - psi_axis)

call density(    xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd, zn,&
                 dn_dpsi, dn_dz, &                                        ! 1st order derivatives
                 dn_dpsi2, dn_dz2, dn_dpsi_dz, &                          ! 2nd order derivatives
                 dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz,  dn_dz3, &           ! 2rd order derivatives
                 dn_dpsi4, dn_dpsi_dz3, dn_dpsi2_dz2, dn_dpsi3_dz, dn_dz4)! 4th order derivatives

if (with_TiTe) then
  
  call temperature_i(xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd, zTi, &
                     dTi_dpsi,  dTi_dz, &                                          ! 1st order derivatives
                     dTi_dpsi2, dTi_dz2, dTi_dpsi_dz, &                            ! 2nd order derivatives
                     dTi_dpsi3, dTi_dpsi_dz2, dTi_dpsi2_dz,  dTi_dz3, &            ! 2rd order derivatives
                     dTi_dpsi4, dTi_dpsi_dz3, dTi_dpsi2_dz2, dTi_dpsi3_dz, dTi_dz4)! 4th order derivatives
  call temperature_e(xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd, zTe, &
                     dTe_dpsi,  dTe_dz, &                                          ! 1st order derivatives
                     dTe_dpsi2, dTe_dz2, dTe_dpsi_dz, &                            ! 2nd order derivatives
                     dTe_dpsi3, dTe_dpsi_dz2, dTe_dpsi2_dz,  dTe_dz3, &            ! 2rd order derivatives
                     dTe_dpsi4, dTe_dpsi_dz3, dTe_dpsi2_dz2, dTe_dpsi3_dz, dTe_dz4)! 4th order derivatives

  zT = zTi + zTe
  dT_dpsi = dTi_dpsi + dTe_dpsi

else
  
  call temperature(xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd, zT, &
                   dT_dpsi,  dT_dz, &                                       ! 1st order derivatives
                   dT_dpsi2, dT_dz2, dT_dpsi_dz, &                          ! 2nd order derivatives
                   dT_dpsi3, dT_dpsi_dz2, dT_dpsi2_dz,  dT_dz3, &           ! 2rd order derivatives
                   dT_dpsi4, dT_dpsi_dz3, dT_dpsi2_dz2, dT_dpsi3_dz, dT_dz4)! 4th order derivatives

endif

call FFprime(    xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd, zFFprime, &
                 dFF_dpsi, dFF_dz, &                                             ! 1st order derivatives
                 dFF_dpsi2, dFF_dz2, dFF_dpsi_dz, &                              ! 2nd order derivatives
                 dFF_dpsi3, dFF_dpsi_dz2, dFF_dpsi2_dz,  dFF_dz3, &              ! 2rd order derivatives
                 dFF_dpsi4, dFF_dpsi_dz3, dFF_dpsi2_dz2, dFF_dpsi3_dz, dFF_dz4, &! 4th order derivatives
                 .true.)

zjz   = zFFprime - R*R * (zn * dT_dpsi + dn_dpsi * zT)

!if ((bootstrap) .and. (restart)) then
!  zjz   = zjz * (0.5d0 - 0.5d0* tanh( (psi_n - (FF_coef(7)-FF_coef(8)))/FF_coef(8) ) )
!endif

return
end subroutine current
