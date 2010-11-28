subroutine current(xpoint2,R,Z,Z_xpoint,psi,psi_axis,psi_bnd,zjz)
!-----------------------------------------------------------------------
! Determine the current at a given position from the density,
! temperature, and FF' input profiles
!-----------------------------------------------------------------------

use phys_module

implicit none

! --- input variables
logical, intent(in)    :: xpoint2
real*8,  intent(in)    :: R, Z
real*8,  intent(in)    :: Z_xpoint
real*8,  intent(in)    :: psi
real*8,  intent(in)    :: psi_axis
real*8,  intent(in)    :: psi_bnd
real*8,  intent(out)   :: zjz        ! Current at the given position.

! --- local variables
real*8  :: psi_n
real*8  :: zn,dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2,dn_dpsi2_dz
real*8  :: zT,dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2,dT_dpsi2_dz
real*8  :: zFFprime, dFFprime_dpsi, dFFprime_dz, dFFprime_dpsi_dz,dFFprime_dpsi2,dFFprime_dz2

psi_n = (psi - psi_axis) / (psi_bnd - psi_axis)

call density(    xpoint2, Z, Z_xpoint, psi,psi_axis,psi_bnd,&
             zn,dn_dpsi,dn_dz,dn_dpsi2, dn_dz2, dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2,dn_dpsi2_dz)

call temperature(xpoint2, Z, Z_xpoint, psi,psi_axis,psi_bnd, &
                 zT,dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2,dT_dpsi2_dz)

call FFprime(    xpoint2, Z, Z_xpoint, psi,psi_axis,psi_bnd, &
             zFFprime,dFFprime_dpsi,dFFprime_dz,dFFprime_dpsi2,dFFprime_dz2, dFFprime_dpsi_dz)

zjz   = zFFprime - R*R * (zn * dT_dpsi + dn_dpsi * zT)

return
end subroutine current
