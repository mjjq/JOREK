!> Calculates the electric and magnetic fields at a specific position
subroutine calc_EB(i_elm,st,phi,E,B,psi,U,delta_fraction)
use parameters
use nodes_elements
use constants
use phys_module, only : F0, central_mass, central_density, tstep_n

implicit none

! Routine parameters
integer, intent(in) :: i_elm
real*8, intent(in)  :: st(2), phi
real*8, intent(in), optional :: delta_fraction

real*8, intent(out) :: E(3), B(3), psi, U

! Internal parameters
integer :: i_var(2)
real*8                    :: P(2), P_s(2), P_t(2), P_phi(2) ! Placeholder for evaluating variables and derivatives locally
! Values
real*8                    :: R, R_s, R_t, Z, Z_s, Z_t
! Deltas
real*8                    :: Pd(2), Pd_s(2), Pd_t(2), Pd_phi(2) ! Placeholder for evaluating variables and derivatives locally
! Others
real*8                    :: inv_st_jac, R_inv
real*8                    :: psi_R, psi_Z, U_R, U_Z, U_phi
real*8                    :: psi_time !< Derivative of psi to time (linearly interpolated from values and deltas)

! Select psi and U
i_var = (/1,2/)

! Interpolate the fields to get psi and U at the current position (and the
! changes u_n - u(n-1))
call       interp_PRZ(node_list,element_list,i_elm,i_var,2,st(1),st(2),phi,P, P_s, P_t, P_phi, R,R_s,R_t,Z,Z_s,Z_t)
call interp_PRZ_delta(node_list,element_list,i_elm,i_var,2,st(1),st(2),phi,Pd,Pd_s,Pd_t,Pd_phi,R,R_s,R_t,Z,Z_s,Z_t)

R_inv = 1.d0/R
inv_st_jac = 1.d0/(R_s * Z_t - R_t * Z_s)

if (present(delta_fraction)) then
  psi_R    = -(  Pd_s(1) * Z_t - Pd_t(1) * Z_s ) * inv_st_jac * delta_fraction
  psi_Z    = -(- Pd_s(1) * R_t + Pd_t(1) * R_s ) * inv_st_jac * delta_fraction
  U_R      = -(  Pd_s(2) * Z_t - Pd_t(2) * Z_s ) * inv_st_jac * delta_fraction
  U_Z      = -(- Pd_s(2) * R_t + Pd_t(2) * R_s ) * inv_st_jac * delta_fraction
  U_phi    = -Pd_phi(2) * delta_fraction
  psi = -Pd(1) * delta_fraction
  U   = -Pd(2) * delta_fraction
else
  psi_R = 0.d0
  psi_Z = 0.d0
  U_R   = 0.d0
  U_Z   = 0.d0
  U_phi = 0.d0
  psi = 0.d0
  U   = 0.d0
endif

! Calculate the derivatives to R and Z (at delta_fraction if presen)
psi_R    = psi_R + (  P_s(1) * Z_t - P_t(1) * Z_s ) * inv_st_jac
psi_Z    = psi_Z + (- P_s(1) * R_t + P_t(1) * R_s ) * inv_st_jac
U_R      = U_R   + (  P_s(2) * Z_t - P_t(2) * Z_s ) * inv_st_jac
U_Z      = U_Z   + (- P_s(2) * R_t + P_t(2) * R_s ) * inv_st_jac

! Update psi and U
psi = psi + P(1)
U   = U   + P(2)

! WARNING: only uses tstep_n(1), does not support variable timestep!
psi_time = Pd(1)/tstep_n(1)

! Calculate the magnetic field (see http://jorek.eu/wiki/doku.php?id=reduced_mhd)
B     = (/ + psi_Z, - psi_R, F0 /) * R_inv
! The local electric field, obtained from E=-Grad (u F0)-\partial_t A
! See http://jorek.eu/wiki/doku.php?id=u_phi
E     = (/ - F0 * U_R, - F0 * U_Z, - F0 * U_phi * R_inv - R * psi_time /)
end subroutine calc_EB
