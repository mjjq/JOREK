!> Module containing base type for field interpolations, with interfaces
!> to implement
module mod_fields
use data_structure
implicit none
private
public fields_base

!> Base type for a field interpolator.
!> Must implement the following interfaces, which are the normal
!> functions and an additional time component (JOREK units)
!> node_list and element_list should be the currently-valid representation of the grid
!> (values themselves should not be used, only for find_RZ etc)
type, abstract :: fields_base
  type(type_node_list)    :: node_list !< Current node list
  type(type_element_list) :: element_list !< Current element list
  contains
    procedure(interp_PRZ), deferred       :: interp_PRZ
    procedure :: calc_EBpsiU
end type fields_base

interface
  !> Interpolate a variable at s, t, phi in i_elm, returning first
  !> derivatives of the variable and of space
  pure subroutine interp_PRZ(this, time, i_elm, i_v, n_v, s, t, phi, P, P_s, P_t, P_phi, P_time, R, R_s, R_t, Z, Z_s, Z_t)
    import fields_base
    class(fields_base),  intent(in)  :: this
    real*8,                   intent(in)  :: time !< Time at which to calculate this variable
    integer,                  intent(in)  :: i_elm
    integer,                  intent(in)  :: n_v, i_v(n_v)
    real*8,                   intent(in)  :: s, t, phi
    real*8,                   intent(out) :: P(n_v), P_s(n_v), P_t(n_v), P_time(n_v)
    real*8,                   intent(out) :: R, R_s, R_t, Z, Z_s, Z_t
    real*8,                   intent(out) :: P_phi(n_v)
  end subroutine interp_PRZ
end interface

contains
!> Calculates the electric and magnetic fields at a specific position
!> in the jorek element `i_elm` at `st`.
pure subroutine calc_EBpsiU(fields, time, i_elm, st, phi, E, B, psi, U)
use phys_module, only: F0, mode, central_mass, central_density
use constants, only: mu_zero, mass_proton

! Routine parameters
class(fields_base), intent(in) :: fields
real*8, intent(in)  :: time
integer, intent(in) :: i_elm !< JOREK element index
real*8, intent(in)  :: st(2) !< element-local coordinates
real*8, intent(in)  :: phi !< toroidal angle
real*8, intent(out) :: E(3) !< Electric field [V/m]
real*8, intent(out) :: B(3) !< Magnetic field [T]
real*8, intent(out) :: psi !< psi in JOREK units (~-0.9 at axis, -0.2 at edge)
real*8, intent(out) :: u !< velocity stream function in m/s

! Internal parameters
integer, parameter :: i_var(2) = [1,2]
real*8             :: P(2), P_s(2), P_t(2), P_phi(2), P_time(2) ! Placeholder for evaluating variables and derivatives locally
! Values
real*8             :: R, R_s, R_t, Z, Z_s, Z_t
! Others
real*8             :: inv_st_jac, R_inv
real*8             :: psi_R, psi_Z, U_R, U_Z, U_phi, t_norm

t_norm  = sqrt(mu_zero * mass_proton * central_mass * central_density * 1.d20) ! 1 jorek time unit in seconds

! Interpolate the fields to get psi and U at the current position (and the
! changes u_n - u(n-1))
call fields%interp_PRZ(time, i_elm, i_var, 2, st(1), st(2), phi, P, P_s, P_t, P_phi, P_time, R, R_s, R_t, Z, Z_s, Z_t)

R_inv = 1.d0/R
inv_st_jac = 1.d0/(R_s * Z_t - R_t * Z_s)

! Calculate the derivatives to R and Z
psi_R    = (  P_s(1) * Z_t - P_t(1) * Z_s ) * inv_st_jac
psi_Z    = (- P_s(1) * R_t + P_t(1) * R_s ) * inv_st_jac
U_R      = (  P_s(2) * Z_t - P_t(2) * Z_s ) * inv_st_jac
U_Z      = (- P_s(2) * R_t + P_t(2) * R_s ) * inv_st_jac
U_phi    = P_phi(2)

! Update psi and U
psi = P(1)
U   = P(2)/t_norm

! Calculate the magnetic field (see http://jorek.eu/wiki/doku.php?id=reduced_mhd)
B     = [+psi_Z, -psi_R, F0] * R_inv
! The local electric field, obtained from E=-Grad (u F0)-\partial_t A
! See http://jorek.eu/wiki/doku.php?id=u_phi
E     = [-F0*U_R, -F0*U_Z, -F0*U_phi*R_inv]/t_norm
E(3)  = E(3) - R*P_time(1) ! because this is not normalized with t_norm
end subroutine calc_EBpsiU
end module mod_fields
