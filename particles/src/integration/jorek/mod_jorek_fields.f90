module mod_jorek_fields
use mod_fields
use data_structure
use mod_action
use mod_particle_sim
implicit none

type, extends(fields) :: jorek_fields
  type(type_node_list), pointer    :: node_list
  type(type_element_list), pointer :: element_list
  real*8 :: F0
  real*8 :: tstep
  ! parameters?
  logical :: static = .false.
  contains
    procedure :: at_particle => at_particle_impl
end type jorek_fields

type, extends(action) :: read_jorek_fields
  contains
    procedure :: do => read_jorek_fields_impl
end type read_jorek_fields

private

contains

!> Read jorek fields from a restart file
subroutine read_jorek_fields_impl(this, sim)
  !use mod_import_export_restart
  class(read_jorek_fields), intent(inout) :: this
  type(particle_sim), intent(inout) :: sim
  write(*,*) "TODO: read JOREK fields"
end subroutine read_jorek_fields_impl

pure subroutine at_particle_impl(this, particle, t, E, B, psi, U)
  use mod_particle_base
  class(jorek_fields),  intent(in)  :: this
  class(particle_base), intent(in)  :: particle
  real*8, intent(in)                :: t
  real*8, dimension(3), intent(out) :: E, B
  real*8, intent(out), optional     :: psi, U

  if (present(psi) .and. present(U)) then
    call calc_EB(particle%i_elm, particle%st, particle%x(3), E, B, psi, U, 0.d0)
  end if
end subroutine at_particle_impl

!> Calculates the electric and magnetic fields at a specific position
!> in the jorek element `i_elm` at `st`.
!> Linear interpolation with element%deltas is performed according to
!> `delta_fraction`, which starts at 1 and goes to 0 for no mixing.
!> If it is 1 we get the fields of the previous timesteps.
pure subroutine calc_EB(i_elm,st,phi,E,B,psi,U,delta_fraction)
use parameters
use constants
use phys_module, only : F0, tstep


interface
  pure subroutine interp_PRZ(node_list, element_list, i_elm, i_v, n_v, &
          s, t, phi, P, P_s, P_t, P_phi, R, R_s, R_t, Z, Z_s, Z_t)
    import :: type_node_list, type_element_list
    type (type_node_list),    intent(in)  :: node_list
    type (type_element_list), intent(in)  :: element_list
    integer,                  intent(in)  :: i_elm
    integer,                  intent(in)  :: n_v, i_v(n_v)
    real*8,                   intent(in)  :: s, t, phi
    real*8,                   intent(out) :: P(n_v), P_s(n_v), P_t(n_v)
    real*8,                   intent(out) :: R, R_s, R_t, Z, Z_s, Z_t
    real*8,                   intent(out) :: P_phi(n_v)
  end subroutine interp_PRZ
  pure subroutine interp_PRZ_delta(node_list, element_list, i_elm, i_v, n_v, &
          s, t, phi, P, P_s, P_t, P_phi, R, R_s, R_t, Z, Z_s, Z_t)
    import :: type_node_list, type_element_list
    type (type_node_list),    intent(in)  :: node_list
    type (type_element_list), intent(in)  :: element_list
    integer,                  intent(in)  :: i_elm
    integer,                  intent(in)  :: n_v, i_v(n_v)
    real*8,                   intent(in)  :: s, t, phi
    real*8,                   intent(out) :: P(n_v), P_s(n_v), P_t(n_v)
    real*8,                   intent(out) :: R, R_s, R_t, Z, Z_s, Z_t
    real*8,                   intent(out) :: P_phi(n_v)
  end subroutine interp_PRZ_delta
end interface

! Routine parameters
integer, intent(in) :: i_elm !< JOREK element index
real*8, intent(in)  :: st(2) !< element-local coordinates
real*8, intent(in)  :: phi !< toroidal angle
real*8, intent(in), optional :: delta_fraction !< linear interpolation factor. goes from 1 to 0 in time

real*8, intent(out) :: E(3), B(3), psi, U !< Fields and potentials

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

! Use the current timestep size (very rough)
psi_time = Pd(1)/tstep

! Calculate the magnetic field (see http://jorek.eu/wiki/doku.php?id=reduced_mhd)
B     = (/ + psi_Z, - psi_R, F0 /) * R_inv
! The local electric field, obtained from E=-Grad (u F0)-\partial_t A
! See http://jorek.eu/wiki/doku.php?id=u_phi
E     = (/ - F0 * U_R, - F0 * U_Z, - F0 * U_phi * R_inv - R * psi_time /)
end subroutine calc_EB
end module mod_jorek_fields
