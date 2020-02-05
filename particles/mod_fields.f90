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
  type(type_node_list), allocatable    :: node_list !< Current node list
  type(type_element_list), allocatable :: element_list !< Current element list
  logical(kind=1) :: flag_zero_dpsidt=.false. !< if true, P_time(1) = dpsi/dt = 0
  contains
    procedure(interp_PRZ),deferred,public   :: interp_PRZ
    procedure(interp_PRZ_2),deferred,public :: interp_PRZ_2
    procedure :: calc_EBpsiU
    procedure :: set_flag_dpsidt !< set the flag_zero_dpsidt
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
  !> interface for interpolation routine returning value, first and second order
  !> derivatives of a variable (excluded second derivatives of phi).
  !> Time derivatives are calculated only for variables and their first order
  !> derivatives on s and t. Values, first and second order derivatives of
  !> the global coordinates R,Z are also returned.
  pure subroutine interp_PRZ_2(this,time,i_elm,i_v,n_v,s,t,phi,P,P_s,P_t,P_phi,&
       P_time,P_ss,P_st,P_tt,P_sphi,P_tphi,P_stime,P_ttime,&
       R,R_s,R_t,R_ss,R_st,R_tt,Z,Z_s,Z_t,Z_ss,Z_st,Z_tt)
    import fields_base
    !> declare input variables
    class(fields_base), intent(in) :: this
    real(kind=8),intent(in) :: time,s,t,phi
    integer,intent(in) :: i_elm,n_v
    integer,dimension(n_v),intent(in) :: i_v
    !> declare ourput variables
    real(kind=8),intent(out) :: R,R_s,R_t,R_ss,R_st,R_tt
    real(kind=8),intent(out) :: Z,Z_s,Z_t,Z_ss,Z_st,Z_tt
    real(kind=8),dimension(n_v),intent(out) :: P,P_s,P_t,P_phi,P_time
    real(kind=8),dimension(n_v),intent(out) :: P_ss,P_st,P_tt,P_sphi,P_tphi
    real(kind=8),dimension(n_v),intent(out) :: P_stime,P_ttime
  end subroutine interp_PRZ_2
end interface

contains
!> Calculates the electric and magnetic fields at a specific position
!> in the jorek element `i_elm` at `st`.
pure subroutine calc_EBpsiU(fields, time, i_elm, st, phi, E, B, psi, U)
use phys_module, only: F0, mode, central_mass, central_density
use constants, only: mu_zero, mass_proton
use mod_pusher_tools, only: transform_derivatives_st_to_RZ
! Routine parameters
class(fields_base), intent(in) :: fields
real*8, intent(in)  :: time
integer, intent(in) :: i_elm !< JOREK element index
real*8, intent(in)  :: st(2) !< element-local coordinates
real*8, intent(in)  :: phi !< toroidal angle
real*8, intent(out) :: E(3) !< Electric field [V/m]
real*8, intent(out) :: B(3) !< Magnetic field [T]
real*8, intent(out) :: psi !< psi in JOREK units
real*8, intent(out) :: u !< velocity stream function in m/s

! Internal parameters
integer, parameter :: i_var(2) = [1,2]
real*8             :: P(2), P_s(2), P_t(2), P_phi(2), P_time(2) ! Placeholder for evaluating variables and derivatives locally
! Values
real*8             :: R, R_s, R_t, Z, Z_s, Z_t
! Others
real*8             :: R_inv,t_norm

t_norm  = sqrt(mu_zero * mass_proton * central_mass * central_density * 1.d20) ! 1 jorek time unit in seconds

! Interpolate the fields to get psi and U at the current position (and the
! changes u_n - u(n-1))
call fields%interp_PRZ(time, i_elm, i_var, 2, st(1), st(2), phi, P, P_s, P_t, P_phi, P_time, R, R_s, R_t, Z, Z_s, Z_t)
R_inv = 1.d0/R

! calculate the derivatives to R and Z
call transform_derivatives_st_to_RZ(P_s,P_t,2,P_s,P_t,R_s,R_t,Z_s,Z_t)

! Update psi and U
psi = P(1)
U   = P(2)/t_norm

! set the poloidal magnetic flux time derivative to zero if true
if(fields%flag_zero_dpsidt) P_time(1) = 0.d0

! Calculate the magnetic field (see http://jorek.eu/wiki/doku.php?id=reduced_mhd)
B     = [P_t(1), -P_s(1), F0] * R_inv
! The local electric field, obtained from E=-Grad (u F0)-\partial_t A
! See http://jorek.eu/wiki/doku.php?id=u_phi
E     = -F0*[P_s(2),P_t(2),P_phi(2)*R_inv]/t_norm
E(3)  = E(3) - R_inv*P_time(1) ! because this is not normalized with t_norm

end subroutine calc_EBpsiU

!> This procedure computes the fields required for resolving
!> the guiding center equations of motions
!> inputs:
!>   fields: (field_base) structure containing methods for computing EM fields
!>   time:   (real8) particle time
!>   i_elm:  (integer) particle mesh element index
!>   st:     (real8) particle position in local mesh coordinates
!>   phi:    (real8) particle toroidal angle
!> outputs:
!>   E:      (real8)(3) electric field in V/m
!>   b:      (real8)(3) magnetic field direction
!>   normB:  (real8)(3) magnetic field intensity in T
!>   gradB:  (real8)(3) gradient of the magnetic field intensity in T/m
!>   curlb:  (real8)(3) curl of the magnetic field direction in 1/m
!>   dbdt:   (real8)(3) magnetic field direction time derivative 1/s
pure subroutine calc_EBNormBGradBCurlbDbdt(fields,time,i_elm,st,phi,E,b,&
     normB,gradB,curlb,dbdt)
  !> load modules
  use phys_module, only: F0, mode, central_mass, central_density
  use constants, only: mu_zero,mass_proton
  use mod_pusher_tools, only: right_handed_cross_product 
  use mod_pusher_tools, only: transform_first_derivatives_st_to_RZ
  use mod_pusher_tools, only: transform_second_derivatives_st_to_RZ
  implicit none

  !> declare input variables
  class(fields_base),intent(in) :: fields
  real(kind=8),intent(in) :: time
  integer,intent(in) :: i_elm
  real(kind=8),dimension(2),intent(in) :: st
  real(kind=8),intent(in) :: phi
  !> declare output variables
  real(kind=8),intent(out) :: normB
  real(kind=8),dimension(3),intent(out) :: E,b,gradB,curlb,dbdt
  !> declare parameters
  !> declare internal variables
  real(kind=8) :: R_inv,normB_inv
  real(kind=8) :: R,R_s,R_t,R_ss,R_st,R_tt
  real(kind=8) :: Z,Z_s,Z_t,Z_ss,Z_st,Z_tt
  real(kind=8),dimension(5) :: U !< stream funciton: [U,U_R,U_Z,U_phi,U_time]
  !> poloidal flux: psi,psi_R,psi_Z,psi_phi,psi_time,psi_RR,psi_RZ,psi_ZZ,
  !>   psi_Rphi, psi_Zphi, psi_Rtime,psi_Ztime
  real(kind=8),dimension(12) :: psi 

  !> interpolate the stream function
  call fields%interp_PRZ(time,i_elm,[2],1,st(1),st(2),phi,U(1),U(2),U(3),&
       U(4),U(5),R,R_s,R_t,Z,Z_s,Z_t)
  !> normalise the electric field in SI units
  U = F0*U/sqrt(mu_zero*mass_proton*central_mass*central_density*1.d20)
  !> transform first U derivatives from st to RZ
  call transform_first_derivatives_st_to_RZ(U(2),U(3),1,[U(2)],[U(3)],R_s,R_t,Z_s,Z_t)
  
  !> interpolate the poloidal flux
  call fields%interp_PRZ_2(time,i_elm,[1],1,st(1),st(2),phi,psi(1),psi(2),&
       psi(3),psi(4),psi(5),psi(6),psi(7),psi(8),psi(9),psi(10),psi(11),&
       psi(12),R,R_s,R_t,R_ss,R_st,R_tt,Z,Z_s,Z_t,Z_ss,Z_st,Z_tt)
  !> set dpsidt to zero if needed
  if(fields%flag_zero_dpsidt) then
     psi(5)  = 0.d0 !< psi_time
     psi(11) = 0.d0 !< psi_stime
     psi(12) = 0.d0 !< psi_ttime
  endif
  R_inv = 1.d0/R !< compute the inverse of R
  !> transform first order psi derivatives from st to RZ
  call transform_first_derivatives_st_to_RZ(psi(2),psi(3),1,psi(2),psi(3),R_s,R_t,Z_s,Z_t)
  call transform_first_derivatives_st_to_RZ(psi(9),psi(10),1,psi(9),psi(10),R_s,R_t,Z_s,Z_t)
  call transform_first_derivatives_st_to_RZ(psi(11),psi(12),1,psi(11),psi(12),R_s,R_t,Z_s,Z_t)
  !> transform second order psi derivatives from st to RZ
  call transform_second_derivatives_st_to_RZ(psi(6),psi(7),psi(8),1,psi(6),psi(7),psi(8),&
       psi(2),psi(3),R_s,R_t,R_ss,R_st,R_tt,Z_s,Z_t,Z_ss,Z_st,Z_tt)
  
  !> compute the electric field
  E = -[U(2),U(3),R_inv*(U(4)+psi(5))] !< V/m
  !> compute the magnetic field
  b = [psi(3),-psi(2),F0]*R_inv !< magnetic field T
  normB = sqrt(b(1)*b(1)+b(2)*b(2)+b(3)*b(3)) !< B field intensity
  normB_inv = 1.d0/normB !< inverse of the B field intensity
  !< direction of the magnetic field
  b = b/normB

  !> compute the dbdt field
  dbdt = ((b(2)*psi(11)-b(1)*psi(12))*b +&
       [psi(12),-psi(11),0.d0])*normb_inv*R_inv

  !> compute the gradB field
  gradB = [psi(2)*psi(6)+psi(3)*psi(7),&
       psi(2)*psi(7)+psi(3)*psi(8),&
       R_inv*(psi(2)*psi(9)+psi(3)*psi(10))]*R_inv*R_inv*normB_inv
  gradB(1) = gradB(1)-normB*R_inv
  
  !> TODO compute the curlb field
  curlb = normB_inv*(right_handed_cross_product(b,gradB) + &
       R_inv*[R_inv*psi(9),R_inv*psi(10),R_inv*psi(2)-psi(6)-psi(8)])
  
end subroutine calc_EBNormBGradBCurlbDbdt

!> This procedure set a flag for setting to zero the poloidal
!> magneticflux time derivative
!> inputs:
!>   this: (fields_base) object of class field
!>   flag_dpsidt_to_zero: (logical1) if true  the poloidal
!>                        magnetic flux time derivative is
!>                        set to zero
!> outputs:
!>   this: (fields_base) object of class field
pure subroutine set_flag_dpsidt(this,flag_dpsidt_to_zero)
  class(fields_base),intent(inout) :: this !< fields object
  logical(kind=1),intent(in) :: flag_dpsidt_to_zero !< flag to set

  !> set the flag_zero_dpsidt flag of this
  this%flag_zero_dpsidt = flag_dpsidt_to_zero
  
end subroutine set_flag_dpsidt

end module mod_fields
