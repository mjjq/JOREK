!> Paritcle pusher module for integrating full orbits of relativistic
!> particles. Available integrators:
!> 5-steps Volume Preserving Integrators: R. Zhang et all, Phys. of Plasmas,
!> vol.22, p.044501 2015
module mod_kinetic_relativistic
use mod_particle_types
! use electric charge, atomic mass unit and speed of ligth
use constants, only: EL_CHG,ATOMIC_MASS_UNIT,SPEED_OF_LIGTH

implicit none

private

public volume_preserving_push_cartesian,relativistic_kinetic_to_gc

contains

!---------------------------------------------------------------------------
!> This subroutine implements a test version of the volume preserving 
!> algorithm: described in R. Zhang, Phys. of Plasmas, vol.22, p.044501, 2015.
!> using the complete Cayley transform in cartesian coordinates.
!> This subroutine has to be used only for tests and comparisons not for
!> production work. WARNING: this pusher works only for constant and
!> uniform E and B.
!> inputs:
!>   particle: (particle_kinetic_relativistic) relativistic particle type
!>   mass:        (real8) particle mass in [AMU]
!>   E:        (real8)(3) electric field in [V/m]
!>   B:        (real8)(3) magnetic field in [T]
!>   dt:       (real8) time step in [s]
!> outputs:
!>   particle: (particle_kinetic_relativistic) relativistic particle type
pure subroutine volume_preserving_push_cartesian(particle,mass,E,B,dt)
  use mod_pusher_tools, only: cayley_transform !< use full Cayley transform
  ! define input output variables
  class(particle_kinetic_relativistic), intent(inout) :: particle !< relativistic particle
  ! define input variables
  real(kind=8), intent(in) :: mass, dt !< mass and time step
  real(kind=8), dimension(3), intent(in) :: E,B !< electric and magnetic field
  ! internal variable
  real(kind=8) :: scaling_factor !< scaling factor [s^2*C/(kg*m)]

  ! compute the dimensional q
  scaling_factor = 5.d-1*dt*particle%q*EL_CHG/(ATOMIC_MASS_UNIT*mass*SPEED_OF_LIGTH)

  ! compute dimensionless momenta
  particle%p = particle%p/(mass*SPEED_OF_LIGTH)

  ! compute position at t_(i+1/2)
  particle%x = particle%x + (5.d-1*dt*SPEED_OF_LIGTH*particle%p)/&
    (sqrt(1.d0+dot_product(particle%p,particle%p)))
  
  ! compute the momenta at t_(i+1/2)
  particle%p = particle%p + scaling_factor*E
  
  ! rotate the momenta with respect to the magnetic field
  particle%p = matmul(cayley_transform(SPEED_OF_LIGTH*scaling_factor/&
    (sqrt(1.d0+dot_product(particle%p,particle%p))),&
    B),particle%p)

  ! update compute momenta at t_(i+1)
  particle%p = particle%p + scaling_factor*E

  ! update particle position at t_(i+1)
  particle%x = particle%x + (5.d-1*dt*SPEED_OF_LIGTH*particle%p)/&
    (sqrt(1.d0+dot_product(particle%p,particle%p)))

  ! compute dimensional momenta
  particle%p = particle%p*mass*SPEED_OF_LIGTH
  
end subroutine volume_preserving_push_cartesian

!---------------------------------------------------------------------------

!> This function transfrom a relativistic kinetic particle in a a
!> particle gc datatype. Phase space particle coordinates are transformed
!> in guiding center position, energy in [eV] and magnetic moment in [eV/T]
!> inputs:
!> outputs:
function relativistic_kinetic_to_gc(node_list,element_list,in,B,mass) result(out)
  use data_structure
  ! declare input variables
  type(type_node_list),intent(in) :: node_list
  type(type_element_list),intent(in) :: element_list
  type(particle_kinetic_relativistic),intent(in) :: in ! input kinetic particle
  real(kind=8),dimension(3),intent(in) :: mass,B !< mass in AMU magnetic field in [T]
  ! delcare output variables
  type(particle_gc) :: out ! output particle guiding center
  ! declare internal variable
  real(kind=8) :: B_norm,p_par !< magnetic intensity, parallel momentum
  real(kind=8),dimension(3) :: B_hat !< magnetic field direction

  ! initialise default variable for particle gc
  out = in
  ! initialise the electric charge
  out%q = in%q

  ! compute magnetic field intensity and direction
  B_norm = norm2(B) !< intensity
  B_hat = B/B_norm  !< direction
  ! compute the parallel momentum
  p_par = dot_product(in%p,B_hat)

  ! compute the guiding center total energy in [eV]
  out%E = ATOMIC_MASS_UNIT*SPEED_OF_LIGTH*(sqrt((mass*SPEED_OF_LIGTH)&
  *(mass*SPEED_OF_LIGTH)+dot_product(in%p,in%p)))/EL_CHG
  ! compute the magnetic moment p_perp^2/(2*B) in [eV/T]
  ! the sign is given by the particle parallel momentum 
  out%mu = sign((ATOMIC_MASS_UNIT*SPEED_OF_LIGTH*dot_product(in%p-p_par*B_hat,&
  in%p-p_par*B_hat))/(2.d0*B_norm*mass*EL_CHG),p_par)

  ! check whether the particle is not a field line
  if(out%q.ne.0) then
    ! compute the gc position
    call relativistic_kinetic_position_to_gc(node_list,element_list,&
    in,B_hat,B_norm,out%x,out%st,out%i_elm)
  endif
  
end function relativistic_kinetic_to_gc

!---------------------------------------------------------------------------

!> This subroutine computes the relativistic gc particle coordinates
!> from relativistic kinetic particle type
!> inputs:
!> outputs:
subroutine relativistic_kinetic_position_to_gc(node_list,element_list,&
in,B_hat,B_norm,x_gc_out,st_gc_out,i_elm_out)
  use data_structure
  use mod_pusher_tools, only: right_handed_cross_product,&
  coordinate_transfrom_RZPHI_to_XYZ,coordinate_transfrom_XYZ_to_RZPHI
  use mod_find_rz_nearby
  ! declare input variables
  type(type_node_list),intent(in) :: node_list
  type(type_element_list),intent(in) :: element_list
  type(particle_kinetic_relativistic),intent(in) :: in
  real(kind=8),dimension(3),intent(in) :: B_hat !< Magnetic field direction B/B_norm
  real(kind=8),intent(in) :: B_norm !< magentic field intensity in [T]
  ! declare output variables
  integer(kind=4),intent(out) :: i_elm_out 
  real(kind=8),dimension(2),intent(out) :: st_gc_out ! local gc position s,t
  real(kind=8),dimension(3),intent(out) :: x_gc_out  !< global position in R,Z,phi
  ! delcare internal variables
  integer :: ifail !< ifail kind not defined in find_RZ_nearby

  ! compute the guiding center position in cartesian reference
  x_gc_out = coordinate_transfrom_RZPHI_to_XYZ(in%x)+&
  (ATOMIC_MASS_UNIT*right_handed_cross_product(in%p,&
  vector_transform_RZPHI_to_XYZ(in%x(3),B_hat)))/&
  (EL_CHG*real(in%q,8)*B_norm)

  ! transform back from a cartesian to a cylindrical coordinate system
  x_gc_out = coordinate_transfrom_XYZ_to_RZPHI(x_gc_out)  

  ! find the guiding center mesh element local coordinates
  call find_RZ_nearby(node_list,element_list,in%x(1),in%x(2),&
  in%st(1),in%st(2),in%i_elm,x_gc_out(1),x_gc_out(2),&
  st_gc_out(1),st_gc_out(2),i_elm_out,ifail)

end subroutine relativistic_kinetic_position_to_gc

!---------------------------------------------------------------------------

end module mod_kinetic_relativistic
