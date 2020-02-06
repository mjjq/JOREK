!> Particle pusher module for integrating full orbits of relativistic
!> particles. Available integrators:
!> -5-steps Volume Preserving Algorithm (VPA): R. Zhang et al., Phys. Plasmas 22 (2015) 044501 
!>                                       (see also C. Sommariva et al., Nucl. Fusion 58 (2018) 016043)
module mod_kinetic_relativistic
use mod_particle_types
use constants, only: EL_CHG,ATOMIC_MASS_UNIT,SPEED_OF_LIGHT

implicit none

private
public volume_preserving_push_cartesian,volume_preserving_push_jorek
public relativistic_kinetic_to_gc
public gc_to_relativistic_kinetic

contains

!---------------------------------------------------------------------------
!> This subroutine computes the first half-step of the JOREK-specific VPA
!> inputs:
!>   half_position: real(8)(3) initial particle position in cartesian coordinates  
!>   particle: (particle_kinetic_relativistic) relativistic particle
!>   mass:     (real8) particle mass in [AMU]
!>   dt:       (real8) time step in [s]
!> outputs:
!>   scaling_factor: (real8)(3) scaling factor to be used in subsequent steps
!>   half_position:  (real8)(2) particle position after half-step 
!>			        in cartesian coordinates
pure subroutine volume_preserving_first_half_step_jorek(particle,half_position,&
       mass,dt,scaling_factor)
  ! input variables
  real(kind=8), intent(in) :: mass, dt !< mass and time step
  ! input/output variables
  real(kind=8),dimension(3),intent(inout) :: half_position
  class(particle_kinetic_relativistic), intent(inout) :: particle !< relativistic particle
  ! output variables
  real(kind=8), intent(out) :: scaling_factor

  scaling_factor = 5.d-1*dt*particle%q*EL_CHG/(ATOMIC_MASS_UNIT*mass*SPEED_OF_LIGHT)
  ! compute dimensionless momentum
  particle%p = particle%p/(mass*SPEED_OF_LIGHT)
  ! compute coordinates at half-step
  half_position = half_position + (5.d-1*dt*SPEED_OF_LIGHT*particle%p) &
                  /(sqrt(1.d0+dot_product(particle%p,particle%p)))
end subroutine volume_preserving_first_half_step_jorek

!---------------------------------------------------------------------------
!> This subroutine computes the second half-step of the JOREK-specific VPA
!> inputs:
!>   particle:       (particle_kinetic_relativistic) relativistic particle
!>   half_position:  (real8)(2) half-position in cartesian coordinates
!>   scaling_factor: (real8) scaling factor for computing momentum
!>   B:		     (real8)(3) magnetic field in cartesian coordinates
!>   E:		     (real8)(3) electric field in cartesian coordinates
!>   mass:	     (real8) particle mass in [AMU]
!>   dt:	     (real8) time step in [s]
!> outputs:
!>   particle: (particle_kinetic_relativistic) relativistic particle type
pure subroutine volume_preserving_second_half_step_jorek(particle,&
                half_position,scaling_factor,E,B,mass,dt)
  ! load methods
  use mod_pusher_tools, only: cayley_transform !< use full Cayley transform
  ! define input/output variables
  class(particle_kinetic_relativistic), intent(inout) :: particle !< relativistic particle
  real(kind=8),dimension(3),intent(inout) :: half_position
  ! define input variables
  real(kind=8),dimension(3),intent(in) :: B, E
  real(kind=8),intent(in) :: scaling_factor, mass, dt

  ! compute momentum at t_(i+1/2)
  particle%p = particle%p + scaling_factor*E
  ! rotate momentum with respect to the magnetic field
  particle%p = matmul(cayley_transform(SPEED_OF_LIGHT*scaling_factor/&
    (sqrt(1.d0+dot_product(particle%p,particle%p))),B),particle%p)
  ! compute momentum at t_(i+1)
  particle%p = particle%p + scaling_factor*E
  ! update position at t_(i+1)
  half_position = half_position + (5.d-1*dt*SPEED_OF_LIGHT*particle%p)/&
    (sqrt(1.d0+dot_product(particle%p,particle%p)))
  ! compute dimensional momentum
  particle%p = particle%p*mass*SPEED_OF_LIGHT
end subroutine volume_preserving_second_half_step_jorek

!---------------------------------------------------------------------------
!> This subroutine integrates a relativistic particle trajectory in JOREK
!> fields using the Volume Preserving Algorithm (VPA)
subroutine volume_preserving_push_jorek(particle,fields,mass,time,timestep,ifail)
  ! load functions
  use mod_coordinate_transforms
  use mod_pusher_tools, only: vector_transform_RZPHI_to_XYZ
  use mod_fields
  use mod_find_rz_nearby
  ! declare input/output variables
  integer(kind=4),intent(inout) :: ifail
  class(particle_kinetic_relativistic), intent(inout) :: particle !< relativistic particle
  ! declare input variables
  real(kind=8),intent(in) :: mass, time, timestep
  class(fields_base), intent(in) :: fields
  ! declare internal variables
  real(kind=8) :: psi, U
  real(kind=8),dimension(3) :: B, E
  ! half_position coordinates: 1:x, 2:y, 3:z, 4:R, 5:Z, 6:phi 
  real(kind=8),dimension(6) :: half_position
  real(kind=8) :: scaling_factor !< in [s^2*C/(kg*m)]

  ! check if the particle is valid
  if(particle%i_elm.eq.0) return
  ! transform the particle position from cylindrical to cartesian coordinates
  half_position(1:3) = cylindrical_to_cartesian(particle%x)
  ! compute first half-step
  call volume_preserving_first_half_step_jorek(particle,half_position(1:3),&
       mass,timestep,scaling_factor)
  ! calculate cylindrical coordinates from cartesian ones
  half_position(4:6) = cartesian_to_cylindrical(half_position(1:3))
  ! find the (i_elm,s,t) coordinates
  call find_RZ_nearby(fields%node_list,fields%element_list,particle%x(1),&
       particle%x(2),particle%st(1),particle%st(2),particle%i_elm,&
       half_position(4),half_position(5),particle%st(1),particle%st(2),&
       particle%i_elm,ifail)
  ! check if the particle is lost, exit if it is the case
  if(particle%i_elm.eq.0) return
  ! copy RZPHI coordinates in particles
  call update_particle_position(particle,half_position(4:6))
  ! compute magnetic and electric fields
  call fields%calc_EBpsiU(time+5.d-1*timestep,particle%i_elm,&
       particle%st,particle%x(3),E,B,psi,U)
  ! get B and E fields in cartesian coordinates
  B = vector_transform_RZPHI_to_XYZ(particle%x(3),B)
  E = vector_transform_RZPHI_to_XYZ(particle%x(3),E)
  ! compute the second half-step  
  call volume_preserving_second_half_step_jorek(particle,half_position(1:3),&
       scaling_factor,E,B,mass,timestep)
  ! calculate cylindrical coordinates from cartesian ones
  half_position(4:6) = cartesian_to_cylindrical(half_position(1:3))
  ! find the (i_elm,s,t) coordinates
  call find_RZ_nearby(fields%node_list,fields%element_list,particle%x(1),&
       particle%x(2),particle%st(1),particle%st(2),particle%i_elm,&
       half_position(4),half_position(5),particle%st(1),particle%st(2),&
       particle%i_elm,ifail)
  ! copy new RZPHI position into particle
  call update_particle_position(particle,half_position(4:6))
end subroutine volume_preserving_push_jorek

!---------------------------------------------------------------------------
!> This subroutine implements a test version of the VPA assuming uniform B and E.
!> This subroutine is to be used for tests and comparisons, not for production.
!> inputs:
!>   particle: (particle_kinetic_relativistic) relativistic particle
!>   mass:     (real8) particle mass in [AMU]
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
  real(kind=8), dimension(3), intent(in) :: E, B !< electric and magnetic fields
  ! internal variable
  real(kind=8) :: scaling_factor !< in [s^2*C/(kg*m)]

  scaling_factor = 5.d-1*dt*particle%q*EL_CHG/(ATOMIC_MASS_UNIT*mass*SPEED_OF_LIGHT)

  ! compute dimensionless momentum
  particle%p = particle%p/(mass*SPEED_OF_LIGHT)

  ! compute position at t_(i+1/2)
  particle%x = particle%x + (5.d-1*dt*SPEED_OF_LIGHT*particle%p)/&
    (sqrt(1.d0+dot_product(particle%p,particle%p)))
  
  ! compute momentum at t_(i+1/2)
  particle%p = particle%p + scaling_factor*E
  
  ! rotate momentum with respect to the magnetic field
  particle%p = matmul(cayley_transform(SPEED_OF_LIGHT*scaling_factor/&
    (sqrt(1.d0+dot_product(particle%p,particle%p))),B),particle%p)

  ! compute momentum at t_(i+1)
  particle%p = particle%p + scaling_factor*E

  ! compute position at t_(i+1)
  particle%x = particle%x + (5.d-1*dt*SPEED_OF_LIGHT*particle%p)/&
    (sqrt(1.d0+dot_product(particle%p,particle%p)))

  ! compute dimensional momentum
  particle%p = particle%p*mass*SPEED_OF_LIGHT  
end subroutine volume_preserving_push_cartesian

!---------------------------------------------------------------------------

!> This function transforms a particle_kinetic_relativistic into a
!> particle_gc. Phase space particle coordinates are transformed
!> in guiding center position, energy in [eV] and magnetic moment in [eV/T]
function relativistic_kinetic_to_gc(node_list,element_list,in,B,mass) result(out)
  use data_structure
  use mod_pusher_tools, only: vector_transform_RZPHI_to_XYZ
  ! declare input variables
  type(type_node_list),intent(in) :: node_list
  type(type_element_list),intent(in) :: element_list
  type(particle_kinetic_relativistic),intent(in) :: in
  real(kind=8),intent(in) :: mass !< particle mass in AMU
  real(kind=8),dimension(3),intent(in) :: B !< magnetic field in [T]
  ! delcare output variables
  type(particle_gc) :: out
  ! declare internal variables
  real(kind=8) :: B_norm, p_par !< magnetic intensity, parallel momentum
  real(kind=8),dimension(3) :: B_hat_cart !< magnetic field direction

  ! initialise default variables for particle_gc
  out = in
  ! initialise the electric charge
  out%q = in%q

  ! compute magnetic field intensity and direction
  B_norm = norm2(B) !< intensity
  B_hat_cart = vector_transform_RZPHI_to_XYZ(in%x(3),B)/B_norm  !< direction
  ! compute the parallel momentum
  p_par = dot_product(in%p,B_hat_cart)

  ! compute the guiding center total (i.e. rest+kinetic) energy in [eV]
  out%E = ATOMIC_MASS_UNIT*SPEED_OF_LIGHT* &
    sqrt((mass*SPEED_OF_LIGHT)*(mass*SPEED_OF_LIGHT)+dot_product(in%p,in%p))/EL_CHG
  ! compute the magnetic moment p_perp^2/(2*B) in [eV/T]
  ! the sign is given by the particle parallel momentum 
  out%mu = sign((ATOMIC_MASS_UNIT*SPEED_OF_LIGHT*dot_product(in%p-p_par*B_hat_cart, &
    in%p-p_par*B_hat_cart))/(2.d0*B_norm*mass*EL_CHG),p_par)

  ! compute the gc position
  if(out%q.ne.0) then 
    call relativistic_kinetic_position_to_gc(node_list,element_list,&
    in%x,in%st,in%i_elm,in%p,in%q,B_hat_cart,B_norm,out%x,out%st,out%i_elm)
  endif  
end function relativistic_kinetic_to_gc

!---------------------------------------------------------------------------

!> This function transforms a particle_gc into a particle_kinetic_relativistic.
!> The gc phase space coordinates have to be provided in [eV] for the energy 
!> and [eV/T] for the magnetic moment.
function gc_to_relativistic_kinetic(node_list,element_list,in,chi,B,mass) result(out)
  use data_structure
  use mod_pusher_tools, only: vector_transform_RZPHI_to_XYZ,get_orthonormals
  ! declare input variables
  type(type_node_list),intent(in) :: node_list
  type(type_element_list),intent(in) :: element_list
  type(particle_gc),intent(in) :: in
  real(kind=8),intent(in) :: mass, chi !< particle mass [AMU], gyroangle [rad]
  real(kind=8),dimension(3),intent(in) :: B !< magnetic field in [T]
  ! declare output variables
  type(particle_kinetic_relativistic) :: out
  ! declare internal variables
  real(kind=8) :: B_norm, p_perp, p_par   !< magnetic intensity, perpendicular/parallel momenta
  real(kind=8),dimension(3) :: B_hat_cart, e1_cart, e2_cart !< B-field-aligned cartesian vector basis

  ! copy the default particle datatype
  out = in
  ! copy the particle charge
  out%q = in%q

  ! compute the magnetic field intensity and direction
  B_norm = norm2(B)
  B_hat_cart = B/B_norm
  ! compute a B-field-aligned orthonormal vector basis
  call get_orthonormals(B_hat_cart,e1_cart,e2_cart)
  ! rotate the basis to a XYZ reference
  B_hat_cart = vector_transform_RZPHI_to_XYZ(in%x(3),B_hat_cart)
  e1_cart = vector_transform_RZPHI_to_XYZ(in%x(3),e1_cart)
  e2_cart = vector_transform_RZPHI_to_XYZ(in%x(3),e2_cart)

  ! compute the perpendicular momentum *squared* in (AMU*m/s)^2
  p_perp = (EL_CHG*2.d0*mass*B_norm*sign(in%mu,1.d0))/ATOMIC_MASS_UNIT
  ! compute the parallel momentum in (AMU*m/s)
  p_par = sign(sqrt((((in%E*EL_CHG)*(in%E*EL_CHG))/ &
  ((ATOMIC_MASS_UNIT*SPEED_OF_LIGHT)*(ATOMIC_MASS_UNIT*SPEED_OF_LIGHT)))-&
  (mass*SPEED_OF_LIGHT)*(mass*SPEED_OF_LIGHT)-p_perp),in%mu)
  ! compute the perpendicular momentum in (AMU*m/s)
  p_perp = sqrt(p_perp)

  ! computing the particle momentum
  out%p = p_par*B_hat_cart + p_perp*(e1_cart*cos(chi)+e2_cart*sin(chi)) 

  ! compute the particle position in R,Z,Phi coordinates
  if(out%q.ne.0) then
    call gc_position_to_relativistic_particle(node_list,element_list,&
    in%x,in%st,in%i_elm,out%p,in%q,B_hat_cart,B_norm,out%x,out%st,out%i_elm)
  endif  
end function gc_to_relativistic_kinetic

!---------------------------------------------------------------------------

!> This subroutine computes the relativistic guiding centre coordinates
!> from the particle position and momentum
subroutine relativistic_kinetic_position_to_gc(node_list,element_list,&
x_in,st_in,i_elm_in,p_in,q_in,B_hat_cart,B_norm,x_gc_out,st_gc_out,i_elm_out)
  use data_structure
  use mod_pusher_tools, only: right_handed_cross_product
  use mod_coordinate_transforms
  use mod_find_rz_nearby
  ! declare input variables
  type(type_node_list),intent(in) :: node_list
  type(type_element_list),intent(in) :: element_list
  integer(kind=1),intent(in) :: q_in !< particle charge
  integer(kind=4),intent(in) :: i_elm_in !< particle element
  real(kind=8),dimension(3),intent(in) :: x_in, p_in !< particle position and momentum
  real(kind=8),dimension(2),intent(in) :: st_in !< particle local coordinates
  real(kind=8),dimension(3),intent(in) :: B_hat_cart !< Magnetic field direction B/B_norm
  real(kind=8),intent(in) :: B_norm !< magnetic field intensity in [T]
  ! declare output variables
  integer(kind=4),intent(out) :: i_elm_out !< gc element
  real(kind=8),dimension(2),intent(out) :: st_gc_out !< local gc position s,t
  real(kind=8),dimension(3),intent(out) :: x_gc_out  !< global position in R,Z,phi
  ! declare internal variables
  integer :: ifail !< ifail kind not defined in find_RZ_nearby

  ! compute the guiding center position in XYZ
  x_gc_out = cylindrical_to_cartesian(x_in) +              &
    (ATOMIC_MASS_UNIT*right_handed_cross_product(p_in,B_hat_cart))/ &
    (EL_CHG*real(q_in,8)*B_norm)

  ! transform from XYZ to RZPHI
  x_gc_out = cartesian_to_cylindrical(x_gc_out)  

  ! find the local coordinates
  call find_RZ_nearby(node_list,element_list,x_in(1),x_in(2),&
  st_in(1),st_in(2),i_elm_in,x_gc_out(1),x_gc_out(2),&
  st_gc_out(1),st_gc_out(2),i_elm_out,ifail)
end subroutine relativistic_kinetic_position_to_gc

!---------------------------------------------------------------------------
!> This subroutine computes the relativistic particle coordinates
!> from relativistic gc particle type
subroutine gc_position_to_relativistic_particle(node_list,element_list,&
x_gc_in,st_gc_in,i_elm_in,p_gc_in,q_gc_in,B_hat_cart,B_norm,x_out,st_out,i_elm_out)
  use data_structure
  use mod_pusher_tools, only: right_handed_cross_product
  use mod_coordinate_transforms
  use mod_find_rz_nearby
  ! declare input variables
  type(type_node_list),intent(in) :: node_list
  type(type_element_list),intent(in) :: element_list
  integer(kind=1),intent(in) :: q_gc_in !< gc charge
  integer(kind=4),intent(in) :: i_elm_in !< gc element
  real(kind=8),dimension(3),intent(in) :: x_gc_in,p_gc_in !< gc position and momentum
  real(kind=8),dimension(2),intent(in) :: st_gc_in !< gc local coordinates
  real(kind=8),dimension(3),intent(in) :: B_hat_cart !< Magnetic field direction B/B_norm
  real(kind=8),intent(in) :: B_norm !< magnetic field intensity in [T]
  ! declare output variables
  integer(kind=4),intent(out) :: i_elm_out !< particle element
  real(kind=8),dimension(2),intent(out) :: st_out !< local particle position s,t
  real(kind=8),dimension(3),intent(out) :: x_out  !< global position in R,Z,phi
  ! declare internal variables
  integer :: ifail !< ifail kind not defined in find_RZ_nearby

  ! compute the particle position in XYZ
  x_out = cylindrical_to_cartesian(x_gc_in) + &
    (ATOMIC_MASS_UNIT*right_handed_cross_product(B_hat_cart,p_gc_in))/ &
    (EL_CHG*real(q_gc_in,8)*B_norm)

  ! transform from XYZ to RZPHI
  x_out = cartesian_to_cylindrical(x_out)

  ! find the local coordinates
  call find_RZ_nearby(node_list,element_list,x_gc_in(1),x_gc_in(2),&
  st_gc_in(1),st_gc_in(2),i_elm_in,x_out(1),x_out(2),&
  st_out(1),st_out(2),i_elm_out,ifail)
end subroutine gc_position_to_relativistic_particle

!---------------------------------------------------------------------------

end module mod_kinetic_relativistic
