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
public relativistic_kinetic_to_particle
public gc_to_relativistic_kinetic
public relativistic_kinetic_to_relativistic_gc

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
  use mod_coordinate_transforms, only: cartesian_to_cylindrical
  use mod_coordinate_transforms, only: cylindrical_to_cartesian
  use mod_coordinate_transforms, only: vector_cylindrical_to_cartesian
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
  ! compute first half step
  call volume_preserving_first_half_step_jorek(particle,half_position(1:3),&
       mass,timestep,scaling_factor)
  ! transform the particle half position from cartesian to cylindrical coordinates
  half_position(4:6) = cartesian_to_cylindrical(half_position(1:3))
  ! call find RZ for identifting the new local particle position
  call find_RZ_nearby(fields%node_list,fields%element_list,particle%x(1),&
       particle%x(2),particle%st(1),particle%st(2),particle%i_elm,&
       half_position(4),half_position(5),particle%st(1),particle%st(2),&
       particle%i_elm,ifail)
  ! check if the particle is lost, exit if it is the case
  if(particle%i_elm.eq.0) return
  ! copy RZPHI coordinates in particles
  particle%x = half_position(4:6)
  ! compute magnetic and electric field
  call fields%calc_EBpsiU(time+5.d-1*timestep,particle%i_elm,&
       particle%st,particle%x(3),E,B,psi,U)
  ! get B and E fields in Cartesian coordinates
  B = vector_cylindrical_to_cartesian(particle%x(3),B)
  E = vector_cylindrical_to_cartesian(particle%x(3),E)
  ! compute the second half-step  
  call volume_preserving_second_half_step_jorek(particle,half_position(1:3),&
       scaling_factor,E,B,mass,timestep)
  ! transform back from cartisian to cylindrical coordinates
  half_position(4:6) = cartesian_to_cylindrical(half_position(1:3))
  ! call the find RZ for tracking the paritcle in local coordinates
  call find_RZ_nearby(fields%node_list,fields%element_list,particle%x(1),&
       particle%x(2),particle%st(1),particle%st(2),particle%i_elm,&
       half_position(4),half_position(5),particle%st(1),particle%st(2),&
       particle%i_elm,ifail)
  ! copy new RZPHI position into particle
  particle%x = half_position(4:6)
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

!--------------------------------------------------------------------------

!> This procedure transforms a particle_kinetic_relativistic into a different
!> particle type
!> inputs:
!>   node_list:    (type_node_list) a jorek node list
!>   element_list: (type_element_list) a jorek element list
!>   particle_in:  (particle_kinetic_relativistic) a relativistic particle
!>   mass:
!>   B:            (real8)(3) magnetic field
!> outputs:
!>   particle_out: (particle_base) output particle
subroutine relativistic_kinetic_to_particle(node_list,element_list,particle_in,&
     particle_out,mass,B)
  !> load modules
  use data_structure
  implicit none
  !> declare input variables
  type(type_node_list), intent(in)                :: node_list
  type(type_element_list), intent(in)             :: element_list
  type(particle_kinetic_relativistic), intent(in) :: particle_in
  real(kind=8), intent(in)                        :: mass
  real(kind=8), dimension(3), intent(in)          :: B
  !> declare output variables
  class(particle_base), intent(out)               :: particle_out

  !> select the type of particle out
  select type (particle_out)
  type is (particle_gc)
     particle_out = relativistic_kinetic_to_gc(node_list,element_list, &
          particle_in,mass,B)
  type is (particle_gc_relativistic)
     particle_out = relativistic_kinetic_to_relativistic_gc(node_list, &
          element_list,particle_in,mass,B)
  end select
  
end subroutine relativistic_kinetic_to_particle

!--------------------------------------------------------------------------

!> This procedure transforms a particle_kinetic_relativistic 
!> into a particle_gc_relativistic
!> inputs:
!>   node_list:    (type_node_list) jorek nodes
!>   element_list: (type_element_list) jorek mesh elements
!>   in:           (particle_kinetic_relativistic) relativistic particle
!>   mass:         (real8) particle mass
!>   B:            (real8)(3) magnetic field
!> outputs:
!>   out: (particle_gc_relativistic) a relativistic gc 
function relativistic_kinetic_to_relativistic_gc(node_list,element_list,&
     in,mass,B) result(out)
  use data_structure
  use mod_coordinate_transforms, only: vector_cartesian_to_cylindrical
  use mod_pusher_tools, only: particle_position_to_gc
  implicit none
  !> declare input variables
  type(type_node_list), intent(in)                :: node_list
  type(type_element_list), intent(in)             :: element_list
  type(particle_kinetic_relativistic), intent(in) :: in
  real(kind=8), intent(in)                        :: mass
  real(kind=8), dimension(3), intent(in)          :: B
  !> declare output variables
  type(particle_gc_relativistic)                  :: out
  !> delcare internal variables
  real(kind=8)                                    :: norm_B
  real(kind=8), dimension(3)                      :: B_hat, p_perp

  !> compute magnetic field direction and intensity
  norm_B = sqrt(B(1)*B(1)+B(2)*B(2)+B(3)*B(3))
  B_hat = B/norm_B
  !> copy base particle
  out = in
  !> copy charge
  out%q = in%q
  !> extract momenta in cylindrical coordinates
  p_perp = vector_cartesian_to_cylindrical(in%x(3),in%p)
  !> compute parallel momentum
  out%p(1) = p_perp(1)*B_hat(1)+p_perp(2)*B_hat(2)+p_perp(3)*B_hat(3)
  !> compute perpendicular momenta
  p_perp = p_perp - out%p(1)*B_hat
  !> compute magnetic moment
  out%p(2) = (p_perp(1)*p_perp(1)+p_perp(2)*p_perp(2)+&
       p_perp(3)*p_perp(3))/(2.d0*norm_B*mass)
  !> check particle validity
  if(out%q.ne.0) then
     call particle_position_to_gc(node_list,element_list,in%x,&
          in%st,in%i_elm,vector_cartesian_to_cylindrical(in%x(3),in%p),&
          in%q,B_hat,norm_B,out%x,out%st,out%i_elm)
  endif
  
end function relativistic_kinetic_to_relativistic_gc

!---------------------------------------------------------------------------

!> This function transforms a particle_kinetic_relativistic into a
!> particle_gc. Phase space particle coordinates are transformed
!> in guiding center position, energy in [eV] and magnetic moment in [eV/T]
!> inputs:
!>   node_list:    (type_node_list) node list
!>   element_list: (type_element_list) element list
!>   in:           (particle_relativistic_kinetic) a relativistic particle
!>   mass:         (real8) particle mass in AMU
!>   B:            (real8)(3)(optional) magnetic field [T]
!> outputs:
!>   out: (particle_gc) a guiding center particle
function relativistic_kinetic_to_gc(node_list,element_list,in,mass,B) result(out)
  use data_structure
  use mod_coordinate_transforms, only: vector_cartesian_to_cylindrical
  use mod_pusher_tools, only: particle_position_to_gc
  ! declare input variables
  type(type_node_list),intent(in) :: node_list
  type(type_element_list),intent(in) :: element_list
  type(particle_kinetic_relativistic),intent(in) :: in ! input kinetic particle
  real(kind=8),intent(in) :: mass !< particle mass in AMU
  real(kind=8),dimension(3),intent(in) :: B !< magnetic field in [T]
  ! delcare output variables
  type(particle_gc) :: out
  ! declare internal variables
  real(kind=8) :: B_norm, p_par !< magnetic intensity, parallel momentum
  real(kind=8),dimension(3) :: p_perp,B_hat !< magnetic field direction

  ! initialise default variables for particle_gc
  out = in
  ! initialise the electric charge
  out%q = in%q

  ! compute the guiding center total (i.e. rest+kinetic) energy in [eV]
  out%E = ATOMIC_MASS_UNIT*SPEED_OF_LIGHT* &
       sqrt((mass*SPEED_OF_LIGHT)*(mass*SPEED_OF_LIGHT)+&
       (in%p(1)*in%p(1)+in%p(2)*in%p(2)+in%p(3)*in%p(3)))/EL_CHG  

  ! compute magnetic field intensity and direction
  B_norm = sqrt(B(1)*B(1)+B(2)*B(2)+B(3)*B(3)) !< intensity
  B_hat = B/B_norm  !< direction
  ! compute the parallel and perpendicular momenta
  p_perp = vector_cartesian_to_cylindrical(in%x(3),in%p)
  p_par = p_perp(1)*B_hat(1)+p_perp(2)*B_hat(2)+p_perp(3)*B_hat(3)
  p_perp = p_perp - p_par*B_hat

  ! compute the magnetic moment p_perp^2/(2*B) in [eV/T]
  ! the sign is given by the particle parallel momentum 
  out%mu = sign((ATOMIC_MASS_UNIT*(p_perp(1)*p_perp(1)+&
       p_perp(2)*p_perp(2)+p_perp(3)*p_perp(3))/&
       (2.d0*B_norm*mass*EL_CHG)),p_par)

  ! compute the gc position
  if(out%q.ne.0) then 
    call particle_position_to_gc(node_list,element_list,&
         in%x,in%st,in%i_elm,&
         vector_cartesian_to_cylindrical(in%x(3),in%p),in%q,&
         B_hat,B_norm,out%x,out%st,out%i_elm)
  endif  
end function relativistic_kinetic_to_gc

!---------------------------------------------------------------------------

!> This function transform a particle gc in a relativistic kinetic particle
!> data type. The gc phase space coordinates have to be provided
!> respectively in [eV] for the energy and [eV/T] for the magnetic moment.
!> inputs:
!>   node_list:    (type_node_list) node list
!>   element_list: (type_element_list) element list
!>   in:           (particle_gc) a guiding center particle
!>   chi:          (real8) gyro-angle in [rad]
!>   time:         (real8) particle time
!>   mass:         (real8) particle mass
!>   B:	           (real8)(3)(optional) magnetic field [T]
!> outputs:
!>   out: (particle_kinetic_relativistic) a kinetic relativistic particle
function gc_to_relativistic_kinetic(node_list,element_list,in,time,mass,chi,B) result(out)
  use data_structure
  use mod_coordinate_transforms, only: vector_cylindrical_to_cartesian
  use mod_pusher_tools, only: get_orthonormals
  use mod_pusher_tools, only: gc_position_to_particle
  ! declare input variables
  type(type_node_list),intent(in) :: node_list
  type(type_element_list),intent(in) :: element_list
  type(particle_gc),intent(in) :: in
  real(kind=8),intent(in) :: time,mass,chi !< time,particle mass [AMU], gyroangle [rad]
  real(kind=8),dimension(3),intent(in) :: B !< mass in AMU magnetic field in [T]
  ! declare output variables
  type(particle_kinetic_relativistic) :: out
  ! declare internal variables
  real(kind=8) :: B_norm, p_perp, p_par   !< magnetic intensity, perpendicular/parallel momenta
  real(kind=8),dimension(3) :: B_hat, e1, e2 !< B-field-aligned cartesian vector basis

  ! copy the default particle datatype
  out = in
  ! copy the particle charge
  out%q = in%q
  
  ! compute the magnetic field intensity and direction
  B_norm = norm2(B)
  B_hat = B/B_norm
  ! compute a B-field-aligned orthonormal vector basis
  call get_orthonormals(B_hat,e1,e2)

  ! compute the perpendicular momentum squared in (AMU*m/s)^2
  p_perp = (EL_CHG*2.d0*mass*B_norm*abs(in%mu))/ATOMIC_MASS_UNIT
  ! compute the parallel momentum in (AMU*m/s)
  p_par = sign(sqrt((((in%E*EL_CHG)*(in%E*EL_CHG))/ &
  ((ATOMIC_MASS_UNIT*SPEED_OF_LIGHT)*(ATOMIC_MASS_UNIT*SPEED_OF_LIGHT)))-&
  (mass*SPEED_OF_LIGHT)*(mass*SPEED_OF_LIGHT)-p_perp),in%mu)
  ! compute the perpendicular momentum in (AMU*m/s)
  p_perp = sqrt(p_perp)

  ! computing the particle momentum
  out%p = p_par*B_hat + p_perp*(e1*cos(chi)+e2*sin(chi)) 

  ! compute the particle position in R,Z,Phi coordinates
  if(out%q.ne.0) then
    call gc_position_to_particle(node_list,element_list,in%x,in%st,&
         in%i_elm,out%p,in%q,B_hat,B_norm,out%x,out%st,out%i_elm)
 endif
 !> transform the momentum in cartesian coodinates
 out%p = vector_cylindrical_to_cartesian(out%x(3),out%p)
end function gc_to_relativistic_kinetic



end module mod_kinetic_relativistic
