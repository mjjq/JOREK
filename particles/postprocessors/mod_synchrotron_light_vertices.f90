!> the mod_synchrotron_light implements variable and
!> procedures defining a synchrotron light source
module mod_synchrotron_light_vertices
use mod_light_vertices, only: light_vertices
implicit none

private
public :: snchrt_light_vertices

!> Variables ---------------------------------------
type,extends(light_vertices) :: snchrt_light_vertices
  contains
  procedure,pass(snchrt_light_vertices) :: init_lights_from_particles =>&
                                           init_snchrt_light_from_particles
  procedure,pass(snchrt_light_vertices) :: directionality_funct => &
                                           snchrt_directionality_funct
end type snchrt_light_vertices
!> Interfaces --------------------------------------

contains

!> Procedures --------------------------------------

!> Tools ------------------------------------------
!> cmpt_snchrt_light_properties computes the
!> synchrotron radiation properties from a
!> kinetic relativistic particle. Variables
!> inputs:
!>   field_size:    (integer) size of the field vectors
!>   property_size: (integer) size of the property vector
!>   particle_in:   (particle_kinetic_relativistic) jorek particle
!>   mass:          (real8) mass of the particle
!>   E_field_cart:  (real8)(field_size) electric field at particle position
!>                  in the cartesian reference system
!>   B_field_cart:  (real8)(field_size) magnetic field at particle position
!>                  in the cartesian reference system
!> outputs:
!>   snchrt_properties: (real8)(property_size) synchrotron radiation properties
!>                      1:3 -> component of the velocity direction (cartesian)
!>                      4   -> beta -> velocity/speed of light
!>                      5   -> relativistic factor
!>                      6   -> orbit curvature (L. Carbakal, PPCF, 2017)
!>                      7   -> total radiation power (L. Carbajal, PPCF, 2017)
subroutine cmpt_snchrt_light_properties(field_size,property_size,particle_in,&
,mass,E_field_cart,B_field_cart,snchrt_properties)
  use constants,          only: PI,EPS_ZERO,EL_CHG,ATOMIC_MASS_UNIT,SPEED_OF_LIGHT
  use mod_math_operators, only: cross_product
  use mod_particle_types, only: particle_kinetic_relativistic
  implicit none
  !> inputs
  type(particle_kinetic_relativistic),intent(in) :: particle_in
  integer                                        :: field_size,property_size
  real*8,intent(in)                              :: mass
  real*8,dimension(field_size),intent(in)        :: B_field,E_field
  !> outputs
  real*8,dimension(property_size),intent(out) :: snchrt_properties
  !> variables
  real*8 :: velocity
  real*8,dimension(3) :: vector_1d_3

  !> compute velocity, velocity direction and relativistic factor
  velocity =  sqrt(particle_in%x(1)*particle_in%x(1)+&
              particle_in%x(2)*particle_in%x(2)+&
              particle_in%x(3)*particle_in%x(3)) 
  snchrt_properties(1:3) = particle_in%p/velocity
  snchrt_properties(4)   = velocity/SPEEF_OF_LIGHT
  snchrt_properties(5)   = sqrt(1.d0 + (snchrt_properties(4)*snchrt_properties(4))/(mass*mass))
  snchrt_properties(4)   = property(4)/(mass*snchrt_properties(5))
  !> compute orbit curvature
  vector_1d_3 = cross_product(snchrt_properties(1:3),E_field+&
                cross_product(particle_in%p/(mass*snchrt_properties(5)),B_field))
  snchrt_properties(6) = (abs(real(particle_in%q,kind=8))*EL_CHG*&
                         sqrt(vector_1d_3(1)*vector_1d_3(1)+&
                         vector_1d_3(2)*vector_1d_3(2)+&
                         vector_1d_3(3)*vector_1d_3(3)))/&
                         (snchrt_properties(5)*mass*ATOMI_MASS_UNIT*&
                         snchrt_properties(4)*snchrt_properties(4)*&
                         SPEED_OF_LIGHT*SPEED_OF_LIGHT)
  !> compute total synchrotron power
  snchrt_properties(7) = (real(particle_in%q*particle_in%q,kind=8)*&
                         EL_CHG*EL_CHG*SPEED_OF_LIGHT*snchrt_properties(4)*&
                         snchrt_properties(4)*snchrt_properties(4)*&
                         snchrt_properties(4)*snchrt_properties(5)*&
                         snchrt_properties(5)*snchrt_properties(5)*&
                         snchrt_properties(5)*snchrt_properties(6)*&
                         snchrt_properties(6))/(6.d0*PI*EPS_ZERO)
end subroutine cmpt_snchrt_light_properties

!>-------------------------------------------------
end module mod_synchrotron_light_vertices
