!> The mod_light_vertices module contains variables
!> and procedures used for defining and defining actions
!> of the light points
module mod_light_vertices
use mod_vertices, only: vertices
implcit none

private
public :: light_vertices

!> Variables --------------------------------------------
type,abstract,extends(vertices) :: light_vertices
  contains
  procedure,pass(light_vertices)                              :: store_light_x_from_particle_id
  procedure(init_lights_parts),deferred,pass(light_vertices)  :: init_lights_from_particles
  procedure(direct_funct),deferred,pass(light_vertices)       :: directionality_funct
end type light_vertices

!> Interfaces -------------------------------------------
interface
  !> computes and store the coordinates and properties of 
  !> lights from particle simulations
  !> inputs:
  !>   light_vert:   (light_vertices) empty light vertices
  !>   sim_particle: (particle_sim) initialised particle simulation
  !> outputs:
  !>   light_vert: (light_vertices) filled light vertices
  subroutine init_lights_parts(light_vert,sim_particle)
    use mod_particle_sim, only: particle_sim
    implicit none
    !> inputs-outputs
    class(light_vertices),intent(inout) :: light_vert
    !> inputs
    class(particle_sim),intent(in) :: sim_particle
  end subroutine init_lights_parts

  !> computes the directionality function for a given point
  !> in space (cartesian coordinate) and a given light for
  !> all spectra wavelengths
  !> inputs:
  !>   light_vert: (light_vertices) initialised light vertices
  !>   spectra:    (spectrum_base) initilises spectra
  !>   light_id:   (integer) id of the light to use
  !>   x_shaded:   (real8)(n_x) point illuminated by the light
  !> outputs:
  !>   light_vert: (light_vertices) initialised light vertices
  !>   light_dstb: (real8)(n_points,n_spectra) light intensity distribution
  !>               from the light_id light to the x_shaded point for all
  !>               spectra points and all spectra
  subroutine direct_funct(light_vert,spectra,light_id,x_shaded,light_dstb)
    use mod_spectra,      only: spectrum_base
    use mod_particle_sim, only: particle_sim
    implicit none
    !> inputs-outpus
    class(light_vertices),intent(inout) :: light_vert
    !> inputs
    class(spectrym_base),intent(in) :: spectra
    integer,intent(in) :: light_id
    real*8,dimension(light_vert%n_x),intent(in) :: x_shaded
    !> outputs
    real*8,dimension(spectra%n_points,spectra%n_spectra),intent(out) :: light_dstb
  end subroutine direct_funct
end interface

contains

!> Procedures -------------------------------------------
!> store the light position in cartesian coordinate give a particle,
!> the light and time index
!> inputs:
!>   light_vert: (light_vertices) light vertices structure
!>   light_id:   (integer) index of the light in the x table
!>   time_id:    (integer) index of the time in the x table
!>   particle:   (particle_base) particle structure
!> outputs:
!>   light_vert: (light_vertices) light vertices with new x entry
subroutine store_light_x_from_particle_id(light_vert,light_id,time_id,particle)
  use mod_particle_types,       only: particle_base
  use mod_coordinate_transform, only: cylindrical_to_cartesian
  implicit none
  !> inputs-outputs
  class(light_vertices),intent(inout) :: light_vertices
  !> inputs
  class(particle_base),intent(in) :: particle
  integer,intent(in)              :: light_id,time_id
  light_vert%x(:,light_id,time_id) = cylindrical_to_cartesian(particle%x)
end subroutine store_light_x_from_particle_id

!> Tools ------------------------------------------------
!>-------------------------------------------------------
end module mod_light_vertices
