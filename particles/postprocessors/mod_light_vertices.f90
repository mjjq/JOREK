!> The mod_light_vertices module contains variables
!> and procedures used for defining and defining actions
!> of the light points
module mod_light_vertices
use mod_vertices, only: vertices
implicit none

private
public :: light_vertices

!> Variables --------------------------------------------
type,abstract,extends(vertices) :: light_vertices
  contains
  procedure,pass(light_vert)                              :: fill_time_vector_particle_sims
  procedure,pass(light_vert)                              :: extract_n_groups_all_particle_sims 
  procedure,pass(light_vert)                              :: extract_n_particles_all_particle_sims
  procedure,pass(light_vert)                              :: extract_particle_types_all_particle_sims
  procedure,pass(light_vert)                              :: store_light_x_from_particle_id
  procedure,pass(light_vert)                              :: find_active_particles_id_time
  procedure(init_lights_parts),pass(light_vert),deferred  :: init_lights_from_particles
  procedure(direct_funct),pass(light_vert),deferred       :: directionality_funct
  procedure(spect_irradiance),pass(light_vert),deferred   :: spectral_irradiance
end type light_vertices

!> Interfaces -------------------------------------------
!> module procedure
!interface fill_time_vector
!  module procedure fill_time_vector_particle_sims
!end interface fill_time_vector

interface store_x_from_id
  module procedure store_light_x_from_particle_id
end interface store_x_from_id

interface
  !> computes and store the coordinates and properties of 
  !> lights from particle simulations
  !> inputs:
  !>   light_vert:       (light_vertices) empty light vertices
  !>   n_times:          (integer) number of times
  !>   sims_particles:   (particle_sim)(n_times) array of particle simulations
  !>   n_sync_lights_in: (integer)(optional) number of requested synchrotron lights
  !> outputs:
  !>   light_vert:     (light_vertices) filled light vertices
  !>   sims_particles: (particle_sim)(n_times) array of particle simulations
  subroutine init_lights_parts(light_vert,n_times,sims_particles,n_sync_lights_in)
    use mod_particle_sim, only: particle_sim
    IMPORT :: light_vertices
    implicit none
    !> inputs-outputs
    class(light_vertices),intent(inout) :: light_vert
    type(particle_sim),dimension(n_times),intent(inout) :: sims_particles
    !> inputs
    integer,intent(in) :: n_times
    integer,intent(in),optional :: n_sync_lights_in
  end subroutine init_lights_parts

  !> computes the directionality function for a given point
  !> in space (cartesian coordinate)-time and a given light
  !> for all spectra wavelengths
  !> inputs:
  !>   light_vert: (light_vertices) initialised light vertices
  !>   spectra:    (spectrum_base) initilises spectra
  !>   light_id:   (integer) id of the light to use
  !>   x_shaded:   (real8)(n_x) point illuminated by the light
  !> outputs:
  !>   light_vert: (light_vertices) initialised light vertices
  !>   spectra:    (spectrum_base) initilises spectra
  !>   light_dstb: (real8)(n_points,n_spectra) light intensity distribution
  !>               from the light_id light to the x_shaded point for all
  !>               spectra points and all spectra
  subroutine direct_funct(light_vert,spectra,time_id,light_id,x_shaded,light_dstb)
    use mod_spectra,  only: spectrum_base
    IMPORT :: light_vertices
    implicit none
    !> inputs-outpus
    class(light_vertices),intent(inout) :: light_vert
    class(spectrum_base),intent(inout)  :: spectra
    !> inputs
    integer,intent(in)                  :: time_id,light_id
    real*8,dimension(light_vert%n_x),intent(in)    :: x_shaded
    !> outputs
    real*8,dimension(spectra%n_points,spectra%n_spectra),intent(out) :: light_dstb
  end subroutine direct_funct

  !> compute the spectral irradiance of a light source for a given point
  !> in spae (cartesian coordinates)-time and a given light source
  !> for all wavelengths and spectra 
  !> inputs:
  !>   light_vert: (light_vertices) initialised light vertices
  !>   spectra:    (spectrum_base) initilises spectra
  !>   light_id:   (integer) id of the light to use
  !>   x_shaded:   (real8)(n_x) point illuminated by the light
  !> outputs:
  !>   light_vert: (light_vertices) initialised light vertices
  !>   spectra:    (spectrum_base) initilises spectra
  !>   light_spec_irradiance: (real8)(n_points,n_spectra) spectral irradiance
  !>                          from the light_id light at ethe time time_id to 
  !>                          the x_shaded point for all spectra points and all spectra
  subroutine spect_irradiance(light_vert,spectra,time_id,light_id,x_shaded,light_spec_irradiance)
    use mod_spectra,  only: spectrum_base
    IMPORT :: light_vertices
    implicit none
    !> inputs-outpus
    class(light_vertices),intent(inout) :: light_vert
    class(spectrum_base),intent(inout)  :: spectra
    !> inputs
    integer,intent(in)                  :: time_id,light_id
    real*8,dimension(light_vert%n_x),intent(in)    :: x_shaded
    !> outputs
    real*8,dimension(spectra%n_points,spectra%n_spectra),intent(out) :: light_spec_irradiance
  end subroutine spect_irradiance
end interface

contains

!> Procedures -------------------------------------------
!> find active particles for all particle lists, particle
!> groups and simulation times
!> inputs:
!>   light_vert: (light_vertices) initialised light vertices
!>   n_groups_max: (integer) maximum number of groups
!>   n_particles_max: (integer) maximum number of particles
!>   n_groups:        (integer)(n_times) number of groups per time
!>   n_particles:     (integer)(n_groups,n_times) number of particles
!>                    per group per time
!>   sims_particles:  (particle_sim)(n_times) array of particle simulations
!> outputs:
!>   light_vert:         (light_vertices) initialised light vertices
!>   sims_particles:     (particle_sim)(n_times) array of particle simulations
!>   n_active_particles: (integer)(n_groups,n_times) number of active
!>                       particles for each group and time
!>   active_particle_id: (integer)(n_particles_max,n_groups_max,n_times)
!>                       particle list index of active particles
subroutine find_active_particles_id_time(light_vert,n_groups_max,n_particles_max,&
n_groups,n_particles,sims_particles,n_active_particles,active_particle_id,p_type)
  use mod_particle_sim,only: particle_sim
  implicit none
  !> inputs-outputs
  class(light_vertices),intent(inout) :: light_vert
  type(particle_sim),dimension(light_vert%n_times),intent(inout) :: sims_particles
  !> inputs
  integer,intent(in) :: n_groups_max,n_particles_max
  integer,dimension(light_vert%n_times),intent(in) :: n_groups
  integer,dimension(n_groups_max,light_vert%n_times),intent(in) :: n_particles
  integer,intent(in),optional :: p_type
  !> outputs
  integer,dimension(n_groups_max,light_vert%n_times),intent(out) :: n_active_particles
  integer,dimension(n_particles_max,n_groups_max,light_vert%n_times),intent(out)::active_particle_id
  !> variables
  integer :: ii
  n_active_particles = 0; active_particle_id = 0;
  if(present(p_type)) then
    do ii=1,light_vert%n_times
      call sims_particles(ii)%find_active_particles_groups(n_groups(ii),n_particles_max,&
      n_particles(:,ii),n_active_particles(:,ii),active_particle_id(:,:,ii),p_type)
      light_vert%n_active_vertices(ii) = sum(n_active_particles(:,ii))
    enddo
  else
    do ii=1,light_vert%n_times
      call sims_particles(ii)%find_active_particles_groups(n_groups(ii),n_particles_max,&
      n_particles(:,ii),n_active_particles(:,ii),active_particle_id(:,:,ii))
      light_vert%n_active_vertices(ii) = sum(n_active_particles(:,ii))    
    enddo
  endif
end subroutine find_active_particles_id_time

!> fill the time vector from particle simulations
!> inputs:
!>   light_vert:     (light_vertices) light vertices type
!>   sims_particles: (particle_sim)(n_times) array of particle simulations
!> outputs:
!>   light_vert:     (light_vertices) light vertices type
!>   sims_particles: (particle_sim)(n_times) array of particle simulations
subroutine fill_time_vector_particle_sims(light_vert,sims_particles)
  use mod_particle_sim, only: particle_sim
  implicit none
  !> inputs-outputs
  class(light_vertices),intent(inout)                            :: light_vert
  type(particle_sim),dimension(light_vert%n_times),intent(inout) :: sims_particles
  !> variables
  integer :: ii
  do ii=1,light_vert%n_times
    light_vert%times(ii) = sims_particles(ii)%time
  enddo
end subroutine fill_time_vector_particle_sims

!> Extract the number of groups for all times (simulations)
!> inputs:
!>   light_vert:     (light_vertices) light vertices type
!>   sims_particles: (particle_sim)(n_times) array of particle simulations
!> outputs:
!>   light_vert:     (light_vertices) light vertices type
!>   sims_particles: (particle_sim)(n_times) array of particle simulations
!>   n_groups:       (integer)(n_times) number of groups
subroutine extract_n_groups_all_particle_sims(light_vert,sims_particles,n_groups)
  use mod_particle_sim, only: particle_sim
  implicit none
  !> inputs-outpus:
  class(light_vertices),intent(inout)                            :: light_vert
  type(particle_sim),dimension(light_vert%n_times),intent(inout) :: sims_particles
  !> outputs:
  integer,dimension(light_vert%n_times),intent(out) :: n_groups
  !> variables
  integer :: ii
  !> extract number of groups for all simulations
  n_groups = 0
  do ii=1,light_vert%n_times
    n_groups(ii) = sims_particles(ii)%compute_group_size()
  enddo
end subroutine extract_n_groups_all_particle_sims

!> extract the number of particles for all groups and times (simulations)
!> inputs:
!>   light_vert:     (light_vertices) light vertices type
!>   sims_particles: (particle_sim)(n_times) array of particle simulations
!>   n_groups_max:   (integer) maximum number of among all simulations
!>   n_groups:       (integer)(n_times) number of groups per simulation
!> outputs:
!>   light_vert:     (light_vertices) light vertices type
!>   sims_particles: (particle_sim)(n_times) array of particle simulations
!>   n_particles:    (integer)(n_groups_max,n_times) number of particles
!>                   per group and per simulation
subroutine extract_n_particles_all_particle_sims(light_vert,sims_particles,&
n_groups_max,n_particles)
  use mod_particle_sim, only: particle_sim
  implicit none
  !> inputs-ouputs:
  class(light_vertices),intent(inout)                            :: light_vert
  type(particle_sim),dimension(light_vert%n_times),intent(inout) :: sims_particles
  !> inputs:
  integer,intent(in)                                             :: n_groups_max
  !> outouts:
  integer,dimension(n_groups_max,light_vert%n_times),intent(out) :: n_particles
  !> variables
  integer :: ii,n_groups
  !> extract number of particles
  n_particles = 0; n_groups = n_groups_max;
  do ii=1,light_vert%n_times
    call sims_particles(ii)%compute_particle_sizes(n_groups,n_particles(:,ii))
  enddo
end subroutine extract_n_particles_all_particle_sims

!> extract the particle types for all groups and simulations
!> inputs:
!>   light_vert:     (light_vertices) light vertices type
!>   sims_particles: (particle_sim)(n_times) array of particle simulations
!>   n_groups_max:   (integer) maximum number of among all simulations
!> outputs:
!>   light_vert:     (light_vertices) light vertices type
!>   sims_particles: (particle_sim)(n_times) array of particle simulations
!>   particle_types: (integer)(n_groups_max,n_times) particle types
!>                   per group and per simulation. CODEX:
!>                   1 -> particle_fieldline
!>                   2 -> particle_gc
!>                   3 -> particle_gc_vpar
!>                   4 -> particle_gc_Qin
!>                   5 -> particle_kinetic
!>                   6 -> particle_kinetic_leapfrog
!>                   7 -> particle_kinetic_relativistic
!>                   8 -> particle_gc_relativistic
subroutine extract_particle_types_all_particle_sims(light_vert,sims_particles,&
n_groups_max,particle_types)
  use mod_particle_sim, only: particle_sim
  implicit none
  !> inputs-outputs:
  class(light_vertices),intent(inout)                            :: light_vert
  type(particle_sim),dimension(light_vert%n_times),intent(inout) :: sims_particles
  !> inputs:
  integer,intent(in)                                          :: n_groups_max
  !> outputs:
  integer,dimension(n_groups_max,light_vert%n_times),intent(out) :: particle_types
  !> variables
  integer :: ii,n_groups
  particle_types = 0; n_groups = n_groups_max;
  do ii=1,light_vert%n_times
    call sims_particles(ii)%find_particle_types(n_groups,particle_types(:,ii))
  enddo
end subroutine extract_particle_types_all_particle_sims

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
  use mod_coordinate_transforms, only: cylindrical_to_cartesian
  implicit none
  !> inputs-outputs
  class(light_vertices),intent(inout) :: light_vert
  !> inputs
  class(particle_base),intent(in) :: particle
  integer,intent(in)              :: light_id,time_id
  light_vert%x(:,light_id,time_id) = cylindrical_to_cartesian(particle%x)
end subroutine store_light_x_from_particle_id

!> Tools ------------------------------------------------
!> compute particle array size per requested type

end module mod_light_vertices
