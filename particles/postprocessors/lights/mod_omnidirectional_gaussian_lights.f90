!> the mod_omnidirectional_gaussian_lights implements variables
!> and procedure defining an omnidirectional light with a
!> gaussian spectrum
module mod_omnidirectional_gaussian_lights
use mod_light_vertices, only: light_vertices
implicit none

private
public :: omnidirectional_gaussian_lights

!> Variables ----------------------------------------------------
type,extends(light_vertices) :: omnidirectional_gaussian_lights
  contains
  procedure,pass(light_vert) :: init_lights_from_particles => &
                                init_omnidir_gaussian_lights_from_particles
  procedure,pass(light_vert) :: directionality_funct => &
                                omnidir_gaussian_directionality_funct
  procedure,pass(light_vert) :: spectral_irradiance => &
                                omnidir_gaussian_spectral_irradiance
  procedure,pass(light_vert),private :: fill_omnidir_gaussian_lights_from_particles
end type omnidirectional_gaussian_lights
!> Interfaces ---------------------------------------------------

contains

!> Procedures ---------------------------------------------------
!> init_omnidir_gaussian_lights_from_particles computes and stores
!> the properties of omnidirectional gaussian lights from particles
!> inputs:
!>   light_vert:     (omnidirectional_gaussian_light) empty
!>                   omnidirectional gaussian light
!>   n_times:        (integer) number of simulation times
!>   sims_particles: (particle_sim)(n_times) array of particle simulations
!>   n_lights_in:    (integer)(optional) number of requested lights
!> outputs:
!>   light_vert:     (omnidirectional_gaussian_light) initialised 
!>                   omnidirectional gaussian light
!>   sims_particles: (particle_sim)(n_times) array of particle simulations
subroutine init_omnidir_gaussian_lights_from_particles(light_vert,&
n_times,sims_particles,n_lights_in)
  use mod_particle_types, only: particle_kinetic_relativistic_id
  use mod_particle_sim,    only: particle_sim
  implicit none
  !> inputs-outputs:
  class(omnidirectional_gaussian_lights),intent(inout) :: light_vert
  type(particle_sim),dimension(n_times),intent(inout)  :: sims_particles
  !> inputs:
  integer,intent(in)                                   :: n_times
  integer,intent(in),optional                          :: n_lights_in
  !> variables:
  integer                              :: ii,n_groups_max,n_particles_max,n_lights
  integer,dimension(n_times)           :: n_groups,n_particles_relativistic
  integer,dimension(:,:),allocatable   :: n_particles,particle_types,n_active_particles
  integer,dimension(:,:,:),allocatable :: active_particle_ids

  !> initialisations
  light_vert%n_property_vertex = 2
  call light_vert%allocate_time_vector(n_times)
  call light_vert%fill_time_vector_particle_sims(sims_particles)
  !> allocate number of particles and particle types
  call light_vert%extract_n_groups_all_particle_sims(sims_particles,n_groups)
  n_groups_max = maxval(n_groups)
  allocate(n_particles(n_groups_max,light_vert%n_times))
  allocate(particle_types(n_groups_max,light_vert%n_times))
  call light_vert%extract_n_particles_all_particle_sims(sims_particles,&
  n_groups_max,n_particles)
  call light_vert%extract_particle_types_all_particle_sims(sims_particles,&
  n_groups_max,particle_types)
  !> compute the number of relativistic particles per each time
  do ii=1,light_vert%n_times
    n_particles_relativistic(ii) = sum(n_particles(:,ii),&
    mask=particle_types(:,ii).eq.particle_kinetic_relativistic_id)
  enddo
  n_particles_max = maxval(n_particles_relativistic)
  n_lights = n_particles_max
  !> allocate active particle arrays and vertices
  allocate(n_active_particles(n_groups_max,light_vert%n_vertices))
  allocate(active_particle_ids(n_particles_max,n_groups_max,light_vert%n_vertices))
  call light_vert%allocate_x_properties(n_lights)

  !> find active particles for all groups and times
  call light_vert%find_active_particles_id_time(n_groups_max,n_particles_max,&
  n_groups,n_particles,sims_particles,n_active_particles,active_particle_ids,&
  particle_kinetic_relativistic_id)
  !> fill the omnidirectional gaussian lights
  call light_vert%fill_omnidir_gaussian_lights_from_particles(sims_particles,&
  n_groups_max,n_particles_max,n_groups,n_active_particles,active_particle_ids)
  !> cleanup
  deallocate(n_particles); deallocate(particle_types);
  deallocate(n_active_particles); deallocate(active_particle_ids);
end subroutine init_omnidir_gaussian_lights_from_particles

!> omnidir_gaussian_spectral_irradiance computes the full spectral anguler
!> power distribution for omnidirectional gaussian lights
!> inputs:
!>   light_vert: (omnidirectional_gaussian_lights) omnidirectional gaussian lights
!>   spectra:    (spectrum_base) spectral intervals and integrators
!>   time_id:    (integer) the time index
!>   light_id:   (integer) the light index
!>   x_shaded:   (real8)(3) shaded point position in cartesian coord
!> outputs:
!>   light_vert:            (omnidirectional_gaussian_lights) omnidirectional gaussian lights
!>   spectra:               (spectrum_base) spectral intervals and integrators  
!>   light_spec_irradiance: (real8)(n_points,n_spectra) omnidirectional gaussian
!>                          full spectral angular distribution
subroutine omnidir_gaussian_spectral_irradiance(light_vert,spectra,time_id,light_id,&
x_shaded,light_spec_irradiance)
  use mod_spectra, only: spectrum_base
  implicit none
  !> inputs-outputs:
  class(omnidirectional_gaussian_lights),intent(inout) :: light_vert
  class(spectrum_base),intent(inout)                   :: spectra
  !> inputs:
  integer,intent(in)                                   :: time_id,light_id
  real*8,dimension(light_vert%n_x),intent(in)          :: x_shaded
  !> outputs:
  real*8,dimension(spectra%n_points,spectra%n_spectra),intent(out) :: light_spec_irradiance
  !> variables
  integer :: ii,jj
  real*8,dimension(spectra%n_spectra) :: spectra_midpoint
  !> initialisation
  spectra_midpoint = 5.d-1*(spectra%points(spectra%n_points,:)-spectra%points(spectra%n_points,:))
  !> compute irradiance
  !$omp parallel do default(shared) private(ii,jj) collapse(2)
  do ii=1,spectra%n_spectra
    do jj=1,spectra%n_points
      light_spec_irradiance(ii,jj) = exp(-((spectra%points(jj,ii)-spectra_midpoint(ii))*&
      (spectra%points(jj,ii)-spectra_midpoint(ii)))/(2.d0*light_vert%properties(1,light_id,time_id)))
    enddo
  enddo
  !$omp end parallel do
end subroutine omnidir_gaussian_spectral_irradiance

!> omnidir_gaussian_directionality_funct computes the directionality function
!> for omnidirectional gaussian lights
!> inputs:
!>   light_vert: (omnidirectional_gaussian_lights) omnidirectional gaussian lights
!>   spectra:    (spectrum_base) spectral intervals and integrators
!>   time_id:    (integer) the time index
!>   light_id:   (integer) the light index
!>   x_shaded:   (real8)(3) shaded point position in cartesian coord
!> outputs:
!>   light_vert: (omnidirectional_gaussian_lights) omnidirectional gaussian lights
!>   spectra:    (spectrum_base) spectral intervals and integrators  
!>   light_dstb: (real8)(n_points,n_spectra) omnidirectional gaussian light
!>               directionality function
subroutine omnidir_gaussian_directionality_funct(light_vert,spectra,time_id,light_id,&
x_shaded,light_dstb)
  use mod_spectra, only: spectrum_base
  implicit none
  !> inputs-outputs:
  class(omnidirectional_gaussian_lights),intent(inout) :: light_vert
  class(spectrum_base),intent(inout)                   :: spectra
  !> inputs:
  integer,intent(in)                                   :: time_id,light_id
  real*8,dimension(light_vert%n_x),intent(in)          :: x_shaded
  !> outputs:
  real*8,dimension(spectra%n_points,spectra%n_spectra),intent(out) :: light_dstb
  !> compute the directionality function
  call light_vert%spectral_irradiance(spectra,time_id,light_id,x_shaded,light_dstb)
  light_dstb = light_dstb/light_vert%properties(2,light_id,time_id)
end subroutine omnidir_gaussian_directionality_funct

!> Tools --------------------------------------------------------
!> fill_omnidir_gaussian_lights_from_particles fill the x and properties 
!> arrays of omnidirectional gaussian lights from particle lists
!> inputs:
!>   light_vert:          (omnidirectional_gaussian_lights) empty
!>                        omnidirectional gaussian lights
!>   sims_particles:      (particle_sim)(n_times) array of particle simulations
!>   n_groups_max:        (integer) maximum size of groups
!>   n_particles_max:     (integer)(n_times) maximum number of particles
!>   n_groups:            (integer)(n_times) size of each group
!>   n_active_particles:  (integer)(n_groups_max,n_times) number of active particles
!>                        per group and per time
!>   active_particle_ids: (integer)(n_particles_max,n_groups_max,n_times) indiced
!>                        of the active particles
!> outputs:
!>   light_vert: (omnidirectional_gaussian_lights) initialised
!>               omnidirectional gaussian lights
subroutine fill_omnidir_gaussian_lights_from_particles(light_vert,sims_particles,&
n_groups_max,n_particles_max,n_groups,n_active_particles,active_particle_ids)
  use constants,          only: PI,SPEED_OF_LIGHT
  use mod_particle_sim,   only: particle_sim
  use mod_particle_types, only: particle_kinetic_relativistic
  implicit none
  !> inputs-outputs:
  class(omnidirectional_gaussian_lights),intent(inout)          :: light_vert
  !> inputs:
  type(particle_sim),dimension(light_vert%n_times),intent(in)   :: sims_particles
  integer,intent(in)                                            :: n_groups_max
  integer,intent(in)                                            :: n_particles_max
  integer,dimension(light_vert%n_times),intent(in)              :: n_groups
  integer,dimension(n_groups_max,light_vert%n_times),intent(in) :: n_active_particles
  integer,dimension(n_particles_max,n_groups_max,light_vert%n_times),intent(in)::active_particle_ids
  real*8                                                        :: rel_factor
  real*8,dimension(light_vert%n_x)                              :: momentum
  !> variables:
  integer :: ii,jj,kk,pp
  !> compute omnidirectional light properties
  do ii=1,light_vert%n_times
    pp=0
    do jj=1,n_groups(ii)
      !$omp parallel default(private) firstprivate(ii,jj,pp,n_active_particles) &
      !$omp shared(sims_particles,active_particle_ids,light_vert)
      select type (p_list=>sims_particles(ii)%groups(jj)%particles)
        type is (particle_kinetic_relativistic)
        !$omp do
        do kk=1,n_active_particles(jj,ii)
          call light_vert%store_light_x_from_particle_id(pp+kk,ii,&
          p_list(active_particle_ids(kk,jj,ii))) !< store position
          !> compute properties
          momentum = p_list(active_particle_ids(kk,jj,ii))%p
          rel_factor = sqrt(1.d0+((momentum(1)*momentum(1) + momentum(2)*momentum(2) + &
                       momentum(3)*momentum(3))/(SPEED_OF_LIGHT*SPEED_OF_LIGHT*&
                       sims_particles(ii)%groups(jj)%mass*sims_particles(ii)%groups(jj)%mass)))
          light_vert%properties(:,pp+kk,ii) = (/rel_factor,1.d0/sqrt(2.d0*PI*rel_factor)/)
        enddo
        !$omp end do
      end select
      !$omp end parallel
      pp = pp + n_active_particles(jj,ii)
    enddo
  enddo
end subroutine fill_omnidir_gaussian_lights_from_particles

!>---------------------------------------------------------------
end module mod_omnidirectional_gaussian_lights

