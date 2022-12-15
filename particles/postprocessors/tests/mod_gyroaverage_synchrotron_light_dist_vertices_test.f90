!> mod_gyroaverage_synchrotron_light_dist_test contains all variables and
!> procedures for testing the gyroaveraged synchrotron light vertices
module mod_gyroaverage_synchrotron_light_dist_vertices_test
use fruit
use mod_particle_types,                              only: particle_gc_relativistic_id
use mod_particle_types,                              only: particle_kinetic_relativistic_id
use mod_particle_types,                              only: particle_kinetic_id
use mod_particle_types,                              only: particle_gc_vpar_id
use mod_particle_sim,                                only: particle_sim
use mod_spectra_deterministic,                       only: spectrum_integrator_2nd
use mod_gyroaverage_synchrotron_light_dist_vertices, only: gyroaverage_synchrotron_light_dist
implicit none

private
public :: run_fruit_gyroaverage_synchrotron_light_dist_vertices

!> Variables ---------------------------------------------------------
!> general parameters
real*8,parameter :: tol_real8=1d-16
real*8,parameter :: tol2_real8=1d-16
real*8,parameter :: mass_RE=5.48579909065d-4
!> parameters for generating synchrotron lights
integer,parameter :: n_mhd_sol=16
integer,parameter :: n_properties=10
integer,parameter :: n_x=3
integer,parameter :: fill_type_base=1 !< use cylindrical initialisation
integer,parameter :: n_times_sol=2
integer,parameter :: n_particle_types_check_sol=1
integer,dimension(n_times_sol),parameter              :: n_groups_per_sim=(/4,3/)
integer,parameter                                     :: n_groups_max=maxval(n_groups_per_sim)
integer,dimension(n_groups_max,n_times_sol),parameter :: n_particles_per_group=&
           reshape((/132,324,42,10,237,143,23,34/),shape(n_particles_per_group))
integer,parameter :: n_particles_max=maxval(n_particles_per_group)
integer,dimension(n_groups_max,n_times_sol),parameter :: particle_types_sol=&
           reshape((/particle_kinetic_relativistic_id,particle_gc_relativistic_id,&
           particle_kinetic_id,particle_gc_vpar_id,particle_gc_relativistic_id,&
           particle_gc_vpar_id,particle_kinetic_id,0/),shape(particle_types_sol))
real*8,parameter                                      :: survival_threshold=0.45
!> parameters for generating spectra
integer,parameter                                     :: n_spectra=2
integer,parameter                                     :: n_lines_per_spectrum=16
real*8,dimension(n_spectra),parameter                 :: min_wlen=(/3d-6,2.5d-7/)
real*8,dimension(n_spectra),parameter                 :: max_wlen=(/3.5d-6,4.2d-7/)
!> parameters for generating shadowed points
integer,parameter                                     :: n_shaded_points=57
integer,parameter                                     :: n_shaded_points_per_particle=7
real*8,dimension(2),parameter                         :: length_shadowed=(/2d-2,7d1/) 
!> variables for generating synchrotron lights
type(gyroaverage_synchrotron_light_dist)              :: vertex_sol
type(particle_sim),dimension(n_times_sol)             :: sims_particles
integer                                               :: n_gc_RE_max
integer,dimension(n_particle_types_check_sol)         :: particle_types_check_sol
integer,dimension(n_times_sol)                        :: n_active_vertices_sol
integer,dimension(n_groups_max,n_times_sol)           :: n_active_particles_sol
integer,dimension(n_particles_max,n_groups_max,n_times_sol) :: active_particle_ids_sol
real*8,dimension(n_times_sol)                         :: time_vector_sol
real*8,dimension(n_x,n_particles_max*n_groups_max,n_times_sol)          :: x_cart_sol
real*8,dimension(n_properties,n_particles_max*n_groups_max,n_times_sol) :: properties_sol
!> variables for generating spectra
type(spectrum_integrator_2nd)                         :: spectrum
!> variables for generating shadowed points
real*8,dimension(:,:,:,:),allocatable                 :: x_shadowed

!> Interfaces --------------------------------------------------------
contains
!> Fruit test basket -------------------------------------------------
!> fruit basket having all set-up, tests and tearing-down procedures
subroutine run_fruit_gyroaverage_synchrotron_light_dist_vertices()
  implicit none
  write(*,'(/A)') "  ... setting-up: gyroaverage synchrotron light vertices tests"
  call setup
  write(*,'(/A)') "  ... running: gyroaverage synchrotron light vertices tests"
  write(*,'(/A)') "  ... tearing-down: gyroaverage synchrotron light vertices tests"
  call teardown
end subroutine run_fruit_gyroaverage_synchrotron_light_dist_vertices

!> Set-up and tear-down procedures------------------------------------
!> allocate and initialise the unit test features
subroutine setup()
  use mod_rng,                        only: type_rng
  use mod_pcg32_rng,                  only: pcg32_rng
  use mod_gnu_rng,                    only: gnu_rng_interval
  use mod_common_test_tools,          only: omp_initialize_rngs
  use mod_particle_common_test_tools, only: sim_time_interval
  use mod_particle_common_test_tools, only: allocate_one_particle_list_type
  use mod_particle_common_test_tools, only: fill_groups,fill_mass_RE 
  use mod_particle_common_test_tools, only: fill_particles_tokamak
  use mod_particle_common_test_tools, only: invalidate_particles
  use mod_particle_common_test_tools, only: obtain_active_particle_ids
  use mod_fields_linear,              only: jorek_fields_interp_linear
  !$ use omp_lib
  implicit none
  !> variables
  integer :: ii,jj,ifail,n_gc_RE_max_loc,n_threads
  class(type_rng),dimension(:),allocatable :: rngs
  !> initialisation
  particle_types_check_sol = (/particle_gc_relativistic_id/)
  vertex_sol%n_property_vertex = n_properties; ifail = 0;
  n_gc_RE_max = 0; n_active_particles_sol = 0; n_threads = 1;
  !$ n_threads = omp_get_max_threads()
  call gnu_rng_interval(n_times_sol,sim_time_interval,time_vector_sol)
  call vertex_sol%allocate_vertices(n_times_sol,n_particles_max*n_groups_max)
  !> allocate and initialise particle lists
  do ii=1,n_times_sol
    sims_particles(ii)%time = time_vector_sol(ii)
    allocate(jorek_fields_interp_linear::sims_particles(ii)%fields)
    allocate(sims_particles(ii)%groups(n_groups_per_sim(ii)))
    call allocate_one_particle_list_type(n_groups_per_sim(ii),&
    n_particles_per_group(1:n_groups_per_sim(ii),ii),&
    particle_types_sol(1:n_groups_per_sim(ii),ii),sims_particles(ii)%groups,ifail)
    call fill_groups(n_groups_per_sim(ii),sims_particles(ii)%groups)
    call fill_mass_RE(n_groups_per_sim(ii),sims_particles(ii)%groups)
    call fill_particles_tokamak(n_groups_per_sim(ii),sims_particles(ii)%groups,fill_type_base)
    call invalidate_particles(n_groups_per_sim(ii),n_particles_max,survival_threshold,&
    n_active_particles_sol(1:n_groups_per_sim(ii),ii),sims_particles(ii)%groups)
    call obtain_active_particle_ids(n_groups_per_sim(ii),n_particles_max,&
    active_particle_ids_sol(:,1:n_groups_per_sim(ii),ii),sims_particles(ii)%groups)
    n_active_vertices_sol(ii) = sum(n_active_particles_sol(:,ii))
    where(particle_types_sol(:,ii).ne.particle_gc_relativistic_id)
      n_active_particles_sol(:,ii) = 0
    endwhere
    n_gc_RE_max_loc = 0
    do jj=1,n_groups_per_sim(ii)
      if(particle_types_sol(jj,ii).ne.particle_gc_relativistic_id) cycle
      n_gc_RE_max_loc = n_gc_RE_max_loc + n_particles_per_group(jj,ii)
    enddo
      n_gc_RE_max = max(n_gc_RE_max,n_gc_RE_max_loc)
  enddo
  !> initialise positions and properties tables
  !TODO
  !> initialise deterministic spectra
  spectrum = spectrum_integrator_2nd(n_lines_per_spectrum,n_spectra,min_wlen,max_wlen)
  call spectrum%generate_spectrum()
  !> generate shadowed point positions
  allocate(x_shadowed(n_x,n_shaded_points_per_particle,n_gc_RE_max,n_times_sol))
  !TODO
end subroutine setup 

!> destroy all test features
subroutine teardown()
  implicit none
  call vertex_sol%deallocate_vertices; call spectrum%deallocate_spectrum;
  if(allocated(x_shadowed)) deallocate(x_shadowed)
end subroutine teardown
!> Tests -------------------------------------------------------------

!> Tools -------------------------------------------------------------

!>--------------------------------------------------------------------
end module mod_gyroaverage_synchrotron_light_dist_vertices_test
 
