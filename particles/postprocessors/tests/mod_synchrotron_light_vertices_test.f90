!> mod_synchrotron_light_vertices_test contains all variables and
!> procedures for testing the synchrotron light vertices
module mod_synchrotron_light_vertices_test
use fruit
use mod_particle_sim,               only: particle_sim
use mod_vertices,                   only: n_x
use mod_synchrotron_light_vertices, only: n_properties
use mod_synchrotron_light_vertices, only: synchrotron_light_vertices
implicit none

private
public :: run_fruit_synchrotron_light_vertices

!> Variables ---------------------------------------------------------
integer,parameter :: n_times_sol=2
integer,dimension(n_times_sol),parameter :: n_groups_per_sim=(/3,2/)
integer,parameter                        :: n_groups_max=maxval(n_groups_per_sim)
integer,dimension(n_groups_max,n_times_sol),parameter :: n_particles_per_group=&
           reshape((/135,247,512,367,413,0/),shape(n_particles_per_group))
integer,parameter                        :: n_particles_max=maxval(n_particles_per_group)
real*8,parameter                         :: survival_threshold=0.33
real*8,parameter                         :: tol_real8=5.d-16
type(synchrotron_light_vertices)            :: vertex_sol
type(particle_sim),dimension(n_times_sol)   :: sims_particles
integer,dimension(n_times_sol)              :: n_active_vertices_sol
integer,dimension(n_groups_max,n_times_sol) :: n_active_particles_sol
integer,dimension(n_particles_max,n_groups_max,n_times_sol) :: active_particle_ids_sol
real*8,dimension(n_times_sol)               :: time_vector_sol
real*8,dimension(n_x,n_particles_max*n_groups_max,n_times_sol) :: x_cart_sol
real*8,dimension(n_properties,n_particles_max*n_groups_max,n_times_sol) :: propeties_sol

!> Interfaces --------------------------------------------------------
contains
!> Fruit test basket -------------------------------------------------
!> fruit basket having all set-up, tests and tearing-down procedures
subroutine run_fruit_synchrotron_light_vertices()
  implicit none
  write(*,'(/A)') "  ... setting-up: synchrotron light vertices tests"
  call setup
  write(*,'(/A)') "  ... running: synchrotron light vertices tests"
  write(*,'(/A)') "  ... tearing-down: synchrotron light vertices tests"
end subroutine run_fruit_synchrotron_light_vertices

!> Set-up and tear-down procedures------------------------------------
!> allocate and initialise the unit test features
subroutine setup()
  use mod_particle_types,                   only: particle_gc_vpar_id
  use mod_particle_types,                   only: particle_kinetic_id
  use mod_particle_types,                   only: particle_kinetic_relativistic_id
  use mod_gnu_rng,                          only: gnu_rng_interval
  use mod_particle_common_test_tools,       only: sim_time_interval
  use mod_particle_common_test_tools,       only: allocate_one_particle_list_type 
  use mod_particle_common_test_tools,       only: fill_particles
  use mod_particle_common_test_tools,       only: invalidate_particles
  use mod_particle_common_test_tools,       only: obtain_active_particle_ids
  use mod_light_vertices_common_test_tools, only: compute_x_cart_particles 
  implicit none
  !> variables
  integer,dimension(n_groups_max,n_times_sol),parameter :: particle_types=&
  reshape((/particle_kinetic_relativistic_id,particle_gc_vpar_id,&
  particle_kinetic_relativistic_id,particle_kinetic_id,&
  particle_kinetic_relativistic_id,-1/),shape(particle_types))
  integer :: ii,ifail
  !> initialisation
  ifail = 0; n_active_particles_sol = 0;
  call gnu_rng_interval(n_times_sol,sim_time_interval,time_vector_sol)
  call vertex_sol%allocate_vertices(n_times_sol,n_particles_max*n_groups_max)

  !> allocate and initialise particle list
  do ii=1,n_times_sol
    sims_particles(ii)%time = time_vector_sol(ii)
    allocate(sims_particles(ii)%groups(n_groups_per_sim(ii)))
    call allocate_one_particle_list_type(n_groups_per_sim(ii),&
    n_particles_per_group(1:n_groups_per_sim(ii),ii),&
    particle_types(1:n_groups_per_sim(ii),ii),sims_particles(ii)%groups,ifail)
    call fill_particles(n_groups_per_sim(ii),sims_particles(ii)%groups)
    call invalidate_particles(n_groups_per_sim(ii),n_particles_max,survival_threshold,&
    n_active_particles_sol(1:n_groups_per_sim(ii),ii),sims_particles(ii)%groups)
    call obtain_active_particle_ids(n_groups_per_sim(ii),n_particles_max,&
    active_particle_ids_sol(:,1:n_groups_per_sim(ii),ii),sims_particles(ii)%groups)
    n_active_vertices_sol(ii) = sum(n_active_particles_sol(:,ii))
  enddo

  !> initialise positions and properties tables
  call compute_x_cart_particles(n_times_sol,n_groups_max,n_particles_max,&
  sims_particles,x_cart_sol)
end subroutine setup

!> Tests -------------------------------------------------------------
!> Tools -------------------------------------------------------------
!>--------------------------------------------------------------------
end module mod_synchrotron_light_vertices_test
