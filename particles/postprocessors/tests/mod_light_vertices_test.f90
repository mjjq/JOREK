!> the mod_light_vertices_test module contains variables
!> and methods for testing all procedures contained in
!> mod_light_vertices.f90. Due to the fact that light
!> vertices is an abstract class, the synchrotron vertex
!> type is used instead.
module mod_light_vertices_test
use fruit
use mod_vertices, only: n_x
use mod_particle_sim, only: particle_sim
use mod_synchrotron_light_vertices, only: synchrotron_light_vertices
use mod_particle_types, only: particle_fieldline_id,particle_gc_id
use mod_particle_types, only: particle_gc_vpar_id,particle_gc_Qin_id
use mod_particle_types, only: particle_kinetic_id,particle_kinetic_leapfrog_id
use mod_particle_types, only: particle_kinetic_relativistic_id
use mod_particle_types, only: particle_gc_relativistic_id
implicit none

private
public :: run_fruit_light_vertices

!> Variables ---------------------------------------------------------
integer,parameter                         :: n_properties_sol=11
integer,parameter                         :: n_times_sol=3
integer,dimension(n_times_sol),parameter  :: n_groups_per_sim=(/3,1,2/)
integer,parameter                         :: n_groups_max=maxval(n_groups_per_sim)
integer,dimension(n_groups_max,n_times_sol),parameter :: n_particles_per_group=&
        reshape((/55,30,44,27,0,0,37,51,0/),shape(n_particles_per_group))
integer,dimension(n_groups_max,n_times_sol),parameter :: particle_types_sol=&
        reshape((/particle_gc_vpar_id,particle_kinetic_relativistic_id,&
        particle_kinetic_relativistic_id,particle_kinetic_relativistic_id,0,0,&
        particle_kinetic_leapfrog_id,particle_kinetic_relativistic_id,0/),&
        shape(particle_types_sol))
real*8,parameter                            :: survival_threshold=0.21
type(synchrotron_light_vertices)            :: vertex_sol
type(particle_sim),dimension(n_times_sol)   :: sims_particles
integer                                     :: n_particles_max
integer,dimension(n_groups_max,n_times_sol) :: n_active_particles_sol
integer,dimension(:,:,:),allocatable        :: active_particle_ids_sol
real*8,dimension(n_times_sol)               :: time_vector_sol
real*8,dimension(:,:,:),allocatable         :: x_cart_sol

!> Interfaces --------------------------------------------------------
contains
!> Fruit basket ------------------------------------------------------
!> fruit basket containing all set-up, test and tear-down methods
subroutine run_fruit_light_vertices()
  implicit none
  write(*,*) "  ... setting-up: light vertices tests"
  call setup()
  write(*,*) "  ... running: light vertices tests"
  write(*,*) "  ... tearing-down: light vertices tests"
  call teardown()
end subroutine run_fruit_light_vertices

!> Set-up and teard-down ---------------------------------------------
!> set-up features common to all unit test
subroutine setup()
  use mod_gnu_rng,                          only: gnu_rng_interval
  use mod_particle_common_test_tools,       only: sim_time_interval
  use mod_particle_common_test_tools,       only: allocate_one_particle_list_type
  use mod_particle_common_test_tools,       only: fill_particles
  use mod_particle_common_test_tools,       only: invalidate_particles
  use mod_particle_common_test_tools,       only: obtain_active_particle_ids
  use mod_light_vertices_common_test_tools, only: compute_store_x_cart_particles
  implicit none
  integer :: ii,ifail
  real*8,dimension(n_times_sol) :: time_vector_sol
  !> initialisation
  ifail = 0; n_particles_max = 0; n_particles_max = maxval(n_particles_per_group); 
  n_active_particles_sol = 0
  call gnu_rng_interval(n_times_sol,sim_time_interval,time_vector_sol)

  !> allocate and initialise particle lists
  allocate(active_particle_ids_sol(1:n_particles_max,n_groups_max,n_times_sol))
  active_particle_ids_sol = 0;
  do ii=1,n_times_sol
    sims_particles(ii)%time = time_vector_sol(ii)
    allocate(sims_particles(ii)%groups(n_groups_per_sim(ii)))
    call allocate_one_particle_list_type(n_groups_per_sim(ii),&
    n_particles_per_group(1:n_groups_per_sim(ii),ii),&
    particle_types_sol(1:n_groups_per_sim(ii),ii),sims_particles(ii)%groups,ifail)
    call fill_particles(n_groups_per_sim(ii),sims_particles(ii)%groups)
    call invalidate_particles(n_groups_per_sim(ii),n_particles_max,survival_threshold,&
    n_active_particles_sol(1:n_groups_per_sim(ii),ii),sims_particles(ii)%groups)
    call obtain_active_particle_ids(n_groups_per_sim(ii),n_particles_max,&
    active_particle_ids_sol(:,1:n_groups_per_sim(ii),ii),sims_particles(ii)%groups)
  enddo

  !> initialise and allocate vertices
  vertex_sol%n_property_vertex = n_properties_sol; 
  call vertex_sol%allocate_vertices(n_times_sol,n_particles_max*n_groups_max)
  !> compute and store the cartesian particle position for all particles
  allocate(x_cart_sol(n_x,n_groups_max*n_particles_max,n_times_sol))
  call compute_store_x_cart_particles(n_times_sol,n_groups_max,n_particles_max,&
  sims_particles,x_cart_sol)
end subroutine setup

subroutine teardown()
  implicit none
  vertex_sol%n_property_vertex=0; deallocate(active_particle_ids_sol); 
  deallocate(x_cart_sol)
end subroutine teardown

!> Tests -------------------------------------------------------------

!> Tools -------------------------------------------------------------
  
!>--------------------------------------------------------------------
end module mod_light_vertices_test
