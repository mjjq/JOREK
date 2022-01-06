!> mod_particle_sim_test contains variables and 
!> procedure for testing the methods coded in
!> mod_particle_sim which do not require MPI
module mod_particle_sim_test
use fruit
use mod_particle_sim, only: particle_sim
use mod_particle_common_test_tools, only: n_particle_types
implicit none

private
public :: run_fruit_particle_sim

!> Variables ------------------------------------
integer,parameter :: n_particles=1000
real*8,parameter  :: survival_prob=6.25d-1
real*8,parameter  :: t_test=2.345d-6
integer,dimension(n_particle_types) :: n_active_particles_sol
integer,dimension(n_particle_types) :: particle_type_list_sol
integer*1,dimension(n_particles,n_particle_types) :: particle_charge_list_sol
integer,dimension(n_particles,n_particle_types)   :: active_particle_ids_sol
type(particle_sim)  :: sim_sol
!> Interfaces -----------------------------------

contains

!> Fruit basket ---------------------------------
!> fruit basket containing all set-up, test and
!> tear-down procedures
subroutine run_fruit_particle_sim()
  implicit none
  write(*,'(/A)') "  ... setting-up: particle sim tests"
  call setup
  write(*,'(/A)') "  ... running: particle sim tests"
  call test_codify_particle_sim
  call test_find_active_particle_id
  call test_find_active_particle_id_type
  write(*,'(/A)') "  ... tearing-up: particle sim tests"
end subroutine run_fruit_particle_sim

!> Set-up and tear-down -------------------------
!> set-up particle types test features
subroutine setup()
  use mod_particle_types, only: particle_fieldline_id,particle_gc_id,particle_gc_vpar_id
  use mod_particle_types, only: particle_kinetic_id,particle_kinetic_leapfrog_id
  use mod_particle_types, only: particle_kinetic_relativistic_id,particle_gc_relativistic_id
  use mod_particle_types, only: particle_gc_Qin_id
  use mod_particle_common_test_tools, only: fill_particles
  use mod_particle_common_test_tools, only: fill_sim_groups
  use mod_particle_common_test_tools, only: invalidate_particles
  use mod_particle_common_test_tools, only: obtain_active_particle_ids
  use mod_particle_common_test_tools, only: obtain_particle_charges
  use mod_particle_common_test_tools, only: allocate_one_particle_list_type
  implicit none
  !> variables
  integer :: ifail
  
  !> allocate the particle lists
  ifail = 0; allocate(sim_sol%groups(n_particle_types))
  call allocate_one_particle_list_type(n_particle_types,n_particles,sim_sol%groups,ifail)
  call assert_true(ifail.eq.0,"Error particle_sim test setup: particle list not allocated!")

  !> fill-up the group and particle base variables
  call fill_sim_groups(n_particle_types,sim_sol%groups)
  call fill_particles(n_particle_types,n_particles,sim_sol%groups)

  !> invalidate particles in particle lists
  call invalidate_particles(n_particle_types,n_particles,survival_prob,&
  n_active_particles_sol,sim_sol%groups)
  call obtain_active_particle_ids(n_particle_types,n_particles,&
  active_particle_ids_sol,sim_sol%groups)

  !> initialise the particle type list
  particle_type_list_sol = (/particle_fieldline_id,particle_gc_id,particle_gc_vpar_id,&
  particle_kinetic_id,particle_kinetic_leapfrog_id,particle_kinetic_relativistic_id,&
  particle_gc_relativistic_id,particle_gc_Qin_id/)
  
  !> get particle charges
  call obtain_particle_charges(n_particle_types,n_particles,&
  particle_charge_list_sol,sim_sol%groups)
end subroutine setup

!> Tests ----------------------------------------
!> test return particle type code for the whole particle list
subroutine test_codify_particle_sim()
  implicit none
  !> variables
  integer :: ii
  integer,dimension(n_particle_types) :: list_code

  !> extract particle list code 
  call sim_sol%find_particle_types(n_particle_types,list_code)
  call assert_equals(list_code,particle_type_list_sol,n_particle_types,&
  "Error in particle_sim codify particle sim: particle types mismatch!")
end subroutine test_codify_particle_sim

!> test find active particle id interface notype
subroutine test_find_active_particle_id()
  implicit none
  !> variables
  integer :: ii
  integer,dimension(n_particle_types) :: n_particles_array,n_active_particles
  integer,dimension(n_particles,n_particle_types) :: active_particle_ids

  !> compute active particles and their id
  n_particles_array = n_particles; n_active_particles = -1; active_particle_ids = -1;
  call sim_sol%find_active_particles_groups(n_particle_types,n_particles,&
  n_particles_array,n_active_particles,active_particle_ids)
  !> check results
  call assert_equals(n_active_particles,n_active_particles_sol,n_particle_types,&
  "Error particle_sim find active particle notype: n_active_particles mismatch!")
  call assert_equals(active_particle_ids,active_particle_ids_sol,n_particles,n_particle_types,&
  "Error particle_sim find active particle notype: active_particle_ids mismatch!")
end subroutine test_find_active_particle_id

!> test find active particle id interface notype
subroutine test_find_active_particle_id_type()
  implicit none
  !> variables
  integer :: ii
  integer,dimension(n_particle_types) :: n_particles_array,n_active_particles
  integer,dimension(n_particle_types) :: n_act_part_sol_loc
  integer,dimension(n_particles,n_particle_types) :: active_particle_ids
  integer,dimension(n_particles,n_particle_types) :: act_part_ids_sol_loc

  !> compute active particles and their id and tests
  n_particles_array = n_particles;
  do ii=1,n_particle_types
    n_active_particles = -1; active_particle_ids = -1;
    n_act_part_sol_loc = 0; act_part_ids_sol_loc = 0;
    n_act_part_sol_loc(ii) = n_active_particles_sol(ii)
    act_part_ids_sol_loc(:,ii) = active_particle_ids_sol(:,ii)
    call sim_sol%find_active_particles_groups(n_particle_types,n_particles,&
    n_particles_array,n_active_particles,active_particle_ids,particle_type_list_sol(ii))
    call assert_equals(n_active_particles,n_act_part_sol_loc,n_particle_types,&
    "Error particle_sim find active particle type: n_active_particles mismatch!")
    call assert_equals(active_particle_ids,act_part_ids_sol_loc,n_particles,n_particle_types,&
    "Error particle_sim find active particle type: active_particle_ids mismatch!")
  enddo 
end subroutine test_find_active_particle_id_type

!>-----------------------------------------------

end module mod_particle_sim_test
