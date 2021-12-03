!> mod_particle_io_mpi_test contains andata and procedures
!> for testing the writing and the reading of particles
!> data to/from HDF5 files (MPI enables).
module mod_particle_io_mpi_test
use mpi
use mod_particle_sim, only: particle_sim
implicit none

private
public :: run_fruit_particle_io_mpi

!> Variables --------------------------------------------
type(particle_sim) :: sim_particles
!> the number of particle groups is set equal to the number
!> of particle types for testing all of them
integer,parameter :: n_groups=8
integer,parameter :: n_particles=5 !< N# of particles per group per task
!> intervals for random number generation
integer,dimension(2),parameter   :: rng_seed_interval=(/-1234,9876/)
integer,dimension(2),parameter   :: q_interval=(/1,10/)
integer,dimension(2),parameter   :: i_elm_interval=(/1,1000000/)
integer,dimension(2),parameter   :: i_life_interval=(/1,1000000/)
real*8,dimension(2),parameter    :: t_birth_interval=(/0.d0,3.45d4/)
real*8,dimension(2),parameter    :: st_interval=(/0.d0,1.d0/)
real*8,dimension(2),parameter    :: dt_interval=(/1.d-12,1.d-2/)
real*8,dimension(2),parameter    :: mass_interval=(/5.485d-4,124.d0/)
real*8,dimension(2),parameter    :: v_interval=(/-6.75d3,8.45d3/)
real*8,dimension(2),parameter    :: Ekin_interval=(/0.d0,1.d7/)
real*8,dimension(2),parameter    :: mu_interval=(/0.d0,1.d-5/)
real*8,dimension(2),parameter    :: Bnorm_interval=(/0.d0,1.4d1/)
real*8,dimension(2),parameter    :: weight_interval=(/0.d0,1.d3/)
real*8,dimension(3),parameter    :: x_lowbnd=(/-5.d2,-1.d2,1.d0/)
real*8,dimension(3),parameter    :: x_uppbnd=(/7.d2,2.d2,4.d2/)
real*8,dimension(3),parameter    :: vp3d_lowbnd=(/-1.25d3,-7.5d2,-8.d1/)
real*8,dimension(3),parameter    :: vp3d_uppbnd=(/7.5d1,2.35d2,4.85d3/)
real*8,dimension(3),parameter    :: ABE_lowbnd=(/-2.67d0,-9.85d0,0.35d0/)
real*8,dimension(3),parameter    :: ABE_uppbnd=(/0.78d0,2.35d0,5.67d0/)
!> Interfaces -------------------------------------------
contains

!> Fruit basket -----------------------------------------
!> run_fruit_particle_io_mpi performs the set-up, 
!> execution and tear-down of test features
!> inputs:
!>   rank:    (integer) mpi task rank
!>   n_tasks: (integer) number of tasks in the commworld
!>   ifail:   (integer) 0 if success
!> outputs:
!>   ifail:   (integer) 0 if success
subroutine run_fruit_particle_io_mpi(rank,n_tasks,ifail)
  implicit none
  !> inputs
  integer,intent(in) :: rank,n_tasks
  !> inputs-outputs
  integer,intent(inout) :: ifail
  if(rank.eq.0) write(*,'(/A)') "  ... setting-up: particle io mpi tests"
  call setup(rank,n_tasks,ifail)
  if(rank.eq.0) write(*,'(/A)') "  ... running: particle io mpi tests"
  if(rank.eq.0) write(*,'(/A)') "  ... tearing-down: particle io mpi tests"
  call teardown(rank,n_tasks,ifail)
end subroutine run_fruit_particle_io_mpi

!> Set-up and tear-down ---------------------------------
!> set-up the test features
!> inputs:
!>   rank:    (integer) mpi task rank
!>   n_tasks: (integer) number of tasks in the commworld
!>   ifail:   (integer) 0 if success
!> outputs:
!>   ifail:   (integer) 0 if success
subroutine setup(rank,n_tasks,ifail)
  use mod_particle_types, only: particle_kinetic,particle_kinetic_leapfrog
  use mod_particle_types, only: particle_gc,particle_fieldline
  use mod_particle_types, only: particle_kinetic_relativistic
  use mod_particle_types, only: particle_gc_relativistic
  use mod_particle_types, only: particle_gc_vpar,particle_gc_Qin
  use mod_gnu_rng,        only: gnu_rng_interval
  use mod_gnu_rng,        only: set_seed_sys_time
  !$ use omp_lib
  implicit none
  !> inputs
  integer,intent(in) :: rank,n_tasks
  !$ integer :: thread_id
  !> inputs-outputs
  integer,intent(inout) :: ifail
  !> variables
  integer :: ii,jj,rng_integer
  real*8              :: rng_real
  real*8,dimension(2) :: rng_real_size2
  real*8,dimension(3) :: rng_real_size3

  !> initialize the particle simulation (requires jorek inputfile)
  call sim_particles%initialize(n_groups,.false.,rank,n_tasks)

  !> allocate particle lists for different particle types
  allocate(particle_kinetic::sim_particles%groups(1)%particles(n_particles))
  allocate(particle_kinetic_leapfrog::sim_particles%groups(2)%particles(n_particles))
  allocate(particle_gc::sim_particles%groups(3)%particles(n_particles))
  allocate(particle_fieldline::sim_particles%groups(4)%particles(n_particles))
  allocate(particle_kinetic_relativistic::sim_particles%groups(5)%particles(n_particles))
  allocate(particle_gc_relativistic::sim_particles%groups(6)%particles(n_particles))
  allocate(particle_gc_vpar::sim_particles%groups(7)%particles(n_particles))
  allocate(particle_gc_Qin::sim_particles%groups(8)%particles(n_particles))

  !> fill-up the particle_base variables for all particles and all groups
  !$omp parallel default(shared) private(ii,jj,rank,rng_integer,rng_real,&
  !$omp rng_real_size2,rng_real_size3,thread_id)
  thread_id = 1
  !$ thread_id = omp_get_thread_num()
  call set_seed_sys_time(rng_seed_interval,rank,thread_id)
  !$omp do collapse(2)
  do jj=1,n_groups
    do ii=1,n_particles
      call gnu_rng_interval(weight_interval,rng_real)
      call gnu_rng_interval(2,st_interval,rng_real_size2)
      call gnu_rng_interval(3,x_lowbnd,x_uppbnd,rng_real_size3)
      sim_particles%groups(jj)%particles(ii)%x       = rng_real_size3
      sim_particles%groups(jj)%particles(ii)%st      = rng_real_size2
      sim_particles%groups(jj)%particles(ii)%weight  = rng_real
      call gnu_rng_interval(t_birth_interval,rng_real)
      sim_particles%groups(jj)%particles(ii)%t_birth = real(rng_real,kind=4)
      call gnu_rng_interval(i_elm_interval,rng_integer)
      sim_particles%groups(jj)%particles(ii)%i_elm   = rng_integer
      call gnu_rng_interval(i_life_interval,rng_integer)
      sim_particles%groups(jj)%particles(ii)%i_life  = rng_integer 
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  !> fill-up the variables for each specific species

end subroutine setup

!> tear-down the test simulation features
!> inputs:
!>   rank:    (integer) mpi task rank
!>   n_tasks: (integer) number of tasks in the commworld
!>   ifail:   (integer) 0 if success
!> outputs:
!>   ifail:   (integer) 0 if success
subroutine teardown(rank,n_tasks,ifail)
  implicit none
  !> inputs
  integer,intent(in) :: rank,n_tasks
  !> inputs-outputs
  integer,intent(inout) :: ifail

  !> deallocate sim_particles and all structures in it
end subroutine teardown
!> Tests ------------------------------------------------
!> Tools ------------------------------------------------
!>-------------------------------------------------------
end module mod_particle_io_mpi_test
