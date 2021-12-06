!> mod_particle_io_mpi_test contains andata and procedures
!> for testing the writing and the reading of particles
!> data to/from HDF5 files (MPI enables).
module mod_particle_io_mpi_test
use fruit
use mpi
use mod_particle_sim, only: particle_sim
implicit none

private
public :: run_fruit_particle_io_mpi

!> Variables --------------------------------------------
type(particle_sim)          :: sim_particles
character(len=28),parameter :: test_filename="test_particle_io_mpi_hdf5.h5"
!> the number of particle groups is set equal to the number
!> of particle types for testing all of them
!> particle_gc_Qin is commented because the I/O for particle_gc_Qin
!> has not been implemented yet
integer,parameter :: n_groups=7!8
integer,parameter :: n_particles=5 !< N# of particles per group per task
real*8,parameter :: tol_real8=1.d-15
!> intervals for random number generation
integer,dimension(2),parameter   :: rng_seed_interval=(/-1234,9876/)
integer,dimension(2),parameter   :: q_interval=(/1,100/)
integer,dimension(2),parameter   :: i_elm_interval=(/1,1000000/)
integer,dimension(2),parameter   :: i_life_interval=(/1,1000000/)
real*8,dimension(2),parameter    :: sim_time_interval=(/0.d0,1.d3/)
real*8,dimension(2),parameter    :: t_birth_interval=(/0.d0,3.45d4/)
real*8,dimension(2),parameter    :: st_interval=(/0.d0,1.d0/)
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
  call test_particle_mpi_io(rank,n_tasks,ifail)
  call test_get_simulation_hdf5_time()
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
  use mod_particle_io,    only: write_simulation_hdf5
  use mod_gnu_rng,        only: gnu_rng_interval
  implicit none
  !> inputs
  integer,intent(in) :: rank,n_tasks
  !> variables
  integer :: ii,jj

  !> inputs-outputs
  integer,intent(inout) :: ifail

  !> initialize the particle simulation (requires jorek inputfile)
  call sim_particles%initialize(n_groups,.false.,rank,n_tasks)

  !> set and broadcast simulation time
  if(rank.eq.0) then
    call gnu_rng_interval(sim_time_interval,sim_particles%time)
  endif
  call MPI_Bcast(sim_particles%time,1,MPI_REAL8,0,MPI_COMM_WORLD,ifail)

  !> allocate particle lists for different particle types
  allocate(particle_fieldline::sim_particles%groups(1)%particles(n_particles))
  allocate(particle_gc::sim_particles%groups(2)%particles(n_particles))
  allocate(particle_gc_vpar::sim_particles%groups(3)%particles(n_particles))
  allocate(particle_kinetic::sim_particles%groups(4)%particles(n_particles))
  allocate(particle_kinetic_leapfrog::sim_particles%groups(5)%particles(n_particles))
  allocate(particle_kinetic_relativistic::sim_particles%groups(6)%particles(n_particles))
  allocate(particle_gc_relativistic::sim_particles%groups(7)%particles(n_particles))
  !allocate(particle_gc_Qin::sim_particles%groups(8)%particles(n_particles)) !< IO not implemented

  !> fill-up the groupd and particle base variables
  call fill_sim_groups(rank,ifail)
  call fill_particle_base(rank)

  !> fill up variables for each species
  call fill_particle_fieldline(rank,sim_particles%groups(1)%particles)
  call fill_particle_gc(rank,sim_particles%groups(2)%particles)
  call fill_particle_gc_vpar(rank,sim_particles%groups(3)%particles)
  call fill_particle_kinetic(rank,sim_particles%groups(4)%particles)
  call fill_particle_kinetic_leapfrog(rank,sim_particles%groups(5)%particles)
  call fill_particle_kinetic_relativistic(rank,sim_particles%groups(6)%particles)
  call fill_particle_gc_relativistic(rank,sim_particles%groups(7)%particles)
  !call fill_particle_gc_Qin(rank,sim_particles%groups(8)%particles) !< I/O not implemented 

  !> write default simulation in file and read it in new simulation
  call write_simulation_hdf5(sim_particles,trim(test_filename))

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

  call MPI_Barrier(MPI_COMM_WORLD,ifail)
  !> remove test file
  if(rank.eq.0) call system("rm "//test_filename)
end subroutine teardown
!> Tests ------------------------------------------------
!> procedure for testing the particle io
subroutine test_particle_mpi_io(rank,n_tasks,ifail)
  use mod_particle_assert_equal, only: assert_equal_particle
  use mod_particle_sim,          only: particle_sim
  use mod_particle_io,           only: read_simulation_hdf5
  implicit none
  !> inputs
  integer,intent(in) :: rank,n_tasks
  !> inputs-outputs
  integer,intent(inout) :: ifail
  !> variables
  type(particle_sim) :: sim_particles_new
  integer :: ii
  real*8 :: comp_real8_1,comp_real8_2

  !> initialize the new particle simulation (requires jorek inputfile)
  call sim_particles_new%initialize(n_groups,.false.,rank,n_tasks)

  !> read default simulation from file and store in new sim
  call read_simulation_hdf5(sim_particles_new,trim(test_filename))

  !> compu variables which are not read from hdf5
  call copy_sim_fieldline_B_hat_prev(sim_particles,sim_particles_new)

  !> check simulation 
  call assert_equals(sim_particles_new%time,sim_particles%time,tol_real8,&
  "Error writing/reading particle simulation: time mismatch!")

  !> check groups
  do ii=1,n_groups
    call assert_equals(sim_particles_new%groups(ii)%mass,sim_particles%groups(ii)%mass,&
    tol_real8,"Error writing/reading particle simulation: mass mismatch!")
    call assert_equal_particle(n_particles,sim_particles_new%groups(ii)%particles,&
    sim_particles%groups(ii)%particles)
  enddo
end subroutine test_particle_mpi_io

!> subroutine for testing the reading of simulation time
subroutine test_get_simulation_hdf5_time()
  use mod_particle_io, only: get_simulation_hdf5_time
  implicit none
  real*8 :: time_new
  time_new = get_simulation_hdf5_time(test_filename)
  call assert_equals(time_new,sim_particles%time,tol_real8,&
  "Error get simulation time hdf5: time mismatch!")
end subroutine test_get_simulation_hdf5_time

!> Tools ------------------------------------------------
!> copy fieldlines B_hat between two simulations 
!> used for IO because it is not stored in hdf5
subroutine copy_sim_fieldline_B_hat_prev(sim_in,sim_out)
  use mod_particle_types, only: particle_fieldline
  use mod_particle_sim,   only: particle_sim
  implicit none
  type(particle_sim),intent(in)    :: sim_in
  type(particle_sim),intent(inout) :: sim_out
  integer :: ii,jj
  !$omp parallel default(private) shared(sim_in,sim_out)
  do ii=1,n_groups
    !$omp do
    do jj=1,n_particles
      select type (p_out=>sim_out%groups(ii)%particles(jj))
      type is (particle_fieldline)
        select type (p_in=>sim_in%groups(ii)%particles(jj))
          type is (particle_fieldline)
          p_out%B_hat_prev = p_in%B_hat_prev
        end select
      end select
    enddo
    !$omp end do
  enddo
  !$omp end parallel
end subroutine copy_sim_fieldline_B_hat_prev

!> generate random values for filling the groups type
!> we do not create random mass for all groups because
!> it seems that only rank 0 is saved in hdf5
subroutine fill_sim_groups(rank,ifail)
  use mod_gnu_rng, only: gnu_rng_interval
  use mod_gnu_rng, only: set_seed_sys_time
  implicit none
  !> inputs
  integer,intent(in) :: rank
  !> input-outputs
  integer,intent(inout) :: ifail
  !> variables
  integer :: ii
  real*8 :: rn_real
  real*8,dimension(n_groups) :: val_to_bcst
  if(rank.eq.0) then
    call set_seed_sys_time(rng_seed_interval,rank)
    do ii=1,n_groups
      call gnu_rng_interval(mass_interval,rn_real)
      sim_particles%groups(ii)%mass = rn_real
      val_to_bcst(ii) = sim_particles%groups(ii)%mass
    enddo
  endif
  call MPI_Bcast(val_to_bcst,n_groups,MPI_REAL8,0,MPI_COMM_WORLD,ifail)
  if(rank.ne.0) then
    do ii=1,n_groups
      sim_particles%groups(ii)%mass = val_to_bcst(ii)
    enddo
  endif
end subroutine fill_sim_groups

!> generate random values for filling the particle base type
subroutine fill_particle_base(rank)
  use mod_gnu_rng, only: gnu_rng_interval
  use mod_gnu_rng, only: set_seed_sys_time
  !$ use omp_lib
  implicit none
  !> inputs
  integer,intent(in) :: rank
  !> variables
  integer :: ii,jj,rn_integer
  !$ integer :: thread_id
  real*8              :: rn_real
  real*8,dimension(2) :: rn_real_size2
  real*8,dimension(3) :: rn_real_size3

  !> fill-up the particle_base variables for all particles and all groups
  !$omp parallel default(shared) private(ii,jj,rank,rn_integer,rn_real,&
  !$omp rn_real_size2,rn_real_size3,thread_id)
  thread_id = 0
  !$ thread_id = omp_get_thread_num()
  call set_seed_sys_time(rng_seed_interval,rank,thread_id)
  !$omp do collapse(2)
  do jj=1,n_groups
    do ii=1,n_particles
      call gnu_rng_interval(weight_interval,rn_real)
      call gnu_rng_interval(2,st_interval,rn_real_size2)
      call gnu_rng_interval(3,x_lowbnd,x_uppbnd,rn_real_size3)
      sim_particles%groups(jj)%particles(ii)%x       = rn_real_size3
      sim_particles%groups(jj)%particles(ii)%st      = rn_real_size2
      sim_particles%groups(jj)%particles(ii)%weight  = rn_real
      call gnu_rng_interval(t_birth_interval,rn_real)
      sim_particles%groups(jj)%particles(ii)%t_birth = real(rn_real,kind=4)
      call gnu_rng_interval(i_elm_interval,rn_integer)
      sim_particles%groups(jj)%particles(ii)%i_elm   = rn_integer
      call gnu_rng_interval(i_life_interval,rn_integer)
      sim_particles%groups(jj)%particles(ii)%i_life  = rn_integer 
    enddo
  enddo
  !$omp end do
  !$omp end parallel
end subroutine fill_particle_base

!> fill up particle_fieldline with random numbers
subroutine fill_particle_fieldline(rank,particles)
  use mod_particle_types, only: particle_base
  use mod_particle_types, only: particle_fieldline
  use mod_gnu_rng,        only: gnu_rng_interval
  use mod_gnu_rng,        only: set_seed_sys_time
  !$ use omp_lib
  implicit none
  !> inputs
  integer,intent(in) :: rank
  class(particle_base),dimension(n_particles),intent(inout) :: particles
  !> variables
  integer :: ii
  !$ integer :: thread_id
  real*8              :: rn_real
  real*8,dimension(3) :: rn_real_size3

  !$omp parallel default(private) shared(particles)
  thread_id = 0
  !$ thread_id = omp_get_thread_num()
  call set_seed_sys_time(rng_seed_interval,rank,thread_id)
  !$omp do
  do ii=1,n_particles
    select type (p=>particles(ii))
      type is (particle_fieldline)
      call gnu_rng_interval(3,ABE_lowbnd,ABE_uppbnd,rn_real_size3)
      call gnu_rng_interval(v_interval,rn_real)
      p%B_hat_prev = rn_real_size3
      p%v = rn_real
    end select
  enddo
  !$omp end do
  !$omp end parallel
end subroutine fill_particle_fieldline

!> fill up particle_gc with random numbers
subroutine fill_particle_gc(rank,particles)
  use mod_particle_types, only: particle_base
  use mod_particle_types, only: particle_gc
  use mod_gnu_rng,        only: gnu_rng_interval
  use mod_gnu_rng,        only: set_seed_sys_time
  !$ use omp_lib
  implicit none
  !> inputs
  integer,intent(in) :: rank
  class(particle_base),dimension(n_particles),intent(inout) :: particles
  !> variables
  integer :: ii,rn_integer
  !$ integer :: thread_id
  real*8 :: rn_real

  !$omp parallel default(private) shared(particles)
  thread_id = 0
  !$ thread_id = omp_get_thread_num()
  call set_seed_sys_time(rng_seed_interval,rank,thread_id)
  !$omp do
  do ii=1,n_particles
    select type (p=>particles(ii))
      type is (particle_gc)
      call gnu_rng_interval(Ekin_interval,rn_real)   
      p%E = rn_real
      call gnu_rng_interval(mu_interval,rn_real)
      p%mu = rn_real
      call gnu_rng_interval(q_interval,rn_integer)
      p%q = int(rn_integer,kind=1)
    end select
  enddo
  !$omp end do
  !$omp end parallel
end subroutine fill_particle_gc

!> fill up particle_gc_vpar with random numbers
subroutine fill_particle_gc_vpar(rank,particles)
  use mod_particle_types, only: particle_base
  use mod_particle_types, only: particle_gc_vpar
  use mod_gnu_rng,        only: gnu_rng_interval
  use mod_gnu_rng,        only: set_seed_sys_time
  !$ use omp_lib
  implicit none
  !> inputs
  integer,intent(in) :: rank
  class(particle_base),dimension(n_particles),intent(inout) :: particles
  !> variables
  integer :: ii,rn_integer
  !$ integer :: thread_id
  real*8 :: rn_real

  !$omp parallel default(private) shared(particles)
  thread_id = 0
  !$ thread_id = omp_get_thread_num()
  call set_seed_sys_time(rng_seed_interval,rank,thread_id)
  !$omp do
  do ii=1,n_particles
    select type (p=>particles(ii))
      type is (particle_gc_vpar)
      call gnu_rng_interval(v_interval,rn_real)   
      p%vpar = rn_real
      call gnu_rng_interval(mu_interval,rn_real)
      p%mu = rn_real
      call gnu_rng_interval(q_interval,rn_integer)
      p%q = int(rn_integer,kind=1)
    end select
  enddo
  !$omp end do
  !$omp end parallel
end subroutine fill_particle_gc_vpar

!> fill up particle_gc_Qin with random numbers
subroutine fill_particle_gc_Qin(rank,particles)
  use mod_particle_types, only: particle_base
  use mod_particle_types, only: particle_gc_Qin
  use mod_gnu_rng,        only: gnu_rng_interval
  use mod_gnu_rng,        only: set_seed_sys_time
  !$ use omp_lib
  implicit none
  !> inputs
  integer,intent(in) :: rank
  class(particle_base),dimension(n_particles),intent(inout) :: particles
  !> variables
  integer :: ii,rn_integer
  !$ integer :: thread_id
  real*8              :: rn_real
  real*8,dimension(3) :: rn_real_size3
  real*8,dimension(3,3) :: rn_real_size33

  !$omp parallel default(private) shared(particles)
  thread_id = 0
  !$ thread_id = omp_get_thread_num()
  call set_seed_sys_time(rng_seed_interval,rank,thread_id)
  !$omp do
  do ii=1,n_particles
    select type (p=>particles(ii))
      type is (particle_gc_Qin)
      call gnu_rng_interval(3,x_lowbnd,x_uppbnd,rn_real_size3)
      p%x_m = rn_real_size3       !< x_m
      call gnu_rng_interval(3,ABE_lowbnd,ABE_uppbnd,rn_real_size3)
      p%Astar_m = rn_real_size3   !< Astar_m
      call gnu_rng_interval(3,ABE_lowbnd,ABE_uppbnd,rn_real_size3)
      p%Astar_k = rn_real_size3   !< Astar_k
      call gnu_rng_interval(3,ABE_lowbnd,ABE_uppbnd,rn_real_size33(:,1))
      call gnu_rng_interval(3,ABE_lowbnd,ABE_uppbnd,rn_real_size33(:,2))
      call gnu_rng_interval(3,ABE_lowbnd,ABE_uppbnd,rn_real_size33(:,3))
      p%dAstar_k = rn_real_size33 !< dAstar_k  
      call gnu_rng_interval(3,ABE_lowbnd,ABE_uppbnd,rn_real_size3)
      p%dBn_k = rn_real_size3     !< dBn_k
      call gnu_rng_interval(3,ABE_lowbnd,ABE_uppbnd,rn_real_size3)
      p%Bnorm_k = rn_real_size3   !< Bnorm_l
      call gnu_rng_interval(3,ABE_lowbnd,ABE_uppbnd,rn_real_size3)
      p%E_k = rn_real_size3       !< E_k
      call gnu_rng_interval(v_interval,rn_real)
      p%vpar_m = rn_real          !< vpar_m
      call gnu_rng_interval(Bnorm_interval,rn_real)
      p%Bn_k = rn_real            !< Bn_k
      call gnu_rng_interval(v_interval,rn_real)
      p%vpar = rn_real            !< vpar
      call gnu_rng_interval(v_interval,rn_real)
      p%mu = rn_real              !< mu
      call gnu_rng_interval(q_interval,rn_integer)
      p%q = int(rn_integer,kind=1)
    end select
  enddo
  !$omp end do
  !$omp end parallel
end subroutine fill_particle_gc_Qin

!> fill up particle_kinetic with random numbers
subroutine fill_particle_kinetic(rank,particles)
  use mod_particle_types, only: particle_base
  use mod_particle_types, only: particle_kinetic
  use mod_gnu_rng,        only: gnu_rng_interval
  use mod_gnu_rng,        only: set_seed_sys_time
  !$ use omp_lib
  implicit none
  !> inputs
  integer,intent(in) :: rank
  class(particle_base),dimension(n_particles),intent(inout) :: particles
  !> variables
  integer :: ii,rn_integer
  !$ integer :: thread_id
  real*8,dimension(3) :: rn_real_size3

  !$omp parallel default(private) private(particles)
  thread_id = 0
  !$ thread_id = omp_get_thread_num()
  call set_seed_sys_time(rng_seed_interval,rank,thread_id)
  !$omp do
  do ii=1,n_particles
    select type (p=>particles(ii))
      type is (particle_kinetic)
      call gnu_rng_interval(3,vp3d_lowbnd,vp3d_uppbnd,rn_real_size3)
      p%v = rn_real_size3
      call gnu_rng_interval(q_interval,rn_integer)
      p%q = int(rn_integer,kind=1)
    end select
  enddo
  !$omp end do
  !$omp end parallel
end subroutine fill_particle_kinetic

!> fill up particle_kinetic_leapfrog with random numbers
subroutine fill_particle_kinetic_leapfrog(rank,particles)
  use mod_particle_types, only: particle_base
  use mod_particle_types, only: particle_kinetic_leapfrog
  use mod_gnu_rng,        only: gnu_rng_interval
  use mod_gnu_rng,        only: set_seed_sys_time
  !$ use omp_lib
  implicit none
  !> inputs
  integer,intent(in) :: rank
  class(particle_base),dimension(n_particles),intent(inout) :: particles
  !> variables
  integer :: ii,rn_integer
  !$ integer :: thread_id
  real*8,dimension(3) :: rn_real_size3

  !$omp parallel default(private) private(particles)
  thread_id = 0
  !$ thread_id = omp_get_thread_num()
  call set_seed_sys_time(rng_seed_interval,rank,thread_id)
  !$omp do
  do ii=1,n_particles
    select type (p=>particles(ii))
      type is (particle_kinetic_leapfrog)
      call gnu_rng_interval(3,vp3d_lowbnd,vp3d_uppbnd,rn_real_size3)
      p%v = rn_real_size3
      call gnu_rng_interval(q_interval,rn_integer)
      p%q = int(rn_integer,kind=1)
    end select
  enddo
  !$omp end do
  !$omp end parallel
end subroutine fill_particle_kinetic_leapfrog

!> fill up particle_kinetic_relativistic with random numbers
subroutine fill_particle_kinetic_relativistic(rank,particles)
  use mod_particle_types, only: particle_base
  use mod_particle_types, only: particle_kinetic_relativistic
  use mod_gnu_rng,        only: gnu_rng_interval
  use mod_gnu_rng,        only: set_seed_sys_time
  !$ use omp_lib
  implicit none
  !> inputs
  integer,intent(in) :: rank
  class(particle_base),dimension(n_particles),intent(inout) :: particles
  !> variables
  integer :: ii,rn_integer
  !$ integer :: thread_id
  real*8,dimension(3) :: rn_real_size3

  !$omp parallel default(private) private(particles)
  thread_id = 0
  !$ thread_id = omp_get_thread_num()
  call set_seed_sys_time(rng_seed_interval,rank,thread_id)
  !$omp do
  do ii=1,n_particles
    select type (p=>particles(ii))
      type is (particle_kinetic_relativistic)
      call gnu_rng_interval(3,vp3d_lowbnd,vp3d_uppbnd,rn_real_size3)
      p%p = rn_real_size3
      call gnu_rng_interval(q_interval,rn_integer)
      p%q = int(rn_integer,kind=1)
    end select
  enddo
  !$omp end do
  !$omp end parallel
end subroutine fill_particle_kinetic_relativistic

!> fill up particle_gc_relativistic with random numbers
subroutine fill_particle_gc_relativistic(rank,particles)
  use mod_particle_types, only: particle_base
  use mod_particle_types, only: particle_gc_relativistic
  use mod_gnu_rng,        only: gnu_rng_interval
  use mod_gnu_rng,        only: set_seed_sys_time
  !$ use omp_lib
  implicit none
  !> inputs
  integer,intent(in) :: rank
  class(particle_base),dimension(n_particles),intent(inout) :: particles
  !> variables
  integer :: ii,rn_integer
  !$ integer :: thread_id
  real*8,dimension(2) :: rn_real_size2

  !$omp parallel default(private) shared(particles)
  thread_id = 0
  !$ thread_id = omp_get_thread_num()
  call set_seed_sys_time(rng_seed_interval,rank,thread_id)
  !$omp do
  do ii=1,n_particles
    select type (p=>particles(ii))
      type is (particle_gc_relativistic)
      call gnu_rng_interval(2,vp3d_lowbnd(1:2),vp3d_uppbnd(1:2),rn_real_size2)
      p%p = rn_real_size2
      call gnu_rng_interval(q_interval,rn_integer)
      p%q = int(rn_integer,kind=1)
    end select
  enddo
  !$omp end do
  !$omp end parallel
end subroutine fill_particle_gc_relativistic

!>-------------------------------------------------------
end module mod_particle_io_mpi_test
