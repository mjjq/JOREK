!>#Example 5
!> Multiple groups with different pushers. (the same here as no other pushers
!> are available yet)
!>
!>* fields: prescribed
!>* pusher: [boris,boris2]
!>* geometry: cartesian
!>
!> Compile with: `make ex4`
!> Run with: `./ex4`
!> See the [annotated source](../sourcefile/ex5.f90.html) for details.
program ex5
use particle_tracer
use mod_diag_print_kinetic_energy
implicit none

! 1. Set up the simulation variables containing
!    sim: particles, time, and io.
!    events: halting points for the pushers and actions to run.
real*8 :: timesteps(2) = [1d-6, 1d-3]
type(particle_sim) :: sim
type(event), dimension(:), allocatable :: events
integer :: i, j, k, n_steps
real*8  :: target_time

! 2. Allocate a group and a particle of type particle_kinetic_leapfrog.
call sim%initialize(num_groups=2)
allocate(particle_kinetic_leapfrog::sim%groups(1)%particles(1))
allocate(particle_kinetic_leapfrog::sim%groups(2)%particles(1))

! 3. Initialize the particles.
!    This should usually be done by a dedicated initialization routine
!    or by reading existing files.
select type (p => sim%groups(1)%particles(1))
type is (particle_kinetic_leapfrog)
  p%x = [0.d0,0.d0,0.d0]
  p%v = [1.d0,0.d0,0.d0]
  p%q = 2
  p%m = 4.d0
end select
select type (p => sim%groups(2)%particles(1))
type is (particle_kinetic_leapfrog)
  p%x = [0.d0,0.d0,0.d0]
  p%v = [10.d0,0.d0,0.d0]
  p%q = 2
  p%m = 4.d0
end select

! 4. Set an event to stop the simulation.
events  = [event(stop_action(), start=1.d0), &
           event(diag_print_kinetic_energy(), step=0.1d0)]

! 5. Check whether all events conform to the requested timestep
call check_and_fix_timesteps(timesteps, events)

! 6. Run first events
call with(sim, events, at=0.d0)

! 7. Loop until we the simulation requests a stop
do while (.not. sim%stop_now)
  ! 7.1 Find out which events are next and when they will run
  call next_event_at(events, sim%time, target_time)

  ! 7.2 Loop over all particle groups
  do i=1,2
    n_steps = nint((target_time - sim%time)/timesteps(i))

    ! 7.3 Select the type of this group once to call the right integrator
    select type (particles => sim%groups(i)%particles)
    type is (particle_kinetic_leapfrog)
      ! 7.4 Loop first over particles, and then over how many steps we can take
      !$omp parallel do default(private) &
      !$omp shared(n_steps, timesteps, i)
      do j=1,size(particles)
        do k=1,n_steps
          call boris_push_cartesian(particles(j), E=[0d0,0d0,0d0], B=[0d0,0d0,1d0], dt=timesteps(i))
        end do ! steps
      end do ! particles
      !$omp end parallel do
    end select
  end do ! groups

  ! 7.5 Update the current time and run events
  sim%time = target_time
  call with(sim, events, at=sim%time)
end do
call sim%finalize
end program ex5
