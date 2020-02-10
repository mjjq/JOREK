!># Example 6
!> Push a single relativistic particle with Cartesian 5-steps Volume Preserving 
!> scheme in static and uniform fields.
!>
!>* Fields:   prescribed: B_x=0, B_y=0, B_z=1 [T] 
!>* Pusher:   5-steps Volume Preserving Scheme
!>* Geometry: cartesian
!>* Particle type: electron
!>* Kinetic energy: 10 MeV
!>* Gyro-period: 3.5723d-13 [s]
!>
!> Note: Only 100 gyro-orbits are considered for fast pusher testing
!>
!> Compile with `make ex6`
!> Run with `./ex6`
program ex6
use particle_tracer
implicit none

! 1. Set up the simulation variables
!    sim: particles, time, io
!    events: pusher halting points and actions to run
real(kind=8)    :: timesteps(1) = [3.5723d-13]
integer(kind=4) :: i,j,k,n_steps
real(kind=8)    :: target_time
type(particle_kinetic_relativistic) :: particle !< define a relativistic particle
type(diag_print_kinetic_energy) :: print_kinetic_energy

! 2. allocate a group and a particles of type particle_kinetic_relativistic
call sim%initialize(num_groups=1)
allocate(particle_kinetic_relativistic::sim%groups(1)%particles(1))

! 3. Hard-coded particle initialisation
sim%groups(1)%mass = 5.4857990907016d-4 !< particle mass in AMU
! select the particle type
select type (p=>sim%groups(1)%particles(1))
type is (particle_kinetic_relativistic)
  p%x = [0.d0,0.d0,0.d0]
  p%p = [3.37886d+6,0.d0,0.d0]
  p%q = -1
end select

! 4. Set an event to stop the simulation at 100 gyro-periods
events = [event(stop_action(), start=3.5723d-9)]

! 5. Check all events conform to the requested timestep
call check_and_fix_timesteps(timesteps, events)

! 6. Run first event
call with(sim, events, at=0.d0)

! 7. Loop until the simulation is stopped
do while (.not. sim%stop_now)
  ! 7.1 extract the next event time
  target_time = next_event_at(sim, events)
  ! 7.2 loop over all particle groups
  do i=1,1
    ! 7.3 compute the number of steps for the particle
    n_steps = nint((target_time - sim%time)/timesteps(i))
    ! 7.4 loop on the particle whithin the i-th group
    !$omp parallel do default(private) shared (i,n_steps,timesteps,sim)
    do j=1,size(sim%groups(i)%particles)
      ! 7.5 copy the i-th particle in particle
      particle = sim%groups(i)%particles(j) 
      ! 7.6 loop on the time steps
      do k=1,n_steps
          ! 7.7 integrate via volume preserving algorithm
          call volume_preserving_push_cartesian(particle,mass=sim%groups(i)%mass,&
               E=[0.d0,0.d0,0.d0],B=[0.d0,0.d0,1.d0],dt=timesteps(i))
      enddo !< time steps
      ! 7.8 copy back the particle into particle groups
      sim%groups(i)%particles(j) = particle
    enddo !< particles
    !$omp end parallel do
  enddo !< groups

  ! 7.9 update current time and run events
  sim%time = target_time
  call with(sim, events, at=sim%time)
enddo !< event

! 8. Print particle information
write(*,*) norm2(sim%groups(1)%particles(1)%x) !< print particle final position
call print_kinetic_energy%do(sim,events(1))  !< print particle kinetic energy

! 9. Finilize the simulation
call sim%finalize

end program ex6

