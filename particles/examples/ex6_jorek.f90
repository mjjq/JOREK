!># Example 6 using JOREK fields
!> Push a single relativistic particle with Volume Preserving Scheme
!> in static JOREK fields, file jorek_restart.h5 or jorek_restart.rst
!>
!>* Particle type: electron
!>* Kinetic energy: 10 MeV
!>* Gyro-period: 3.5723d-13 [s]
!>
!> Compile with `make ex6_jorek`
!> Run with `./ex6_jorek`
program ex6_jorek
use particle_tracer
use mod_particle_io
implicit none

! Set up the simulation variables
real(kind=8)    :: timesteps(1) = [3.5723d-13]
integer(kind=4) :: i, j, k, n_steps, i_elm_old, ifail, n_lost
real(kind=8)    :: target_time, t
real*8 :: E(3), B(3), rz_old(2), st_old(2), psi, U
type(read_jorek_fields_interp_linear) :: fieldreader
type(particle_kinetic_relativistic)   :: particle !< define a relativistic particle
type(diag_print_kinetic_energy)       :: print_kinetic_energy

! Allocate a group and a particle of type particle_kinetic_relativistic
call sim%initialize(num_groups=1)
allocate(particle_kinetic_relativistic::sim%groups(1)%particles(1))
!allocate(particle_kinetic_leapfrog::sim%groups(1)%particles(1))

sim%groups(1)%mass = 5.4857990907016d-4 !< particle mass in AMU

! Set up the field reader
fieldreader = read_jorek_fields_interp_linear(basename='jorek', i=-1)
call with(sim, fieldreader)

select type (p=>sim%groups(1)%particles(1))
type is (particle_kinetic_relativistic)
!type is (particle_kinetic_leapfrog)
  p%x = [3.d0,0.d0,0.d0]
  call find_RZ(sim%fields%node_list, sim%fields%element_list, &
               p%x(1), p%x(2), & ! inputs
               p%x(1), p%x(2), p%i_elm, p%st(1), p%st(2), ifail) ! outputs
  p%p = [0.d0,0.d0,3.37886d+6]
!  p%v = [0.d0,0.d0,1.0d+7]
  p%q = -1
end select

! Set an event to stop the simulation
events = [event(stop_action(), start=3.5723d-7)]

! Check all events conform to the requested timestep
call check_and_fix_timesteps(timesteps, events)

! Run first event
call with(sim, events, at=0.d0)

! Loop until the simulation is stopped
do while (.not. sim%stop_now)
  ! Extract the next event time
  target_time = next_event_at(sim, events)
  ! Loop over all particle groups
  do i=1,1
    ! Compute the number of steps for the particle
    n_steps = nint((target_time - sim%time)/timesteps(i))
    n_lost = 0

    select type (particles => sim%groups(i)%particles)
!    type is (particle_kinetic_leapfrog)
    type is (particle_kinetic_relativistic)	
      ! Loop on the particle whithin the i-th group
      !$omp parallel do default(private) &
	  !$omp shared (i, n_steps, timesteps, sim) &
      !$omp reduction(+:n_lost)	
      do j=1,size(particles,1)
        do k=1,n_steps
          if (particles(j)%i_elm .eq. 0) exit
          t = sim%time + k*timesteps(i)

!          call sim%fields%calc_EBpsiU(t, particles(j)%i_elm, particles(j)%st, particles(j)%x(3), E, B, psi, U)
!          rz_old    = particles(j)%x(1:2)
!          st_old    = particles(j)%st
!          i_elm_old = particles(j)%i_elm
!          call boris_push_cylindrical(particles(j), sim%groups(i)%mass, E, B, timesteps(i))
!          call find_RZ_nearby(sim%fields%node_list, sim%fields%element_list, rz_old(1), rz_old(2), st_old(1), st_old(2), i_elm_old, &
!              particles(j)%x(1), particles(j)%x(2), particles(j)%st(1), particles(j)%st(2), particles(j)%i_elm, ifail)

          call volume_preserving_push_jorek(particles(j),sim%fields,sim%groups(i)%mass,t,timesteps(i),ifail)

		  if (particles(j)%i_elm .eq. 0) n_lost = n_lost + 1		
        end do !< time steps
      end do !< particles
      !$omp end parallel do
    end select
    write(*,*) "number of lost particles: ", n_lost	  
  enddo !< groups

  ! Update current time and run events
  sim%time = target_time
  call with(sim, events, at=sim%time)
enddo !< event

! Print particle information
write(*,*) 'Final x,y,z: ', sim%groups(1)%particles(1)%x(1), sim%groups(1)%particles(1)%x(2), sim%groups(1)%particles(1)%x(3)

call print_kinetic_energy%do(sim,events(1))  !< print particle kinetic energy

call write_simulation_hdf5(sim, 'part_restart.h5')

! Finalize the simulation
call sim%finalize

end program ex6_jorek

