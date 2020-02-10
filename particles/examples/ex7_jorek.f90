!> Push relativistic guiding centre(s) in JOREK fields
!>
!> Compile with `make ex6_jorek`
!> Run with `./ex6_jorek < JOREK_namelist`

program ex7_jorek

use particle_tracer
use mod_particle_io
use mod_particle_diagnostics
use mod_fields_linear   
use mod_fields_hermite_birkhoff                  
implicit none

! Set up the simulation variables
real(kind=8)                      :: timesteps(1) = [3.5723d-13]
real(kind=8)                      :: target_time, t
integer(kind=4)                   :: n_part, i, j, k, n_steps, ifail, n_lost
logical                           :: restart
type(diag_print_kinetic_energy)   :: print_kinetic_energy
type(write_particle_diagnostics)  :: diag

call sim%initialize(num_groups=1)

! Set up the diagnostics output
diag = write_particle_diagnostics(filename='diag.h5',only=[1,2,6,12,13,14]) ! store total and kinetic energies, p_phi, phi, R, Z

restart = .false. !

if (.not. restart) then
  sim%time = 1.d-7  ! start time 

  ! Allocate a group and particle(s) of type particle_gc_relativistic
  n_part = 1
  allocate(particle_gc_relativistic::sim%groups(1)%particles(n_part))
  sim%groups(1)%mass = 5.4857990907016d-4 !< particle mass in AMU

  ! Set events to write output data and stop the simulation.
  ! One can use read_jorek_fields_interp_linear or read_jorek_fields_interp_hermite_birkhoff,
  ! and i=-1 (to read jorek_restart.h5 and keep this field at all time) or i=last_file_before_time(sim%time)
  ! (to read a sequel of jorekXXXXX.h5 files and use time-evolving fields)
  events = [event(read_jorek_fields_interp_linear(i=-1)), & 
            event(diag,start=sim%time,step=1d-8),         &
            event(stop_action(),start=sim%time+5.d-8)]

  ! Run first event to read the JOREK fields
  call with(sim, events, at=0.d0)

  select type (p=>sim%groups(1)%particles(1))
  type is (particle_gc_relativistic)
    p%x = [3.d0,0.d0,0.d0]
    call find_RZ(sim%fields%node_list, sim%fields%element_list, &
                 p%x(1), p%x(2), & ! inputs
                 p%x(1), p%x(2), p%i_elm, p%st(1), p%st(2), ifail) ! outputs
    p%p = [1.6d6,0.]
    p%q = -1
  end select
else
  write(*,*) '*******************************'
  write(*,*) 'Reading particle restart file'
  call read_simulation_hdf5(sim, 'part_restart.h5')
  write(*,*) 'Time = ', sim%time
  write(*,*) '*******************************' 
  ! Set events to write output data and stop the simulation.
  ! One can use read_jorek_fields_interp_linear or read_jorek_fields_interp_hermite_birkhoff,
  ! and i=-1 (to read jorek_restart.h5 and keep this field at all time) or i=last_file_before_time(sim%time)
  ! (to read a sequel of jorekXXXXX.h5 files and use time-evolving fields)
  events = [event(read_jorek_fields_interp_linear(i=-1)), & 
            event(diag,start=sim%time,step=1d-8),         &
	    event(stop_action(),start=sim%time+5.d-8)]

  ! Run first event to read the JOREK fields
  call with(sim, events, at=0.d0)   
endif

! Check all events conform to the requested timestep
call check_and_fix_timesteps(timesteps, events)

! Loop until the simulation is stopped
do while (.not. sim%stop_now)
  ! Extract the next event time
  target_time = next_event_at(sim, events)
  ! Loop over all particle groups
  do i=1,1
    ! Compute the number of steps for the particle
    n_steps = nint((target_time - sim%time)/timesteps(i))
	write(*,*) 'n_steps', n_steps
    n_lost = 0

    select type (particles => sim%groups(i)%particles)
    type is (particle_gc_relativistic)	
!      !$omp parallel do default(private) &
!      !$omp shared (i, n_steps, timesteps, sim) &
!      !$omp reduction(+:n_lost)	
      do j=1,size(particles,1)
        do k=1,n_steps
          if (particles(j)%i_elm .eq. 0) exit
          t = sim%time + k*timesteps(i)
!         PUT PUSHER HERE... TO BE COMPLETED	  
          if (particles(j)%i_elm .eq. 0) n_lost = n_lost + 1		
        end do !< time steps
      end do !< particles
!      !$omp end parallel do
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

end program ex7_jorek

