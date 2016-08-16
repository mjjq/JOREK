program particles

use mod_particle_sim
use mod_event
use mod_action
use mod_jorek_fields
implicit none

type(particle_sim) :: sim
type(event), dimension(:), allocatable :: events

integer :: i

! events
!call read_parameters(sim)
!call allocate(particle_boris::particles(sim%n_particles))
!call read_adas
!call initialize_rng(rng)
!call initialize_particles(sim)

events = [ &
  event(stop_action(), start=1.d0) & ! Stop the sim after 1 second
!, event(read_all()) & ! Read an existing state at t=0
!, event(write_state, step=1.d-1) & ! Write output every 0.1 second
, event(read_jorek_fields(), step=1.d-2) & ! Read new fields every 0.01 second
]



end program particles
