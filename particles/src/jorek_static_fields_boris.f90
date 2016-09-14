!> Calculate the particle trajectories with the Boris method in
!> static JOREK fields, taking from the first input file argument
program jorek_static_fields_boris
use mod_particle_sim
use mod_event
use mod_action
use mod_main_loop
use mod_boris
use mod_pusher, only: pusher_container
use mod_jorek_fields
use mpi
implicit none

type(particle_sim) :: sim
type(event), dimension(:), allocatable :: events
type(pusher_container), dimension(:), allocatable :: pushers

allocate(sim%fields, source=jorek_fields_static())
allocate(sim%groups(1))
pushers = [pusher_container(pusher_boris(fixed_timestep=1d-9))]

events = [ &
  event(read_jorek_fields(basename='jorek_restart')), &
  event(seed_particles()), &
  event(stop_action(), start=1d-3) & ! Stop the sim after 1 second
]

call main_loop(sim, pushers, events)
end program jorek_static_fields_boris
