program static_field_example
use mod_constants, only: CARTESIAN
use mod_particle_sim, only: particle_sim
use mod_event, only: event
use mod_action, only: stop_action
use mod_prescribed_fields, only: prescribed_fields
use mod_main_loop, only: main_loop
use mod_boris, only: pusher_boris, new_pusher_boris
use mod_pusher, only: pusher_container
use mpi
implicit none

type(particle_sim) :: sim
type(pusher_container), dimension(:), allocatable :: pushers
type(event), dimension(:), allocatable :: events

allocate(sim%fields, source=prescribed_fields(CARTESIAN, E, B))
allocate(sim%groups(1))

pushers = [ &
  pusher_container(pusher_boris(fixed_timestep=0.1d0)) &
]

events = [ &
!  event(seed_particles()), &
  event(stop_action(), start=1.d0) & ! Stop the sim after 1 second
]

call main_loop(sim, pushers, events)

contains
pure function E(x, t)
  real*8, dimension(3), intent(in) :: x
  real*8, intent(in) :: t
  real*8, dimension(3) :: E
  E = (/0.d0, 0.d0, 0.d0/)
end function E
pure function B(x, t)
  real*8, dimension(3), intent(in) :: x
  real*8, intent(in) :: t
  real*8, dimension(3) :: B
  B = (/0.d0, 0.d0, 1.d0/)
end function B
end program static_field_example
