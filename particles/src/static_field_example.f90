program static_field_example
use particle_tracer
use mpi
implicit none

type(particle_sim) :: sim
type(pusher_container), dimension(:), allocatable :: pushers
type(event), dimension(:), allocatable :: events

allocate(sim%fields, source=prescribed_fields(CARTESIAN, E, B))
allocate(sim%groups(1))
allocate(particle_boris::sim%groups(1)%particles(1))

pushers = [ &
  pusher_container(pusher_boris(fixed_timestep=0.1d0)) &
]

events = [ &
  event(stop_action(), start=1.d0) & ! Stop the sim after 1 second
]

call main_loop(sim, pushers, events)
write(*,*) sim%groups(1)%particles(1)%x

contains
pure function E(x, t)
  real*8, intent(in) :: x(3), t
  real*8 :: E(3)
  E = [0,0,0]
end function E
pure function B(x, t)
  real*8, intent(in) :: x(3), t
  real*8 :: B(3)
  B = [0,0,1]
end function B
subroutine init_particle(p)
  type(particle_boris), intent(inout) :: p
  p%x = [0,0,0]
  p%v = [1,0,0]
  p%q = 2
  p%m = 4.d0
  p%lost = .false.
end subroutine init_particle
end program static_field_example
