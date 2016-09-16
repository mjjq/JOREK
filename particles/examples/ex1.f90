!>#Example 1
!> 
!>* fields: prescribed
!>* pusher: boris
!>* geometry: cartesian
!>
!> Compile with: `make ex1`
!> 
!> Sample runs: `./ex1`
!>
!>## Description
!> This example follows a particle in a static, uniform magnetic field
!> in the z-direction of strength 1 Tesla.
program ex1
use particle_tracer
implicit none

! 1. Set up the simulation variables: sim, pushers, events
type(particle_sim) :: sim
type(pusher_container), dimension(:), allocatable :: pushers
type(event), dimension(:), allocatable :: events

! 2. Set up the fields to be used in the simulation (E and B defined below)
allocate(sim%fields, source=prescribed_fields(CARTESIAN, E, B))

! 3. Allocate groups and particles
allocate(sim%groups(1))
allocate(particle_boris::sim%groups(1)%particles(1))

! 4. Initialize a particle
select type (p => sim%groups(1)%particles(1))
type is (particle_boris)
  p%x = [0.d0,0.d0,0.d0]
  p%v = [1.d0,0.d0,0.d0]
  p%q = 2
  p%m = 4.d0
end select

! 5. Set up pushers
pushers = [pusher_container(pusher_boris(fixed_timestep=0.1d0))]

! 6. Set an event to stop the simulation
events  = [event(stop_action(), start=1.d0)]

! 7. Run the main loop
call main_loop(sim, pushers, events)


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
end program ex1
