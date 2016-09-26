!>#Example 2
!> Demonstrate input/output capabilities
!>
!>* fields: prescribed
!>* pusher: boris
!>* geometry: cartesian
!>
!> Compile with: `make ex2`
!> Run with: `./ex2`
!> See the [annotated source](../sourcefile/ex2.f90.html) for details.
!>
!>## Description
!> This example follows a particle in a static, uniform magnetic field
!> in the z-direction of strength 1 Tesla,
!> and writes it out to a file every 0.1 second.
program ex2
use particle_tracer
implicit none


! 1. Set up the simulation variables containing
!    sim: particles, time, and io.
!    pushers: information on timestepping (and hooks).
!    events: halting points for the pushers and actions to run.
type(particle_sim) :: sim
type(pusher_container), dimension(:), allocatable :: pushers
type(event), dimension(:), allocatable :: events
integer :: ierr

! 2. Set up MPI IO, needed for writing
call MPI_Init(ierr)

! 3. Set up the fields to be used in the simulation. (E and B defined below)
allocate(sim%fields, source=prescribed_fields(CARTESIAN, E_zero, B_z))

! 4. Allocate a group and a particle of type particle_boris.
allocate(sim%groups(1))
allocate(particle_kinetic_leapfrog::sim%groups(1)%particles(1))

! 5. Initialize the particle.
!    This should usually be done by a dedicated initialization routine
!    or by reading existing files.
select type (p => sim%groups(1)%particles(1))
type is (particle_kinetic_leapfrog)
  p%x = [0.d0,0.d0,0.d0]
  p%v = [1.d0,0.d0,0.d0]
  p%q = 2
  p%m = 4.d0
end select

! 6. Set up one pushers to be used for all groups.
pushers = [pusher_container(pusher_boris(fixed_timestep=0.1d0))]

! 7. Set an event to stop the simulation.
events  = [event(stop_action(), start=1.d0), &
    event(write_action(), step=0.1d0)]

! 8. Run the main loop.
call main_loop(sim, pushers, events)

! 9. Finalize MPI
call MPI_Finalize(ierr)

contains
!> A field that is zero everywhere.
pure function E_zero(x, t) result(E)
  real*8, intent(in) :: x(3), t
  real*8 :: E(3) !< Electric field in V/M
  E = [0,0,0]
end function E_zero
!> A field that is one in the z-direction and zero in others.
pure function B_z(x, t) result(B)
  real*8, intent(in) :: x(3), t
  real*8 :: B(3) !< Magnetic field in Tesla
  B = [0,0,1]
end function B_z
end program ex2
