!>#Example 5
!> Using pusher hooks
!>
!>* fields: prescribed
!>* pusher: boris
!>* geometry: cartesian
!>
!> Compile with: `make ex5`
!> Run with: `./ex5`
!> See the [annotated source](../sourcefile/ex5.f90.html) for details.
!>
!>## Description
!> This example outputs statistics on the kinetic energy of all particles.
!> A custom hook has been implemented, to reduce the particle speed by a fraction dt/1s every timestep.
program ex5
use particle_tracer
use mpi
use mod_diag_print_kinetic_energy
implicit none

! 1. Set up the simulation data type and actions
type(particle_sim) :: sim
type(pusher_container), dimension(:), allocatable :: pushers
type(event), dimension(:), allocatable :: events
type(diag_print_kinetic_energy) :: print_ke
integer :: i, ierr

! 2. Set up MPI IO, needed for reading/writing
call MPI_Init(ierr)

! 3. Set up the fields to be used in the simulation. (E and B defined below)
!    This is mostly needed because it sets the geometry of the fields, which the event might use.
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
!    Here we add a hook that does nothing with the particle
pushers = [pusher_container(pusher_boris(fixed_timestep=0.1d0, &
    hooks=[hook_base(particle_action_noop())] &
    ))]

! 7. Set an event to stop the simulation.
events  = [event(stop_action(), start=1.d0), &
           event(print_ke, step=0.1d0)]

! 8. Run the main loop.
call main_loop(sim, pushers, events)
  
! 9. Close the MPI subsystem
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
end program ex5
