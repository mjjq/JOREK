!>#Example 3
!> Postprocess existing data. Requires the output of `ex2`
!>
!>* geometry: cartesian
!>
!> Compile with: `make ex2 ex3`
!> Run with: `./ex2 && ./ex3`
!> See the [annotated source](../sourcefile/ex3.f90.html) for details.
!>
!>## Description
!> This example outputs statistics on the kinetic energy of all particles
program ex3
use particle_tracer
use mpi
use mod_diag_print_kinetic_energy
implicit none

! 1. Set up the simulation data type and actions
type(particle_sim) :: sim
type(read_action)  :: reader
type(diag_print_kinetic_energy) :: print_ke
integer :: i, ierr

! 2. Set up MPI IO, needed for reading/writing
call MPI_Init(ierr)

! 3. Set up the fields to be used in the simulation. (E and B defined below)
!    This is mostly needed because it sets the geometry of the fields, which the event might use.
allocate(sim%fields, source=prescribed_fields(CARTESIAN, E_zero, B_z))

do i=1,command_argument_count()
  ! 4. Read files from commandline
  call get_command_argument(i, reader%filename)
  call reader%run(sim)

  ! 5. Perform an analysis on this sim
  call print_ke%run(sim)
end do
  
! 6. Close the MPI subsystem
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
end program ex3
