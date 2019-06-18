!>#Example 3
!> Postprocess existing data. Requires the output of `ex2`
!>
!>* geometry: cartesian
!>
!> Compile with: `make ex2 ex3`
!> Run with: `./ex2 && ./ex3 part*h5`
!> See the [annotated source](../sourcefile/ex3.f90.html) for details.
!>
!>## Description
!> This example outputs statistics on the kinetic energy of all particles
program ex3
use particle_tracer
use mod_diag_print_kinetic_energy
implicit none

! 1. Set up the simulation data type and actions
type(read_action)  :: reader
type(diag_print_kinetic_energy) :: print_ke
integer :: i

! 1. Set up MPI
call sim%initialize(num_groups=0,skip_jorek2help=.true.)
! 2. Loop over command arguments
do i=1,command_argument_count()
  ! 3. Read files from commandline
  call get_command_argument(i, reader%filename)
  call reader%run(sim)

  ! 4. Perform an analysis on this sim
  call print_ke%run(sim)
end do

! 5. Close the MPI subsystem
call sim%finalize
end program ex3
