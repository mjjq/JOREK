!> Module containing a datatype for simulation parameters
module mod_particle_sim
use mod_particle_types
use mod_constants, only: CYLINDRICAL
implicit none
private
public particle_group, particle_sim

!> A group of particles, implemented as an allocatable array
type :: particle_group
  class(particle_base), dimension(:), allocatable :: particles
end type particle_group

!> Particle simulation type, containing all variables pertaining to a simulation.
type :: particle_sim
  real*8                                          :: time = 0.d0 !< time of the simulation. Only accurate when in events with sync or at
  !< the start of the simulation
  integer*1                                       :: geometry = CYLINDRICAL
  logical                                         :: stop_now = .false.
  type(particle_group), dimension(:), allocatable :: groups
contains
  procedure :: finalize
  procedure :: initialize
end type particle_sim

contains
!> Actions to perform when setting up a simulation
subroutine initialize(sim, num_groups)
  use mpi
  class(particle_sim), intent(inout) :: sim !< why is this class() and not type()?
  integer, intent(in) :: num_groups
  integer :: provided, ierr
  call MPI_Init_thread(MPI_THREAD_MULTIPLE, provided, ierr)
  if (ierr .ne. 0) write(*,*) "Error ", ierr, " in MPI_Init_thread"
  if (provided .ne. MPI_THREAD_MULTIPLE) write(*,*) "WARNING: provided(", provided, ") != MPI_THREAD_MULTIPLE"

  allocate(sim%groups(num_groups))
end subroutine
!> Actions to perform when stopping the simulation.
subroutine finalize(sim)
  class(particle_sim), intent(in) :: sim !< why is this class() and not type()?
  integer :: ierr
  if (sim%stop_now) then
    write(*,"(A,g12.6,A)") "INFO: Stop requested at ", sim%time, " , exiting"
  else
    write(*,"(A,g12.6,A)") "INFO: End of events at ", sim%time, " , exiting"
  end if
  call MPI_Finalize(ierr)
end subroutine
end module mod_particle_sim
