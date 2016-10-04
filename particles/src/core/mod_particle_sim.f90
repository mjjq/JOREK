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

  real*8, private :: time !< time at which this group is currently
  real*8, private :: t0, t1
contains
  procedure :: before_push
  procedure :: after_push
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
  if (sim%stop_now) then
    write(*,"(A,g12.6,A)") "INFO: Stop requested at ", sim%time, " , exiting"
  else
    write(*,"(A,g12.6,A)") "INFO: End of events at ", sim%time, " , exiting"
  end if
end subroutine


!> Calculate push time 
subroutine before_push(this)
  Class(particle_group), intent(inout) :: this
  call cpu_time(this%t0)
end subroutine before_push

!> Calculate and print push time, with MPI support
subroutine after_push(this, i, time)
  class(particle_group), intent(inout) :: this
  real*8, intent(in) :: time
  integer, intent(in) :: i
  call cpu_time(this%t1)
  write(*,"(A,i3,A,g12.6,A,g12.6,A,g12.6,A)") "INFO: Group ", i, " pushed from ", this%time, " to ", time, &
      " in ", this%t1-this%t0, " s"
  this%time = time
end subroutine after_push
end module mod_particle_sim
