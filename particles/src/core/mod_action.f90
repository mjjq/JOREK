!> Module to formalize performing an action every now and then
module mod_action
use mod_particle_sim
implicit none

!> Action abstract type, representing anything that can be done to a simulation
type, abstract, public :: action
  !> Logging variable, set this in an initializer
  character(len=30) :: name = "unset action" !< Event name for logging
  logical :: log = .false. !< Output event duration

  !> Timing variables
  real*8, private :: t0 = 0.0
contains
  procedure, pass, public :: run
  procedure(do), deferred, pass, private :: do
end type action
interface
  subroutine do(this, sim)
    import :: action, particle_sim
    class(action), intent(inout)      :: this
    type(particle_sim), intent(inout) :: sim
  end subroutine do
end interface


!> Example action (stops the simulation using MPI_ABORT)
type, extends(action), public :: stop_action
contains
  procedure :: do => do_stop_action
end type stop_action
interface stop_action
  module procedure new_stop_action
end interface

private
contains
!> Constructor for stop_action
function new_stop_action()
  type(stop_action) :: new_stop_action
  new_stop_action%name = "Stop"
  new_stop_action%log  = .false.
end function new_stop_action

!> Perform the stop action
subroutine do_stop_action(this, sim)
  use mpi
  class(stop_action), intent(inout) :: this
  type(particle_sim), intent(inout) :: sim
  integer :: ierr
  call exit(0)
end subroutine do_stop_action




!> Run an action
subroutine run(this, sim)
  class(action), intent(inout)      :: this
  type(particle_sim), intent(inout) :: sim
  real*8 :: t1

  call cpu_time(this%t0)
  call this%do(sim)
  call cpu_time(t1)
  if (this%log) write(*,"(A,A,f7.4,A)") trim(this%name), " finished in ", t1-this%t0, "s"
end subroutine run
end module mod_action
