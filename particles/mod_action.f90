!> Module to formalize performing an action every now and then
module mod_action
use mod_particle_sim
implicit none
private
public action, stop_action

!> Action abstract type, representing anything that can be done to a simulation
type, abstract :: action
  !> Logging variable, set this in an initializer
  character(len=30) :: name = "unset action" !< Event name for logging
  logical :: log = .false. !< Output event duration

  !> Timing variables
  real*8, private :: t0 = 0.d0, w0 = 0.d0
contains
  procedure, pass, public :: run
  procedure(do_interface), deferred, pass, private :: do
end type action
interface
  subroutine do_interface(this, sim)
    import :: action, particle_sim
    class(action), intent(inout)      :: this
    type(particle_sim), intent(inout) :: sim
  end subroutine do_interface
end interface


!> Example action (stops the simulation)
type, extends(action) :: stop_action
contains
  procedure :: do => do_stop_action
end type stop_action
interface stop_action
  module procedure new_stop_action
end interface

contains
!> Constructor for stop_action
function new_stop_action()
  type(stop_action) :: new_stop_action
  new_stop_action%name = "Stop"
  new_stop_action%log  = .false.
end function new_stop_action

!> Perform the stop action
subroutine do_stop_action(this, sim)
  class(stop_action), intent(inout) :: this
  type(particle_sim), intent(inout) :: sim
  sim%stop_now = .true.
end subroutine do_stop_action




!> Run an action
subroutine run(this, sim)
  !$ use omp_lib
  class(action), intent(inout)      :: this
  type(particle_sim), intent(inout) :: sim
  real*8 :: t1, w1
  logical :: has_omp
  has_omp = .false.
  !$ has_omp = .true.

  call cpu_time(this%t0)
  !$ this%w0 = omp_get_wtime()
  call this%do(sim)
  call cpu_time(t1)
  !$ w1 = omp_get_wtime()

  ! this is only on node 0 of course
  if (this%log) then
    if (.not. has_omp) write(*,"(A,A,f7.4,A)") trim(this%name), " finished in ", t1-this%t0, "s"
    !$ write(*,"(A,A,f7.4,A,f7.4,A)") trim(this%name), " finished in ", w1-this%w0, &
    !$ "s (cpu time: ", t1-this%t0, ")"
  end if
end subroutine run
end module mod_action
