!> Module to formalize performing an action every now and then
!> Events are objects (derived types) containing information
!> about when to run, and an action they should run.
module mod_event
use mod_action
use mod_constants
implicit none
private
public event

!> Event type
type :: event
  real*8  :: start    = 0.d0       !< Physical starting time
  real*8  :: step     = huge(0.d0) !< Step every how long?
  real*8  :: end      = huge(0.d0) !< Stop after time end

  integer, dimension(:), allocatable :: sync_groups !< which groups to require at a full-timestep (default = all, empty array = none)

  !> Action to perform when this event runs
  class(action), allocatable :: action
end type event
interface event
  module procedure new_event
end interface

contains
!> Constructor for an event
!> This is needed to allow changing default values
function new_event(act, start, step, end)
  type(event) :: new_event
  class(action), intent(in)     :: act
  real*8, intent(in), optional  :: start, step, end
  if (present(start))    new_event%start    = start
  if (present(step))     new_event%step     = step
  if (present(end))      new_event%end      = end
  allocate(new_event%action, source=act) ! because assignment is not yet supported in gfortran 6.1.1
end function new_event
end module mod_event
