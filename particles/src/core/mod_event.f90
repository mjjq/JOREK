!> Module to formalize performing an action every now and then
!> Events are objects (derived types) containing information
!> about when to run, and an action they should run.
module mod_event
use mod_action
implicit none

!> Event type
type, public :: event
  real*8  :: start !< Physical starting time (default 0)
  real*8  :: step  !< Step every how long?
  real*8  :: end   !< Stop at time tend
  logical :: sync  !< Should group(s) be at a full-timestep?

  !> Action
  class(action), allocatable :: action
end type event
interface event
  module procedure new_event
end interface

private
contains
!> Constructor for an event
!> This is so extensive to allow default values and overloading
function new_event(act, start, step, end, sync)
  use mod_action
  type(event) :: new_event
  class(action), intent(in)     :: act
  real*8, intent(in), optional  :: start, step, end
  logical, intent(in), optional :: sync
  new_event%start = 0.d0
  if (present(start)) new_event%start = start
  new_event%step = huge(0.d0)
  if (present(step)) new_event%step = step
  new_event%end = huge(0.d0)
  if (present(end)) new_event%end = end
  new_event%sync = .true.
  if (present(sync)) new_event%sync = sync
  allocate(new_event%action, source=act) ! because assignment is not yet supported in gfortran 6.1.1
end function new_event
end module mod_event
