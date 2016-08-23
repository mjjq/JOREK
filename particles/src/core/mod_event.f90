!> Module to formalize performing an action every now and then
!> Events are objects (derived types) containing information
!> about when to run, and an action they should run.
module mod_event
use mod_action
implicit none

!> Event type
type, public :: event
  real*8  :: start    = 0.d0       !< Physical starting time
  real*8  :: step     = huge(0.d0) !< Step every how long?
  real*8  :: end      = huge(0.d0) !< Stop at time end
  logical :: sync     = .false.    !< Should all groups be at a full-timestep? (if not, they may be slightly before or after this time
  logical :: mpi_sync = .false.    !< Should all processes do this event simultaneously?

  !> Action
  class(action), allocatable :: action
end type event
interface event
  module procedure new_event
end interface

contains
!> Constructor for an event
!> This is needed to allow changing default values
function new_event(act, start, step, end, sync, mpi_sync)
  type(event) :: new_event
  class(action), intent(in)     :: act
  real*8, intent(in), optional  :: start, step, end
  logical, intent(in), optional :: sync, mpi_sync
  if (present(start))    new_event%start    = start
  if (present(step))     new_event%step     = step
  if (present(end))      new_event%end      = end
  if (present(sync))     new_event%sync     = sync
  if (present(mpi_sync)) new_event%mpi_sync = mpi_sync
  allocate(new_event%action, source=act) ! because assignment is not yet supported in gfortran 6.1.1
end function new_event
end module mod_event
