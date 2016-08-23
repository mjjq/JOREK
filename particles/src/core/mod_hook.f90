!> Module to formalize performing an action on a single particle every (n) timestep(s)
!> Hooks are objects (derived types) containing information
!> about when to run, and an action to run.
module mod_hook
use mod_particle_group
use mod_particle_action
implicit none

!> Hook type
type, public :: hook
  integer :: istep !< Every how many steps should this hook run? (takes precedence over step)
  real*8  :: step  !< Every how often should this hook run? (rounded to the nearest istep > 1)

  !> Action to execute
  class(particle_action), allocatable :: action
end type hook
interface hook
  module procedure new_hook
end interface

private
contains
!> Constructor for an event
!> This is extensive to allow default values and overloading
function new_hook(act, istep, step)
  use mod_action
  use mpi
  type(hook) :: new_hook
  class(action), intent(in)     :: act
  integer, intent(in), optional :: istep !< Step every istep (takes precedence over step)
  real*8, intent(in), optional  :: step !< Step every step time (specify this or istep)
  integer :: ierr

  if (.not. present(istep) .and. .not. present(step)) then
    write(*,*) "Error: istep or step is required"
    call MPI_ABORT(MPI_COMM_WORLD, 0, ierr)
  else
    if (present(istep)) then
      new_hook%istep = istep
    else
      new_hook%step = step
    end if
  end if
  allocate(new_hook%action, source=act) ! because assignment is not yet supported in gfortran 6.1.1
end function new_hook
end module mod_hook
