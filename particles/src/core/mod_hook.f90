!> Module to formalize performing an action on a single particle every (n) timestep(s)
!> Hooks are objects (derived types) containing information
!> about when to run. These are subclassed for the different particle pushers, with hooks
!> that are aware of the particle type.
module mod_hook
use mod_particle_types
implicit none
private
public particle_action, particle_action_noop, hook_base

!> Particle action abstract type, representing anything that can be done to a single particle
type, abstract :: particle_action
contains
  procedure(do), deferred, pass :: do
end type particle_action
interface
  pure subroutine do(this, particle)
    import :: particle_action, particle_base
    class(particle_action), intent(in) :: this
    class(particle_base), intent(inout) :: particle
  end subroutine do
end interface
type, extends(particle_action) :: particle_action_noop
  contains
    procedure :: do => do_nothing
end type

!> Hook type
type :: hook_base
  !not implemented yet: integer :: istep = 1 !< Every how many steps should this hook run? (takes precedence over step if > 0)
  !not implemented yet: real*8  :: step  = 0.d0 !< Every how often should this hook run? (rounded to the nearest istep > 1)

  !> Action to perform when this event runs
  class(particle_action), allocatable :: action
end type hook_base
interface hook_base
  module procedure new_hook_base
end interface hook_base

contains
function new_hook_base(action)
  type(hook_base) :: new_hook_base
  class(particle_action), intent(in) :: action
  allocate(new_hook_base%action, source=action)
end function new_hook_base

pure subroutine do_nothing(this, particle)
  class(particle_action_noop), intent(in) :: this
  class(particle_base), intent(inout) :: particle
end subroutine do_nothing
end module mod_hook
