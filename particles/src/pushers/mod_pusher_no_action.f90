!> Example particle pusher which does nothing, used for tests.
!> This module consists of a single type, [[pusher_no_action]], extending
!> [[pusher_base]]. Required elements to implement are the type, along with
!> any time-step control fields you would like to have (in addition to the
!> fixed_timestep field in pusher_base), and a subroutine implementing [[push_single]]
!> (in this case [[push_no_action]].)
!>
!> You also need to create an interface for the object, allocating groups
!> and hooks if present (and your custom fields). This is [[new_pusher_no_action]].
module mod_pusher_no_action
use mod_pusher
implicit none

!> Example pusher doing nothing. Must implement push_single according to the
!> interface in [[mod_pusher]].
type, extends(pusher_base) :: pusher_no_action
contains
  procedure :: push_single => push_no_action
end type pusher_no_action
!> Interface for creation of new `pusher_no_action`. This allows initialisation
!> with or without allocating components.
interface pusher_no_action
  module procedure new_pusher_no_action
end interface pusher_no_action

contains
!> The actual implementation of the pusher. This does not alter the particles,
!> except for the action of any hooks set.
!> There should be a small speed benefit of being pure, but it is not
!> yet required by the interface.
!>
!> It must push a single particle from time_start to time_end (or as close as possible
!> in the case of fixed_timesteps), and will be called from within an openmp block.
pure subroutine push_no_action(this, fields, particle, time_start, time_end)
  use mod_fields
  use mod_particle_types
  class(pusher_no_action), intent(in) :: this !< The current pusher
  class(fields_base), intent(in)      :: fields !< Fields used in the simulation, such as [[mod_prescribed_fields]]
  class(particle_base), intent(inout) :: particle !< The particle to push now
  real*8, intent(in) :: time_start, time_end !< How long to push for

  integer :: n, i, k

  n = 1
  if (allocated(this%fixed_timestep)) n = nint((time_end-time_start)/this%fixed_timestep)
  do i=1,n
    ! Run hooks if there are any
    if (allocated(this%hooks)) then
      do k=1,size(this%hooks,1)
        ! Test if there is any action allocated
        call this%hooks(k)%action%do(particle)
      end do
    end if
  end do
end subroutine push_no_action

!> This constructor function is provided as a
!> workaround for https://gcc.gnu.org/bugzilla/show_bug.cgi?id=77412
!> and to allow optional arguments.
pure function new_pusher_no_action(groups, fixed_timestep, hooks)
  type(pusher_no_action) :: new_pusher_no_action !< The constructor result
  integer, dimension(:), intent(in), optional         :: groups !< Which groups to push (such as [1,3] or unallocated to push all
  !< remaining groups)
  real*8, intent(in), optional                        :: fixed_timestep
  type(hook_base), intent(in), optional, dimension(:) :: hooks !< Hooks to add to the pusher
  if (present(groups)) allocate(new_pusher_no_action%groups, source=groups)
  if (present(fixed_timestep)) allocate(new_pusher_no_action%fixed_timestep, source=fixed_timestep)
  if (present(hooks)) allocate(new_pusher_no_action%hooks, source=hooks)
end function new_pusher_no_action
end module mod_pusher_no_action
