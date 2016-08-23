module mod_pusher_no_action
use mod_pusher, only: pusher_base
implicit none

type, extends(pusher_base) :: pusher_no_action
contains
  procedure :: push_single => push_no_action
end type pusher_no_action

contains
pure subroutine push_no_action(this, fields, particle, time_start, time_end)
  use mod_fields
  use mod_particle_base
  class(pusher_no_action), intent(inout)   :: this
  class(fields_base), intent(in)      :: fields
  class(particle_base), intent(inout) :: particle
  real*8, intent(in) :: time_start, time_end
end subroutine push_no_action
end module mod_pusher_no_action
