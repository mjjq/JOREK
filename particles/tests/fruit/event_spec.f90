!> This module contains testcases for the event system
module event_spec
use mod_event
use mod_action
use fruit
implicit none

contains

subroutine test_create_event_with_stop_action
  type(event), dimension(:), allocatable :: events
  events = [event(stop_action(), start=1.d0)]
  call assert_equals(events(1)%start, 1.d0)
end subroutine test_create_event_with_stop_action
end module event_spec
