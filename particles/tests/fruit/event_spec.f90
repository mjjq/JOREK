!> This module contains testcases for the event system
module event_spec
use mod_event, only: event
use mod_action, only: action, stop_action
use fruit
implicit none

contains

subroutine test_create_event_with_stop_action
  type(event), dimension(:), allocatable :: events
  class(action), allocatable :: a

  allocate(a, source=stop_action())

  events = [ &
    event(stop_action(), start=1.d0) &
  ]

  call assert_equals(events(1)%start, 1.d0)
end subroutine test_create_event_with_stop_action
end module event_spec
