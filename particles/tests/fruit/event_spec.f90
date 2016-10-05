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

subroutine test_event_run_at
  use mod_constants, only: tick
  type(event) :: e
  e = event(stop_action(), start=1.d0)
  call assert_false(e%run_at(0.d0), 'must not run before start')
  call assert_true(e%run_at(1.d0), 'must run at start')
  call assert_false(e%run_at(2.d0), 'must not run after start if no step is set')
  e = event(stop_action(), start=1.d0, step=2.d0)
  call assert_true(e%run_at(3.d0), 'must run at start + 1 step')
  call assert_false(e%run_at(3.1d0), 'must not run at start + 1.05 step')
  call assert_false(e%run_at(3.0001d0), 'must not run at start + 1.00005 step')
  call assert_false(e%run_at(-1.d0), 'must not run at start - 1 step')
  e = event(stop_action(), start=1.d0, step=2.d0, end=3.d0)
  call assert_true(e%run_at(3.d0), 'must run at start + 1 step == end')
  call assert_true(e%run_at(3.d0+tick/2), 'must run at start + 1 step == end')
  call assert_false(e%run_at(5.d0), 'must not run > end')
  e = event(stop_action(), start=1.d0, end=2.d0)
  call assert_true(e%run_at(1.d0), 'must run at start')
  call assert_false(e%run_at(2.d0), 'must not run after start if no step is set')
end subroutine test_event_run_at
end module event_spec
