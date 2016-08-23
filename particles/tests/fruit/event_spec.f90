!> This module contains testcases for the event system
module event_spec
use mod_event, only: event
use mod_action, only: action, stop_action
use mod_particle_sim, only: particle_sim
use mod_main_loop, only: main_loop
use mod_pusher, only: pusher_base, pusher_container
use mod_pusher_no_action, only: pusher_no_action
use fruit
implicit none

type, extends(action) :: increment_action
  integer :: counter = 0
contains
  procedure :: do => do_increment_action
end type increment_action

contains

subroutine test_create_event_with_stop_action
  type(event), dimension(:), allocatable :: events
  events = [event(stop_action(), start=1.d0)]
  call assert_equals(events(1)%start, 1.d0)
end subroutine test_create_event_with_stop_action


subroutine test_main_loop_with_three_increments
  type(event), dimension(:), allocatable :: events
  type(pusher_container), dimension(:), allocatable :: pushers
  type(particle_sim) :: sim
  type(increment_action) :: inc
  pushers = [pusher_container(pusher_no_action())]
  events = [ &
    event(inc, step=0.33d0), &
    event(stop_action(), start=1.d0) &
  ]
  call main_loop(sim, pushers, events)
  call assert_equals(3, inc%counter, "increment_action should run three times")
end subroutine test_main_loop_with_three_increments

subroutine test_main_loop_with_no_increments
  type(event), dimension(:), allocatable :: events
  type(pusher_container), dimension(:), allocatable :: pushers
  type(particle_sim) :: sim
  type(increment_action) :: inc

  pushers = [pusher_container(pusher_no_action())]
  events = [ &
    event(inc, step=1.33d0), &
    event(stop_action(), start=1.d0) &
  ]
  call main_loop(sim, pushers, events)
  call assert_equals(0, inc%counter, "increment_action should run zero times")
end subroutine test_main_loop_with_no_increments

subroutine do_increment_action(this, sim)
  use mod_particle_sim
  class(increment_action), intent(inout) :: this
  type(particle_sim), intent(inout) :: sim
  this%counter = this%counter + 1
end subroutine do_increment_action
end module event_spec
