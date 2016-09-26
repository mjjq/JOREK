!> This module contains testcases for the event system
module event_spec
use mod_event
use mod_particle_sim
use mod_main_loop
use mod_pusher
use mod_pusher_no_action
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

subroutine test_main_loop_with_five_increments
  type(event), dimension(:), allocatable :: events
  type(pusher_container), dimension(:), allocatable :: pushers
  type(particle_sim) :: sim
  pushers = [pusher_container(pusher_no_action())]
  events = [event(increment_action(), step=0.21d0, end=1.d0)]
  call main_loop(sim, pushers, events)
  select type (a => events(1)%action)
  type is (increment_action)
    call assert_equals(5, a%counter, "increment_action should run five times")
  end select
end subroutine test_main_loop_with_five_increments

subroutine test_main_loop_with_eleven_increments
  type(event), dimension(:), allocatable :: events
  type(pusher_container), dimension(:), allocatable :: pushers
  type(particle_sim) :: sim
  pushers = [pusher_container(pusher_no_action())]
  events = [event(increment_action(), step=0.1d0, end=1.d0)] ! this works or not depending on fp issues
  call main_loop(sim, pushers, events)
  select type (a => events(1)%action)
  type is (increment_action)
    call assert_equals(11, a%counter, "increment_action should run eleven times")
  end select
end subroutine test_main_loop_with_eleven_increments

subroutine test_main_loop_without_increments
  type(event), dimension(:), allocatable :: events
  type(pusher_container), dimension(:), allocatable :: pushers
  type(particle_sim) :: sim

  pushers = [pusher_container(pusher_no_action())]
  events = [event(increment_action(), start=1.33d0, end=1.d0)]
  call main_loop(sim, pushers, events)
  select type (a => events(1)%action)
  type is (increment_action)
    call assert_equals(0, a%counter, "increment_action should run zero times")
  end select
end subroutine test_main_loop_without_increments

subroutine do_increment_action(this, sim)
  use mod_particle_sim
  class(increment_action), intent(inout) :: this
  type(particle_sim), intent(inout) :: sim
  this%counter = this%counter + 1
  write(*,*) "incrementing counter to ", this%counter, ' at t=', sim%time
end subroutine do_increment_action
end module event_spec
