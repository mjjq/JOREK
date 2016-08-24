!> The main loop calculates the required timesteps to pass all events
!> (exactly or not, depending on the value of sync for the events)
module mod_main_loop
use mod_particle_sim, only: particle_sim
use mod_event, only: event
use mod_pusher, only: pusher_container
implicit none

public :: main_loop
private
contains
subroutine main_loop(sim, pushers, events)
  type(particle_sim), intent(inout)         :: sim
  type(pusher_container), intent(inout), dimension(:) :: pushers
  type(event), intent(inout), dimension(:)  :: events
  integer :: ierr, i
  integer, dimension(:), allocatable :: next_events
  real*8 :: time, next_event_time
  
  ! Check if pushers contain enough for all groups
  call check_pusher_groups(pushers, size(sim%groups), ierr)
  if (ierr .ne. 0) return

  ! Calculate required timesteps
  call calculate_timesteps(pushers, events, ierr)
  if (ierr .ne. 0) return

  time = 0.d0
  do
    ! Calculate which of the events is next, or stop if there are no more
    call next_event_index(events, time, next_events, next_event_time)
    if (size(next_events) .eq. 0) then
      write(*,*) "INFO: end of events, exiting"
      exit ! stop this loop
    end if

    ! push particles until that time
    !$omp parallel default(none) &
    !$omp shared(sim, pushers, events) &
    !$omp private(ierr, next_event_time, next_events, time)
    ! TODO
    !$omp end parallel

    time = next_event_time
    ! run event(s) on the master thread
    do i=1,size(next_events)
      call events(next_events(i))%action%run(sim)
    end do
  end do
end subroutine

!> Return the number of the next event(s) to run
!> If the time is > event_start, calculate the time from
!> ```
!> |--DT--|              dt
!> |------|------|------|--|---|
!> |T0                     t   te
!> ```
!> where DT is the event%step, T0=event%start,
!> dt = mod(t-T0, DT) and te = t + DT - dt
subroutine next_event_index(events, current_time, next_events, event_time)
  type(event), intent(inout), dimension(:) :: events
  real*8, intent(in) :: current_time
  integer, dimension(:), allocatable, intent(out) :: next_events
  real*8, intent(out) :: event_time
  real*8 :: event_run
  logical, dimension(size(events)) :: event_first
  integer :: i

  event_time = huge(0.d0)
  event_first(:) = .false.
  do i=1,size(events)
    if (events(i)%start .gt. current_time) then
      event_run = events(i)%start
    else
      event_run = current_time + events(i)%step - mod(current_time - events(i)%start, events(i)%step)
      ! Do not run any events that are really close to now (to fix floating-point issues)
      if (event_run - current_time .le. 1d-14) cycle
    end if
    if (event_run .le. event_time .and. event_run .le. events(i)%end) then
      event_first(i) = .true.
      event_time = event_run
    end if
  end do

  next_events = pack([(i, i=1, size(events))], event_first)
end subroutine next_event_index

!> Check whether there is exactly one pusher for each group
!> At most 1 pusher can have the groups array unallocated, and it will be used for all
!> unlisted groups
subroutine check_pusher_groups(pushers, num_groups, ierr)
  use mpi
  type(pusher_container), intent(inout), dimension(:) :: pushers !< all the requested pushers
  integer, intent(in) :: num_groups !< number of groups in total (1..num_groups)
  integer, intent(out) :: ierr !< if nonzero we cannot run the simulation with this config
  integer :: num_unallocated, index_unallocated, i, j
  logical, dimension(num_groups) :: group_pushed

  ierr = 0
  num_unallocated = 0
  index_unallocated = 0
  group_pushed = .false.
  do i=1,size(pushers)
    if (.not. allocated(pushers(i)%pusher)) then
      write(*,*) "ERROR: unallocated pusher"
      ierr = 1
      return
    else
      if (.not. allocated(pushers(i)%pusher%groups)) then
        num_unallocated = num_unallocated + 1
        index_unallocated = i
      else
        do j=1,size(pushers(i)%pusher%groups)
          if (group_pushed(j)) then
            write(*,*) "ERROR: multiple pushers for group", j
            ierr = 3
            return
          else
            group_pushed(j) = .true.
          end if
        end do
      end if
    end if
  end do
  if (num_unallocated > 1) then
    write(*,*) "ERROR: too many unallocated groups in pusher list"
    ierr = 2
    return
  end if

  ! set the unallocated pusher to do all missing groups (select from implied do-loop from 1..num_groups)
  if (index_unallocated .gt. 0) then
    pushers(index_unallocated)%pusher%groups = pack([(i, i=1, num_groups)], .not. group_pushed)
  else
    if (.not. all(group_pushed)) then
      write(*,*) "ERROR: not all groups have pushers. Missing: ", pack([(i, i=1, num_groups)], .not. group_pushed)
      ierr = 4
      return
    end if
  end if
end subroutine check_pusher_groups



!> Calculate whether we need to change any of the fixed timesteps or events to match
subroutine calculate_timesteps(pushers, events, ierr)
  type(pusher_container), intent(inout), dimension(:) :: pushers !< all the requested pushers
  type(event), intent(inout), dimension(:)  :: events
  integer, intent(out) :: ierr !< if nonzero we cannot run the simulation with this config

  ! Test whether step > 0
end subroutine calculate_timesteps
end module mod_main_loop
