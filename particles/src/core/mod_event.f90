!> Module to formalize performing an action every now and then
!> Events are objects (derived types) containing information
!> about when to run, and an action they should run.
module mod_event
use mod_action
use mod_particle_sim
use mod_constants, only: tick
implicit none
private
public event, with, next_event_at, check_and_fix_timesteps

!> Event type
type :: event
  real*8  :: start    = 0.d0       !< Physical starting time
  real*8  :: step     = huge(0.d0) !< Step every how long?
  real*8  :: end      = huge(0.d0) !< Stop after time end

  integer, dimension(:), allocatable :: sync_groups !< which groups to require at a full-timestep (default = all, empty array = none)

  !> Action to perform when this event runs
  class(action), allocatable :: action
end type event
interface event
  module procedure new_event
end interface

interface with
  module procedure with_event_0D, with_action_0D, &
        with_event_1D, with_event_1D_mask, &
        with_action_1D, with_action_1D_mask
  end interface
contains
!> Constructor for an event
!> This is needed to allow changing default values
function new_event(act, start, step, end)
  type(event) :: new_event
  class(action), intent(in)     :: act
  real*8, intent(in), optional  :: start, step, end
  if (present(start))    new_event%start    = start
  if (present(step))     new_event%step     = step
  if (present(end))      new_event%end      = end
  allocate(new_event%action, source=act) ! because assignment is not yet supported in gfortran 6.1.1
end function new_event


subroutine with_event_0D(sim, single_event)
  type(particle_sim), intent(inout) :: sim
  type(event), intent(inout) :: single_event
  call single_event%action%run(sim)
end subroutine with_event_0D
subroutine with_action_0D(sim, single_action)
  type(particle_sim), intent(inout) :: sim
  class(action), intent(inout) :: single_action
  call single_action%run(sim)
end subroutine with_action_0D
subroutine with_event_1D(sim, events)
  type(particle_sim), intent(inout) :: sim
  type(event), intent(inout), dimension(:) :: events
  integer :: i
  do i=1,size(events)
    call events(i)%action%run(sim)
  end do
end subroutine with_event_1D
subroutine with_action_1D(sim, actions)
  type(particle_sim), intent(inout) :: sim
  class(action), intent(inout), dimension(:) :: actions
  integer :: i
  do i=1,size(actions)
    call actions(i)%run(sim)
  end do
end subroutine with_action_1D
subroutine with_event_1D_mask(sim, events, mask)
  type(particle_sim), intent(inout) :: sim
  type(event), intent(inout), dimension(:) :: events
  !logical, dimension(size(events,1)), intent(in) :: mask ! internal compiler error in gfortran, workaround below
  logical, dimension(:), intent(in) :: mask
  integer :: i
  do i=1,size(events)
    if (mask(i)) call events(i)%action%run(sim)
  end do
end subroutine with_event_1D_mask
subroutine with_action_1D_mask(sim, actions, mask)
  type(particle_sim), intent(inout) :: sim
  class(action), intent(inout), dimension(:) :: actions
  !logical, dimension(size(actions,1)), intent(in) :: mask
  logical, dimension(:), intent(in) :: mask
  integer :: i
  do i=1,size(actions)
    if (mask(i)) call actions(i)%run(sim)
  end do
end subroutine with_action_1D_mask

!> Return the number of the next event(s) to run
!> If the time is > event_start, calculate the time from
!> ```
!>  <-DT->               dt
!> |------|------|------|--|---|
!> |T0                     t   te
!> ```
!> where DT is the event%step, T0=event%start,
!> dt = mod(t-T0, DT) and te = t + DT - dt
!>
!> Any events that are within 1d-14 of the current time will not run
!> (to prevent double events due to floating-point issues)
subroutine next_event_at(events, current_time, run_event, event_time)
  type(event), intent(inout), dimension(:) :: events
  real*8, intent(in) :: current_time
  logical, dimension(size(events)), intent(out) :: run_event
  real*8, intent(out) :: event_time
  real*8 :: event_run !< when events(i) is to run (at the soonest)
  integer :: i

  event_time = huge(0.d0)
  run_event(:) = .false.
  do i=1,size(events)
    ! next event needs to be at least tick in the future
    if (events(i)%start .gt. current_time + tick) then
      event_run = events(i)%start
    else
      event_run = current_time + events(i)%step - mod(current_time - events(i)%start, events(i)%step)
      if (abs(event_run - current_time) .le. tick) event_run = current_time + events(i)%step
    end if

    ! If this event has ended already
    if (event_run .gt. events(i)%end + tick) cycle

    ! if this event occurs faster than the previously fastest (event_time)
    if (event_run .lt. event_time - tick) then
      run_event(:) = .false.
      run_event(i) = .true.
      event_time = event_run
    else if (event_run .le. event_time + tick) then ! if it is equally fast
      run_event(i) = .true.
    end if
  end do
end subroutine next_event_at


!> Calculate whether we need to change any of the fixed timesteps or events to match
!> For each of the pushers with a fixed timestep
subroutine check_and_fix_timesteps(pusher_timestep, events)
  use mod_event_timestep
  real*8,      intent(inout), dimension(:) :: pusher_timestep
  type(event), intent(inout), dimension(:) :: events

  real*8, dimension(:), allocatable    :: event_start, event_step, pusher_timestep_work
  logical, dimension(:,:), allocatable :: constraints
  integer :: i, j, pusher, ierr

  ! select start and step times of all events
  event_start = [(events(i)%start, i=1, size(events))]
  event_step  = [(events(i)%step,  i=1, size(events))]
  pusher_timestep_work = pusher_timestep

  ! find constraints
  allocate(constraints(size(events),size(pusher_timestep)))
  constraints = .false.
  do i=1,size(events)
    if (allocated(events(i)%sync_groups)) then ! if we have specified specific groups to sync only
      ! for each of the groups, check if it needs to sync to this event, and set it for that pusher
      do j=1,size(events(i)%sync_groups)
        if ((events(i)%sync_groups(j) .lt. 1) .or. &
            (events(i)%sync_groups(j) .gt. size(pusher_timestep))) then
          write(*,*) "ERROR: cannot find group ", events(i)%sync_groups(j)
        else
          constraints(i,events(i)%sync_groups(j)) = .true.
        end if
      end do
    else ! sync all groups
      constraints(i,:) = .true. ! add all the pushers for this event
    end if
  end do

  call fix_event_timestep(pusher_timestep_work, event_start, event_step, constraints, ierr)
  if (ierr .ne. 0) then
    write(*,*) "error in fix_event_timestep"
    call exit(1)
  end if

  ! show changes in timesteps and set them
  do i=1,size(pusher_timestep)
    if (abs(pusher_timestep_work(i) - pusher_timestep(i)) .gt. TICK) then
      write(*,*) "INFO: changing timestep of pusher", i, " from ", pusher_timestep(i), " to ", pusher_timestep_work(i)
    end if
    ! always update, but notify only for significant changes
    pusher_timestep(i) = pusher_timestep_work(i)
  end do
  do i=1,size(events)
    if (abs(event_start(i) - events(i)%start) .gt. TICK) then
      write(*,"(A,i3,A,A,A,g14.8,A,g14.8)") "INFO: changing start time of event ", i, &
          "(", trim(events(i)%action%name), ") from ", events(i)%start, " to ", event_start(i)
    end if
    if (abs(event_step(i) - events(i)%step) .gt. TICK) then
      write(*,"(A,i3,A,A,A,g14.8,A,g14.8)") "INFO: changing timestep of event ", i, &
          "(", trim(events(i)%action%name), ") from ", events(i)%step, " to ", event_step(i)
    end if
    ! always update, notify only for significant changes
    events(i)%start = event_start(i)
    events(i)%step = event_step(i)
  end do
end subroutine check_and_fix_timesteps
end module mod_event
