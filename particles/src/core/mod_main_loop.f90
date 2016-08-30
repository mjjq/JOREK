!> The main loop calculates the required timesteps to pass all events
!> (exactly or not, depending on the value of sync for the events)
module mod_main_loop
use mod_particle_sim, only: particle_sim
use mod_event, only: event
use mod_pusher, only: pusher_container
use mod_constants, only: tick
implicit none

public :: main_loop
private
contains
subroutine main_loop(sim, pushers, events)
  type(particle_sim), intent(inout)         :: sim
  type(pusher_container), intent(inout), dimension(:) :: pushers
  type(event), intent(inout), dimension(:)  :: events
  integer :: ierr, i, j
  integer, dimension(:), allocatable :: next_events
  real*8 :: next_event_time

  
  ! Check if pushers contain enough for all groups
  call check_set_pusher_groups(sim, pushers, ierr)
  if (ierr .ne. 0) return

  ! Check if we have any fixed_timestep pushers and reschedule events to fit
  call check_and_fix_timesteps(sim, pushers, events, ierr)
  if (ierr .ne. 0) return

  ! if we are at 0.d0 (to within one tick) run all of the start-events
  if (abs(sim%time) .lt. tick) then
    call next_event_index(events, sim%time, next_events, next_event_time, include_now=.true.)
    do i=1,size(next_events)
      call events(next_events(i))%action%run(sim)
    end do
  end if

  do
    ! Calculate which of the events is next, or stop if there are no more
    call next_event_index(events, sim%time, next_events, next_event_time)
    if (size(next_events) .eq. 0) then
      write(*,*) "INFO: end of events, exiting"
      exit ! stop this loop
    end if

    !$omp parallel default(none) &
    !$omp shared(sim, pushers, events, next_event_time) &
    !$omp private(ierr, i, j)
    ! note that this is not a parallel do loop, just the start of this parallel region
    do i=1,size(sim%groups)
      !$omp do
      do j=1,size(sim%groups(i)%particles)
        ! pusher guarantees to push the particle until at least next_event_time (or a little bit further)
        call pushers(sim%groups(i)%pusher)%pusher%push_single(sim%fields, sim%groups(i)%particles(j), &
            sim%groups(i)%time, next_event_time)
      end do
      !$omp end do
      sim%groups%time = next_event_time
    end do
    !$omp end parallel

    sim%time = next_event_time
    ! run event(s) on the master thread
    do i=1,size(next_events)
      call events(next_events(i))%action%run(sim)
    end do
  end do
end subroutine


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
subroutine next_event_index(events, current_time, next_events, event_time, include_now)
  type(event), intent(inout), dimension(:) :: events
  real*8, intent(in) :: current_time
  integer, dimension(:), allocatable, intent(out) :: next_events
  real*8, intent(out) :: event_time
  logical, optional :: include_now !< if set to .true., do not remove events occurring at current_time (+- tolerance)
  real*8 :: event_run !< when events(i) is to run (at the soonest)
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
      if (.not. (present(include_now) .and. include_now) &
          .and. abs(event_run - current_time) .le. tick) event_run = current_time + events(i)%step
    end if

    ! If this event has ended already
    if (event_run .gt. events(i)%end + tick) cycle

    ! if this event occurs faster than the previously fastest (event_time)
    if (event_run .lt. event_time - tick) then
      event_first(:) = .false.
      event_first(i) = .true.
      event_time = event_run
    else if (event_run .le. event_time + tick) then
      event_first(i) = .true.
    end if
  end do

  ! Select all indices which have true in event_first
  next_events = pack([(i, i=1, size(events))], event_first)
end subroutine next_event_index


!> Check whether there is exactly one pusher for each group
!> At most 1 pusher can have the groups array unallocated, and it will be used for all
!> unlisted groups
subroutine check_set_pusher_groups(sim, pushers, ierr)
  type(particle_sim), intent(inout)                   :: sim
  type(pusher_container), intent(inout), dimension(:) :: pushers !< all the requested pushers
  integer, intent(out) :: ierr !< if nonzero we cannot run the simulation with this config
  integer :: num_unallocated, index_unallocated, i, j
  logical, dimension(size(sim%groups)) :: group_pushed

  if (.not. allocated(sim%groups)) then
    write(*,*) "WARNING: no groups allocated"
    allocate(sim%groups(1:0)) ! if we have no groups, allow this for testing purposes
  end if
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
            sim%groups%pusher = i
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
    pushers(index_unallocated)%pusher%groups = pack([(i, i=1, size(sim%groups))], .not. group_pushed)
    do i=1,size(sim%groups)
      if (.not. group_pushed(i)) then
        sim%groups(i)%pusher = index_unallocated
      end if
    end do
  else
    if (.not. all(group_pushed)) then
      write(*,*) "ERROR: not all groups have pushers. Missing: ", pack([(i, i=1, size(sim%groups))], .not. group_pushed)
      ierr = 4
      return
    end if
  end if
end subroutine check_set_pusher_groups


!> Calculate whether we need to change any of the fixed timesteps or events to match
!> For each of the pushers with a fixed timestep 
subroutine check_and_fix_timesteps(sim, pushers, events, ierr)
  use mod_event_timestep
  type(particle_sim),     intent(in)                  :: sim
  type(pusher_container), intent(inout), dimension(:) :: pushers !< all the requested pushers
  type(event),            intent(inout), dimension(:) :: events
  integer, intent(out) :: ierr !< if nonzero we cannot run the simulation with this config

  logical, dimension(size(pushers))    :: pusher_timestep_fixed
  real*8, dimension(:), allocatable    :: pusher_timestep
  real*8, dimension(:), allocatable    :: event_start, event_step
  logical, dimension(:,:), allocatable :: constraints
  integer :: i, j, pusher

  ierr = 0

  ! check if we have any pushers with a fixed timestep
  pusher_timestep_fixed = .false.
  do i=1,size(pushers)
    ! whether pushers(i)%pusher is allocated was tested before, in [[check_pusher_groups]]
    if (allocated(pushers(i)%pusher%fixed_timestep)) pusher_timestep_fixed(i) = .true.
  end do

  ! if so, create constraints and call fix_event_timestep
  if (any(pusher_timestep_fixed)) then
    ! select fixed_timesteps of all pushers
    pusher_timestep = pack([(pushers(i)%pusher%fixed_timestep, i=1, size(pushers))], pusher_timestep_fixed)
    ! select start and step times of all events
    event_start = [(events(i)%start, i=1, size(events))]
    event_step  = [(events(i)%step,  i=1, size(events))]

    ! find constraints
    allocate(constraints(size(events),size(pusher_timestep)))
    constraints = .false.
    do i=1,size(events)
      if (allocated(events(i)%sync_groups)) then ! if we have specified specific groups to sync only
        ! for each of the groups, check if it needs to sync to this event, and set it for that pusher
        do j=1,size(events(i)%sync_groups)
          if ((events(i)%sync_groups(j) .lt. lbound(sim%groups,1)) .or. &
              (events(i)%sync_groups(j) .gt. ubound(sim%groups,1))) then
            write(*,*) "ERROR: invalid sync group"
          else
            ! get the index of this pusher in pusher_timestep
            pusher = count(pusher_timestep_fixed(1:sim%groups(events(i)%sync_groups(j))%pusher))
            constraints(i,pusher) = .true.
          end if
        end do
      else ! sync all groups
        constraints(i,:) = .true. ! add all the pushers for this event
      end if
    end do

    call fix_event_timestep(pusher_timestep, event_start, event_step, constraints, ierr)
    if (ierr .ne. 0) return

    ! show changes in timesteps and set them
    j = 0
    do i=1,size(pushers)
      if (pusher_timestep_fixed(i)) then
        j = j + 1
        if (abs(pusher_timestep(j) - pushers(i)%pusher%fixed_timestep) .gt. TICK) then
          write(*,*) "INFO: changing timestep of pusher", i, " from ", pushers(i)%pusher%fixed_timestep, " to ", pusher_timestep(j)
        end if
        ! always update, but notify only for significant changes
        pushers(i)%pusher%fixed_timestep = pusher_timestep(j)
      end if
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
  end if
end subroutine check_and_fix_timesteps
end module mod_main_loop
