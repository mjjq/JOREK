!> Module to formalize performing an action every now and then
!> Events are objects (derived types) containing information
!> about when to run, and an action they should run.
module mod_event
use mod_particle_sim
use mod_event_timestep
implicit none
private
public action, stop_action, cycle_time_action
public event, with, next_event_at, check_and_fix_timesteps
public mpi_minmeanmax
public :: TICK

!> @name Time precision
real*8,  parameter :: TICK             = 1d-12                    !< Time precision for events [s]

!> Action abstract type, representing anything that can be done to a simulation
type, abstract :: action
  !> Logging variable, set this in an initializer
  character(len=50) :: name = "unset action" !< Event name for logging
  logical :: log = .false. !< Output event duration

  !> Timing variables
  real*8, private :: t0 = 0.d0, w0 = 0.d0
contains
  procedure, pass, public :: run
  procedure(do_interface), deferred, pass, private :: do
end type action
!> Example action (stops the simulation)
type, extends(action) :: stop_action
contains
  procedure :: do => do_stop_action
end type stop_action
interface stop_action
  module procedure new_stop_action
end interface
!> Timing output action (prints time since previous round)
type, extends(action) :: cycle_time_action
  real*8 :: last_t = 0.d0 !< Time at previous iteration
contains
  procedure :: do => do_cycle_time_action
end type cycle_time_action
interface cycle_time_action
  module procedure new_cycle_time_action
end interface






!> Event type
type :: event
  real*8  :: start    = 0.d0       !< Physical starting time
  real*8  :: step     = huge(0.d0) !< Step every how long?
  real*8  :: end      = huge(0.d0) !< Stop after time end. If equal to start, runs once

  integer, dimension(:), allocatable :: sync_groups !< which groups to require at a full-timestep (unallocated = all, empty array = none)

  !> Action to perform when this event runs
  class(action), allocatable :: action
contains
  procedure run_at
end type event
interface event
  module procedure new_event
end interface

interface with
  module procedure with_event_0D, with_action_0D, &
        with_event_1D, with_event_1D_at, with_event_1D_mask, &
        with_action_1D, with_action_1D_mask
end interface
interface
  subroutine do_interface(this, sim, ev)
    import :: action, particle_sim, event
    class(action), intent(inout)      :: this
    type(particle_sim), intent(inout) :: sim
    type(event), intent(inout), optional :: ev
  end subroutine do_interface
end interface


contains
!> Constructor for an event
!> This is needed to allow changing default values
function new_event(act, start, step, end, sync_groups, sync_groups_none)
  type(event), target           :: new_event
  class(action), intent(in)     :: act
  real*8, intent(in), optional  :: start, step, end
  integer, dimension(:), intent(in), optional :: sync_groups !< Groups to sync
  logical, intent(in), optional :: sync_groups_none !< Sync no groups. If sync_groups is present that takes precedence
  if (present(start))    new_event%start    = start
  if (present(step))     new_event%step     = step
  if (present(end))      new_event%end      = end
  allocate(new_event%action, source=act) ! because assignment is not yet supported in gfortran 6.1.1
  if (present(sync_groups)) then
    new_event%sync_groups = sync_groups
  else
    if (present(sync_groups_none)) allocate(new_event%sync_groups(0))
  end if
end function new_event

!> Should this event run at this time?
function run_at(this, time)
  class(event), intent(in) :: this
  real*8, intent(in) :: time
  logical :: run_at
  integer :: closest_iteration
  closest_iteration = nint((time - this%start)/this%step)
  run_at = .false.
  if (closest_iteration .ge. 0 &
    .and. abs(time - (this%start + closest_iteration*this%step)) .le. tick &
    .and. (time - this%end .le. tick)) then
    run_at = .true.
  end if
end function

subroutine with_event_0D(sim, single_event)
  type(particle_sim), intent(inout) :: sim
  type(event), intent(inout) :: single_event
  call single_event%action%run(sim, single_event)
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
    call events(i)%action%run(sim, events(i))
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
subroutine with_event_1D_at(sim, events, at)
  type(particle_sim), intent(inout) :: sim
  type(event), intent(inout), dimension(:) :: events
  real*8, intent(in) :: at
  integer :: i
  do i=1,size(events)
    if (events(i)%run_at(at)) call events(i)%action%run(sim, events(i))
  end do
end subroutine with_event_1D_at
subroutine with_event_1D_mask(sim, events, mask)
  type(particle_sim), intent(inout) :: sim
  type(event), intent(inout), dimension(:) :: events
  !logical, dimension(size(events,1)), intent(in) :: mask ! internal compiler error in gfortran, workaround below
  logical, dimension(:), intent(in) :: mask
  integer :: i
  do i=1,size(events)
    if (mask(i)) call events(i)%action%run(sim, events(i))
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
function next_event_at(sim, events) result(at)
  type(particle_sim), intent(inout) :: sim
  type(event), intent(in), dimension(:) :: events
  real*8 :: at
  logical, dimension(size(events)) :: run_event
  real*8 :: event_run !< when events(i) is to run (at the soonest)
  integer :: i

  at = huge(0.d0)
  run_event(:) = .false.
  do i=1,size(events)
    ! next event needs to be at least tick in the future
    if (events(i)%start .gt. sim%time + tick) then
      event_run = events(i)%start
    else
      event_run = sim%time + events(i)%step - mod(sim%time - events(i)%start, events(i)%step)
      if (abs(event_run - sim%time) .le. tick) event_run = sim%time + events(i)%step
    end if

    ! If this event has ended already
    if (event_run .gt. events(i)%end + tick) cycle

    ! if this event occurs faster than the previously fastest (event_time)
    if (event_run .lt. at - tick) then
      run_event(:) = .false.
      run_event(i) = .true.
      at = event_run
    else if (event_run .le. at + tick) then ! if it is equally fast
      run_event(i) = .true.
    end if
  end do

  ! Exit the simulation if there are no more events to do
  if (at .ge. maxval(events(:)%end)) then
    sim%stop_now = .true.
    if (at .eq. huge(0.d0)) at = sim%time ! if the next event is not occurring or at infinity keep the current time
  end if
end function next_event_at


!> Calculate whether we need to change any of the fixed timesteps or events to match
!> For each of the pushers with a fixed timestep
subroutine check_and_fix_timesteps(pusher_timestep, events)
  use mpi
  real*8,      intent(inout), dimension(:) :: pusher_timestep
  type(event), intent(inout), dimension(:) :: events

  real*8, dimension(:), allocatable    :: event_start, event_step, pusher_timestep_work
  logical, dimension(:,:), allocatable :: constraints
  integer :: i, j, ierr, my_id
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)

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
      if (my_id .eq. 0) write(*,*) "INFO: changing timestep of pusher", i, " from ", pusher_timestep(i), " to ", pusher_timestep_work(i)
    end if
    ! always update, but notify only for significant changes
    pusher_timestep(i) = pusher_timestep_work(i)
  end do
  do i=1,size(events)
    if (abs(event_start(i) - events(i)%start) .gt. TICK) then
      if (my_id .eq. 0) write(*,"(A,i3,A,A,A,g16.8,A,g16.8,A,g14.6,A)") "INFO: changing start time of event ", i, &
          "(", trim(events(i)%action%name), ") from ", events(i)%start, " to ", event_start(i), ' by ', event_start(i)-events(i)%start
      events(i)%start = event_start(i)
    end if
    if (abs(event_step(i) - events(i)%step) .gt. TICK) then
      if (my_id .eq. 0) write(*,"(A,i3,A,A,A,g16.8,A,g16.8,A,g14.6,A)") "INFO: changing timestep of event ", i, &
          "(", trim(events(i)%action%name), ") from ", events(i)%step, " to ", event_step(i), ' by ', event_step(i)-events(i)%step
      events(i)%step = event_step(i)
    end if
  end do
end subroutine check_and_fix_timesteps





!> Constructor for stop_action
function new_stop_action()
  type(stop_action) :: new_stop_action
  new_stop_action%name = "Stop"
  new_stop_action%log  = .false.
end function new_stop_action

!> Perform the stop action
subroutine do_stop_action(this, sim, ev)
  class(stop_action), intent(inout) :: this
  type(particle_sim), intent(inout) :: sim
  type(event), intent(inout), optional :: ev
  sim%stop_now = .true.
end subroutine do_stop_action




!> Perform the cycle_time action
subroutine do_cycle_time_action(this, sim, ev)
  use mpi
  class(cycle_time_action), intent(inout) :: this
  type(particle_sim), intent(inout) :: sim
  type(event), intent(inout), optional :: ev
  real*8 :: t1, m4(4)
  t1 = MPI_WTIME()

  if (this%last_t .eq. 0.d0) then
    this%last_t = t1
  else
    m4 = mpi_minmeanmedmax(t1 - this%last_t)
    if (sim%my_id .eq. 0) write(*,"(A,4f8.3,A)") "Cycle time (min/mean/median/max): ", m4, "s"
  end if
end subroutine do_cycle_time_action
!> Constructor for stop_action
function new_cycle_time_action() result(new)
  type(cycle_time_action) :: new
  new%name = "CycleTime"
  new%log  = .false.
end function new_cycle_time_action



!> Run an action
subroutine run(this, sim, ev)
  !$ use omp_lib
  use mpi
  class(action), intent(inout)      :: this
  type(particle_sim), intent(inout) :: sim
  type(event), intent(inout), optional :: ev !< If run from an event, get a pointer to it here
  real*8 :: t1, w1, mmm(3), mmm2(3), m4(4)
  logical :: has_omp
  integer :: ierr
  has_omp = .false.
  !$ has_omp = .true.

  if (this%log) then
    ! Calculate imbalance between MPI processes
    if (MPI_WTIME_IS_GLOBAL .eq. 1) then
      t1 = MPI_WTIME()
    else
      ! Be aware that clock skew happens! That will distort the results over time
      ! The load balancing action has an explicit mpi_barrier and could act as a
      ! reset point for the clocks. Note that mpi_barrier only implies execution
      ! synchronisation, not time synchronisation, which means that sim%wtime_start
      ! is probably only accurate up to the network latency.
      !
      ! If MPI_WTIME_IS_GLOBAL is set the MPI implementation will attempt to
      ! synchronise clocks and it is actually better not to use sim%wtime_start.
      t1 = MPI_WTIME() - sim%wtime_start
    end if
    m4 = mpi_minmeanmedmax(t1)
    m4 = m4 - m4(1) ! Measure times from first process to reach (i.e. min = 0)
    if (m4(4) .gt. 0.1) then ! only write if there is a significant imbalance
      if (sim%my_id .eq. 0) write(*,"(A,A,3f8.3,A)") trim(this%name), " MPI entry imbalance (mean/median/max): ", m4(2:4), "s"
    end if
  end if


  call cpu_time(this%t0)
  !$ this%w0 = omp_get_wtime()
  if (present(ev)) then
    call this%do(sim, ev)
  else
    call this%do(sim)
  end if
  call cpu_time(t1)
  !$ w1 = omp_get_wtime()

  if (this%log) then
    if (sim%n_cpu .gt. 1) then
      mmm = mpi_minmeanmax(t1-this%t0)
      !$ mmm2 = mpi_minmeanmax(w1-this%w0)
      if (sim%my_id .eq. 0) then
        if (.not. has_omp) write(*,"(A,A,3f9.3,A)") trim(this%name), " finished in (min/mean/max): ", &
            mmm, "s"
        !$ write(*,"(A,A,3f8.3,A,3f8.3,A)") trim(this%name), " finished in (min/mean/max): ", &
        !$   mmm2, &
        !$ "s (cpu time: ", mmm, ")"
      end if
    else
      if (.not. has_omp) write(*,"(A,A,f8.3,A)") trim(this%name), " finished in: ", &
         t1-this%t0, "s"
      !$ write(*,"(A,A,f8.3,A,f8.3,A)") trim(this%name), " finished in: ", &
      !$   w1-this%w0, &
      !$ "s (cpu time: ", t1-this%t0, ")"
    end if

    ! Recalculate imbalance between MPI processes
    if (MPI_WTIME_IS_GLOBAL .eq. 1) then
      t1 = MPI_WTIME()
    else
      t1 = MPI_WTIME() - sim%wtime_start
    end if
    m4 = mpi_minmeanmedmax(t1)
    m4 = m4 - m4(1) ! Measure times from first process to reach (i.e. min = 0)
    if (m4(4) .gt. 0.1) then ! only write if there is a significant imbalance
      if (sim%my_id .eq. 0) write(*,"(A,A,3f8.3,A)") trim(this%name), " MPI exit imbalance (mean/median/max): ", m4(2:4), "s"
    end if
  end if
end subroutine run

!> Calculate min, mean and max values over all nodes
!> returns something only on root process.
function mpi_minmeanmax(in) result(out)
  use mpi
  real*8, intent(in) :: in
  real*8 :: out(3)
  integer :: n_cpu, my_id, ierr
  real*8, allocatable :: in_all(:)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
  allocate(in_all(n_cpu))
  call MPI_Gather(in, 1, MPI_REAL8, in_all, 1, MPI_REAL8, 0, MPI_COMM_WORLd, ierr)
  if (my_id .eq. 0) then
    out(1) = minval(in_all,1)
    out(2) = sum(in_all,1)/real(n_cpu)
    out(3) = maxval(in_all,1)
  else
    out = 0.d0
  end if
end function mpi_minmeanmax

!> Calculate min, mean and max values over all nodes
!> returns something only on root process.
function mpi_minmeanmedmax(in) result(out)
  use mpi
  use mod_quicksort
  real*8, intent(in) :: in
  real*8 :: out(4)
  integer :: n_cpu, my_id, ierr
  real*8, allocatable :: in_all(:)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
  allocate(in_all(n_cpu))
  call MPI_Gather(in, 1, MPI_REAL8, in_all, 1, MPI_REAL8, 0, MPI_COMM_WORLd, ierr)
  if (my_id .eq. 0) then
    call quicksort(in_all)
    out(1) = in_all(1)
    out(2) = sum(in_all,1)/real(n_cpu)
    out(3) = in_all((n_cpu+1)/2)
    out(4) = in_all(n_cpu)
  else
    out = 0.d0
  end if
end function mpi_minmeanmedmax
end module mod_event
