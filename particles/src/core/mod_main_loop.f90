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
  integer :: ierr
  
  ! Check if pushers contain enough for all groups
  call check_pusher_groups(pushers, size(sim%group), ierr)
  if (ierr .ne. 0) return

  ! Calculate required timesteps
  call calculate_timesteps(pushers, events, ierr)
  if (ierr .ne. 0) return

  !$omp parallel default(none) &
  !$omp shared(sim, pushers, events) &
  !$omp private(ierr)

  do while (.true.)
  ! Calculate which of the events is next, or stop if there are no more
  ! push particles until that time
  ! run event on the master thread
  end do

  !$omp end parallel
end subroutine


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
end subroutine calculate_timesteps
end module mod_main_loop
