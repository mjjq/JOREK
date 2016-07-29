!> Module to formalize performing an action every now and then
!> Events are objects (derived types) containing information
!> about when to run, and a function they should run.
!> Event is the base class, meant to be extended.
!> Each extending class can give a different signature for the function to run
module events
implicit none

!> Base type for all events
type, abstract :: type_event
  real*8  :: tstart = 0.d0 !< Physical starting time (default 0)
  integer :: istart = 0 !< Time step starting time (default 0)
  integer :: istep  = huge(0) !< Repeat every istep
  real*8  :: tend   = huge(0.d0)  !< Stop at time tend
  integer :: iend   = huge(0) !< Stop at istep iend
contains
  private
  procedure, private :: can_start
end type type_event

contains
  !> Function to test whether an event can run now
  logical function can_start(this, index_now, t_now)
    implicit none
    class(type_event), intent(in) :: this
    integer, intent(in)           :: index_now
    real*8 , intent(in)           :: t_now

    can_start = (t_now .gt. this%tstart .and. index_now .gt. this%istart .and. &
                 t_now .lt. this%tend   .and. index_now .lt. this%iend   .and. &
                 mod(index_now - this%istart, this%istep) .eq. 0)
  end function can_start
end module events
