!> Particle input-output event module
module mod_io_actions
use mod_particle_io
use mod_action
use mod_particle_sim
implicit none

type, private, extends(action), abstract :: io_action
  character(len=80) :: basename = 'part'
  integer           :: decimal_digits = 3
  integer           :: fractional_digits = 8
  character(len=5)  :: extension = '.h5'
  contains
    procedure :: get_filename
end type io_action

type, public, extends(io_action) :: read_action
  real*8 :: time
contains
  procedure :: do => do_read_action
end type read_action
interface read_action
  module procedure new_read_action
end interface read_action

type, public, extends(io_action) :: write_action
contains
  procedure :: do => do_write_action
end type write_action
interface write_action
  module procedure new_write_action
end interface write_action

private
contains
!> Write a filename consisting of this%basename and time as floating-point
!> number, without spaces (tricky), with this%decimal_digits before the `.` and
!> this%fractional_digits behind the `.`. Set these two to 0 to just set the name.
function get_filename(this, time) result(filename)
  class(io_action), intent(in) :: this
  real*8, intent(in)           :: time
  character(len=120)           :: filename
  character(len=20)            :: format
  if (this%decimal_digits .eq. 0 .and. this%fractional_digits .eq. 0) then
    write(filename,'(A,A)') trim(this%basename), this%extension
  else
    write(format,'(A,I0,A,I0,A,I0,A)') '(a,i', this%decimal_digits, '.', this%decimal_digits, &
        ',f0.', this%fractional_digits, ',a)'
    write(filename,trim(format)) trim(this%basename), int(time), time-int(time), this%extension
  end if
end function get_filename

  



!> Constructor for read_action
function new_read_action()
  type(read_action) :: new_read_action
  new_read_action%name = "ReadAction"
end function new_read_action

!> Action for reading the simulation
subroutine do_read_action(this, sim)
  class(read_action), intent(inout) :: this
  type(particle_sim), intent(inout) :: sim
  call read_simulation_hdf5(sim, this%get_filename(this%time))
end subroutine do_read_action


!> Constructor for write_action
function new_write_action()
  type(write_action) :: new_write_action
  new_write_action%name = "WriteAction"
end function new_write_action

!> Action for writing the simulation
subroutine do_write_action(this, sim)
  class(write_action), intent(inout) :: this
  type(particle_sim), intent(inout)  :: sim
  call write_simulation_hdf5(sim, this%get_filename(sim%time))
end subroutine do_write_action
end module mod_io_actions
