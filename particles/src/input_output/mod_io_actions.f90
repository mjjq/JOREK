!> Particle input-output event module
module mod_io_actions
use mod_particle_io
use mod_action
use mod_particle_sim
implicit none
private
public read_action, write_action, &
    get_filename ! public because we test it externally

type, extends(action), abstract :: io_action
  character(len=120), public :: filename = '' !< Filename to use
  character(len=80), public  :: basename = 'part' !< If no filename, use basename + digits
  integer, public            :: decimal_digits = 3 !< Number of decimals before the point in timestamp
  integer, public            :: fractional_digits = 8 !< Number of decimals after the point
  character(len=5), public   :: extension = '.h5'
  contains
    procedure :: get_filename
end type io_action

type, extends(io_action) :: read_action
  real*8 :: time !< used with the formats from io_action if filename is unset
contains
  procedure :: do => do_read_action
end type read_action
interface read_action
  module procedure new_read_action
end interface read_action

type, extends(io_action) :: write_action
contains
  procedure :: do => do_write_action
end type write_action
interface write_action
  module procedure new_write_action
end interface write_action

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
  else if (this%fractional_digits .eq. 0) then
    write(format,'(A,I0,A)') '(a,i', this%decimal_digits, 'a)'
    write(filename,trim(format)) trim(this%basename), int(time), this%extension
  else
    write(format,'(A,I0,A,I0,A,I0,A)') '(a,i', this%decimal_digits, '.', this%decimal_digits, &
        ',f0.', this%fractional_digits, ',a)'
    write(filename,trim(format)) trim(this%basename), int(time), time-int(time), this%extension
  end if
end function get_filename

!> Constructor for read_action.
!> Be sure to use keyword arguments when initializing, to avoid confusion
function new_read_action(filename, basename, decimal_digits, fractional_digits, extension)
  type(read_action) :: new_read_action
  character(len=*), intent(in), optional :: filename
  character(len=*), intent(in), optional :: basename
  integer, intent(in), optional          :: decimal_digits
  integer, intent(in), optional          :: fractional_digits
  character(len=*), intent(in), optional :: extension
  if (present(filename)) new_read_action%filename = filename
  if (present(basename)) new_read_action%basename = basename
  if (present(decimal_digits)) new_read_action%decimal_digits = decimal_digits
  if (present(fractional_digits)) new_read_action%fractional_digits = fractional_digits
  if (present(extension)) new_read_action%extension = extension
  new_read_action%name = "ReadAction"
end function new_read_action

!> Action for reading the simulation
subroutine do_read_action(this, sim)
  class(read_action), intent(inout) :: this
  type(particle_sim), intent(inout) :: sim
  if (len_trim(this%filename) .eq. 0) then
    call read_simulation_hdf5(sim, this%get_filename(this%time))
  else
    call read_simulation_hdf5(sim, trim(this%filename))
  end if
end subroutine do_read_action


!> Constructor for write_action
!> Be sure to use keyword arguments when initializing, to avoid confusion
function new_write_action(filename, basename, decimal_digits, fractional_digits, extension)
  type(write_action) :: new_write_action
  character(len=*), intent(in), optional :: filename
  character(len=*), intent(in), optional :: basename
  integer, intent(in), optional          :: decimal_digits
  integer, intent(in), optional          :: fractional_digits
  character(len=*), intent(in), optional :: extension
  new_write_action%name = "WriteAction"
  if (present(filename)) new_write_action%filename = filename
  if (present(basename)) new_write_action%basename = basename
  if (present(decimal_digits)) new_write_action%decimal_digits = decimal_digits
  if (present(fractional_digits)) new_write_action%fractional_digits = fractional_digits
  if (present(extension)) new_write_action%extension = extension
end function new_write_action

!> Action for writing the simulation
subroutine do_write_action(this, sim)
  class(write_action), intent(inout) :: this
  type(particle_sim), intent(inout)  :: sim
  if (len_trim(this%filename) .eq. 0) then
    call write_simulation_hdf5(sim, this%get_filename(sim%time))
  else
    call write_simulation_hdf5(sim, trim(this%filename))
  end if
end subroutine do_write_action
end module mod_io_actions
