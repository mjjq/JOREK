!> Particle input-output event module
module mod_io_events
use mod_particle_io
use mod_events
implicit none

type, private, extends(event), abstract :: io_event
  character(len=80) :: filename
end type io_event

type, public, extends(io_event) :: read_state
contains
  procedure :: do => do_read_state
end type read_state
interface read_state
  module procedure new_read_state
end interface read_state

type, public, extends(io_event) :: write_state
contains
  procedure :: do => do_write_state
end type write_state
interface write_state
  module procedure new_write_state
end interface write_state

public :: new_read_state, new_write_state
private ! everything else is private
contains


!> Constructor for read_state
function new_read_state()
  type(read_state) :: new_read_state
  new_read_state%sync = .true.
  new_read_state%name = "ReadState"
end function new_read_state

!> Action for reading the state
subroutine do_read_state(this, sim)
  class(read_state), intent(inout)  :: this
  type(particle_sim), intent(inout) :: sim
end subroutine do_read_state


!> Constructor for write_state
function new_write_state()
  type(write_state) :: new_write_state
  new_write_state%sync = .true.
  new_write_state%name = "writeState"
end function new_write_state

!> Action for writing the state
subroutine do_write_state(this, sim)
  class(write_state), intent(inout) :: this
  type(particle_sim), intent(inout) :: sim
end subroutine do_write_state
end module mod_io_events
