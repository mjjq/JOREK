!> Particle pusher base module.
!>
!>## How to create new pushers
!> The best starting point is an existing pusher. Copy the file and decide which
!> base particle to use (or create a new one, see [[mod_particle_types]] for this).
!>
!> Modifications then need to be made in the following files
!>* pushers/mod_pusher_NAME.f90 (see [[mod_pusher_no_action]] for an annotated example)
!>* [[particle_tracer]] (add your module to the imports list)
module mod_pusher
use mod_hook
implicit none
private
public pusher_base, pusher_container

!> Abstract type to be extended by new pushers
type, abstract :: pusher_base
  integer, dimension(:), allocatable :: groups !< groups from which to push particles, or all groups if unallocated
  real*8, allocatable :: fixed_timestep !< Fixed timestep, if changing the timestep is not supported

  type(hook_base), dimension(:), allocatable :: hooks
contains
  procedure(push_single), deferred :: push_single
end type pusher_base

!> Container to allow having an array of different pushers
type :: pusher_container
  class(pusher_base), allocatable :: pusher
end type pusher_container
interface pusher_container
  module procedure new_pusher_container
end interface pusher_container


interface
  !> Push a single particle from time_start to time_end (or a bit further if the timestep is fixed)
  subroutine push_single(this, fields, particle, time_start, time_end)
    use mod_particle_types
    use mod_fields
    import :: pusher_base
    class(pusher_base), intent(in)      :: this
    class(fields_base), intent(in)      :: fields
    class(particle_base), intent(inout) :: particle
    real*8, intent(in) :: time_start, time_end
  end subroutine push_single
end interface

contains
!> Constructor for a pusher-container
!> contains an allocate statement to convert type(pusher_*) to class(pusher_base)
pure function new_pusher_container(pusher)
  type(pusher_container) :: new_pusher_container
  class(pusher_base), intent(in) :: pusher
  allocate(new_pusher_container%pusher, source=pusher)
end function new_pusher_container
end module mod_pusher
