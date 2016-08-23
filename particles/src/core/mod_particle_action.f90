!> Module to formalize performing an action on a particle
module mod_particle_action
use mod_particle_base
implicit none

!> Action abstract type, representing anything that can be done to a simulation
type, abstract, public :: particle_action
contains
  procedure(do), deferred, pass :: do
end type particle_action
interface
  pure subroutine do(this, particle)
    import :: particle_action, particle_base
    class(particle_action), intent(in)  :: this
    class(particle_base), intent(inout) :: particle
  end subroutine do
end interface
end module mod_particle_action
