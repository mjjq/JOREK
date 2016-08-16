!> Module to formalize performing an action on a particle
module mod_particle_action
use mod_particle_base
use mod_particle_sim
implicit none

!> Action abstract type, representing anything that can be done to a simulation
type, abstract, public :: particle_action
contains
  procedure(do), deferred, pass :: do
end type particle_action
interface
  pure subroutine do(this, sim, particle)
    import :: particle_action, particle_base, particle_sim
    class(particle_action), intent(in)  :: this
    type(particle_sim), intent(in)      :: sim
    class(particle_base), intent(inout) :: particle
  end subroutine do
end interface
end module mod_particle_action
