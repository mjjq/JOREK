!> Fields as used in the simulation
module mod_fields
use mod_constants
implicit none

!> Base type for different field
type, abstract :: fields_base
  integer*1 :: geometry
contains
  procedure(at_particle), deferred :: at_particle
end type fields_base

interface
  pure subroutine at_particle(this, particle, t, E, B, psi, U)
    use mod_particle_base
    import :: fields_base
    implicit none
    class(fields_base), intent(in)    :: this
    class(particle_base), intent(in)  :: particle
    real*8, intent(in)                :: t !< The current time
    real*8, dimension(3), intent(out) :: E, B
    real*8, intent(out)               :: psi, U
  end subroutine at_particle
end interface
end module mod_fields
