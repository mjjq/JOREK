!> Module for using predefined fields
module mod_prescribed_fields
use mod_fields
implicit none
private
public prescribed_fields

!> Use a function of x and t for the electric and magnetic field
type, extends(fields_base) :: prescribed_fields
  procedure(position_dependent_field), pointer, public, nopass :: electric_field
  procedure(position_dependent_field), pointer, public, nopass :: magnetic_field
contains
  procedure :: at_particle => at_particle_impl
end type prescribed_fields

interface
  pure function position_dependent_field(x, t)
    real*8, dimension(3), intent(in) :: x
    real*8, intent(in) :: t
    real*8, dimension(3) :: position_dependent_field
  end function position_dependent_field
end interface

contains
pure subroutine at_particle_impl(this, particle, t, E, B, psi, U)
  use mod_particle_types
  class(prescribed_fields), intent(in) :: this
  class(particle_base), intent(in)     :: particle
  real*8, intent(in)                   :: t !< The current time
  real*8, dimension(3), intent(out)    :: E, B
  real*8, intent(out), optional        :: psi, U

  if (present(psi)) psi = 0.d0
  if (present(U))   U   = 0.d0 ! artefact of JOREK integration

  E = this%electric_field(particle%x, t)
  B = this%magnetic_field(particle%x, t)
end subroutine at_particle_impl
end module mod_prescribed_fields
