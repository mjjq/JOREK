!> Module for using predefined fields
module mod_prescribed_fields
  use mod_fields
  implicit none

  type, extends(type_fields) :: type_prescribed_fields
    procedure(position_dependent_field), pointer, private, nopass :: electric_field
    procedure(position_dependent_field), pointer, private, nopass :: magnetic_field
    contains
      procedure :: at_particle => at_particle_impl
  end type type_prescribed_fields

  interface
    pure function position_dependent_field(x)
      implicit none
      real*8, dimension(3), intent(in) :: x
      real*8, dimension(3) :: position_dependent_field
    end function position_dependent_field
  end interface

contains
  pure subroutine at_particle_impl(this, particle, E, B, psi, U)
    use mod_particle_type
    implicit none
    class(type_prescribed_fields), intent(in) :: this
    class(type_particle), intent(in)      :: particle
    real*8, dimension(3), intent(out)     :: E, B
    real*8, intent(out)                   :: psi, U

    psi = 0.d0
    U = 0.d0 ! artefact of JOREK integration

    E = this%electric_field(particle%x)
    B = this%magnetic_field(particle%x)
  end subroutine at_particle_impl
end module mod_prescribed_fields
