!> Base module for a testcase, such as [[mod_penning_case]], [[mod_gradb_case]]
module mod_case
  use mod_particle_types
  implicit none

  !> A case, defining an end time, a run function and requiring an 
  !> initialization routine and error calculation routine
  type, abstract :: case
  contains
    procedure(field), nopass, public, deferred :: E, B
    procedure(initialize_particle), deferred, pass, public :: initialize_particle
    procedure(calculate_error), deferred, pass, private :: calculate_error
  end type
  interface
    pure function field(x, t)
      real*8, dimension(3), intent(in) :: x
      real*8, intent(in) :: t
      real*8, dimension(3) :: field 
    end function field
    subroutine initialize_particle(this, particle)
      import :: case, particle_base
      class(case), intent(in)             :: this
      class(particle_base), intent(inout) :: particle
    end subroutine initialize_particle
    function calculate_error(this, particle)
      import :: case, particle_base
      class(case), intent(in)          :: this
      class(particle_base), intent(in) :: particle
      real*8 :: calculate_error
    end function calculate_error
  end interface
contains
end module mod_case
