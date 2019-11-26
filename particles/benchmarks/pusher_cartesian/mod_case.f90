!> Base module for a testcase, such as [[mod_penning_case]], [[mod_gradb_case]]
module mod_case
  use mod_particle_types, only: particle_base
  implicit none

  private
  public :: case

  !> A case, defining an end time, a run function and requiring an
  !> initialization routine and error calculation routine
  type, abstract :: case
  contains
    procedure(field), nopass, public, deferred :: E
    procedure(field), nopass, public, deferred :: B
    procedure(initialize), deferred, pass, public :: initialize
    procedure(calc_error), deferred, pass, public :: calc_error
  end type
  interface
    pure function field(x, t)
      real*8, dimension(3), intent(in) :: x
      real*8, intent(in) :: t
      real*8, dimension(3) :: field
    end function field
    subroutine initialize(this, particle)
      import :: case, particle_base
      class(case), intent(in)             :: this
      class(particle_base), intent(inout) :: particle
    end subroutine initialize
    function calc_error(this, particle)
      import :: case, particle_base
      class(case), intent(in)          :: this
      class(particle_base), intent(in) :: particle
      real*8 :: calc_error
    end function calc_error
  end interface
end module mod_case
