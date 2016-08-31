!> Base module for a testcase, such as the penning trap, a tokamak equilibrium
!> or the grad-B drift
module mod_case
  use mod_particle_base
  use mod_pusher
  use mod_prescribed_fields
  implicit none

  !> 
  type, abstract :: case
    real*8                  :: time_end !< How long to simulate for
    type(prescribed_fields) :: fields !< The fields used in this case
  contains
    procedure, pass, public :: run
    procedure(initialize_particle), deferred, pass, private :: initialize_particle
    procedure(calculate_error), deferred, pass, private     :: calculate_error
  end type
  interface
    subroutine initialize_particle(this, particle, fixed_timestep)
      import :: case, particle_base
      class(case), intent(in)             :: this
      class(particle_base), intent(inout) :: particle
      real*8, intent(in), optional        :: fixed_timestep
    end subroutine initialize_particle
    function calculate_error(this, particle)
      import :: case, particle_base
      class(case), intent(in)          :: this
      class(particle_base), intent(in) :: particle
      real*8 :: calculate_error
    end function calculate_error
  end interface
contains

!> Run a single testcase with the given particle and pusher and report the results
subroutine run(this, pusher, particle, err, runtime)
  class(case), intent(in)        :: this
  class(pusher_base), intent(in) :: pusher
  class(particle_base), intent(inout) :: particle
  real*8, intent(out) :: err
  real*8, intent(out) :: runtime

  real*8 :: time_start
  call cpu_time(time_start)
  call pusher%push_single(this%fields, particle, 0.d0, this%time_end)
  call cpu_time(runtime)
  runtime = runtime - time_start

  err = this%calculate_error(particle)
end subroutine run
end module mod_case
