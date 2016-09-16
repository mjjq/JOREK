!> Base module for a testcase, such as [[mod_penning_case]], [[mod_gradb_case]]
module mod_case
  use mod_particle_base, only: particle_base
  use mod_pusher, only: pusher_base
  use mod_prescribed_fields, only: prescribed_fields
  implicit none

  !> A case, defining an end time, a run function and requiring an 
  !> initialization routine and error calculation routine
  type, abstract :: case
    real*8                  :: time_end !< How long to simulate for
    type(prescribed_fields) :: fields !< The fields used in this case
  contains
    procedure, pass, public :: run
    procedure, pass, public :: initialize_particle_all
    procedure(initialize_particle), deferred, pass, public :: initialize_particle
    procedure(calculate_error), deferred, pass, private :: calculate_error
  end type
  interface
    subroutine initialize_particle(this, particle, pusher)
      import :: case, particle_base, pusher_base
      class(case), intent(in)             :: this
      class(particle_base), intent(inout) :: particle
      class(pusher_base), intent(in)      :: pusher
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
subroutine run(this, pusher, particle, err, runtime, output_file)
  class(case), intent(in)        :: this
  class(pusher_base), intent(in) :: pusher !< the particle pusher to use
  class(particle_base), intent(inout) :: particle !< the particle to push
  real*8, intent(out) :: err !< The error as defined by [[calculate_error]]
  real*8, intent(out) :: runtime !< the runtime as reported by the cpu_time intrinsic
  character(len=*), optional :: output_file
  real*8 :: time_start
  integer :: u, i
  integer, parameter :: requested_output_points = 10000

  if (present(output_file)) then
    open(newunit=u, file=output_file)
    if (allocated(pusher%fixed_timestep)) then
      do i=1,nint(this%time_end/pusher%fixed_timestep)
        call pusher%push_single(this%fields, particle, 0.d0, pusher%fixed_timestep)
        write(u,*) particle%x
      end do
    else
      do i=1,requested_output_points
        call pusher%push_single(this%fields, particle, 0.d0, this%time_end/requested_output_points)
        write(u,*) particle%x
      end do
    end if
    runtime = 0.d0
    close(u)
  else
    call cpu_time(time_start)
    call pusher%push_single(this%fields, particle, 0.d0, this%time_end)
    call cpu_time(runtime)
    runtime = runtime - time_start
  end if

  err = this%calculate_error(particle)
end subroutine run

!> Initialize a particle for a testcase (perhaps move this to the particle init modules?)
subroutine initialize_particle_all(this, particle, pusher, x0, v0, charge, mass)
  use mod_coordinate_transforms
  use mod_boris
  use mod_constants, only: CARTESIAN, CYLINDRICAL
  class(case), intent(in)             :: this
  class(particle_base), intent(inout) :: particle
  class(pusher_base), intent(in)      :: pusher
  real*8, dimension(3), intent(in)    :: x0, v0
  integer*1, intent(in)               :: charge
  real*4, intent(in)                  :: mass
  if (this%fields%geometry .eq. CARTESIAN)   particle%x = x0
  if (this%fields%geometry .eq. CYLINDRICAL) particle%x = cartesian_to_cylindrical(x0)
  particle%q = charge
  particle%m = mass
  particle%lost = .false.
  select type (particle)
  type is (particle_boris)
    ! set up the velocity
    if (this%fields%geometry .eq. CARTESIAN)   particle%v = v0
    if (this%fields%geometry .eq. CYLINDRICAL) particle%v = vector_rotation(v0, particle%x(2))
    select type (pusher)
    type is (pusher_boris)
      call pusher%initial_half_step_backwards(this%fields, particle)
    class default
      write(*,*) "ERROR: passed incompatible particle and pusher"
    end select
  class default
    write(*,*) "ERROR: particle type not implemented in initialize_particle_all"
  end select
end subroutine initialize_particle_all
end module mod_case
