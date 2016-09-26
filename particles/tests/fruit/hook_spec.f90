!> This module contains some trivial tests of the boris pusher
module hook_spec
use mod_constants, only: CARTESIAN
use mod_prescribed_fields
use mod_pusher_no_action
use mod_boris, only: particle_boris
use mod_hook
use fruit
implicit none

type, extends(particle_action) :: particle_action_half_velocity
contains
  procedure :: do => half_particle_velocity
end type particle_action_half_velocity
contains

pure subroutine half_particle_velocity(this, particle)
  class(particle_action_half_velocity), intent(in) :: this
  class(particle_base), intent(inout) :: particle
  select type (particle)
    type is (particle_boris)
      particle%v = particle%v / 2.d0
  end select
end subroutine half_particle_velocity

!> Create a hook to half the velocity of the particle every time it runs.
subroutine test_hook_half_velocity
  type(particle_boris)    :: particle
  type(pusher_no_action)  :: pusher
  type(prescribed_fields) :: fields
  real*8, parameter    :: timestep = 1d-1, time = 1d0
  real*8 :: x(2)

  fields = prescribed_fields(CARTESIAN, E_zero, B_minus_z)
  pusher = pusher_no_action(fixed_timestep=timestep, hooks=[hook_base(particle_action_half_velocity())])
  particle%v(:)  = [0d0, 1d0, 0d0]

  call pusher%push_single(fields, particle, 0.d0, time)
  call assert_equals(2.d0**(-10), norm2(particle%v), 'particle velocity must be 1/2**10')
end subroutine test_hook_half_velocity

pure function E_zero(x, t) result(E)
  real*8, dimension(3), intent(in) :: x
  real*8, intent(in) :: t
  real*8, dimension(3) :: E
  E = (/0.d0, 0.d0, 0.d0/)
end function E_zero
pure function B_minus_z(x, t) result(B)
  real*8, dimension(3), intent(in) :: x
  real*8, intent(in) :: t
  real*8, dimension(3) :: B
  B = (/0.d0, 0.d0, -1.d0/)
end function B_minus_z
end module hook_spec
