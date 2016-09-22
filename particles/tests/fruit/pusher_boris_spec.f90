!> This module contains some trivial tests of the boris pusher
module pusher_boris_spec
use mod_constants, only: CARTESIAN, TWOPI, EL_CHG, ATOMIC_MASS_UNIT
use mod_prescribed_fields, only: prescribed_fields
use mod_boris, only: pusher_boris, particle_boris
use fruit
implicit none
contains

!> Test the gyromotion of a single particle, orbiting with radius 1 in 1 second
subroutine test_single_gyro_orbit_dt0_01
  real*8 :: x(2)
  x = x_orbit(0.01d0, 0.25d0)
  call assert_equals(0.d0, x(1), 1d-1, "x(0.25) = 0")
  call assert_equals(1.d0, x(2), 1d-1, "y(0.25) = 1")
  x = x_orbit(0.01d0, 0.5d0)
  call assert_equals(-1.d0, x(1), 1d-3, "x(0.50) = -1")
  call assert_equals(0.d0,  x(2), 1d-1, "y(0.50) = 0")
  x = x_orbit(0.01d0, 0.75d0)
  call assert_equals(0.d0,  x(1), 1d-1, "x(0.75) = 0")
  call assert_equals(-1.d0, x(2), 1d-1, "y(0.75) = -1")
  x = x_orbit(0.01d0, 1.0d0)
  call assert_equals(1.d0, x(1), 1d-7, "x(1.00) = 1")
  call assert_equals(0.d0, x(2), 1d-7, "y(1.00) = 0") ! why does this converge so nicely?
end subroutine test_single_gyro_orbit_dt0_01

!> TODO: we still have linear convergence for the y-component (which is also the quickly-varying one)
subroutine test_half_gyro_orbit_convergence
  real*8 :: x(2)
  x = x_orbit(1d-2, 0.5d0)
  call assert_equals(0.d0, x(1)+1.d0, 1d-3, "dt=1d-2 x(0.50) = -1")
  call assert_equals(0.d0, x(2)     , 1d-1, "dt=1d-2 y(0.50) = 0")
  x = x_orbit(1d-3, 0.5d0)
  call assert_equals(0.d0, x(1)+1.d0, 1d-5, "dt=1d-3 x(0.50) = -1")
  call assert_equals(0.d0, x(2)     , 1d-2, "dt=1d-3 y(0.50) = 0")
  x = x_orbit(1d-4, 0.5d0)
  call assert_equals(0.d0, x(1)+1.d0, 1d-7, "dt=1d-4 x(0.50) = -1")
  call assert_equals(0.d0, x(2)     , 1d-3, "dt=1d-4 y(0.50) = 0")
  x = x_orbit(1d-5, 0.5d0)
  call assert_equals(0.d0, x(1)+1.d0, 1d-7, "dt=1d-5 x(0.50) = -1")
  call assert_equals(0.d0, x(2)     , 1d-4, "dt=1d-5 y(0.50) = 0")
  x = x_orbit(1d-6, 0.5d0)
  call assert_equals(0.d0, x(1)+1.d0, 1d-7, "dt=1d-6 x(0.50) = -1")
  call assert_equals(0.d0, x(2)     , 1d-5, "dt=1d-6 y(0.50) = 0")
end subroutine test_half_gyro_orbit_convergence

function x_orbit(timestep, time) result(x)
  type(particle_boris) :: particle
  type(pusher_boris)   :: pusher
  type(prescribed_fields) :: fields
  real*8, intent(in)   :: timestep, time
  real*8 :: x(2)

  fields = prescribed_fields(CARTESIAN, E_zero, B_minus_z)
  pusher = pusher_boris(fixed_timestep=timestep)

  particle%x(:)  = [1.d0, 0.d0, 0.d0]
  particle%v(:)  = [0.d0, TWOPI, 0.d0] ! v^(-1/2), see Delzanno, JCP (2013)
  particle%q     = 1 ! +1 e
  particle%m     = EL_CHG/TWOPI/ATOMIC_MASS_UNIT ! mass in unified atomic mass units to have f=1Hz in a field of 1 Tesla
  particle%lost = .false.

  call pusher%push_single(fields, particle, 0.d0, time)
  x = particle%x(1:2)
end function

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
end module pusher_boris_spec
