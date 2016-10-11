!> This module contains some trivial tests of the boris pusher
module pusher_boris_spec
use constants, only: CARTESIAN, TWOPI, EL_CHG, ATOMIC_MASS_UNIT
use mod_particle_types
use mod_boris
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
  type(particle_kinetic_leapfrog) :: particle
  real*8, intent(in)   :: timestep, time
  real*8 :: x(2)
  integer :: i

  particle%x(:)  = [1.d0, 0.d0, 0.d0]
  particle%v(:)  = [0.d0, TWOPI, 0.d0] ! TODO get accurate v^(-1/2), see Delzanno, JCP (2013) and pusher_test
  particle%q     = 1 ! +1 e
  particle%m     = EL_CHG/TWOPI/ATOMIC_MASS_UNIT ! mass in unified atomic mass units to have f=1Hz in a field of 1 Tesla
  particle%lost = .false.

  do i=1,nint(time/timestep)
    call boris_push_cartesian(particle, [0d0,0d0,0d0], [0d0,0d0,-1d0], timestep)
  end do
  x = particle%x(1:2)
end function
end module pusher_boris_spec
