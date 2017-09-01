!> This module contains some simple tests of the boris pusher in cartesian geometry.
module pusher_boris_penning_spec
use mod_particle_types
use mod_boris, only: boris_push_cartesian, boris_push_cylindrical, boris_initial_half_step_backwards_XYZ, boris_initial_half_step_backwards_RZPhi
use mod_penning_case
use fruit
implicit none
contains

subroutine test_penning_case_cartesian_4
  type(case_penning_cartesian)  :: case
  class(particle_base), allocatable :: particle
  real*8, parameter :: dt = 1d-4
  real*8 :: err
  integer :: k, n_steps
  allocate(particle_kinetic_leapfrog::particle)

  n_steps = nint(case%time_end/dt)
  call case%initialize(particle)
  select type (particle)
  type is (particle_kinetic_leapfrog)
    call assert_equals(1, int(particle%q,4), "Charge must be set correctly")
    call boris_initial_half_step_backwards_XYZ(particle, case%mass, case%E(particle%x, 0.d0), &
      case%B(particle%x, 0.d0), dt)
    do k=1,n_steps
      call boris_push_cartesian(particle, case%mass, case%E(particle%x, 0.d0), case%B(particle%x, 0.d0), dt)
    end do
    err = case%calc_error(particle)
  end select
  call assert_equals(0.d0, err, 1.1d-4, "Error must be below 1.1d-4")
end subroutine test_penning_case_cartesian_4
subroutine test_penning_case_cartesian_5
  type(case_penning_cartesian)  :: case
  type(particle_kinetic_leapfrog) :: particle
  real*8, parameter :: dt = 1d-5
  real*8 :: err
  integer :: k, n_steps

  n_steps = nint(case%time_end/dt)
  call case%initialize(particle)
  call assert_equals(1, int(particle%q,4), "Charge must be set correctly")
  call boris_initial_half_step_backwards_XYZ(particle, case%mass, case%E(particle%x, 0.d0), &
    case%B(particle%x, 0.d0), dt)
  do k=1,n_steps
    call boris_push_cartesian(particle, case%mass, case%E(particle%x, 0.d0), case%B(particle%x, 0.d0), dt)
  end do
  err = case%calc_error(particle)
  call assert_equals(0.d0, err, 1.1d-6, "Error must be below 1.1d-6")
end subroutine test_penning_case_cartesian_5


subroutine test_penning_case_cylindrical_4
  type(case_penning_cylindrical)  :: case
  type(particle_kinetic_leapfrog) :: particle
  real*8, parameter :: dt = 1d-4
  real*8 :: err
  integer :: k, n_steps

  n_steps = nint(case%time_end/dt)
  call case%initialize(particle)
  call assert_equals(1, int(particle%q,4), "Charge must be set correctly")
  call boris_initial_half_step_backwards_RZPhi(particle, case%mass, case%E(particle%x, 0.d0), &
    case%B(particle%x, 0.d0), dt)
  do k=1,n_steps
    call boris_push_cylindrical(particle, case%mass, case%E(particle%x, 0.d0), case%B(particle%x, 0.d0), dt)
  end do
  err = case%calc_error(particle)
  call assert_equals(0.d0, err, 1.1d-4, "Error must be below 1.1d-4")
end subroutine test_penning_case_cylindrical_4

subroutine test_penning_case_cylindrical_5
  type(case_penning_cylindrical)  :: case
  type(particle_kinetic_leapfrog) :: particle
  real*8, parameter :: dt = 1d-5
  real*8 :: err
  integer :: k, n_steps

  n_steps = nint(case%time_end/dt)
  call case%initialize(particle)
  call assert_equals(1, int(particle%q,4), "Charge must be set correctly")
  call boris_initial_half_step_backwards_RZPhi(particle, case%mass, case%E(particle%x, 0.d0), &
    case%B(particle%x, 0.d0), dt)
  do k=1,n_steps
    call boris_push_cylindrical(particle, case%mass, case%E(particle%x, 0.d0), case%B(particle%x, 0.d0), dt)
  end do
  err = case%calc_error(particle)
  call assert_equals(0.d0, err, 1.1d-6, "Error must be below 1.1d-6")
end subroutine test_penning_case_cylindrical_5
end module pusher_boris_penning_spec
