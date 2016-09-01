!# Cases for testing the trajectory of a particle in a strongly inhomogeneous
! magnetic field of the form \[ \mathbf{B} = \left(0, 0, B_0 \exp(\lambda x)\right)\]
! 
! The general solution of the particle trajectory is calculated in [[gradB_solution]]
! 
!## Results
! 
! The trajectory of a particle is followed from \(\mathbf{x}_0\) at \(t=0 s\) to \(t=10 s\).
! An example trajectory, calculated with the [[mod_boris]] method and relatively large timesteps (\(\delta t = 0.01\)) is shown below.
! 
! ![xy-trajectory-boris](|media|/tests/gradB/gradB_xy_boris.png)
!
!### Comparing pushers
! The different pushers in the [[pusher_test]] perform as follows on this case
! ![pusher-test-gradB](|media|/tests/all_pushers/gradB.png)
module mod_gradB_case
  use mod_constants, only: CARTESIAN, CYLINDRICAL, ATOMIC_MASS_UNIT, EL_CHG
  use mod_case
  use mod_boris
  use mod_pusher
  use mod_coordinate_transforms
  implicit none

  !> Case for a penning trap in cartesian coordinates
  type, extends(case), public :: case_gradB
    contains
      procedure :: initialize_particle => initialize_particle_gradB
      procedure :: calculate_error => calculate_error_gradB
  end type
  interface case_gradB
    module procedure new_case_gradB
  end interface case_gradB

  ! gradB parameters
  real*8, parameter :: B0 = 1d0 !< Tesla
  real*8, parameter :: lambda = 1d-1 !< 1/m
  real*8, parameter :: theta_zero = 0.d0 !< radians
  real*8, parameter :: v_perp = 1d0 !< m/s
  real*8, parameter :: x0(3)   = [0d0,0d0,0d0] !< cartesian, meters
  real*8, parameter :: v0(3)   = v_perp*[cos(theta_zero),sin(theta_zero),0.d0] !< cartesian, meters
  real*4, parameter :: mass    = 1.d7 !< atomic mass units
  integer*1, parameter :: charge = 1 !< electron charges
  real*8, parameter :: time_end = 10 !< s
  public :: gradB_solution
  private
contains

!> Interface exists because we need to set prescribed_fields to the right value.
!> it is important to use this! otherwise you need to set the fields manually.
pure function new_case_gradB(geometry) result(new)
  integer*1, intent(in) :: geometry
  type(case_gradB) :: new
  new%time_end = time_end
  if (geometry .eq. CARTESIAN)   new%fields = prescribed_fields(geometry, E_zero, B_cartesian)
  if (geometry .eq. CYLINDRICAL) new%fields = prescribed_fields(geometry, E_zero, B_cylindrical)
end function new_case_gradB

!> Initialize a particle for the gradB test case
subroutine initialize_particle_gradB(this, particle, pusher)
  class(case_gradB), intent(in)       :: this
  class(particle_base), intent(inout) :: particle
  class(pusher_base), intent(in)      :: pusher
  call this%initialize_particle_all(particle, pusher, x0, v0, charge, mass)
end subroutine initialize_particle_gradB

!> Calculate the error as the difference between particle posiition vectors
pure function calculate_error_gradB(this, particle) result(err)
  class(case_gradB), intent(in) :: this
  class(particle_base), intent(in)         :: particle
  real*8 :: err
  if (this%fields%geometry .eq. CARTESIAN)   err = norm2(gradB_solution(this%time_end) - particle%x)
  if (this%fields%geometry .eq. CYLINDRICAL) err = norm2(gradB_solution(this%time_end) - cylindrical_to_cartesian(particle%x))
end function calculate_error_gradB

pure function B_cartesian(x, t) result(B)
  real*8, dimension(3), intent(in) :: x
  real*8, intent(in) :: t
  real*8, dimension(3) :: B
  B = [0.d0, 0.d0, B0*exp(lambda*x(1))]
end function B_cartesian
pure function B_cylindrical(x, t) result(B)
  real*8, dimension(3), intent(in) :: x
  real*8, intent(in) :: t
  real*8, dimension(3) :: B
  B = [0.d0, 0.d0, B0*exp(lambda*x(1)*cos(x(2)))]
end function B_cylindrical
pure function E_zero(x, t) result(E)
  real*8, dimension(3), intent(in) :: x
  real*8, intent(in) :: t
  real*8, dimension(3) :: E
  E = [0.d0, 0.d0, 0.d0]
end function E_zero


!> Particle position at time \(t\) for motion with the parameters defined in [[mod_gradb_case]]
pure function gradB_solution(t) result(x)
  use mod_constants, only: PI, TWOPI
  real*8, intent(in)   :: t !< time [s]
  real*8, dimension(3) :: x !< position vector [m]

  real*8 :: B_at_zero(3)
  real*8 :: alpha, big_gamma, phi_zero, theta, g ! g=small_gamma

  B_at_zero   = B_cartesian(x0, 0.d0)
  big_gamma   = sin(theta_zero) + (charge*EL_CHG)/(mass*ATOMIC_MASS_UNIT)*B_at_zero(3)/(v_perp*lambda)
  g           = sqrt(big_gamma**2-1)
  phi_zero    = atan(1.d0/g * (-1.d0 + big_gamma * tan(theta_zero*0.5d0)))
  alpha       = 0.5d0 * v_perp * lambda * g * t - phi_zero
  theta       = 2.d0 * atan(1.d0/big_gamma - g/big_gamma * tan(alpha)) - TWOPI*real(nint(alpha/PI),8)

  x = x0 + [1.d0/lambda * (log((g**4 + g**2)/(g**2 - g*sin(2.d0*alpha) + cos(2.d0*alpha) + 1.d0)) &
                           - log(big_gamma*(big_gamma - sin(theta_zero)))) &
          , big_gamma * v_perp * t + (theta - theta_zero)/lambda &
          , 0.d0]
end function gradB_solution
end module mod_gradB_case
