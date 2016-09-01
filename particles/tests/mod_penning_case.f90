!> {!tests/mod_penning_case.md!}
module mod_penning_case
  use mod_constants, only: CARTESIAN, CYLINDRICAL, ATOMIC_MASS_UNIT, EL_CHG
  use mod_coordinate_transforms
  use mod_case
  use mod_pusher
  use mod_boris
  implicit none

  type, extends(case), public :: case_penning
    contains
      procedure :: initialize_particle => initialize_particle_penning
      procedure :: calculate_error => calculate_error_norm
  end type
  interface case_penning
    module procedure new_case_penning
  end interface case_penning

  public :: penning_trajectory
  ! Penning trap parameters (in SI units)
  real*8, parameter :: omega_e = 4.9d0 !< rad/s
  real*8, parameter :: omega_b = 25.d0 !< rad/s
  real*8, parameter :: epsilon = -1.d0
  real*8, parameter :: x0(3)   = [10.d0,0.d0,0.d0] !< m (xyz)
  real*8, parameter :: v0(3)   = [50.d0,0.d0,20.d0] !< m (xyz)
  real*4, parameter :: mass    = 1.d0 !< atomic mass units
  integer*1, parameter :: charge = 1 !< electron charges
  real*8, parameter :: time_end = 16.d0 !< s
  private
contains

!> Interface exists because we need to set prescribed_fields to the right value.
!> it is important to use this! otherwise you need to set the fields manually.
pure function new_case_penning(geometry) result(new)
  type(case_penning) :: new
  integer*1, intent(in) :: geometry !< one of [[mod_constants:CARTESIAN]] or [[mod_constants:CYLINDRICAL]]
  new%time_end = time_end
  if (geometry .eq. CARTESIAN)   new%fields = prescribed_fields(geometry, E_cartesian, B_z)
  if (geometry .eq. CYLINDRICAL) new%fields = prescribed_fields(geometry, E_cylindrical, B_z)
end function new_case_penning

!> Wrapper to the global initialization routine with the parameters for this case.
subroutine initialize_particle_penning(this, particle, pusher)
  class(case_penning), intent(in)     :: this
  class(particle_base), intent(inout) :: particle
  class(pusher_base), intent(in)      :: pusher
  call this%initialize_particle_all(particle, pusher, x0, v0, charge, mass)
end subroutine initialize_particle_penning


!> Calculate the error as the difference between particle position vectors
pure function calculate_error_norm(this, particle) result(err)
  class(case_penning), intent(in)  :: this
  class(particle_base), intent(in) :: particle
  real*8 :: err
  if (this%fields%geometry .eq. CARTESIAN)   err = norm2(penning_trajectory(this%time_end) - particle%x)
  if (this%fields%geometry .eq. CYLINDRICAL) err = norm2(penning_trajectory(this%time_end) - cylindrical_to_cartesian(particle%x))
end function calculate_error_norm


!> Magnetic field in the penning trap (valid for both cylindrical and cartesian
!> cases because there is only a z-component)
pure function B_z(x, t) result(B)
  real*8, dimension(3), intent(in) :: x
  real*8, intent(in) :: t
  real*8, dimension(3) :: B
  B = [0.d0, 0.d0, 1.d0]*omega_b*mass*ATOMIC_MASS_UNIT/(real(charge)*EL_CHG)
end function B_z

pure function E_cartesian(x, t) result(E)
  real*8, dimension(3), intent(in) :: x
  real*8, intent(in) :: t
  real*8, dimension(3) :: E
  E = epsilon*omega_e**2/(real(charge)*el_chg)*mass*atomic_mass_unit * &
      [-x(1), -x(2), 2.d0*x(3)]
end function E_cartesian
pure function E_cylindrical(x, t) result(E)
  real*8, dimension(3), intent(in) :: x
  real*8, intent(in) :: t
  real*8, dimension(3) :: E
  E = epsilon*omega_e**2/(real(charge)*el_chg)*mass*atomic_mass_unit * &
      [-x(1), 0.d0, 2.d0*x(3)]
end function E_cylindrical

!> Calculate the position of a particle in the penning trap, released at 
!> \(x_0\) with speed \(v_0\), at time \(t\)
pure function penning_trajectory(t) result(x)
  real*8, intent(in) :: t !< The time at which to calculate the solution value
  real*8             :: x(3) !< The position of the particle at time \(t\) in cartesian coordinates

  ! Internal variables
  real*8 :: omega_plus, omega_minus
  real*8 :: R_plus, R_minus
  real*8 :: T_plus, T_minus
  real*8 :: omega
  complex(kind=8) :: w

  ! Some initialization
  omega = sqrt(-2.d0*epsilon)*omega_e
  omega_plus  = 0.5d0*(omega_b + sqrt(omega_b**2 + 4.d0*epsilon*omega_e**2))
  omega_minus = 0.5d0*(omega_b - sqrt(omega_b**2 + 4.d0*epsilon*omega_e**2))
  R_minus = (omega_plus * x0(1) + v0(2))/(omega_plus - omega_minus)
  R_plus  = x0(1) - R_minus
  T_minus = (omega_plus * x0(2) - v0(1))/(omega_plus - omega_minus)
  T_plus  = x0(2) - T_minus

  ! Calculate the result in the x-y plane in terms of w = x + iy
  w = cmplx(R_plus ,T_plus , 8)*exp(cmplx(0.d0,-omega_plus*t, 8)) + &
      cmplx(R_minus,T_minus, 8)*exp(cmplx(0.d0,-omega_minus*t, 8))

  x(1) = real(real(w),8)
  x(2) = real(aimag(w),8)
  x(3) = x0(3)*cos(omega*t)+v0(3)*sin(omega*t)/omega
end function penning_trajectory
end module mod_penning_case
