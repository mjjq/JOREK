!> {!particles/benchmarks/pusher/mod_penning_case.md!}
module mod_penning_case
  use constants, only: ATOMIC_MASS_UNIT, EL_CHG
  use mod_coordinate_transforms
  use mod_case
  use mod_boris ! Is due to a gfortran bug with derived types and modules
  use mod_particle_types
  implicit none

  private
  public :: penning_trajectory
  public :: case_penning_cartesian, case_penning_cylindrical
  public :: omega_e, omega_b, epsilon, charge, mass
  public :: x0, v0
  public :: time_end

  ! Penning trap parameters (in SI units)
  real*8, parameter :: omega_e = 4.9d0 !< rad/s
  real*8, parameter :: omega_b = 25.d0 !< rad/s
  real*8, parameter :: epsilon = -1.d0
  real*8, parameter :: x0(3)   = [10.d0,0.d0,0.d0] !< m (xyz), (RZphi)
  real*8, parameter :: v0(3)   = [50.d0,0.d0,20.d0] !< m (xyz), [50,20,0] in RZphi
  real*4, parameter :: mass    = 1.d0 !< atomic mass units
  integer*1, parameter :: charge = 1 !< charge units
  real*8, parameter :: time_end = 16.d0 !< s
  ! The particle remains between +- 3 in Z and 8 and 13 in R for these parameters. set this in an input file

  !> See [[mod_penning_case]] for a description of the testcase.
  type, extends(case), abstract :: case_penning
    real*8 :: time_end = time_end
    real*8 :: mass = mass
    contains
      procedure :: initialize => initialize_particle_penning
      procedure :: calc_error => calculate_error_norm
  end type
  type, extends(case_penning) :: case_penning_cartesian
    contains
      procedure, nopass :: E => E_cartesian
      procedure, nopass :: B => B_z
  end type
  !> See [[mod_penning_case]] for a description of the testcase.
  type, extends(case_penning) :: case_penning_cylindrical
    contains
      procedure, nopass :: E => E_cylindrical
      procedure, nopass :: B => B_z_cyl
  end type
contains

subroutine initialize_particle_penning(this, particle)
  use mod_particle_types
  class(case_penning), intent(in)   :: this
  class(particle_base), intent(inout) :: particle
  select type (this)
  type is (case_penning_cartesian)
    particle%x = x0
    select type (p => particle)
    type is (particle_kinetic_leapfrog)
      p%q = charge
      p%v = v0
      write(*,*) "OK"
    class default
      write(*,*) "WARNING: Unknown type in initialize_particle_penning"
    end select
  type is (case_penning_cylindrical)
    particle%x = cartesian_to_cylindrical(x0)
    select type (p => particle)
    type is (particle_kinetic_leapfrog)
      p%q = charge
      p%v = vector_rotation(cartesian_to_cylindrical(v0), p%x(3))
    class default
      write(*,*) "WARNING: Unknown type in initialize_particle_penning"
    end select
  end select
end subroutine

!> Calculate the error as the difference between particle position vectors
pure function calculate_error_norm(this, particle) result(err)
  class(case_penning), intent(in)  :: this
  class(particle_base), intent(in) :: particle
  real*8 :: err
  err = 0.d0
  select type (this)
  type is (case_penning_cartesian)
    err = norm2(penning_trajectory(this%time_end) - particle%x)
  type is (case_penning_cylindrical)
    err = norm2(penning_trajectory(this%time_end) - cylindrical_to_cartesian(particle%x))
  end select
end function calculate_error_norm


!> Magnetic field in the penning trap
pure function B_z(x, t) result(B)
  real*8, dimension(3), intent(in) :: x
  real*8, intent(in) :: t
  real*8, dimension(3) :: B
  B = [0.d0, 0.d0, 1.d0]*omega_b*mass*ATOMIC_MASS_UNIT/(real(charge)*EL_CHG)
end function B_z
!> Magnetic field in the penning trap in RZPhi coordinates
pure function B_z_cyl(x, t) result(B)
  real*8, dimension(3), intent(in) :: x
  real*8, intent(in) :: t
  real*8, dimension(3) :: B
  B = [0.d0, 1.d0, 0.d0]*omega_b*mass*ATOMIC_MASS_UNIT/(real(charge)*EL_CHG)
end function B_z_cyl

!> Electric field in the penning trap
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
      [-x(1), 2.d0*x(2), 0.d0]
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
