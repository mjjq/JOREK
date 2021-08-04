!< Radiation reaction force a.k.a. synchrotron loss operator
!!
!! This module contains operators for relativistic particle and guiding center
!! to take into account the radiation reaction force. This is only relevant
!! for relativistic electrons.
!!
!! The implemented operators are taken from here:
!! https://arxiv.org/pdf/1412.1966.pdf
module mod_radreactforce
  use constants, only: EL_CHG, SPEED_OF_LIGHT, ATOMIC_MASS_UNIT, PI, EPS_ZERO, MASS_ELECTRON
  use mod_fields, only: fields_base
  use mod_particle_types

  implicit none
  
contains

  !< Calculate the characteristic time for the radiation reaction force
  !!
  !! The characteristic time is calculated assuming test particle is an
  !! electron.
  function radreactforce_chartime(B, gamma)
    implicit none
    real*8 :: radreactforce_chartime
    real*8, intent(in) :: B     !< Magnetic field strength [T]
    real*8, intent(in) :: gamma !< Particle Lorentz factor

    radreactforce_chartime = 6 * PI * EPS_ZERO * gamma * ( MASS_ELECTRON * SPEED_OF_LIGHT )**3 &
                           / ( EL_CHG**4 * B**2 )

  end function radreactforce_chartime

  !< Apply radiation reaction force to a particle
  !!
  !! Radiation reaction force is applied with forward Euler method. The term
  !! containing explicit magnetic field time-dependence is omitted for now.
  subroutine radreactforce_kinetic(fields, t, dt, mass, particle)

    use mod_coordinate_transforms, only: vector_cylindrical_to_cartesian

    implicit none
    
    class(fields_base), intent(in) :: fields !< The field data
    real(kind=8), intent(in)       :: t      !< Current time
    real(kind=8), intent(in)       :: dt     !< Time step
    real(kind=8), intent(in)       :: mass   !< Particle mass [AMU]
    type(particle_kinetic_relativistic), intent(inout) :: particle !< Particle to be operated

    real*8 :: tau
    real*8 :: pperp(3), p(3) !< Perpendicular and total momentum in cartesian coordinates

    real*8 :: E(3), B(3), Bxyz(3), psi, U
    real*8 :: tau, Bnorm, ppar, mu, gamma

    ! Convert to SI units
    mass = mass * ATOMIC_MASS_UNIT
    p    = particle%p * ATOMIC_MASS_UNIT

    ! Calculate the characteristic time and pperp
    call fields%calc_EBpsiU(t, particle%i_elm, particle%st, particle%x(3), E, B, psi, U)

    Bxyz  = vector_cylindrical_to_cartesian(particle%x(3), B)
    Bnorm = norm2(B)
    gamma = sqrt(1.0 + ( norm2(p) / ( mass * SPEED_OF_LIGHT ) )**2 )
    tau   = radreactforce_chartime(B, gamma)
    pperp = p * ( 1.0 - dot_product(Bxyz,p) / ( Bnorm * norm2(p) ) )

    ! Apply RR-force and convert back to JOREK units
    particle%p = -dt * ( pperp + p * norm2(pperp)**2 / ( mass * SPEED_OF_LIGHT )**2 ) / tau
    particle%p = particle%p / ATOMIC_MASS_UNIT

  end subroutine radreactforce_kinetic


  !< Apply radiation reaction force to a guiding center
  !!
  !! Radiation reaction force is applied with forward Euler method. Only the zeroth
  !! order terms are included for now (magnetic field non-uniformity is neglegted).
  subroutine radreactforce_gc(fields, t, dt, mass, particle)
    implicit none

    class(fields_base), intent(in) :: fields !< The field data
    real(kind=8), intent(in)       :: t      !< Current time
    real(kind=8), intent(in)       :: dt     !< Time step
    real(kind=8), intent(in)       :: mass   !< Particle mass [AMU]
    type(particle_gc_relativistic), intent(inout) :: particle !< Particle to be operated
    
    real*8 :: E(3), B(3), psi, U
    real*8 :: tau, Bnorm, ppar, mu, gamma

    ! Convert to SI units
    mass = mass * ATOMIC_MASS_UNIT
    ppar = particle%p(1) * ATOMIC_MASS_UNIT
    mu   = particle%p(2) * ATOMIC_MASS_UNIT

    ! Calculate the characteristic time
    call fields%calc_EBpsiU(t, particle%i_elm, particle%st, particle%x(3), E, B, psi, U)

    Bnorm = norm2(B)
    gamma = sqrt( 1.0 + ppar**2 / ( mass * SPEED_OF_LIGHT )**2 + 2 * Bnorm * mass * mu )
    tau   = radreactforce_chartime(Bnorm, gamma)

    ! Apply RR-force
    particle%p(1) = ppar - dt * 2 * ppar * mu * Bnorm / ( mass * SPEED_OF_LIGHT**2 * tau )
    particle%p(2) = mu - dt * 2 * mu * ( 1.0 + 2 * mu * Bnorm / ( mass * SPEED_OF_LIGHT**2 ) ) / tau

    ! Convert back to JOREK units
    particle%p(1) = particle%p(1) / ATOMIC_MASS_UNIT
    particle%p(2) = particle%p(2) / ATOMIC_MASS_UNIT

  end subroutine radreactforce_gc

end module mod_radreactforce
