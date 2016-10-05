!> Particle pusher module with the Boris scheme
module mod_boris
  use mod_particle_types
  implicit none
  private

  public boris_push_cylindrical, boris_push_cartesian, boris_initial_half_step_backwards
contains

!> Push a single particle for some timesteps with the boris method
!> See G.L. Delzanno, E. Camporeale / JCP 253 (2013) 259-277 for details
pure subroutine boris_push_cylindrical(particle, E, B, dt)
  type(particle_kinetic_leapfrog), intent(inout)  :: particle
  real*8, dimension(3), intent(in) :: E, B
  real*8, intent(in) :: dt
  real*8 :: R, Rphi

  ! update the velocity from v^(n-1/2) to v^(n+1/2)
  call boris_method_v_only(particle, E, B, dt)
  ! update the position from v^n to v^(n+1)
  ! Calculate the new R and RPhi
  R    = particle%x(1) + particle%v(1) * dt
  RPhi = particle%v(2) * dt

  ! Calculate the new R, Phi, Z
  particle%x(1) = sqrt(R**2 + RPhi**2)
  particle%x(2) = particle%x(2) + asin(RPhi / particle%x(1))
  particle%x(3) = particle%x(3) + dt * particle%v(3)

  ! Adjust R and Phi velocities to the new reference frame (z component stays the same)
  particle%v(1:2) = [R     * particle%v(1) + RPhi * particle%v(2), &
                     -RPhi * particle%v(1) + R    * particle%v(2)] / particle%x(1)
end subroutine boris_push_cylindrical

!> Push a single particle for some timesteps with the boris method
!> See G.L. Delzanno, E. Camporeale / JCP 253 (2013) 259-277 for details
pure subroutine boris_push_cartesian(particle, E, B, dt)
  class(particle_kinetic_leapfrog), intent(inout)  :: particle
  real*8, dimension(3), intent(in) :: E, B
  real*8, intent(in) :: dt

  ! update the velocity from v^(n-1/2) to v^(n+1/2)
  call boris_method_v_only(particle, E, B, dt)
  particle%x = particle%x + particle%v * dt
end subroutine boris_push_cartesian



pure subroutine boris_method_v_only(particle, E, B, dt)
  use mod_constants, only: EL_CHG, ATOMIC_MASS_UNIT
  real*8, dimension(3), intent(in) :: E, B
  class(particle_kinetic_leapfrog), intent(inout) :: particle
  real*8, intent(in) :: dt

  real*8  :: fE, fB, eom
  real*8  :: psi, U, B2, Bnorm

  eom = EL_CHG / (particle%m * ATOMIC_MASS_UNIT)

  B2    = dot_product(B,B)
  Bnorm = sqrt(B2)

  ! Calculate the geometric factor f = tan(q/m delta_t/2 |B|)/|B|
  fE = particle%q*eom * dt * 0.5d0
  fB = tan(particle%q*eom * dt * 0.5d0 * Bnorm) / Bnorm

  ! Calculate the electric field update (v^n-1/2 -> v-) with the Boris method
  particle%v = particle%v + fE * E
  ! Calculate the rotation
  particle%v = (particle%v + 2.d0*fB/(1.d0+fB*fB*B2)*( &
    left_handed_cross_product(particle%v,B) &
    - fB * particle%v * B2 &
    + fB * B * dot_product(particle%v,B)))
  ! Calculate the next electric field update (v+ -> v^n+1/2)
  particle%v = particle%v + fE * E
end subroutine boris_method_v_only

!> Given a particle with position x and velocity v at time t=0 (t^0), calculate v^-1/2
!> in cylindrical coordinates this is the same as in cartesian coordinates
!> (see G.L. Delzanno, E. Camporeale / JCP 253 (2013) 259-277 for details)
pure subroutine boris_initial_half_step_backwards(particle, E, B, dt)
  use mod_constants, only: EL_CHG, ATOMIC_MASS_UNIT
  class(particle_kinetic_leapfrog), intent(inout) :: particle
  real*8, dimension(3), intent(in) :: E, B 
  real*8, intent(in) :: dt
  real*8, dimension(3) :: v !< for calculating the initial half-step
  real*8 :: f, B2
  f = - (EL_CHG * real(particle%q)) / (ATOMIC_MASS_UNIT * particle%m) * dt * 0.25d0
  B2 = dot_product(B, B)
  v = particle%v + f*E
  v = (v + 2.d0*f/(1.d0+f**2*B2) &
      * (left_handed_cross_product(v, B) - f*v*B2 + f*B*dot_product(v,B)))
  v = v + f*E
  particle%v = v
end subroutine boris_initial_half_step_backwards

!> The cross product in a left-handed coordinate system (e.g. XYZ or RPhiZ)
pure function left_handed_cross_product(a, b)
  real*8, dimension(3) :: left_handed_cross_product
  real*8, dimension(3), intent(in) :: a, b

  left_handed_cross_product(1) = a(2) * b(3) - a(3) * b(2)
  left_handed_cross_product(2) = a(3) * b(1) - a(1) * b(3)
  left_handed_cross_product(3) = a(1) * b(2) - a(2) * b(1)
end function left_handed_cross_product
end module mod_boris
