!> Calculate an approximation to the guiding center position of a particle
!> being time-stepped with the Boris method (or another leapfrog method)
pure function guiding_center_position(particle, dt) result(x_gc)

use mod_particles
use constants
use data_structure
use phys_module, only: central_density, central_mass

implicit none

! Routine parameters
type (type_particle), intent(in)  :: particle !< The particle to calculate the guiding centre position of
real*8,               intent(in)  :: dt !< The timestep size (used to calculate v at the current time)
real*8 :: x_gc(3) !< The guiding centre position of the particle

! Internal parameters
real*8 :: E(3), B(3), psi, U, B2, fE, fB, qom, t_norm
real*8 :: v(3), R_new, R_update, RPhi_update, v_tmp(3), v_n(3)

!! Obtain the B and E fields
call calc_EB(particle%i_elm,particle%st,particle%x(3),E,B,psi,U)
B2 = dot_product(B,B)

! q/m in JOREK units, including correction for E and B having semi-SI units
t_norm = sqrt(mu_zero * mass_proton * central_mass * central_density * 1.d20)    ! jorek time normalisation
qom = real(particle%q) * EL_CHG / (particle%mass * ATOMIC_MASS_UNIT) * t_norm

!! Calculate the velocity at t_n+1/2 (boris method step, velocity only)
! Calculate the geometric factor f = tan(q/m delta_t/2 |B|)/|B|
fB = tan(qom * dt * 0.5d0 * sqrt(B2)) / sqrt(B2)
fE = qom * dt * 0.5d0

! Calculate the electric field update (v^n-1/2 -> v-) with the Boris method
v = particle%v + fE * E
! Calculate the rotation
v = (v + 2.d0*fB/(1.d0+fB*fB*B2)*( cross_product(v,B) &
  - fB * v * B2 &
  + fB * B * dot_product(v,B)))
! Calculate the next electric field update (v+ -> v^n+1/2)
v = v + fE * E

R_update = particle%x(1) + v(1) * dt
RPhi_update = v(3) * dt

! Calculate the new R and Z
R_new = sqrt(R_update**2 + RPhi_update**2)

! Adjust R and Phi velocities to the new reference frame (z stays the same)
v_tmp = v
v(1) =  R_update/R_new    * v_tmp(1) + RPhi_update/R_new * v_tmp(3)
v(3) = -RPhi_update/R_new * v_tmp(1) + R_update/R_new    * v_tmp(3)


!! Calculate the gyroradius and direction of the guiding center
v_n = (v + particle%v)/2 ! Calculate v at t^n as the average of v^n-1/2 and v^n+1/2
x_gc = particle%x + cross_product(v_n,B) / qom / B2

end function guiding_center_position
