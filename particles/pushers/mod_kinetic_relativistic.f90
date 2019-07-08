!> Paritcle pusher module for integrating full orbits of relativistic
!> particles. Available integrators:
!> 5-steps Volume Preserving Integrators: R. Zhang et all, Phys. of Plasmas,
!> vol.22, p.044501 2015
module mod_kinetic_relativistic
use mod_particle_types
! use electric charge, atomic mass unit and speed of ligth
use constants, only: EL_CHG,ATOMIC_MASS_UNIT,SPEED_OF_LIGTH

implicit none

private

public relativistic_volume_preserving_push_cartesian

contains

!---------------------------------------------------------------------------
!> This subroutine implements a test version of the volume preserving 
!> algorithm: described in R. Zhang, Phys. of Plasmas, vol.22, p.044501, 2015.
!> using the complete Cayley transform in cartesian coordinates.
!> This subroutine has to be used only for tests and comparisons not for
!> production work. WARNING: this pusher works only for constant and
!> uniform E and B.
!> inputs:
!>   particle: (particle_kinetic_relativistic) relativistic particle type
!>   m:        (real8) particle mass in [AMU]
!>   E:        (real8)(3) electric field in [V/m]
!>   B:        (real8)(3) magnetic field in [T]
!>   dt:       (real8) time step in [s]
!> outputs:
!>   particle: (particle_kinetic_relativistic) relativistic particle type
pure subroutine volume_preserving_push_cartesian(particle,m,E,B,dt)
  use mod_pusher_tools, only: cayley_transform !< use full Cayley transform
  ! define input output variables
  class(particle_kinetic_relativistic), intent(inout) :: particle !< relativistic particle
  ! define input variables
  real(kind=8), intent(in) :: m, dt !< mass and time step
  real(kind=8), dimension(3), intent(in) :: E,B !< electric and magnetic field
  ! internal variable
  real(kind=8) :: scaling_factor !< scaling factor [s^2*C/(kg*m)]

  ! compute the dimensional q
  scaling_factor = 5.d-1*dt*particle%q*EL_CHG/(ATOMIC_MASS_UNIT*mass*SPEED_OF_LIGTH)

  ! compute dimensionless momenta
  particle%p = particle%p/(mass*SPEED_OF_LIGTH)

  ! compute position at t_(i+1/2)
  particle%x = particle%x + (5.d-1*dt*SPEED_OF_LIGTH*particle%p)/&
    (sqrt(1.d0+dot_product(particle%p,particle%p)))
  
  ! compute the momenta at t_(i+1/2)
  particle%p = particle%p + scaling_factor*E
  
  ! rotate the momenta with respect to the magnetic field
  particle%p = matmul(cayley_transform(SPEED_OF_LIGTH*scaling_factor/&
    (sqrt(1.d0+dot_product(particle%p,particle%p))),&
    B),particle%p)

  ! update compute momenta at t_(i+1)
  particle%p = particle%p + scaling_factor*E

  ! update particle position at t_(i+1)
  particle%x = particle%x + (5.d-1*dt*SPEED_OF_LIGTH*particle%p)/&
    (sqrt(1.d0+dot_product(particle%p,particle%p)))

  ! compute dimensional momenta
  particle%p = particle%p*mass*SPEED_OF_LIGTH
  
end subroutine volume_preserving_push_cartesian

!---------------------------------------------------------------------------

end module mod_kinetic_relativistic
