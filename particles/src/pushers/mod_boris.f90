!> Particle pusher module with the Boris scheme
!> See [[mod_pusher_no_action]] for an annotated example of the components below.
module mod_boris
  use mod_pusher
  use mod_particle_types

  type, extends(pusher_base) :: pusher_boris
  contains
    procedure :: push_single => boris_push_single_particle
    procedure, private :: boris_method_v_only => boris_method_v_only
    procedure, private :: boris_method_cylindrical_correction => boris_method_cylindrical_correction
    procedure :: initial_half_step_backwards
  end type pusher_boris
  interface pusher_boris
    module procedure new_pusher_boris
  end interface pusher_boris

  public :: pusher_boris
  public :: particle_base
  private
contains

!> Workaround for https://gcc.gnu.org/bugzilla/show_bug.cgi?id=77412
pure function new_pusher_boris(groups, fixed_timestep, hooks)
  type(pusher_boris) :: new_pusher_boris
  integer, dimension(:), intent(in), optional         :: groups
  real*8, intent(in), optional                        :: fixed_timestep
  type(hook_base), intent(in), optional, dimension(:) :: hooks
  if (present(groups)) allocate(new_pusher_boris%groups, source=groups)
  if (present(fixed_timestep)) allocate(new_pusher_boris%fixed_timestep, source=fixed_timestep)
  if (present(hooks)) allocate(new_pusher_boris%hooks, source=hooks)
end function new_pusher_boris

!> Push a single particle for some timesteps with the boris method
pure subroutine boris_push_single_particle(this, fields, particle, time_start, time_end)
  use mod_fields
  use mod_constants, only: TICK, EL_CHG, ATOMIC_MASS_UNIT
  class(pusher_boris), intent(in)      :: this
  class(fields_base), intent(in)       :: fields
  class(particle_base), intent(inout)  :: particle
  real*8, intent(in) :: time_start, time_end
  real*8 :: eom
  integer :: j, k, n_step

  ! Skip this particle if it left the domain
  if (particle%lost) return
  eom = EL_CHG / (particle%m * ATOMIC_MASS_UNIT)
  n_step = ceiling((time_end - time_start - TICK)/this%fixed_timestep)

  select type(particle)
    type is (particle_kinetic_leapfrog)
      do j=1,n_step
        ! update the velocity from v^(n-1/2) to v^(n+1/2)
        call this%boris_method_v_only(fields, particle, eom, time_start + (j-1)*this%fixed_timestep)
        ! update the position from v^n to v^(n+1)
        select case (fields%geometry)
          case (CARTESIAN)
            particle%x = particle%x + particle%v * this%fixed_timestep
          case (CYLINDRICAL)
            call this%boris_method_cylindrical_correction(particle)
        end select
        ! Run hooks if there are any
        if (allocated(this%hooks)) then
          do k=1,size(this%hooks,1)
            call this%hooks(k)%action%do(particle)
          end do
        end if
      end do
  end select
end subroutine boris_push_single_particle


!> The cylindrical velocity correction and position update for the Boris method
!> See G.L. Delzanno, E. Camporeale / JCP 253 (2013) 259-277 for details
pure subroutine boris_method_cylindrical_correction(this, particle)
  class(pusher_boris), intent(in)      :: this
  class(particle_kinetic_leapfrog), intent(inout) :: particle
  real*8 :: R, Rphi

  ! Calculate the new R and RPhi
  R    = particle%x(1) + particle%v(1) * this%fixed_timestep
  RPhi = particle%v(2) * this%fixed_timestep

  ! Calculate the new R, Phi, Z
  particle%x(1) = sqrt(R**2 + RPhi**2)
  particle%x(2) = particle%x(2) + asin(RPhi / particle%x(1))
  particle%x(3) = particle%x(3) + this%fixed_timestep * particle%v(3)

  ! Adjust R and Phi velocities to the new reference frame (z component stays the same)
  particle%v(1:2) = [R     * particle%v(1) + RPhi * particle%v(2), &
                     -RPhi * particle%v(1) + R    * particle%v(2)] / particle%x(1)
end subroutine boris_method_cylindrical_correction


pure subroutine boris_method_v_only(this, fields, particle, eom, t)
  use mod_fields
  class(pusher_boris), intent(in)      :: this
  class(fields_base), intent(in)       :: fields
  class(particle_kinetic_leapfrog), intent(inout) :: particle
  real*8, intent(in) :: eom
  real*8, intent(in) :: t

  real*8  :: fE, fB
  real*8  :: psi, U, B2, Bnorm
  real*8, dimension(3) :: E, B

  call fields%at_particle(particle, t, E, B, psi, U)
  B2    = dot_product(B,B)
  Bnorm = sqrt(B2)

  ! Calculate the geometric factor f = tan(q/m delta_t/2 |B|)/|B|
  fE = particle%q*eom * this%fixed_timestep * 0.5d0
  fB = tan(particle%q*eom * this%fixed_timestep * 0.5d0 * Bnorm) / Bnorm

  ! Calculate the electric field update (v^n-1/2 -> v-) with the Boris method
  particle%v = particle%v + fE * E
  ! Calculate the rotation
  particle%v = (particle%v + 2.d0*fB/(1.d0+fB*fB*B2)*( &
    cross_product(particle%v,B) &
    - fB * particle%v * B2 &
    + fB * B * dot_product(particle%v,B)))
  ! Calculate the next electric field update (v+ -> v^n+1/2)
  particle%v = particle%v + fE * E
  ! above steps are the same for any left-handed coordinate system
end subroutine boris_method_v_only

!> Given a particle with position x and velocity v at time t=0 (t^0), calculate v^-1/2
!> in cylindrical coordinates this is the same as in cartesian coordinates
!> (see G.L. Delzanno, E. Camporeale / JCP 253 (2013) 259-277 for details)
pure subroutine initial_half_step_backwards(this, fields, particle)
  use mod_fields
  use mod_constants, only: EL_CHG, ATOMIC_MASS_UNIT
  class(pusher_boris), intent(in)      :: this
  class(fields_base), intent(in)       :: fields
  class(particle_kinetic_leapfrog), intent(inout) :: particle
  real*8, dimension(3) :: v, E, B !< for calculating the initial half-step
  real*8 :: f, B2
  f = - (EL_CHG * real(particle%q)) / (ATOMIC_MASS_UNIT * particle%m) * this%fixed_timestep * 0.25d0
  call fields%at_particle(particle, 0.d0, E, B)
  B2 = dot_product(B, B)
  v = particle%v + f*E
  v = (v + 2.d0*f/(1.d0+f**2*B2) &
      * (cross_product(v, B) - f*v*B2 + f*B*dot_product(v,B)))
  v = v + f*E
  particle%v = v
end subroutine initial_half_step_backwards

!> The cross product in a left-handed coordinate system (e.g. XYZ or RPhiZ)
pure function cross_product(a, b)
  real*8, dimension(3) :: cross_product
  real*8, dimension(3), intent(in) :: a, b

  cross_product(1) = a(2) * b(3) - a(3) * b(2)
  cross_product(2) = a(3) * b(1) - a(1) * b(3)
  cross_product(3) = a(1) * b(2) - a(2) * b(1)
end function cross_product
end module mod_boris
