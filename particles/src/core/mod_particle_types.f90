!> Mod_particle_type contains the default particle type, type_particle
!> If you wish to add a new particle type, create one here, inheriting from
!> [[particle_base]], and update [[mod_particle_io]] (search for `particle_kinetic`
!> and add your particle at each spot)
module mod_particle_types
  implicit none
  private
  public particle_base, particle_kinetic, particle_kinetic_leapfrog

  !> The base type for all other particles. Includes everything but velocity
  !> as this is different in different pushers.
  !> Integration in a 2D finite element method is included in the form of 2 coordinates
  !> and an element index.
  type, abstract :: particle_base
    real*8    :: x(3)             !< particle position in real space
    real*4    :: m                !< mass [atomic mass units]
    real*4    :: weight = 1.0     !< weight (i.e. number of particles)
    real*8    :: st(2)            !< JOREK integration: particle position in the element
    integer*4 :: i_elm            !< JOREK integration: index in element_list
    integer*2 :: label            !< Particle type number (i in species(i))
    integer*1 :: q                !< charge [e]
    logical*1 :: lost = .false.   !< particle is active or lost
  end type particle_base

  !> For most kinetic methods the velocity is required at time \(t\)
  type, extends(particle_base) :: particle_kinetic
    real*8, dimension(3) :: v !< Velocity [m/s]
  end type particle_kinetic

  !> Leapfrog methods define the particle velocity at time \(t^{n-1/2}\)
  !> and are therefore incompatible with normal kinetic methods (but a conversion
  !> function should not be too difficult)
  type, extends(particle_base) :: particle_kinetic_leapfrog
    real*8, dimension(3) :: v !< Velocity [m/s] at t=t^(n-1/2) (where the position is known at t^n)
  end type particle_kinetic_leapfrog
end module mod_particle_types
