!> Mod_particle_type contains the default particle type, type_particle
!> If you wish to add a new particle type, create one here, inheriting from
!> [[particle_base]], and update [[mod_particle_io]] (search for `particle_kinetic`
!> and add your particle at each spot)
module mod_particle_types
  implicit none
  private
  public particle_base, particle_kinetic, particle_kinetic_leapfrog, particle_gc, particle_fieldline, copy_particle_base

  !> The base type for all other particles. Includes only the position and weight elements
  !> Integration in a 2D finite element method is included in the form of 2 coordinates
  !> and an element index.
  type, abstract :: particle_base
    real*8    :: x(3)             !< particle position in real space
    real*8    :: st(2)            !< JOREK integration: particle position in the element
    real*4    :: weight = 1.0     !< weight (i.e. number of particles)
    integer*4 :: i_elm            !< JOREK integration: index in element_list
  end type particle_base

  !> A simple type just for fieldline tracing
  type, extends(particle_base) :: particle_fieldline
  end type particle_fieldline

  !> A simple guiding-center particle type.
  type, extends(particle_base) :: particle_gc
    real*8    :: E !< The particle energy
    real*8    :: mu !< The magnetic moment
    integer*1 :: q !< Charge [e]
  end type particle_gc

  !> For most kinetic methods the velocity is required at time \(t\)
  type, extends(particle_base) :: particle_kinetic
    real*8, dimension(3) :: v !< Velocity [m/s]
    integer*1            :: q !< charge [e]
  end type particle_kinetic

  !> Leapfrog methods define the particle velocity at time \(t^{n-1/2}\)
  !> and are therefore incompatible with normal kinetic methods (but a conversion
  !> function should not be too difficult)
  type, extends(particle_base) :: particle_kinetic_leapfrog
    real*8, dimension(3) :: v !< Velocity [m/s] at t=t^(n-1/2) (where the position is known at t^n)
    integer*1            :: q !< charge [e]
  end type particle_kinetic_leapfrog

contains
  !> Copy the base variables from one particle of class(particle_base) to another
  pure subroutine copy_particle_base(in, out)
    class(particle_base), intent(in)    :: in
    class(particle_base), intent(inout) :: out
    out%x      = in%x
    out%st     = in%st
    out%weight = in%weight
    out%i_elm  = in%i_elm
  end subroutine copy_particle_base
end module mod_particle_types
