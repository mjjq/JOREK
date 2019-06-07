!> Mod_particle_type contains the default particle type, type_particle
!> If you wish to add a new particle type, create one here, inheriting from
!> [[particle_base]], and update [[mod_particle_io]] (search for `particle_kinetic`
!> and add your particle at each spot)
module mod_particle_types
  implicit none
  private
  public particle_base, particle_kinetic, particle_kinetic_leapfrog, particle_gc, particle_fieldline, copy_particle_base, copy_particle_derived
  public particle_get_q

  !> The base type for all other particles. Includes only the position and weight elements
  !> Integration in a 2D finite element method is included in the form of 2 coordinates
  !> and an element index.
  type, abstract :: particle_base
    real*8    :: x(3)             !< particle position in real space
    real*8    :: st(2)            !< particle position in the element
    real*4    :: weight = 1.0     !< weight (i.e. number of particles)
    integer*4 :: i_elm            !< index in element_list
  end type particle_base

  !> A simple type just for fieldline tracing in two-step methods (Adams Bashforth) or for forward euler
  type, extends(particle_base) :: particle_fieldline
    real*8    :: B_hat_prev(3) !< Field direction at previous timestep
    real*8    :: v !< Parallel velocity along the fieldline
  end type particle_fieldline

  !> A simple guiding-center particle type.
  type, extends(particle_base) :: particle_gc
    real*8    :: E !< The particle energy [eV]
    real*8    :: mu !< The magnetic moment [eV/T]. Sign determines sign of v_par
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

  !> Convenience function to obtain q if it exists, or 0 otherwise
  !> Here also because of https://gcc.gnu.org/bugzilla/show_bug.cgi?id=82064
  !> which means that we cannot use the same derived type in too many modules which will be
  !> imported in the main program (roughly)
  pure function particle_get_q(in) result(q)
    class(particle_base), intent(in) :: in
    integer*1 :: q
    select type (p => in)
    type is (particle_kinetic)
      q = p%q
    type is (particle_kinetic_leapfrog)
      q = p%q
    type is (particle_gc)
      q = p%q
    class default
      q = 0
    end select
  end function particle_get_q

  !> This subroutine copies derived particles type.
  !> we do not use a different subroutine for each derived particle type
  !> because fortran for avoiding problems when type selectors are given 
  !> in functions / subroutines within OpenMP regions
  !> inputs:
  !>   in: (type particle) particle to be copied
  !> outputs:
  !>   out: (type particle) copy destination
  !> Author: C. Sommariva, 07/06/2019, email: cristian.sommariva[at]epfl.ch
  pure subroutine copy_particle_derived(in,out)
    !< definitions
    class(particle_base),intent(in) :: in !< define the particle to be copied
    class(particle_base),intent(out) :: out !< define the destination particle

    !< actions
    call copy_particle_base(in, out) !< first, copy the particle base
    !< select the "in" particle type
    select type (p_in => in)
    type is (particle_fieldline) !< in particle is particle_fieldline type
      select type (p_out => out) !< the particle out type is selected
          type is (particle_fieldline) !< in particle is particle_fieldline type
          p_out%B_hat_prev = p_in%B_hat_prev !< copy the magnetic field
          p_out%v = p_in%v                   !< copy the parallel velocity along a magnetic field line
      end select !< end the particle "out" select
    type is (particle_gc)        !< in particle is particle_gc type
      select type (p_out => out) !< the particle out type is selected
        type is (particle_gc)    !< in particle is particle_gc type
        p_out%E = p_in%E   !< copy the energy
        p_out%mu = p_in%mu !< copy the magnetic moment
        p_out%q = p_in%q   !< copy the charge
      end select !< end the particle "out" select
    type is (particle_kinetic) !< in particle is particle_kinetic type
      select type (p_out => out)   !< the particle out type is selected
        type is (particle_kinetic) !< in particle is particle_kinetic type
        p_out%v = p_in%v !< copy the velocity
        p_out%q = p_in%q !< copy the charge
      end select!< end the particle "out" select
    type is (particle_kinetic_leapfrog) !< in particle is particle_kinetic type
      select type (p_out => out)            !< the particle out type is selected
        type is (particle_kinetic_leapfrog) !< in particle is particle_kinetic type
        p_out%v = p_in%v !< copy the velocity
        p_out%q = p_in%q !< copy charge
      end select !< end the particle "out" select
    end select !< end the particle "in" select
  end subroutine copy_particle_derived
end module mod_particle_types
