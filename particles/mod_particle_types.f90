!> Mod_particle_type contains the default particle type, type_particle
!> If you wish to add a new particle type, create one here, inheriting from
!> [[particle_base]], and update [[mod_particle_io]] (search for `particle_kinetic`
!> and add your particle at each spot)
module mod_particle_types
  implicit none
  private
  public particle_base, particle_kinetic, particle_kinetic_leapfrog, particle_gc, particle_fieldline, copy_particle_base
  public particle_get_q
  !< access to derived type get subroutines
  public get_particle_fieldline_attributes,get_particle_gc_attributes
  !< access to derived type get subroutines
  public get_particle_kinetic_attributes,get_particle_kinetic_leapfrog_attributes
  !< access to derived type set subroutines
  public set_particle_fieldline_attributes,set_particle_gc_attributes
  !< access to derived type set subroutines
  public set_particle_kinetic_attributes,set_particle_kinetic_leapfrog_attributes

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
  !pure subroutine copy_particle_base(in, out)
  subroutine copy_particle_base(in, out)
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

  !> Get attribute particle_fieldline from particle_base datatype. This data are copied 
  !> into a particle_fieldline type.
  !> inputs:
  !>   particle_in:   (particle_base) particle from which data are copied
  !> ouputs:
  !>    particle_out: (particle_fieldline) particle in which data are stored
  !> Author: C. Sommariva, 11/06/2019, email: cristian.sommariva[at]epfl.ch
  pure subroutine get_particle_fieldline_attributes(particle_in,particle_out)
    class(particle_base),intent(in) :: particle_in          !< define input particle as particle base
    class(particle_fieldline),intent(inout) :: particle_out !< define particle output as particle_fieldline

    particle_out%B_hat_prev = (/0.d0,0.d0,0.d0/)!< initialise magnetic field to zero
    particle_out%v = 0.d0                       !< initialise parallel velocity to zero

    select type (p_in => particle_in) !< select data type for particle_in
    type is (particle_fieldline)      !< define the type as field line for particle_in
      particle_out%B_hat_prev = p_in%B_hat_prev !< copy magnetic field in particle out 
      particle_out%v = p_in%v                   !< copy the field line parallel velocity in particle out 
    end select                        !< end select data type for particle_in
  end subroutine get_particle_fieldline_attributes !< end of subroutine get_particle_fieldline_attributes

  !> Set attribute from particle_fieldline in particle_base datatype.
  !> inputs:
  !>    particle_in:  (particle_fieldline) particle from which data are copied
  !> outputs:
  !>    particle_out: (particle_base) particle in which data are stored
  !> Author: C. Sommariva, 11/06/2019, email: cristian.sommariva[at]epfl.ch
  pure subroutine set_particle_fieldline_attributes(particle_in,particle_out)
    class(particle_fieldline),intent(in)  :: particle_in !< define particle input as particle_fieldline
    class(particle_base),intent(inout) :: particle_out   !< define output particle as particle base

    select type (p_out => particle_out) !< select data type for particle_out
    type is (particle_fieldline)        !< define the type as field line for particle_out
      p_out%B_hat_prev = particle_in%B_hat_prev !< copy magnetic field in particle out
      p_out%v = particle_in%v                   !< copy the field line parallel velocity in particle out
    end select                          !< end select data type for particle_out
  end subroutine set_particle_fieldline_attributes !< end of subroutine set_particle_fieldline_attributes


  !> Get attribute particle_gc from particle_base datatype. This data are copied 
  !> into a particle_gc type.
  !> inputs:
  !>   particle_in:   (particle_base) particle from which data are copied
  !> ouputs:
  !>    particle_out: (particle_gc) particle in which data are stored
  !> Author: C. Sommariva, 11/06/2019, email: cristian.sommariva[at]epfl.ch
  pure subroutine get_particle_gc_attributes(particle_in,particle_out)
    class(particle_base),intent(in) :: particle_in   !< define input particle as particle base
    class(particle_gc),intent(inout) :: particle_out !< define particle output as particle_gc

    particle_out%E =  0.d0 !< initialise gc energy to 0
    particle_out%mu = 0.d0 !< initialise gc magnetic moment to 0
    particle_out%q =  0    !< initialise gc charge to 0

    select type (p_in => particle_in) !< select data type for particle_in
    type is (particle_gc)             !< define the type as gc for particle_in
      particle_out%E = p_in%E   !< copy the gc energy in particle out 
      particle_out%mu = p_in%mu !< copy the gc magnetic moment in particle out
      particle_out%q = p_in%q   !< copy the gc charge in particle out
    end select                        !< end select data type for particle_in
  end subroutine get_particle_gc_attributes !< end of subroutine get_particle_gc_attributes

  !> Set attribute from particle_gc in particle_base datatype.
  !> inputs:
  !>    particle_in:  (particle_gc) particle from which data are copied
  !> outputs:
  !>    particle_out: (particle_base) particle in which data are stored
  !> Author: C. Sommariva, 11/06/2019, email: cristian.sommariva[at]epfl.ch
  pure subroutine set_particle_gc_attributes(particle_in,particle_out)
    class(particle_gc),intent(in)  :: particle_in      !< define particle input as particle_gc
    class(particle_base),intent(inout) :: particle_out !< define output particle as particle base

    select type (p_out => particle_out) !< select data type for particle_out
    type is (particle_gc)               !< define the type as field line for particle_out
      p_out%E = particle_in%E   !< copy the gc energy in particle out
      p_out%mu = particle_in%mu !< copy the gc magnetic moment in particle out
      p_out%q = particle_in%q   !< copy the gc charge in particle out 
    end select                          !< end select data type for particle_out
  end subroutine set_particle_gc_attributes !< end of subroutine set_particle_gc_attributes

  !> Get attribute particle_kinetic from particle_base datatype. This data are copied 
  !> into a particle_kinetic type.
  !> inputs:
  !>   particle_in:   (particle_base) particle from which data are copied
  !> ouputs:
  !>    particle_out: (particle_kinetic) particle in which data are stored
  !> Author: C. Sommariva, 11/06/2019, email: cristian.sommariva[at]epfl.ch
  pure subroutine get_particle_kinetic_attributes(particle_in,particle_out)
    class(particle_base),intent(in) :: particle_in        !< define input particle as particle base
    class(particle_kinetic),intent(inout) :: particle_out !< define particle output as particle_kinetic

    particle_out%v = (/0.d0,0.d0,0.d0/) !< initialise particle out velocity to zero
    particle_out%q = 0                  !< initialise particle out charge to zero

    select type (p_in => particle_in) !< select data type for particle_in
    type is (particle_kinetic)        !< define the type as particle kinetic for particle_in
      particle_out%v = p_in%v !< copy the velocity in particle out 
      particle_out%q = p_in%q !< copy the charge in particle out 
    end select                        !< end select data type for particle_in
  end subroutine get_particle_kinetic_attributes !< end of subroutine get_particle_kinetic_attributes

  !> Set attribute from particle_kinetic in particle_base datatype.
  !> inputs:
  !>    particle_in:  (particle_kinetic) particle from which data are copied
  !> outputs:
  !>    particle_out: (particle_base) particle in which data are stored
  !> Author: C. Sommariva, 11/06/2019, email: cristian.sommariva[at]epfl.ch
  pure subroutine set_particle_kinetic_attributes(particle_in,particle_out)
    class(particle_kinetic),intent(in)  :: particle_in !< define particle input as particle_kinetic
    class(particle_base),intent(inout) :: particle_out !< define output particle as particle base

    select type (p_out => particle_out) !< select data type for particle_out
    type is (particle_kinetic)          !< define the type as particle kinetic for particle_out
      p_out%v = particle_in%v !< copy the velocity in particle out
      p_out%q = particle_in%q !< copy the charge parallel velocity in particle out
    end select                          !< end select data type for particle_out
  end subroutine set_particle_kinetic_attributes !< end of subroutine set_particle_kinetic_attributes

  !> Get attribute particle_kinetic_leapfrog from particle_base datatype. This data are copied 
  !> into a particle_kinetic_leapfrog type.
  !> inputs:
  !>   particle_in:   (particle_base) particle from which data are copied
  !> ouputs:
  !>    particle_out: (particle_kinetic_leapfrog) particle in which data are stored
  !> Author: C. Sommariva, 11/06/2019, email: cristian.sommariva[at]epfl.ch
  pure subroutine get_particle_kinetic_leapfrog_attributes(particle_in,particle_out)
    class(particle_base),intent(in) :: particle_in                 !< define input particle as particle base
    class(particle_kinetic_leapfrog),intent(inout) :: particle_out !< define particle output as
                                                                   !< particle_kinetic_leapfrog
    particle_out%v = (/0.d0,0.d0,0.d0/) !< initialise particle out velocity to zero
    particle_out%q = 0                  !< initialise particle out charge to zero

    select type (p_in => particle_in)   !< select data type for particle_in
    type is (particle_kinetic_leapfrog) !< define the type as particle kinetic leapfrog for particle_in
      particle_out%v = p_in%v !< copy the velocity in particle out 
      particle_out%q = p_in%q !< copy the charge in particle out 
    end select                          !< end select data type for particle_in
  end subroutine get_particle_kinetic_leapfrog_attributes !< end of subroutine get_particle_kinetic_leapfrog_attributes

  !> Set attribute from particle_kinetic_leapfrog in particle_base datatype.
  !> inputs:
  !>    particle_in:  (particle_kinetic_leapfrog) particle from which data are copied
  !> outputs:
  !>    particle_out: (particle_base) particle in which data are stored
  !> Author: C. Sommariva, 11/06/2019, email: cristian.sommariva[at]epfl.ch
  pure subroutine set_particle_kinetic_leapfrog_attributes(particle_in,particle_out)
    class(particle_kinetic_leapfrog),intent(in)  :: particle_in !< define particle input as
                                                                !< particle_kinetic_leapfrog
    class(particle_base),intent(inout) :: particle_out          !< define output particle as particle base
    select type (p_out => particle_out) !< select data type for particle_out
    type is (particle_kinetic_leapfrog) !< define the type as particle kinetic leapfrog for particle_out
      p_out%v = particle_in%v !< copy the velocity in particle out
      p_out%q = particle_in%q !< copy the charge parallel velocity in particle out
    end select                          !< end select data type for particle_out
  end subroutine set_particle_kinetic_leapfrog_attributes !< end of subroutine set_particle_kinetic_leapfrog_attributes

end module mod_particle_types !< end-of-module mod_particle_types

