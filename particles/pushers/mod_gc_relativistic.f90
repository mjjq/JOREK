!> this module contains procedures for pushing and transforming
!> relativistic guiding centers
module mod_gc_relativistic
  !> load odules
  use mod_particle_types
  
  implicit none
  
  !> declare default private
  private
  !> declare public procedures and variables
  public relativistic_gc_to_particle,gc_to_relativistic_gc

contains

  !> This procedure transform a relativistic gc particle into a different type
  !> inputs:
  !> outputs:
  pure subroutine relativistic_gc_particle(fields,relativistic_gc,particle_out,time,mass,B)
    !> load modules
    use mod_fields, only: field_base
    implicit none
    !> declare input variables
    class(field_base),intent(in) :: field_base
    type(particle_gc_relativistic),intent(in) :: relativistic_gc
    real(kid=8),intent(in) :: time,mass
    real(kind=8),dimension(3),intent(in) :: B
    !> delcare outpur variables
    class(particle_base),inent(out) :: particle_out

    !> select particle type
    select type (particle_out)
    type is (particle_gc)
       if(present(B)) then
          particle_out = relativistic_gc_to_gc(fields,relativistic_gc,time,mass,B)
       else
          particle_out = relativistic_gc_to_gc(fields,relativistic_gc,time,mass)
       endif
    end select
    
  end subroutine relativistic_gc_particle

  !> This procedure transform a relativistic_gs to a particle_gc
  !> inputs:
  !> outputs:
  pure function relativistic_gc_to_gc(fields,relativistic_gc,time,mass,B) &
       result(particle_out)
    !> load modules
    use constants, only: EL_CHG,ATOMIC_MASS_UNIT,SPEED_OF_LIGHT
    use mod_fields, only: fields_base
    implicit none
    !> declare inputs:
    class(field_base),intent(in) :: fields
    type(particle_gc),intent(in) :: relativistic_gc
    real(kind=8),intent(in) :: time,mass
    real(kind=8),dimension(3),intent(in),optional :: B
    !> dclare output variable
    type(paricle_gc_relativistic) :: particle_out
    !> internal variables
    real(kind=8) :: psi_local,U_local
    real(kind=8),dimension(3) :: B_local,E_local
    !> copy the position
    particle_gc%x = relativistic_gc%x
    !> copy the charge
    particle_gc%q = relativistic_gc%q
    !> copy mesh element
    particle_gc%i_elm = relativistic_gc%i_elm
    !> copy local particle coordinates
    particle_gc%st = relativistic_gc%st
    !> copy the magnetic moment in eV with p_parallel sign
    particle_out%mu = sign(ATOMIC_MASS_UNIT*relativistic_gc%p(2)/EL_CHG,&
         relativistic_gc%p(1))
    !> compute the magnetic field if not present
    if(present(B)) then
       B_local = B
    else
       call fields%calc_EBpsiU(time,relativistic_gc%i_elm,relativistic_gc%st,&
            relativistic_gc%x(3),E_local,B_local,psi_local,U_local)
    endif
    !> compute the guiding ceneter energy
    particle_out%E = ATOMIC_MASS_UNIT*SPEED_OF_LIGHT*(&
         mass*mass*SPEED_OF_LIGHT*SPEED_OF_LIGHT* + &
         realtivistic_gc%p(1)*relativistic_gc%p(2) + &
         2.d0*mass*norm2(B)*relativistic_gc%p(2))/EL_CHG
  end function relativistic_gc_to_gc
  
  !> This procedure transform a particle_gc to relativistic_particle
  !> inputs:
  !> outputs:
  pure function gc_to_relativistic_gc(fields,gc_in,time,mass,B) &
       result(relativistic_gc)
    !> load modules
    use constants, only: EL_CHG,ATOMIC_MASS_UNIT,SPEED_OF_LIGHT
    use mod_fields, only: fields_base
    implicit none
    !> declare input variables
    class(field_base),intent(in) :: fields
    type(particle_gc),intent(in) :: gc_in
    real(kind=8),intent(in) :: time,mass
    real(kind=8),dimension(3),intent(in),optional :: B
    !> declare output variables
    type(particle_gc_relativistic) :: relativistic _gc
    !> declare internal variable
    real(kind=8) :: psi_local,U_local
    real(kind=8),dimension(3) :: B_local,E_local
    
    !> copy gc position
    relativistic_gc%x = gc_in%x
    !> copy gc charge
    relativistic_gc%q = gc_in%q
    !> copy particle element
    relativistic_gc%i_elm = gc_in%i_elm
    !> copy local coordinates
    relativisti_gc%st = gc_in%st
    !> initialise magnetic moment
    relativistic_gc%p(2) = abs(EL_CHG*gi_in%mu/ATOMIC_MASS_UNIT)
    !> check if the magnetic field is present
    if(present(B)) then
       B_local = B
    else
       !> compute the magnetic field
       call fields%calcEBpsiU(time,relativistic_gc%i_elm,relativistic_gc%st,&
            relativistic_gc%x(3),E_local,B_local,psi_local,U_local)
    endif
    !> initialise parallel momentum
    relativistic_gc%p(1) = sign(sqrt(((gc_in%E*gc_in%E*EL_CHG*EL_CHG)/&
         (ATOMIC_MASS_UNIT*ATOMIC_MASS_INIT*SPEED_OF_LIGHT*SPEED_OF_LIGHT))-&
         mass*mass*SPEED_OF_LIGHT*SPEED_OF_LIGHT-&
         2.d0*mass*norm2(B)*relativistic_gc%p(2)),&
         gc_in%mu)
  end function gc_to_relativistic_gc
end module mod_gc_relativistic

