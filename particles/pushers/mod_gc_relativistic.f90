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
  !>   node_list:       (type_node_list) jorek node list
  !>   element_list:    (type_element_list) jorek element list  
  !>   relativsitci_gc: (particle_gc_relativistic) a particle gc relativistic
  !>   time:            (real8) particle time
  !>   mass:            (real8) paritcle mass
  !>   B:               (real8)(3) magnetic field 
  !> outputs:
  !>   particle_out: (particle_base) the output particle
  subroutine relativistic_gc_to_particle(node_list,element_list,&
       relativistic_gc,particle_out,time,mass,B)
    !> load modules
    use data_structure
    implicit none
    !> declare input variables
    type(type_node_list),intent(in) :: node_list
    type(type_element_list),intent(in) :: element_list
    type(particle_gc_relativistic),intent(in) :: relativistic_gc
    real(kind=8),intent(in) :: time,mass
    real(kind=8),dimension(3),intent(in) :: B
    !> delcare outpur variables
    class(particle_base),intent(out) :: particle_out

    !> select particle type
    select type (particle_out)
    type is (particle_gc)
       particle_out = relativistic_gc_to_gc(relativistic_gc,time,mass,B)
    end select
    
  end subroutine relativistic_gc_to_particle

  !> This procedure transform a relativistic_gs to a particle_gc
  !> inputs:
  !>   relativistic_gc: (particle_gc_relativistic) relativistic gc
  !>   time:            (real8) particle time
  !>   mass:            (real8) particle mass
  !>   B:               (real8)(3) magnetic field
  !> outputs:
  !>   paericle_out:    (particle_gc) gc in energy and magnetic moment
  pure function relativistic_gc_to_gc(relativistic_gc,time,mass,B) &
       result(particle_out)
    !> load modules
    use constants, only: EL_CHG,ATOMIC_MASS_UNIT,SPEED_OF_LIGHT
    implicit none
    !> declare inputs:
    type(particle_gc_relativistic),intent(in) :: relativistic_gc
    real(kind=8),intent(in) :: time,mass
    real(kind=8),dimension(3),intent(in) :: B
    !> dclare output variable
    type(particle_gc) :: particle_out

    !> copy the position
    particle_out%x = relativistic_gc%x
    !> copy the charge
    particle_out%q = relativistic_gc%q
    !> copy mesh element
    particle_out%i_elm = relativistic_gc%i_elm
    !> copy local particle coordinates
    particle_out%st = relativistic_gc%st
    !> copy the magnetic moment in eV with p_parallel sign
    particle_out%mu = sign(ATOMIC_MASS_UNIT*relativistic_gc%p(2)/EL_CHG,&
         relativistic_gc%p(1))
    !> compute the guiding ceneter energy
    particle_out%E = ATOMIC_MASS_UNIT*SPEED_OF_LIGHT*(&
         mass*mass*SPEED_OF_LIGHT*SPEED_OF_LIGHT* + &
         relativistic_gc%p(1)*relativistic_gc%p(2) + &
         2.d0*mass*norm2(B)*relativistic_gc%p(2))/EL_CHG
  end function relativistic_gc_to_gc
  
  !> This procedure transform a particle_gc to relativistic_particle
  !> inputs:
  !>   gc_in: (particle_gc) guiding ceneter in energy momentum
  !>   time:  (real8) particle time
  !>   mass:  (real8) particle mass
  !>   B:     (real8)(3) magnetic field
  !> outputs:
  !>   relativistic_gc: (particle_gc_relativistic) relativistic gc
  pure function gc_to_relativistic_gc(gc_in,time,mass,B) result(relativistic_gc)
    !> load modules
    use constants, only: EL_CHG,ATOMIC_MASS_UNIT,SPEED_OF_LIGHT
    implicit none
    !> declare input variables
    type(particle_gc),intent(in) :: gc_in
    real(kind=8),intent(in) :: time,mass
    real(kind=8),dimension(3),intent(in) :: B
    !> declare output variables
    type(particle_gc_relativistic) :: relativistic_gc
    
    !> copy gc position
    relativistic_gc%x = gc_in%x
    !> copy gc charge
    relativistic_gc%q = gc_in%q
    !> copy particle element
    relativistic_gc%i_elm = gc_in%i_elm
    !> copy local coordinates
    relativistic_gc%st = gc_in%st
    !> initialise magnetic moment
    relativistic_gc%p(2) = abs(EL_CHG*gc_in%mu/ATOMIC_MASS_UNIT)
    !> check if the magnetic field is present
    !> initialise parallel momentum
    relativistic_gc%p(1) = sign(sqrt(((gc_in%E*gc_in%E*EL_CHG*EL_CHG)/&
         (ATOMIC_MASS_UNIT*ATOMIC_MASS_UNIT*SPEED_OF_LIGHT*SPEED_OF_LIGHT))-&
         mass*mass*SPEED_OF_LIGHT*SPEED_OF_LIGHT-&
         2.d0*mass*norm2(B)*relativistic_gc%p(2)),&
         gc_in%mu)
  end function gc_to_relativistic_gc
end module mod_gc_relativistic

