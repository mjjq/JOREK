!> this module contains procedures for pushing and transforming
!> relativistic guiding centers
module mod_gc_relativistic
  !> load modules
  use mod_particle_types
  
  implicit none
  
  !> declare default private
  private
  !> declare public procedures and variables
  public relativistic_gc_to_particle, gc_to_relativistic_gc

contains

  !> This procedure computes the guiding ceneter derivatives required
  !> for the runge_kutta integration.
  !> inputs:
  !>   fields:
  !>   n_variables:       (integer) number of variables=4:
  !>   n_int_parameters:  (integer) number of integer parameters=2
  !>   n_real_parameters: (integer) number of real parameters
  !>   t:
  !>   solution_old:     (real8)(n_variables) initial solution:
  !>                     1:R, 2:Z, 3:phi, 4:parallel momentum
  !>   solution:         (real8)(n_variables) solution at a runge-kutta stage
  !>   int_parameters:   (integer)(n_int_parameters) integer parameters
  !>                     1:old i_elm, 2: charge
  !>   real_parameters:  (real8)(n_real_parameters) real parameters:
  !>                     1:mass 2:magnetic moment
  !>   derivatives:      (real8)(n_variables) runge-kutta derivatives
  !>   ifail:            (integer) if 0 calculation failed
  subroutine compute_derivatives_runge_kutta(fields,n_variables,&
       n_int_parameters,n_real_parameters,t,solution_old,solution,&
       int_parameters,real_parameters,derivatives,ifail)
    !> declare input variables
    
    
  end subroutine compute_derivatives_runge_kutta
  
  !> This procedure transforms a relativistic gc particle into a different type
  !> inputs:
  !>   node_list:       (type_node_list) jorek node list
  !>   element_list:    (type_element_list) jorek element list  
  !>   relativistic_gc: (particle_gc_relativistic) a particle gc relativistic
  !>   time:            (real8) particle time
  !>   mass:            (real8) particle mass
  !>   B:               (real8)(3) magnetic field
  !>   gyro_angle:      (real8)(optional) the gyro-angle, defaut=0
  !> outputs:
  !>   particle_out: (particle_base) the output particle
  !>   gyro_angle:   (real8)(optiona) the gyro-angle, default=0
  subroutine relativistic_gc_to_particle(node_list,element_list,&
       relativistic_gc,particle_out,mass,B,gyro_angle)
    !> load modules
    use data_structure
    implicit none
    !> declare input variables
    type(type_node_list), intent(in)           :: node_list
    type(type_element_list), intent(in)        :: element_list
    type(particle_gc_relativistic), intent(in) :: relativistic_gc
    real(kind=8), intent(in)                   :: mass
    real(kind=8), dimension(3), intent(in)     :: B
    !> declare output variables
    class(particle_base), intent(out)          :: particle_out
    !> declare input output variables
    real(kind=8), intent(inout), optional      :: gyro_angle

    !< if the gyroangle is not present set it to zero
    if(.not.present(gyro_angle)) gyro_angle = 0.d0
    
    !> select particle type
    select type (particle_out)
    type is (particle_gc)
       particle_out = relativistic_gc_to_gc(relativistic_gc,mass,B)
    type is (particle_kinetic_relativistic)
       particle_out = relativistic_gc_to_relativistic_kinetic(&
            node_list,element_list,relativistic_gc,mass,B,gyro_angle)
    end select
    
  end subroutine relativistic_gc_to_particle

  !> This procedure transforms a particle_gc_relativistic 
  !> into a particle_kinetic_relativistic
  !> inputs:
  !>   node_list:       (type_node_list) jorek node list
  !>   element_list:    (type_element_list) jorek element list
  !>   relativistic_gc: (particle_gc_relativistic) a relativistic gc
  !>   mass:            (real8) particle mass 
  !>   B:               (real8)(3) magnetic field
  !>   gyro_angle:      (real8) gyro-angle for initialising a particle
  !> outputs:
  !>   relativistic_particle: (particle_kinetic_relativistic) a
  !>                           relativistic kinetic particle
  function relativistic_gc_to_relativistic_kinetic(node_list,element_list, &
       in,mass,B,gyro_angle) result(out)
    !> load modules
    use data_structure
    use mod_coordinate_transforms, only: vector_cylindrical_to_cartesian
    use mod_pusher_tools, only: gc_position_to_particle,get_orthonormals
    implicit none
    !> declare input varibales
    type(type_node_list), intent(in)           :: node_list
    type(type_element_list), intent(in)        :: element_list
    type(particle_gc_relativistic), intent(in) :: in
    real(kind=8), intent(in)                   :: mass, gyro_angle
    real(kind=8), dimension(3), intent(in)     :: B
    !> declare output variables
    type(particle_kinetic_relativistic)        :: out
    !> declare internal varibales
    real(kind=8)                               :: B_norm, p_perp
    real(kind=8), dimension(3)                 :: B_hat, e1, e2
    
    out = in !< copy the default particle type
    out%q = in%q !< copy the charge
    
    !> compute the magnetic field intensity and direction
    B_norm = norm2(B)
    B_hat = B/B_norm
    !> construct orthogonal basis
    call get_orthonormals(B_hat,e1,e2)
    !> compute the perpendicular momentum
    p_perp = sqrt(2.d0*mass*B_norm*in%p(2))
    !> compute particle momenta
    out%p = in%p(1)*B_hat + p_perp*(e1*cos(gyro_angle)+e2*sin(gyro_angle))
    !> check if the particle is valid
    if(out%q.ne.0) then
       call gc_position_to_particle(node_list,element_list,in%x,&
            in%st,in%i_elm,out%p,out%q,B_hat,B_norm,out%x,&
            out%st,out%i_elm)
    endif
    !> transform the momenta to cartesian coordinates
    out%p = vector_cylindrical_to_cartesian(out%x(3),out%p)
  end function relativistic_gc_to_relativistic_kinetic 

  !> This procedure transforms a particle_gc_relativistic
  !> into a particle_gc
  !> inputs:
  !>   relativistic_gc: (particle_gc_relativistic) relativistic gc
  !>   mass:            (real8) particle mass
  !>   B:               (real8)(3) magnetic field
  !> outputs:
  !>   particle_out:    (particle_gc) gc in energy and magnetic moment
  pure function relativistic_gc_to_gc(relativistic_gc,mass,B) &
       result(particle_out)
    !> load modules
    use constants, only: EL_CHG, ATOMIC_MASS_UNIT, SPEED_OF_LIGHT
    implicit none
    !> declare inputs:
    type(particle_gc_relativistic), intent(in) :: relativistic_gc
    real(kind=8), intent(in)                   :: mass
    real(kind=8), dimension(3), intent(in)     :: B
    !> declare output variable
    type(particle_gc)                          :: particle_out

    !> copy the position
    particle_out%x = relativistic_gc%x
    !> copy the charge
    particle_out%q = relativistic_gc%q
    !> copy mesh element
    particle_out%i_elm = relativistic_gc%i_elm
    !> copy local particle coordinates
    particle_out%st = relativistic_gc%st
    !> copy the magnetic moment in eV/T with p_parallel sign
    particle_out%mu = sign(ATOMIC_MASS_UNIT*relativistic_gc%p(2)/EL_CHG, &
         relativistic_gc%p(1))
    !> compute the guiding center energy in eV
    particle_out%E = ATOMIC_MASS_UNIT*SPEED_OF_LIGHT          &
         * sqrt(mass*mass*SPEED_OF_LIGHT*SPEED_OF_LIGHT +     &
                relativistic_gc%p(1)*relativistic_gc%p(1) +   &
                2.d0*mass*norm2(B)*relativistic_gc%p(2)) /EL_CHG
  end function relativistic_gc_to_gc
  
  !> This procedure transforms a particle_gc into particle_gc_relativistic
  !> inputs:
  !>   gc_in: (particle_gc) guiding center in energy momentum
  !>   mass:  (real8) particle mass
  !>   B:     (real8)(3) magnetic field
  !> outputs:
  !>   relativistic_gc: (particle_gc_relativistic) relativistic gc
  pure function gc_to_relativistic_gc(gc_in,mass,B) result(relativistic_gc)
    !> load modules
    use constants, only: EL_CHG, ATOMIC_MASS_UNIT, SPEED_OF_LIGHT
    implicit none
    !> declare input variables
    type(particle_gc), intent(in)          :: gc_in
    real(kind=8), intent(in)               :: mass
    real(kind=8), dimension(3), intent(in) :: B
    !> declare output variables
    type(particle_gc_relativistic)         :: relativistic_gc
    
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
    !> initialise parallel momentum
    relativistic_gc%p(1) = sign(sqrt(((gc_in%E*gc_in%E*EL_CHG*EL_CHG)/         &
         (ATOMIC_MASS_UNIT*ATOMIC_MASS_UNIT*SPEED_OF_LIGHT*SPEED_OF_LIGHT)) -  &
         mass*mass*SPEED_OF_LIGHT*SPEED_OF_LIGHT -                             &
         2.d0*mass*norm2(B)*relativistic_gc%p(2)),gc_in%mu)
  end function gc_to_relativistic_gc
end module mod_gc_relativistic

