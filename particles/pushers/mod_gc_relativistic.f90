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
  public relativistic_gc_to_relativistic_kinetic
  public relativistic_gc_momenta_from_E_cospitch
  public runge_kutta_fixed_dt_gc_push_jorek
  public runge_kutta_fixed_dt_gc_push

contains

  !> This procedure push a relativistic guiding center particle in JOREK
  !> fields using standard runge kutta integrator withou time step control
  !> inputs:
  !>   fields:   (fields_base) jorek fields
  !>   t:        (real8) simulation time
  !>   dt:       (real8) simulation time step
  !>   mass:     (real8) the gc mass in AMU
  !>   particle: (particle_gc_relativistic) the gc to integrate
  !> outputs:
  !>   particle: (particle_gc_relativistic) the integrated gc
  subroutine runge_kutta_fixed_dt_gc_push_jorek(fields,t,dt,&
       mass,particle)
    !> load modules
    use mod_fields, only: fields_base
    use mod_find_rz_nearby
    use mod_runge_kutta, only: runge_kutta_fixed_dt
    implicit none
    !> declare input output variables
    type(particle_gc_relativistic),intent(inout) :: particle
    !> declare input variables
    class(fields_base),intent(in) :: fields
    real(kind=8),intent(in) :: t,dt,mass
    !> internal variables
    integer :: ifail,i_elm_new !< particle new element
    !> new particle local coordinates: 1:s,2:t
    real(kind=8),dimension(2) :: st_new
    !> integrate global coordinates: 1:R, 2:Z, 3:phi, 4:p_parallel
    real(kind=8),dimension(4) :: solution_new


    !> compute runge kutta differentials
    call runge_kutta_fixed_dt(compute_relativistic_gc_derivatives_jorek,&
         fields,4,2,4,t,dt,[particle%x(1),particle%x(2),particle%x(3),&
         particle%p(1)],[particle%i_elm,int(particle%q)],[particle%st(1),&
         particle%st(2),mass,particle%p(2)],solution_new,i_elm_new)
    
    !> compute the new local coordinates
    if(i_elm_new.ne.0) call find_rz_nearby(fields%node_list,&
         fields%element_list,particle%x(1),particle%x(2),&
         particle%st(1),particle%st(2),particle%i_elm,&
         solution_new(1),solution_new(2),st_new(1),st_new(2),&
         i_elm_new,ifail)
    
    !> overwrite particle fields
    particle%x = solution_new(1:3)
    particle%p(1) = solution_new(4)
    particle%st = st_new
    particle%i_elm = i_elm_new
    
  end subroutine runge_kutta_fixed_dt_gc_push_jorek

  !> This subroutine integrates a relativistic gc using fixed time step runge
  !> kutta method on analytical fields.
  !> inputs:
  !>   fields:   (fields_base) jorke fields
  !>   t:        (real8) integration time
  !>   dt:       (real8) time step
  !>   mass:     (real8) particle am
  !>   particle: (particle_gc_relativistic) particle to push
  !> outputs:
  !>   particle: (particle_gc_relativistic) pushed particle
  subroutine runge_kutta_fixed_dt_gc_push(fields,t,dt,mass,particle)
    !> load modules
    use mod_fields, only: fields_base
    use mod_runge_kutta, only: runge_kutta_fixed_dt
    implicit none
    !> declare input output varibales
    type(particle_gc_relativistic),intent(inout) :: particle
    !> delcare input variables
    class(fields_base),intent(in) :: fields
    real(kind=8),intent(in) :: t,dt,mass
    !> delcare internal variables
    integer :: ifail
    real(kind=8),dimension(4) :: solution_new

    !> integrate particle trajectory
    call runge_kutta_fixed_dt(compute_relativistic_gc_derivatives,&
         fields,4,1,2,t,dt,[particle%x(1),particle%x(2),particle%x(3),&
         particle%p(1)],[int(particle%q)],[mass,particle%p(2)],&
         solution_new,ifail)
    !> overwrite the new particle position
    particle%x = solution_new(1:3)
    particle%p(1) = solution_new(4)
    
  end subroutine runge_kutta_fixed_dt_gc_push
  
  !> This procedure computes the guiding ceneter derivatives required
  !> for the runge_kutta integration.
  !> inputs:
  !>   fields:
  !>   n_variables:       (integer) number of variables=4:
  !>   n_int_parameters:  (integer) number of integer parameters=2
  !>   n_real_parameters: (integer) number of real parameters=4
  !>   t:
  !>   solution_old:     (real8)(n_variables) initial solution:
  !>                     1:R, 2:Z, 3:phi, 4:parallel momentum
  !>   solution:         (real8)(n_variables) solution at a runge-kutta stage
  !>   int_parameters:   (integer)(n_int_parameters) integer parameters
  !>                     1:old i_elm, 2: charge
  !>   real_parameters:  (real8)(n_real_parameters) real parameters:
  !>                     1:s_old, 2:t_old, 3:mass 4:magnetic moment
  !> outputs:
  !>   derivatives:      (real8)(n_variables) runge-kutta derivatives
  !>   ifail:            (integer) if 0 calculation failed
  subroutine compute_relativistic_gc_derivatives_jorek(fields,n_variables,&
       n_int_parameters,n_real_parameters,t,solution_old,solution,&
       int_parameters,real_parameters,derivatives,ifail)
    !> load modules
    use constants, only: SPEED_OF_LIGHT,EL_CHG,ATOMIC_MASS_UNIT
    use mod_fields, only: fields_base
    use mod_math_operators, only: cross_product
    use mod_find_rz_nearby
    implicit none
    !> declare input variables
    class(fields_base),intent(in) :: fields
    integer,intent(in) :: n_variables,n_int_parameters,n_real_parameters
    integer,dimension(n_variables),intent(in) :: int_parameters
    real(kind=8),intent(in) :: t
    real(kind=8),dimension(n_variables),intent(in) :: solution,solution_old
    real(kind=8),dimension(n_real_parameters),intent(in) :: real_parameters
    !> declare output variables
    integer,intent(out) :: ifail
    real(kind=8),dimension(n_variables),intent(out) :: derivatives
    !> internal variables
    integer :: ierr !< error for find_rz nearby
    real(kind=8) :: normB,gamma !< magnetic field intensity and relativistic factor
    real(kind=8),dimension(2) :: st_new !< local postion particle at stage
    !> define guiding center fields
    real(kind=8),dimension(3) :: E,b,gradB,curlb,dbdt,B_star

    !> find the gc at stage local position
    call find_RZ_nearby(fields%node_list,fields%element_list,solution_old(1),&
         solution_old(2),real_parameters(1),real_parameters(2),&
         int_parameters(1),solution(1),solution(2),st_new(1),&
         st_new(2),ifail,ierr)

    !> compute the guiding center fields
    if(ifail.ne.0) call fields%calc_EBNormBGradBCurlbDbdt(t,ifail,st_new,&
         solution(3),E,b,normB,gradB,curlb,dbdt)

    !> compute the guiding center rhs
    derivatives = compute_relativistic_gc_rhs(int_parameters(2),&
         real_parameters(3),real_parameters(4),solution(1),&
         solution(4),normB,E,b,gradB,curlb,dbdt)
  end subroutine compute_relativistic_gc_derivatives_jorek

  !> This procedure compute the guiding center derivatives in analytical
  !> fields. This is mainly used for testing models.
  !> inputs:
  !>   fields:            (fields_base) jorek fields
  !>   n_variables:       (n_variables) number of variables
  !>   n_int_parameters:  (integer) number of integer parameters = 1
  !>   n_real_parameters: (integer) number of real parameters = 2
  !>   t:                 (real8) time at stage
  !>   solution_old:      (real8)(n_variables) old particle state
  !>   solution:          (real8)(n_variables) new particle state at stage
  !>   int_parameters:    (integer)(n_int_parameters) 1: charge
  !>   real_parameters:   (real8)(n_real_parameters) 1:mass, 2:magnetic moment
  !> outputs:
  !>   ifail:       (integer) if 0 the integration failed
  !>   derivatives: (real8)(n_variables) guiding center right field side
  pure subroutine compute_relativistic_gc_derivatives(fields,n_variables,&
       n_int_parameters,n_real_parameters,t,solution_old,solution,&
       int_parameters,real_parameters,derivatives,ifail)
    !> load modules
    use mod_fields, only: fields_base
    implicit none
    !> declare input varibales
    class(fields_base),intent(in) :: fields
    integer,intent(in) :: n_variables,n_int_parameters,n_real_parameters
    real(kind=8),intent(in) :: t
    real(kind=8),dimension(n_variables),intent(in) :: solution_old,solution
    integer,dimension(n_int_parameters),intent(in) :: int_parameters
    real(kind=8),dimension(n_real_parameters),intent(in) :: real_parameters
    !> declare output variables
    integer,intent(out) :: ifail
    real(kind=8),dimension(n_variables),intent(out) :: derivatives
    !> internal variables
    real(kind=8) :: normB
    real(kind=8),dimension(3) :: E,b,gradB,curlb,dbdt

    !> compute the new electromagnetic fields
    call fields%calc_analytical_EBNormBGradBCurlbDbdt(solution(1:2),E,b,normB,&
         gradB,curlb,dbdt)
    !> compute gc right hand side
    derivatives = compute_relativistic_gc_rhs(int_parameters(1),real_parameters(1),&
         real_parameters(2),solution(1),solution(4),normB,E,b,gradB,curlb,dbdt)
    !> set ifail to true
    ifail = 1
    
  end subroutine compute_relativistic_gc_derivatives
  
  !> This procedure computes the guiding center equaction right hand side
  !> as reported in J.R Cary, A.J. Brizard, Rev. Mod. Phys, vol.81, p.693 ,2009
  !> inputs:
  !> outputs:
  !>   derivatives: (real8)(4) gc right hand side: 1:R_dot,2:Z_dot,
  !>   3:phi_dot,4:p_parallel_dot
  pure function compute_relativistic_gc_rhs(charge,mass,magnetic_moment,&
       R,p_parallel,normB,E,b,gradB,curlb,dbdt) result(derivatives)
    use constants, only: EL_CHG,SPEED_OF_LIGHT,ATOMIC_MASS_UNIT
    use mod_math_operators, only: cross_product
    implicit none
    !> declare input variables
    integer,intent(in) :: charge
    real(kind=8),intent(in) :: mass,R,p_parallel,magnetic_moment,normB
    real(kind=8),dimension(3),intent(in) ::b,E,gradB,curlb,dbdt
    !> declare output variable
    real(kind=8),dimension(4) :: derivatives
    !> declare input variables
    real(kind=8) :: gamma !< relativistic factor
    real(kind=8),dimension(3) :: B_star,E_star

    !> compute the relativistic factor
    gamma = sqrt(mass*mass*SPEED_OF_LIGHT*SPEED_OF_LIGHT + &
         p_parallel*p_parallel+2.d0*mass*normB*magnetic_moment)/&
         (mass*SPEED_OF_LIGHT)

    !> compute B_star and E_star
    B_star = p_parallel*curlb + &
         ((EL_CHG*real(charge,kind=8)*normB)/ATOMIC_MASS_UNIT)*b !< B_star
    E_star = (EL_CHG*real(charge,kind=8)*E/ATOMIC_MASS_UNIT) - &
         p_parallel*dbdt - ((magnetic_moment*gradB)/gamma) !< E_star
    
    !> compute the guiding center position derivatives
    derivatives(1:3) = (cross_product(E_star,b) +&
         ((p_parallel*B_star)/(mass*gamma)))
    derivatives(3) = derivatives(3)/R
    derivatives(4) = B_star(1)*E_star(1)+B_star(2)*E_star(2)+B_star(3)*E_star(3)
    derivatives = derivatives/(B_star(1)*b(1)+B_star(2)*b(2)+B_star(3)*b(3))
    
  end function compute_relativistic_gc_rhs
  
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
    !> declare internal variables
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

  !> This function fills in the p(1) and p(2) fields of a particle_gc_relativistic
  !> particle from its energy and (cosine of) pitch-angle
  function relativistic_gc_momenta_from_E_cospitch(rel_gc_in,energy,ksi,mass,fields,time) result(rel_gc_out)

    use constants, only: ATOMIC_MASS_UNIT, EL_CHG, SPEED_OF_LIGHT
    use mod_fields

    implicit none

    !> input variables
    type(particle_gc_relativistic), intent(in) :: rel_gc_in
    real*8, intent(in)                         :: energy !< Particle energy in eV
    real*8, intent(in)                         :: ksi    !< Cosine of particle pitch-angle 
    real*8, intent(in)                         :: mass   !< Particle mass in AMU
    class(fields_base), intent(in)             :: fields
    real*8, intent(in)                         :: time

    !> output variables
    type(particle_gc_relativistic) :: rel_gc_out

    !> internal variables
    real*8               :: p_norm_sq, energy_at_rest
    real*8               :: psi, U, B_norm
    real*8, dimension(3) :: E, B

    !> Copy existing fields from input to output
    rel_gc_out%q     = rel_gc_in%q
    rel_gc_out%x     = rel_gc_in%x
    rel_gc_out%i_elm = rel_gc_in%i_elm
    rel_gc_out%st    = rel_gc_in%st
    
    if (abs(ksi)>1) then
      write(*,*) 'Error in relativistic_gc_momenta_from_E_cospitch: abs(cos(pitch-angle)) should be <1, exiting'
      stop
    end if
    
    energy_at_rest = ATOMIC_MASS_UNIT*mass*SPEED_OF_LIGHT**2
    
    p_norm_sq = ((energy*EL_CHG)**2-energy_at_rest**2)/SPEED_OF_LIGHT**2

    !> mu [AMU (m/s)**2 /T]
    call fields%calc_EBpsiU(time, rel_gc_out%i_elm, rel_gc_out%st, rel_gc_out%x(3), E, B, psi, U)
    B_norm = norm2(B)
    rel_gc_out%p(2) = p_norm_sq*(1.d0-ksi**2)/(2.d0*B_norm*mass*ATOMIC_MASS_UNIT**2)

    !> Parallel momentum [AMU m/s]
    rel_gc_out%p(1) = sqrt(p_norm_sq)*ksi/ATOMIC_MASS_UNIT		       

  end function relativistic_gc_momenta_from_E_cospitch

end module mod_gc_relativistic

