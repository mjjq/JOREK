!> Field line tracing using Runge-Kutta scheme
module mod_fieldline_runge_kutta

  use mod_particle_types
  
  implicit none
  
  private
  public runge_kutta_fixed_dt_fieldline_push_jorek

contains

  !> Computes the RHS of the field line tracing equations
  !!
  !! This routine is called by the Runge-Kutta integrator and hence is
  !! declared private.
  subroutine compute_fieldline_rhs(fields,n_variables, &
       n_int_parameters,n_real_parameters,t,solution_old,solution,            &
       int_parameters,real_parameters,derivatives,ifail)
  
    !> load modules
    use mod_fields, only: fields_base
    use mod_find_rz_nearby
    implicit none
  
    class(fields_base), intent(in)                         :: fields
    integer, intent(in)                                    :: n_variables, n_int_parameters, n_real_parameters
    integer, dimension(n_variables), intent(in)            :: int_parameters
    real(kind=8), intent(in)                               :: t
    real(kind=8), dimension(n_variables), intent(in)       :: solution, solution_old
    real(kind=8), dimension(n_real_parameters), intent(in) :: real_parameters
  
    integer, intent(out)                              :: ifail
    real(kind=8), dimension(n_variables), intent(out) :: derivatives

    real(kind=8) :: E(3), B(3), psi, U, Bnorm, st_new(2)
    integer :: ierr

    ! Evaluate magnetic field vector at the marker position
    call find_RZ_nearby(fields%node_list,fields%element_list,solution_old(1), &
         solution_old(2),real_parameters(1),real_parameters(2),               &
         int_parameters(1),solution(1),solution(2),st_new(1),                 &
         st_new(2),ifail,ierr)

    if(ifail .ne. 0) then
       call fields%calc_EBpsiU(t, int_parameters(1), st_new, solution(3), E, B, psi, U)
    end if

    ! Marker is pushed along the field line vector. Note 1/R in the toroidal component.
    Bnorm = norm2(B)
    derivatives(1) = real_parameters(3) * B(1) / Bnorm
    derivatives(2) = real_parameters(3) * B(2) / Bnorm
    derivatives(3) = real_parameters(3) * B(3) / ( Bnorm * solution(1) )

  end subroutine compute_fieldline_rhs


  !> Solves marker trajectory along a field line using RK4
  !!
  !! This is the public routine in this module which is used to integrate field line
  !! trajectory for a single step. The distance (in meters) marker is pushed is particle%v * dt 
  !! so choose those parameters accordingly (particle%v is just a parameter that is not evolved).
  subroutine fieldline_runge_kutta_fixed_dt_push_jorek(fields, particle, t, dt)
    !> modules                                                                                                                                                                                                                                                                           
    use mod_fields, only: fields_base
    use mod_find_rz_nearby
    use mod_runge_kutta, only: runge_kutta_fixed_dt
    implicit none
    !> input/output variables                                                                                                                                                                                                                                                            
    type(particle_fieldline), intent(inout) :: particle !< Marker that is advanced
    real(kind=8), intent(in)                :: t        !< Time at which the field is evaluated
    class(fields_base), intent(in)          :: fields   !< The field data
    real(kind=8), intent(in)                :: dt       !< Step size (in seconds if particle%v is in m/s)
    !> internal variables
    integer                    :: ifail, i_elm_new 
    real(kind=8), dimension(2) :: st_new
    
    !> Solved marker (R,Z,Phi) position
    real(kind=8), dimension(3) :: solution_new

    !> Compute RK4 step
    call runge_kutta_fixed_dt(compute_fieldline_rhs,                    &
         fields,3,1,2,t,dt,[particle%x(1),particle%x(2),particle%x(3)], &
         [particle%i_elm],[particle%st(1),particle%st(2),particle%v],   &
         solution_new, i_elm_new)

    !> Compute the new local coordinates
    if(i_elm_new.ne.0) call find_rz_nearby(fields%node_list,  &
         fields%element_list,particle%x(1),particle%x(2),     &
         particle%st(1),particle%st(2),particle%i_elm,        &
         solution_new(1),solution_new(2),st_new(1),st_new(2), &
         i_elm_new,ifail)

    !> Update fields
    particle%x     = solution_new(1:3)
    particle%st    = st_new
    particle%i_elm = i_elm_new
  
  end subroutine fieldline_runge_kutta_fixed_dt_push_jorek

end module mod_fieldline_runge_kutta
