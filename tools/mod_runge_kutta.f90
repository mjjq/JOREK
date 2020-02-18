!> This module contains all procedures and variables required
!> for solving a seto of ordinary differential equatons using
!>a Runge-Kutta scheme. By default, the Runge-Kutta 4(5)
!> Cash-Karp is implemented.
module mod_runge_kutta

  private !< set all as private
  public runge_kutta_fixed_dt,runge_kutta_order_error_control_dt
  public runge_kutta_adaptative_dt

  !> declare module parameters
  integer,parameter :: n_stages=6 !< number of stages

  !> interfaces

  !> this is is a common interface for all procedures computing
  !> solutions from runge-kutta differentials
  interface compute_runge_kutta_solution
     module procedure compute_runge_kutta_solution_1,&
          compute_runge_kutta_solution_2,compute_runge_kutta_solution_5
  end interface compute_runge_kutta_solution
  
  abstract interface
     !> This is the interface for procedures computing the ODE
     !> right hanf side to be used in compute_runge_kutta_derivatives
     !> inputs:
     !>   n_variables:       (integer) number of variables
     !>   n_int_parameters:  (integer) number of integer parameters
     !>   n_real_parameters: (integer) number of real parameters
     !>   t:                 (real8) integration coordinate
     !>   solution_old:      (real8)(n_variables) old solution
     !>   solution:          (real8)(n_variables) solution at a give stage
     !>   int_parameters:    (integer)(n_int_parameters) integer parameters
     !>   real_parameters:   (real8)(n_real_parameters) real parameters
     !>   fields:            (field_base)(optional) fields for particle pushing
     !> outputs
     !>   derivatives:  (real8)(n_variables) new derivatives
     !>   ifail:        (integer) if 0 derivative calculation failed
     subroutine compute_runge_kutta_rhs(fields,n_variables,&
          n_int_parameters,n_real_parameters,t,solution_old,solution,&
          int_parameters,real_parameters,derivatives,ifail)
       !> load module
       use mod_fields, only: fields_base
       implicit none
       !> declare input variables
       integer,intent(in) :: n_variables,n_int_parameters,n_real_parameters
       integer,dimension(n_int_parameters),intent(in) :: int_parameters
       real(kind=8),intent(in) :: t
       real(kind=8),dimension(n_variables),intent(in) :: solution_old,solution
       real(kind=8),dimension(n_real_parameters),intent(in) :: real_parameters
       class(fields_base),intent(in) :: fields
       !> declare output variables
       integer,intent(out) :: ifail
       real(kind=8),dimension(n_variables),intent(out) :: derivatives
     end subroutine compute_runge_kutta_rhs

     !> interface of the function computing user defined errors
     !> inputs:
     !>   fields:            (field_base)(optional) fields for particle pushing
     !>   n_variables:       (integer) number of variables
     !>   n_int_parameters:  (integer) number of integer parameters
     !>   n_real_parameters: (integer) number of real parameters
     !>   t:                 (real8) integration coordinate
     !>   solution_1:        (real8)(n_variables) old / low-order solution
     !>   solution_2:        (real8)(n_variables) new / high-order solution
     !>   int_parameters:    (integer)(n_int_parameters) integer parameters
     !>   real_parameters:   (real8)(n_real_parameters) real parameters
     !> outputs:
     !>   maximum_norm_error: (real8) maximum normalised error error/tolerance
     function compute_user_error(fields,n_variables,n_int_parameters,&
          n_real_parameters,t,solution_1,solution_2,int_parameters,&
          real_parameters) result(maximum_norm_error)
       !> load modules
       use mod_fields, only: fields_base
       implicit none
       !> declare inputs
       class(fields_base),intent(in) :: fields
       integer,intent(in) :: n_variables,n_int_parameters,n_real_parameters
       real(kind=8),intent(in) :: t
       real(kind=8),dimension(n_variables),intent(in) :: solution_1,solution_2
       integer,dimension(n_int_parameters),intent(in) :: int_parameters
       real(kind=8),dimension(n_real_parameters),intent(in) :: real_parameters
       !> declare outputs
       real(kind=8) :: maximum_norm_error
     end function compute_user_error

     !> interface of the function adapting time steps using
     !> user defined rules.
     !> inputs:
     !>   fields:            (field_base)(optional) fields for particle pushing
     !>   n_variables:       (integer) number of variables
     !>   n_int_parameters:  (integer) number of integer parameters
     !>   n_real_parameters: (integer) number of real parameters
     !>   t:                 (real8) integration coordinate
     !>   dt:                (real8) old time step
     !>   solution:        (real8)(n_variables) solution
     !>   int_parameters:    (integer)(n_int_parameters) integer parameters
     !>   real_parameters:   (real8)(n_real_parameters) real parameters
     !> outputs:
     !>   dt_new: (real8) new time step
     function adapt_time_step(fields,n_variables,n_int_parameters,&
          n_real_parameters,t,dt,solution,int_parameters,&
          real_parameters) result(dt_new)
       !> load modules
       use mod_fields, only: fields_base
       implicit none
       !> declare inputs
       class(fields_base),intent(in) :: fields
       integer,intent(in) :: n_variables,n_int_parameters,n_real_parameters
       real(kind=8),intent(in) :: t,dt
       real(kind=8),dimension(n_variables),intent(in) :: solution
       integer,dimension(n_int_parameters),intent(in) :: int_parameters
       real(kind=8),dimension(n_real_parameters),intent(in) :: real_parameters
       !> declare outputs
       real(kind=8) :: dt_new
     end function adapt_time_step   
  end interface
  
contains

  !> this procedure implement a runge-kutta integrator with
  !> time step control with error estimated from multiple
  !> orders.This is basically a feedback integration step
  !> cotroller.
  !> inputs:
  !>   compute_rhs:       (procedure) subroutine for computing
  !>                      the ODE(s) right hand side
  !>   fields:            (fields_base) jorek fields structure
  !>   n_variables:       (integer) number of varibales
  !>   n_int_parameters:  (integer) number of integer parameters
  !>   n_real_parameters: (integer) number of real parameters
  !>   t:                 (real8) integration variables
  !>   t_stop:            (real8) ccomputation stop time
  !>   dt:                (real8) integration step
  !>   solution_old:      (real8)(n_variables) old solution
  !>   int_parameters:    (integer)(n_integer_parameters)
  !>                      integer parameters
  !>   real_parameters:   (real8)(n_real_parameters)
  !>                      real parameters
  !>   tolerances:        (real8)(n_variables) base error tolerances
  !>   compute_user_err:  (real8)(n_errors)(optional) external function
  !>                       allowing the user to compute an error based
  !>                       on case specific error
  !> outputs:
  !>   dt:       real(8) new time step
  !>   solution: (n_variables) runge kutta step solution
  !>   ifail:    (integer) if 0 the integration failed
  subroutine runge_kutta_order_error_control_dt(compute_rhs,fields,n_variables,&
       n_int_parameters,n_real_parameters,t,t_stop,dt,solution_old,&
       int_parameters,real_parameters,tolerances,solution,ifail,compute_user_err)
    !> load modules
    use mod_fields, only: fields_base
    implicit none
    !> declare parameters
    integer,parameter :: maximum_iteration=100
    !> shampine error estimate parameters:
    !>   1: safety factor
    !>   2: time step reduction
    !>   3: time step increment
    real(kind=8),dimension(3),parameter :: error_parameters=[&
         9.d-1,1.d0/4.d0,1.d0/5.d0]
    !> delcare input / outputs
    real(kind=8),intent(inout) :: dt
    !> delcare inputs
    procedure(compute_runge_kutta_rhs) :: compute_rhs
    procedure(compute_user_error),optional :: compute_user_err
    class(fields_base),intent(in) :: fields
    integer,intent(in) :: n_variables,n_int_parameters,n_real_parameters
    real(kind=8),intent(in) :: t,t_stop
    real(kind=8),dimension(n_variables),intent(in) :: solution_old,tolerances
    integer,dimension(n_int_parameters),intent(in) :: int_parameters
    real(kind=8),dimension(n_real_parameters),intent(in) :: real_parameters
    !> declare outputs
    integer,intent(out) :: ifail
    real(kind=8),dimension(n_variables),intent(out) :: solution
    !> declare internal variables
    real(kind=8),dimension(n_variables*n_stages) :: differentials
    !> internal variables
    integer :: iteration !< number of iterations
    real(kind=8) :: error,dt_new !< error and new time step
    real(kind=8),dimension(n_variables) :: solution_low_order

    !> initialise counter to zero, error to 2 and copy time step
    iteration = 0
    dt_new = dt

    !> compute first step
    !> compute runge-kutta differential
    call compute_runge_kutta_differentials(compute_rhs,fields,n_variables,&
         n_int_parameters,n_real_parameters,t,dt,solution_old,&
         int_parameters,real_parameters,differentials,ifail)
    !> compute solution
    call compute_runge_kutta_solution(n_variables,solution_old,&
       differentials,solution,solution_low_order)
    !> compute error
    error = compute_base_error(n_variables,tolerances,&
         solution,solution_low_order)!< base error
    !> check if a user selected method has to be used
    if(present(compute_user_err)) error = max(error,&
         abs(compute_user_err(fields,n_variables,n_int_parameters,&
         n_real_parameters,t,solution,solution_low_order,&
         int_parameters,real_parameters)))
    !> if good error compute larger time step
    if(error.lt.1.d0) call compute_time_step_shampine(dt_new,&
         error,[error_parameters(1),error_parameters(3)])

    !> loop on the error control
    do while(error.ge.1.d0 .and. iteration.le.maximum_iteration)
       !> compute runge-kutta differentials
       call compute_runge_kutta_differentials(compute_rhs,fields,n_variables,&
            n_int_parameters,n_real_parameters,t,dt,solution_old,&
            int_parameters,real_parameters,differentials,ifail)
       !> compute solution
       call compute_runge_kutta_solution(n_variables,solution_old,&
          differentials,solution,solution_low_order)
       !> compute error
       error = compute_base_error(n_variables,tolerances,&
            solution,solution_low_order)!< base error
       !> check if a user selected method has to be used
       if(present(compute_user_err)) error = max(error,&
            abs(compute_user_err(fields,n_variables,n_int_parameters,&
            n_real_parameters,t,solution,solution_low_order,&
            int_parameters,real_parameters)))
       !> compute new time step
       call compute_time_step_shampine(dt_new,error,&
              error_parameters(1:2))
    enddo

    !> check of the time step is larger than the stop time
    if((t+dt_new).gt.t_stop) dt_new = t_stop-t
    !> overwrite time step
    dt = dt_new

  end subroutine runge_kutta_order_error_control_dt

  !> this subroutine implement a runge kutta integrator with
  !> adaptative integration step. Rules for adaptation step
  !> integration has to be provided by the user. This is
  !> basically a feed-forward integration step controller.
  !> inputs:
  !>   compute_rhs:       (procedure) subroutine for computing
  !>                      the ODE(s) right hand side
  !>   compute_dt:        (procedure) function for computing
  !>                      the new time step
  !>   fields:            (fields_base) jorek fields structure
  !>   n_variables:       (integer) number of varibales
  !>   n_int_parameters:  (integer) number of integer parameters
  !>   n_real_parameters: (integer) number of real parameters
  !>   t:                 (real8) integration variables
  !>   t_stop:            (real8) stop time
  !>   dt:                (real8) integration step
  !>   solution_old:      (real8)(n_variables) old solution
  !>   int_parameters:    (integer)(n_integer_parameters)
  !>                      integer parameters
  !>   real_parameters:   (real8)(n_real_parameters)
  !>                      real parameters
  !> outputs:
  !>   dt:       (real8) new integration step
  !>   solution: (n_variables) runge kutta step solution
  !>   ifail:    (integer) if 0 the integration failed
  subroutine runge_kutta_adaptative_dt(compute_rhs,compute_dt,&
       fields,n_variables,n_int_parameters,n_real_parameters,&
       t,t_stop,dt,solution_old,int_parameters,real_parameters,&
       solution,ifail)
    !> load modules
    use mod_fields, only: fields_base
    implicit none
    !> declare input-output variables
    real(kind=8),intent(inout) :: dt
    !> delcare inputs
    procedure(compute_runge_kutta_rhs) :: compute_rhs
    procedure(adapt_time_step) :: compute_dt
    class(fields_base),intent(in) :: fields
    integer,intent(in) :: n_variables,n_int_parameters,n_real_parameters
    real(kind=8),intent(in) :: t,t_stop
    real(kind=8),dimension(n_variables),intent(in) :: solution_old
    integer,dimension(n_int_parameters),intent(in) :: int_parameters
    real(kind=8),dimension(n_real_parameters),intent(in) :: real_parameters
    !> declare outputs
    integer,intent(out) :: ifail
    real(kind=8),dimension(n_variables),intent(out) :: solution
    !> declare internal variables
    real(kind=8),dimension(n_variables*n_stages) :: differentials

    !> compute the new time step
    dt = compute_dt(fields,n_variables,n_int_parameters,&
          n_real_parameters,t,dt,solution,int_parameters,&
          real_parameters)

    !> check if the time step is not too large
    if((t+dt) .gt. t_stop) dt = t_stop - dt

    !> compute runge-kutta differentials
    call compute_runge_kutta_differentials(compute_rhs,fields,n_variables,&
         n_int_parameters,n_real_parameters,t,dt,solution_old,&
         int_parameters,real_parameters,differentials,ifail)

    !> compute runge-kutta solutions
    call compute_runge_kutta_solution(n_variables,solution_old,&
         differentials,solution)

  end subroutine runge_kutta_adaptative_dt  
  
  !> this subroutine implement a runge kutta integrator with
  !> fixed integration step
  !> inputs:
  !>   compute_rhs:       (procedure) subroutine for computing
  !>                      the ODE(s) right hand side
  !>   fields:            (fields_base) jorek fields structure
  !>   n_variables:       (integer) number of varibales
  !>   n_int_parameters:  (integer) number of integer parameters
  !>   n_real_parameters: (integer) number of real parameters
  !>   t:                 (real8) integration variables
  !>   dt:                (real8) integration step
  !>   solution_old:      (real8)(n_variables) old solution
  !>   int_parameters:    (integer)(n_integer_parameters)
  !>                      integer parameters
  !>   real_parameters:   (real8)(n_real_parameters)
  !>                      real parameters
  !> outputs:
  !>   solution: (n_variables) runge kutta step solution
  !>   ifail:    (integer) if 0 the integration failed
  subroutine runge_kutta_fixed_dt(compute_rhs,fields,n_variables,&
       n_int_parameters,n_real_parameters,t,dt,solution_old,&
       int_parameters,real_parameters,solution,ifail)
    !> load modules
    use mod_fields, only: fields_base
    implicit none
    !> delcare inputs
    procedure(compute_runge_kutta_rhs) :: compute_rhs
    class(fields_base),intent(in) :: fields
    integer,intent(in) :: n_variables,n_int_parameters,n_real_parameters
    real(kind=8),intent(in) :: t,dt
    real(kind=8),dimension(n_variables),intent(in) :: solution_old
    integer,dimension(n_int_parameters),intent(in) :: int_parameters
    real(kind=8),dimension(n_real_parameters),intent(in) :: real_parameters
    !> declare outputs
    integer,intent(out) :: ifail
    real(kind=8),dimension(n_variables),intent(out) :: solution
    !> declare internal variables
    real(kind=8),dimension(n_variables*n_stages) :: differentials

    !> compute runge-kutta differentials
    call compute_runge_kutta_differentials(compute_rhs,fields,n_variables,&
         n_int_parameters,n_real_parameters,t,dt,solution_old,&
         int_parameters,real_parameters,differentials,ifail)

    !> compute runge-kutta solutions
    call compute_runge_kutta_solution(n_variables,solution_old,&
         differentials,solution)
    
  end subroutine runge_kutta_fixed_dt

  !> This procedure computes the Runge-Kutta differentials
  !> for each stage. Default Runge-Kutta: Cash-Karp 4(5)
  !> inputs:
  !>   compute_rhs:       (procedure) subroutine for computing
  !>                      the ODEs right hand side
  !>   n_variables:       (integer) number of variables (ODEs)
  !>   n_int_parameters:  (integer) number of integer parameters
  !>   n_real_parameters: (integer) number of real parameters
  !>   t:                 (real8) coordinate of solution old
  !>   dt:                (real8) integration step 
  !>   solution_old:      (real8)(n_variables) initial solution
  !>   int_parameters:    (integer)(n_int_parameters) integer parameters
  !>   real_parameters:   (real8)(n_real_parameters) real parameters
  !>   fields:            (fields_base) fields for particle tracking
  !> outputs:
  !>   differentials: (real8)(n_varibales*(n_stages+1)) differentials
  !>   ifail:               (integer) if 0 integration failed
  subroutine compute_runge_kutta_differentials(compute_rhs,&
       fields,n_variables,n_int_parameters,n_real_parameters,&
       t,dt,solution_old,int_parameters,real_parameters,&
       differentials,ifail)
    !> load module
    use mod_fields, only: fields_base
    implicit none
    !> step coefficients
    real(kind=8),dimension(5),parameter :: A_vect=[2.d-1,3.d-1,6.d-1,1.d0,8.75d-1]
    !> derivatives coefficients
    real(kind=8),dimension(15),parameter :: B_vect=[2.d-1,7.5d-2,2.25d-1,3.d-1,&
         -9.d-1,1.2d0,-2.037037037037037d-1,2.5d0,-2.592592592592593d0,&
         1.296296296296296d0,2.949580439814815d-2,3.41796875d-1,4.159432870370371d-2,&
         4.003454137731481d-1,6.1767578125d-2]
    !> declare input variables
    procedure(compute_runge_kutta_rhs) :: compute_rhs
    integer,intent(in) :: n_variables,n_int_parameters,n_real_parameters
    integer,dimension(n_int_parameters),intent(in) :: int_parameters
    real(kind=8),intent(in) :: t,dt
    real(kind=8),dimension(n_variables),intent(in) :: solution_old
    real(kind=8),dimension(n_real_parameters),intent(in) :: real_parameters
    class(fields_base),intent(in) :: fields
    !> declare output variable
    integer,intent(out) :: ifail
    real(kind=8),dimension(n_variables*n_stages),intent(out) :: differentials
    !> declare internal variables
    integer :: i,j !< indexes
    integer :: counter !< stage counter
    real(kind=8) :: t_new !< new coordinate
    real(kind=8),dimension(n_variables) :: deltas !< reduction of derivative stages
    !> interface of the subroutine compute derivatives

    counter = 0 !< initialise counter to 0    
    deltas = 0.d0  !< initialise deltas to zero

    !> compute the first derivatives
    call compute_rhs(fields,n_variables,n_int_parameters,n_real_parameters,&
         t,solution_old,solution_old,int_parameters,real_parameters,&
         differentials(1:n_variables),ifail)
    
    !> loop on the number of stages 
    do i=1,n_stages-1
       if(ifail.eq.0) exit !< exit if integration failed
       t_new = t + A_vect(i)*dt !< update the time
       !> loop on the previous stages for reduction
       do j=1,i
          deltas = deltas + B_vect(counter+j)*&
               differentials((j-1)*n_variables+1:j*n_variables)
       enddo
       counter = counter + i !< update counter
       !> computing the new derivatives

       call compute_rhs(fields,n_variables,n_int_parameters,&
            n_real_parameters,t_new,solution_old,solution_old+dt*deltas,&
            int_parameters,real_parameters,&
            differentials(i*n_variables+1:(i+1)*n_variables),ifail)
       deltas = 0.d0 !< re-initialise the deltas variable
    enddo
    differentials = dt*differentials !< comute differentials from derivatives

  end subroutine compute_runge_kutta_differentials

  !> This procedure computes the runge kutta solution of five different orders
  !> for embedded runge-kutta. This procedure allows to control both
  !> integration step and order. Default: runge-kutta 4(5) cash-karp
  !> inputs:
  !>   n_variables:   (integer) number of variables
  !>   solution_old:  (real8)(n_varibales) old solution
  !>   differentials: (real8)(n_variables*(n_stages+1)) differentials
  !> outputs:
  !>   solution_1: (real8)(n_variables) highest order solution
  !>   solution_2: (real8)(n_variables) middle-high order solution
  !>   solution_3: (real8)(n_variables) middle order solution
  !>   solution_4: (real8)(n_variables) middle-low order solution
  !>   solution_5: (real8)(n_variables) lowest order solution
  pure subroutine compute_runge_kutta_solution_5(n_variables,solution_old,&
       differentials,solution_1,solution_2,solution_3,solution_4,solution_5)
    !> coefficients for computing solutions
    real(kind=8),dimension(30),parameter :: C_vect=[9.788359788359788d-2,&
         0.d0,4.025764895330113d-1,2.104377104377105d-1,0.d0,2.891022021456804d-1,&
         1.021773726851852d-1,0.d0,3.839079034391534d-1,2.445927372685185d-1,&
         1.932198660714286d-2,2.5d-1,3.518518518518519d-1,0.d0,-3.703703703703703d-1,&
         1.018518518518519d0,0.d0,0.d0,-1.5d0,2.5d0,0.d0,0.d0,0.d0,0.d0,&
         1.d0,0.d0,0.d0,0.d0,0.d0,0.d0]
    !> declare input variables
    integer,intent(in) :: n_variables
    real(kind=8),dimension(n_variables),intent(in) :: solution_old
    real(kind=8),dimension(n_variables*n_stages),intent(in) :: differentials
    !> declare output variables
    real(kind=8),dimension(n_variables),intent(out) :: solution_1,solution_2,&
         solution_3,solution_4,solution_5
    !> delcare internal variables
    integer :: i !< indexes

    !> initialise solutions
    solution_1 = solution_old
    solution_2 = solution_old
    solution_3 = solution_old
    solution_4 = solution_old
    solution_5 = solution_old

    !> loop on the number of stages
    do i=1,n_stages
       !> compute solutions
       solution_1 = solution_1 + &
            C_vect(i)*differentials(n_variables*(i-1)+1:i*n_variables)
       solution_2 = solution_2 + &
            C_vect(n_stages+i)*differentials(n_variables*(i-1)+1:i*n_variables)
       solution_3 = solution_3 + &
            C_vect(2*n_stages+i)*differentials(n_variables*(i-1)+1:i*n_variables)
       solution_4 = solution_4 + &
            C_vect(3*n_stages+i)*differentials(n_variables*(i-1)+1:i*n_variables)
       solution_5 = solution_5 + &
            C_vect(4*n_stages+1)*differentials(n_variables*(i-1)+1:i*n_variables)
    enddo
    
  end subroutine compute_runge_kutta_solution_5

  !> This procedure computes the runge kutta solution for two different orders
  !> for embedded runge-kutta. This procedure allows integration step control.
  !> Defaul: runge-kutta 4(5) cash-karp
  !> inputs:
  !>   n_variables:   (integer) number of variables
  !>   solution_old:  (real8)(n_variables) old solution
  !>   differentials: (real8)(n_variables*(n_stages+1)) differentials
  !> outputs:
  !>   solution_1: (real8)(n_variables) highest order solution
  !>   solution_2: (real8)(n_variables) lowest order solution
  pure subroutine compute_runge_kutta_solution_2(n_variables,solution_old,&
       differentials,solution_1,solution_2)
    !> coefficients for computing solutions
    real(kind=8),dimension(12),parameter :: C_vect=[9.788359788359788d-2,&
         0.d0,4.025764895330113d-1,2.104377104377105d-1,0.d0,2.891022021456804d-1,&
         1.021773726851852d-1,0.d0,3.839079034391534d-1,2.445927372685185d-1,&
         1.932198660714286d-2,2.5d-01]
    !> declare input variables
    integer,intent(in) :: n_variables
    real(kind=8),dimension(n_variables),intent(in) :: solution_old
    real(kind=8),dimension(n_variables*n_stages),intent(in) :: differentials
    !> declare output variables
    real(kind=8),dimension(n_variables),intent(out) :: solution_1,solution_2
    !> delcare internal variables
    integer :: i !< index

    !> initialise solution
    solution_1 = solution_old
    solution_2 = solution_old
    
    !> loop on the number of stages
    do i=1,n_stages
       !> compute solution
       solution_1 = solution_1 + &
            C_vect(i)*differentials(n_variables*(i-1)+1:i*n_variables)
       solution_2 = solution_2 + &
            C_vect(n_stages+1)*differentials(n_variables*(i-1)+1:i*n_variables)
    enddo
    
  end subroutine compute_runge_kutta_solution_2
  
  !> This procedure computes the runge kutta solution for one specific order.
  !> Default: runge-kutta 4(5) cash-karp
  !> inputs:
  !>   n_variables:   (integer) number of variables
  !>   solution_old:  (real8)(n_variables) old solution
  !>   differentials: (real8)(n_variables*(n_stages+1)) differentials
  !> outputs:
  !>   solution: (real8)(n_variables) runge-kutta solution
  pure subroutine compute_runge_kutta_solution_1(n_variables,solution_old,&
       differentials,solution)
    !> coefficients for computing solutions
    real(kind=8),dimension(6),parameter :: C_vect=[9.788359788359788d-2,&
         0.d0,4.025764895330113d-1,2.104377104377105d-01,0.d0,2.891022021456804d-1]
    !> declare input variables
    integer,intent(in) :: n_variables
    real(kind=8),dimension(n_variables),intent(in) :: solution_old
    real(kind=8),dimension(n_stages*n_variables),intent(in) :: differentials
    !> declare output variables
    real(kind=8),dimension(n_variables),intent(out) :: solution
    !> declare internal variables
    integer :: i !< index
    
    solution = solution_old !< initialise solution
    !> loop on the stages
    do i=1,n_stages
       solution = solution + C_vect(i)*differentials(n_variables*(i-1)+1:i*n_variables)
    enddo
       
  end subroutine compute_runge_kutta_solution_1

  !> This procedure compute the base error between two solutions as
  !> max(abs). WARNING: tolerance must have the size of the variables
  !> inputs:
  !>   n_variables: (integer) number of variables
  !>   tolerances:  (real8)(n_variables) tolerances
  !>   solution_1:  (real8)(n_variables) first solution
  !>   solution_2:  (real8)(n_variables) second solution
  !> outputs:
  !>   error: (real8) error estimate
  pure function compute_base_error(n_variables,tolerances,solution_1,&
       solution_2) result(error)
    !> declare intput varibales
    integer,intent(in) :: n_variables
    real(kind=8),dimension(n_variables),intent(in) :: solution_1,solution_2
    real(kind=8),dimension(n_variables),intent(in) :: tolerances
    !> declare output_variables
    real(kind=8) :: error

    !> compute the error in infinite norm (most restrictive)
    error = maxval(abs((solution_1-solution_2)/tolerances))
    
  end function compute_base_error
  
  !> This procedure estimate the new time step using the shimpine method
  !> inputs:
  !>   dt:         (real8) old time step
  !>   error:      (real8) integration error
  !<   parameters: (real8)(2) 1:safety factor, 2:exponent
  !> outputs:
  !>   dt: (real8) new time step
  pure subroutine compute_time_step_shampine(dt,error,parameters)
    !> declare input/ourputs
    real(kind=8),intent(inout) :: dt
    !> declare inputs:
    real(kind=8),intent(in) :: error
    real(kind=8),dimension(2),intent(in) :: parameters

    !> compute new time step
    dt = dt*parameters(1)*(error**(-1.d0*parameters(2)))
    
  end subroutine compute_time_step_shampine
  
end module mod_runge_kutta

