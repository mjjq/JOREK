!> This module contains all procedures and variables required
!> for solving a seto of ordinary differential equatons using
!>a Runge-Kutta scheme. By default, the Runge-Kutta 4(5)
!> Cash-Karp is implemented.
module mod_runge_kutta

  private !< set all as private
  public n_stages
  public compute_runge_kutta_derivatives
  public compute_runge_kutta_solution

  !> declare module parameters
  integer,parameter :: n_stages=6 !< number of stages

  !> interfaces

  !> this is is a common interface for all procedures computing
  !> solutions from runge-kutta differentials
  interface compute_runge_kutta_solution
     module procedure compute_runge_kutta_solution_1,&
          compute_runge_kutta_solution_2,compute_runge_kutta_solution_5
  end interface compute_runge_kutta_solution
  
  interface
     !> This is the interface for procedures computing the ODE
     !> right hanf side to be used in compute_runge_kutta_derivatives
     !> inputs:
     !>   n_variables:  (integer) number of variables
     !>   t:            (real8) integration coordinate
     !>   solution_old: (real8)(n_variables) old solution
     !>   deltas:       (real8)(n_variables) sum of previous satages
     !> outputs
     !>   derivatives:  (real8)(n_variables) new derivatives
     !>   ifail:        (integer) if 0 derivative calculation failed
     subroutine compute_derivtives_runge_kutta(n_variables,t,&
          solution_old,deltas,derivatives,ifail)
       !> declare input variables
       integer,intent(in) :: n_variables
       real(kind=8),intent(in) :: t
       real(kind=8),dimension(n_variables),intent(in) :: solution_old,deltas
       !> declare output variables
       integer,intent(out) :: ifail
       real(kind=8),dimension(n_variables),intent(out) :: derivatives
     end subroutine compute_derivtives_runge_kutta
  end interface
  
contains

  !> This procedure computes the Runge-Kutta differentials
  !> for each stage. Default Runge-Kutta: Cash-Karp 4(5)
  !> inputs:
  !>   n_variables:         (integer) number of variables (ODEs)
  !>   t:                   (real8) coordinate of solution old
  !>   dt:                  (real8) integration step 
  !>   solution_old:        (real8)(n_variables) initial solution
  !> outputs:
  !>   differentials: (real8)(n_varibales*(n_stages+1)) differentials
  !>   ifail:               (integer) if 0 integration failed
  pure subroutine compute_runge_kutta_derivatives(n_variables,&
       t,dt,solution_old,differentials,ifail)
    !> step coefficients
    real(kind=8),dimension(6),parameter :: A_vect=[2.d-1,3.d-1,6.d-1,1.d0,8.75d-1]
    !> derivatives coefficients
    real(kind=8),dimension(15),parameter :: B_vect=[2.d-1,7.5d-2,2.25d-1,3.d-1,&
         -9.d-1,1.2d0,-2.037037037037037d-1,2.5d0,-2.592592592592593d0,&
         1.296296296296296d0,2.949580439814815d-2,3.41796875d-1,4.159432870370371d-2,&
         4.003454137731481d-1,6.1767578125d-2]
    !> declare input variables
    integer,intent(in) :: n_variables
    real(kind=8),intent(in) :: t,dt
    real(kind=8),dimension(n_variables),intent(in) :: solution_old
    !> declare output variable
    integer,intent(out) :: ifail
    real(kind=8),dimension(n_variables*n_stages),intent(out) :: differentials
    !> declare internal variables
    integer :: i,j !< indexes
    integer :: counter=0 !< stage counter
    real(kind=8) :: t_new !< new coordinate
    real(kind=8),dimension(n_variables) :: deltas=0.d0 !< reduction of derivative stages
    
    !> compute the first derivatives
    call compute_derivatives_runge_kutta(n_variables,t,solution_old,&
         deltas,differentials(1:n_variables),ifail)
    
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
       call compute_derivatives_runge_kutta(n_variables,t_new,solution_old,&
            dt*deltas,differentials(i*n_variables+1:(i+1)*n_variables),ifail)
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
  pure function compute_runge_kutta_solution_5(n_variables,solution_old,&
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
            C_vect(i)*differentials(n_variabls*(i-1)+1:i*n_variables)
       solution_2 = solution_2 + &
            C_vect(n_stages+i)*differentials(n_variables*(i-1)+1:i*n_variables)
       solution_3 = solution_3 + &
            C_vect(2*n_stages+i)*differentials(n_variables*(i-1)+1:i*n_variables)
       solution_4 = solution_4 + &
            C_vect(3*n_stages+i)*differentials(n_variables*(i-1)+1:i*n_variables)
       solution_5 = solution_5 + &
            C_vect(4*n_stages+1)*differentials(n_variables*(i-1)+1:i*n_variables)
    enddo
    
  end function compute_runge_kutta_solution_5

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
  pure function compute_runge_kutta_solution_2(n_variables,solution_old,&
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
    do i=1:n_stages
       !> compute solution
       solution_1 = solution_1 + &
            C_vect(i)*differentials(n_variables*(i-1)+1:i*n_variables)
       solution_2 = solution_2 + &
            C_vect(n_stages+1)*differentials(n_variables*(i-1)+1:i*n_variables)
    enddo
    
  end function compute_runge_kutta_solution_2
  
  !> This procedure computes the runge kutta solution for one specific order.
  !> Default: runge-kutta 4(5) cash-karp
  !> inputs:
  !>   n_variables:   (integer) number of variables
  !>   solution_old:  (real8)(n_variables) old solution
  !>   differentials: (real8)(n_variables*(n_stages+1)) differentials
  !> outputs:
  !>   solution: (real8)(n_variables) runge-kutta solution
  pure function compute_runge_kutta_solution_1(n_variables,solution_old,&
       differentials,solution)
    !> coefficients for computing solutions
    real(kind=6),dimension(6),parameter :: C_vect=[9.788359788359788d-2,&
         0.d0,4.025764895330113d-1,2.104377104377105d-01,0.d0,2.891022021456804d-1]
    !> declare input variables
    integer,intent(in) :: n_variables    integer,intent(in) :: n_variables
    real(kind=8),dimension(n_variables),intent(in) :: solution_old
    real(kind=8),dimension(n_stages*n_variables),intent(in) :: differentials
    !> declare output variables
    real(kind=8),dimension(n_variables),intent(out) :: solution
    !> declare internal variables
    integer :: i !< index
    
    solution = solution_old !< initialise solution
    !> loop on the stages
    do i=1:n_stages
       solution = solution + C_vect(i)*differentials(n_variables*(i-1)+1:i*n_variables)
    enddo
       
  end function compute_runge_kutta_solution_1
  
end module mod_runge_kutta

