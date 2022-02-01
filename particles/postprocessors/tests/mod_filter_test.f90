!> mod_filter_test contains variables an procedures used for
!> testing the filter class. The filter class is an abstract
!> class hence the filter_unity class is used instead.
module mod_filter_test
use fruit
implicit none

private
public :: run_fruit_filter

!> Variables and and parameters -------------------------------
integer,parameter                                  :: n_dimensions_sol=2
integer,parameter                                  :: n_positions_sol=13
real*8,parameter                                   :: tol_real8=5.d-16
integer,dimension(n_dimensions_sol),parameter      :: filter_shape_sol=(/5,10/)
real*8,dimension(n_dimensions_sol),parameter       :: pos_lowbnd=(/-3.d1,1.d-2/)
real*8,dimension(n_dimensions_sol),parameter       :: pos_uppbnd=(/3.d1,5.d0/)
real*8,dimension(n_dimensions_sol,n_positions_sol) :: positions_sol
real*8,dimension(n_positions_sol)                  :: weights_sol

!> Interfaces -------------------------------------------------
contains
!> Fruit basket -----------------------------------------------
!> fruit basket containing all set-up, test and tear-down 
!> procedure to be executed
subroutine run_fruit_filter()
  implicit none
  write(*,*) "  ... setting-up filter tests"
  call setup
  write(*,*) "  ... running: filter tests"
  call test_filter_allocation_deallocation
  call test_initialisation_filter_unity
  call test_compute_weights_filter_unity
  call test_compute_weights_filter_unity_vectorial
  write(*,*) "  ... tearing-down: filter tests"
end subroutine run_fruit_filter

!> Set-up and tear-down ---------------------------------------
!> initiliase test features
subroutine setup()
  use mod_gnu_rng, only: gnu_rng_interval
  implicit none
  !> variables
  integer :: ii
  !> initilaise positions and weights
  do ii=1,n_positions_sol
    call gnu_rng_interval(n_dimensions_sol,pos_lowbnd,pos_uppbnd,positions_sol(:,ii))
  enddo 
  weights_sol = 1.d0
end subroutine setup

!> Tests ------------------------------------------------------
!> test the filter allocation and deallocation methods
subroutine test_filter_allocation_deallocation()
  use mod_assert_equals_tools, only: assert_equals_allocatable_arrays
  use mod_filter_unity,        only: filter_unity
  implicit none
  !> variables
  type(filter_unity) :: filter_test
  !> tests allocation
  call filter_test%allocate_filter(n_dimensions_sol)
  call assert_equals(filter_test%n_dimensions,n_dimensions_sol,&
  "Error allocate filter: n_dimensions mistmatech!")
  call assert_equals_allocatable_arrays(n_dimensions_sol,&
  filter_test%stencil_shape,0,"Error allocate filter: stencil_shape")
  !> test deallocation
  call filter_test%deallocate_filter
  call assert_equals(filter_test%n_dimensions,0,&
  "Error deallocate filter: n_dimensions not reset to 0!")
  call assert_false(allocated(filter_test%stencil_shape),&
  "Error deallocate filter: stencil shape not deallocated!")
end subroutine test_filter_allocation_deallocation

!> test initialisation of filter unity
subroutine test_initialisation_filter_unity()
  use mod_filter_unity, only: filter_unity
  implicit none
  !> variables
  type(filter_unity) :: filter_test
  !> test initialisation without stencil shape
  call filter_test%init_filter(n_dimensions_sol)
  call assert_equals(filter_test%n_dimensions,n_dimensions_sol,&
  "Error filter unity initialisation: n dimensions mismatch")
  call filter_test%deallocate_filter
  !> test initialisation with stencil shape
  call filter_test%init_filter(n_dimensions_sol,filter_shape_sol)
  call assert_equals(filter_test%n_dimensions,n_dimensions_sol,&
  "Error filter unity initialisation stencil shape: n dimensions mismatch")
end subroutine test_initialisation_filter_unity 

!> test compute weights filter unity
subroutine test_compute_weights_filter_unity()
  use mod_filter_unity, only: filter_unity
  implicit none
  !> variables
  integer :: ii
  real*8,dimension(n_positions_sol) :: weights
  type(filter_unity) :: filter_test
  !> initialisation
  call filter_test%init_filter(n_dimensions_sol)
  !> test
  do ii=1,n_positions_sol
    call filter_test%compute_filter_from_position(positions_sol(:,ii),weights(ii))
  enddo
  call assert_equals(weights,weights_sol,n_positions_sol,tol_real8,&
  "Error filter unity compute filter form positions: filter weights mistmatch!")
  !> deallocation
  call filter_test%deallocate_filter
end subroutine test_compute_weights_filter_unity

!> test compute weights filter unity: vectorial version
subroutine test_compute_weights_filter_unity_vectorial()
  use mod_filter_unity, only: filter_unity
  implicit none
  !> variables
  real*8,dimension(n_positions_sol) :: weights
  type(filter_unity) :: filter_test
  !> initialisation
  call filter_test%init_filter(n_dimensions_sol)
  !> test
  call filter_test%compute_filter_from_position_vectorial(n_positions_sol,positions_sol,weights)
  call assert_equals(weights,weights_sol,n_positions_sol,tol_real8,&
  "Error filter unity compute filter form positions (vectorial): filter weights mistmatch!")
  !> deallocation
  call filter_test%deallocate_filter
end subroutine test_compute_weights_filter_unity_vectorial

!> Tools ------------------------------------------------------
!>-------------------------------------------------------------
end module mod_filter_test
