!> module mod_math_operators_test contains variables and
!> and procedures for testing the procedures in 
!> mod_math_operators
module mod_math_operators_test
use fruit
implicit none

private
public :: run_fruit_math_operators

!> Variables --------------------------------------------
integer,parameter :: zero_r4=real(0.d0,kind=4)
integer,parameter :: n_vectors=4
real*4,parameter  :: tol_r4=real(1.0d-1,kind=4)
real*8,parameter  :: tol_r8=2.5d-10
!> vectors for testing the vector product
real*8,dimension(2) :: a_interval_r8=(/-3.5d1,5.4d1/)
real*8,dimension(2) :: b_interval_r8=(/5.4d1,1.15d2/)
real*4,dimension(3,n_vectors) :: vec_prod_a_r4,vec_prod_b_r4
real*8,dimension(3,n_vectors) :: vec_prod_a_r8,vec_prod_b_r8

contains

!> Test basket ------------------------------------------
!> test basket for setting-up, executing and tearing-down
!> all test feature
subroutine run_fruit_math_operators()
  implicit none
  write(*,'(/A)') "  ... setting-up: math operators tests"
  call setup
  write(*,'(/A)') "  ... running: math operators tests"
  call test_cross_product
  write(*,'(/A)') "  ... tearing-down: math operators tests"
  call teardown
end subroutine run_fruit_math_operators

!> Set-up and tear-down ---------------------------------
!> Set-up of the math operator test features
subroutine setup()
  use mod_gnu_rng, only: gnu_rng_interval
  implicit none
  integer :: ii

  !> generate vectors for cross product tests
  call gnu_rng_interval(3,n_vectors,a_interval_r8,vec_prod_a_r8)
  call gnu_rng_interval(3,n_vectors,b_interval_r8,vec_prod_b_r8)
  vec_prod_a_r4 = real(vec_prod_a_r8,kind=4)
  vec_prod_b_r4 = real(vec_prod_b_r8,kind=4)
!  do ii=1,n_vectors
!    vec_prod_b_r8(:,ii) = vector_rotation(angles_r8(ii),vec_prod_a_r8(:,ii))
!    vec_prod_b_r4(:,ii) = vector_rotation(angles_r4(ii),vec_prod_a_r4(:,ii))
!  enddo
end subroutine setup

!> tear-down the math operator test features
subroutine teardown()
  implicit none
  !> set all static vectors to zero
  vec_prod_a_r4 = zero_r4; vec_prod_b_r4 = zero_r4;
  vec_prod_a_r8 = 0.d0; vec_prod_b_r8 = 0.d0;
end subroutine teardown

!> Tests ------------------------------------------------
!> Test vector product both single and double precision
subroutine test_cross_product()
  use mod_math_operators, only: cross_product
  implicit none
  integer :: ii
  real*4,dimension(3) :: vec_c_r4
  real*8,dimension(3) :: vec_c_r8

  !> test single precision cross product
  do ii=1,n_vectors
    vec_c_r4 = cross_product(vec_prod_a_r4(:,ii),vec_prod_b_r4(:,ii))
    call assert_true((abs(dot_product(vec_prod_a_r4(:,ii),vec_c_r4)).lt.tol_r4),&
    "Error math operators cross product (float): vectors a and c not orthogonal")
    call assert_true((abs(dot_product(vec_prod_b_r4(:,ii),vec_c_r4)).lt.tol_r4),&
    "Error math operators cross product (float): vectors b and c not orthogonal")
  enddo
  do ii=1,n_vectors
    vec_c_r8 = cross_product(vec_prod_a_r8(:,ii),vec_prod_b_r8(:,ii))
    call assert_equals(dot_product(vec_prod_a_r8(:,ii),vec_c_r8),0.d0,tol_r8,&
    "Error math operators cross product (double): vectors a and c not orthogonal")
    call assert_equals(dot_product(vec_prod_b_r8(:,ii),vec_c_r8),0.d0,tol_r8,&
    "Error math operators cross product (double): vectors b and c not orthogonal")
  enddo
end subroutine test_cross_product

!>-------------------------------------------------------

end module mod_math_operators_test
