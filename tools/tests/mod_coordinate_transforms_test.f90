!> mod_coordinate_transformrs_test contains variables and
!> procedures for testing the mod_coordinate_transforms
!> module procedures
module mod_coordinate_transforms_test
use fruit
implicit none

private
public :: run_fruit_coordinate_transforms

!> Variables--------------------------------------------
integer,parameter :: n_points=4                !< number of test positions
integer,parameter :: n_origins=4               !< number of sphere origins
real*8,parameter  :: tol_r8=1.d-14             !< tolerance double
real*4,parameter  :: tol_r4=real(1.d-6,kind=4) !< tolerance float
!> tolerance for calculatons
real*4,parameter  :: tol_calc_r4=real(5.d-5,kind=4)
real*8,parameter  :: tol_calc_r8=7.5d-14
real*4,parameter  :: zero_r4=real(0.d0,kind=4)
real*4,parameter  :: one_r4=real(1.d0,kind=4)
real*4,dimension(3),parameter :: zeros_r4=real((/0.d0,0.d0,0.d0/),kind=4)
real*8,dimension(3),parameter :: zeros_r8=(/0.d0,0.d0,0.d0/)
!> intervals for randomly chosing the first, second and third
!> position components
real*8,dimension(3),parameter :: x_lowbnd=(/-2.3d1,-3.2d2,-9.d-1/)
real*8,dimension(3),parameter :: x_uppbnd=(/4.23d2,1.45d1,7.50d1/)
!> intervals for randomly chosing the first, second and third
!> vector components
real*8,dimension(3),parameter :: a_lowbnd=(/-3.41d2,-4.67d1,-9.35d1/)
real*8,dimension(3),parameter :: a_uppbnd=(/6.75d1,8.70d1,2.43d2/)
real*4,dimension(3,n_points)  :: x_r4           !< set o positions
real*4,dimension(3,n_origins) :: origin_r4      !< set of origins
real*4,dimension(3,n_origins) :: T_r4,N_r4,B_r4 !< sphere directions
real*8,dimension(3,n_points)  :: x_r8           !< set o positions
real*8,dimension(3,n_origins) :: origin_r8      !< set of origins
real*8,dimension(3,n_origins) :: T_r8,N_r8,B_r8 !< sphere directions

!> Interfaces ------------------------------------------
!> function for testing basis orthonormality
interface test_orthonormality_basis
  module procedure test_orthonormality_basis_r4
  module procedure test_orthonormality_basis_r8
end interface test_orthonormality_basis

contains

!> Fruit test basket -----------------------------------
!> Test basket for executing set-up, tests and tear-down
subroutine run_fruit_coordinate_transforms()
  implicit none
  write(*,'(/A)') "  ... setting-up: coordinate transfroms tests"
  call setup
  write(*,'(/A)') "  ... running: coordinate transforms tests"
  call test_cartisian_tofrom_cylindrical_transform
  write(*,'(/A)') "  ... tearing-down: coordinate transforms tests"
  call teardown
end subroutine run_fruit_coordinate_transforms

!> Set-up and tear-down --------------------------------
!> Set-up test features common to all tests
subroutine setup()
  use mod_gnu_rng,        only: gnu_rng_interval
  use mod_math_operators, only: cross_product
  implicit none
  !> variables
  integer :: ii
  
  !> generate random positions (assume cartesian coord.)
  do ii=1,n_points
    call gnu_rng_interval(n_points,x_lowbnd,x_uppbnd,x_r8(:,ii))
  enddo
  !> generate random origins (assume cartesian coord.)
  do ii=1,n_origins
    call gnu_rng_interval(n_points,x_lowbnd,x_uppbnd,origin_r8(:,ii))
    !> generate random orthonormal directions (assume cartesian coord.)
    call gnu_rng_interval(n_points,a_lowbnd,a_uppbnd,T_r8(:,ii))
    call gnu_rng_interval(n_points,a_lowbnd,a_uppbnd,N_r8(:,ii))
    !> compute the vector for the ith origin
    T_r8(:,ii) = T_r8(:,ii)-origin_r8(:,ii)
    N_r8(:,ii) = N_r8(:,ii)-origin_r8(:,ii)
    !> compute basis
    T_r8(:,ii) = T_r8(:,ii)/norm2(T_r8(:,ii))
    N_r8(:,ii) = N_r8(:,ii) - (dot_product(T_r8(:,ii),N_r8(:,ii)))*T_r8(:,ii)
    N_r8(:,ii) = N_r8(:,ii)/norm2(N_r8(:,ii))
    B_r8(:,ii) = cross_product(T_r8(:,ii),N_r8(:,ii))
    B_r8(:,ii) = B_r8(:,ii)/norm2(B_r8(:,ii))
    !> test the correct generation of the orthonormal basis
    call test_orthonormality_basis(T_r8(:,ii),N_r8(:,ii),B_r8(:,ii))
  enddo

  !> convert to float precision
  x_r4 = real(x_r8,kind=4); origin_r4 = real(origin_r8,kind=8); 
  T_r4 = real(T_r8,kind=4); N_r4 = real(N_r8,kind=4); B_r4 = real(B_r8,kind=4);
  do ii=1,n_origins
    call test_orthonormality_basis(T_r4(:,ii),N_r4(:,ii),B_r4(:,ii))
  enddo
end subroutine setup

!> Clean-up all common test features
subroutine teardown()
  implicit none
  !> set all variables to 0
  x_r4 = zero_r4; origin_r4 = zero_r4;
  T_r4 = zero_r4; N_r4 = zero_r4; B_r4 = zero_r4;
  x_r8 = 0.d0; origin_r8 = 0.d0;
  T_r8 = 0.d0; N_r8 = 0.d0; B_r8 = 0.d0;
end subroutine teardown
!> Tests -----------------------------------------------
!> Test cartesian to cylindrical and cylindrical to 
!> cartesian transformations for single and double
!> precision function
subroutine test_cartisian_tofrom_cylindrical_transform()
  use mod_coordinate_transforms, only: cartesian_to_cylindrical
  use mod_coordinate_transforms, only: cylindrical_to_cartesian
  implicit none
  !> variables
  integer             :: ii
  real*4,dimension(3) :: x_cart_new_r4,x_cyl_r4
  real*8,dimension(3) :: x_cart_new_r8,x_cyl_r8
  !> test single procision
  do ii=1,n_points
    x_cyl_r4 = cartesian_to_cylindrical(x_r4(:,ii))
    x_cart_new_r4 = cylindrical_to_cartesian(x_cyl_r4)   
    call assert_equals(x_r4(:,ii)-x_cart_new_r4,zeros_r4,3,tol_calc_r4,&
    "Error test cartesian to cylindrical (float): x-cartesian mismatch!")
  enddo
  !> test double procision
  do ii=1,n_points
    x_cyl_r8 = cartesian_to_cylindrical(x_r8(:,ii))
    x_cart_new_r8 = cylindrical_to_cartesian(x_cyl_r8)
    call assert_equals(x_r8(:,ii)-x_cart_new_r8,zeros_r8,3,tol_calc_r8,&
    "Error test cartesian to cylindrical (double): x-cartesian mismatch!")
  enddo
end subroutine test_cartisian_tofrom_cylindrical_transform

!> Tools -----------------------------------------------
!> method for testing the orthonormality 
!> of a basis function. Single precision.
subroutine test_orthonormality_basis_r4(v1,v2,v3)
  implicit none
  real*4,dimension(3),intent(in) :: v1,v2,v3
  call assert_equals(dot_product(v1,v1),one_r4,tol_r4,&
  "Error basis orthonormality (float): v1 is not normalized!")
  call assert_equals(dot_product(v2,v2),one_r4,tol_r4,&
  "Error basis orthonormality (float): v2 is not normalized!")
  call assert_equals(dot_product(v3,v3),one_r4,tol_r4,&
  "Error basis orthonormality (float): v3 is not normalized!")
  call assert_equals(dot_product(v1,v2),zero_r4,tol_r4,&
  "Error basis orthonormality (float): v1 and v2 are not orthogonal!")
  call assert_equals(dot_product(v1,v3),zero_r4,tol_r4,&
  "Error basis orthonormality (float): v1 and v3 are not orthogonal!")
  call assert_equals(dot_product(v2,v3),zero_r4,tol_r4,&
  "Error basis orthonormality (float): v2 and v3 are not orthogonal!")
end subroutine test_orthonormality_basis_r4

!> method for testing the orthonormality 
!> of a basis function. Double precision.
subroutine test_orthonormality_basis_r8(v1,v2,v3)
  implicit none
  real*8,dimension(3),intent(in) :: v1,v2,v3
  call assert_equals(dot_product(v1,v1),1.d0,tol_r8,&
  "Error basis orthonormality (double): v1 is not normalized!")
  call assert_equals(dot_product(v2,v2),1.d0,tol_r8,&
  "Error basis orthonormality (double): v2 is not normalized!")
  call assert_equals(dot_product(v3,v3),1.d0,tol_r8,&
  "Error basis orthonormality (double): v3 is not normalized!")
  call assert_equals(dot_product(v1,v2),0.d0,tol_r8,&
  "Error basis orthonormality (double): v1 and v2 are not orthogonal!")
  call assert_equals(dot_product(v1,v3),0.d0,tol_r8,&
  "Error basis orthonormality (double): v1 and v3 are not orthogonal!")
  call assert_equals(dot_product(v2,v3),0.d0,tol_r8,&
  "Error basis orthonormality (double): v2 and v3 are not orthogonal!")
end subroutine test_orthonormality_basis_r8

!>------------------------------------------------------
end module mod_coordinate_transforms_test
