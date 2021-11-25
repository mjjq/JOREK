!> mod_coordinate_transformrs_test contains variables and
!> procedures for testing the mod_coordinate_transforms
!> module procedures
module mod_coordinate_transforms_test
use fruit
implicit none

private
public :: run_fruit_coordinate_transforms

!> Variables--------------------------------------------
integer,parameter           :: n_points=4
integer,parameter           :: n_origins=4
!> intervals for randomly chosing the first, second and third
!> position components
real*8,dimension(2),parameter :: x1_interval=(/-2.3d1,4.23d2/)
real*8,dimension(2),parameter :: x2_interval=(/-3.2d2,1.45d1/)
real*8,dimension(2),parameter :: x3_interval=(/-9.d-1,7.50d1/)
!> intervals for randomly chosing the first, second and third
!> vector components
real*8,dimension(2),parameter :: a1_interval=(/-3.41d2,6.75d1/)
real*8,dimension(2),parameter :: a2_interval=(/-4.67d1,8.70d1/)
real*8,dimension(2),parameter :: a3_interval=(/-9.35d1,2.43d2/)
real*8,dimension(n_points,3)  :: x       !< set o positions
real*8,dimension(n_origins,3) :: origins !< set of origins
real*8,dimension(n_origins,3) :: T,N,B   !< sphere directions
!>------------------------------------------------------

contains

!> Fruit test basket -----------------------------------
!> Test basket for executing set-up, tests and tear-down
subroutine run_fruit_coordinate_transforms()
  implicit none
  write(*,'(/A)') "  ... setting-up: coordinate transfroms tests"
  call setup
  write(*,'(/A)') "  ... running: coordinate transforms tests"
  write(*,'(/A)') "  ... tearing-down: coordinate transforms tests"
  call teardown
end subroutine run_fruit_coordinate_transforms

!> Set-up and tear-down --------------------------------
!> Set-up test features common to all tests
subroutine setup()
  implicit none
  x=0.d0
end subroutine setup

!> Clean-up all common test features
subroutine teardown()
  implicit none
  !> set all variables to 0
  x = 0.d0; origins = 0.d0;
  T = 0.d0; N = 0.d0; B = 0.d0;
end subroutine teardown
!> Procedure -------------------------------------------
!>------------------------------------------------------
end module mod_coordinate_transforms_test
