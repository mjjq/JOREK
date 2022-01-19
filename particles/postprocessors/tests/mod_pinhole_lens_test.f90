!> mod_pinhole_lens_test contains all variable and procedures required
!> for testing the pinhole lens class.
module mod_pinhole_lens_test
  use fruit
  implicit none

  private
  public :: run_fruit_pinhole_lens

!> Variables --------------------------------------------------
integer,parameter :: n_x=3
integer,parameter :: n_samples=12453
real*8,parameter  :: tol_r8=5.d-16
real*8,dimension(n_samples),parameter :: ones_r8=1.d0
real*8,dimension(n_x),parameter :: center_uppbnd=(/4.5d3,4.67d2,-2.34d-1/)
real*8,dimension(n_x),parameter :: center_lowbnd=(/-1.d2,2.45d1,-2.3d2/)
real*8,dimension(n_x)           :: center_sol
real*8,dimension(n_x,n_samples) :: x_sol

!> Interfaces -------------------------------------------------
contains
!> Fruit basket -----------------------------------------------
!> basket running the set-up, test and teard-down procedures
subroutine run_fruit_pinhole_lens()
  implicit none
  write(*,'(/A)') "  ... setting-up: lens tests"
  call setup
  write(*,'(/A)') "  ... running: lens tests"
  call test_pinhole_lens_sampling
  call test_pinhole_initialisation
  call test_pinhole_lens_pdf
  write(*,'(/A)') "  ... tearing-down: lens tests"
  call teardown
end subroutine run_fruit_pinhole_lens

!> Set-up and teard-down---------------------------------------
!> set-up test features
subroutine setup()
  use mod_gnu_rng, only: gnu_rng_interval
  implicit none
  integer :: ii
  !> initialize center array
  call gnu_rng_interval(n_x,center_lowbnd,center_uppbnd,center_sol)
  do ii=1,n_samples
    x_sol(:,ii) = center_sol
  enddo
end subroutine setup

!> tear-down test features
subroutine teardown()
  implicit none
  center_sol = 0.d0; x_sol = 0.d0;
end subroutine teardown

!> Tests ------------------------------------------------------
!> Test pinhole initialisation
subroutine test_pinhole_initialisation()
  use mod_pinhole_lens, only: pinhole_lens
  implicit none
  type(pinhole_lens) :: pinhole
  call pinhole%init_pinhole(n_x,center_sol)
  call assert_equals(pinhole%n_x,n_x,&
  "Error pinhole lens initialisation: n_x mismatch!")
  call assert_true(allocated(pinhole%center),&
  "Error pinhole lens initialisation: center not allocated!")
  call assert_equals(pinhole%center,center_sol,n_x,tol_r8,&
  "Error pinhole lens initialisation: center mismatch!")
  call pinhole%deallocate_lens
end subroutine test_pinhole_initialisation
 
!> Test pinhole sampling routine
subroutine test_pinhole_lens_sampling()
  use mod_pinhole_lens, only: pinhole_lens
  implicit none
  type(pinhole_lens)  :: pinhole
  real*8,dimension(n_x,n_samples) :: x_loc
  call pinhole%init_pinhole(n_x,center_sol)
  call pinhole%lens_sampling(n_samples,x_loc)
  call assert_equals(x_loc,x_sol,n_x,n_samples,tol_r8,&
  "Error pinhole lens sampling: positions mismatch!")
  call pinhole%deallocate_lens
end subroutine test_pinhole_lens_sampling

!> Test pinhole pdf routine
subroutine test_pinhole_lens_pdf()
  use mod_pinhole_lens, only: pinhole_lens
  implicit none
  type(pinhole_lens)  :: pinhole
  real*8,dimension(n_samples) :: pdf_loc
  call pinhole%init_pinhole(n_x,center_sol)
  call pinhole%lens_pdf(n_samples,x_sol,pdf_loc)
  call assert_equals(pdf_loc,ones_r8,n_samples,tol_r8,&
  "Error pinhole lens sampling: positions mismatch!")
  call pinhole%deallocate_lens
end subroutine test_pinhole_lens_pdf

!> Tools ------------------------------------------------------
!>-------------------------------------------------------------
end module mod_pinhole_lens_test
