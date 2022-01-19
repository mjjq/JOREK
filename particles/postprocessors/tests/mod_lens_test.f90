!> mod_lens_test contains all variable and procedures required
!> for testing the lens class. Due to the fact that lens is an 
!> abstract class, the pinhole_lens class instead
module mod_lens_test
  use fruit
  use mod_pinhole_lens, only: pinhole_lens
  implicit none

  private
  public :: run_fruit_lens

!> Variables --------------------------------------------------
integer,parameter :: n_x=3
real*8,parameter :: tol_real8=5.d-16
real*8,dimension(n_x),parameter :: center_uppbnd=(/4.5d3,4.67d2,-2.34d-1/)
real*8,dimension(n_x),parameter :: center_lowbnd=(/-1.d2,2.45d1,-2.3d2/)
real*8,dimension(n_x) :: center_sol

!> Interfaces -------------------------------------------------
contains
!> Fruit basket -----------------------------------------------
!> basket running the set-up, test and teard-down procedures
subroutine run_fruit_lens()
  implicit none
  write(*,'(/A)') "  ... setting-up: lens tests"
  call setup()
  write(*,'(/A)') "  ... running: lens tests"
  write(*,'(/A)') "  ... tearing-down: lens tests"
  call teardown()
end subroutine run_fruit_lens

!> Set-up and teard-down---------------------------------------
!> set-up test features
subroutine setup()
  use mod_gnu_rng, only: gnu_rng_interval
  implicit none
  !> initialize center array
  call gnu_rng_interval(n_x,center_lowbnd,center_uppbnd,center_sol)
end subroutine setup

!> tear-down test features
subroutine teardown()
  implicit none
  center_sol = 0.d0;
end subroutine teardown

!> Tests ------------------------------------------------------

!> Tools ------------------------------------------------------
!>-------------------------------------------------------------
end module mod_lens_test
