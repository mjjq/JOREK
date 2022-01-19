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
!> Interfaces -------------------------------------------------
contains
!> Fruit basket -----------------------------------------------
!> basket running the set-up, test and teard-down procedures
subroutine run_fruit_lens()
  implicit none
  write(*,'(/A)') "  ... setting-up: lens tests"
  write(*,'(/A)') "  ... running: lens tests"
  write(*,'(/A)') "  ... tearing-down: lens tests"
end subroutine run_fruit_lens

!> Set-up and teard-down---------------------------------------

!> Tests ------------------------------------------------------
!> Tools ------------------------------------------------------
!>-------------------------------------------------------------
end module mod_lens_test
