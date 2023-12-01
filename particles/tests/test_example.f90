module 
use fruit
implicit none
private
public :: run_fruit_
contains
!> Fruit basket -----------------------------------
subroutine run_fruit_
  implicit none
  write(*,'(/A)') "  ... setting-up: "
  write(*,'(/A)') "  ... running: "
  call run_test_case(,'')
  write(*,'(/A)') "  ... tearing-down: "
end subroutine run_fruit_

!> Tests ------------------------------------------
!> ------------------------------------------------
end module 
