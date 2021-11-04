! This module contains the fruit basket for running
! all besselik test cases
module mod_besselik_test_basket
use fruit
use mod_besselik_test
implicit none
contains

! run_all_basselik_tests performs the setup,
! run the tests and performs the teardown of
! the besselik functions
subroutine run_all_besselik_tests()
  implicit none

  write(*,'(/A)') "  ... set-up: besselik tests" 
	call setup !< setup test variables
  write(*,'(/A)') "  ... running: besselik tests"
	call test_split_array_value
  write(*,'(/A)') "  ... tearing-down: besselik tests" 
	call teardown !< cleanup test variables

end subroutine run_all_besselik_tests

! define the fruit basket 
subroutine fruit_basket
  implicit none
	call run_all_besselik_tests
end subroutine fruit_basket

end module mod_besselik_test_basket
