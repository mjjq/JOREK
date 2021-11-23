!> mod_spectra_deterministic_test contains variables and procedures
!> for testing the deterministic spectrum integrators
module mod_spectra_deterministic_test
use fruit
implicit none

private
public :: run_fruit_spectra_deterministic_test

!> Variables --------------------------------------------------------
!> Interfaces -------------------------------------------------------

contains

!> Fruit basket -----------------------------------------------------
!> run_fruit_spectra_deterministic executes the set-up, runs the tests
!> and clean-up all test features
subroutine run_fruit_spectra_deterministic_test()
  implicit none
  write(*,'(/A)') "  ... setting-up: spectra deterministic integrator tests"   
  write(*,'(/A)') "  ... running: spectra deterministic integrator tests"
  write(*,'(/A)') "  ... tearing-down: spectra deterministic integrator tests"
end subroutine run_fruit_spectra_deterministic_test

!> Set-up and tear-down ---------------------------------------------
!> Tests ------------------------------------------------------------

end module mod_spectra_deterministic_test
