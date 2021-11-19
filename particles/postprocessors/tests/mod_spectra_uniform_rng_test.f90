! mod_spectra_uniform_rng_test contains all setup, teadown and test function
! for verifying the correctness of the procedure defining spectra from uniform 
! random number distributions 
module mod_spectra_uniform_rng_test
use fruit
use mod_rng
implicit none

private
public :: run_fruit_spectra_uniform_rng_test

!> Variables -----------------------------------------------------------------
integer,parameter :: n_points=1000
integer,parameter :: n_spectra=5
!>----------------------------------------------------------------------------

contains
!> Test basket ---------------------------------------------------------------
!> test basket for executing the simulation set-up, tests and tear-down
subroutine run_fruit_spectra_uniform_rng_test()
implicit none
  write(*,'(/A)') "  ... settin-up: spectra uniform rng tests"
  write(*,'(/A)') "  ... running: spectra uniform rng tests"
  call test_spectrum_rng_uniform_construction_noinit
  write(*,'(/A)') "  ... tearing-down: spectra uniform rng tests"
end subroutine run_fruit_spectra_uniform_rng_test

!> Set-up and tear-down ------------------------------------------------------
!> Set-up test variables
subroutine setup()
end subroutine setup

!> clean up all test variables
subroutine teardown()
end subroutine teardown

!> Tests ---------------------------------------------------------------------
!> test allocation, deallocation and construction of spectrum_base class
subroutine test_spectrum_rng_uniform_construction_noinit()
  use mod_spectra, only: spectrum_rng_uniform
  implicit none
  !> variables 
  type(spectrum_rng_uniform) :: spectrum

  !> try allocation
  call spectrum%allocate_spectrum(n_points,n_spectra)
  call assert_equals(spectrum%n_points,n_points,&
  "Error spectrum base allocation: wrong n_points value!")
  call assert_equals(spectrum%n_spectra,n_spectra,&
  "Error spectrum base allocation: wrong n_spectra value!")
  call assert_true(allocated(spectrum%points),&
  "Error spectrum base allocation: points not allocated!")
  call assert_equals(shape(spectrum%points),(/n_points,n_spectra/),2,&
  "Error spectrum base allocation: points array has a wrong shape")

  !> try deallocate
  call spectrum%deallocate_spectrum
  call assert_equals(spectrum%n_points,-1,&
  "Error spectrum base deallocation: failed cleaning n_points!")
  call assert_equals(spectrum%n_spectra,-1,&
  "Error spectrum base deallocation: failed cleaning n_spectra!")
  call assert_false(allocated(spectrum%points),&
  "Error spectrum base deallocation: points still allocated!")

  !> try construction
  spectrum = spectrum_rng_uniform(n_points,n_spectra) 
  call assert_equals(spectrum%n_points,n_points,&
  "Error spectrum base construction: wrong n_points value!")
  call assert_equals(spectrum%n_spectra,n_spectra,&
  "Error spectrum base construction: wrong n_spectra value!")
  call assert_true(allocated(spectrum%points),&
  "Error spectrum base construction: points not allocated!")
  call assert_equals(shape(spectrum%points),(/n_points,n_spectra/),2,&
  "Error spectrum base construction: points array has a wrong shape")
  call spectrum%deallocate_spectrum
  
end subroutine test_spectrum_rng_uniform_construction_noinit
!>----------------------------------------------------------------------------


end module mod_spectra_uniform_rng_test
