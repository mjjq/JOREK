!> mod_spectra_deterministic_test contains variables and procedures
!> for testing the deterministic spectrum integrators
module mod_spectra_deterministic_test
use fruit
implicit none

private
public :: run_fruit_spectra_deterministic_test

!> Variables --------------------------------------------------------
integer,parameter :: n_convergence=5 !< number of points for convergence
integer,parameter :: n_points=500000 !< number of points
integer,parameter :: n_spectra=2     !< number of spectra
!> n_points for convergence study
integer,dimension(n_convergence) :: n_points_conv=(/1000,10000,100000,1000000,10000000/)
real*8,dimension(2),parameter :: min_wlen=(/3.0d-6,3.0d-7/) !< minimum wavelength
real*8,dimension(2),parameter :: max_wlen=(/3.5d-6,4.0d-7/) !< maximum wavelength
real*8,dimension(2),parameter :: min_angle=(/6.0d-1,2.3d0/) !< minimum angle for integration
real*8,dimension(2),parameter :: max_angle=(/3.6d0,3.4d0/) !< maximum angle for integration
real*8,dimension(n_spectra)   :: wbin_size !> size of the wavelength integration interval

!> Interfaces -------------------------------------------------------

contains

!> Fruit basket -----------------------------------------------------
!> run_fruit_spectra_deterministic executes the set-up, runs the tests
!> and clean-up all test features
subroutine run_fruit_spectra_deterministic_test()
  implicit none
  write(*,'(/A)') "  ... setting-up: spectra deterministic integrator tests"
  call setup
  write(*,'(/A)') "  ... running: spectra deterministic integrator tests"
  write(*,'(/A)') "  ... tearing-down: spectra deterministic integrator tests"
  call teardown
end subroutine run_fruit_spectra_deterministic_test

!> Set-up and tear-down ---------------------------------------------
!> Set-up the test variables
subroutine setup()
  implicit none
  !> set spectrum integration interval
  wbin_size = (max_wlen-min_wlen)/n_points
end subroutine setup

!> Tear-down the test variables
subroutine teardown()
  implicit none
  wbin_size = 0.d0
end subroutine teardown

!> Tests ------------------------------------------------------------
!> test the allocation and deallocation of the spectrum_integrator_1st
!> without initialisation
subroutine test_deterministic_allocation_noinit()
  use mod_spectra_deterministic, only: spectrum_integrator_1st
  implicit none
  !> variables
  type(spectrum_integrator_1st) :: spectrum

  !> test allocation and deallocation
  call spectrum%allocate_spectrum(n_points,n_spectra)
  call assert_equals(spectrum%n_points,n_points,&
  "Error spectrum integration 1st allocation: n_points do not match!")
  call assert_equals(spectrum%n_spectra,n_spectra,&
  "Error spectrum integration 1st allocation: n_spectra do not match!")
  call assert_true(allocated(spectrum%points),&
  "Error spectrum integration 1st allocation: points not allocated!")
  call spectrum%deallocate_spectrum
  call assert_equals(spectrum%n_points,-1,&
  "Error spectrum integration 1st deallocation: n_points not set to default!")
  call assert_equals(spectrum%n_spectra,-1,&
  "Error spectrum integration 1st deallocation: n_spectra not set to default!")
  call assert_fail(allocated(spectrum%points),&
  "Error spectrum integration 1st allocation: points not deallocated!")

  !> test constructor
  spectrum = spectrum_integrator_1st(n_points,n_spectra)
   call assert_equals(spectrum%n_points,n_points,&
  "Error spectrum integration 1st construction: n_points do not match!")
  call assert_equals(spectrum%n_spectra,n_spectra,&
  "Error spectrum integration 1st construction: n_spectra do not match!")
  call assert_true(allocated(spectrum%points),&
  "Error spectrum integration 1st construction: points not allocated!") 
  call spectrum%deallocate_spectrum

end subroutine test_deterministic_allocation_noinit

end module mod_spectra_deterministic_test
