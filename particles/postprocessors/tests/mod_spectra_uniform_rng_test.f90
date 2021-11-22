! mod_spectra_uniform_rng_test contains all setup, teadown and test function
! for verifying the correctness of the procedure defining spectra from uniform 
! random number distributions 
module mod_spectra_uniform_rng_test
use fruit
use mod_rng, only: type_rng
implicit none

private
public :: run_fruit_spectra_uniform_rng_test

!> Variables -----------------------------------------------------------------
integer,parameter :: n_points=1000000
integer,parameter :: n_spectra=2
real*8,parameter  :: tol_real8=1.d-16 !< tolerance for assert
real*8,dimension(2),parameter :: min_wlen=(/3.d-6,3.0d-7/) !< minimum wavelength
real*8,dimension(2),parameter :: max_wlen=(/3.5d-6,4.d-7/) !< maximum wavelength
class(type_rng),dimension(:),allocatable :: rngs !< random number generators
integer                                    :: n_threads,thread_id !< N# and id omp threads
real*8,dimension(2)                        :: i_pdf !< 1/pdf=min_wlen,max_wlen
!>----------------------------------------------------------------------------

contains
!> Test basket ---------------------------------------------------------------
!> test basket for executing the simulation set-up, tests and tear-down
subroutine run_fruit_spectra_uniform_rng_test()
implicit none
  write(*,'(/A)') "  ... settin-up: spectra uniform rng tests"
  call setup
  write(*,'(/A)') "  ... running: spectra uniform rng tests"
  call test_spectrum_rng_uniform_construction_noinit
  call test_set_uniform_spectrum_interval
  call test_spectrum_rng_uniform_construction_init
  call test_spectrum_generation_rng_uniform
  write(*,'(/A)') "  ... tearing-down: spectra uniform rng tests"
  call teardown
end subroutine run_fruit_spectra_uniform_rng_test

!> Set-up and tear-down ------------------------------------------------------
!> Set-up test variables
subroutine setup()
  use mod_random_seed,only: random_seed
  use mod_pcg32_rng,  only: pcg32_rng
  !$ use omp_lib
  implicit none
  !> variables
  integer :: n_streams,thread_id,ifail

  !> set the inverse of the pdf
  i_pdf = max_wlen-min_wlen
  !> initialise the rngs using the pcg32
  n_threads = 1
  thread_id = 1
  !$ n_threads = omp_get_max_threads()
  allocate(pcg32_rng::rngs(n_threads))
  !$omp parallel default(private) shared(rngs,n_threads,ifail)
  !$ thread_id = omp_get_thread_num()+1
  call rngs(thread_id)%initialize(n_dims=n_spectra,seed=random_seed(),&
  n_streams=n_threads*n_spectra,i_stream=thread_id,ierr=ifail)
  !$omp end parallel 
end subroutine setup

!> clean up all test variables
subroutine teardown()
  implicit none
  i_pdf = 0.d0
  deallocate(rngs)
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
  "Error spectrum base allocation: n_points mismatch!")
  call assert_equals(spectrum%n_spectra,n_spectra,&
  "Error spectrum base allocation: n_spectra mismatch!")
  call assert_true(allocated(spectrum%points),&
  "Error spectrum base allocation: points not allocated!")
  call assert_equals(shape(spectrum%points),(/n_points,n_spectra/),2,&
  "Error spectrum base allocation: points array shape mismatch!")

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
  "Error spectrum base construction: n_points mismatch!")
  call assert_equals(spectrum%n_spectra,n_spectra,&
  "Error spectrum base construction: n_spectra mismatch!")
  call assert_true(allocated(spectrum%points),&
  "Error spectrum base construction: points not allocated!")
  call assert_equals(shape(spectrum%points),(/n_points,n_spectra/),2,&
  "Error spectrum base construction: points array shape mismatch!")
  call spectrum%deallocate_spectrum
  
end subroutine test_spectrum_rng_uniform_construction_noinit

!> test set and change of the spectrum interval
subroutine test_set_uniform_spectrum_interval()
  use mod_spectra, only: spectrum_rng_uniform
  implicit none
  !> variables
  type(spectrum_rng_uniform) :: spectrum
  real*8,dimension(n_spectra),parameter :: min_wlen_2=(/1.d-5,5.d-9/)
  real*8,dimension(n_spectra),parameter :: max_wlen_2=(7.4d-5,1.25d-8)
  real*8,dimension(2*n_spectra)         :: min_wlen_3,max_wlen_3

  !> test setting of min_wlen,max_wlen
  spectrum = spectrum_rng_uniform(n_points,n_spectra)
  call spectrum%set_spectrum_interval(n_spectra,min_wlen,max_wlen)
  call assert_equals(spectrum%n_spectra,n_spectra,"Error set uniform interval: n_spectra mismatch!")
  call assert_equals(spectrum%min_wlen,min_wlen,n_spectra,tol_real8,&
  "Error set uniform interval: min_wlen mismatch")
  call assert_equals(spectrum%i_pdf,i_pdf,n_spectra,tol_real8,&
  "Error set uniform interval: i_pdf mismatch")

  !> test change with equal number of spectral intervals
  call spectrum%set_spectrum_interval(n_spectra,min_wlen_2,max_wlen_2)
  call assert_equals(spectrum%min_wlen,min_wlen_2,n_spectra,tol_real8,&
  "Error set uniform interval: min_wlen_2 mismatch")
  call assert_equals(spectrum%i_pdf,max_wlen_2-min_wlen_2,n_spectra,tol_real8,&
  "Error set uniform interval: i_pdf_2 mismatch")

  !> test change with different number of intervals
  min_wlen_3(1:n_spectra) = min_wlen; min_wlen_3(n_spectra+1:2*n_spectra) = min_wlen_2;
  max_wlen_3(1:n_spectra) = max_wlen; max_wlen_3(n_spectra+1:2*n_spectra) = max_wlen_2;
  call spectrum%set_spectrum_interval(2*n_spectra,min_wlen_3,max_wlen_3)
  call assert_equals(spectrum%n_spectra,2*n_spectra,&
  "Error set uniform interval: n_spectra_2 mismatch!")
  call assert_equals(spectrum%min_wlen,min_wlen_3,2*n_spectra,tol_real8,&
  "Error set uniform interval: min_wlen_3 mismatch")
  call assert_equals(spectrum%i_pdf,max_wlen_3-min_wlen_3,2*n_spectra,tol_real8,&
  "Error set uniform interval: i_pdf_3 mismatch")  

  !> deallocate variables
  call spectrum%deallocate_spectrum()
end subroutine test_set_uniform_spectrum_interval

!> test allocation and construction with initialisation
subroutine test_spectrum_rng_uniform_construction_init()
  use mod_spectra, only: spectrum_rng_uniform
  implicit none
  !> variables
  type(spectrum_rng_uniform) :: spectrum
  real*8,dimension(2*n_spectra) :: real8_param

  !> test allocation with initialisation
  real8_param(1:n_spectra) = min_wlen
  real8_param(n_spectra+1:2*n_spectra) = max_wlen
  spectrum = spectrum_rng_uniform(n_points,n_spectra)
  call spectrum%allocate_spectrum(n_points,n_spectra,real8_param)
  call assert_equals(spectrum%n_points,n_points,&
  "Error spectrum rng uniform allocation: n_points mismatch!")
  call assert_equals(spectrum%n_spectra,n_spectra,&
  "Error spectrum rng uniform allocation: n_spectra mismatch!")
  call assert_true(allocated(spectrum%points),&
  "Error spectrum rng uniform allocation: points not allocated!")
  call assert_equals(shape(spectrum%points),(/n_points,n_spectra/),2,&
  "Error spectrum rng uniform allocation: points array shape mismatch!")
  call assert_equals(spectrum%min_wlen,min_wlen,n_spectra,tol_real8,&
  "Error spectrum rng uniform allocation: min_wlen mismatch")
  call assert_equals(spectrum%i_pdf,i_pdf,n_spectra,tol_real8,&
  "Error spectrum rng uniform allocation: i_pdf mismatch")
  call spectrum%deallocate_spectrum

  !> test constructor with initialisation
  spectrum = spectrum_rng_uniform(n_points,n_spectra,min_wlen,max_wlen)
  call assert_equals(spectrum%n_points,n_points,&
  "Error spectrum rng uniform constructor: n_points mismatch!")
  call assert_equals(spectrum%n_spectra,n_spectra,&
  "Error spectrum rng uniform constructor: n_spectra mismatch!")
  call assert_true(allocated(spectrum%points),&
  "Error spectrum rng uniform constructor: points not allocated!")
  call assert_equals(shape(spectrum%points),(/n_points,n_spectra/),2,&
  "Error spectrum rng uniform constructor: points array shape mismatch!")
  call assert_equals(spectrum%min_wlen,min_wlen,n_spectra,tol_real8,&
  "Error spectrum rng uniform constructor: min_wlen mismatch")
  call assert_equals(spectrum%i_pdf,i_pdf,n_spectra,tol_real8,&
  "Error spectrum rng uniform constructor: i_pdf mismatch")
  call spectrum%deallocate_spectrum

end subroutine test_spectrum_rng_uniform_construction_init

!> test the random generation of uniform spectra within interval
subroutine test_spectrum_generation_rng_uniform()
  use mod_spectra, only: spectrum_rng_uniform
  implicit none
  !> variables
  type(spectrum_rng_uniform) :: spectrum
  integer :: ii

  !> construct the spectrum class
  spectrum = spectrum_rng_uniform(n_points,n_spectra,min_wlen,max_wlen)
  !> generate new random spectra, check intervals bounds
  call spectrum%generate_spectrum(rngs)
  do ii=1,n_spectra
    call assert_true(minval(spectrum%points(:,ii)).ge.min_wlen(ii),&
    "Error generate uniform random spectrum: minmum point value out-of-bound")
    call assert_true(maxval(spectrum%points(:,ii)).le.max_wlen(ii),&
    "Error generate uniform random spectrum: maximum point value out-of-bound")
  enddo
  !> clean-up the spectrum class
  call spectrum%deallocate_spectrum

end subroutine test_spectrum_generation_rng_uniform

!> test the integration via Monte-Carlo method (uniform distribution)

!>----------------------------------------------------------------------------


end module mod_spectra_uniform_rng_test
