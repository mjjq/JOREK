!> mod_spectra_deterministic_test contains variables and procedures
!> for testing the deterministic spectrum integrators
module mod_spectra_deterministic_test
use fruit
implicit none

private
public :: run_fruit_spectra_deterministic_test

!> Variables --------------------------------------------------------
integer,parameter :: n_convergence=5      !< number of points for convergence
integer,parameter :: n_points=512437      !< number of points
integer,parameter :: n_spectra=2          !< number of spectra
real*8,parameter  :: tol_grid=3.d-16      !< tolerance for grid check
real*8,parameter  :: accuracy_order=-2.d0  !< accuracy order
real*8,parameter  :: tol_accuracy=5.d-2   !< tolerance on the accuray order
real*8,parameter  :: tol_int_error=5.d-12 !< tolerance on the minim integratio error
!> n_points for convergence study
integer,dimension(n_convergence) :: n_points_conv=(/997,10725,100000,1003757,10023947/)
real*8,dimension(2),parameter :: min_wlen=(/3.0d-6,3.0d-7/) !< minimum wavelength
real*8,dimension(2),parameter :: max_wlen=(/3.5d-6,4.0d-7/) !< maximum wavelength
real*8,dimension(2),parameter :: min_angle=(/1.75d-1,8.4d1/) !< minimum angle for integration
real*8,dimension(2),parameter :: max_angle=(/3.25d1,1.75d2/) !< maximum angle for integration
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
  call test_deterministic_allocation_noinit
  call test_deterministic_allocation_init
  call test_set_spectrum_int_2nd_properties
  call test_generate_midpoint_spectra
  call test_2nd_order_rectangle_integrator
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
!> test the allocation, deallocation and construction of the
!> spectrum_integrator_2nd without initialisation
subroutine test_deterministic_allocation_noinit()
  use mod_spectra_deterministic, only: spectrum_integrator_2nd
  implicit none
  !> variables
  type(spectrum_integrator_2nd) :: spectrum

  !> test allocation and deallocation
  call spectrum%allocate_spectrum(n_points,n_spectra)
  call assert_equals(spectrum%n_points,n_points,&
  "Error spectrum integration 2nd allocation: n_points do not match!")
  call assert_equals(spectrum%n_spectra,n_spectra,&
  "Error spectrum integration 2nd allocation: n_spectra do not match!")
  call assert_true(allocated(spectrum%points),&
  "Error spectrum integration 2nd allocation: points not allocated!")
  call assert_equals(shape(spectrum%points),(/n_points,n_spectra/),2,&
  "Error spectrum integration 2nd allocation: points size mismatch!")
  call spectrum%deallocate_spectrum
  call assert_equals(spectrum%n_points,-1,&
  "Error spectrum integration 2nd deallocation: n_points not set to default!")
  call assert_equals(spectrum%n_spectra,-1,&
  "Error spectrum integration 2nd deallocation: n_spectra not set to default!")
  call assert_false(allocated(spectrum%points),&
  "Error spectrum integration 2nd deallocation: points not deallocated!")

  !> test constructor
  spectrum = spectrum_integrator_2nd(n_points,n_spectra)
   call assert_equals(spectrum%n_points,n_points,&
  "Error spectrum integration 2nd construction: n_points do not match!")
  call assert_equals(spectrum%n_spectra,n_spectra,&
  "Error spectrum integration 2nd construction: n_spectra do not match!")
  call assert_true(allocated(spectrum%points),&
  "Error spectrum integration 2nd construction: points not allocated!") 
  call assert_equals(shape(spectrum%points),(/n_points,n_spectra/),2,&
  "Error spectrum integration 2nd construction: points size mismatch!")
  call spectrum%deallocate_spectrum

end subroutine test_deterministic_allocation_noinit

!> test allocation, deallocation and construction of the
!> spectrum_integrator_2nd with initialisation
subroutine test_deterministic_allocation_init()
  use mod_spectra_deterministic, only: spectrum_integrator_2nd
  implicit none
  !> variables
  type(spectrum_integrator_2nd) :: spectrum
  real*8,dimension(2*n_spectra) :: real8_param

  !> initialisation
  real8_param(1:n_spectra) = min_wlen
  real8_param(n_spectra+1:2*n_spectra) = max_wlen
  !> test allocation and deallocation with initialisation
  call spectrum%allocate_spectrum(n_points,n_spectra,real8_param)
  call assert_equals(spectrum%n_points,n_points,&
  "Error spectrum integration 2nd allocation init: n_points do not match!")
  call assert_equals(spectrum%n_spectra,n_spectra,&
  "Error spectrum integration 2nd allocation init: n_spectra do not match!")
  call assert_equals(spectrum%min_wlen,min_wlen,n_spectra,&
  "Error spectrum integration 2nd allocation init: min wavelengths do not match!")
  call assert_equals(spectrum%wbin_size,wbin_size,n_spectra,&
  "Error spectrum integration 2nd allocation init: wavelengths bin size do not match!")
  call assert_true(allocated(spectrum%points),&
  "Error spectrum integration 2nd allocation init: points not allocated!") 
  call assert_equals(shape(spectrum%points),(/n_points,n_spectra/),2,&
  "Error spectrum integration 2nd allocation init: points size mismatch!")
  call spectrum%deallocate_spectrum
  call assert_equals(spectrum%n_points,-1,&
  "Error spectrum integration 2nd deallocation init: n_points not set to default!")
  call assert_equals(spectrum%n_spectra,-1,&
  "Error spectrum integration 2nd deallocation init: n_spectra not set to default!")
  call assert_false(allocated(spectrum%points),&
  "Error spectrum integration 2nd deallocation init: points not deallocated!")
  call assert_false(allocated(spectrum%min_wlen),&
  "Error spectrum integration 2nd deallocation init: min_wlen not deallocated!")
  call assert_false(allocated(spectrum%wbin_size),&
  "Error spectrum integration 2nd deallocation init: wbin_size not deallocated!")

  !> test construction with initialisation
  spectrum = spectrum_integrator_2nd(n_points,n_spectra,min_wlen,max_wlen)
  call assert_equals(spectrum%n_points,n_points,&
  "Error spectrum integration 2nd construction init: n_points do not match!")
  call assert_equals(spectrum%n_spectra,n_spectra,&
  "Error spectrum integration 2nd construction init: n_spectra do not match!")
  call assert_equals(spectrum%min_wlen,min_wlen,n_spectra,&
  "Error spectrum integration 2nd construction init: min wavelengths do not match!")
  call assert_equals(spectrum%wbin_size,wbin_size,n_spectra,&
  "Error spectrum integration 2nd construction init: wavelengths bin size do not match!")
  call assert_true(allocated(spectrum%points),&
  "Error spectrum integration 2nd construction init: points not allocated!") 
  call assert_equals(shape(spectrum%points),(/n_points,n_spectra/),2,&
  "Error spectrum integration 2nd construction init: points size mismatch!")
  call spectrum%deallocate_spectrum
end subroutine test_deterministic_allocation_init

!> test set properties deterministic integrator 2nd order 
subroutine test_set_spectrum_int_2nd_properties()
  use mod_spectra_deterministic, only: spectrum_integrator_2nd
  implicit none
  !> variables
  type(spectrum_integrator_2nd) :: spectrum
  real*8,dimension(n_spectra),parameter :: min_wlen_2=(/1.d-5,5.d-9/)
  real*8,dimension(n_spectra),parameter :: max_wlen_2=(/7.4d-5,1.25d-8/)
  real*8,dimension(n_spectra)           :: wbin_size_2
  real*8,dimension(2*n_spectra)         :: min_wlen_3,max_wlen_3,wbin_size_3

  !> initialisation
  wbin_size_2 = (max_wlen_2-min_wlen_2)/n_points
  min_wlen_3(1:n_spectra) = min_wlen; min_wlen_3(n_spectra+1:2*n_spectra) = min_wlen_2;
  max_wlen_3(1:n_spectra) = max_wlen; max_wlen_3(n_spectra+1:2*n_spectra) = max_wlen_2;
  wbin_size_3 = (max_wlen_3-min_wlen_3)/n_points

  !> test settings on unallocated properties
  spectrum = spectrum_integrator_2nd(n_points,n_spectra)
  call spectrum%set_spectrum_interval(n_spectra,min_wlen,max_wlen)
  call assert_equals(spectrum%n_spectra,n_spectra,&
  "Error spectrum integration set interval not allocated: n_spectra do not match!")
  call assert_equals(spectrum%min_wlen,min_wlen,n_spectra,&
  "Error spectrum integration set interval not allocated: min_wlen mismatch!")
  call assert_equals(spectrum%wbin_size,wbin_size,n_spectra,&
  "Error spectrum integration set interval not allocated: wbin_size mismatch!")

  !> Test change of interval alrady allocated
  call spectrum%set_spectrum_interval(n_spectra,min_wlen_2,max_wlen_2)
  call assert_equals(spectrum%n_spectra,n_spectra,&
  "Error spectrum integration set interval allocated: n_spectra do not match!")
  call assert_equals(spectrum%min_wlen,min_wlen_2,n_spectra,&
  "Error spectrum integration set interval allocated: min_wlen mismatch!")
  call assert_equals(spectrum%wbin_size,wbin_size_2,n_spectra,&
  "Error spectrum integration set interval allocated: wbin_size mismatch!")

  !> Test change of interval alrady with re-allocated
  call spectrum%set_spectrum_interval(2*n_spectra,min_wlen_3,max_wlen_3)
  call assert_equals(spectrum%n_spectra,2*n_spectra,&
  "Error spectrum integration set interval reallocated: n_spectra do not match!")
  call assert_equals(spectrum%min_wlen,min_wlen_3,2*n_spectra,&
  "Error spectrum integration set interval reallocated: min_wlen mismatch!")
  call assert_equals(spectrum%wbin_size,wbin_size_3,2*n_spectra,&
  "Error spectrum integration set interval reallocated: wbin_size mismatch!")

  !> cleanup
  call spectrum%deallocate_spectrum
end subroutine test_set_spectrum_int_2nd_properties

!> test the generation of midpoint grids
subroutine test_generate_midpoint_spectra()
  use mod_spectra_deterministic, only: spectrum_integrator_2nd
  implicit none
  !> variables
  type(spectrum_integrator_2nd) :: spectrum
  integer :: ii,jj
  real*8,dimension(n_points+1) :: interval_nodes,grid_nodes

  !> initialise variables
  spectrum = spectrum_integrator_2nd(n_points,n_spectra,min_wlen,max_wlen)
  !> generate grids and check correctness
  call spectrum%generate_spectrum
  do jj=1,spectrum%n_spectra
    do ii=1,spectrum%n_points+1
      interval_nodes(ii) = min_wlen(jj) +  wbin_size(jj)*real(ii-1,kind=8)
    enddo
    grid_nodes(1:n_points) = (spectrum%points(:,jj)-5.d-1*spectrum%wbin_size(jj))
    grid_nodes(n_points+1) = (spectrum%points(n_points,jj)+5.d-1*spectrum%wbin_size(jj))
    grid_nodes = grid_nodes/interval_nodes
    interval_nodes = 1.d0
    call assert_equals(grid_nodes,interval_nodes,n_points,tol_grid,&
    "Error spectrum integration generate spectrum: spectral grid mismatch!")
  enddo
  !> cleanaup
  call spectrum%deallocate_spectrum

end subroutine test_generate_midpoint_spectra

!> test the 2nd order rectangle method integrator. The reltive error 
!> is used due to the large value of the integral.
subroutine test_2nd_order_rectangle_integrator()
  use omp_lib
  use constants,                 only: PI
  use mod_test_functions,        only: expxsin2x,int_expxsin2x
  use mod_linear_reg,            only: linear_regression
  use mod_spectra_deterministic, only: spectrum_integrator_2nd 
  implicit none
  !> variables
  type(spectrum_integrator_2nd) :: spectrum
  integer                       :: ii,kk
  !$ integer                    :: jj
  real*8,dimension(2)           :: conv_coeff !< convergence coeff. from linear regression
  real*8,dimension(n_spectra)   :: integral
  real*8,dimension(:,:),allocatable         :: integrands
  real*8,dimension(n_spectra,n_convergence) :: rel_int_error !< integretion error

  !> convergence loop
  do kk=1,n_convergence
    !> initialise
    allocate(integrands(n_points_conv(kk),n_spectra))
    spectrum = spectrum_integrator_2nd(n_points_conv(kk),n_spectra,min_angle,max_angle)
    call spectrum%generate_spectrum
#ifdef _OPENMP
   !$omp parallel do default(private) shared(spectrum,integrands) collapse(2)
   do jj=1,spectrum%n_spectra
     do ii=1,spectrum%n_points
       integrands(ii,jj) = expxsin2x(spectrum%points(ii,jj))
     enddo
   enddo
   !$omp end parallel do
#else
    do ii=1,spectrum%n_spectra
      integrands(:,ii) = expxsin2x(spectrum%n_points,spectrum%points(:,ii))
    enddo
#endif
    call spectrum%integrate_data(integrands,integral) !< integrate values
    rel_int_error(:,kk) = abs((integral - (int_expxsin2x(n_spectra,max_angle)-&
    int_expxsin2x(n_spectra,min_angle)))/(int_expxsin2x(n_spectra,max_angle)-&
    int_expxsin2x(n_spectra,min_angle))) !< compute the error  
    deallocate(integrands)
    call spectrum%deallocate_spectrum
  enddo

  !> compute convergence rate and check for error
  do ii=1,n_spectra
    call linear_regression(n_convergence,log10(real(n_points_conv,kind=8)),&
    log10(rel_int_error(ii,:)),conv_coeff)
    call assert_equals(conv_coeff(1),accuracy_order,tol_accuracy,&
    "Error spectrum integration: expected accuracy order not matched!")
    call assert_true(rel_int_error(ii,n_convergence).lt.tol_int_error,&
    "Error spectrum integration: expected minimum error not achieved!")
  enddo
end subroutine test_2nd_order_rectangle_integrator


!>-------------------------------------------------------------------
end module mod_spectra_deterministic_test
