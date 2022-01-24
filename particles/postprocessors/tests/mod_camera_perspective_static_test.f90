!> mod_camera_perspective_static_test contains all 
!> variables and procedures used for unit testing 
!> the camera perspective static model
module mod_camera_perspective_static_test
use fruit
use mod_assert_equals_tools,   only: assert_equals_allocatable_arrays
use mod_pinhole_lens,          only: pinhole_lens
use mod_spectra_deterministic, only: spectrum_integrator_2nd
use mod_camera_perspective_static, only: camera_perspective_static
implicit none

private
public :: run_fruit_camera_perspective_static

!> Variable and data types -------------------------
integer,parameter :: n_properties=1
integer,parameter :: n_x_sol=3
integer,parameter :: n_times_sol=1
integer,parameter :: n_plane_vertices=3
integer,parameter :: n_pixels_x=512
integer,parameter :: n_pixels_y=512
integer,parameter :: n_points_on_lens_sol=2345
integer,parameter :: n_lines_per_spectrum=13
integer,parameter :: n_spectra=2
real*8,parameter  :: tol_real8=5.d-16 
real*8,dimension(3),parameter :: center_pos_lowbnd=(/-2.d-1,4.d1,-7.d0/)
real*8,dimension(3),parameter :: center_pos_uppbnd=(/3.4d1,3.d2,5.d0/)
real*8,dimension(n_spectra),parameter :: min_wlen=(/3.d0-6,2.5d-7/)
real*8,dimension(n_spectra),parameter :: max_wlen=(/3.5d0-6,4.2d-7/)
real*8,dimension(:,:,:),allocatable :: points_on_lens
real*8,dimension(:,:,:),allocatable :: pdf_points_on_lens
type(camera_perspective_static) :: camera_sol
type(pinhole_lens)              :: pinhole_sol
type(spectrum_integrator_2nd)   :: spectrum_sol
!> Interfaces --------------------------------------
contains
!> Fruit basket ------------------------------------
!> fruit basket executes all set-up, test and 
!> tear-down procedures
subroutine run_fruit_camera_perspective_static
  implicit none
  write(*,'(/A)') "  ... setting-up: camera perspective static tests"
  call setup
  write(*,'(/A)') "  ... running: camera perspective static tests"
  call test_points_on_lens_pdf_pinhole
  call test_de_allocation_camera_perspective_static
  write(*,'(/A)') "  ... tearing-up: camera perspective static tests"
  call teardown
end subroutine run_fruit_camera_perspective_static

!> Set-up and tear-down ----------------------------
!> set-up unit test features
subroutine setup()
  use mod_gnu_rng, only: gnu_rng_interval
  implicit none
  real*8,dimension(3) :: center
  !> initialise the pinhole lens
  call gnu_rng_interval(3,center_pos_lowbnd,center_pos_uppbnd,center)
  call pinhole_sol%init_pinhole(n_x_sol,center)
  !> initialise camera spectrum
  spectrum_sol = spectrum_integrator_2nd(n_lines_per_spectrum,&
  n_spectra,min_wlen,max_wlen)
  !> initialise camera
  camera_sol%n_plane_points = n_plane_vertices
end subroutine setup

!> tearing-down unit test features
subroutine teardown()
  implicit none
  call pinhole_sol%deallocate_lens; call spectrum_sol%deallocate_spectrum;
  if(allocated(points_on_lens)) deallocate(points_on_lens)
  if(allocated(pdf_points_on_lens)) deallocate(pdf_points_on_lens)
end subroutine teardown

!> Tests -------------------------------------------
!> test allocation and deallocation of camera_perspective_static
!> attributs (only)
subroutine test_de_allocation_camera_perspective_static()
  use mod_assert_equals_tools,       only: assert_equals_allocatable_arrays
  implicit none
  !> allocate camera perspective static
  call camera_sol%allocate_camera_perspective_static(spectrum_sol,3,&
  (/n_points_on_lens_sol,n_pixels_x,n_pixels_y/))
  !> test allocation
  call assert_equals(camera_sol%n_vertices,n_points_on_lens_sol,&
  "Error allocate camera perspective static: n vertices mismatch")
  call assert_equals_allocatable_arrays(n_times_sol,camera_sol%n_active_vertices,&
  "Error allocate camera perspective static: n active vertices")
  call assert_equals_allocatable_arrays(n_times_sol,camera_sol%times,&
  "Error allocate camera perspective static: times")
  call assert_equals_allocatable_arrays(n_x_sol,n_points_on_lens_sol,&
  n_times_sol,camera_sol%x,"Error allocate camera perspective static: x")
  call assert_equals_allocatable_arrays(n_properties,n_points_on_lens_sol,n_times_sol,&
  camera_sol%properties,"Error allocate camera perspective static: properties")
  call assert_equals(camera_sol%n_pixels_spectra,(/n_spectra,n_pixels_x,n_pixels_y/),&
  3,"Error allocate camera perspective static: n pixels / spectra")
  call assert_equals_allocatable_arrays(spectrum_sol%n_spectra,n_pixels_x,&
  n_pixels_y,n_times_sol,camera_sol%pixel_intensities,&
  "Error allocate camera perspective static: pixel intensities")
  call assert_equals_allocatable_arrays(n_x_sol,n_plane_vertices,&
  camera_sol%image_plane,"Error allocate_camera perspective static: image plane") 
  !> deallocate camera perspective static
  call camera_sol%deallocate_camera_perspective_static
  call assert_equals(camera_sol%n_vertices,0,&
  "Error deallocate camera perspective static: n vertices not 0")
  call assert_equals(camera_sol%n_pixels_spectra,(/0,0,0/),3,&
  "Error deallocate camera perspective static: n pixels spectra not 0")
  call assert_false(allocated(camera_sol%n_active_vertices),&
  "Error deallocate camera perspective static: n active vertices not deallocated")
  call assert_false(allocated(camera_sol%times),&
  "Error deallocate camera perspective static: times not deallocated")
  call assert_false(allocated(camera_sol%x),&
  "Error deallocate camera perspective static: x not deallocated")
  call assert_false(allocated(camera_sol%properties),&
  "Error deallocate camera perspective static: properties not deallocated")
  call assert_false(allocated(camera_sol%pixel_intensities),&
  "Error deallocate camera perspective static: pixel intensities not deallocated")
  call assert_false(allocated(camera_sol%image_plane),&
  "Error deallocate camera perspective static: image plane not deallocated")
end subroutine test_de_allocation_camera_perspective_static

!> test generation of points and pdf on lens using
!> a pinhole lens
subroutine test_points_on_lens_pdf_pinhole()
  implicit none
  !> allocate and initialise variables
  camera_sol%n_property_vertex = n_properties
  call camera_sol%allocate_vertices(1,1)
  allocate(points_on_lens(n_x_sol,1,1))
  allocate(pdf_points_on_lens(1,1,1))
  call pinhole_sol%sampling(1,points_on_lens)
  call pinhole_sol%pdf(1,points_on_lens(:,:,1),pdf_points_on_lens(:,1,1))
  !> test generation without input number of points
  call camera_sol%generate_points_on_lens_pdf(pinhole_sol)
  call assert_equals(camera_sol%n_vertices,1,&
  "Error camera perspective static generate points on pinhole: n vertices not 1!")
  call assert_equals_allocatable_arrays(n_x_sol,camera_sol%n_vertices,&
  camera_sol%n_times,camera_sol%x,points_on_lens,tol_real8,&
  "Error camera perspective static generate points on pinhole: x")
  call assert_equals_allocatable_arrays(n_properties,camera_sol%n_vertices,&
  camera_sol%n_times,camera_sol%properties,pdf_points_on_lens,tol_real8,&
  "Error camera perspective static generate points on pinhole: pdf points on lens")
  !> test generation with input number points
  call camera_sol%generate_points_on_lens_pdf(pinhole_sol,n_points_on_lens_sol)
  call assert_equals(camera_sol%n_vertices,1,&
  "Error camera perspective static generate points on pinhole: n vertices not 1!")
  call assert_equals_allocatable_arrays(n_x_sol,camera_sol%n_vertices,&
  camera_sol%n_times,camera_sol%x,points_on_lens,tol_real8,&
  "Error camera perspective static generate points on pinhole: points on lens")
  call assert_equals_allocatable_arrays(n_properties,camera_sol%n_vertices,&
  camera_sol%n_times,camera_sol%properties,pdf_points_on_lens,tol_real8,&
  "Error camera perspective static generate points on pinhole: pdf points on lens") 
  !> deallocate variable
  call camera_sol%deallocate_vertices
  deallocate(points_on_lens); deallocate(pdf_points_on_lens);
end subroutine test_points_on_lens_pdf_pinhole

!> Tools -------------------------------------------
!>--------------------------------------------------
end module mod_camera_perspective_static_test

