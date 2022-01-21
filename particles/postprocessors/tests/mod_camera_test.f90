!> mod_camera_test implements all variables and procedures
!> used for testing the camera model. Camera is an abstract
!> type => the type camera_perspective_static is used instead.
module mod_camera_test
use fruit
use mod_camera_perspective_static, only: camera_perspective_static
implicit none

private
public :: run_fruit_camera

!> Variables --------------------------------------------------------
integer,parameter :: n_x_sol=3
integer,parameter :: n_properties_sol=11
integer,parameter :: n_times_sol=12
integer,parameter :: n_pixels_x_sol=256
integer,parameter :: n_pixels_y_sol=512
integer,parameter :: n_spectra_sol=3
integer,parameter :: n_vertices_sol=1245
integer,parameter :: n_times_2_sol=33
integer,parameter :: n_pixels_x_2_sol=1024
integer,parameter :: n_pixels_y_2_sol=128
integer,parameter :: n_spectra_2_sol=12
integer,parameter :: n_vertices_2_sol=234
type(camera_perspective_static) :: camera_sol
!> Interfaces -------------------------------------------------------
contains
!> Fruit basket -----------------------------------------------------
!> fruit basket runs all the set-up, test and tear-down procedures
subroutine run_fruit_camera
  implicit none
  write(*,'(/A)') "  ... setting-up: camera tests"
  call setup()
  write(*,'(/A)') "  ... running: camera tests"
  call test_de_allocate_camera
  write(*,'(/A)') "  ... tearing-down: camera tests"
  call teardown()
end subroutine run_fruit_camera

!> Set-up and tear-down ---------------------------------------------
!> set-up test features
subroutine setup()
  implicit none
  !> set n_properties
  camera_sol%n_property_vertex = n_properties_sol
end subroutine setup

!> tear-down test features
subroutine teardown()
  implicit none
  camera_sol%n_property_vertex=0
end subroutine teardown 

!> Tests ------------------------------------------------------------
!> test allocation and deallocation of the camera classes
subroutine test_de_allocate_camera()
  implicit none
  !> test allocation from allocated
  call camera_sol%allocate_camera(n_times_sol,n_vertices_sol,&
  n_spectra_sol,n_pixels_x_sol,n_pixels_y_sol)
  call assert_equals(camera_sol%n_times,n_times_sol,&
  "Error allocate camera from unallocated: n_times mismatch!")
  call assert_equals(camera_sol%n_vertices,n_vertices_sol,&
  "Error allocate camera from unallocated: n_vertices mismatch!")
  call assert_equals(camera_sol%n_pixels_spectra,(/n_spectra_sol,n_pixels_x_sol,&
  n_pixels_y_sol/),3,"Error allocate camera from unallocated: n_pixels_spectra mismatch!")
  call assert_equals(size(camera_sol%n_active_vertices),n_times_sol,&
  "Error allocate camera from unallocated: n_active_vertices size mismatch!")
  call assert_equals(size(camera_sol%times),n_times_sol,&
  "Error allocate camera from unallocated: times size mismatch!")
  call assert_equals(shape(camera_sol%x),(/n_x_sol,n_vertices_sol,n_times_sol/),3,&
  "Error allocate camera from unallocated: x shape mismatch!")
  call assert_equals(shape(camera_sol%properties),(/n_properties_sol,n_vertices_sol,&
  n_times_sol/),3,"Error allocate camera from unallocated: properties shape mismatch!")
  call assert_equals(shape(camera_sol%pixel_intensities),(/n_spectra_sol,n_pixels_x_sol,&
  n_pixels_y_sol,n_times_sol/),4,&
  "Error allocate camera from unallocated: pixel intensities shape mismatch!")
  call assert_true(all(camera_sol%times.eq.0.d0),&
  "Error allocate camera from unallocated: pixel intensities not zero!")
  call assert_true(all(camera_sol%x.eq.0.d0),&
  "Error allocate camera from unallocated: x not zero!")
  call assert_true(all(camera_sol%properties.eq.0.d0),&
  "Error allocate camera from unallocated: properties not zero!")
  call assert_true(camera_sol%exposure_time.eq.0.d0,&
  "Error allocate camera from unallocated: exposure time not zero!")
  call assert_true(all(camera_sol%pixel_intensities.eq.0.d0),&
  "Error allocate camera from unallocated: pixel intensities not zero!")
  !> test allocation from allocated
  call camera_sol%allocate_camera(n_times_2_sol,n_vertices_2_sol,&
  n_spectra_2_sol,n_pixels_x_2_sol,n_pixels_y_2_sol)
  call assert_equals(camera_sol%n_times,n_times_2_sol,&
  "Error allocate camera from allocated: n_times mismatch!")
  call assert_equals(camera_sol%n_vertices,n_vertices_2_sol,&
  "Error allocate camera from allocated: n_vertices mismatch!")
  call assert_equals(camera_sol%n_pixels_spectra,(/n_spectra_2_sol,n_pixels_x_2_sol,&
  n_pixels_y_2_sol/),3,"Error allocate camera from allocated: n_pixels_spectra mismatch!")
  call assert_equals(size(camera_sol%n_active_vertices),n_times_2_sol,&
  "Error allocate camera from allocated: n_active_vertices size mismatch!")
  call assert_equals(size(camera_sol%times),n_times_2_sol,&
  "Error allocate camera from allocated: times size mismatch!")
  call assert_equals(shape(camera_sol%x),(/n_x_sol,n_vertices_2_sol,n_times_2_sol/),3,&
  "Error allocate camera from allocated: x shape mismatch!")
  call assert_equals(shape(camera_sol%properties),(/n_properties_sol,n_vertices_2_sol,&
  n_times_2_sol/),3,"Error allocate camera from allocated: properties shape mismatch!")
  call assert_equals(shape(camera_sol%pixel_intensities),(/n_spectra_2_sol,n_pixels_x_2_sol,&
  n_pixels_y_2_sol,n_times_2_sol/),4,&
  "Error allocate camera from allocated: pixel intensities shape mismatch!")
  !> test deallocation
  call camera_sol%deallocate_camera
  call assert_equals(camera_sol%n_times,0,&
  "Error deallocate camera: n_times not zero!")
  call assert_equals(camera_sol%n_vertices,0,&
  "Error deallocate camera: n_vertices not zero!")
  call assert_equals(camera_sol%n_pixels_spectra,(/0,0,0/),3,&
  "Error deallocate camera: n_pixels_spectra not zero!")
  call assert_true(.not.allocated(camera_sol%n_active_vertices),&
  "Error deallocate camera: n_active_vertices allocated!")
  call assert_true(.not.allocated(camera_sol%times),&
  "Error deallocate camera: times size allocated!")
  call assert_true(.not.allocated(camera_sol%x),&
  "Error deallocate camera: x allocated!")
  call assert_true(.not.allocated(camera_sol%properties),&
  "Error deallocate camera: properties allocated!")
  call assert_true(.not.allocated(camera_sol%pixel_intensities),&
  "Error deallocate camera: pixel intensities allocated!") 
   call assert_true(camera_sol%exposure_time.eq.0.d0,&
  "Error deallocate camera: exposure time not zero!")
end subroutine test_de_allocate_camera

!> Tools ------------------------------------------------------------
!>-------------------------------------------------------------------
end module mod_camera_test
