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
integer,parameter :: n_plane_points_sol=3
integer,parameter :: n_pixels_x_sol=256
integer,parameter :: n_pixels_y_sol=512
integer,parameter :: n_spectra_sol=3
integer,parameter :: n_vertices_sol=1245
integer,parameter :: n_times_2_sol=33
integer,parameter :: n_plane_points_2_sol=3
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
  !> set the number of points of a plane
  camera_sol%n_plane_points = n_plane_points_sol
end subroutine setup

!> tear-down test features
subroutine teardown()
  implicit none
  camera_sol%n_property_vertex=0; camera_sol%n_plane_points=0;
end subroutine teardown 

!> Tests ------------------------------------------------------------
!> test allocation and deallocation of the camera classes
subroutine test_de_allocate_camera()
  use mod_assert_equals_tools, only: assert_equals_allocatable_arrays
  implicit none
  !> test allocation from allocated
  call camera_sol%allocate_camera(n_times_sol,n_vertices_sol,&
  n_pixels_x_sol,n_pixels_y_sol,n_spectra_sol)
  call assert_equals(camera_sol%n_times,n_times_sol,&
  "Error allocate camera from unallocated: n_times mismatch!")
  call assert_equals(camera_sol%n_vertices,n_vertices_sol,&
  "Error allocate camera from unallocated: n_vertices mismatch!")
  call assert_equals(camera_sol%n_pixels_spectra,(/n_spectra_sol,n_pixels_x_sol,&
  n_pixels_y_sol/),3,"Error allocate camera from unallocated: n_pixels_spectra mismatch!")
  call assert_equals_allocatable_arrays(n_times_sol,camera_sol%n_active_vertices,0,&
  "Error allocate camera from unallocated: n_active_vertices")  
  call assert_equals_allocatable_arrays(n_times_sol,camera_sol%times,0.d0,&
  "Error allocate camera from unallocated: times") 
  call assert_equals_allocatable_arrays(n_x_sol,n_vertices_sol,n_times_sol,&
  camera_sol%x,0.d0,"Error allocate camera from unallocated: x")
  call assert_equals_allocatable_arrays(n_properties_sol,n_vertices_sol,n_times_sol,&
  camera_sol%properties,0.d0,"Error allocate camera from unallocated: properties")
  call assert_true(camera_sol%exposure_time.eq.0.d0,&
  "Error allocate camera from unallocated: exposure_time not zero!")
  call assert_equals_allocatable_arrays(2,n_times_sol,camera_sol%plane_edge_length,&
  0.d0,"Error allocate camera from unallocated: plane_edge_length")
  call assert_equals_allocatable_arrays(n_x_sol,n_times_sol,camera_sol%image_plane_direction,&
  0.d0,"Error allocate camera from unallocated: image_plane_direction")
  call assert_equals_allocatable_arrays(n_x_sol,n_plane_points_sol,n_times_sol,&
  camera_sol%image_plane,0.d0,"Error allocate camera from unallocated: image_plane_direction")
  !> test allocation from allocated
  camera_sol%n_plane_points = n_plane_points_2_sol
  call camera_sol%allocate_camera(n_times_2_sol,n_vertices_2_sol,&
  n_pixels_x_2_sol,n_pixels_y_2_sol,n_spectra_2_sol)
  call assert_equals(camera_sol%n_times,n_times_2_sol,&
  "Error allocate camera from allocated: n_times mismatch!")
  call assert_equals(camera_sol%n_vertices,n_vertices_2_sol,&
  "Error allocate camera from allocated: n_vertices mismatch!")
  call assert_equals(camera_sol%n_pixels_spectra,(/n_spectra_2_sol,n_pixels_x_2_sol,&
  n_pixels_y_2_sol/),3,"Error allocate camera from allocated: n_pixels_spectra mismatch!")
  call assert_equals_allocatable_arrays(n_times_2_sol,camera_sol%n_active_vertices,0,&
  "Error allocate camera from allocated: n_active_vertices")  
  call assert_equals_allocatable_arrays(n_times_2_sol,camera_sol%times,0.d0,&
  "Error allocate camera from allocated: times") 
  call assert_equals_allocatable_arrays(n_x_sol,n_vertices_2_sol,n_times_2_sol,&
  camera_sol%x,0.d0,"Error allocate camera from allocated: x")
  call assert_equals_allocatable_arrays(n_properties_sol,n_vertices_2_sol,n_times_2_sol,&
  camera_sol%properties,0.d0,"Error allocate camera from allocated: properties")
  call assert_equals_allocatable_arrays(2,n_times_2_sol,camera_sol%plane_edge_length,&
  0.d0,"Error allocate camera from allocated: plane_edge_length")
  call assert_equals_allocatable_arrays(n_x_sol,n_times_2_sol,camera_sol%image_plane_direction,&
  0.d0,"Error allocate camera from allocated: image_plane_direction")
  call assert_equals_allocatable_arrays(n_x_sol,n_plane_points_2_sol,n_times_2_sol,&
  camera_sol%image_plane,0.d0,"Error allocate camera from allocated: image_plane_direction")
  !> test deallocation
  call camera_sol%deallocate_camera
  call assert_equals(camera_sol%n_times,0,&
  "Error deallocate camera: n_times not zero!")
  call assert_equals(camera_sol%n_vertices,0,&
  "Error deallocate camera: n_vertices not zero!")
  call assert_equals(camera_sol%n_plane_points,0,&
  "Error deallocate camera: n_plane_points not zero!")
  call assert_equals(camera_sol%n_pixels_spectra,(/0,0,0/),3,&
  "Error deallocate camera: n_pixels_spectra not zero!")
  call assert_equals(camera_sol%pixel_size,(/0d0,0d0/),2,&
  "Error deallocate camera: pixels_size not zero!")
  call assert_false(allocated(camera_sol%n_active_vertices),&
  "Error deallocate camera: n_active_vertices allocated!")
  call assert_false(allocated(camera_sol%times),&
  "Error deallocate camera: times size allocated!")
  call assert_false(allocated(camera_sol%x),&
  "Error deallocate camera: x allocated!")
  call assert_false(allocated(camera_sol%properties),&
  "Error deallocate camera: properties allocated!")
  call assert_true(camera_sol%exposure_time.eq.0.d0,&
  "Error deallocate camera: exposure time not zero!")
  call assert_false(allocated(camera_sol%plane_edge_length),&
  "Error deallocate camera: plane_edge_length allocated!")
  call assert_false(allocated(camera_sol%image_plane_direction),&
  "Error deallocate camera: image_plane_direction allocated!")
  call assert_false(allocated(camera_sol%image_plane),&
  "Error deallocate camera: image_plane allocated!")
end subroutine test_de_allocate_camera

!> Tools ------------------------------------------------------------
!>-------------------------------------------------------------------
end module mod_camera_test
