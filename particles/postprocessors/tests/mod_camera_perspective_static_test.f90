!> mod_camera_perspective_static_test contains all 
!> variables and procedures used for unit testing 
!> the camera perspective static model
module mod_camera_perspective_static_test
use fruit
use mod_assert_equals_tools, only: assert_equals_allocatable_arrays
use mod_pinhole_lens, only: pinhole_lens
implicit none

private
public :: run_fruit_camera_perspective_static

!> Variable and data types -------------------------
integer,parameter :: n_properties=1
integer,parameter :: n_x_sol=3
integer,parameter :: n_points_on_lens_sol=2345
real*8,parameter  :: tol_real8=5.d-16 
real*8,dimension(3),parameter :: center_pos_lowbnd=(/-2.d-1,4.d1,-7.d0/)
real*8,dimension(3),parameter :: center_pos_uppbnd=(/3.4d1,3.d2,5.d0/)
real*8,dimension(:,:,:),allocatable :: points_on_lens
real*8,dimension(:,:,:),allocatable :: pdf_points_on_lens
type(pinhole_lens) :: pinhole
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
  call pinhole%init_pinhole(n_x_sol,center)
end subroutine setup

!> tearing-down unit test features
subroutine teardown()
  implicit none
  call pinhole%deallocate_lens;
  if(allocated(points_on_lens)) deallocate(points_on_lens)
  if(allocated(pdf_points_on_lens)) deallocate(pdf_points_on_lens)
end subroutine teardown

!> Tests -------------------------------------------
!> test generation of points and pdf on lens using
!> a pinhole lens
subroutine test_points_on_lens_pdf_pinhole()
  use mod_camera_perspective_static, only: camera_perspective_static
  implicit none
  !> variables
  type(camera_perspective_static) :: camera
  !> allocate and initialise variables
  camera%n_property_vertex = n_properties
  call camera%allocate_vertices(1,1)
  allocate(points_on_lens(n_x_sol,1,1))
  allocate(pdf_points_on_lens(1,1,1))
  call pinhole%sampling(1,points_on_lens)
  call pinhole%pdf(1,points_on_lens(:,:,1),pdf_points_on_lens(:,1,1))
  !> test generation without input number of points
  call camera%generate_points_on_lens_pdf(pinhole)
  call assert_equals(camera%n_vertices,1,&
  "Error camera perspective static generate points on pinhole: n vertices not 1!")
  call assert_equals_allocatable_arrays(n_x_sol,camera%n_vertices,&
  camera%n_times,camera%x,points_on_lens,tol_real8,&
  "Error camera perspective static generate points on pinhole: x")
  call assert_equals_allocatable_arrays(n_properties,camera%n_vertices,&
  camera%n_times,camera%properties,pdf_points_on_lens,tol_real8,&
  "Error camera perspective static generate points on pinhole: pdf points on lens")
  !> test generation with input number points
  call camera%generate_points_on_lens_pdf(pinhole,n_points_on_lens_sol)
  call assert_equals(camera%n_vertices,1,&
  "Error camera perspective static generate points on pinhole: n vertices not 1!")
  call assert_equals_allocatable_arrays(n_x_sol,camera%n_vertices,&
  camera%n_times,camera%x,points_on_lens,tol_real8,&
  "Error camera perspective static generate points on pinhole: points on lens")
  call assert_equals_allocatable_arrays(n_properties,camera%n_vertices,&
  camera%n_times,camera%properties,pdf_points_on_lens,tol_real8,&
  "Error camera perspective static generate points on pinhole: pdf points on lens") 
  !> deallocate variable
  call camera%deallocate_vertices
  deallocate(points_on_lens); deallocate(pdf_points_on_lens);
end subroutine test_points_on_lens_pdf_pinhole

!> Tools -------------------------------------------
!>--------------------------------------------------
end module mod_camera_perspective_static_test

