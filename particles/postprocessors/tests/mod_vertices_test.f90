!> the mod_vertices_test module contains variables
!> and methods for testing all procedures contained
!> in mod_vertices.f90. Due to the fact that vertices
!> is an abstract class, the synchrotron vertex type
!> is used instead
module mod_vertices_test
use fruit
use mod_vertices, only: n_x
use mod_synchrotron_light_vertices, only: synchrotron_light_vertices
implicit none

private
public :: run_fruit_vertices

!> Variables ---------------------------------------------------------
integer,parameter :: n_times_sol=3
integer,parameter :: n_times_2_sol=5
integer,parameter :: n_vertices_sol=111
integer,parameter :: n_vertices_2_sol=222
integer,parameter :: n_properties_sol=11
type(synchrotron_light_vertices) :: vertex_sol

!> Interfaces --------------------------------------------------------

contains
!> Fruit basket ------------------------------------------------------
!> fruit basked containing all set-up, test and teard-down methods
subroutine run_fruit_vertices()
  implicit none
  write(*,'(/A)') "  ... setting-up: vertices tests"
  call setup
  write(*,'(/A)') "  ... running: vertices tests"
  call test_de_allocate_time_vector 
  call test_de_allocates_vertex_x_properties
  call test_de_allocates_vertices
  write(*,'(/A)') "  ... tearing-down: vertices tests"
  call teardown
end subroutine run_fruit_vertices

!> Set-up and tear-down ----------------------------------------------
!> set-up unit test features
subroutine setup()
  implicit none
  vertex_sol%n_property_vertex = n_properties_sol !< set "fake" number of properties
end subroutine setup

!> tear-down all unit test features
subroutine teardown()
  implicit none
  vertex_sol%n_property_vertex = 0 
end subroutine teardown

!> Tests -------------------------------------------------------------
!> test allocate and deallocate time vector
subroutine test_de_allocate_time_vector()
  implicit none
  !> test allocate time vector from unallocated from allocated vector
  call vertex_sol%allocate_time_vector(n_times_sol)
  call assert_equals(vertex_sol%n_times,n_times_sol,&
  "Error allocate vertices time vector from unallocated: n_times mismatch!")
  call assert_true(allocated(vertex_sol%n_active_vertices),&
  "Error allocate vertices time vector from unallocated: n_active_vertices not allocated!")
  call assert_true(allocated(vertex_sol%times),&
  "Error allocate vertices time vector from unallocated: times not allocated!")
  call assert_equals(size(vertex_sol%n_active_vertices),n_times_sol,&
  "Error allocate vertices time vector from unallocated: n_active_vertices size mismatch!")
  call assert_equals(size(vertex_sol%times),n_times_sol,&
  "Error allocate vertices time vector from unallocated: times size mismatch!")
  call assert_true(all(vertex_sol%n_active_vertices.eq.0),&
  "Error allocate vertices time vector from unallocated: n_active_vertices not zero!")
  call assert_true(all(vertex_sol%times.eq.0.d0),&
  "Error allocate vertices time vector from unallocated: times not zero!")
  !> test allocated time vector from allocated vector
  call vertex_sol%allocate_time_vector(n_times_2_sol) 
  call assert_equals(vertex_sol%n_times,n_times_2_sol,&
  "Error allocate vertices time vector from allocated: n_times mismatch!")
  call assert_true(allocated(vertex_sol%n_active_vertices),&
  "Error allocate vertices time vector from allocated: n_active_vertices not allocated!")
  call assert_true(allocated(vertex_sol%times),&
  "Error allocate vertices time vector from allocated: times not allocated!")
  call assert_equals(size(vertex_sol%n_active_vertices),n_times_2_sol,&
  "Error allocate vertices time vector from allocated: n_active_vertices size mismatch!")
  call assert_equals(size(vertex_sol%times),n_times_2_sol,&
  "Error allocate vertices time vector from allocated: times size mismatch!")
  call assert_true(all(vertex_sol%n_active_vertices.eq.0),&
  "Error allocate vertices time vector from allocated: n_active_vertices not zero!")
  call assert_true(all(vertex_sol%times.eq.0.d0),&
  "Error allocate vertices time vector from allocated: times not zero!")
  !> test deallocate time vector
  call vertex_sol%deallocate_time_vector
  call assert_equals(vertex_sol%n_times,0,&
  "Error deallocate vertices time vector: n_times not reset!")
  call assert_true(.not.allocated(vertex_sol%n_active_vertices),&
  "Error deallocate vertices time vector: n_active_vertices allocated!")
  call assert_true(.not.allocated(vertex_sol%times),&
  "Error deallocate vertices time vector: times allocated!")
end subroutine test_de_allocate_time_vector

!> test the allocation / deallocation of vertex position and properties
subroutine test_de_allocates_vertex_x_properties()
  implicit none
  !> initialise the time structure first
  call vertex_sol%allocate_time_vector(n_times_sol) 
  !> allocate x and properties from unallocated vertex
  call vertex_sol%allocate_x_properties(n_vertices_sol)
  call assert_equals(vertex_sol%n_vertices,n_vertices_sol,&
  "Error allocate vertices x-property vector from unallocated: n_vertices mismatch!")
  call assert_true(allocated(vertex_sol%x),&
  "Error allocate vertices  x-property vector from unallocated: x not allocated!")
  call assert_true(allocated(vertex_sol%properties),&
  "Error allocate vertices x-property vector from unallocated: properties not allocated!")
  call assert_equals(shape(vertex_sol%x),(/n_x,n_vertices_sol,n_times_sol/),3,&
  "Error allocate vertices x-property vectior from unallocated: x size mismatch!")
  call assert_equals(shape(vertex_sol%properties),&
  (/n_properties_sol,n_vertices_sol,n_times_sol/),3,&
  "Error allocate vertices x-property vector from unallocated: properties size mismatch!")
  call assert_true(all(vertex_sol%x.eq.0),&
  "Error allocate vertices x-property vector from unallocated: x not zero!")
  call assert_true(all(vertex_sol%properties.eq.0.d0),&
  "Error allocate vertices x-property vector from unallocated: properties not zero!")
  !> allocate x and properties from allocated vertes
  call vertex_sol%allocate_x_properties(n_vertices_2_sol)
  call assert_equals(vertex_sol%n_vertices,n_vertices_2_sol,&
  "Error allocate vertices x-property vector from allocated: n_vertices mismatch!")
  call assert_true(allocated(vertex_sol%x),&
  "Error allocate vertices  x-property vector from allocated: x not allocated!")
  call assert_true(allocated(vertex_sol%properties),&
  "Error allocate vertices x-property vector from allocated: properties not allocated!")
  call assert_equals(shape(vertex_sol%x),(/n_x,n_vertices_2_sol,n_times_sol/),3,&
  "Error allocate vertices x-property vector from allocated: x size mismatch!")
  call assert_equals(shape(vertex_sol%properties),&
  (/n_properties_sol,n_vertices_2_sol,n_times_sol/),3,&
  "Error allocate vertices x-property vector from allocated: properties size mismatch!")
  call assert_true(all(vertex_sol%x.eq.0),&
  "Error allocate vertices x-property vector from allocated: x not zero!")
  call assert_true(all(vertex_sol%properties.eq.0.d0),&
  "Error allocate vertices x-property vector from allocated: properties not zero!")
  !> deallocate x and properties
  call vertex_sol%deallocate_x_properties
  call assert_equals(vertex_sol%n_vertices,0,&
  "Error deallocate vertices x-property vector: n_vertices not reset!")
  call assert_true(.not.allocated(vertex_sol%x),&
  "Error deallocate vertices x-property vector: x allocated!")
  call assert_true(.not.allocated(vertex_sol%properties),&
  "Error deallocate vertices x-property vector: properties allocated!") 
end subroutine test_de_allocates_vertex_x_properties

!> check allocation and deallocation of vertices
subroutine test_de_allocates_vertices()
  implicit none

  !> test vertex allocation from unallocated
  call vertex_sol%allocate_vertices(n_times_sol,n_vertices_sol)
  call assert_equals(vertex_sol%n_times,n_times_sol,&
  "Error allocate vertices from unallocated: n_times mismatch!")
  call assert_equals(vertex_sol%n_vertices,n_vertices_sol,&
  "Error allocate vertices from unallocated: n_vertices mismatch!")
  call assert_equals(size(vertex_sol%n_active_vertices),n_times_sol,&
  "Error allocate vertices from unallocated: n_active_vertices size mismatch!")
  call assert_equals(size(vertex_sol%times),n_times_sol,&
  "Error allocate vertices from unallocated: times size mismatch!")
  call assert_equals(shape(vertex_sol%x),(/n_x,n_vertices_sol,n_times_sol/),3,&
  "Error allocate vertices from unallocated: x size mismatch!")
  call assert_equals(shape(vertex_sol%properties),(/n_properties_sol,n_vertices_sol,n_times_sol/),3,&
  "Error allocate vertices from unallocated: properties size mismatch!")
  !> test vertex allocation from allocated
  call vertex_sol%allocate_vertices(n_times_2_sol,n_vertices_2_sol)
  call assert_equals(vertex_sol%n_times,n_times_2_sol,&
  "Error allocate vertices from allocated: n_times mismatch!")
  call assert_equals(vertex_sol%n_vertices,n_vertices_2_sol,&
  "Error allocate vertices from allocated: n_vertices mismatch!")
  call assert_equals(size(vertex_sol%n_active_vertices),n_times_2_sol,&
  "Error allocate vertices from allocated: n_active_vertices size mismatch!")
  call assert_equals(size(vertex_sol%times),n_times_2_sol,&
  "Error allocate vertices from allocated: times size mismatch!")
  call assert_equals(shape(vertex_sol%x),(/n_x,n_vertices_2_sol,n_times_2_sol/),3,&
  "Error allocate vertices from allocated: x size mismatch!")
  call assert_equals(shape(vertex_sol%properties),&
  (/n_properties_sol,n_vertices_2_sol,n_times_2_sol/),3,&
  "Error allocate vertices from allocated: properties size mismatch!")
  !> test vertex deallocation
  call vertex_sol%deallocate_vertices
  call assert_equals(vertex_sol%n_times,0,&
  "Error deallocate vertices: n_times not reset!")
  call assert_equals(vertex_sol%n_vertices,0,&
  "Error deallocate vertices: n_vertices not reset!")
  call assert_true(.not.allocated(vertex_sol%n_active_vertices),&
  "Error deallocate vertices: n_active_vertices allocated!")
  call assert_true(.not.allocated(vertex_sol%times),&
  "Error deallocate vertices: times allocated!")
  call assert_true(.not.allocated(vertex_sol%x),&
  "Error deallocate vertices: x allocated!")
  call assert_true(.not.allocated(vertex_sol%properties),&
  "Error deallocate vertices: properties allocated!")
end subroutine test_de_allocates_vertices

!> Tools -------------------------------------------------------------
!>--------------------------------------------------------------------
end module mod_vertices_test

