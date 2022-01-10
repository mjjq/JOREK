!> the mod_light_vertices_test module contains variables
!> and methods for testing all procedures contained in
!> mod_light_vertices.f90. Due to the fact that light
!> vertices is an abstract class, the synchrotron vertex
!> type is used instead.
module mod_light_vertices_test
use fruit
use mod_vertices, only: n_x
use mod_synchrotron_light_vertices, only: synchrotron_light_vertices
implicit none

private
public :: run_fruit_light_vertices

!> Variables ---------------------------------------------------------
!> Interfaces --------------------------------------------------------
contains
!> Fruit basket ------------------------------------------------------
!> fruit basket containing all set-up, test and tear-down methods
subroutine run_fruit_light_vertices()
  implicit none
  write(*,*) "  ... setting-up: light vertices tests"
  write(*,*) "  ... running: light vertices tests"
  write(*,*) "  ... tearing-down: light vertices tests"
end subroutine run_fruit_light_vertices

!> Set-up and teard-down ---------------------------------------------
!> Tests -------------------------------------------------------------
!> Tools -------------------------------------------------------------
!>--------------------------------------------------------------------
end module mod_light_vertices_test
