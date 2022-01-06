!> the mod_vertices_test module contains variables
!> and methods for testing all procedures contained
!> in mod_vertices.f90. Due to the fact that vertices
!> is an abstract class, the synchrotron vertex type
!> is used instead
module mod_vertices_test
use fruit
use mod_synchrotron_light_vertices, only: synchrotron_light_vertices
implicit none

private
public :: run_fruit_vertices

!> Variables ---------------------------------------------------------
type(synchrotron_light_vertices) :: sync_light_sol

!> Interfaces --------------------------------------------------------

contains
!> Fruit basket ------------------------------------------------------
!> fruit basked containing all set-up, test and teard-down methods
subroutine run_fruit_vertices()
  implicit none
  write(*,'(/A)') "  ... setting-up: vertices tests"
  write(*,'(/A)') "  ... running: vertices tests"
  write(*,'(/A)') "  ... tearing-down: vertices tests"
end subroutine run_fruit_vertices
!> Set-up and tear-down ----------------------------------------------
!> Tests -------------------------------------------------------------
!> Tools -------------------------------------------------------------
!>--------------------------------------------------------------------
end module mod_vertices_test

