!> mod_synchrotron_light_vertices_test contains all variables and
!> procedures for testing the synchrotron light vertices
module mod_synchrotron_light_vertices_test
use fruit
use mod_particle_sim,               only: particle_sim
use mod_vertices,                   only: n_x
use mod_synchrotron_light_vertices, only: synchrotron_light_vertices
implicit none

private
public :: run_fruit_synchrotron_light_vertices

!> Variables ---------------------------------------------------------
!> Interfaces --------------------------------------------------------
contains
!> Fruit test basket -------------------------------------------------
!> fruit basket having all set-up, tests and tearing-down procedures
subroutine run_fruit_synchrotron_light_vertices()
  implicit none
  write(*,'(/A)') "  ... setting-up: synchrotron light vertices tests"
  write(*,'(/A)') "  ... running: synchrotron light vertices tests"
  write(*,'(/A)') "  ... tearing-down: synchrotron light vertices tests"
end subroutine run_fruit_synchrotron_light_vertices

!> Set-up and tear-down procedures------------------------------------
!> Tests -------------------------------------------------------------
!> Tools -------------------------------------------------------------
!>--------------------------------------------------------------------
end module mod_synchrotron_light_vertices_test
