!> mod_camera_perspective_static_mpi_test cantains all variables
!> and procedures used for testing mpi enabled procedures of the
!> camera_perspective_static class.
module mod_camera_perspective_static_mpi_test
use fruit
use fruit_mpi
use mod_omnidirectional_gaussian_lights, only: omnidirectional_gaussian_lights
implicit none
private
public :: run_fruit_camera_perspective_static_mpi

!> Variable and data types ---------------------------------------------------
type(omnidirectional_gaussian_lights) :: lights_sol

!> Interfaces ----------------------------------------------------------------
contains
!> Fruit basket --------------------------------------------------------------
!> fruit basket executes all set-up, test and tear-down procedures
subroutine run_fruit_camera_perspective_static_mpi(rank,n_tasks,ifail)
  implicit none
  !> inputs-outputs:
  integer,intent(inout) :: ifail
  !> inputs:
  integer,intent(in) :: rank,n_tasks
  if(rank.eq.0) write(*,*) "  ... setting-up: camera perspective static mpi tests"
  if(rank.eq.0) write(*,*) "  ... running: camera perspective static mpi tests"
  if(rank.eq.0) write(*,*) "  ... tearing-down: camera perspective static mpi tests"
end subroutine run_fruit_camera_perspective_static_mpi

!> Set-up and tear-down ------------------------------------------------------
!> Tests ---------------------------------------------------------------------
!> Tools ---------------------------------------------------------------------
!>----------------------------------------------------------------------------
end module mod_camera_perspective_static_mpi_test

