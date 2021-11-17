! The module mod_mpi_tools_test contains variables and procedures 
! used for testing the routines of the modules mod_mpi_tools
module mod_mpi_tools_init_test
use fruit
implicit none

private
public :: run_fruit_mpi_tools_init

! Variable --------------------------------------------------------
integer :: my_id,n_tasks
real*8 :: start_time
contains

! Test basket ----------------------------------------------------
! run_fruit_mpi_tools executes the mpi_tools test set-up,
! tear-down and run the tests
subroutine run_fruit_mpi_tools_init
  implicit none

  ! execute setup -> tests -> teardown
  write(*,"(/A)") "  ... setting-up: mpi_tools tests"
  call setup
  write(*,"(/A)") "  ... running: mpi_tools tests"
  call test_init_mpi_threads
  write(*,"(/A)") "  ... tearing-down: mpi_tools tests" 
  call teardown

end subroutine run_fruit_mpi_tools_init

! Set-up and tear-down -------------------------------------------

! the setup procedure, initiliase the common test variables
subroutine setup()
  implicit none

  ! variables initialise to 0
  my_id = -999; n_tasks = -999; start_time = -9.d2;
end subroutine setup

! tear-down procedure, clean-up all variables and closes
! the mpi communicator if initialised
subroutine teardown()
  implicit none
  logical :: initialised
  integer :: ierr
  call MPI_Initialized(initialised,ierr)
  if(initialised) call MPI_Finalize(ierr)
end subroutine teardown

! Tests ----------------------------------------------------------

! test_mpi_thread_init tests the mpi thread initialisation
subroutine test_init_mpi_threads()
  use mod_mpi_tools, only: init_mpi_threads
  implicit none
  ! variables
  integer :: my_id_2,n_tasks_2,ierr

  ! initialisation
  ierr = 0; my_id_2 = -999; n_tasks_2 = -999;

  ! initialise the mpi threads
  call init_mpi_threads(my_id,n_tasks,ierr)
  
  ! checks
  call assert_equals(ierr,0,"Error: init mpi threads: en error occurred")
  call assert_true(my_id.ge.0,"Error: my_id must be non-negative")
  call assert_true(n_tasks.gt.0,"Error: number of tasks must be positive")

  ! extract time
  call init_mpi_threads(my_id_2,n_tasks_2,ierr,start_time)

  ! checks
  call assert_equals(ierr,0,"Error: init mpi threads: en error occurred")
  call assert_true(my_id_2.ge.0,"Error: my_id must be non-negative")
  call assert_true(n_tasks_2.gt.0,"Error: number of tasks must be positive")
  call assert_equals(my_id,my_id_2,"Error: tasks id must be equal")
  call assert_equals(n_tasks,n_tasks_2,"Error: number of tasks must be equal")
  call assert_true(start_time.ge.0.d0,"Error: the start time must be non-negative")
  
end subroutine test_init_mpi_threads

!-----------------------------------------------------------------

end module mod_mpi_tools_init_test
