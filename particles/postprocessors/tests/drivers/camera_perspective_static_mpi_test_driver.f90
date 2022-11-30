! This program is the driver of the camera_perspective_static_mpi tests
program camera_perspective_static_mpi_test_driver
use fruit
use fruit_mpi
use mod_mpi_tools,                          only: init_mpi_threads
use mod_mpi_tools,                          only: finalize_mpi_threads
use mod_camera_perspective_static_mpi_test, only: run_fruit_camera_perspective_static_mpi
  implicit none
  integer :: rank,n_tasks,ifail

  ! initialize the mpi comm_world
  call init_mpi_threads(rank,n_tasks,ifail)

  ! init fruit suite
  call init_fruit
  call fruit_init_mpi_xml(rank)

  ! run the camera perspective static mpi tests
  call run_fruit_camera_perspective_static_mpi(rank,n_tasks,ifail)

  ! write test summary and finilize test suit
  call fruit_summary_mpi(n_tasks,rank)
  call fruit_summary_mpi_xml(n_tasks,rank)
  call fruit_finalize_mpi(n_tasks,rank)

  ! close the mpi communicator
  call finalize_mpi_threads(ifail)
end program camera_perspective_static_mpi_test_driver
