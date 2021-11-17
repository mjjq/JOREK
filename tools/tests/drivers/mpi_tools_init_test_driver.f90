! This program is the driver of the boost besselk tests
program mpi_tools_init_test_driver
use fruit
use mod_mpi_tools_init_test, only: run_fruit_mpi_tools_init
  implicit none

  ! init fruit suite
  call init_fruit
  call init_fruit_xml

  ! run then mpi_tools test basket
  call run_fruit_mpi_tools_init

  ! write test summary and finilize test suit
  call fruit_summary
  call fruit_summary_xml
  call fruit_finalize

end program mpi_tools_init_test_driver
