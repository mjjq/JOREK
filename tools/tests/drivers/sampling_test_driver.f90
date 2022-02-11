!> This program is the driver of the sampling tests
program sampling_test_driver
use fruit
use mod_sampling_test, only: run_fruit_sampling
  implicit none

  ! init fruit suite
  call init_fruit
  call init_fruit_xml

  ! run the sampling test basket
  call run_fruit_sampling

  ! write test summary and finilize test suit
  call fruit_summary
  call fruit_summary_xml
  call fruit_finalize

end program sampling_test_driver
