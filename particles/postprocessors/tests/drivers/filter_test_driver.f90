! This program is the driver of the filter tests
program filter_test
use fruit
use mod_filter_test, only: run_fruit_filter
  implicit none

  ! init fruit suite
  call init_fruit
  call init_fruit_xml

  ! run the filter test basket
  call run_fruit_filter

  ! write test summary and finilize test suit
  call fruit_summary
  call fruit_summary_xml
  call fruit_finalize

end program filter_test
