! This program is the driver of the besselk tests
program besselk_test_driver
use fruit
use mod_besselk_test, only: run_fruit_besselk
  implicit none

  ! init fruit suite
  call init_fruit
  call init_fruit_xml

  ! run the besselk test basket
  call run_fruit_besselk

  ! write test summary and finilize test suit
  call fruit_summary
  call fruit_summary_xml
  call fruit_finalize

end program besselk_test_driver
