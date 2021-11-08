! This program is the driver of the boost besselk tests
program boost_besselk_test_driver
use fruit
use mod_boost_besselk_test, only: run_fruit_boost_besselk
  implicit none

  ! init fruit suite
  call init_fruit
  call init_fruit_xml

  ! run then boost_besselk test basket
  call run_fruit_boost_besselk

  ! write test summary and finilize test suit
  call fruit_summary
  call fruit_summary_xml
  call fruit_finalize

end program boost_besselk_test_driver
