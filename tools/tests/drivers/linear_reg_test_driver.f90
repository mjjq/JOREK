! This program is the driver of the linear regression tests
program linear_reg_test_driver
use fruit
use mod_linear_reg_test, only: run_fruit_linear_reg
  implicit none

  ! init fruit suite
  call init_fruit
  call init_fruit_xml

  ! run the linear regression test basket
  call run_fruit_linear_reg

  ! write test summary and finilize test suit
  call fruit_summary
  call fruit_summary_xml
  call fruit_finalize

end program linear_reg_test_driver
