! This program is the driver of the math operators tests
program math_operators_test_driver
use fruit
use mod_math_operators_test, only: run_fruit_math_operators
  implicit none

  ! init fruit suite
  call init_fruit
  call init_fruit_xml

  ! run the math operators test basket
  call run_fruit_math_operators

  ! write test summary and finilize test suit
  call fruit_summary
  call fruit_summary_xml
  call fruit_finalize

end program math_operators_test_driver
