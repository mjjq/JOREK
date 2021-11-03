! This program is the driver of the besselik tests
program besselik_test_driver
use fruit
use mod_besselik_test_basket
  implicit none

  ! init fruit suite
	call init_fruit
	call init_fruit_xml

  ! run the test basket
	call fruit_basket

	! write test summary and finilize test suit
	call fruit_summary
	call fruit_summary_xml
	call fruit_finalize

end program besselik_test_driver
