! This program is the driver of the pinhole lens tests
program pinhole_lens_test
use fruit
use mod_pinhole_lens_test, only: run_fruit_pinhole_lens
  implicit none

  ! init fruit suite
  call init_fruit
  call init_fruit_xml

  ! run the pinhole lens test basket
  call run_fruit_pinhole_lens

  ! write test summary and finilize test suit
  call fruit_summary
  call fruit_summary_xml
  call fruit_finalize

end program pinhole_lens_test
