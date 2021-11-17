! This program is the driver of the tiles tests
program tiles_test_driver
use fruit
use mod_tiles_test, only: run_fruit_tiles
  implicit none

  ! init fruit suite
  call init_fruit
  call init_fruit_xml

  ! run the tiles test basket
  call run_fruit_tiles

  ! write test summary and finilize test suit
  call fruit_summary
  call fruit_summary_xml
  call fruit_finalize

end program tiles_test_driver
