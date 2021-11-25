!> This program is the driver of the coordinate
!> transform tests
program coordinate_transforms_test_driver
use fruit
use mod_coordinate_transforms_test, only: run_fruit_coordinate_transforms
  implicit none

  ! init fruit suite
  call init_fruit
  call init_fruit_xml

  ! run the coordinate transform test basket
  call run_fruit_coordinate_transforms

  ! write test summary and finilize test suit
  call fruit_summary
  call fruit_summary_xml
  call fruit_finalize

end program coordinate_transforms_test_driver
