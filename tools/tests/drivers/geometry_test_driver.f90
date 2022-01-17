! This program is the driver of the geometry tests
program geometry_test_driver
use fruit
use mod_geometry_test, only: run_fruit_geometry
  implicit none

  ! init fruit suite
  call init_fruit
  call init_fruit_xml

  ! run the geometry test basket
  call run_fruit_geometry

  ! write test summary and finilize test suit
  call fruit_summary
  call fruit_summary_xml
  call fruit_finalize

end program geometry_test_driver
