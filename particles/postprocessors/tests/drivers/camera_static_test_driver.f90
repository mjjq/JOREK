! This program is the driver of the camera static tests
program camera_static_test
use fruit
use mod_camera_static_test, only: run_fruit_camera_static
  implicit none

  ! init fruit suite
  call init_fruit
  call init_fruit_xml

  ! run the camera static test basket
  call run_fruit_camera_static

  ! write test summary and finilize test suit
  call fruit_summary
  call fruit_summary_xml
  call fruit_finalize

end program camera_static_test
