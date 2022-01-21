! This program is the driver of the camera tests
program camera_test
use fruit
use mod_camera_test, only: run_fruit_camera
  implicit none

  ! init fruit suite
  call init_fruit
  call init_fruit_xml

  ! run the camera test basket
  call run_fruit_camera

  ! write test summary and finilize test suit
  call fruit_summary
  call fruit_summary_xml
  call fruit_finalize

end program camera_test
