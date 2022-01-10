! This program is the driver of the light vertices tests
program light_vertices_test
use fruit
use mod_light_vertices_test, only: run_fruit_light_vertices
  implicit none

  ! init fruit suite
  call init_fruit
  call init_fruit_xml

  ! run the light vertices test basket
  call run_fruit_light_vertices

  ! write test summary and finilize test suit
  call fruit_summary
  call fruit_summary_xml
  call fruit_finalize

end program light_vertices_test
