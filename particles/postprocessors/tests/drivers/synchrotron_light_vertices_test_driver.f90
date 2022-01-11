! This program is the driver of the synchrotron light vertices tests
program synchrotron_light_vertices_test
use fruit
use mod_synchrotron_light_vertices_test, only: run_fruit_synchrotron_light_vertices
  implicit none

  ! init fruit suite
  call init_fruit
  call init_fruit_xml

  ! run the synchrotron light vertices test basket
  call run_fruit_synchrotron_light_vertices

  ! write test summary and finilize test suit
  call fruit_summary
  call fruit_summary_xml
  call fruit_finalize

end program synchrotron_light_vertices_test
