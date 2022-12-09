! This program is the driver of the synchrotron light vertices tests
program full_synchrotron_light_dist_vertices_test
use fruit
use mod_full_synchrotron_light_dist_vertices_test, only: run_fruit_full_synchrotron_light_dist_vertices
  implicit none

  ! init fruit suite
  call init_fruit
  call init_fruit_xml

  ! run the full synchrotron light distribution vertices test basket
  call run_fruit_full_synchrotron_light_dist_vertices

  ! write test summary and finilize test suit
  call fruit_summary
  call fruit_summary_xml
  call fruit_finalize

end program full_synchrotron_light_dist_vertices_test
